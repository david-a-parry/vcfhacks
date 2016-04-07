#!/usr/bin/env perl
use warnings;
use strict;
use Parallel::ForkManager;
use Getopt::Long;
use Sys::CPU;
use Pod::Usage;
use Term::ProgressBar;
use Data::Dumper;
use List::Util qw (sum);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use POSIX qw/strftime/;
use FindBin qw($RealBin);
use lib "$RealBin/lib";
use VcfReader;
use VcfhacksUtils;

my @samples;
my @evs;
my @evs_pop = qw ( EVS_EA_AF EVS_AA_AF EVS_ALL_AF ) ;
my $freq;
my $forks = 1;
my $buffer_size;
my $cpus = Sys::CPU::cpu_count();
my %opts = (samples => \@samples, esp_file => \@evs, freq => \$freq, cache => \$buffer_size, forks => \$forks);

GetOptions(\%opts,
        "output=s",
        "input=s",
        "samples=s{,}",
        "esp_file=s{,}",
        "dir=s",
        "help",
        "manual",
        "f|freq=f" => \$freq,
        "t|forks=i" => \$forks,
        "cache=i"      => \$buffer_size,
        "Progress",
        ) or pod2usage(-exitval => 2, -message => "Syntax error") ;
pod2usage (-verbose => 2) if $opts{manual};
pod2usage (-verbose => 1) if $opts{help};

pod2usage(-exitval => 2, -message => "Syntax error") if not $opts{input} or (not $opts{dir} and not @{$opts{esp_file}});
if (defined $freq){
        pod2usage(-exitval => 2, -message => "value for --freq argument must be greater than 0 and less than 50") if $freq > 50 or $freq <= 0;
        print STDERR "Filtering variant alleles with an allele frequency of $freq percent or above.\n";
}else{
    print STDERR "WARNING - no allele frequency (--freq) specified, matching alleles will be labelled but not filtered.\n";
}
$freq /= 100 if ($freq);
if ( $forks < 2 ) {
    $forks = 0;    #no point having overhead of forks for one fork
}else{
    if ($forks > $cpus){
        print STDERR "[Warning]: Number of forks ($forks) exceeds number of CPUs on this machine ($cpus)\n";
    }
    if ( not $buffer_size ) {
        $buffer_size = 10000 > $forks * 1000 ? 10000 : $forks * 1000;
    }
    print STDERR
"[INFO] Processing in batches of $buffer_size variants split among $forks forks.\n";
}


my $OUT;
if ($opts{output}){
        open ($OUT, ">$opts{output}") || die "Can't open $opts{output} for writing: $!";
}else{
        $OUT = *STDOUT;
}
my $time = strftime("%H:%M:%S", localtime);
my $total_vcf = -1;
print STDERR "[$time] Initializing input VCF...\n";
my ($header, $first_var, $VCF)  = VcfReader::getHeaderAndFirstVariant($opts{input});
die "Header not ok for input ($opts{input}) "
    if not VcfReader::checkHeader( header => $header );
if ( defined $opts{Progress} ) {
    $time = strftime( "%H:%M:%S", localtime );
    if (-p $opts{input} or $opts{input} eq "-" ) {
        print STDERR "\n[$time] INFO - Input is from STDIN or pipe - will report progress per 10000 variants. ";
    }else {
        $total_vcf= VcfReader::countVariants($opts{input});
        print STDERR "[$time] INFO - $opts{input} has $total_vcf variants.\n";
    }
}
my %sample_to_col = ();
if (@samples) {
    %sample_to_col = VcfReader::getSamples(
        get_columns => 1,
        header      => $header,
    );
}

$time = strftime( "%H:%M:%S", localtime );
print STDERR "[$time] Finished initializing input VCF\n";

if ($opts{dir}){
    opendir (my $DIR, $opts{dir}) or die "Can't read directory $opts{dir}: $!\n";
    my @files = grep {/\.vcf(\.gz)*$/}readdir($DIR);
    die "No vcf files found in $opts{dir}\n" if not @files;
    foreach my $file (@files){
        push @evs, "$opts{dir}/$file";
    }
}


#initialize EVS VCFs in parallel
my %evs_to_info = ();
my %evs_to_index = ();
my $dbpm = Parallel::ForkManager->new($forks);
for ( my $i = 0 ; $i < @evs ; $i++ ) {
    $dbpm->run_on_finish(    # called BEFORE the first call to start()
        sub {
            my ( $pid, $exit_code, $ident, $exit_signal, $core_dump,
                $data_structure_reference )
              = @_;

            if ( defined($data_structure_reference) ){                
                if ( ref $data_structure_reference eq 'ARRAY' ) {
                    $evs_to_info{$data_structure_reference->[2]} = $data_structure_reference->[0];
                    $evs_to_index{ $data_structure_reference->[2] } =
                      $data_structure_reference->[1];
                }else {
                    die
                      "Unexpected return from dbSNP reference initialization:\n"
                      . Dumper $data_structure_reference;
                }
            }else{ # problems occuring during storage or retrieval 
                die "No message received from child process $pid!\n";
            }
        }
    );
    $time = strftime( "%H:%M:%S", localtime );
    print STDERR "[$time] Initializing $evs[$i] EVS reference VCF "
      . ( $i + 1 ) . " of "
      . scalar(@evs) . "\n";
    $dbpm->start() and next;
    my @info_and_index = initializeDbsnpVcfs( $evs[$i] );
    push @info_and_index, $evs[$i];
    $time = strftime( "%H:%M:%S", localtime );
    print STDERR
      "[$time] Finished initializing $evs[$i] EVS reference VCF.\n";
    $dbpm->finish( 0, \@info_and_index );
}
print STDERR "Waiting for children...\n" if $opts{VERBOSE};
$dbpm->wait_all_children;

checkEvsInfoFields();

writeHeaders(); 

$time = strftime("%H:%M:%S", localtime);
print STDERR "[$time] EVS filter starting\n";

my $prog_total;
my $progressbar;
my $next_update = 0;
my $prev_percent = 0;
if ( defined $opts{Progress} and $total_vcf ) {
    $prog_total = $total_vcf * 3;
    $progressbar = Term::ProgressBar->new(
        #{ name => "Filtering", count => ($prog_total), ETA => "linear", } );
        { name => "Filtering", 
          count => ($prog_total), 
          ETA => "linear",
        } 
    );
}
my $kept = 0; #variants not filtered
my $filtered = 0; #variants filtered
my $found = 0; #variants that match a known SNP

my $n                = 0;
my $vars             = 0;
my @lines_to_process = ();
my %no_fork_args = ();
if ($forks < 2){
    foreach my $e (@evs) {
        my %s = VcfReader::getSearchArguments( $e, $evs_to_index{$e} );
        $no_fork_args{$e} = \%s;
    }
}
processLine($first_var);
VAR: while (my $line = <$VCF> ) {
    processLine($line);
}

process_buffer() if $forks > 1;
close $VCF;
close $OUT;
if ($progressbar){
    $time = strftime("%H:%M:%S", localtime);
    $progressbar->message( "[INFO - $time] $vars variants processed" );
    $progressbar->update($prog_total) if $prog_total >= $next_update;
}

$time = strftime("%H:%M:%S", localtime);
print STDERR "Time finished: $time\n";
print STDERR "$found matching variants identified.\n";
print STDERR "$filtered variants filtered, $kept variants retained.\n";


################################################
#####################SUBS#######################
################################################
sub processLine{
    my $line = shift;
    return if $line =~ /^#/;
    $n++;
    $vars++;
    checkProgress(1);
    if ($forks > 1){
        push @lines_to_process, $line;
        if ( @lines_to_process >= $buffer_size ) {
            process_buffer();
            @lines_to_process = ();
        }
    }else{
        chomp $line;
        #our VcfReader methods should be more efficient on pre-split lines
        my @split_line = split( "\t", $line );
        my %res = filter_on_evs_maf( \@split_line, \%no_fork_args);
        if ($res{keep}){
            print $OUT join("\t", @{$res{keep}}) ."\n"; 
            $kept++;
        }
        $kept++ if $res{keep};
        $found++ if $res{found};
        $filtered++ if $res{filter};
        $n += 2;
        checkProgress();
    }
}

################################################
sub checkProgress{
    return if not $progressbar;
    my $do_count_check = shift;
    if ($prog_total > 0){
        $next_update = $progressbar->update($n) if $n >= $next_update;
    }elsif($do_count_check){#input from STDIN/pipe
        if (not $vars % 10000) {
            my $time = strftime( "%H:%M:%S", localtime );
            $progressbar->message( "[INFO - $time] $vars variants read" );
        }
    }
}

################################################
sub process_buffer {

    #this arrays is an array of refs to batches
    # so that we can quickly sort our batches rather than
    # performing a big sort on all lines
    return if not @lines_to_process;
    my @lines_to_print;
    my $i               = 0;
    my $t               = 0;
    my $lines_per_slice = @lines_to_process;
    if ( $forks > 0 ) {
        $lines_per_slice = int( @lines_to_process / $forks ) > 1 ? int( @lines_to_process / $forks ) : 1;
    }
    my @batch = ();

    #get a batch for each thread
    for ( my $i = 0 ; $i < @lines_to_process ; $i += $lines_per_slice ) {
        my $last =
          ( $i + $lines_per_slice - 1 ) < $#lines_to_process
          ? $i + $lines_per_slice - 1
          : $#lines_to_process;
        if ($i + $lines_per_slice >= @lines_to_process){
            $last = $#lines_to_process;
        }
        my @temp = @lines_to_process[ $i .. $last ];
        push @batch, \@temp;
    }
    my $pm = Parallel::ForkManager->new($forks);
    $pm->run_on_finish(    # called BEFORE the first call to start()
        sub {
            my ( $pid, $exit_code, $ident, $exit_signal, $core_dump,
                $data_structure_reference )
              = @_;

            if ( defined($data_structure_reference) )
            {
                my %res = %{$data_structure_reference}
                  ;        
                if ( ref $res{keep} eq 'ARRAY' ) {
                    $lines_to_print[$res{order}] = \@{ $res{keep} } if @{ $res{keep} };
                }
                if ( $res{found} ) {
                    $found += $res{found} ;
                }
                if ( $res{filter} ) {
                    $filtered += $res{filter} ;
                }
                $n += $res{batch_size};
                checkProgress();
            }
            else
            { 
                die "ERROR: no message received from child process $pid!\n";
            }
        }
    );
    my $order = -1;
    foreach my $b (@batch) {
        $order++;
        $pm->start() and next;
        my %results = process_batch($b, $order);
        $pm->finish( 0, \%results );
    }
    $pm->wait_all_children;

    #print them
    if (@lines_to_print){
        my $incr_per_batch = @lines_to_process / @lines_to_print;
        foreach my $batch (@lines_to_print) {
            if (not defined $batch){
                $n += $incr_per_batch;
                checkProgress();
                next;
            }
            my $incr_per_line = $incr_per_batch / @$batch;
            foreach my $l (@$batch) {
                if ( ref $l eq 'ARRAY' ) {
                    print $OUT join( "\t", @$l ) . "\n";
                }
                else {
                    print $OUT "$l\n";
                }
                $kept++;
                $n += $incr_per_line;
                checkProgress();
            }
        }
    }else{
        if ($progressbar) {
            $n += @lines_to_process;
            checkProgress();
        }

    }
}

################################################
sub process_batch {

    #filter a set of lines
    my ($batch, $order) = @_;
    my %results = 
    ( 
        batch_size => scalar(@$batch),
        order      => $order,
    );
    my %sargs;
    foreach my $e (@evs) {

        #WE'VE FORKED AND COULD HAVE A RACE CONDITION HERE -
        # WHICH IS WHY WE DO A FIRST PASS WITH OUR initializeDbsnpVcfs
        # METHOD TOWARDS THE START OF THE PROGRAM
        my %s = VcfReader::getSearchArguments( $e, $evs_to_index{$e} );
        $sargs{$e} = \%s;
    }
    foreach my $line ( @{$batch} ) {
        chomp $line;
        my @split_line = split( "\t", $line );
        my %res = filter_on_evs_maf( \@split_line, \%sargs );
        push @{ $results{keep} },   $res{keep}   if $res{keep};
        $results{filter}++  if $res{filter};
        $results{found}++  if $res{found};
    }
    foreach my $d (keys %sargs){
        if (exists $sargs{$d}->{file_handle}){
            close $sargs{$d}->{file_handle} if $sargs{$d}->{file_handle};
        }
    }
    return %results;
}


################################################
sub filter_on_evs_maf{
    my ( $vcf_line, $search_args ) = @_;
    my %r = (keep => undef, found => 0, filter => 0);
    if (defined $opts{Progress}){
           $next_update = $progressbar->update($n) if $n >= $next_update;
    }
    my %min_vars       = VcfReader::minimizeAlleles($vcf_line);
    my %sample_alleles = ();
    if (@samples){
        %sample_alleles = map { $_ => undef } VcfReader::getSampleCall(
            multiple            => \@samples,
            sample_to_columns   => \%sample_to_col,
            return_alleles_only => 1,
            line                => $vcf_line
        );
    }
ALLELE: foreach my $allele ( sort { $a <=> $b } keys %min_vars ) {
        if ($min_vars{$allele}->{ORIGINAL_ALT} eq '*'){
            $min_vars{$allele}->{filter_snp} = 1;
            next;
        }
        $min_vars{$allele}->{filter_snp} = 0;
        if (@samples){
            #doesn't exist in any sample...
            next if not exists $sample_alleles{$allele};
        }
        foreach my $k ( keys %{$search_args} ) {
            if (
                my @snp_hits = VcfReader::searchForPosition(
                    %{ $search_args->{$k} },
                    chrom => $min_vars{$allele}->{CHROM},
                    pos   => $min_vars{$allele}->{POS}
                )
              )
            {
                foreach my $snp_line (@snp_hits) {

                    #check whether the snp line(s) match our variant
                    my @snp_split = split( "\t", $snp_line );
                    if ( my $match = checkVarMatches( $min_vars{$allele}, \@snp_split ) ) {
                        $r{found} = 1;
                        annotateEvsAlleleFrequencies($min_vars{$allele}, \@snp_split, $match);#this annotates allele frequencies per allele and adds {filter_snp} annotation if >= $freq
                    }
                }
            }
        }
    }
    #filter if ALL ALT alleles have been found to have MAF above $freq
    if ($freq){
        my $al_filter = 0;
FILTER: foreach my $allele (keys %min_vars){
            if (not $min_vars{$allele}->{filter_snp}){
                last FILTER;
            }
            $al_filter++;
        }
        if ($al_filter == keys %min_vars){
            $r{filter} = 1;
            return %r;
        }
    }
    #not all alleles' MAF above $freq, annotate information and return
    my %evs_info = ();
    #create hash of MAF values, one value per allele
    foreach my $allele ( sort { $a <=> $b } keys %min_vars){
        foreach my $pop (@evs_pop){
            if (exists $min_vars{$allele}->{$pop}){
                push @{$evs_info{$pop}}, $min_vars{$allele}->{$pop};
            }else{
                push @{$evs_info{$pop}}, '.';
            }
        }
    }
    #add EVS MAF values to INFO field
    foreach my $pop (@evs_pop){
        $vcf_line = VcfReader::addVariantInfoField 
            (
                line => $vcf_line,
                id   => $pop,
                value => join(",", @{ $evs_info{$pop} } ), 
            );
    }
    $r{keep} = $vcf_line;
    return %r;
}

################################################
sub annotateEvsAlleleFrequencies{
    my ($min_allele, $snp_line, $snp_allele) = @_;
    
    $min_allele->{EVS_EA_AF}  = getAf($snp_line, $snp_allele, "EA_AC");
    $min_allele->{EVS_AA_AF}  = getAf($snp_line, $snp_allele, "AA_AC");
    $min_allele->{EVS_ALL_AF} = getAf($snp_line, $snp_allele, "TAC");
    if ($freq){
        foreach my $af (qw / EVS_EA_AF EVS_AA_AF EVS_ALL_AF / ){
            next if $min_allele->{$af} eq '.';
            $min_allele->{filter_snp}++ if $min_allele->{$af} >= $freq;
        }
    }
}
################################################
sub getAf{
    my ($snp_line, $alt, $field) = @_;
    my $counts = VcfReader::getVariantInfoField( $snp_line, $field);
    my @af = split(",", $counts); 
    if (@af < $alt){
        die "Not enough allele frequencies found for allele counts!\n";
    }
    my $total = sum (@af);
    return '.' if not $total;
    #EVS VCF ACs have the alt alleles first, then the REF allele
    return sprintf("%g", $af[$alt -1] / $total) ;
}
################################################
sub checkVarMatches {
    my ( $min_allele, $snp_line ) = @_;
    my $matches = 0;
    my %snp_min = VcfReader::minimizeAlleles($snp_line);
    foreach my $snp_allele ( keys %snp_min ) {
        $min_allele->{CHROM} =~ s/^chr//;
        $snp_min{$snp_allele}->{CHROM} =~ s/^chr//;
        next if $min_allele->{CHROM} ne $snp_min{$snp_allele}->{CHROM};
        next if $min_allele->{POS} ne $snp_min{$snp_allele}->{POS};
        next if $min_allele->{REF} ne $snp_min{$snp_allele}->{REF};
        next if $min_allele->{ALT} ne $snp_min{$snp_allele}->{ALT};
        return $snp_allele;
    }
    return 0;
}
        
################################################
sub initializeDbsnpVcfs {
    my ($snpfile) = @_;
    #we getSearchArguments here simply to prevent a race condition later
    my @head = VcfReader::getHeader($snpfile);
    die "Header not ok for $snpfile " 
        if not VcfReader::checkHeader( header => \@head );
    #my %sargs = VcfReader::getSearchArguments($snpfile);
    my %info = VcfReader::getInfoFields( header => \@head);
    my %index = VcfReader::readIndex($snpfile);
    return ( \%info, \%index );
}
################################################
sub checkEvsInfoFields{
    foreach my $field ( qw / EA_AC AA_AC TAC / ){ 
        foreach my $d (@evs){
            $time = strftime( "%H:%M:%S", localtime );
            if ( not $evs_to_info{$d}->{$field} ){
                die "Could not find required EVS field ($field) in file: $d\n";
            }
        }
    }
}

#####################################
sub writeHeaders{
    my $meta_head = join("\n", grep {/^##/} @$header);
    print $OUT "$meta_head\n";
    my %ea_af = 
    (
        ID          => "EVS_EA_AF",
        Number      => "A",
        Type        => "Float",
        Description => "European American allele frequencies calculated ".
                       "from NHLBI EVS VCFs"
    );
    
    my %aa_af = 
    (
        ID          => "EVS_AA_AF",
        Number      => "A",
        Type        => "Float",
        Description => "African American allele frequencies calculated ".
                       "from NHLBI EVS VCFs",
    );

    my %all_af = 
    (
        ID          => "EVS_ALL_AF",
        Number      => "A",
        Type        => "Float",
        Description => "Total allele frequencies calculated from NHLBI ".
                       "EVS VCFs",
    );
    foreach my $inf_hash( \%ea_af, \%aa_af, \%all_af){
        print $OUT VcfhacksUtils::getInfoHeader(%{$inf_hash}) . "\n";
    }
    my $optstring = VcfhacksUtils::getOptsVcfHeader(%opts);
    print $OUT "$optstring\n$header->[-1]\n"; 
}

#####################################

=head1 NAME

filterOnEvsMaf.pl - annotate/filter NHLBI ESP variants from a VCF file 

=head1 SYNOPSIS

        filterOnEvsMaf.pl -i [vcf file] -d [directory containing ESP VCF files] [options]
        filterOnEvsMaf.pl -i [vcf file] -e [ESP VCF file(s)] [options]
        filterOnEvsMaf.pl -h (display help message)
        filterOnEvsMaf.pl -m (display manual page)

=cut

=head1 ARGUMENTS

=over 8

=item B<-i    --input>

Input VCF file 

=item B<-o    --output>

Output snp filtered file. Optional - default is STDOUT.

=item B<-d    --dir>

Directory containing ESP VCF files. These can be downloaded from http://evs.gs.washington.edu/EVS/.

=item B<-e    --esp_file>

One or more ESP VCF files to use. You may use this in addition to or instead of --dir argument, but for most uses it is envisaged that you will simply use the --dir argument to specify the location of a directory containing ESP files for each chromosome.

=item B<-s    --samples>

One or more samples to check variants for.  Default is to check all variants specified by the ALT field. If specified, variants will not be filtered unless variant is present in at least one of the samples specified by this argument.

=item B<-f    --freq>

Percent minor allele frequency to filter from. Only variants with minor alleles equal to or over this frequency will be removed. Default is to not filter variants, only annotate.

=item B<-t    --forks>

Number of forks to create for parallelising your analysis. By default no forking is done. To speed up your analysis you may specify the number of parallel processes to use here. (N.B. forking only occurs if a value of 2 or more is given here as creating 1 fork only results in increased overhead with no performance benefit).

=item B<-c    --cache>

Cache size. Variants are processed in batches to allow for efficient parallelisation. When forks are used the default is to process up to 10,000 variants at once or 1,000 x no. forks if more than 10 forks are used. If you find this program comsumes too much memory when forking you may want to set a lower number here. When using forks you may get improved performance by specifying a higher cache size, however the increase in memory usage is proportional to your cache size multiplied by the number of forks.

=item B<-p    --progress>

Use this flag to show a progress bar while this program is running.

=item B<-h    --help>

Display help message.

=item B<-m    --manual>

Show manual page.

=back 

=cut


=head1 DESCRIPTION

This program will annotate variants in a VCF file with NHLBI ESP allele frequencies and optionally filter variants with an allele frequency equal to or greater than a user supplied allele frequency (as specified with --freq). The following INFO fields are added to variants: 

    EVS_EA_AF  (European-American allele frequency)
    EVS_AA_AF  (African-American allele frequency)
    EVS_ALL_AF (Allele frequency for all EVS samples)
    

ESP VCF files are available from http://evs.gs.washington.edu/EVS/.

=cut

=head1 AUTHOR

David A. Parry

=head1 COPYRIGHT AND LICENSE

Copyright 2013,2014,2015  David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut

