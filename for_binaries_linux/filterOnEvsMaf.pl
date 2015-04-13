#!/usr/bin/perl
use warnings;
use strict;
use Parallel::ForkManager;
use Getopt::Long;
use Sys::CPU;
use Pod::Usage;
use Term::ProgressBar;
use Data::Dumper;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use POSIX qw/strftime/;
use FindBin;
use lib "$FindBin::Bin/lib";
use VcfReader;
my @samples;
my @evs;
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
    print STDERR "WARNING - no allele frequency (--freq) specified, ANY matching alleles will be filtered.\n";
}

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
my $total_vcf;
print STDERR "[$time] Initializing input VCF... ";
die "Header not ok for input ($opts{input}) "
    if not VcfReader::checkHeader( vcf => $opts{input} );
if ( defined $opts{Progress} ) {
    if ( $opts{input} eq "-" ) {
        print STDERR "Can't use --progress option when input is from STDIN\n";
        delete $opts{Progress};
    }
    else {
        $total_vcf= VcfReader::countVariants($opts{input});
        print STDERR "$opts{input} has $total_vcf variants. ";
    }
}
my %contigs = ();
if ($forks){
    %contigs       = VcfReader::getContigOrder( $opts{input} );
}
my %sample_to_col = ();
if (@samples) {
    %sample_to_col = VcfReader::getSamples(
        vcf         => $opts{input},
        get_columns => 1,
    );
}

$time = strftime( "%H:%M:%S", localtime );
print STDERR "\n[$time] Finished initializing input VCF\n";

if ($opts{dir}){
    opendir (my $DIR, $opts{dir}) or die "Can't read directory $opts{dir}: $!\n";
    my @files = grep {/\.vcf(\.gz)*$/}readdir($DIR);
    die "No vcf files found in $opts{dir}\n" if not @files;
    foreach my $file (@files){
        push @evs, "$opts{dir}/$file";
    }
}


#initialize EVS VCFs in parallel
my @evs_headers = ();
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
                    push @evs_headers, @{ $data_structure_reference->[0] };
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
    my @head_and_index = initializeDbsnpVcfs( $evs[$i] );
    push @head_and_index, $evs[$i];
    $time = strftime( "%H:%M:%S", localtime );
    print STDERR
      "[$time] Finished initializing $evs[$i] EVS reference VCF.\n";
    $dbpm->finish( 0, \@head_and_index );
}
print STDERR "Waiting for children...\n" if $opts{VERBOSE};
$dbpm->wait_all_children;


if ($freq){
    if (not grep {/##INFO=<ID=MAF/} @evs_headers ){
        die  "Can't find MAF fields in evs file headers.\n";
    }
}
my $meta_head = VcfReader::getMetaHeader( $opts{input} );
print $OUT "$meta_head\n";
print $OUT  "##filterOnEvsMaf=\"";
my @opt_string = ();
foreach my $k (sort keys %opts){
    if (not ref $opts{$k}){
        push @opt_string, "$k=$opts{$k}";
    }elsif (ref $opts{$k} eq 'SCALAR'){
        if (defined ${$opts{$k}}){
            push @opt_string, "$k=${$opts{$k}}";
        }else{
            push @opt_string, "$k=undef";
        }
    }elsif (ref $opts{$k} eq 'ARRAY'){
        if (@{$opts{$k}}){
            push @opt_string, "$k=" .join(",", @{$opts{$k}});
        }else{
            push @opt_string, "$k=undef";
        }
    }
}
my $col_header = VcfReader::getColumnHeader( $opts{input} );
print $OUT join( " ", @opt_string ) . "\"\n" . $col_header . "\n";

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
my @lines_to_process = ();
my $VCF              = VcfReader::_openFileHandle( $opts{input} )
  ;    # or die "Can't open $opts{input} for reading: $!\n";
my %no_fork_args = ();
if ($forks < 2){
    foreach my $e (@evs) {
        my %s = VcfReader::getSearchArguments( $e, $evs_to_index{$e} );
        $no_fork_args{$e} = \%s;
    }
}
VAR: while ( my $line = <$VCF> ) {
    next if $line =~ /^#/;
    $n++;
    if ($progressbar) {
        $next_update = $progressbar->update($n) if $n >= $next_update;
    }
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
        print $OUT "$line\n" if $res{keep};
        $kept++ if $res{keep};
        $found++ if $res{found};
        $filtered++ if $res{filter};
        $n += 2;
        if ($progressbar) {
            $next_update = $progressbar->update($n) if $n >= $next_update;
        }
    }
}
close $VCF;
process_buffer() if $forks > 1;
close $OUT;
if ( defined $opts{Progress} ) {
    $progressbar->update($prog_total) if $prog_total >= $next_update;
}

$time = strftime("%H:%M:%S", localtime);
print STDERR "Time finished: $time\n";
print STDERR "$found matching variants identified.\n";
print STDERR "$filtered variants filtered, $kept variants retained.\n";


################################################
#####################SUBS#######################
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
                    push @lines_to_print, \@{ $res{keep} } if @{ $res{keep} };
                }
                if ( $res{found} ) {
                    $found += $res{found} ;
                }
                if ( $res{filter} ) {
                    $filtered += $res{filter} ;
                }
                if ($progressbar) {
                    $n += $res{batch_size};
                    $next_update = $progressbar->update($n)
                      if $n >= $next_update;
                }
            }
            else
            { 
                die "ERROR: no message received from child process $pid!\n";
            }
        }
    );
    foreach my $b (@batch) {
        $pm->start() and next;
        my %results = process_batch($b);
        $pm->finish( 0, \%results );
    }
    $pm->wait_all_children;

    #print them
    @lines_to_print = sort {
        VcfReader::by_first_last_line($a, $b, \%contigs) 
        } @lines_to_print;


    if (@lines_to_print){
        my $incr_per_batch = @lines_to_process / @lines_to_print;
        foreach my $batch (@lines_to_print) {
            my $incr_per_line = $incr_per_batch / @$batch;
            foreach my $l (@$batch) {
                if ( ref $l eq 'ARRAY' ) {
                    print $OUT join( "\t", @$l ) . "\n";
                }
                else {
                    print $OUT "$l\n";
                }
                $kept++;
                if ($progressbar) {
                    $n += $incr_per_line;
                    $next_update = $progressbar->update($n) if $n >= $next_update;
                }
            }
        }
    }else{
        if ($progressbar) {
            $n += @lines_to_process;
            $next_update = $progressbar->update($n) if $n >= $next_update;
        }

    }
}

################################################
sub process_batch {

    #filter a set of lines
    my ($batch) = @_;
    my %results = ( batch_size => scalar(@$batch) );
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
                    if ( checkVarMatches( $min_vars{$allele}, \@snp_split ) ) {
                        $r{found} = 1;
                        #perform filtering if freq is set
                        if ($freq){
                            my $maf_field = VcfReader::getVariantInfoField( \@snp_split, 'MAF');
                            die "Expected MAF field to be a comma separated list of 3 numbers for $k line: $snp_line " 
                                if $maf_field !~/(\d+\.\d+,){2,}\d+\.\d+/;
                            my @mafs = split(',', $maf_field);
                            foreach my $maf (@mafs){
                                if ($maf >= $freq){
                                    $min_vars{$allele}->{filter_snp}++ ;
                                    next ALLELE;
                                }
                            }
                        }else{#default behaviour is to filter everything that matches
                            $min_vars{$allele}->{filter_snp}++ ;
                            next ALLELE;
                        }
                    }
                }
            }
        }
    }
    foreach my $allele (keys %min_vars){
        if (not $min_vars{$allele}->{filter_snp}){
            #print if ANY of the alleles haven't met 
            #our filter criteria
            $r{keep} = $vcf_line;
            return %r;
        }
    }
    #if we're here all alleles must have met filter criteria
    $r{filter} = 1;
    return %r;
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
        $matches++;
    }
    return $matches;
}
        
################################################
sub initializeDbsnpVcfs {
    my ($snpfile) = @_;
    #we getSearchArguments here simply to prevent a race condition later
    my @head = VcfReader::getHeader($snpfile);
    die "Header not ok for $snpfile " 
        if not VcfReader::checkHeader( header => \@head );
    #my %sargs = VcfReader::getSearchArguments($snpfile);
    my %index = VcfReader::readIndex($snpfile);
    return ( \@head, \%index );
}
################################################
=head1 NAME

filterOnEvsMaf - filter NHLBI ESP variants from a VCF file 

=head1 SYNOPSIS

        filterOnEvsMaf -i [vcf file] -d [directory containing ESP VCF files] [options]
        filterOnEvsMaf -i [vcf file] -e [ESP VCF file(s)] [options]
        filterOnEvsMaf -h (display help message)
        filterOnEvsMaf -m (display manual page)

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

Percent minor allele frequency to filter from. Only variants with minor alleles equal to or over this frequency will be removed. Default is to remove any matching variant.

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

This program will filter a VCF file to remove variants matching variants in VCF files from NHLBI ESP.  If --freq is set only variants in ESP with a MAF equal to or greater than --freq will be filtered. ESP VCF files are available from http://evs.gs.washington.edu/EVS/.

=cut

=head1 AUTHOR

David A. Parry

University of Leeds

=head1 COPYRIGHT AND LICENSE

Copyright 2013,2014  David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut

