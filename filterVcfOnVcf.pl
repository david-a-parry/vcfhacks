#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Parallel::ForkManager;
use Sys::CPU;
use Pod::Usage;
use Term::ProgressBar;
use Data::Dumper;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use POSIX qw(strftime);
use FindBin;
use lib "$FindBin::Bin/lib";
use VcfReader;

=head1 NAME

filterVcfOnVcf.pl - filter a VCF file using variants from one or more other VCF files.

=head1 SYNOPSIS

filterVcfOnVcf.pl -i <input vcf file> -f <vcfs to use to filter input> [options]

filterVcfOnVcf.pl -i <input vcf file> -d <directory containing vcfs to use to filter input> [options]

filterVcfOnVcf.pl -h (show help message)

filterVcfOnVcf.pl -m (show manual page)


=head1 ARGUMENTS

=over 8 

=item B<-i    --input>

Input VCF file.

=item B<-o    --output>

Output file. Optional - defaults to STDOUT.

=item B<-f    --filter>

One or more VCF files to use as to filter input. Variants that match from the input VCF will be filtered from the output. By default, if any sample matches a variant it will be filtered but combinations of --reject, --not_samples, --allele_frequency_filter, --threshold and --filter_homozygotes can be used to modify this behaviour.

=item B<-d    --directories>

One or more directories containing VCF files to use to filter input. Only files with a ".vcf" or ".vcf.gz" extension will be used.

=item B<-e    --expression>

Perl style regular expression to use to match VCF files when using --directories argument.

=item B<-s    --samples>

Samples from input file to check. If specified only variant alleles from these samples will be used to compare against variants in filter VCFs. Any line in your input without a variant call in these samples will be filtered. Default is to look at all alleles for each variant.

=item B<-r    --reject>

Samples from filter VCF files to use for filtering.  If specified only variant alleles from these samples will be used to compare with the input VCF. Default is to look at all alleles.

=item B<-x    --not_samples>

Samples from filter VCF files to ignore.  If specified variant alleles from all samples except these will be used to compare with the input VCF. Default is to look at all alleles.

=item B<-w    --info_filter>

Use this flag to use metrics written to the INFO field of your filter VCFs for filtering on frequency, homozygosity or threshold counts rather than checking sample genotypes. Filtering on allele frequency (see --allele_frequency_filter option below) can be performed as long as your VCF contains allele frequency (AF) INFO tags and optionally population specific allele count (AC) and allele number (AN) tags (see --population_ids option below). This option is particularly useful for filtering against VCFs from ExAC using global and population allele frequencies. Alternatively, you can pre-process your own filter VCFs with the sampleCallsToInfo.pl script to reduce genotype information to AF, PGTS and GTC fields. Converting filter VCFs containing many hundreds/thousands of samples  with sampleCallsToInfo.pl for use with this option will speed up your analyses by many orders of magnitude.

=item B<-y    --allele_frequency_filter>

Reject variants if the allele frequency in all --filter VCFs is equal to or greater than this value. All samples will be used to calculate allele frequency regardless of --reject or --not_samples settings, although variants will only be counted if the genotype quality is greater than the value given for --un_quality. The value must be a float between 0 and 1. 

=item B<--population_ids>

If using the --info_filter option, --allele_frequency_filter option and only one filter VCF, by default the filter VCF will be scanned for allele counts (e.g. AFR_AC) and allele number (e.g. AFR_AN) INFO tags for the population codes as used by EXAC (AFR, AMR, EAS, FIN, NFE, SAS). An alternative list of population IDs to search can be entered here. To disable this feature use this option specifying 'disable' as an argument.
    
=item B<--minimum_allele_number>

If filtering using population allele counts (see --population_ids above) use this value to specify a minimum number of alleles to be called for a population in order to filter a variant on that population's allele frequency. Default = 100.

=item B<-l    --threshold>

Filter lines if the number of samples containing at least one matching variant allele is equal to or greater than this value.

=item B<-z    --filter_homozygotes>

Filter lines if a sample is homozygous for a variant. If --reject or --not_samples options are used only the relevant samples will be checked. If this option is used with --allele_frequency_filter or --threshold options lines will be filtered if the variant is homozygous in a sample regardless of allele frequency; if there is no homozygous sample lines will be filtered n frequency as normal. If used without --allele_frequency_filter or --threshold options lines will only be filtered if the variant is homozygous in at least one sample.

=item B<-q    --quality>

Minimum variant quality as defined by the QUAL field.  Variants in the input will only be printed if they have an equal or higher variant quality score. Variants in filter VCFs will only be used for filtering if they have an equal or higher variant quality score.

=item B<-g    --genotype_quality>

Minimum genotype quality.  Only applicable if --samples,  --reject or --threshold arguments are used. Sample alleles will only be counted if they have a genotype quality score equal to or higher than this value. Default is 0 (i.e. no filtering on genotype quality).

=item B<-a    --aff_quality>

Minimum genotype qualities to consider for samples specified by --samples argument only. Any sample call with a genotype quality below this threshold will be considered a no call. Default is 0.

=item B<-u    --un_quality>

Minimum genotype qualities to consider for samples specified by --reject argument only. Any sample call with a genotype quality below this threshold will be considered a no call. Default is 0.

=item B<-p    --print_matching>

Use this flag to reverse the script so that only matching lines are printed.

=item B<-t    --forks>

Number of forks to create for parallelising your analysis. By default no forking is done. To speed up your analysis you may specify the number of parallel processes to use here. (N.B. forking only occurs if a value of 2 or more is given here as creating 1 fork only results in increased overhead with no performance benefit).

=item B<-c    --cache>

Cache size. Variants are processed in batches to allow for efficient parallelisation. When forks are used the default is to process up to 10,000 variants at once or 1,000 x no. forks if more than 10 forks are used. If you find this program comsumes too much memory when forking you may want to set a lower number here. When using forks you may get improved performance by specifying a higher cache size, however the increase in memory usage is proportional to your cache size multiplied by the number of forks.

=item B<-b    --progress>

Show a progress bar.

=item B<-h    --help>

Show help message.

=item B<-m    --manual>

Show manual page.

=back

=head1 DESCRIPTION

Filter variants from a VCF using one or more other VCFs.  You may specify samples or a number of samples with matching variants to filter with.  You can also use minimum genotype or variant qualities to filter with. Alternatively, filtering can be performed on certain annotations in the INFO field of each variant using the --info_filter option.

By default, each variant from the given --input file will be assessed and if all alternative/variant alleles for a given variant are represented in the --filter VCFs provided, the variant will be filtered. You may use --reject or --not_samples arguments to only filter variant alleles if present in specific samples in the --filter VCFs. You may also use the --samples argument to only compare variants from specific samples in your --input VCF and filter any variants that do not have alleles represented by these samples. See above for details of all available options. 

=head1 EXAMPLES

    filterVcfOnVcf.pl -i input.vcf -f controls.vcf -o filtered.vcf
    (filter variants in input.vcf if present in controls.vcf)
    
    filterVcfOnVcf.pl -i input.vcf -f controls.vcf -o filtered.vcf -y 0.01
    (filter variants in input.vcf if present at an allele frequency of 1% or higher in controls.vcf)

    filterVcfOnVcf.pl -i input.vcf -f ExAC.r0.3.sites.vep.vcf.gz -o filtered.vcf -y 0.01  -w 
    (filter variants in input.vcf if present at an allele frequency of 1% or higher in ExAC using INFO fields for filtering)


=head1 AUTHOR

David A. Parry

=head1 COPYRIGHT AND LICENSE

Copyright 2013, 2014, 2015  David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut

my $vcf;
my @filter_vcfs;
my @dirs;
my @samples
  ; #samples in vcf input to check calls for - filter only if they contain same allele as in filter_vcf. Default is to simply check alleles in ALT field and ignore samples.
my @reject
  ;    #if specified will only check alleles for these samples in filter_vcfs
my @ignore_samples;    #these samples will be ignored in either VCF
my $filter_with_info
  ; #use allele counts/genotype counts in INFO field to filter with, not samples
my $threshold =
  0;    #only filter if we see the allele this many times in filter_vcfs
my $filter_homozygotes
  ;     #flag to filter if any of the reject samples are homozygous
my $maf            = 0;
my $print_matching = 0
  ; #flag telling the script to invert so we print matching lines and filter non-matching lines
my $out;
my $min_qual = 0;
my $minGQ    = 0;
my $aff_quality;      #will convert to $minGQ value if not specified
my $unaff_quality;    #will convert to $minGQ value if not specified
my $help;
my $man;
my $progress;
my $regex;            #match this regex if looking in dir
my @pop_acs = ();
    #if we have only one vcf and find population allele counts and allele numbers
    #we'll use those to calculate population specific frequencies and filter on MAF 
    #if --allele_frequency_filter argument is specified
my $MIN_AN = 100; 
    # min number of alleles required before @pop_acs values are used for filtering
my @pop_ids =  ();
    # population IDs to look for for @pop_acs
my $forks = 0;
my $buffer_size;
my %opts = (
    not_samples             => \@ignore_samples,
    expression              => \$regex,
    input                   => \$vcf,
    output                  => \$out,
    filter                  => \@filter_vcfs,
    directories             => \@dirs,
    samples                 => \@samples,
    reject                  => \@reject,
    info_filter             => \$filter_with_info,
    allele_frequency_filter => \$maf,
    threshold               => \$threshold,
    filter_homozygotes      => \$filter_homozygotes,
    quality                 => \$min_qual,
    genotype_quality        => \$minGQ,
    aff_quality             => \$aff_quality,
    un_quality              => \$unaff_quality,
    progress                => \$progress,
    population_ids          => \@pop_ids,
    cache                   => \$buffer_size,
    forks                   => \$forks,
    help                    => \$help,
    manual                  => \$man,
    print_matching          => \$print_matching
);

GetOptions(
    \%opts,
    'x|not_samples=s{,}'          => \@ignore_samples,
    'expression=s'                => \$regex,
    'i|input=s'                   => \$vcf,
    'output=s'                    => \$out,
    'f|filter=s{,}'               => \@filter_vcfs,
    'directories=s{,}'            => \@dirs,
    'samples=s{,}'                => \@samples,
    'reject=s{,}'                 => \@reject,
    'w|info_filter'               => \$filter_with_info,
    'y|allele_frequency_filter=f' => \$maf,
    'l|threshold=i'               => \$threshold,
    'z|filter_homozygotes'        => \$filter_homozygotes,
    'quality=f'                   => \$min_qual,
    'genotype_quality=f'          => \$minGQ,
    'a|aff_quality=i'             => \$aff_quality,
    'un_quality=i'                => \$unaff_quality,
    'population_ids=s{,}'         => \@pop_ids,
    'minimum_allele_number=i'     => \$MIN_AN,
    'bar|progress'                => \$progress,
    "t|forks=i"                   => \$forks,
    "cache=i"                     => \$buffer_size,
    'help'                        => \$help,
    'm|manual'                    => \$man,
    'p|print_matching'            => \$print_matching
  )
  or pod2usage(
    -message => "Syntax error.",
    -exitval => 2
  );

pod2usage( -verbose => 2 ) if $man;
pod2usage( -verbose => 1 ) if $help;
pod2usage( -message => "Syntax error", exitval => 2 )
  if ( not $vcf or ( not @filter_vcfs and not @dirs ) );
pod2usage(
    -message =>
      "Syntax error - cannot use --reject and --not_samples argument together",
    exitval => 2
) if ( @reject and @ignore_samples );
pod2usage(
    -message =>
"Syntax error - cannot use --reject or --not_samples arguments with --info_filter option",
    exitval => 2
  )
  if ( @reject or @ignore_samples )
  and $filter_with_info;
pod2usage(
    -message => "Variant quality scores must be 0 or greater.\n",
    -exitval => 2
) if ( $min_qual < 0 );
pod2usage(
    -message => "Genotype quality scores must be 0 or greater.\n",
    -exitval => 2
) if ( $minGQ < 0 );
pod2usage(
    -message =>
      "--allele_frequency_filter (-y) value must be between 0 and 1.\n",
    -exitval => 2
) if ( $maf < 0 or $maf > 1 );

if ( defined $aff_quality ) {
    pod2usage(
        -message => "Genotype quality scores must be 0 or greater.\n",
        -exitval => 2
    ) if $aff_quality < 0;
}
else {
    $aff_quality = $minGQ;
}
if ( defined $unaff_quality ) {
    pod2usage(
        -message => "Genotype quality scores must be 0 or greater.\n",
        -exitval => 2
    ) if $unaff_quality < 0;
}
else {
    $unaff_quality = $minGQ;
}
if (not @pop_ids){
    @pop_ids =  qw / AFR AMR EAS FIN NFE SAS / ; #default is ExAC pop codes
} 
elsif ( grep{ /^disable$/i } @pop_ids ){
    @pop_ids = ();
    print STDERR "[INFO] Filtering on population allele count / allele number disabled.\n";
}

my $cpus = Sys::CPU::cpu_count();
if ( $forks < 2 ) {
    $forks = 0;    #no point having overhead of forks for one fork
}
else {
    if ( $forks > $cpus ) {
        print STDERR
"[Warning]: Number of forks ($forks) exceeds number of CPUs on this machine ($cpus)\n";
    }
    if ( not $buffer_size ) {
        $buffer_size = 10000 > $forks * 1000 ? 10000 : $forks * 1000;
    }
    print STDERR
"[INFO] Processing in batches of $buffer_size variants split among $forks forks.\n";
}
push @filter_vcfs, get_vcfs_from_directories( \@dirs, $regex ) if @dirs;
die "No VCFs found to use as filters.\n" if not @filter_vcfs;

my $time = strftime( "%H:%M:%S", localtime );
my $total_variants;
print STDERR "[$time] Initializing input VCF... ";
die "Header not ok for input ($vcf) "
  if not VcfReader::checkHeader( vcf => $vcf );
if ( defined $progress ) {
    if ( $vcf eq "-" ) {
        print STDERR "Can't use --progress option when input is from STDIN\n";
        $progress = 0;
    }
    else {
        $total_variants = VcfReader::countVariants($vcf);
        print STDERR "$vcf has $total_variants variants. ";
    }
}
my %contigs = ();
if ($forks){
    %contigs       = VcfReader::getContigOrder( $vcf );
}
my %sample_to_col = ();
if (@samples) {
    %sample_to_col = VcfReader::getSamples(
        vcf         => $vcf,
        get_columns => 1,
    );
}

$time = strftime( "%H:%M:%S", localtime );
print STDERR "\n[$time] Finished initializing input VCF\n";

my %filter_vcf_to_index = ();
my %filter_vcf_samples  = ();
my %filter_vcf_info     = ();

my $dbpm = Parallel::ForkManager->new($forks);
for ( my $i = 0 ; $i < @filter_vcfs ; $i++ ) {
    $dbpm->run_on_finish(    # called BEFORE the first call to start()
        sub {
            my ( $pid, $exit_code, $ident, $exit_signal, $core_dump,
                $data_structure_reference )
              = @_;

            if ( defined($data_structure_reference) ) {
                if ( ref $data_structure_reference eq 'ARRAY' ) {
                    $filter_vcf_to_index{ $data_structure_reference->[3] } =
                      $data_structure_reference->[0];
                    $filter_vcf_samples{ $data_structure_reference->[3] } =
                      $data_structure_reference->[1];
                    $filter_vcf_info{ $data_structure_reference->[3] } =
                      $data_structure_reference->[2];
                }
                else {
                    die
                      "Unexpected return from dbSNP reference initialization:\n"
                      . Dumper $data_structure_reference;
                }
            }
            else {    # problems occuring during storage or retrieval
                die "No message received from child process $pid!\n";
            }
        }
    );
    $time = strftime( "%H:%M:%S", localtime );
    print STDERR "[$time] Initializing $filter_vcfs[$i] filter reference VCF "
      . ( $i + 1 ) . " of "
      . scalar(@filter_vcfs) . "\n";
    $dbpm->start() and next;
    my @index_samp_and_info = initializeFilterVcfs( $filter_vcfs[$i] );
    push @index_samp_and_info, $filter_vcfs[$i];
    $time = strftime( "%H:%M:%S", localtime );
    print STDERR "[$time] Finished initializing $filter_vcfs[$i] filter VCF.\n";
    $dbpm->finish( 0, \@index_samp_and_info );
}
$dbpm->wait_all_children;
if ($filter_with_info) {
    check_filter_vcf_info_fields();
}

#only retain samples from filter vcfs if not specified by --not_samples argument
#if --reject argument is used only keep samples specified by --reject
#if neither argument is used keep all
foreach my $f ( keys %filter_vcf_samples ) {
    if ( keys %{ $filter_vcf_samples{$f} } == 0 ) {
        if (@reject) {
            die
"Filter VCF $f has no samples - cannot run with --reject argument.\n"
              if keys %{ $filter_vcf_samples{$f} } == 0;
        }
        if (@ignore_samples) {
            die
"Filter VCF $f has no samples - cannot run with --not_samples argument.\n";
        }
        if ( $threshold or $maf or $filter_homozygotes ) {
            if ( not $filter_with_info ) {
                die
"Filter VCF $f has no samples. Use of --allele_frequency_filter, --threshold".
" or --filter_homozygotes options is not allowed with filter VCFs without samples ".
"unless used with the --info_filter option.\n";
            }
        }
    }

    foreach my $ignore (@ignore_samples) {
        if ( exists $filter_vcf_samples{$f}->{$ignore} ) {
            delete $filter_vcf_samples{$f}->{$ignore};
        }
    }
    if (@reject) {
        foreach my $samp ( keys %{ $filter_vcf_samples{$f} } ) {
            if ( not grep { $_ eq $samp } @reject ) {
                delete $filter_vcf_samples{$f}->{$samp};
            }
        }
    }
}
my $OUT;
if ($out) {
    open( $OUT, ">$out" ) or die "Can't open $out for writing: $!\n";
}
else {
    $OUT = \*STDOUT;
}

my $prev_chrom = 0;
my $progressbar;
my $next_update = 0;
if ($progress) {
    $progressbar = Term::ProgressBar->new(
        {
            name  => "Filtering",
            count => $total_variants * 3,

            #ETA   => "linear",
        }
    );
}else{
    $time = strftime( "%H:%M:%S", localtime );
    print STDERR "[$time] Filtering started.\n";
}

my $meta_head = VcfReader::getMetaHeader($vcf);
print $OUT "$meta_head\n";
print $OUT "##filterVcfOnVcf.pl=\"";
my @opt_string = ();
foreach my $k ( sort keys %opts ) {
    if ( not ref $opts{$k} ) {
        push @opt_string, "$k=$opts{$k}";
    }
    elsif ( ref $opts{$k} eq 'SCALAR' ) {
        if ( defined ${ $opts{$k} } ) {
            push @opt_string, "$k=${$opts{$k}}";
        }
        else {
            push @opt_string, "$k=undef";
        }
    }
    elsif ( ref $opts{$k} eq 'ARRAY' ) {
        if ( @{ $opts{$k} } ) {
            push @opt_string, "$k=" . join( ",", @{ $opts{$k} } );
        }
        else {
            push @opt_string, "$k=undef";
        }
    }
}
my $col_header = VcfReader::getColumnHeader($vcf);
print $OUT join( " ", @opt_string ) . "\"\n" . $col_header . "\n";
my $kept             = 0;
my $filtered         = 0;
my $n                = 0;
my $variants_done    = 0;
my @lines_to_process = ();
my $VCF              = VcfReader::_openFileHandle($vcf);
my %no_fork_args     = ();

if ( $forks < 2 ) {
    foreach my $f (@filter_vcfs) {
        my %s = VcfReader::getSearchArguments( $f, $filter_vcf_to_index{$f} );
        $no_fork_args{$f} = \%s;
    }
}
LINE: while ( my $line = <$VCF> ) {
    next if $line =~ /^#/;
    $variants_done++;
    $n++;
    if ($progress) {
        $next_update = $progressbar->update($n) if $n >= $next_update;
    }
    if ( $forks > 1 ) {
        push @lines_to_process, $line;
        if ($progressbar) {
            $next_update = $progressbar->update($n) if $n >= $next_update;
        }
        if ( @lines_to_process >= $buffer_size ) {
            process_buffer();
            @lines_to_process = ();
        }
    }
    else {
        chomp $line;
        my @split = split( "\t", $line );
        my $l = filter_on_vcf( \@split, \%no_fork_args );
        if ($l) {
            print $OUT "$line\n";
            $kept++;
        }
        else {
            $filtered++;
        }
        $n += 2;
        if ($progressbar) {
            $next_update = $progressbar->update($n) if $n >= $next_update;
        }
    }
}
process_buffer() if $forks > 1;
if ($progressbar) {
    $progressbar->update( $total_variants * 3 )
      if $total_variants * 3 >= $next_update;
}
close $VCF;
close $OUT;
$time = strftime( "%H:%M:%S", localtime );
print STDERR "\n[$time] $filtered variants filtered, $kept printed ";
print STDERR "($total_variants total)" if $total_variants;
print STDERR "\n";

################################################
#####################SUBS#######################
################################################
sub process_buffer {
    return if not @lines_to_process;
    my @lines_to_print;
    my $lines_per_slice = @lines_to_process;
    if ( $forks > 0 ) {
        $lines_per_slice =
            int( @lines_to_process / $forks ) > 1
          ? int( @lines_to_process / $forks )
          : 1;
    }
    my @batch = ();

    #get a batch for each thread
    for ( my $i = 0 ; $i < @lines_to_process ; $i += $lines_per_slice ) {
        my $last =
          ( $i + $lines_per_slice - 1 ) < $#lines_to_process
          ? $i + $lines_per_slice - 1
          : $#lines_to_process;
        if ( $i + $lines_per_slice >= @lines_to_process ) {
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

            if ( defined($data_structure_reference) ) {
                my %res = %{$data_structure_reference};
                if ( ref $res{keep} eq 'ARRAY' ) {
                    push @lines_to_print, \@{ $res{keep} } if @{ $res{keep} };
                }
                if ( $res{filter} ) {
                    $filtered += $res{filter};
                }
                if ($progressbar) {
                    $n += $res{batch_size};
                    $next_update = $progressbar->update($n)
                      if $n >= $next_update;
                }
            }
            else {
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
    @lines_to_print =
      sort { VcfReader::by_first_last_line( $a, $b, \%contigs ) }
      @lines_to_print;
    if (@lines_to_print) {
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
                    $next_update = $progressbar->update($n)
                      if $n >= $next_update;
                }
            }
        }
    }
    else {
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
    my %sargs;
    foreach my $f (@filter_vcfs) {

        #WE'VE FORKED AND COULD HAVE A RACE CONDITION HERE -
        # WHICH IS WHY WE DO A FIRST PASS WITH OUR initializeDbsnpVcfs
        # METHOD TOWARDS THE START OF THE PROGRAM
        my %s = VcfReader::getSearchArguments( $f, $filter_vcf_to_index{$f} );
        $sargs{$f} = \%s;
    }
    my %results = ( batch_size => scalar(@$batch) );
    foreach my $line ( @{$batch} ) {
        chomp $line;
        my @split = split( "\t", $line );
        my $l = filter_on_vcf( \@split, \%sargs );
        if ($l) {
            push @{ $results{keep} }, $l;
        }
        else {
            $results{filter}++;
        }
    }
    return %results;
}

################################################
sub filter_on_info_fields {
    my ( $vcf_line, $search_args ) = @_;
    my $qual  = VcfReader::getVariantField( $vcf_line, 'QUAL' );
    my $chrom = VcfReader::getVariantField( $vcf_line, 'CHROM' );
    if ($min_qual &&  $qual < $min_qual ) {
        return;
    }

#process each allele separately in case we have MNVs/deletions that need to be simplified

    my %min_vars       = VcfReader::minimizeAlleles( $vcf_line, );
    my %sample_alleles = ();
    if (@samples) {
        %sample_alleles = map { $_ => undef } VcfReader::getSampleCall(
            line                => $vcf_line,
            multiple            => \@samples,
            return_alleles_only => 1,
            minGQ               => $aff_quality,
            sample_to_columns   => \%sample_to_col
        );
        delete $sample_alleles{0};
        delete $sample_alleles{'.'};
        return
          if not keys
          %sample_alleles;  #filter if we don't have any variants in our samples
    }
    else {
        #if no samples specified use all alleles for %sample_alleles
        %sample_alleles = map { $_ => undef } keys %min_vars;
    }

    my %sample_matches = ()
      ; #check each allele matches in all filter_vcfs but don't reset after each file
    my %thresh_counts =
      ();    #count samples in all filter_vcfs i.e. don't reset after each file
    my %af_counts
      ; #count allele occurences and total alleles to calculate allele frequency and don't reset after each file
    my %f_genos = (); #store genotypes as keys if we're using $filter_homozygote
    my %alleles_over_maf = ();    #store alleles that have exceeded maf in here

  FILTER: foreach my $f (@filter_vcfs) {
      ALLELE: foreach my $allele ( keys %sample_alleles ) {
            if (
                my @snp_hits = VcfReader::searchForPosition(
                    %{ $search_args->{$f} },
                    chrom => $min_vars{$allele}->{CHROM},
                    pos   => $min_vars{$allele}->{POS}
                )
              )
            {

                #get genotype call codes (0, 1, 2 etc.) for filter samples
              FILTER_LINE: foreach my $snp_line (@snp_hits) {
                    my @f_alts      = ();
                    my %geno_counts = ();
                    my @snp_split   = split( "\t", $snp_line );
                    if ($min_qual) {
                        my $filter_qual =
                          VcfReader::getVariantField( \@snp_split, 'QUAL' );
                        next FILTER_LINE if $filter_qual < $min_qual;
                    }
                    my %filter_min = VcfReader::minimizeAlleles( \@snp_split, );

#%f_alts = map {$_ => undef} $filter_obj->getSampleActualGenotypes(multiple => \@temp_reject, return_alleles_only => 1, minGQ => $unaff_quality);
                    my @alts = VcfReader::readAlleles( line => \@snp_split, );
                    for ( my $i = 0 ; $i < @alts ; $i++ ) {
                        push @f_alts, $i;
                    }
                    if ( $filter_homozygotes or ( $maf and @filter_vcfs > 1 ) )
                    {
                        my @gtcs = split(
                            ",",
                            VcfReader::getVariantInfoField( \@snp_split,
                                "GTC", )
                        );
                        my @pgts = split(
                            ",",
                            VcfReader::getVariantInfoField(
                                \@snp_split, "PGTS",
                            )
                        );
                        die
"Genotype count (GTC) does not have the same number of alleles as possible genotypes (PGTS) field for line:\n$snp_line\n"
                          if @gtcs != @pgts;
                        for ( my $i = 0 ; $i < @gtcs ; $i++ ) {
                            $geno_counts{ $pgts[$i] } = $gtcs[$i];
                        }
                    }
                    my $filter_match = ''
                      ; #if one of the filter's ALTs matches store the ALT allele code here
                  ALT: foreach my $alt (@f_alts) {
                        next ALT if $alt eq '.';
                        next ALT if $alt == 0;
                        next ALT
                          if $min_vars{$allele}->{POS} ne
                          $filter_min{$alt}->{POS};
                        next ALT
                          if $min_vars{$allele}->{REF} ne
                          $filter_min{$alt}->{REF};
                        next ALT
                          if $min_vars{$allele}->{ALT} ne
                          $filter_min{$alt}->{ALT};
                        $min_vars{$allele}->{CHROM} =~ s/^chr//;
                        $filter_min{$alt}->{CHROM} =~ s/^chr//;
                        next ALT
                          if $min_vars{$allele}->{CHROM} ne
                          $filter_min{$alt}->{CHROM};
                        #filter ALT matches input allele
                        if (    $filter_homozygotes
                            and not $threshold
                            and not $maf )
                        {
           #if using $filter_homozygotes on its own we only consider something a
           #'match' if it's homozygous
                            if (   exists $geno_counts{"$alt/$alt"}
                                or exists $geno_counts{"$alt|$alt"} )
                            {
                                $filter_match = $alt;
                                $sample_matches{$allele}++;
                                last ALT;
                            }
                        }
                        else {
                            $filter_match = $alt;
                            $sample_matches{$allele}++;
                            last ALT;
                        }
                    }
                    if ( not $filter_match ) {
                        next FILTER_LINE;
                    }

                    if ($threshold) {
                        foreach my $k ( keys %geno_counts ) {
                            my @g_alleles = split( /[\/\|]/, $k );
                            if ( grep { $_ eq $filter_match } @g_alleles ) {
                                $thresh_counts{$allele} += $geno_counts{$k};
                            }
                        }
                    }

                    if ($maf) {
                        if (    @filter_vcfs == 1
                            and exists $filter_vcf_info{$f}->{AF})
                        {
                            my @afs = split(
                                ",",
                                VcfReader::getVariantInfoField(
                                    \@snp_split, "AF",
                                )
                            );
                            my $freq = $afs[ $filter_match - 1 ];
                            if ( $freq >= $maf ) {
                                $alleles_over_maf{$allele}++;
                            }elsif(@pop_acs){
                                if ( pop_freq_over_maf(\@snp_split, $filter_match) ) {
                                    $alleles_over_maf{$allele}++;
                                }
                            }
                        }
                        else {
                            foreach my $k ( keys %geno_counts ) {
                                my @g_alleles = split( /[\/\|]/, $k );
                                foreach my $g (@g_alleles) {
                                    if ( $g eq $filter_match ) {
                                        $af_counts{$allele}->{counts} +=
                                          $geno_counts{$k};
                                    }
                                    $af_counts{$allele}->{total} +=
                                      $geno_counts{$k};
                                }
                            }
                            if ( @filter_vcfs == 1 ) {
                                my $freq = 0;
                                eval {
                                    $freq =
                                      $af_counts{$allele}->{counts} /
                                      $af_counts{$allele}->{total};
                                };
                                $alleles_over_maf{$allele}++ if $freq >= $maf;
                            }
                        }
                    }
                }    #read pos
            }    #search
        }   #foreach allele
            #done each allele - see if we've got enough data to filter this line
            #and can skip other filter vcfs for win
        next FILTER
          if keys %sample_alleles !=
          keys %sample_matches;    #can't skip - not all alleles accounted for
        my $homozygous_alleles = 0;
        if ($filter_homozygotes) {
            foreach my $allele ( keys %sample_alleles ) {
                if (   exists $f_genos{"$allele/$allele"}
                    or exists $f_genos{"$allele|$allele"} )
                {
                    $homozygous_alleles++;
                }
            }
            if ( $homozygous_alleles == keys %sample_alleles ) {
                if ($print_matching) {
                    return $vcf_line;
                }
                else {
                    return;
                }
            }
        }

        foreach my $allele ( keys %sample_alleles ) {
            if ($threshold) {
                next FILTER if $thresh_counts{$allele} < $threshold;
            }
        }    #if we haven't gone to next FILTER and
             #are not filtering on allele frequency we can filter this line
             #no need to look at other filter_vcfs
        if ( not $maf ) {
            if ($print_matching) {
                return $vcf_line;
            }
            else {
                return;
            }
        }
    }    #foreach filter vcf
         #no line meeting criteria for allele in any filter vcf
         #now check allele frequency if given
    if ($maf) {
        if ( @filter_vcfs != 1 ) {
          COUNTS: foreach my $allele ( keys %sample_alleles ) {
                last COUNTS if not exists $af_counts{$allele};
                if (
                    (
                        $af_counts{$allele}->{counts} /
                        $af_counts{$allele}->{total}
                    ) >= $maf
                  )
                {
                    $alleles_over_maf{$allele}++;
                }
                else {
                    last COUNTS;
                }
            }
        }
        if ( keys %alleles_over_maf == keys %sample_alleles ) {
            if ($print_matching) {
                return $vcf_line;
            }
            else {
                return;
            }
        }
    }
    if ($print_matching) {
        return;
    }
    else {
        return $vcf_line;
    }
}
#################################################
sub pop_freq_over_maf{
    my ($split, $allele_code) = @_;
    # calculate AF for each @pops pop 
    # and return 1 if > $maf
    foreach my $pop (@pop_acs){
        my $an = VcfReader::getVariantInfoField($split, "AN_$pop");
        next if not $an;
        next if $an < $MIN_AN;
        my $acs = VcfReader::getVariantInfoField($split, "AC_$pop");
        my $count = (split ",", $acs)[$allele_code -1];
        if ( $count / $an > $maf){
            return 1;
        }
    }
    return 0;
}
################################################
sub filter_on_vcf {
    my ( $vcf_line, $search_args ) = @_;
    if ($filter_with_info) {
        return filter_on_info_fields( $vcf_line, $search_args );
    }
    else {
        return filter_on_vcf_samples( $vcf_line, $search_args );
    }
}
################################################
sub filter_on_vcf_samples {
    my ( $vcf_line, $search_args ) = @_;
    my $qual  = VcfReader::getVariantField( $vcf_line, 'QUAL' );
    my $chrom = VcfReader::getVariantField( $vcf_line, 'CHROM' );
    if ($min_qual &&  $qual < $min_qual ) {
        return;
    }

#process each allele separately in case we have MNVs/deletions that need to be simplified

    my %min_vars       = VcfReader::minimizeAlleles( $vcf_line, );
    my %sample_alleles = ();
    if (@samples) {
        %sample_alleles = map { $_ => undef } VcfReader::getSampleCall(
            line                => $vcf_line,
            multiple            => \@samples,
            return_alleles_only => 1,
            minGQ               => $aff_quality,
            sample_to_columns   => \%sample_to_col
        );
        delete $sample_alleles{0};
        delete $sample_alleles{'.'};
        return
          if not keys
          %sample_alleles;  #filter if we don't have any variants in our samples
    }
    else {
        #if no samples specified use all alleles for %sample_alleles
        %sample_alleles = map { $_ => undef } keys %min_vars;
    }
    my %sample_matches = ()
      ; #check each allele matches in all filter_vcfs but don't reset after each file
    my %thresh_counts =
      ();    #count samples in all filter_vcfs i.e. don't reset after each file
    my %af_counts
      ; #count allele occurences and total alleles to calculate allele frequency
    my %f_genos =
      ();    #store genotypes as keys if we're using $filter_homozygotes

  FILTER: foreach my $f (@filter_vcfs) {
        my @temp_reject = keys %{ $filter_vcf_samples{$f} };
      ALLELE: foreach my $allele ( keys %sample_alleles ) {
            if (
                my @snp_hits = VcfReader::searchForPosition(
                    %{ $search_args->{$f} },
                    chrom => $min_vars{$allele}->{CHROM},
                    pos   => $min_vars{$allele}->{POS}
                )
              )
            {

                #get genotype call codes (0, 1, 2 etc.) for filter samples
              FILTER_LINE: foreach my $snp_line (@snp_hits) {
                    my @snp_split = split( "\t", $snp_line );
                    if ($min_qual) {
                        my $filter_qual =
                          VcfReader::getVariantField( \@snp_split, 'QUAL' );
                        next FILTER_LINE if $filter_qual < $min_qual;
                    }
                    my %filter_min = VcfReader::minimizeAlleles( \@snp_split, );

#%f_alts = map {$_ => undef} $filter_obj->getSampleActualGenotypes(multiple => \@temp_reject, return_alleles_only => 1, minGQ => $unaff_quality);
                    my @f_alts = VcfReader::getSampleCall(
                        line                => \@snp_split,
                        multiple            => \@temp_reject,
                        sample_to_columns   => $filter_vcf_samples{$f},
                        return_alleles_only => 1,
                        minGQ               => $unaff_quality
                    );
                    if ($filter_homozygotes) {
                        my %genos = VcfReader::getSampleCall(
                            line              => \@snp_split,
                            multiple          => \@temp_reject,
                            minGQ             => $unaff_quality,
                            sample_to_columns => $filter_vcf_samples{$f},
                        );
                        %genos = reverse %genos;
                        foreach my $k ( keys %genos ) {
                            $f_genos{$k}++;
                        }
                    }
                    my $filter_match = ''
                      ; #if one of the filter's ALTs matches store the ALT allele code here
                  ALT: foreach my $alt (@f_alts) {
                        next ALT if $alt eq '.';
                        next ALT if $alt == 0;
                        next ALT
                          if $min_vars{$allele}->{POS} ne
                          $filter_min{$alt}->{POS};
                        next ALT
                          if $min_vars{$allele}->{REF} ne
                          $filter_min{$alt}->{REF};
                        next ALT
                          if $min_vars{$allele}->{ALT} ne
                          $filter_min{$alt}->{ALT};
                        $min_vars{$allele}->{CHROM} =~ s/^chr//;
                        $filter_min{$alt}->{CHROM} =~ s/^chr//;
                        next ALT
                          if $min_vars{$allele}->{CHROM} ne
                          $filter_min{$alt}->{CHROM};

                        if (    $filter_homozygotes
                            and not $threshold
                            and not $maf )
                        {
           #is using $filter_homozygotes on its own we only consider something a
           #'match' if it's homozygous
                            if (   exists $f_genos{"$alt/$alt"}
                                or exists $f_genos{"$alt|$alt"} )
                            {
                                $filter_match = $alt;
                                $sample_matches{$allele}++;
                                last ALT;
                            }
                        }
                        else {
                            $filter_match = $alt;
                            $sample_matches{$allele}++;
                            last ALT;
                        }
                    }
                    if ( not $filter_match ) {
                        next FILTER_LINE;
                    }

                    if ($threshold) {
                        foreach my $t_samp (@temp_reject) {
                            my %t_alleles =
                              map { $_ => undef } VcfReader::getSampleCall(
                                line                => \@snp_split,
                                sample              => $t_samp,
                                return_alleles_only => 1,
                                minGQ               => $unaff_quality,
                                sample_to_columns   => $filter_vcf_samples{$f},
                              );
                            if ( exists $t_alleles{$filter_match} ) {
                                $thresh_counts{$allele}++;
                            }
                        }

                       #foreach my $alt (@alts){
                       #   next FILTER_LINE if not exists $thresh_counts{$alt};
                       #  next FILTER_LINE if $thresh_counts{$alt} < $threshold;
                       #}
                    }
                    if ($maf) {
                        my %f_allele_counts = VcfReader::countAlleles(
                            line  => \@snp_split,
                            minGQ => $unaff_quality
                        );
                        foreach my $f_al ( keys %f_allele_counts ) {
                            if ( $f_al eq $filter_match ) {
                                $af_counts{$allele}->{counts} +=
                                  $f_allele_counts{$f_al};
                            }
                            $af_counts{$allele}->{total} +=
                              $f_allele_counts{$f_al};
                        }
                    }
                }    #read pos
            }    #search
        }   #foreach allele
            #done each allele - see if we've got enough data to filter this line
            #and can skip other filter vcfs for win
        next FILTER
          if keys %sample_alleles !=
          keys %sample_matches;    #can't skip - not all alleles accounted for
        my $homozygous_alleles = 0;
        if ($filter_homozygotes) {
            foreach my $allele ( keys %sample_alleles ) {
                if (   exists $f_genos{"$allele/$allele"}
                    or exists $f_genos{"$allele|$allele"} )
                {
                    $homozygous_alleles++;
                }
            }
            if ( $homozygous_alleles == keys %sample_alleles ) {
                if ($print_matching) {
                    return $vcf_line;
                }
                else {
                    return;
                }
            }
        }

        foreach my $allele ( keys %sample_alleles ) {
            if ($threshold) {
                next FILTER if $thresh_counts{$allele} < $threshold;
            }
        }    #if we haven't gone to next FILTER and
             #are not filtering on allele frequency we can filter this line
             #no need to look at other filter_vcfs
        if ( not $maf ) {
            if ($print_matching) {
                return $vcf_line;
            }
            else {
                return;
            }
        }
    }    #foreach filter vcf
         #no line meeting criteria for allele in any filter vcf
         #now check allele frequency if given
    if ($maf) {
        my $alleles_over_maf = 0;
      COUNTS: foreach my $allele ( keys %sample_alleles ) {
            last COUNTS if not exists $af_counts{$allele};
            if (
                (
                    $af_counts{$allele}->{counts} / $af_counts{$allele}->{total}
                ) >= $maf
              )
            {
                $alleles_over_maf++;
            }
            else {
                last COUNTS;
            }
        }
        if ( $alleles_over_maf == keys %sample_alleles ) {
            if ($print_matching) {
                return $vcf_line;
            }
            else {
                return;
            }
        }
    }
    if ($print_matching) {
        return;
    }
    else {
        return $vcf_line;
    }
}
#################################################
sub initializeFilterVcfs {
    my ($snpfile) = @_;

    #we getSearchArguments here simply to prevent a race condition later
    my @head = VcfReader::getHeader($snpfile);
    die "Header not ok for $snpfile "
      if not VcfReader::checkHeader( header => \@head );

    #my %sargs = VcfReader::getSearchArguments($snpfile);
    my %index = VcfReader::readIndex($snpfile);
    my %samp  = VcfReader::getSamples(
        vcf         => $snpfile,
        get_columns => 1
    );
    my %info_fields = VcfReader::getInfoFields( header => \@head );
    return ( \%index, \%samp, \%info_fields );
}

#################################################
sub check_filter_vcf_info_fields {
    foreach my $f (@filter_vcfs) {
        if (   $threshold
            or $filter_homozygotes
            or ( $maf and @filter_vcfs > 1 ))
        {
            check_pgts_gtc($f);
        }
        if ( $maf and @filter_vcfs == 1 ){
            if ( not exists $filter_vcf_info{$f}->{AF} ) {
                check_pgts_gtc($f);
            }
            @pop_acs = check_pop_ac_info($f);
            if (@pop_acs){
                print STDERR "[INFO] Found " . scalar(@pop_acs) . " population " . 
                    "allele counts in VCF INFO fields. These will be used for per ".
                    "population allele-frequency filtering. To disable use " .
                    "'--population_ids disable' when running this program.\n";
            }
        }
    }
}

#################################################
sub check_pop_ac_info{
    my $f = shift;
    my @pop_codes = ();
    foreach my $pop ( @pop_ids ){
        if (exists $filter_vcf_info{$f}->{"AC_$pop"} and 
            exists $filter_vcf_info{$f}->{"AN_$pop"}){
            next if ($filter_vcf_info{$f}->{"AC_$pop"}->{Number} ne 'A');
            next if ($filter_vcf_info{$f}->{"AC_$pop"}->{Type} ne 'Integer');
            next if ($filter_vcf_info{$f}->{"AN_$pop"}->{Number} ne '1');
            next if ($filter_vcf_info{$f}->{"AN_$pop"}->{Type} ne 'Integer');
            push @pop_codes, $pop;
        }
    }
    return @pop_codes;
}
#################################################
sub check_pgts_gtc {
    my $f = shift;
    if ( not exists $filter_vcf_info{$f}->{PGTS} ) {
        die
"Cannot filter on genotype information with the --info_filter option without " . 
"possible genotypes (PGTS) info field. Filter VCF $f is missing this field from " .
"its header. This INFO field is available by processing a VCF containing multiple" . 
" samples using the sampleCallsToInfo.pl script.\n";
    }
    elsif ( $filter_vcf_info{$f}->{PGTS}->{Description} ne
        "\"Possible Genotype Call Codes (for ease of reference).\"" )
    {
        die
"Filter VCF $f possible genotypes (PGTS) field does not match expected INFO field " .
"description. This INFO field is available by processing a VCF containing multiple "
."samples using the sampleCallsToInfo.pl script.\n";
    }

    if ( not exists $filter_vcf_info{$f}->{GTC} ) {
        die
"Cannot filter on genotype information with the --info_filter option without ".
"genotype counts (GTC) info field. Filter VCF $f is missing this field from ".
"its header. This INFO field is available by processing a VCF containing multiple "
."samples using the sampleCallsToInfo.pl script.\n";
    }
    elsif ( $filter_vcf_info{$f}->{GTC}->{Description} ne
        "\"Genotype counts in order of genotypes listed by PGTS field.\"" )
    {
        die
"Filter VCF $f genotype counts (GTC) field does not match expected INFO field ".
"description. This INFO field is available by processing a VCF containing multiple ".
"samples using the sampleCallsToInfo.pl script.\n";
    }
}

#################################################
sub get_vcfs_from_directories {
    my ( $dirs, $regex ) = @_;
    my @vcfs = ();
    foreach my $d (@$dirs) {
        $d =~ s/\/$//;
        opendir( my $DIR, $d ) or die "Can't open directory $d: $!\n";

        #my @dir_vcfs = grep {/\.vcf(\.gz)*$/i} readdir $DIR;
        my @dir_vcfs =
          grep { /\.vcf$/i }
          readdir $DIR;   #can't use gzipped vcfs with ParseVCF search functions
        if ( not @dir_vcfs ) {
            print STDERR "WARNING - no VCF files in directory $d\n";
        }
        elsif ($regex) {
            @dir_vcfs = grep { /$regex/ } @dir_vcfs;
            if ( not @dir_vcfs ) {
                print STDERR
"WARNING - no VCF files matching regex /$regex/ in directory $d\n"
                  if not @dir_vcfs;
            }
            else {
                foreach my $v (@dir_vcfs) {
                    push @vcfs, "$d/$v";
                }
            }
        }
        else {
            foreach my $v (@dir_vcfs) {
                push @vcfs, "$d/$v";
            }
        }
    }
    return @vcfs;
}
