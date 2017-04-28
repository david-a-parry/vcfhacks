#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use Parallel::ForkManager;
use Sys::CPU;
use Pod::Usage;
use Data::Dumper;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use POSIX qw(strftime);
use FindBin qw($RealBin);
use lib "$RealBin/lib";
use lib "$RealBin/lib/dapPerlGenomicLib";
use VcfReader 0.3;
use VcfhacksUtils;

=head1 NAME

filterVcfOnVcf.pl - filter a VCF file using variants from one or more other VCF files.

=head1 SYNOPSIS

filterVcfOnVcf.pl -i <input vcf file> -f <vcf to use to filter input> [options]

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

A VCF file to use to filter input. Variants from the input VCF that match a variant in this file will be filtered from the output. By default, if any sample matches a variant it will be filtered but combinations of --reject, --not_samples, --allele_frequency_filter, --threshold and --filter_homozygotes can be used to modify this behaviour. However, if more than one allele is present for a variant in your input VCF it will only be filtered if ALL alleles are matched in your filter VCF. 

Alternatively, if you want to filter ANY matching variant regardless of genotypes in your --filter VCF, simply add the -w/--info_filter argument without specifying any values for --allele_frequency_filter, --threshold or --filter_homozygotes.

=item B<-s    --samples>

Samples from input file to check. If specified only variant alleles from these samples will be used to compare against variants in filter VCFs. Any line in your input without a variant call in these samples will be filtered. Default is to look at all alleles for each variant.

=item B<-r    --reject>

Samples from filter VCF files to use for filtering.  If specified only variant alleles from these samples will be used to compare with the input VCF. Default is to look at all alleles.

=item B<-x    --not_samples>

Samples from filter VCF files to ignore.  If specified variant alleles from all samples except these will be used to compare with the input VCF. Default is to look at all alleles.

=item B<-w    --info_filter>

Use this flag to use metrics written to the INFO field of your filter VCF for filtering on frequency, homozygosity or threshold counts rather than checking sample genotypes. Filtering on allele frequency (see --allele_frequency_filter option below) can be performed as long as your VCF contains allele frequency (AF) INFO tags and optionally population specific allele count (AC) and allele number (AN) tags (see --population_ids option below). This option is particularly useful for filtering against VCFs from ExAC using global and population allele frequencies. Alternatively, you can pre-process your own filter VCFs with the sampleCallsToInfo.pl script to reduce genotype information to AF, PGTS and GTC fields. Converting filter VCFs containing many hundreds/thousands of samples  with sampleCallsToInfo.pl for use with this option will speed up your analyses by many orders of magnitude.

This option can also be used without using --allele_frequency_filter/--threshold/--filter_homozygotes arguments in order to filter any matching variant regardless of sample genotypes. 

=item B<-y    --allele_frequency_filter>

Reject variants if the allele frequency in your --filter VCF is equal to or greater than this value. If samples from the filter VCF have been specified using the --reject or --not_samples argument, only the relevant samples will be used for calculating allele frequency. Variants will only be counted if the genotype quality is greater than the value given for --un_quality. The value must be a float between 0 and 1. 

=item B<--population_ids>

If using the --info_filter option and --allele_frequency_filter option, by default the filter VCF will be scanned for allele counts (e.g. AFR_AC) and allele number (e.g. AFR_AN) INFO tags for the population codes as used by EXAC (AFR, AMR, EAS, FIN, NFE, SAS). An alternative list of population IDs to search can be entered here. To disable this feature use this option specifying 'disable' as an argument.
    
=item B<-n    --minimum_allele_number>

Minimum number of allele counts for a variant found in a --filter VCF before it will be used for filtering on allele frequency. If filtering using population allele counts (see --population_ids above) this sets the minimum number of alleles that have to be called in a population (default = 100) before the frequency value will be used for filtering. If filtering using samples in a VCF the default is 0 (i.e. no minimum allele number).

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

=item B<--annotation>

Use this option to specify an annotation for allele frequencies and allele numbers calculated from samples/INFO fields in your --filter VCF. These annotations will be added to the INFO field of the output as 'annotation_AF' and 'annotation_AN' where 'annotation' is the value you supply to this option. For example, if specifying '--annotation controls' it would add 'controls_AF' and 'controls_AN' INFO fields giving the allele frequency and allele numbers respectively for the samples used in your --filter VCF. These annotations can be used by programs such as findBiallelic.pl and getFunctionalVariants.pl for filtering using custom allele frequencies.

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

Filter variants from a VCF using another VCFs.  You may specify samples or a number of samples with matching variants to filter with.  You can also use minimum genotype or variant qualities to filter with. Alternatively, filtering can be performed on certain annotations in the INFO field of each variant using the --info_filter option.

By default, each variant from the given --input file will be assessed and if all alternative/variant alleles for a given variant are represented in the --filter VCF provided, the variant will be filtered. You may use --reject or --not_samples arguments to only filter variant alleles if present in specific samples in the --filter VCF. You may also use the --samples argument to only compare variants from specific samples in your --input VCF and filter any variants that do not have alleles represented by these samples. See above for details of all available options. 

=head1 EXAMPLES

    filterVcfOnVcf.pl -i input.vcf -f controls.vcf -o filtered.vcf
    (filter variants in input.vcf if present in controls.vcf - at least one sample in controls.vcf must carry the variant for it to be filtered)
    
    filterVcfOnVcf.pl -i input.vcf -f controls.vcf -o filtered.vcf -y 0.01
    (filter variants in input.vcf if present at an allele frequency of 1% or higher in controls.vcf)

    filterVcfOnVcf.pl -i input.vcf -f controls.vcf -o filtered.vcf -w
    (filter variants in input.vcf if present in controls.vcf - genotypes of samples in controls.vcf will not be checked)
    
    filterVcfOnVcf.pl -i input.vcf -f ExAC.r0.3.sites.vep.vcf.gz -o filtered.vcf -y 0.01  -w 
    (filter variants in input.vcf if present at an allele frequency of 1% or higher in ExAC using INFO fields for filtering)


=head1 AUTHOR

David A. Parry

=head1 COPYRIGHT AND LICENSE

Copyright 2013, 2014, 2015  David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut

my $vcf;
my $filter_vcf;
my @dirs;
my @samples
  ; #samples in vcf input to check calls for - filter only if they contain same allele as in filter_vcf. Default is to simply check alleles in ALT field and ignore samples.
my @reject
  ;    #if specified will only check alleles for these samples in filter_vcfs
my @ignore_samples;    #these samples will be ignored in either VCF
my $filter_with_info
  ; #use allele counts/genotype counts in INFO field to filter with, not samples
my $threshold =
  0;    #only filter if we see the allele this many times in filter_vcf
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
my $MIN_AN ; 
    # min number of alleles required before @pop_acs values are used for filtering if using info
    # or min number of alleles for freq filtering on samples
my @pop_ids =  ();
    # population IDs to look for for @pop_acs
my $annotate_af = '';
my $forks = 0;
my $buffer_size;
my %opts = (
    a              => \$aff_quality,
    annotation     => \$annotate_af,
    b              => \$progress,
    c              => \$buffer_size,
    f              => \$filter_vcf,
    g              => \$minGQ,
    h              => \$help,
    i              => \$vcf,
    l              => \$threshold,
    m              => \$man,
    n              => \$MIN_AN,
    o              => \$out,
    p              => \$print_matching,
    population_ids => \@pop_ids,
    q              => \$min_qual,
    r              => \@reject,
    s              => \@samples,
    t              => \$forks,
    u              => \$unaff_quality,
    w              => \$filter_with_info,
    x              => \@ignore_samples,
    y              => \$maf,
    z              => \$filter_homozygotes,
);

GetOptions(
    \%opts,
    'annotation=s',
    'a|aff_quality=i',
    'b|progress',
    'c|cache=i',
    'f|filter=s{,}',
    'g|genotype_quality=f',
    'h|?|help',
    'i|input=s',
    'l|threshold=i',
    'n|minimum_allele_number=i',
    'm|manual',
    'o|output=s',
    'population_ids=s{,}',
    'p|print_matching',
    'q|quality=f',
    'r|reject=s{,}',
    's|samples=s{,}',
    't|forks=i',
    'u|un_quality=i',
    'w|info_filter',
    'x|not_samples=s{,}',
    'y|allele_frequency_filter=f',
    'z|filter_homozygotes',
) or pod2usage
(
    -message => "Syntax error.",
    -exitval => 2
);

pod2usage( -verbose => 2, -exitval => 0 ) if $man;
pod2usage( -verbose => 1, -exitval => 0 ) if $help;

if ( not $vcf or ( not $filter_vcf ) ){
    pod2usage( -message => "Syntax error", exitval => 2 );
}

if ( @reject and @ignore_samples ){
    pod2usage
    (
        -message => "Syntax error - cannot use --reject ".
                    "and --not_samples argument together",
        exitval  => 2,
    );
} 

if (( @reject or @ignore_samples ) and $filter_with_info){
    pod2usage
    (
        -message => "Syntax error - cannot use --reject or --not_samples ".
                    "arguments with --info_filter option",
        exitval  => 2,
    );
}

if ( $min_qual < 0 ){
    pod2usage
    (
        -message => "Variant quality scores must be 0 or greater.\n",
        -exitval => 2
    );
}
if ( $minGQ < 0 ){
    pod2usage
    (
        -message => "Genotype quality scores must be 0 or greater.\n",
        -exitval => 2
    );
}

if ( $maf < 0 or $maf > 1 ){
    pod2usage
    (
        -message =>  "--allele_frequency_filter (-y) value ".
                     "must be between 0 and 1.\n",
        -exitval => 2
    );
}

if ( defined $aff_quality ) {
    pod2usage
    (
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
    @pop_ids =  qw /POPMAX AFR AMR EAS FIN NFE SAS/ 
    ;#default are gnomAD/ExAC pop codes except for OTH and ASJ
} 
elsif ( grep{ /^disable$/i } @pop_ids ){
    @pop_ids = ();
    my $time = strftime( "%H:%M:%S", localtime );
    print STDERR "[$time] INFO - Filtering on population allele count / allele".
                 " number disabled.\n";
}

if (not defined $MIN_AN){
#if using --info_filter we use require at least 100 alleles by default 
# before filtering on freq
    if ($filter_with_info) {
        $MIN_AN = 100;
    }else{
        $MIN_AN = 0;
    }
}

my $cpus = Sys::CPU::cpu_count();
if ( $forks < 2 ) {
    $forks = 0;    #no point having overhead of forks for one fork
}
else {
    if ( $forks > $cpus ) {
        print STDERR
        my $time = strftime( "%H:%M:%S", localtime );
        print STDERR "[$time] WARNING - Number of forks ($forks) exceeds " .
                     "number of CPUs on this machine ($cpus)\n";
    }
    if ( not $buffer_size ) {
        $buffer_size = 10000 > $forks * 1000 ? 10000 : $forks * 1000;
    }
    my $time = strftime( "%H:%M:%S", localtime );
    print STDERR "[$time] INFO - Processing in batches of $buffer_size ".
                 "variants split among $forks forks.\n";
}

my $time = strftime( "%H:%M:%S", localtime );
my $total_variants = 0;
print STDERR "[$time] INFO - Initializing input VCF...\n";
my ($header, $first_var, $VCF)  = VcfReader::getHeaderAndFirstVariant($vcf);
die "[ERROR] Header not ok for input ($vcf) "
    if not VcfReader::checkHeader( header => $header );
if ( defined $progress ) {
    $time = strftime( "%H:%M:%S", localtime );
    if (-p $vcf or $vcf eq "-" ) {
        print STDERR "\n[$time] INFO - Input is from STDIN or pipe - ".
                     "will report progress per 10000 variants. ";
    }else {
        $total_variants = VcfReader::countVariants($vcf);
        print STDERR "[$time] - INFO $vcf has $total_variants variants.\n";
    }
}
my %sample_to_col = ();
if (@samples) {
    %sample_to_col = VcfReader::getSamples(
        get_columns => 1,
        header      => $header,
    );
    foreach my $samp (@samples){ 
        if (not exists $sample_to_col{$samp}){
            die "ERROR - Sample $samp was not found in input VCF ($vcf)\n"; 
        }
    }
}

$time = strftime( "%H:%M:%S", localtime );
print STDERR "[$time] INFO - Finished initializing input VCF\n";

print STDERR "[$time] INFO - Initializing filter VCF ($filter_vcf)...\n";
my (
    $filter_vcf_index, 
    $filter_vcf_samples, 
    $filter_info
) = initializeFilterVcfs( $filter_vcf );
print STDERR "[$time] INFO - Finished initializing filter VCF.\n";
if ($filter_with_info) {
    check_filter_vcf_info_fields();
}

#only retain samples from filter vcfs if not specified by --not_samples argument
#if --reject argument is used only keep samples specified by --reject
#if neither argument is used keep all
if ( keys %{ $filter_vcf_samples } == 0 ) {
    if (@reject) {
        die
"Filter VCF $filter_vcf has no samples - cannot run with --reject argument.\n";
    }
    if (@ignore_samples) {
        die "Filter VCF $filter_vcf has no samples - cannot run with ".
            "--not_samples argument.\n";
    }
    if ( $threshold or $maf or $filter_homozygotes ) {
        if ( not $filter_with_info ) {
            die "Filter VCF $filter_vcf has no samples. Use of ".
                "--allele_frequency_filter, --threshold or ".
                "--filter_homozygotes options is not allowed with ".
                "filter VCFs without samples unless used with the ".
                "--info_filter option.\n";
        }
    }
}

foreach my $ignore (@ignore_samples) {
    if ( exists $filter_vcf_samples->{$ignore} ) {
        delete $filter_vcf_samples->{$ignore};
    }else{
        die "[ERROR] Sample $ignore not found in --filter VCF ($filter_vcf).\n";
    }
}

foreach my $rej (@reject){
    if (not exists $filter_vcf_samples->{$rej} ) {
        die "[ERROR] Sample $rej not found in --filter VCF ($filter_vcf).\n";
    }
}

if (@reject) {
    foreach my $samp ( keys %{ $filter_vcf_samples } ) {
        if ( not grep { $_ eq $samp } @reject ) {
            delete $filter_vcf_samples->{$samp};
        }
    }
}
my $OUT;
if ($out) {
    open( $OUT, ">$out" ) or die "[ERROR] Can't open $out for writing: $!\n";
}
else {
    $OUT = \*STDOUT;
}

printHeader();

my $prev_chrom = 0;
my $progressbar;
my $next_update = 0;
if ($progress) {
    ($progressbar, $total_variants) = VcfhacksUtils::getProgressBar(
        input  => $vcf,
        name   => "Filtering",
        factor => 3,
    );
}else{
    $time = strftime( "%H:%M:%S", localtime );
    print STDERR "[$time] INFO - Filtering started.\n";
}

my $kept             = 0;
my $filtered         = 0;
my $n                = 0;
my $variants_done    = 0;
my @lines_to_process = ();
my %no_fork_args     = ();

if ( $forks < 2 ) {
    %no_fork_args = VcfReader::getSearchArguments
    ( 
        $filter_vcf, 
        $filter_vcf_index 
    );
}

processLine($first_var);
LINE: while ( my $line = <$VCF> ) {
    processLine($line);
}
process_buffer() if $forks > 1;
if ($progressbar) {
    $time = strftime( "%H:%M:%S", localtime );
    $progressbar->update( $total_variants * 3 )
      if $total_variants * 3 >= $next_update;
    $progressbar->message( "[INFO - $time] $variants_done variants processed" );
}
$time = strftime( "%H:%M:%S", localtime );
print STDERR "[$time] INFO - $filtered variants filtered, $kept printed ";
print STDERR "($total_variants total)" if $total_variants;
print STDERR "\n";

################################################
#####################SUBS#######################
################################################
sub processLine{
    my $line = shift;
    return if $line =~ /^#/;
    $variants_done++;
    $n++;
    checkProgress(1);
    if ( $forks > 1 ) {
        push @lines_to_process, $line;
        checkProgress();
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
            print $OUT join("\t", @$l) . "\n";
            $kept++;
        }
        else {
            $filtered++;
        }
        $n += 2;
        checkProgress();
    }
}

################################################
sub checkProgress{
    return if not $progress;
    my $do_count_check = shift;
    if ($progressbar) {
        $next_update = $progressbar->update($n) if $n >= $next_update;
    }elsif($do_count_check){#input from STDIN/pipe
        VcfhacksUtils::simpleProgress($variants_done, 0, " variants read" );
    }
}

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
                    if (@{ $res{keep} }){
                        $lines_to_print[$res{order}] = \@{ $res{keep} } ;
                    }
                }
                if ( $res{filter} ) {
                    $filtered += $res{filter};
                }
                $n += $res{batch_size};
                checkProgress();
            }
            else {
                die "[ERROR] no message received from child process $pid!\n";
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
    if (@lines_to_print) {
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
    }
    else {
        $n += @lines_to_process;
        checkProgress();
    }
}

################################################
sub process_batch {

    #filter a set of lines
    my ($batch, $order) = @_;
    my %sargs = VcfReader::getSearchArguments( $filter_vcf, $filter_vcf_index );
    my %results = 
    ( 
        batch_size => scalar(@$batch),
        order      => $order,
    );
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
sub printHeader{
    my $meta_head = join("\n", grep {/^##/} @$header);
    print $OUT "$meta_head\n";
    
    foreach my $pop (@pop_acs){
        my $desc = $filter_info->{"AC_$pop"}->{Description} ;
        $desc =~ s/\"/\'/g;
        my $an_desc = $filter_info->{"AN_$pop"}->{Description} ;
        $an_desc =~ s/\"/\'/g;
        my %af_info = 
        (
            ID          => "FVOV_AF_$pop",
            Number      => "A",
            Type        => "Float",
            Description => "Putative population allele number from $filter_vcf.".
                           " Description of original AC_$pop was as follows:" .
                           " $desc",
        );
        my %ac_info = 
        (
            ID          => "FVOV_AN_$pop",
            Number      => "A",
            Type        => "Integer",
            Description => "Putative population allele number from $filter_vcf. ".
                           "Description of original AN_$pop was as follows:".
                           " $an_desc",
        );  
        print $OUT VcfhacksUtils::getInfoHeader(%af_info) . "\n"; 
        print $OUT VcfhacksUtils::getInfoHeader(%ac_info) . "\n"; 
        
    }
    if ($annotate_af){
        $annotate_af =~ s/\s/_/g;#no white space in INFO field
        my %af_info = 
        (
            ID          => "$annotate_af" ."_AF",
            Number      => "A",
            Type        => "Float",
            Description => "Allele frequency calculated by ".
                           "filterVcfOnVcf.pl from $filter_vcf."
        );
        my %ac_info = 
        (
            ID          => "$annotate_af" ."_AN",
            Number      => "A",
            Type        => "Integer",
            Description => "Allele number annotated by ".
                           "filterVcfOnVcf.pl from $filter_vcf.",
        );  
        my $time = strftime( "%H:%M:%S", localtime );
        print STDERR "[$time] INFO - Adding INFO fields $annotate_af" .
                     "_AF and $annotate_af" ."_AN to output for allele ".
                     "frequency and allele numbers calculated from ".
                     "$filter_vcf...\n";
        print $OUT VcfhacksUtils::getInfoHeader(%af_info) . "\n"; 
        print $OUT VcfhacksUtils::getInfoHeader(%ac_info) . "\n"; 
    }
    print $OUT VcfhacksUtils::getOptsVcfHeader(%opts) . "\n"; 
    print $OUT "$header->[-1]\n";
}
################################################
sub filter_on_info_fields {
    my ( $vcf_line, $search_args ) = @_;
    my $qual  = VcfReader::getVariantField( $vcf_line, 'QUAL' );
    my $chrom = VcfReader::getVariantField( $vcf_line, 'CHROM' );
    if ($min_qual &&  $qual < $min_qual ) {
        return;
    }

#process each allele separately in case we have MNVs/deletions 
#that need to be simplified

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
        deleteRefAndNoCallAlleles(\%sample_alleles, \%min_vars);
        return
          if not keys
          %sample_alleles;  #filter if we don't have any variants in our samples
    }
    else {
        #if no samples specified use all alleles for %sample_alleles
        %sample_alleles = map { $_ => undef } keys %min_vars;
        deleteRefAndNoCallAlleles(\%sample_alleles, \%min_vars);
    }

    my %sample_matches   = (); #check each allele matches in all filter_vcfs 
                               # but don't reset after each file
    my %thresh_counts    = (); #count samples in all filter_vcfs 
                               # i.e. don't reset after each file
    my %af_counts;             #count allele occurences and total alleles to calc 
                               # allele frequency and don't reset after each file
    my %f_genos          = (); #store genotypes as keys if using $filter_homozygote
    my %alleles_over_maf = (); #store alleles that have exceeded maf in here

  ALLELE: foreach my $allele ( keys %sample_alleles ) {
        if (
            my @snp_hits = VcfReader::searchForPosition(
                %{$search_args},
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

                my @alts = VcfReader::readAlleles( line => \@snp_split, );
                for ( my $i = 0 ; $i < @alts ; $i++ ) {
                    push @f_alts, $i;
                }
                if ( $filter_homozygotes )
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
                    die "Genotype count (GTC) does not have the same number ".
                        "of alleles as possible genotypes (PGTS) field for ".
                        "line:\n$snp_line\n"
                      if @gtcs != @pgts;
                    for ( my $i = 0 ; $i < @gtcs ; $i++ ) {
                        $geno_counts{ $pgts[$i] } = $gtcs[$i];
                    }
                }
                my $filter_match = ''; #if one of the filter's ALTs matches 
                                       # store the ALT allele code here
ALT:            foreach my $alt (@f_alts) {
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
                    $filter_match = $alt;
                    $sample_matches{$allele}++;
                    last ALT;
                }
                if ( not $filter_match ) {
                    next FILTER_LINE;
                }
                if(@pop_acs){
                    add_pop_freqs_to_allele
                    (
                        $min_vars{$allele}, 
                        \@snp_split, 
                        $filter_match
                    );
                }
                if ($annotate_af){
                    add_annotated_af_to_allele
                    (
                        $min_vars{$allele}, 
                        \@snp_split, 
                        $filter_match
                    );
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
                    if ( exists $filter_info->{AF}) {
                        my @afs = split(
                            ",",
                            VcfReader::getVariantInfoField(
                                \@snp_split, "AF",
                            )
                        );
                        my $freq = $afs[ $filter_match - 1 ];
                        if ( exists $filter_info->{AN}) {
                            my $an = VcfReader::getVariantInfoField
                            (
                                \@snp_split, 
                                "AN"
                            );
                            if ($an >= $MIN_AN){
                                if ( $freq >= $maf ) {
                                    $alleles_over_maf{$allele}++;
                                }
                            }
                        }elsif ( $freq >= $maf ) {
                            $alleles_over_maf{$allele}++;
                        }
                        if(@pop_acs and not $alleles_over_maf{$allele}){
                        #if total AF not over maf check individual populations
                            if ( pop_freq_over_maf($min_vars{$allele}) ) {
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
                        my $freq = 0;
                        eval {
                            $freq =
                              $af_counts{$allele}->{counts} /
                              $af_counts{$allele}->{total};
                        };
                        $alleles_over_maf{$allele}++ if $freq >= $maf;
                    }
                }
            }    #read pos
        }    #search
    }   #foreach allele
        #done each allele - see if we've got enough data to filter this line
    if ( not $maf and not $filter_homozygotes ) { 
    # no filters - default is filter if we found a match
        if (keys %sample_matches == keys %sample_alleles ){#found all alleles
            if ($print_matching) {
                $vcf_line = annotate_pop_freqs(\%min_vars, $vcf_line);
                return $vcf_line;
            }else {
                return;
            }
        }else{#not all alleles found so don't filter
            if (not $print_matching) {
                $vcf_line = annotate_pop_freqs(\%min_vars, $vcf_line);
                return $vcf_line;
            }else {
                return;
            }
        }
    }
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

     #check allele frequency if given
    if ($maf) {
        if ( keys %alleles_over_maf == keys %sample_alleles ) {
            if ($print_matching) {
                $vcf_line = annotate_pop_freqs(\%min_vars, $vcf_line);
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
        $vcf_line = annotate_pop_freqs(\%min_vars, $vcf_line);
        return $vcf_line;
    }
}

#################################################
sub deleteRefAndNoCallAlleles{
    my $al_hash = shift;
    my $min = shift;
    my @del = ();
    foreach my $k (keys %$al_hash){
        if ($k eq '0' or $k eq '.'){
            push @del, $k;
        }elsif($min->{$k}->{ORIGINAL_ALT} eq '*'){
            push @del, $k;
        }
    }
    foreach my $d (@del){
        delete $al_hash->{$d};
    }
}


#################################################
sub pop_freq_over_maf{
    my $min_allele = shift;
    foreach my $pop (@pop_acs){
        next if not $min_allele->{"FVOV_AN_$pop"} ;
        next if $min_allele->{"FVOV_AN_$pop"} eq '.';
        next if $min_allele->{"FVOV_AN_$pop"} < $MIN_AN;
        if ( $min_allele->{"FVOV_AF_$pop"} >= $maf){
            return 1;
        }
    }
    return 0;
}

#################################################
sub annotate_pop_freqs{
    my $m_var = shift;
    my $vcf_line = shift;
    foreach my $pop (@pop_acs){
        my @an = ();
        my @af = ();
        foreach my $allele (sort {$a <=> $b} keys %{$m_var}){
            push @an, 
              exists $m_var->{$allele}->{"FVOV_AN_$pop"} 
              ? $m_var->{$allele}->{"FVOV_AN_$pop"} : '.';
            push @af, 
              exists $m_var->{$allele}->{"FVOV_AF_$pop"} 
              ? $m_var->{$allele}->{"FVOV_AF_$pop"} : '.';
        }
        
        $vcf_line = VcfReader::addVariantInfoField 
            (
                line => $vcf_line,
                id   => "FVOV_AN_$pop",
                value => join(",", @an ), 
            );
        $vcf_line = VcfReader::addVariantInfoField 
            (
                line => $vcf_line,
                id   => "FVOV_AF_$pop",
                value => join(",", @af ), 
            );
    }
    if ($annotate_af){
        my @an = ();
        my @af = ();
        foreach my $allele (sort {$a <=> $b} keys %{$m_var}){
            push @an, 
              exists $m_var->{$allele}->{$annotate_af . "_AN"} 
              ? $m_var->{$allele}->{$annotate_af . "_AN"} : '.';
            push @af, 
              exists $m_var->{$allele}->{$annotate_af . "_AF"} 
              ? $m_var->{$allele}->{$annotate_af . "_AF"} : '.';
        }
        
        $vcf_line = VcfReader::addVariantInfoField 
            (
                line => $vcf_line,
                id   => $annotate_af . "_AN",
                value => join(",", @an ), 
            );
        $vcf_line = VcfReader::addVariantInfoField 
            (
                line => $vcf_line,
                id   => $annotate_af . "_AF",
                value => join(",", @af ), 
            );

    }
    return $vcf_line;
}

#################################################
sub add_annotated_af_to_allele{
    my ($min_allele, $split, $allele_code) = @_;
    # calculate AF for each @pops pop 
    # and return 1 if > $maf
    my $an = VcfReader::getVariantInfoField($split, "AN");
    my $ac = VcfReader::getVariantInfoField($split, "AC");
    my $af = VcfReader::getVariantInfoField($split, "AF");
    if ($an){
        $min_allele->{$annotate_af . "_AN"} = $an;
    }else{
        $min_allele->{"$annotate_af"."_AN"} = '.';
    }
    if (defined $ac){
        my $count = (split ",", $ac)[$allele_code -1];
        $min_allele->{$annotate_af . "_AC"} = $count;
    }else{
        $min_allele->{"$annotate_af"."_AC"} = '.';
    }
    if ($af){
        my $f = (split ",", $af)[$allele_code -1];
        $min_allele->{$annotate_af."_AF"} = $f;
    }elsif($an and defined $ac){
        my $count = (split ",", $ac)[$allele_code -1];
        $min_allele->{$annotate_af . "_AF"} = sprintf("%g", $count/$an);
    }else{
        $min_allele->{"$annotate_af"."_AF"} = '.';
    }
}

#################################################
sub add_pop_freqs_to_allele{
    my ($min_allele, $split, $allele_code) = @_;
    # calculate AF for each @pops pop 
    # and return 1 if > $maf
    foreach my $pop (@pop_acs){
        my $an = VcfReader::getVariantInfoField($split, "AN_$pop");
        if (not $an or $an eq '.'){
            $min_allele->{"FVOV_AN_$pop"} = '.';
            $min_allele->{"FVOV_AC_$pop"} = '.';
            next;
        }
        my @ans = split(",", $an);#some AN (e.g.POPMAX) have Number=A, not 1
        if (@ans > 1){
            $an = $ans[$allele_code -1]
        }
        next if $an eq '.';
        $min_allele->{"FVOV_AN_$pop"} = $an;
        my $acs = VcfReader::getVariantInfoField($split, "AC_$pop");
        my $count = (split ",", $acs)[$allele_code -1];
        $min_allele->{"FVOV_AF_$pop"} = sprintf("%g", $count/$an);
    }
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

#process each allele separately in case MNVs/deletions need simplified

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
        deleteRefAndNoCallAlleles(\%sample_alleles, \%min_vars);
        return
          if not keys
          %sample_alleles;  #filter if we don't have any variants in our samples
    }
    else {
        #if no samples specified use all alleles for %sample_alleles
        %sample_alleles = map { $_ => undef } keys %min_vars;
        deleteRefAndNoCallAlleles(\%sample_alleles, \%min_vars);
    }
    my %sample_matches = (); #check each allele matches in all filter_vcfs 
                             # but don't reset after each file
    my %thresh_counts  = (); #count samples in all filter_vcfs i.e. 
                             # don't reset after each file
    my %f_genos        = (); #store genotypes as keys if using $filter_homozygotes

    my @temp_reject = keys %{ $filter_vcf_samples };
ALLELE: foreach my $allele ( keys %sample_alleles ) {
        if (
            my @snp_hits = VcfReader::searchForPosition(
                %{$search_args},
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

                my @f_alts = ();
                if (@temp_reject){
                    @f_alts = VcfReader::getSampleCall(
                        line                => \@snp_split,
                        multiple            => \@temp_reject,
                        sample_to_columns   => $filter_vcf_samples,
                        return_alleles_only => 1,
                        minGQ               => $unaff_quality
                    );
                }else{#no samples specified - check all alt alleles
                    my @alts = VcfReader::readAlleles( line => \@snp_split, );
                    for ( my $i = 0 ; $i < @alts ; $i++ ) {
                        push @f_alts, $i;
                    }
                }
                if ($filter_homozygotes) {
                    my %genos = VcfReader::getSampleCall(
                        line              => \@snp_split,
                        multiple          => \@temp_reject,
                        minGQ             => $unaff_quality,
                        sample_to_columns => $filter_vcf_samples,
                    );
                    %genos = reverse %genos;
                    foreach my $k ( keys %genos ) {
                        $f_genos{$k}++;
                    }
                }
                my $filter_match = ''; #if one of the filter's ALTs matches 
                                       # store the ALT allele code here
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

                    $filter_match = $alt;
                    $sample_matches{$allele}++;
                    last ALT;
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
                            sample_to_columns   => $filter_vcf_samples,
                          );
                        if ( exists $t_alleles{$filter_match} ) {
                            $thresh_counts{$allele}++;
                        }
                    }

                }
                if ($maf or $annotate_af) {
                    my %f_allele_counts = VcfReader::countAlleles(
                        line                => \@snp_split,
                        minGQ               => $unaff_quality,
                        samples             => \@temp_reject,
                        sample_to_columns   => $filter_vcf_samples,
                    );
                    foreach my $f_al ( keys %f_allele_counts ) {
                        if ( $f_al eq $filter_match ) {
                            $min_vars{$allele}->{AC} = $f_allele_counts{$f_al};
                        }
                        $min_vars{$allele}->{AN} += $f_allele_counts{$f_al};
                    }
                    if ($min_vars{$allele}->{AN} > 0){
                        $min_vars{$allele}->{AF} = sprintf
                        (
                            "%g", 
                            $min_vars{$allele}->{AC} / $min_vars{$allele}->{AN}
                        ); 
                    }else{
                        $min_vars{$allele}->{AN} = 0;
                        $min_vars{$allele}->{AF} = 0;
                    }
                }
            }    #read pos
        }    #search
    }   #foreach allele
    #done each allele 
    if (  not $filter_homozygotes and not $threshold and not $maf ){
        if (keys %sample_matches == keys %sample_alleles ){#found
            return filter_line(\%min_vars, $vcf_line);
        }else{#not found
            if ($print_matching){
                return;
            }else{
                if ($annotate_af){
                    $vcf_line = annotate_sample_freqs
                    (
                        \%min_vars, 
                        $vcf_line, 
                        $annotate_af
                    );
                }
                return $vcf_line ;
            }
        }
    }
        
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
            return filter_line(\%min_vars, $vcf_line);
        }
    }

    if ($threshold) {
        my $alleles_above_threshold = 0;
        foreach my $allele ( keys %sample_alleles ) {
            last if $thresh_counts{$allele} < $threshold;
            $alleles_above_threshold++;
        }
        if ($alleles_above_threshold == keys %sample_alleles ){
            return filter_line(\%min_vars, $vcf_line);
        }
    }
    #check allele frequency if given
    if ($maf) {
        my $alleles_over_maf = 0;
COUNTS: foreach my $allele ( keys %sample_alleles ) {
            last COUNTS if not exists $min_vars{$allele}->{AF};
            if (
                $min_vars{$allele}->{AF} >= $maf and 
                $min_vars{$allele}->{AN} >= $MIN_AN
              )
            {
                $alleles_over_maf++;
            }
            else {
                last COUNTS;
            }
        }
        if ( $alleles_over_maf == keys %sample_alleles ) {
            return filter_line(\%min_vars, $vcf_line);


        }
    }
    if ($print_matching) {
        return;
    }
    else {
        if ($annotate_af){
            $vcf_line = annotate_sample_freqs
            (
                \%min_vars, 
                $vcf_line, 
                $annotate_af
            );
        }
        return $vcf_line;
    }
}

#################################################
sub filter_line{
    if (not $print_matching){
        return;
    }#if using print matching we return the annotated vcf line
    my $m_var = shift;
    my $vcf_line = shift;
    if ($annotate_af){
        $vcf_line = annotate_sample_freqs($m_var, $vcf_line, $annotate_af);
    }
    return $vcf_line;
}
#################################################
sub annotate_sample_freqs{
    my $m_var = shift;
    my $vcf_line = shift;
    my $annot = shift;
    my @an = ();
    my @af = ();
    foreach my $allele (sort {$a <=> $b} keys %{$m_var}){
        if (exists  $m_var->{$allele}->{AN}){
            push @an, $m_var->{$allele}->{AN};
        }else{
            push @an, '.';
        }
        if (exists  $m_var->{$allele}->{AF}){
            push @af, $m_var->{$allele}->{AF};
        }else{
            push @af, '.';
        }
    }
    $vcf_line = VcfReader::addVariantInfoField 
        (
            line => $vcf_line,
            id   => $annot."_AN",
            value => join(",", @an ), 
        );
    $vcf_line = VcfReader::addVariantInfoField 
        (
            line => $vcf_line,
            id   => $annot."_AF",
            value => join(",", @af ), 
        );
    return $vcf_line;
}
#################################################
sub initializeFilterVcfs {
    my ($snpfile) = @_;

    #we getSearchArguments here simply to prevent a race condition later
    my @head = VcfReader::getHeader($snpfile);
    die "[ERROR] Header not ok for $snpfile "
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
    if (   $threshold or $filter_homozygotes)
    {
        check_pgts_gtc($filter_vcf);
    }
    if ( $maf ){
        if ( not exists $filter_info->{AF} ) {
            check_pgts_gtc($filter_vcf);
        }
        @pop_acs = check_pop_ac_info();
        if (@pop_acs){
            my $time = strftime( "%H:%M:%S", localtime );
            print STDERR "[$time] INFO - Found " . scalar(@pop_acs) . 
                " population allele counts in VCF INFO fields. These will be ".
                "used for per population allele-frequency filtering. To " .
                "disable use '--population_ids disable' when running this ".
                "program.\n";
        }
    }
}

#################################################
sub check_pop_ac_info{
    my @pop_codes = ();
    foreach my $pop ( @pop_ids ){
        if (exists $filter_info->{"AC_$pop"} and 
            exists $filter_info->{"AN_$pop"}){
            next if ($filter_info->{"AC_$pop"}->{Number} ne 'A');
            next if ($filter_info->{"AC_$pop"}->{Type} ne 'Integer');
            next if ($filter_info->{"AN_$pop"}->{Number} ne '1' and 
                     $filter_info->{"AN_$pop"}->{Number} ne 'A' );
            next if ($filter_info->{"AN_$pop"}->{Type} ne 'Integer');
            push @pop_codes, $pop;
        }
    }
    return @pop_codes;
}
#################################################
sub check_pgts_gtc {
    if ( not exists $filter_info->{PGTS} ) {
        die <<EOT
Cannot filter on genotype information with the --info_filter option without 
possible genotypes (PGTS) info field. Filter VCF $filter_vcf is missing this 
field from its header. This INFO field is available by processing a VCF 
containing multiple samples using the sampleCallsToInfo.pl script.
EOT
        ;
    }
    elsif ( $filter_info->{PGTS}->{Description} ne
        "\"Possible Genotype Call Codes (for ease of reference).\"" )
    {
        die <<EOT
Filter VCF $filter_vcf possible genotypes (PGTS) field does not match expected 
INFO field description. This INFO field is available by processing a VCF 
containing multiple samples using the sampleCallsToInfo.pl script.
EOT
        ;
    }

    if ( not exists $filter_info->{GTC} ) {
        die <<EOT
Cannot filter on genotype information with the --info_filter option without 
genotype counts (GTC) info field. Filter VCF $filter_vcf is missing this field 
from its header. This INFO field is available by processing a VCF containing 
multiple samples using the sampleCallsToInfo.pl script.
EOT
        ;
    }
    elsif ( $filter_info->{GTC}->{Description} ne
        "\"Genotype counts in order of genotypes listed by PGTS field.\"" )
    {
        die <<EOT
Filter VCF $filter_vcf genotype counts (GTC) field does not match expected INFO 
field description. This INFO field is available by processing a VCF containing 
multiple samples using the sampleCallsToInfo.pl script.
EOT
        ;
    }
}

