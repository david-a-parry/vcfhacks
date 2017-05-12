#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Term::ProgressBar;
use File::Temp qw/ tmpnam / ;
use Data::Dumper;
use FindBin qw($RealBin);
my $script_prefix = "perl $RealBin/..";

my $loc = tmpnam(); 
my $func = tmpnam(); 
my @genes = ();
my @classes = ();
my %opts = 
(
    g => \@genes,
    b => 'GRCh37',
    a => 0.01,
    v => \@classes,
);
GetOptions(
        \%opts,
        'i|input=s',
        'g|genes=s{,}',
        'o|output=s',
        'v|add_variant_classes=s{,}',
        'f|force',
        'b|build=s',
        'a|allele_frequency=f',
        'h|?|help',
) or die "Syntax error in option spec!\n";

usage() if $opts{h};
usage("--input argument is required!\n") if not $opts{i};
usage("--genes argument is required!\n") if not @genes;

my $useHg38 = checkBuild();
if ($useHg38){
    print STDERR "INFO: Using build 38 of the human genome\n";
}else{
    print STDERR "INFO: Using build 37 of the human genome\n";
}

my $output;
if ($opts{o}){
    $output = $opts{o};
}else{
    $output = "geneSummary.xlsx";
}
$output .= ".xlsx" if $output !~ /\.xlsx$/;
if (-e $output and not $opts{f}){
    die "$output already exists - use the --force flag to overwrite\n";
}


my $cmd = "$script_prefix/getVariantsByLocation.pl -g " . join(" ", @genes) . " -i $opts{i} -o $loc";
if ($useHg38){
    $cmd .= " -d";
}
print STDERR "Executing command: $cmd\n";
system($cmd);
die "getVariantsByLocation.pl Error $? / $@\n" if $?;

$cmd = "$script_prefix/getFunctionalVariants.pl -i $loc -o $func -add splice_region_variant";
if ($opts{a}){
     $cmd .= " -a $opts{a}";
}
if (@classes){
    $cmd .= " --add_classes " . join(" ", @classes) ;
}
print STDERR "Executing command: $cmd\n";
system($cmd);
die "getFunctionalVariants.pl Error $? / $@\n" if $?;

$cmd = "$script_prefix/annovcfToSimple.pl -v -u  all --contains_variant --functional -i $func -o $output -n CaddPhredScore ";
foreach my $af (qw /
    FVOV_AF_AFR
    FVOV_AF_AMR
    FVOV_AF_EAS
    FVOV_AF_FIN
    FVOV_AF_NFE
    FVOV_AF_SAS
    AS_CLNDBN
    AS_CLNSIG
    AS_COMMON
    AS_dbSNPBuildID
    AS_G5
    AS_G5A
    AS_CAF
    AS_MAF
/){
    $cmd .= " $af";
}

print STDERR "Executing command: $cmd\n";
system($cmd);
die "annovcfToSimple.pl Error $? / $@\n" if $?;

END{
    for my $f ($loc, $func){
        if (-e $f){
           unlink($f);
        }
    }
    print STDERR "Done\n";
}



##################################################
sub checkBuild{
    if (uc($opts{b}) eq 'GRCH38' or uc($opts{b}) eq 'HG38'){
        return 1;
    }elsif(uc($opts{b}) eq 'GRCH37' or uc($opts{b}) eq 'HG19'){
        return 0;
    }else{
        die "Did not recognise genome build '$opts{b}' supplied by -b/--build"
            ."argument. Valid values are GRCh37, GRCh38, hg19 or hg38.\n";
    }
}

##################################################
sub usage{
    my $msg = shift;
    print STDERR "\n$msg\n" if $msg;
    print <<EOT

Usage: $0 -i input.vcf -o output.xlsx -g gene1 [gene2 ... geneN] 

Options:
    -i,--input FILE
        Input VCF file.  Required.

    -g,--genes STRING(s)
        One or more genes to find overlapping variants for. Required.
    
    -b,--build STRING
        Use this human genome build. Valid values are GRCh37, GRCh38, hg19, hg38.

    -a,--allele_frequency FREQ
        Value between 0.00 and 1.00 for filtering variants on allele frequency.
        Variants with a frequency greater than this value will be filtered. 
        Default = 0.01

    -o,--output FILE
        Output XLSX file. Default = geneSummary.xlsx.

    -v,--add_variant_classes
        Specify additional VEP variant classes to retrieve (e.g. intron_variant) 
        in addition to the default. By default the following classes are 
        retrieved:
                    frameshift_variant
                    inframe_deletion
                    inframe_insertion
                    initiator_codon_variant
                    missense_variant
                    protein_altering_variant
                    regulatory_region_ablation
                    regulatory_region_amplification
                    splice_acceptor_variant
                    splice_donor_variant
                    splice_region_variant
                    stop_gained
                    stop_lost
                    transcript_ablation
                    transcript_amplification
                    TFBS_ablation
                    TFBS_amplification


    -f,--force
        Use this flag to overwrite existing output files.

    -h,--help
        Show this message and exit

EOT
;
    exit 1 if $msg;
    exit;
}





