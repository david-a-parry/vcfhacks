#!/usr/bin/perl
use warnings;
use strict;
use FileHandle;
use Getopt::Long;
use Pod::Usage;
use Term::ProgressBar;
use Data::Dumper;
use POSIX qw/strftime/;
use FindBin;
use lib "$FindBin::Bin/lib";
use ParseVCF;
my @samples;
my @filters;
my $input; 
my %opts = (samples => \@samples, filters => \@filters);

GetOptions(\%opts,
        "output=s",
        "input=s",
        "known_snps", # i.e. only process lines with rs ID
        "samples=s{,}",
        "quality=i",#phred like genotype quality required
        "variant_quality=f",
        "genome_build=s",
        "filters=s{,}",
        "force",
        "help",
        "manual",
        "progress",
        ) or pod2usage(-exitval => 2, -message => "Syntax error") ;
pod2usage(-verbose => 2) if $opts{manual};
pod2usage(-verbose => 1) if $opts{help};
pod2usage(-exitval => 2, -message => "Syntax error: --input argument is required") if not $opts{input};

my @acceptable_chromosomes = (1 .. 22, "X", "Y", "MT"); 
my %chrom_check = map {$_ => 1} @acceptable_chromosomes;
my %filter_check = map {$_ => 1} @filters;

#open input vcf with ParseVCF
print STDERR "Initialising VCF...\n";
my $vcf_obj = ParseVCF->new(file=> $opts{input});
my $total_vcf = $vcf_obj->countLines("variants") if defined $opts{progress};
print STDERR "$opts{input} has $total_vcf variants\n" if defined $opts{progress};
my $time = strftime("%H:%M:%S", localtime);

if (not @{$opts{samples}}){
#default is to output all samples
    @samples = $vcf_obj->getSampleNames();
    die "No samples found in VCF input $opts{input}!\n" if not @samples;
}else{
#for each sample check it exists in input vcf 
    foreach my $sample (@samples){
    if (not grep{$sample eq $_} $vcf_obj->getSampleNames()){
            die "Can't find sample $sample in $opts{input}. Found these samples: "
            .join(", ", $vcf_obj->getSampleNames()) . "\n";
        }
    }
}
#create an output file per sample, open filehandle and  
#put filehandles in a hash (key = sample). Write header. 
my %filehandles = ();
$opts{genome_build} = 'hg19' if not $opts{genome_build};
foreach my $sample (@samples){
    my $out;
    if ($opts{output}){
    $out = $opts{output} . "_$sample.snpview.txt";
    }else{
    $out = "$sample.snpview.txt";
    }
    if (-e $out and not $opts{force}){
        die "$out already exists - choose a new name or use --force to overwrite\n";
    }
    my $fh = FileHandle->new("> $out");
    $fh->print("#VCF File=$opts{input}\n");
    $fh->print("#Converted $0 " .join(" ", @ARGV) ."\n");
    $fh->print("#%genome-version-ucsc=$opts{genome_build}\n");
    $fh->print("Probe Set ID\tCall Codes\tConfidence\tdbSNP RS ID\tChromosome\tChromosomal Position\n");
    $filehandles{$sample} = $fh;
}

#set up progressbar 
my $progressbar;
my $n = 0;
my $next_update = 0;
if (defined $opts{progress} and $total_vcf){
    $progressbar = Term::ProgressBar->new(
    {
        name => "Converting", 
        count => $total_vcf, 
        ETA => "linear", 
    }
    );
}

#for each line of vcf file
#check there are only 2 alleles
#for each sample get genotype call and write line to sample's output file
my $converted = 0;
VAR: while (my $line = $vcf_obj->readLine){
    $n++;
    if ($opts{known_snps}){
        if (defined $progressbar){
            $next_update = $progressbar->update($n) if $n >= $next_update;
        }
        next VAR if $vcf_obj->getVariantField("ID") !~ /rs\d+/i;
    }
    if (defined $opts{variant_quality}){
        next VAR if $vcf_obj->getVariantField("QUAL") < $opts{variant_quality};
    }
    if (@filters){
       my $f = $vcf_obj->getVariantField("FILTER");
       next VAR if not exists $filter_check{$f}; 
    }
    next VAR if $vcf_obj->readAlleles() > 2;#skip line if not biallelic
    my $chrom = $vcf_obj->getVariantField("CHROM");
    $chrom =~ s/^chr//i;
    next VAR if not exists $chrom_check{$chrom};
    my $pos = $vcf_obj->getVariantField("POS");
    my $line_id = $chrom ."_" . $pos . "_" 
                . $vcf_obj->getVariantField("REF") . "/" 
                . $vcf_obj->getVariantField("ALT") ;
    
SAMPLE: foreach my $sample (@samples){
        my $call;
        my @alleles = $vcf_obj->getSampleActualGenotypes
                            (
                                sample => $sample, 
                                return_alleles_only => 1
                            );
        next SAMPLE if not @alleles;#skip if no call
        if (@alleles > 2){
            die "Internal Error - More than 2 alleles, should have skipped. At line:\n$line\n";
        }elsif (@alleles == 2){
            $call = 'AB';
        }else{
            if ($alleles[0] eq $vcf_obj->getVariantField("REF")){
                $call = 'AA';
            }elsif($alleles[0] eq $vcf_obj->getVariantField("ALT")){
                $call = 'BB';
            }else{
                die "Couldn't determine allele for $sample line:\n$line\n";
            }
        }
        my $gq =  $vcf_obj->getSampleGenotypeField
                        (
                            field => 'GQ',
                            sample => $sample
                        );
        if (defined $opts{quality}){#filter on genotype phred like quality score
            next SAMPLE if not defined $gq;
            next SAMPLE if $gq < $opts{quality};
        }
        my $genotype_quality;
        if (defined $gq){
            $genotype_quality = 10**(-$gq/10);
        }else{
            $genotype_quality = 1.000;
        }
        my $id = $vcf_obj->getVariantField("ID") ;
        if ($id eq '.'){
            $id = $line_id;
        }
        $filehandles{$sample}->print("$line_id\t$call\t$genotype_quality\t$id\t$chrom\t$pos\n");
        if (defined $progressbar){
            $next_update = $progressbar->update($n) if $n >= $next_update;
        }
    }
}

if (defined $progressbar){
    $progressbar->update($total_vcf) if $total_vcf >= $next_update;
}
#tidy up
foreach my $sample (@samples){
    $filehandles{$sample}->close;
}

$vcf_obj->DESTROY();

=head1 NAME

vcfToSnpViewer.pl - convert variants from a VCF file to per sample input files for SnpViewer (http://sourceforge.net/projects/snpviewer/)

=head1 SYNOPSIS

        vcfToSnpViewer.pl -i [vcf file]  [options]
        vcfToSnpViewer.pl -h (display help message)
        vcfToSnpViewer.pl -m (display manual page)

=cut

=head1 ARGUMENTS

=over 8

=item B<-i    --input>

Input VCF file 

=item B<-o    --output>

Prefix for output files. Files will be named this prefix_samplename.snpviewer.txt or if this argument isn't used output files will simply be samplename.snpviewer.txt

=item B<-k    --known_snps>

Only output variants with a dbSNP rs identifier in the ID field.

=item B<-s    --samples>

One or more samples to ourput variants for. Default is to output all samples present in input VCF. 

=item B<-q    --quality>

Minimum phred-like sample genotype quality to output.

=item B<-v    --variant_quality>

Minimum phred-like variant quality (QUAL field) to output.

=item B<-g    --genome_build>

Genome build to annotate output files with.  Defaults to hg19.

=item B<--filter>

Value in FILTER field to select variants with. Only output variants with a FILTER field exactly matching one of these terms.

=item B<--force>

Use this flag to overwrite existing output files.

=item B<-p    --progress>

Use this flag to display progress bar.

=item B<-h    --help>

Display help message.

=item B<-m    --manual>

Show manual page.

=back 

=cut


=head1 DESCRIPTION

Reads a VCF file and outputs variants in a format readable by my autozygosity mapping tool SnpViewer (http://sourceforge.net/projects/snpviewer/). If you want fine control of the types of SNPs you output you may want to run annotateSnps.pl to get a list of known SNPs matching your criteria. This program allows filtering on the basis of genotype quality scores, variant quality scores, filters, or SNP rs ID.


=cut

=head1 AUTHOR

David A. Parry

=head1 COPYRIGHT AND LICENSE

Copyright 2013  David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut

