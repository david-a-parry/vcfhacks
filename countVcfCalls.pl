#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
#use Tie::File;
use Net::FTP;
use File::Temp qw/ tempfile tempdir /;
my $vcf;
my $help;
my $depth;
my $total;
my $quality;
my $out;
my $get_sample_names;
my $exclude_non_variants;
GetOptions("exclude_non_variants" => \$exclude_non_variants, "ouput=s" => \$out, "input=s" => \$vcf, "depth=i" => \$depth, "quality=i" => \$quality, "help" => \$help, "retrieve_samples" => \$get_sample_names)
    or usage("Syntax error");
usage($help) if ($help);
usage() unless($vcf);
open (my $IN, $vcf) || die "can't open $vcf:$!\n";
my $OUT;
if ($out){
    open ($OUT, ">$out") || die "can't open $out for writing:$!\n";
}else{
    $OUT = \*STDOUT;
}
my ($tmp, $TEMP);
my $tmp_prefix = 'count';
$tmp_prefix = $out if $out;
($TEMP, $tmp) = tempfile("$tmp_prefix.tmp_countVcfXXXX", UNLINK => 1) or die "Can't create temporary output file\n";
my $dp_field = 0;
my $gq_field = 0;
my $gt_field = 0;
my $format = 8;
my @samples = ();
my @header = ();
my $chrom = 0;
my $num_sample_genotype_fields = 0; 
print STDERR "Reading input, counting variants...\n";
while (<$IN>){
    next if (/^##/);
    chomp;
    if (/^#\w/){
        if (@header){
            die "Multiple header lines found - please check your input!\n";
        }
        @header = split("\t");
        $chrom++ until $header[$chrom] =~ /CHROM/ or $chrom > $#header;
        die "Can't find CHROM field in header\n" if $chrom > $#header;
        $format = 0;
        $format++ until $header[$format] eq "FORMAT" or $format > $#header;
        die "Can't find format field in header.\n" if ($format > $#header);
        if ($get_sample_names){
            @samples = @header[($format+1)..$#header];
            die "No samples in header!\n" if not @samples;
        }
        print $TEMP join("\t", @header[0..$chrom+4]);
        print $TEMP "\tAllele Counts\tGenotype Counts";
        print $TEMP "\n";
    }else{
        if ($get_sample_names and not @header){
            print STDERR "WARNING - can't retrieve sample names - no header found.\n" if not @header;
            $get_sample_names = 0;
        }
        my @line = split(/\t/, $_);
        my %variant_samples;
        $dp_field = 0;
        $gq_field = 0;
        $gt_field = 0;
        my @fields = split(/:/, $line[$format]);
        $dp_field++ until $fields[$dp_field] =~ "DP" or $dp_field >= $#fields;
        if ($depth and $fields[$dp_field] ne "DP"){
            print STDERR "Warning - no depth field found\n";
            $depth = 0;
        }
        $gq_field++ until $fields[$gq_field] eq "GQ" or $gq_field >= $#fields;
        if ($quality and $fields[$gq_field] ne "GQ"){
            print STDERR "Warning - no genotype quality field found\n";
            $quality = 0;
        }
        $gt_field++ until $fields[$gt_field] eq "GT" or $gt_field >= $#fields;
        die "No genotype field\n$_" if $fields[$gt_field] ne "GT";
        
        my $ref = $line[$chrom + 3];
                my @alts = split(",", $line[$chrom + 4]);
        @alts = grep {! /\./} @alts;
                my @all_alleles = ($ref, @alts);
                my @calls = @line[($format+1)..$#line];        
        my $refs = 0;
        my $no_calls = 0;
                my %other_alleles = ();
                my %genotypes = ();
                for (my $i = 0; $i < @calls; $i++){
            my @sample_call = split(/:/, $calls[$i]);
                    my $gt = $sample_call[$gt_field];
            if ($gt =~ /\.[\/\|]\./){
                $no_calls++;
                next;
            }
            if ($depth){
                next if ($sample_call[$dp_field] < $depth);
            }
            if ($quality){
                next if ($sample_call[$gq_field] < $quality);
            }
            
                        my @alleles = split(/[\/\|]/, $gt);
            my $sample_has_variant = 0;#we'll collect all samples with variants if $get_sample_names flag is true
                        my @sorted_genotype = ();#the order of alleles varies in 1000 genomes vcf files
                        foreach my $allele (sort {$a<=>$b} @alleles){
                             if ($allele == 0){
                                    push (@sorted_genotype, $all_alleles[0]);
                                        $refs++;
                                }elsif ($allele =~ /(\d+)/){
                                        push (@sorted_genotype, $all_alleles[$1]);
                                        $other_alleles{$all_alleles[$1]}++;
                    $sample_has_variant++;
                                }
                        }
                        $genotypes{join("/", @sorted_genotype)}++; 
            if ($sample_has_variant and $get_sample_names){
                push(@{$variant_samples{join("/", @sorted_genotype)}}, $samples[$i]);
            }
                }
        if ($exclude_non_variants){
            if (keys %genotypes == 1){
                next if (keys %genotypes)[0] eq "$all_alleles[0]/$all_alleles[0]";
            }
        }
                print $TEMP join("\t", @line[0..$chrom+4]);
                print $TEMP "\t$refs reference";
                foreach my $key (keys %other_alleles){
                    print $TEMP ", $other_alleles{$key} $key";
                }
                print $TEMP "\t";
                my @possible_genotypes = get_possible_genotypes(\@all_alleles);
                my @genotype_string = ();
        foreach my $genotype (@possible_genotypes){
                    if (exists $genotypes{$genotype}){
                        push (@genotype_string, "$genotype $genotypes{$genotype}");
                        }else{
                        push (@genotype_string, "$genotype 0");
                        }
                }
                print $TEMP join(", ", @genotype_string);
        my @sample_string = ();
        if ($get_sample_names){
            foreach my $genotype (@possible_genotypes[1..$#possible_genotypes]){#skipping the first (hom reference) genotype
                if (exists $variant_samples{$genotype}){
                    my $string = "Samples with $genotype: " . join(", ", @{$variant_samples{$genotype}});
                    push(@sample_string, $string);
                }else{
                    push(@sample_string, "-");
                }
            }    
            print $TEMP "\t". join("\t", @sample_string) ;
            $num_sample_genotype_fields  = @sample_string if $num_sample_genotype_fields < @sample_string; 
        }
        print $TEMP "\n"; 
    }
}
close $TEMP;
print STDERR "Done counting variants - writing output...\n";
open ($TEMP, $tmp) or die "Can't open temp file $tmp for reading and writing output.\n";
my $total_fields = $chrom + 7 + $num_sample_genotype_fields;
while (<$TEMP>){
    if (/^#\w/){
        chomp;
        print $OUT $_;
        if ($get_sample_names){
            for (my $n = 1; $n <= $num_sample_genotype_fields; $n++){
                print $OUT "\tAlt Genotype $n";
            }
        }
        print $OUT "\n";
    }else{
        my @split = split("\t");
        if (@split < $total_fields){
            chomp;
            my $extra = $total_fields - @split;
            print $OUT $_;
            print $OUT "\t-" x $extra;
            print $OUT "\n";
        }else{
            print $OUT $_;
        }    
    }
}
        
print STDERR "Done writing output\n";

my $prob = 10**(-$quality/10)  if ($quality);

print STDERR "Parameters: ";
print STDERR "Minimum depth = $depth " if ($depth);
print STDERR "Minimum per sample phred genotype quality = $quality (P = $prob)" if ($quality);
print STDERR "No minimum depth or quality specified" unless ($depth or $quality);
print STDERR "\n";

close $IN;


####
sub get_possible_genotypes{
        my ($allele_ref) = @_;
        my @combinations = ();
        for (my $n = 0; $n < @$allele_ref; $n++){
                for (my $m = 0; $m <= $n; $m++){
                        push (@combinations, "$$allele_ref[$m]/$$allele_ref[$n]");
                }
        }
        return @combinations;
}


####
sub usage{
my $help = @_;
if ($help){
    print <<EOT

    Reads VCF files and reports the number of observed genotypes for each variant.
    
    Usage:

    $0 -i input.vcf [options]

    Options:

    -i,--input                [vcf file]
    -d,--depth                [minimum depth of reads in order to count variant - optional]
    -q,--quality              [minimum phred quality score for genotype - optional]
    -e,--exclude_non_variants [only include lines with at least one sample with a variant call]
    --r,--retrieve_samples    [print columns with sample names per genotype]

EOT
;
    exit;
}
else{
    print "--input argument is required.  For help use --help\n";
    exit 1;
}
}
