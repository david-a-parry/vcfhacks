#!/usr/bin/perl
#writen by  David Parry May 2011
#
##This program is free software: you can redistribute it and/or modify
##it under the terms of the GNU General Public License as published by
##the Free Software Foundation, either version 3 of the License, or
##(at your option) any later version.
##
##This program is distributed in the hope that it will be useful,
##but WITHOUT ANY WARRANTY; without even the implied warranty of
##MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
##GNU General Public License for more details.
##
##You should have received a copy of the GNU General Public License
##along with this program.If not, see <http://www.gnu.org/licenses/>.
##
#
#
sub usage{
    my ($help) = @_;
    print "\nUsage:\t";
    print "findDeNovo.pl -i [vcf input] -s [sample id to search for denovos] -p [parents or other samples to compare] [-a -q [int]]\n\n";
    if ($help){
        print "Identifies potential denovos in vcf file in a specified sample using parents or other samples as filters.\n\n";
        print "--input            -i    [vcf file(s) to analyze]\n";
        print "--sample           -s    [sample to analyze]\n";
        print "--parents          -p    [multiple samples to use to filter out non-denovo mutations]\n";
        print "--all              -a    [output everything that can't be ruled out as de novo rather than just variants which are demonstrably not in both parents]\n";
        print "--quality          -q    [minimum genotype quality confidence score - defaults to 20]\n";
        print "--reject           -r    [reject any variants present in this many additional samples]\n";
        print "--column           -c    [column offset of vcf file (i.e. if additional fields precede the CHROM field)]\n";
        print "--lose_homozygotes -l    [skip any homozygous calls even if not present in parents]\n";
        print "--help             -h    [show this help and exit]\n\n";
    }
    exit;
}
    
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
my @vcf;
my $out;
my $sample;
my %genes;
my %geneline;
my @head;
my $offset = 0;
my $reject = 0;
my $skip_homozygotes;
my @parents;
my $all;
my $help;
my $quality = 20.0; #min genotype quality
my $possible; #flag if we want to keep anything we can't rule out from parents 
GetOptions('output=s' => \$out, 'input=s{,}' =>\@vcf, 'sample=s' => \$sample, 'parents=s{,}' => \@parents, 'column_offset=s' => \$offset, 'reject=i' => \$reject,  'help' => \$help, 'all' => \$possible, "quality=f" => \$quality, "lose_homozygotes" => \$skip_homozygotes);
usage($help) if $help;
usage() if (not @vcf);
my %var_fields;
my $format = 0;
my $header_done = 0;
my $OUT;
if ($out){
    open ($OUT, ">$out") or die "Can't open $out for writing: $!\n";
}else{
    $OUT = \*STDOUT;
}
my @other_samples;#indexes of positions of other samples in vcf file for use with reject option
foreach my $vcf (@vcf){
    open (my $IN, $vcf) || die "Can't open $vcf\n";
    LINE: while (<$IN>){
        if (/^##/){
            print $OUT $_;
            next LINE;
        }
        if (/^#[A-Z]/){
            print $OUT $_ if not $header_done;
            $header_done++;
            chomp (my @header = split("\t"));
            if ($sample){
                $var_fields{$sample} = 0;
                $var_fields{$sample}++ until $header[$var_fields{$sample}] eq $sample or $var_fields{$sample} > $#header;
                die "Can't find $sample in $vcf header\nHeader:$_" if $var_fields{$sample} > $#header;
            }
            if (@parents){
                foreach my $parent (@parents){
                    $var_fields{$parent} = 0;
                    $var_fields{$parent}++ until $header[$var_fields{$parent}] eq $parent or $var_fields{$parent} > $#header;
                    die "Can't find $parent in $vcf header\nHeader: $_" if $var_fields{$parent} > $#header;
                }
            }
            $format++ until $header[$format] eq "FORMAT" or $format > $#header;
            die "Can't find FORMAT field in $vcf header\nHeader: $_" if $format > $#header;
            my $field = $format+1;
OTHER_SAMPLES:        while ($field < @header){
                my $other_sample = $header[$field];
                $field++;
                next OTHER_SAMPLES if $other_sample eq $sample;
                foreach my $parent (@parents){
                    next OTHER_SAMPLES if $other_sample eq $parent;
                }
                #$other_sample is not parent or our test sample if we exit this loop
                push (@other_samples, $field-1);#index of positions of other samples
            }
            next LINE;
        }
        my %gene;
        chomp (my @line = split("\t"));
        foreach my $line (@line){
            $line=~ s/[\r\s\n]//;
        }
        next LINE if $line[$var_fields{$sample}] =~ /^\.[\/\|]\./;
        my @format_field = split(":", $line[$format]);
        my $gt = 0;
        $gt++ until $format_field[$gt] eq "GT" or $gt > $#format_field;
        die "Can't find genotype field in vcf line:\n$_" if $gt > $#format_field;
        my $gq = 0;
        $gq++ until $format_field[$gq] eq "GQ" or $gq > $#format_field;
        die "Can't find genotype quality field in vcf line:\n$_" if $gq > $#format_field;
        #SEE IF SAMPLE IS HET AND IF SO SEE IF PARENTS ARE WT
        my @sample_gt = split(":", $line[$var_fields{$sample}]);
        my $no_call = 0;
        if ($sample_gt[$gt] =~ /(\d+)[\/\|](\d+)/){#catch alleles with $1 and $2
            next LINE if $1 == 0 and $2 == 0;#if sample is hom reference then skip
            if ($skip_homozygotes){
                next LINE if $1 == $2;
            }
            my %sample_alleles;
            $sample_alleles{$sample}->{$1}++;
            $sample_alleles{$sample}->{$2}++;
            foreach my $parent (@parents){
                my @par_gt = split(":", $line[$var_fields{$parent}]);
                my @par_alleles = split(/[\/\|]/, $par_gt[$gt]);
                if ($par_gt[$gt] =~ /\.[\/\|]\./){
                    next LINE if not $possible;#skip no calls unless $possible flag is set
                }else{#there will only be a quality field if there is a genotype call
                    next LINE if $par_gt[$gq] < $quality;#skip low quality
                }
                foreach my $par_allele (@par_alleles){
                    $sample_alleles{$parent}->{$par_allele}++;
                }
            }#test alleles caught in our %sample_alleles hash
ALLELE:            foreach my $key (keys %{$sample_alleles{$sample}}){#iterate through our samples alleles
                next ALLELE if $key == 0;#skip reference calls
                foreach my $parent (@parents){
                    next ALLELE if (exists $sample_alleles{$parent}->{$key});
                }
                #if we've exited the foreach @parent loop then the allele is not present in any of the parent samples 
                if ($reject){
                    my $match = 0;
REJECT:                    foreach my $other_sample (@other_samples){
                        my @other_sample_gt = split(":", $line[$other_sample]);
                        next REJECT if ($other_sample_gt[$gt] =~ /\.[\/\|]\./);
                        next REJECT if $other_sample_gt[$gq] < $quality;#skip low quality
                        my @other_sample_alleles = split(/[\/\|]/, $other_sample_gt[$gt]);
                        foreach my $other_sample_allele (@other_sample_alleles){
                            $match++ if $key == $other_sample_allele;
                        }
                    }
                    next ALLELE if $match >= $reject;#skip if allele is present in $reject samples or more
                }
                #if we've exited or skipped this block then allele doesn't exceed $reject threshold in other samples
                print $OUT  $_;
                next LINE;
            }
        }
    }
}


###########
