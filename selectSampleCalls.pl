#!/usr/bin/perl
#David Parry June 2015
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/lib";
use VcfReader;

my @samples = ();
my %opts = (s => \@samples);
GetOptions(
    \%opts,
    'i|input=s', 
    'o|output=s',
    'f|file=s', #list of names
    's|sample=s{,}', 
    'h|?|help' ,
) or usage("Syntax error");
usage() if $opts{h};
usage( "--input is required.") if not $opts{i};
usage( "--sample or --file argument is required.") if not $opts{s} and not $opts{f};

if ($opts{f}){
    open (my $SAMP, $opts{f}) or die "Can't open $opts{f} for reading: $!\n";
    while (my $line = <$SAMP>){
        chomp $line; 
        push @samples, $line;
    }
}

my @header  = VcfReader::getHeader($opts{i});
die "Bad header!\n" if (not VcfReader::checkHeader(header => \@header));
my %samples_to_col = VcfReader::getSamples(
        header => \@header,
        get_columns => 1
);
foreach my $s (@samples){
    if (not exists $samples_to_col{$s}){
        die "ERROR: Sample $s not found in VCF!\n";
    }
}

my $OUT = \*STDOUT;
if ($opts{o}){
    open ($OUT, ">$opts{o}") or die "Couldn't open $opts{o} for writing: $!\n";
}

print $OUT join("\n", @header[0..$#header-1]) . "\n";
print $OUT "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";
print $OUT join("\t", @samples) . "\n";

my $IN = VcfReader::_openFileHandle($opts{i});
while (my $line = <$IN>){
    next if $line =~ /^#/;
    chomp $line;
    my @split = split("\t", $line);
    my @output;
    foreach my $f (qw / CHROM POS ID REF ALT QUAL FILTER INFO FORMAT /){
        push @output, VcfReader::getVariantField(\@split, $f);
    }
    foreach my $s (@samples){
        push @output, VcfReader::getSampleVariant(\@split, $samples_to_col{$s});
    }
    print $OUT join("\t", @output) . "\n";
}
close $IN;
close $OUT;


#################################################
sub usage{
    my $error = shift;
    if ($error){
        print "ERROR: $error\n";
    }
    print <<EOT

    Retrieve variant calls for select samples from a VCF.

    Usage: 

    $0 -i input.vcf -s sample1 [sample2... ] [-o output.vcf]

    Options:

    -i,--input       [vcf input]
    -o,--output      [optional output file. Defaults to printing to STDOUT]
    -s,--sample      [one or more samples to retrieve]
    -f,--file        [text file with one sample per line]
    -h,-?,--help     [show this help message and exit]

    This program will output the calls for samples specified by --sample or --file arguments from a VCF in VCF format. The idea is just to quickly output the call information for a handful of samples from a multisample VCF. However INFO fields will not be modified so INFO fields such as allele counts and frequencies will be inaccurate after processing with this program.

EOT
;
    exit 1 if $error;
    exit;
}


