#!/usr/bin/env perl
#David Parry June 2015
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use FindBin qw($RealBin);
use lib "$RealBin/lib/dapPerlGenomicLib";
use VcfReader;

my %opts = ();
GetOptions(
    \%opts,
    'i|input=s', 
    'o|output=s',
    's|sample_file=s', #old name -> new name
    'd|delimiter=s',
    'r|rename_duplicates',
    'h|?|help' ,
) or usage("Syntax error");
usage() if $opts{h};
usage( "--input is required.") if not $opts{i};
usage( "--sample_file is required.") if not $opts{s};

my %sample_to_change = readSampleFile($opts{s});

my @header  = VcfReader::getHeader($opts{i});
die "Bad header!\n" if (not VcfReader::checkHeader(header => \@header));
my @samples = VcfReader::getSamples(header => \@header);

my $changed = 0;
my @new_samples = (); 
foreach my $s (@samples){
    if (exists $sample_to_change{$s}){
        push @new_samples, $sample_to_change{$s};
        $changed++;
    }else{
        push @new_samples, $s;
    }
}

checkDuplicates(\@new_samples);


print STDERR "Renamed $changed samples out of " . scalar@samples . "\nWriting output.\n";
my $OUT = \*STDOUT;
if ($opts{o}){
    open ($OUT, ">$opts{o}") or die "Couldn't open $opts{o} for writing: $!\n";
}

print $OUT join("\n", @header[0..$#header-1]) . "\n";
print $OUT "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";
print $OUT join("\t", @new_samples) . "\n";

my $IN = VcfReader::_openFileHandle($opts{i});
while (my $line = <$IN>){
    next if $line =~ /^#/;
    print $OUT $line;
}
close $IN;
close $OUT;

print STDERR "Finished renaming samples.\n";


#################################################
sub readSampleFile{
    my $f = shift;
    my $delimiter = "\t";
    $delimiter = $opts{d} if $opts{d};
    open (my $SAMP, $f) or die "Can't open $f for reading: $!\n";
    my %samps = ();
    while (my $line = <$SAMP>){
        chomp $line;
        next if $line =~ /^#/;
        my @s = split(/$delimiter/, $line);
        if (@s < 2){
            die "\nSample file error - not enough columns for line:\n\n$line\n" .
                "Do you Need to specify your delimiter?\nExiting\n";
        }
        if (exists $samps{$s[0]}){
            die "ERROR: Sample $s[0] appears more than once in $f!\n";
        }
        $samps{$s[0]} = $s[1];
    }
    die "No sample mappings identified in input file $f\n" if not %samps;
    return %samps; 
}

#################################################
sub checkDuplicates{
    my $ar = shift;
    my %seen = ();
    if ($opts{r}){
        renameDuplicates($ar);
    }else{
        my @dups = grep { $seen{$_}++ } @$ar;
        if (@dups){
            die "ERROR: Detected duplicate samples after conversion:\n" 
                . join("\n", @dups) . 
                "\nExiting\nUse the -r/--rename_duplicates option to rename ".
                "duplicate samples rather than exiting.";
        }
    }
}

#################################################a
sub renameDuplicates{
    my $ar = shift;
    my %dups = ();
    foreach my $s (@$ar){
        if (exists $dups{$s}){
            print STDERR "Renaming duplicate sample $s -> ";
            $s = renameSample($s, \%dups);
            print STDERR "$s\n";
        }
        $dups{$s}++;
    }
}

#################################################
sub renameSample{
    my $s = shift;
    my $d = shift;
    my $r = $s . "_$d->{$s}";
    if (exists $d->{$r}){
        $d->{$s}++;
        return renameSample($s, $d);
    }
    return $r;
}
    
#################################################
sub usage{
    my $error = shift;
    if ($error){
        print STDERR "ERROR: $error\n";
    }
    print STDERR <<EOT

    Renames samples in a VCF according to user-supplied sample mapping file.

    Usage: 

    $0 -i input.vcf -s sample_mapping.txt [-o output.vcf]

    Options:

    -i,--input              [vcf input]
    -o,--output             [optional output file]
    -s,--sample_file        [text file with one column for old names and another for new names]
    -d,--delimiter          [string separator for columns in --sample_file. Default is any whitespace]
    -r,--rename_duplicates  [use this flag if duplicate samples should be renamed (e.g. for duplicate Sample_A this would be renamed Sample_A_1)]
    -h,-?,--help            [show this help message and exit]

    Sample mapping file must contain two columns separated by whitespace or a delimiter as specified by the --delimiter option. The first column must be the old sample name and the second the new sample name. Lines beginning with # will be ignored. Extra columns will also be ignored. Any samples in the VCF not present in the sample file will remain unchanged. 

EOT
;
    exit 1 if $error;
    exit;
}


