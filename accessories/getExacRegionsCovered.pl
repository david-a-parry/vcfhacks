#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

my @input = ();
my %opts = 
(
    i => \@input, 
    c => 20,
    f => 0.85,
);

GetOptions
(
    \%opts,
    'c|covered=i',
    'f|fraction_over=f',
    'i|input=s{,}',
    'h|help',
) or usage("Syntax error in option spec\n");
usage() if $opts{h};
usage("-i/--input argument is required\n") if not @input; 
my %cov_columns = ();
foreach my $f (@input){
    my $IN = openCovFile($f);
    my $head = <$IN>;
    chomp $head;
    close $IN;
    if ($head !~ /^#chrom\tpos\tmean\tmedian(\t(\d+))+$/){
        die "Invalid header for file '$f'!\n";
    }
    my @split = split("\t", $head); 
    my $i = 0;   
    $i++ until $split[$i] eq $opts{c} or $i > $#split;
    if ( $i > $#split ){
        die "Could not find coverage column '$opts{c}' for file $f\n";
    }
    $cov_columns{$f} = $i;
}

foreach my $f (@input){
    my %region = (); 
    my $IN = openCovFile($f);
    my $prev_chrom = '';
    while (my $line = <$IN>){
        next if $line =~ /^#/;
        chomp $line;
        my @split = split("\t", $line); 
        my $fract = $split[ $cov_columns{$f} ];
        if ($region{chrom} and $region{chrom} ne $split[0]){
            outputRegion(\%region, \@split);
        }
        if ($fract >= $opts{f}){
            if (not keys %region){
                $region{start} = $split[1];
                $region{last} = $split[1];
                $region{chrom} = $split[0];
            }elsif($region{last} == $split[1] - 1){
                $region{last} = $split[1];
            }else{
                outputRegion(\%region, \@split);
            }
        }

    }
    close $IN;
    outputRegion(\%region);
}

sub outputRegion{
    my $region = shift;
    my $split = shift;
    if (keys %{$region}){
        my $bedstart = $region->{start} - 1;
        print "$region->{chrom}\t$bedstart\t$region->{last}\n";
    }
    if ($split){
        $region->{start} = $split->[1];
        $region->{last} = $split->[1];
        $region->{chrom} = $split->[0];
    }
}

sub openCovFile{
    my $f = shift;
    if ($f =~ /\.gz$/){
        open (my $IN, "gzip -dc $f |") or die "Can't open input '$f' via gzip: $!\n";
        return $IN;
    }else{
        open (my $IN, "<", $f) or die "Can't open input '$f': $!\n";
        return $IN;
    }
}

sub usage{
    my $msg = shift;
    print STDERR "\n$msg\n" if $msg;
    print STDERR <<EOT

Retrieves regions from ExAC coverage files where a given fraction of samples are covered by more N reads or more.

Usage: $0 -i [one or more ExAC coverage files] [options]

Options:
    
    -i,--input FILE(s)
        One or more coverage files from ExAC

    -c,--covered INT
        Minimum coverage for output. Default = 20X.

    -f,--fraction_over FLOAT
        Fraction of samples to be covered creater than the minimum coverage.
        Default = 0.85.
    
    -h,--help
        Show this message and exit

EOT
    ;
    exit 1 if $msg;
    exit;
}


