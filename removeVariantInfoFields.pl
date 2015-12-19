#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use FindBin qw($RealBin);
use lib "$RealBin/lib";
use VcfReader;
my @remove = ();
my %opts = (r => \@remove);
GetOptions(
    \%opts,
    'i|input=s', 
    'o|output=s',
    'r|remove_field=s{,}', 
    'h|?|help' ,
) or usage("Syntax error");
usage() if $opts{h};
usage( "--input is required.") if not $opts{i};
usage( "at least one field to remove must be specified using the -r/--remove_field option.") if not @remove;

my @head = VcfReader::getHeader($opts{i});
die "Header not ok for input ($opts{i}) "
    if not VcfReader::checkHeader( header => \@head );
my %info = VcfReader::getInfoFields( header => \@head);

my $OUT = \*STDOUT;
if ($opts{o}){
    open ($OUT, ">$opts{o}") or die "Couldn't open $opts{o} for writing: $!\n";
}
foreach my $h (@head){
    if ($h =~ /^##INFO=<ID=(\w+),Number=([\.\w+]),Type=(\w+),Description=(.+)>$/){
        my $id = $1;
        next if grep {$id eq $_ } @remove;
    }
    print $OUT "$h\n";
}


my $VCF = VcfReader::_openFileHandle( $opts{i} );

while (my $line = <$VCF>){
    next if $line =~ /^#/;
    chomp $line;
    my @s = split("\t", $line);
    my @info = split(";", VcfReader::getVariantField(\@s, 'INFO') ); 
    my @new_info = ();
    foreach my $inf (@info){
        next if grep {$inf eq $_} @remove;
        next if grep {$inf =~ /^$_=.+/} @remove;
        push @new_info, $inf;
    }
    my $replaced = VcfReader::replaceVariantField(\@s, 'INFO', join(";", @new_info));
    print $OUT join("\t", @$replaced) . "\n";
}

#################################################
sub usage{
    
    my $error = shift;
    if ($error){
        print STDERR "ERROR: $error\n";
    }

    print STDERR <<EOT

    Removes specified INFO headers and INFO fields from a VCF

    Usage: 

    $0 -i input.vcf -r INF1 [INF2 INF3 ...]
    Options:

    -i,--input          [vcf input]
    -o,--output         [optional output file]
    -r,--remove_field   [one or more INFO fields to remove]
    -h,-?,--help        [show this help message and exit]
 

EOT
;
    exit 1 if $error;
    exit;
}


