#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use List::Util qw ( first ) ;
use POSIX qw/strftime/;
use FindBin;
use Term::ProgressBar;
use lib "$FindBin::Bin/../lib";
use lib "$FindBin::Bin/../lib/dapPerlGenomicLib";
use VcfReader 0.3;
use VcfhacksUtils;

my %opts = (f => 'INFO');
GetOptions(
    \%opts,
    'i|input=s',
    'o|output=s',
    'f|field=s',
    'h|?|help',
);
usage() if $opts{h};
usage("ERROR: -i/--input argument is required") if not $opts{i};

my $OUT = \*STDOUT;
if ($opts{o}){
    open ($OUT, ">", $opts{o}) or die "cannot open $opts{o} for writing: $!\n";
}

my @header = VcfReader::getHeader($opts{i});
print $OUT join("\n", grep { /^##/ } @header) . "\n";
print $OUT "##removeInfoFieldWhitespace.pl -i $opts{i}\n";
print $OUT "$header[-1]\n";

my $VCF = VcfReader::openVcf($opts{i}); 
LINE: while (my $line = <$VCF>){
    next if $line =~ /^#/;#skip header
    chomp $line;
    my @split = split("\t", $line); 
    my $inf = VcfReader::getVariantField(\@split, $opts{f});
    my $r = $inf =~ s/\s/_/g;#replace all white space with underscore
    my $line_ref = \@split;
    if ($r){
        $line_ref = VcfReader::replaceVariantField(\@split, $opts{f}, $inf);
    }
    print $OUT join("\t", @$line_ref) . "\n";
}
close $VCF;


#####################################################################

sub usage{
    my $msg = shift;

    print STDERR "\n$msg\n" if $msg;

    print <<EOT

Replaces any whitespace found in a field of VCF with underscores.
By default replaces whitespace in INFO field, but you may specify the field to change with the --field option.

Usage: $0 -i input.vcf [-f INFO] [-o output.vcf]

Options: 
        -i,--input
            Input vcf
        
        -o,--output
            Optional output file. Default = STDOUT
        
        -f,--field
            VCF field to remove whitespace from. Default = INFO

        -h,--help
            Show this message and exit

EOT
;
    exit 1 if $msg;
    exit;
}
 


