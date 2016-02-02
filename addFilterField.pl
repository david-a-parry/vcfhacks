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
use lib "$FindBin::Bin/lib";
use VcfReader;
use VcfhacksUtils;

my %opts = (r => 0); 

GetOptions(
    \%opts,
    'i|input=s',
    'o|output=s',
    'f|filter=s',
    'd|description=s',
    'r|replace',
    'm|match_replace=s',
    'h|?|help',
);
usage() if $opts{h};
usage("ERROR: -i/--input argument is required") if not $opts{i};
usage("ERROR: -f/--filter argument is required") if not $opts{f};
usage("ERROR: -d/--description argument is required") if not $opts{d};

my $total_vars = VcfReader::countVariants($opts{i});
print STDERR "$opts{i} has $total_vars variants.\n";
my $OUT = \*STDOUT;
if ($opts{o}){
    open ($OUT, ">", $opts{o}) or die "cannot open $opts{o} for writing: $!\n";
}

my $progressbar = Term::ProgressBar->new(
    {
        name  => "Annotating",
        count => $total_vars,
        ETA   => "linear",
    }
);
my @header = VcfReader::getHeader($opts{i});

print $OUT join("\n", grep { /^##/ } @header) . "\n";
my %filter = 
(
    ID          => $opts{f}, 
    Description => $opts{d},
); 
print $OUT VcfhacksUtils::getFilterHeader(%filter) . "\n"; 
print $OUT "$header[-1]\n";

my $next_update = 0;
my $var_count = 0;
my $VCF = VcfReader::openVcf($opts{i}); 
LINE: while (my $line = <$VCF>){
    next if $line =~ /^#/;#skip header
    chomp $line;
    my @split = split("\t", $line); 
    
    my $line_ref = VcfReader::addVariantFilterField
    (
        line  => \@split,
        id    => $opts{f},
        replace => $opts{r},
        match_replace => $opts{m},
    );
    print $OUT join("\t", @$line_ref) . "\n";
    $var_count++;#for progress bar
    if ($var_count >= $next_update){
        $next_update = $progressbar->update( $var_count )
    }
}
close $VCF;


#####################################################################

sub usage{
    my $msg = shift;

    print STDERR "\n$msg\n" if $msg;

    print <<EOT

Add or replace FILTER field for ALL variants in a VCF - e.g. if you've retrieved 
off target variants and want to label them all as such.

Usage: $0 -i input.vcf -f 'SomeFilter' 'Description of SomeFilter FILTER field'

Options: 
        -i,--input
            Input vcf
        
        -o,--output
            Optional output file. Default = STDOUT
        
        -f,--filter
            Name of FILTER to add. Required

        -d,--description
            Description of FILTER to add. Required

        -r,--replace
            Replace existing FILTER fields rather than prepending to them
        
        -m,--match_replace
            Replace if existing FILTER field matches this string

        -h,--help
            Show this message and exit

Examples:

    $0 -i input.vcf -f OffTarget -d "Off target regions" -o output.vcf
    #prepend 'OffTarget' to FILTER field of all variants in input.vcf

    $0 -i input.vcf -f OffTarget -d "Off target regions" -o output.vcf
    #replace FILTER field with 'OffTarget' for all variants in input.vcf

    $0 -i input.vcf -f OffTarget -d "Off target regions" -o output.vcf -m PASS
    #prepend 'OffTarget' to FILTER field unless existing FILTER field is 'PASS',
    #in which case replace 'PASS' with 'OffTarget'


EOT
;
    exit 1 if $msg;
    exit;
}
 


