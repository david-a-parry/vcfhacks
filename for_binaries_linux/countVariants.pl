#!/usr/bin/perl
use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/lib";
use VcfReader;

die "Usage: countVariants input1.vcf input2.vcf ...\n" if @ARGV <  1;
my $total = 0;
while (my $vcf = shift){ 
    my $variants = VcfReader::countVariants($vcf) ;
    print "$vcf:\t$variants variants\n";
    $total += $variants;
}
print "Total:\t$total variants\n";
