#!/usr/bin/perl
use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin";
use ParseVCF;

die "Usage: countVariants.pl input1.vcf input2.vcf ...\n" if @ARGV <  1;
my $total = 0;
while (my $vcf = shift){ 
    my $obj = ParseVCF->new(file => $vcf);
    my $variants = $obj->countLines() ;
    print "$vcf:\t$variants variants\n";
    $total += $variants;
}
print "Total:\t$total variants\n";
