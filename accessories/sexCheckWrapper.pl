#!/usr/bin/env perl
use strict; 
use warnings;
use FindBin qw($RealBin);
my $script_prefix = "perl $RealBin/..";

die "Usage: $0 input.vcf\n" if @ARGV != 1;

my $vcf = shift;

my @samples = split("\n", `$script_prefix/getSampleNames.pl -i $vcf -l  -n`); 
checkExit($?);
my @contigs = split("\n", `tabix -l $vcf `);
checkExit($?);
my $x = 'X';
if (grep  {/^chrX$/} @contigs){
    $x = "chrX";
}elsif(not grep  {/^X$/} @contigs){
    die "Could not find 'X' or 'chrX' in list of contigs in $vcf\n";
}
print join("\t", qw ( Sample Het Hom Prop.Het Gender ) ). "\n";
foreach my $s (@samples){
    my $cmd = "tabix -h $vcf $x | $script_prefix/getHetVariants.pl -i - -s $s -gq 30 | $script_prefix/countVariants.pl - ";
    my $het = execute($cmd);
    next if not defined $het;
    $cmd = "tabix -h $vcf $x | $script_prefix/getHetVariants.pl -i - -s $s  -r -gq 30 | $script_prefix/countVariants.pl - ";
    my $hom = execute($cmd);
    next if not defined $hom;
    my $total = $het + $hom;
    my $prop = 0;
    my $gender = '?';
    if ($total){
        $prop = $het/$total;
        $gender = $prop > 0.05 ? 'F' : 'M';
    }
    print join("\t", $s, $het, $hom, $prop, $gender) . "\n";
}
    
    

##################################################################################
sub execute{
    my $cmd = shift;
    print STDERR "Executing command: $cmd\n";
    my $out = `$cmd`;
    checkExit($?);
    if ($out =~ /Total:\s+(\d+)\s+variants/){
        return $1;
    }else{
        warn "Could not parse output from $cmd. Output was:\n$out\n";
    }
    return undef;
}

##################################################################################
sub checkExit{
    my $e = shift;
    if($e == -1) {
        print "Failed to execute: $!\n";
    }elsif($e & 127) {
        printf "Child died with signal %d, %s coredump\n",
        ($e & 127),  ($e & 128) ? 'with' : 'without';
    }elsif($e != 0){
        printf "Child exited with value %d\n", $e >> 8;
    }
}
