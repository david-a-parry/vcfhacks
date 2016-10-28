#!/usr/bin/env perl
#David Parry August 2016
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use FindBin qw ( $RealBin ) ;
use lib "$FindBin::Bin/../lib";
use ParsePedfile;

my %opts = ();

GetOptions
(
    \%opts,
    "p|ped=s",
    "r2=s",
    "r=s",
    "h",
) or usage("Syntax error!");

usage("--ped option is required") if not $opts{p};
usage("--r or --r2 option is required") if not $opts{r} and not $opts{r2};
usage("please only use one of --r or --r2 options") if  $opts{r} and  $opts{r2};

my $p = ParsePedfile->new(file => $opts{p});
my $rel = $opts{r} ? $opts{r} : $opts{r2}; 
open (my $REL, $rel) or die "Error opening $rel: $!\n";
my %scores = ();
while(my $line = <$REL>){
    chomp $line;
    my @split = split /\t/, $line; 
    $scores{$split[0]}->{$split[1]} = $split[-1];
    if ($opts{r}){
        $scores{$split[1]}->{$split[0]} = $split[-1];
    }
}

print join("\t", 
qw/
    Family 
    Sample
    Father
    Mother
    Siblings
    FatherScore
    MotherScore
    SiblingScore
    ValidTrio?
    Culprits?
/
) . "\n";

my $p_cutoff = $opts{r} ? 0.35 : 0.225;
my $s_cutoff = $opts{r} ? 0.3 : 0.2;

foreach my $f ( $p->getAllFamilies() ){
    foreach my $s ( $p->getSamplesFromFamily($f) ){
        if ($p->getParents($s) != 2){
            #print join("\t", $f, $s, split("", "-" x 7), ) . "\n"; 
            next; 
        }
        my @not_ok = ();
        my $father = $p->getFather($s);
        my $mother = $p->getMother($s);
        my $f_score = $scores{$s}->{$father};
        my $m_score = $scores{$s}->{$mother};
        if ($f_score < $p_cutoff){
            push @not_ok, $father;
        }
        if ($m_score < $p_cutoff){
            push @not_ok, $mother;
        }
        my @sib_scores = ();
        my $siblings = join(",", $p->getFullSiblings($s)); 
        $siblings ||= "-";
        foreach my $sib ($p->getFullSiblings($s)){
            push @sib_scores, $scores{$s}->{$sib};
            if ($scores{$s}->{$sib} < $s_cutoff){
                push @not_ok, $sib;
            }
        }
        my $s_score = @sib_scores ? join(",", @sib_scores) : '-';
        my $culprits = @not_ok ? join(",", @not_ok) : "-";
        my $ok = @not_ok ? "FAIL" : "OK";
        print join
        (
            "\t",
            $f,
            $s,
            $father,
            $mother,
            $siblings,
            $f_score,
            $m_score, 
            $s_score,
            $ok,
            $culprits,
        ) . "\n";

    }
}

