#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long qw(:config no_ignore_case);
use Data::Dumper;
use POSIX qw/strftime/;
use FindBin qw($RealBin);
use lib "$RealBin/../lib";
use lib "$RealBin/../lib/dapPerlGenomicLib";
use VcfReader 0.3;
use VcfhacksUtils;

my $minGq = 20;
my %impacts = 
( 
    HIGH     => 0,
    MODERATE => 1,
    LOW      => 2,
    MODIFIER => 3,
);

die "Usage: $0 input.vcf\n" if @ARGV !=1;
my $input = shift;
my $time = strftime( "%H:%M:%S", localtime );
print STDERR "[$time] INFO - Initializing input VCF...\n";
my ($header, $first_var, $VCF)  = VcfReader::getHeaderAndFirstVariant($input);
die "Header not ok for input ($input) "
    if not VcfReader::checkHeader( header => $header );
my %sample_to_col = ();
%sample_to_col = VcfReader::getSamples
(
    get_columns => 1,
    header      => $header,
);
my %csq_header = getAndCheckCsqHeader();
$time = strftime( "%H:%M:%S", localtime );
print STDERR "[$time] INFO - Finished initializing input VCF\n";

my %counts = ();
my @columns = 
qw / 
    Sample
    Total
    Pass
    Filtered
    dbSNP
    Novel
    Het
    Hom
    Ref
/;
if (%csq_header){
    push @columns,
    qw/
        HIGH
        MODERATE
        LOW
        MODIFIER
        HOM_HIGH
    /; 
}
    
print join("\t", @columns) . "\n";
my $n = 0;
countCalls($first_var); 
while (my $line = <$VCF>){
    countCalls($line);
}
close $VCF;

outputCounts();

#################################################
sub outputCounts{
    foreach my $sample (sort keys %counts){
        print join("\t", $sample, 
            map { $counts{$sample}->{$_} ? 
                  $counts{$sample}->{$_} : 
                  0 
            } @columns[1..$#columns] ) .  "\n";
    }
    informUser("Done processing $n autosomal variants\n");
}

#################################################
sub countCalls{
    my $var = shift; 
    chomp $var;
    my @split = split("\t", $var);
    if ($split[0] =~ /^(chr)*[XYM]/){
        return; #only do autosomes
    }
    $n++;
    my %var_info = map { $_ => 0 } qw/ Filtered Pass dbSNP Novel /;
    
    my ($filter, $id) = VcfReader::getMultipleVariantFields(\@split, 'FILTER', 'ID');
    if ($filter eq 'PASS'){
        $var_info{Pass} = 1;
    }else{
        $var_info{Filtered} = 1;
    }
    if ($id =~ /^rs\d/){
        $var_info{dbSNP} = 1;
    }else{
        $var_info{Novel} = 1;
    }
    my @alleles = VcfReader::readAlleles(line => \@split);
    my %allele_impacts = ();
    if (%csq_header){
        my @csq = VcfReader::getVepFields
        ( 
            line        => \@split,
            field       => [ qw / allele impact canonical / ],
            vep_header  => \%csq_header,
        );
        my @alts_to_vep_allele = VcfReader::altsToVepAllele
        (
            ref => $alleles[0],
            alts => [ @alleles[1..$#alleles] ], 
        );
        my %vep_to_allele = map { $alts_to_vep_allele[$_] => $_ + 1 } 0..$#alts_to_vep_allele;
        my @canon = grep { $_->{canonical} eq 'YES' } @csq;
        foreach my $c (@canon){
            my $alt = $vep_to_allele{$c->{allele}};
            if (moreDamaging( $c->{impact}, $allele_impacts{$alt}) ) {
                $allele_impacts{$alt} = $c->{impact};
            }
        }
        foreach (my $i = 1; $i < @alleles; $i++){
            $allele_impacts{$i} ||= "MODIFIER";#e.g. intergenic var with no canonical transcript
        }
    }
    my %samp_to_gt = VcfReader::getSampleCall
    (
        line              => \@split, 
        sample_to_columns => \%sample_to_col,
        minGQ             => $minGq,
        all               => 1,
    ); 
SAMP: foreach my $s (keys %samp_to_gt){#we only support non-phased GTs (i.e. separated by '/').
        next SAMP if $samp_to_gt{$s} eq './.';
        if ($samp_to_gt{$s} eq '0/0'){
            $counts{$s}->{Ref}++;
            next SAMP;
        }
        my @gt = split(/\//, $samp_to_gt{$s});
        if ($gt[0] == $gt[1]){
            next SAMP if $alleles[$gt[0]] eq '*';
            $counts{$s}->{Hom}++;
            if ($allele_impacts{$gt[0]} eq "HIGH"){
                $counts{$s}->{HOM_HIGH}++;
            }
        }else{
            if ( grep { $_ ne 0 and $alleles[$_] ne '*' } @gt){
                $counts{$s}->{Het}++;
            }else{
                next SAMP; # skip REF/* calls
            }
        }
        if (%impacts){
            my $imp = 'MODIFIER';
            foreach my $g (@gt){
                next if $g == 0;
                if (moreDamaging( $allele_impacts{$g}, $imp) ){
                    $imp = $allele_impacts{$g};
                }
            }
            $counts{$s}->{$imp}++;
        }
        $counts{$s}->{Total}++;
        map { $counts{$s}->{$_} += $var_info{$_} } keys %var_info;
    }
}

#################################################
sub moreDamaging{
    my ($new, $old) = @_;
    if (not $old){
        return 1;
    }
    return $impacts{$new} < $impacts{$old};
}

#################################################
sub getAndCheckCsqHeader{
    my %csq_head = ();
    eval { 
        %csq_head = VcfReader::readVepHeader
        (
            header => $header
        ); 
    } ;
    if (not $@){
        informUser("Found VEP header - using VEP 'CSQ' annotations.\n");
    }else{
        informUser
        (
            "Could not find VEP headers in input. ".
            "Will not produce stats for variant consequences.\n"
        );
    }
    return %csq_head;
}

#################################################
sub informUser{
    my $msg = shift;
    my $time = strftime( "%H:%M:%S", localtime );
    print STDERR "[INFO - $time] $msg";
}
