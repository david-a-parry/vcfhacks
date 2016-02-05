use strict;
use warnings;
use Test::More tests => 11;
use FindBin qw($RealBin $Script);
use lib "$RealBin/../lib";
use lib "$RealBin/../lib/Bioperl";
use lib "$RealBin/../lib/BioASN1EntrezGene/lib";
BEGIN 
{ 
    use_ok("VcfhacksUtils");
}
my $n_tests = 1;

my @ar = qw ( a b c c c d c e f e g h ) ;
my @br = qw ( a b c d e f g h );
@ar = VcfhacksUtils::removeDups(@ar); 
is_deeply(
    \@ar,
    \@br,
    "remove duplicate array entries"
);
$n_tests++; 

ok(
    VcfhacksUtils::getAndCheckInSilicoPred('vep', ['all']),
    "getAndCheckInSilicoPred VEP",
);
$n_tests++; 

ok(
    VcfhacksUtils::getAndCheckInSilicoPred('snpeff', ['all']),
    "getAndCheckInSilicoPred SnpEff",
);
$n_tests++; 

ok(
    VcfhacksUtils::readBiotypesFile(),
    "read Biotypes data",
);
$n_tests++; 

ok(
    VcfhacksUtils::readVepClassesFile(),
    "read VEP classes"
);
$n_tests++; 

ok(
    VcfhacksUtils::readSnpEffClassesFile(),
    "read SnpEff classes"
);
$n_tests++; 

my %filter_field = 
(
    ID          => "AFilterField",
    Description => "A made up VCF FILTER field.",
);
my $f_string = '##FILTER=<ID=AFilterField,Description="A made up VCF FILTER field.">';
is(
    VcfhacksUtils::getFilterHeader(%filter_field),
    $f_string,
    "create a FILTER header line"
);
$n_tests++; 

my %format_field = 
(
 ID          => "AFormatField",
 Number      => "A",
 Type        => "String",
 Description => "A made up VCF INFO field.",
);
$f_string = '##FORMAT=<ID=AFormatField,Number=A,Type=String,Description="A made up VCF INFO field.">';
is(
    VcfhacksUtils::getFormatHeader(%format_field),
    $f_string,
    "create a FORMAT header line"
);
$n_tests++; 

my %info_field = 
(
 ID          => "AnInfoField",
 Number      => "A",
 Type        => "String",
 Description => "A made up VCF INFO field.",
);
my $inf_string = '##INFO=<ID=AnInfoField,Number=A,Type=String,Description="A made up VCF INFO field.">';
is(
    VcfhacksUtils::getInfoHeader(%info_field),
    $inf_string,
    "create an INFO header line"
);
$n_tests++; 

my %opts = ( a => "option", b => 1, x => '', y => ["some", "more", "options"]);
my $optstring = "##$Script" .'"a=option b=1 x= y=some,more,options"';
is(
    VcfhacksUtils::getOptsVcfHeader(%opts),
    $optstring,
    "create an Options header line"
);
$n_tests++; 

done_testing($n_tests);
