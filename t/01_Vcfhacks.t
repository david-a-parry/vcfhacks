use strict;
use warnings;
use Test::More; # tests => 44;
use FindBin qw($RealBin $Script);
use lib "$RealBin/../lib";
use lib "$RealBin/../lib/Bioperl";
use lib "$RealBin/../lib/BioASN1EntrezGene/lib";
BEGIN 
{ 
    #test loading of 9 local modules
    use_ok("VcfReader");
    use_ok("ClinVarReader");
    use_ok("ParsePedfile");
    use_ok("VcfhacksUtils");
    use_ok("TextToExcel");
    use_ok("Tabix");
    use_ok("SortGenomicCoordinates");
    use_ok("IdParser");
    use_ok("Bio::SeqIO::entrezgene");
}
my $n_tests = 9;

#Set up our input files
my $vcf = "$RealBin/test_data/test1.vcf";
my $index = "$vcf.vridx";
my $shuf = "$RealBin/test_data/test1_shuffled.vcf";
my $tmpsort = "$shuf.tmpsort";
my $gz = "$RealBin/test_data/test1.vcf.gz";
my $tbi = "$gz.tbi";

foreach my $f ($tmpsort, $index, $tbi){
    if (-e $f){
        unlink($f) 
        or die "Could not remove pre-existing test output ($f) - tests aborted: $! "; 
    }
}


#start function tests
ok 
(
    my $FH = VcfReader::openVcf($vcf), 
    "use openVcf to return a filehandle"
);
$n_tests++; 

ok 
(
    VcfReader::checkHeader( vcf => $vcf ),
    "check header for $vcf"
);
$n_tests++; 

ok
(
    my %sargs = VcfReader::getSearchArguments($vcf),
    "get search parameters for uncompressed vcf '$vcf'"
);
$n_tests++; 

ok
(
    my @hits = VcfReader::searchByRegion(
        %sargs,
        chrom => '1',
        start => 3544327, 
        end   => 3544327,
    ),
    "search for coordinate on uncompressed vcf '$vcf'"
);
$n_tests++; 

is
(
    @hits, 
    1,
    "get one hit for search 1:3544327-3544327"
); 
$n_tests++; 

ok
(
    my %tsargs = VcfReader::getSearchArguments($gz),
    "get search parameters for compressed vcf '$gz'"
);
$n_tests++; 

ok
(
    my @thits = VcfReader::searchByRegion(
        %tsargs,
        chrom => '1',
        start => 3544327, 
        end   => 3544327,
    ),
    "search for coordinate on compressed vcf '$gz'"
);
$n_tests++; 

is
(
    @thits, 
    1,
    "get one hit for search 1:3544327-3544327"
); 
$n_tests++; 

   
ok
(
    @hits = VcfReader::searchByRegion(
        %sargs,
        chrom => '6',
        start => 24658070, 
        end   => 25581425,
    ),
    "search for region on uncompressed vcf '$vcf'"
);
$n_tests++; 

is
(
    @hits, 
    4,#4
    "get four hits for search 6:24658070-25581425"
);
$n_tests++; 

ok
(
    @thits = VcfReader::searchByRegion(
        %tsargs,
        chrom => '6',
        start => 24658070, 
        end   => 25581425,
    ),
    "search for region on compressed vcf '$gz'"
);
$n_tests++; 

is
(
    @thits, 
    4,
    "get four hits for search 6:24658070-25581425"
);
$n_tests++; 

ok
(
    @hits = VcfReader::searchByRegion(
        %sargs,
        chrom => 'X',
        start => 153247733, 
        end   => 153357614,
    ),
    "search for non-autosomal region on uncompressed vcf '$vcf'"
);
$n_tests++; 

is
(
    @hits, 
    2,
    "get two hits for search X:153247733-153357614"
);
$n_tests++; 

ok
(
    @thits = VcfReader::searchByRegion(
        %tsargs,
        chrom => 'X',
        start => 153247733, 
        end   => 153357614,
    ),
    "search for non-autosomal region on compressed vcf '$gz'"
);
$n_tests++; 

is
(
    @thits, 
    2,
    "get two hits for search X:153247733-153357614"
);
$n_tests++; 

is
(
    VcfReader::getLineCount($vcf),
    2706,
    "get number of lines from vcf",
); 
$n_tests++; 

is
(
    VcfReader::getFileLengthFromIndex($vcf),
    2706,
    "get number of lines from index",
); 
$n_tests++; 

is
(
    VcfReader::countVariants($vcf),
    2638,
    "get number of variants from vcf",
); 
$n_tests++; 

is
(
    VcfReader::getLineCount($gz),
    2706,
    "get number of lines from compressed vcf",
); 
$n_tests++; 

is
(
    VcfReader::countVariants($gz),
    2638,
    "get number of variants from compressed vcf",
); 
$n_tests++; 

is
(
    VcfReader::checkCoordinateSorted($vcf),
    1,
    "correctly assert vcf is sorted"
);
$n_tests++; 

is
(
    VcfReader::checkCoordinateSorted($shuf),
    0,
    "correctly assert vcf is not sorted"
);
$n_tests++; 

ok
(
    VcfReader::sortVcf(vcf => $shuf, output => "$tmpsort"),
    "use sortVcf on $shuf"
);
$n_tests++; 

is
(
    VcfReader::checkCoordinateSorted($tmpsort),
    1,
    "check sortVcf output is sorted"
);
$n_tests++; 

my @ar = qw ( a b c c c d c e f e g h ) ;
my @br = qw ( a b c d e f g h );
@ar = VcfhacksUtils::removeDups(@ar); 
is(
    @ar,
    @br,
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
