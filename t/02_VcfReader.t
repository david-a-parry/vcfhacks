use strict;
use warnings;
use Test::More tests => 32;
use FindBin qw($RealBin);
use lib "$RealBin/../lib";
use lib "$RealBin/../lib/dapPerlGenomicLib";
use lib "$RealBin/../lib/Bioperl";
use lib "$RealBin/../lib/BioASN1EntrezGene/lib";
BEGIN 
{ 
    use_ok("VcfReader");
    use_ok("Bio::DB::HTS::Tabix");
}
my $n_tests = 2;

#Set up our input files
my $vcf = "$RealBin/test_data/test1.vcf";
my $index = "$vcf.vridx";
my $shuf = "$RealBin/test_data/test1_shuffled.vcf";
my $tmpsort = "$shuf.tmpsort";
my $gz = "$RealBin/test_data/test1.vcf.gz";
my $tbi = "$gz.tbi";
my $dbsnp = "$RealBin/test_data/dbSnpTestRef.vcf.gz";

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
    my $BH = VcfReader::openVcf($gz), 
    "use openVcf to return a filehandle from a compressed VCF"
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



my %db_sargs = VcfReader::getSearchArguments($dbsnp);
is
(
    VcfReader::searchForPosition(
        %db_sargs,
        chrom => '22',
        pos   => 42991511,
    ),
    1,
    "searchForPosition returns single hit for snp with a snp immediately after",
);
$n_tests++; 
   
is
(
    VcfReader::searchForPosition(
        %db_sargs,
        chrom => '22',
        pos   => 42991512,
    ),
    1,
    "searchForPosition returns single hit for snp with a snp immediately before",
);
$n_tests++; 


is
(
    VcfReader::searchByRegion(
        %db_sargs,
        chrom => '22',
        start => 42991511, 
        end   => 42991511,
    ),
    1,
    "searchByRegion returns single hit for snp with a snp immediately after",
);
$n_tests++; 
   
is
(
    VcfReader::searchByRegion(
        %db_sargs,
        chrom => '22',
        start => 42991512, 
        end   => 42991512,
    ),
    1,
    "searchByRegion returns single hit for snp with a snp immediately before",
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

done_testing($n_tests);
foreach my $f ($tmpsort, $index, $tbi){
    if (-e $f){
        unlink($f) 
         or warn "Could not remove test output ($f): $! "; 
    }
}


