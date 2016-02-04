use strict;
use warnings;
use Test::More; #tests => 10;
use FindBin qw($RealBin);
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

ok 
(
    VcfReader::checkHeader( vcf => $vcf ),
    "check header for $vcf"
);

ok
(
    my %sargs = VcfReader::getSearchArguments($vcf),
    "get search parameters for uncompressed vcf '$vcf'"
);

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

is
(
    @hits, 
    1,
    "get one hit for search 1:3544327-3544327"
); 

ok
(
    my %tsargs = VcfReader::getSearchArguments($gz),
    "get search parameters for compressed vcf '$gz'"
);

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

is
(
    @thits, 
    1,
    "get one hit for search 1:3544327-3544327"
); 

   
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

is
(
    @hits, 
    4,#4
    "get four hits for search 6:24658070-25581425"
);

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

is
(
    @thits, 
    4,
    "get four hits for search 6:24658070-25581425"
);

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

is
(
    @hits, 
    2,
    "get two hits for search X:153247733-153357614"
);

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

is
(
    @thits, 
    2,
    "get two hits for search X:153247733-153357614"
);

is
(
    VcfReader::getLineCount($vcf),
    2706,
    "get number of lines from vcf",
); 

is
(
    VcfReader::getFileLengthFromIndex($vcf),
    2706,
    "get number of lines from index",
); 

is
(
    VcfReader::countVariants($vcf),
    2638,
    "get number of variants from vcf",
); 

is
(
    VcfReader::getLineCount($gz),
    2706,
    "get number of lines from compressed vcf",
); 

is
(
    VcfReader::countVariants($gz),
    2638,
    "get number of variants from compressed vcf",
); 

is
(
    VcfReader::checkCoordinateSorted($vcf),
    1,
    "correctly assert vcf is sorted"
);

is
(
    VcfReader::checkCoordinateSorted($shuf),
    0,
    "correctly assert vcf is not sorted"
);

ok
(
    VcfReader::sortVcf(vcf => $shuf, output => "$tmpsort"),
    "use sortVcf on $shuf"
);

is
(
    VcfReader::checkCoordinateSorted($tmpsort),
    1,
    "check sortVcf output is sorted"
);



done_testing();
