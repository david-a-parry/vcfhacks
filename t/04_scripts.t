use strict;
use warnings;
#use Test::More;
use Test::More tests => 14;
use FindBin qw($RealBin);
use lib "$RealBin/../lib";
use lib "$RealBin/../lib/Bioperl";
use lib "$RealBin/../lib/BioASN1EntrezGene/lib";
BEGIN 
{ 
    use_ok("VcfReader");
    use_ok("Tabix");
}
my $n_tests = 2;

#Set up our input files
my $vcf = "$RealBin/test_data/test1.vcf";
my $index = "$vcf.vridx";
my $shuf = "$RealBin/test_data/test1_shuffled.vcf";
my $tmpsort = "$shuf.tmpsort2";
my $gz = "$RealBin/test_data/test1.vcf.gz";
my $tbi = "$gz.tbi";
my $dref = "$RealBin/test_data/dbSnpTestRef.vcf.gz";
my $dtest = "$RealBin/test_data/dbSnpTestCommonAndRare.vcf";
my $biallelic = "$RealBin/test_data/biallelic_test1.vcf";
my $annoremote = "$RealBin/test_data/annoTestRemote.vcf";

my $script_prefix = "perl $RealBin/..";
my $output = `$script_prefix/findBiallelic.pl -s Sample_133 Sample_76 -r Sample_139  -i $vcf | grep -v '#'`;
my $expected = getBiallelicData();
is
(
    $output,
    $expected,
    "find biallelic variants",
);
$n_tests++; 

$expected = getLocationData();

$output = `$script_prefix/getVariantsByLocation.pl -i $vcf -r 6:24658070-25581425  | cut -sf 1-9`;
is(
    $output,
    $expected,
    "retrieve variants with getVariantsByLocation.pl with an uncompressed vcf '$vcf'"
);
$n_tests++; 

$output = `$script_prefix/getVariantsByLocation.pl -i $gz -r 6:24658070-25581425  | cut -sf 1-9`;
is(
    $output,
    $expected,
    "retrieve variants with getVariantsByLocation.pl with an compressed vcf '$gz'"
);
$n_tests++; 

$expected = getLocatioGeneData();
$output = `$script_prefix/getVariantsByLocation.pl -i $gz -g COL3A1 COL13A1  | cut -sf 1-9`;
is(
    $output,
    $expected,
    "retrieve variants with getVariantsByLocation.pl using gene symbols"
);
$n_tests++; 

$expected = 647;
chomp($output = `$script_prefix/getFunctionalVariants.pl -i $vcf | wc -l `);
is(
    $output,
    $expected,
    "get expected number of functional variants with getFunctionalVariants.pl"
);
$n_tests++; 

$expected = 114;
chomp($output = `$script_prefix/filterOnSample.pl -s Sample_100 -i $vcf | wc -l `);
is(
    $output,
    $expected,
    "get expected number of variants with filterOnSample.pl"
);
$n_tests++; 

`$script_prefix/sortVcf.pl -i $shuf -o $tmpsort`;
is
(
    VcfReader::checkCoordinateSorted($tmpsort),
    1,
    "check sortVcf.pl output is sorted"
);
$n_tests++; 

chomp ($output = `$script_prefix/annotateSnps.pl -i $dtest -d $dref -f 1 | grep -v '#' | wc -l `);
is
(
    $output,
    12,
    "annotateSnps.pl filters on frequency "
);
$n_tests++; 


chomp ($output = `$script_prefix/annotateSnps.pl -i $dtest -d $dref -b 129 | grep -v '#' | wc -l `);
is
(
    $output,
    22,
    "annotateSnps.pl filters on build"
);
$n_tests++; 

chomp ($output = `$script_prefix/annotateSnps.pl -i $dtest -d $dref -f 1 -b 129 | grep -v '#' | wc -l `);
is
(
    $output,
    10,
    "annotateSnps.pl filters on frequency and build"
);
$n_tests++; 

$expected = getAnnoData();
$output = `$script_prefix/geneAnnotator.pl -i $biallelic | cut -sf 1-9`;
like
(
    $output,
    $expected,
    "geneAnnotator.pl retrieves information from local database",
);
$n_tests++; 


$expected = getAnnoDataRemote();
$output = `$script_prefix/geneAnnotator.pl -i $annoremote | cut -sf 1-9`;
like
(
    $output,
    $expected,
    "geneAnnotator.pl retrieves information via REST server",
);
$n_tests++; 

done_testing($n_tests);
cleanup();

###############################################################################
sub cleanup{
    foreach my $f ($tmpsort, $index, $tbi){
        if (-e $f){
            unlink($f) 
             or warn "Could not remove test output ($f): $! "; 
        }
    }
}

###############################################################################
sub getAnnoDataRemote{
#exact data retrieved will change on database updates so just use a regex to
#check this is working
    return qr/;GeneAnno=ENSG00000158109\|\d+\|\w+\|\S+\|\S+\|\S+\|\S*\|\S*\|\S*\|\S*\tGT:AD:DP:GQ:PL/;
}

###############################################################################
sub getAnnoData{
#exact data retrieved will change on database updates so just use a regex to
#check this is working
    return qr/;GeneAnno=ENSG00000197467\|\d+\|\w+\|\S+\|\S+\|\S+\|\S*\|\S+\|\S+\tGT:AD:DP:GQ:PL/;
}

###############################################################################
sub getBiallelicData{
    open (my $IN, "<", $biallelic) or die "Can't open findBiallelic.pl test data: $!\n";
    my $data;
    while (<$IN>){
        $data .= $_ unless /^#/;
    }
    close $IN;
    return $data;
}

###############################################################################
sub getLocationData{
    return <<EOT
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
6	24658076	.	G	C	810.90	PASS	AC=2;AF=4.202e-03;AN=476;BaseQRankSum=-2.062e+00;ClippingRankSum=0.380;DP=4919;FS=0.000;GQ_MEAN=58.63;GQ_STDDEV=40.44;InbreedingCoeff=-0.0073;MLEAC=2;MLEAF=4.202e-03;MQ=70.00;MQ0=0;MQRankSum=1.08;NCC=0;QD=13.08;ReadPosRankSum=0.462;SOR=0.708;VQSLOD=10.34;culprit=MQ;set=variant19;SnpAnnotation=[.];CSQ=C|ENSG00000111802|ENST00000545995|Transcript|missense_variant|571|571|191|L/V|Cta/Gta|||-1|TDP2|HGNC|17768|protein_coding||||ENSP00000437637|TYDP2_HUMAN|I3QNV0_HUMAN|UPI00001584B7|tolerated(0.1)|possibly_damaging(0.477)|4/7||Pfam_domain:PF03372&Superfamily_domains:SSF56219|ENST00000545995.1:c.571C>G|ENSP00000437637.1:p.Leu191Val|||||||||||||||||||neutral(0.444),C|ENSG00000111802|ENST00000478507|Transcript|intron_variant&non_coding_transcript_variant||||||||-1|TDP2|HGNC|17768|processed_transcript|||||||||||2/3||ENST00000478507.1:n.320-4695C>G||||||||||||||||||||,C|ENSG00000111802|ENST00000478285|Transcript|non_coding_transcript_exon_variant&non_coding_transcript_variant|668|||||||-1|TDP2|HGNC|17768|processed_transcript||||||||||2/3|||ENST00000478285.1:n.668C>G||||||||||||||||||||,C|ENSG00000111802|ENST00000341060|Transcript|missense_variant|703|307|103|L/V|Cta/Gta|||-1|TDP2|HGNC|17768|protein_coding||||ENSP00000345345|||UPI00001A3A86|tolerated(0.06)|possibly_damaging(0.499)|3/6||Pfam_domain:PF03372&Superfamily_domains:SSF56219|ENST00000341060.3:c.307C>G|ENSP00000345345.3:p.Leu103Val|||||||||||||||||||deleterious(0.489),C|ENSG00000111802|ENST00000378198|Transcript|missense_variant|652|481|161|L/V|Cta/Gta|||-1|TDP2|HGNC|17768|protein_coding|YES||CCDS4557.1|ENSP00000367440|TYDP2_HUMAN||UPI0000032018|deleterious(0.01)|probably_damaging(0.934)|4/7||Pfam_domain:PF03372&Superfamily_domains:SSF56219|ENST00000378198.4:c.481C>G|ENSP00000367440.4:p.Leu161Val|||||||||||||||||||deleterious(0.786)	GT:AD:DP:GQ:PL
6	25109829	.	CCAA	C	958.97	PASS	AC=2;AF=4.237e-03;AN=472;BaseQRankSum=1.95;ClippingRankSum=-4.120e-01;DP=4944;FS=2.777;GQ_MEAN=60.15;GQ_STDDEV=45.84;InbreedingCoeff=-0.0070;MLEAC=2;MLEAF=4.237e-03;MQ=70.00;MQ0=0;MQRankSum=-7.660e-01;NCC=2;QD=19.57;ReadPosRankSum=0.450;SOR=1.347;VQSLOD=3.58;culprit=FS;set=variant19;SnpAnnotation=[.];CSQ=-|ENSG00000168405|ENST00000493257|Transcript|upstream_gene_variant|||||||3152|-1|CMAHP|HGNC|2098|processed_transcript|||||||||||||||||||||||||||||||||,-|ENSG00000168405|ENST00000436589|Transcript|intron_variant&non_coding_transcript_variant||||||||-1|CMAHP|HGNC|2098|unitary_pseudogene|||||||||||4/12||ENST00000436589.2:n.632-23_632-21delTTG||||||||||||||||||||,-|ENSG00000168405|ENST00000377989|Transcript|intron_variant&non_coding_transcript_variant||||||||-1|CMAHP|HGNC|2098|processed_transcript|YES||||||||||5/13||ENST00000377989.4:n.1174-23_1174-21delTTG||||||||||||||||||||,-|ENSG00000168405|ENST00000490939|Transcript|upstream_gene_variant|||||||920|-1|CMAHP|HGNC|2098|retained_intron|||||||||||||||||||||||||||||||||,-|ENSG00000168405|ENST00000377993|Transcript|intron_variant&non_coding_transcript_variant||||||||-1|CMAHP|HGNC|2098|unitary_pseudogene|||||||||||5/13||ENST00000377993.2:n.662-23_662-21delTTG||||||||||||||||||||,-|ENSG00000168405|ENST00000493981|Transcript|upstream_gene_variant|||||||76|-1|CMAHP|HGNC|2098|processed_transcript|||||||||||||||||||||||||||||||||,-|ENSG00000168405|ENST00000458373|Transcript|intron_variant&non_coding_transcript_variant||||||||-1|CMAHP|HGNC|2098|processed_transcript|||||||||||5/5||ENST00000458373.1:n.920-23_920-21delTTG||||||||||||||||||||,-|ENSG00000168405|ENST00000424282|Transcript|intron_variant&non_coding_transcript_variant||||||||-1|CMAHP|HGNC|2098|unitary_pseudogene|||||||||||5/9||ENST00000424282.2:n.930-23_930-21delTTG||||||||||||||||||||	GT:AD:DP:GQ:PL
6	25471559	rs144894905;rs374781233	GT	G,GTT,GTTT	21522.17	PASS	AC=13,143,15;AF=0.028,0.303,0.032;AN=472;BaseQRankSum=-1.980e-01;ClippingRankSum=0.340;DP=7057;FS=2.897;GQ_MEAN=47.68;GQ_STDDEV=62.40;InbreedingCoeff=-0.3509;MLEAC=7,147,12;MLEAF=0.015,0.311,0.025;MQ=70.00;MQ0=0;MQRankSum=0.063;NCC=2;QD=6.63;ReadPosRankSum=-5.360e-01;SOR=0.356;VQSLOD=3.28;culprit=FS;set=variant19;SnpAnnotation=[.,dbSNPBuildID=134/138,.];CSQ=-|ENSG00000079691|ENST00000377969|Transcript|intron_variant||||||||1|LRRC16A|HGNC|21581|protein_coding||||ENSP00000367206|LR16A_HUMAN||UPI0000D61400||||10/11||ENST00000377969.3:c.296+75delT||||||||||||||||||||,TTT|ENSG00000079691|ENST00000377969|Transcript|intron_variant||||||||1|LRRC16A|HGNC|21581|protein_coding||||ENSP00000367206|LR16A_HUMAN||UPI0000D61400||||10/11||ENST00000377969.3:c.296+75[3]T||||||||||||||||||||,TT|ENSG00000079691|ENST00000377969|Transcript|intron_variant||||||||1|LRRC16A|HGNC|21581|protein_coding||||ENSP00000367206|LR16A_HUMAN||UPI0000D61400||||10/11||ENST00000377969.3:c.296+75dupT||||||||||||||||||||,-|ENSG00000079691|ENST00000497227|Transcript|intron_variant&non_coding_transcript_variant||||||||1|LRRC16A|HGNC|21581|processed_transcript|||||||||||2/2||ENST00000497227.1:n.502+75delT||||||||||||||||||||,TTT|ENSG00000079691|ENST00000497227|Transcript|intron_variant&non_coding_transcript_variant||||||||1|LRRC16A|HGNC|21581|processed_transcript|||||||||||2/2||ENST00000497227.1:n.502+75[3]T||||||||||||||||||||,TT|ENSG00000079691|ENST00000497227|Transcript|intron_variant&non_coding_transcript_variant||||||||1|LRRC16A|HGNC|21581|processed_transcript|||||||||||2/2||ENST00000497227.1:n.502+75dupT||||||||||||||||||||,-|ENSG00000079691|ENST00000329474|Transcript|intron_variant||||||||1|LRRC16A|HGNC|21581|protein_coding|YES||CCDS54973.1|ENSP00000331983|LR16A_HUMAN||UPI00004588AB||||10/36||ENST00000329474.6:c.779+75delT||||||||||||||||||||,TTT|ENSG00000079691|ENST00000329474|Transcript|intron_variant||||||||1|LRRC16A|HGNC|21581|protein_coding|YES||CCDS54973.1|ENSP00000331983|LR16A_HUMAN||UPI00004588AB||||10/36||ENST00000329474.6:c.779+75[3]T||||||||||||||||||||,TT|ENSG00000079691|ENST00000329474|Transcript|intron_variant||||||||1|LRRC16A|HGNC|21581|protein_coding|YES||CCDS54973.1|ENSP00000331983|LR16A_HUMAN||UPI00004588AB||||10/36||ENST00000329474.6:c.779+75dupT||||||||||||||||||||	GT:AD:DP:GQ:PGT:PID:PL
6	25581415	.	T	C	1832.82	PASS	AC=2;AF=4.202e-03;AN=476;BaseQRankSum=4.33;ClippingRankSum=-2.120e-01;DP=6052;FS=0.000;GQ_MEAN=73.79;GQ_STDDEV=71.95;InbreedingCoeff=-0.0046;MLEAC=2;MLEAF=4.202e-03;MQ=70.00;MQ0=0;MQRankSum=0.804;NCC=0;QD=18.33;ReadPosRankSum=1.69;SOR=0.607;VQSLOD=9.95;culprit=MQ;set=variant19;SnpAnnotation=[.];CSQ=C|ENSG00000079691|ENST00000329474|Transcript|intron_variant||||||||1|LRRC16A|HGNC|21581|protein_coding|YES||CCDS54973.1|ENSP00000331983|LR16A_HUMAN||UPI00004588AB||||30/36||ENST00000329474.6:c.2810-56T>C||||||||||||||||||||	GT:AD:DP:GQ:PL
EOT
;
}

###############################################################################
sub getLocatioGeneData{
    return <<EOT
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
2	189868932	.	A	G	1643.81	PASS	AC=2;AF=4.202e-03;AN=476;BaseQRankSum=2.94;ClippingRankSum=0.817;DP=6280;FS=0.000;GQ_MEAN=74.84;GQ_STDDEV=68.40;InbreedingCoeff=-0.0042;MLEAC=2;MLEAF=4.202e-03;MQ=70.00;MQ0=0;MQRankSum=-1.400e-02;NCC=0;QD=17.12;ReadPosRankSum=-4.660e-01;SOR=0.749;VQSLOD=10.90;culprit=MQ;set=variant15;SnpAnnotation=[.];CSQ=G|ENSG00000168542|ENST00000317840|Transcript|intron_variant||||||||1|COL3A1|HGNC|2201|protein_coding||||ENSP00000315243||Q6LBY7_HUMAN&E7ENY8_HUMAN|UPI0000D61281||||36/41||ENST00000317840.5:c.2520+1177A>G||||||||||||||||||||,G|ENSG00000168542|ENST00000304636|Transcript|intron_variant||||||||1|COL3A1|HGNC|2201|protein_coding|YES||CCDS2297.1|ENSP00000304408|CO3A1_HUMAN|Q6LBY7_HUMAN&D2JYH5_HUMAN|UPI0000456EBA||||39/50||ENST00000304636.3:c.2824-51A>G||||||||||||||||||||,G|ENSG00000168542|ENST00000467886|Transcript|downstream_gene_variant|||||||289|1|COL3A1|HGNC|2201|retained_intron|||||||||||||||||||||||||||||||||,G|ENSG00000168542|ENST00000487010|Transcript|upstream_gene_variant|||||||4096|1|COL3A1|HGNC|2201|retained_intron|||||||||||||||||||||||||||||||||	GT:AD:DP:GQ:PL
10	71648058	.	AG	A	1415.82	LOW_VQSLOD	AC=2;AF=4.219e-03;AN=474;DP=6179;FS=0.000;GQ_MEAN=68.46;GQ_STDDEV=17.15;InbreedingCoeff=0.9954;MLEAC=2;MLEAF=4.219e-03;MQ=60.64;MQ0=0;NCC=0;NEGATIVE_TRAIN_SITE;QD=33.71;SOR=0.788;VQSLOD=-9.568e-01;culprit=MQ;set=variant;SnpAnnotation=[.];CSQ=-|ENSG00000197467|ENST00000517713|Transcript|splice_acceptor_variant||||||||1|COL13A1|HGNC|2190|protein_coding|||CCDS44423.2|ENSP00000430061|CODA1_HUMAN|Q9UP45_HUMAN|UPI0000F6E6D3||||7/36||ENST00000517713.1:c.523-1delG|||||||||||||||||||HC|,-|ENSG00000197467|ENST00000354547|Transcript|splice_acceptor_variant||||||||1|COL13A1|HGNC|2190|protein_coding|||CCDS44424.2|ENSP00000346553|CODA1_HUMAN|Q9UP45_HUMAN|UPI0000EE047D||||7/38||ENST00000354547.3:c.523-1delG|||||||||||||||||||HC|,-|ENSG00000197467|ENST00000398978|Transcript|splice_acceptor_variant||||||||1|COL13A1|HGNC|2190|protein_coding|YES||CCDS44419.1|ENSP00000381949|CODA1_HUMAN|Q9UP45_HUMAN|UPI000046FD72||||7/39||ENST00000398978.3:c.523-1delG|||||||||||||||||||HC|,-|ENSG00000197467|ENST00000398964|Transcript|splice_acceptor_variant||||||||1|COL13A1|HGNC|2190|protein_coding||||ENSP00000381936||Q9UP45_HUMAN&E7ES56_HUMAN|UPI0000160933||||4/36||ENST00000398964.3:c.436-1delG|||||||||||||||||||HC|,-|ENSG00000197467|ENST00000398971|Transcript|splice_acceptor_variant||||||||1|COL13A1|HGNC|2190|protein_coding||||ENSP00000381943||Q9UP45_HUMAN&E7ES49_HUMAN|UPI0000160930||||6/37||ENST00000398971.3:c.523-1delG|||||||||||||||||||HC|,-|ENSG00000197467|ENST00000398973|Transcript|splice_acceptor_variant||||||||1|COL13A1|HGNC|2190|protein_coding||||ENSP00000381945||Q9UP45_HUMAN&E9PEG9_HUMAN|UPI0002065296||||7/37||ENST00000398973.3:c.523-1delG|||||||||||||||||||HC|,-|ENSG00000197467|ENST00000398966|Transcript|splice_acceptor_variant||||||||1|COL13A1|HGNC|2190|protein_coding||||ENSP00000381938||Q9UP45_HUMAN&E7ES55_HUMAN|UPI000016092F||||6/37||ENST00000398966.3:c.523-1delG|||||||||||||||||||HC|,-|ENSG00000197467|ENST00000398969|Transcript|splice_acceptor_variant||||||||1|COL13A1|HGNC|2190|protein_coding||||ENSP00000381941||Q9UP45_HUMAN&E7ES50_HUMAN|UPI0002065294||||4/34||ENST00000398969.3:c.409-1delG|||||||||||||||||||HC|,-|ENSG00000197467|ENST00000520133|Transcript|splice_acceptor_variant||||||||1|COL13A1|HGNC|2190|protein_coding|||CCDS44428.2|ENSP00000430173|CODA1_HUMAN|Q9UP45_HUMAN|UPI0001CA7E68||||5/33||ENST00000520133.1:c.436-1delG|||||||||||||||||||HC|,-|ENSG00000197467|ENST00000520267|Transcript|splice_acceptor_variant||||||||1|COL13A1|HGNC|2190|protein_coding|||CCDS44427.2|ENSP00000428057|CODA1_HUMAN|Q9UP45_HUMAN|UPI000192C3E9||||4/34||ENST00000520267.1:c.409-1delG|||||||||||||||||||HC|,-|ENSG00000197467|ENST00000479733|Transcript|splice_acceptor_variant&NMD_transcript_variant||||||||1|COL13A1|HGNC|2190|nonsense_mediated_decay||||ENSP00000430089||Q9UP45_HUMAN&E7EX21_HUMAN|UPI0001E8F0C8||||8/40||ENST00000479733.1:c.550-1delG||||||||||||||||||||,-|ENSG00000197467|ENST00000522165|Transcript|splice_acceptor_variant||||||||1|COL13A1|HGNC|2190|protein_coding|||CCDS44425.2|ENSP00000428342|CODA1_HUMAN|Q9UP45_HUMAN|UPI0000F6E6D7||||7/37||ENST00000522165.1:c.523-1delG|||||||||||||||||||HC|,-|ENSG00000197467|ENST00000398968|Transcript|splice_acceptor_variant||||||||1|COL13A1|HGNC|2190|protein_coding||||ENSP00000381940||Q9UP45_HUMAN&E7ES51_HUMAN|UPI000016092D||||6/37||ENST00000398968.3:c.523-1delG|||||||||||||||||||HC|,-|ENSG00000197467|ENST00000398974|Transcript|splice_acceptor_variant||||||||1|COL13A1|HGNC|2190|protein_coding||||ENSP00000381946||Q9UP45_HUMAN&E7ES46_HUMAN|UPI0000160926||||5/37||ENST00000398974.3:c.487-1delG|||||||||||||||||||HC|,-|ENSG00000197467|ENST00000356340|Transcript|splice_acceptor_variant||||||||1|COL13A1|HGNC|2190|protein_coding||||ENSP00000348695||Q9UP45_HUMAN&G5E987_HUMAN|UPI0000167B85||||6/38||ENST00000356340.3:c.523-1delG|||||||||||||||||||HC|,-|ENSG00000197467|ENST00000398972|Transcript|splice_acceptor_variant||||||||1|COL13A1|HGNC|2190|protein_coding||||ENSP00000381944||Q9UP45_HUMAN&E7ES47_HUMAN|UPI0002065295||||7/38||ENST00000398972.3:c.523-1delG|||||||||||||||||||HC|,-|ENSG00000197467|ENST00000357811|Transcript|splice_acceptor_variant||||||||1|COL13A1|HGNC|2190|protein_coding||||ENSP00000350463|CODA1_HUMAN|Q9UP45_HUMAN|UPI000046FD74||||7/37||ENST00000357811.3:c.523-1delG|||||||||||||||||||HC|	GT:AD:DP:GQ:PL
10	71682523	.	AG	A	1150.87	PASS	AC=3;AF=6.329e-03;AN=474;BaseQRankSum=-1.170e-01;ClippingRankSum=1.29;DP=6195;FS=1.808;GQ_MEAN=65.26;GQ_STDDEV=21.48;InbreedingCoeff=0.6411;MLEAC=3;MLEAF=6.329e-03;MQ=70.00;MQ0=0;MQRankSum=-4.280e-01;NCC=0;QD=22.13;ReadPosRankSum=0.371;SOR=1.162;VQSLOD=3.34;culprit=SOR;set=variant;SnpAnnotation=[.];CSQ=-|ENSG00000197467|ENST00000517713|Transcript|frameshift_variant|1105|1105|369|G/X|Ggg/gg|||1|COL13A1|HGNC|2190|protein_coding|||CCDS44423.2|ENSP00000430061|CODA1_HUMAN|Q9UP45_HUMAN|UPI0000F6E6D3|||21/37||Low_complexity_(Seg):Seg&Pfam_domain:PF01391&PROSITE_profiles:PS50315&PROSITE_profiles:PS50099|ENST00000517713.1:c.1105delG|ENSP00000430061.1:p.Leu370SerfsTer71||||||||||||||||||HC|,-|ENSG00000197467|ENST00000354547|Transcript|frameshift_variant|1597|1105|369|G/X|Ggg/gg|||1|COL13A1|HGNC|2190|protein_coding|||CCDS44424.2|ENSP00000346553|CODA1_HUMAN|Q9UP45_HUMAN|UPI0000EE047D|||21/39||Low_complexity_(Seg):Seg&Pfam_domain:PF01391&PROSITE_profiles:PS50315&PROSITE_profiles:PS50099|ENST00000354547.3:c.1105delG|ENSP00000346553.3:p.Leu370SerfsTer71||||||||||||||||||HC|,-|ENSG00000197467|ENST00000398975|Transcript|upstream_gene_variant|||||||4320|1|COL13A1|HGNC|2190|protein_coding||||ENSP00000381947|||UPI000059D163||||||||||||||||||||||||||,-|ENSG00000197467|ENST00000398978|Transcript|frameshift_variant|1663|1171|391|G/X|Ggg/gg|||1|COL13A1|HGNC|2190|protein_coding|YES||CCDS44419.1|ENSP00000381949|CODA1_HUMAN|Q9UP45_HUMAN|UPI000046FD72|||22/40||Low_complexity_(Seg):Seg&Pfam_domain:PF01391&PROSITE_profiles:PS50315&PROSITE_profiles:PS50099|ENST00000398978.3:c.1171delG|ENSP00000381949.3:p.Leu392SerfsTer71||||||||||||||||||HC|,-|ENSG00000197467|ENST00000398964|Transcript|frameshift_variant|1620|1084|362|G/X|Ggg/gg|||1|COL13A1|HGNC|2190|protein_coding||||ENSP00000381936||Q9UP45_HUMAN&E7ES56_HUMAN|UPI0000160933|||19/37||Low_complexity_(Seg):Seg&Pfam_domain:PF01391&PROSITE_profiles:PS50315&PROSITE_profiles:PS50099|ENST00000398964.3:c.1084delG|ENSP00000381936.3:p.Leu363SerfsTer71||||||||||||||||||HC|,-|ENSG00000197467|ENST00000398971|Transcript|frameshift_variant|1707|1171|391|G/X|Ggg/gg|||1|COL13A1|HGNC|2190|protein_coding||||ENSP00000381943||Q9UP45_HUMAN&E7ES49_HUMAN|UPI0000160930|||21/38||Low_complexity_(Seg):Seg&Pfam_domain:PF01391&PROSITE_profiles:PS50315&PROSITE_profiles:PS50099|ENST00000398971.3:c.1171delG|ENSP00000381943.3:p.Leu392SerfsTer71||||||||||||||||||HC|,-|ENSG00000197467|ENST00000398973|Transcript|frameshift_variant|1707|1171|391|G/X|Ggg/gg|||1|COL13A1|HGNC|2190|protein_coding||||ENSP00000381945||Q9UP45_HUMAN&E9PEG9_HUMAN|UPI0002065296|||22/38||Low_complexity_(Seg):Seg&Pfam_domain:PF01391&PROSITE_profiles:PS50315&PROSITE_profiles:PS50099|ENST00000398973.3:c.1171delG|ENSP00000381945.3:p.Leu392SerfsTer134||||||||||||||||||HC|,-|ENSG00000197467|ENST00000398966|Transcript|frameshift_variant|1641|1105|369|G/X|Ggg/gg|||1|COL13A1|HGNC|2190|protein_coding||||ENSP00000381938||Q9UP45_HUMAN&E7ES55_HUMAN|UPI000016092F|||20/38||Low_complexity_(Seg):Seg&Pfam_domain:PF01391&PROSITE_profiles:PS50315&PROSITE_profiles:PS50099|ENST00000398966.3:c.1105delG|ENSP00000381938.3:p.Leu370SerfsTer71||||||||||||||||||HC|,-|ENSG00000197467|ENST00000398969|Transcript|frameshift_variant|1536|1000|334|G/X|Ggg/gg|||1|COL13A1|HGNC|2190|protein_coding||||ENSP00000381941||Q9UP45_HUMAN&E7ES50_HUMAN|UPI0002065294|||18/35||Low_complexity_(Seg):Seg&Pfam_domain:PF01391&PROSITE_profiles:PS50315&PROSITE_profiles:PS50099|ENST00000398969.3:c.1000delG|ENSP00000381941.3:p.Leu335SerfsTer134||||||||||||||||||HC|,-|ENSG00000197467|ENST00000520133|Transcript|frameshift_variant|1018|1018|340|G/X|Ggg/gg|||1|COL13A1|HGNC|2190|protein_coding|||CCDS44428.2|ENSP00000430173|CODA1_HUMAN|Q9UP45_HUMAN|UPI0001CA7E68|||19/34||Low_complexity_(Seg):Seg&Pfam_domain:PF01391&PROSITE_profiles:PS50315&PROSITE_profiles:PS50099|ENST00000520133.1:c.1018delG|ENSP00000430173.1:p.Leu341SerfsTer71||||||||||||||||||HC|,-|ENSG00000197467|ENST00000520267|Transcript|frameshift_variant|1306|1000|334|G/X|Ggg/gg|||1|COL13A1|HGNC|2190|protein_coding|||CCDS44427.2|ENSP00000428057|CODA1_HUMAN|Q9UP45_HUMAN|UPI000192C3E9|||18/35||Low_complexity_(Seg):Seg&Pfam_domain:PF01391&PROSITE_profiles:PS50315&PROSITE_profiles:PS50099|ENST00000520267.1:c.1000delG|ENSP00000428057.1:p.Leu335SerfsTer71||||||||||||||||||HC|,-|ENSG00000197467|ENST00000479733|Transcript|frameshift_variant&NMD_transcript_variant|1691|1198|400|G/X|Ggg/gg|||1|COL13A1|HGNC|2190|nonsense_mediated_decay||||ENSP00000430089||Q9UP45_HUMAN&E7EX21_HUMAN|UPI0001E8F0C8|||23/41||Low_complexity_(Seg):Seg&Pfam_domain:PF01391&PROSITE_profiles:PS50315&PROSITE_profiles:PS50099|ENST00000479733.1:c.1198delG|ENSP00000430089.1:p.Leu401SerfsTer71|||||||||||||||||||,-|ENSG00000197467|ENST00000522165|Transcript|frameshift_variant|1114|1114|372|G/X|Ggg/gg|||1|COL13A1|HGNC|2190|protein_coding|||CCDS44425.2|ENSP00000428342|CODA1_HUMAN|Q9UP45_HUMAN|UPI0000F6E6D7|||21/38||Low_complexity_(Seg):Seg&Pfam_domain:PF01391&PROSITE_profiles:PS50315&PROSITE_profiles:PS50099|ENST00000522165.1:c.1114delG|ENSP00000428342.1:p.Leu373SerfsTer71||||||||||||||||||HC|,-|ENSG00000197467|ENST00000398968|Transcript|frameshift_variant|1650|1114|372|G/X|Ggg/gg|||1|COL13A1|HGNC|2190|protein_coding||||ENSP00000381940||Q9UP45_HUMAN&E7ES51_HUMAN|UPI000016092D|||20/38||Low_complexity_(Seg):Seg&Pfam_domain:PF01391&PROSITE_profiles:PS50315&PROSITE_profiles:PS50099|ENST00000398968.3:c.1114delG|ENSP00000381940.3:p.Leu373SerfsTer71||||||||||||||||||HC|,-|ENSG00000197467|ENST00000398974|Transcript|frameshift_variant|1671|1135|379|G/X|Ggg/gg|||1|COL13A1|HGNC|2190|protein_coding||||ENSP00000381946||Q9UP45_HUMAN&E7ES46_HUMAN|UPI0000160926|||20/38||Low_complexity_(Seg):Seg&Pfam_domain:PF01391&PROSITE_profiles:PS50315&PROSITE_profiles:PS50099|ENST00000398974.3:c.1135delG|ENSP00000381946.3:p.Leu380SerfsTer71||||||||||||||||||HC|,-|ENSG00000197467|ENST00000356340|Transcript|frameshift_variant|1707|1171|391|G/X|Ggg/gg|||1|COL13A1|HGNC|2190|protein_coding||||ENSP00000348695||Q9UP45_HUMAN&G5E987_HUMAN|UPI0000167B85|||21/39||Low_complexity_(Seg):Seg&Pfam_domain:PF01391&PROSITE_profiles:PS50315&PROSITE_profiles:PS50099|ENST00000356340.3:c.1171delG|ENSP00000348695.3:p.Leu392SerfsTer71||||||||||||||||||HC|,-|ENSG00000197467|ENST00000398972|Transcript|frameshift_variant|1707|1171|391|G/X|Ggg/gg|||1|COL13A1|HGNC|2190|protein_coding||||ENSP00000381944||Q9UP45_HUMAN&E7ES47_HUMAN|UPI0002065295|||22/39||Low_complexity_(Seg):Seg&Pfam_domain:PF01391&PROSITE_profiles:PS50315&PROSITE_profiles:PS50099|ENST00000398972.3:c.1171delG|ENSP00000381944.3:p.Leu392SerfsTer134||||||||||||||||||HC|,-|ENSG00000197467|ENST00000357811|Transcript|frameshift_variant|1597|1105|369|G/X|Ggg/gg|||1|COL13A1|HGNC|2190|protein_coding||||ENSP00000350463|CODA1_HUMAN|Q9UP45_HUMAN|UPI000046FD74|||21/38||Low_complexity_(Seg):Seg&Pfam_domain:PF01391&PROSITE_profiles:PS50315&PROSITE_profiles:PS50099|ENST00000357811.3:c.1105delG|ENSP00000350463.3:p.Leu370SerfsTer71||||||||||||||||||HC|	GT:AD:DP:GQ:PL
EOT
;
}
