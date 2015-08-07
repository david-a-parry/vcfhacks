#!/usr/bin/env perl 
use strict;
use warnings;
use Getopt::Long;
use Bio::DB::Sam;
use Data::Dumper;
use File::Temp qw/ tempfile /;
use FindBin;
use lib "$FindBin::Bin/lib";
use VcfReader;

my %opts = ();
GetOptions(
    \%opts,
    'i|input=s',        #hgmd input
    'o|output=s',       #optional output file
    'f|fasta=s',        #genome fasta
    'u|unconverted=s',  #optional file for unconverted variants
    's|sort',           #flag to sort final output in coordinate order
    'h|?|help',
) or pod2usage(-exitval => 2, -message => "Syntax error.\n"); 

usage() if $opts{h};
usage("-i/--input is required" ) if (not $opts{i});
usage("-f/--fasta is required" ) if (not $opts{f});

#OPEN AND CHECK HGMD INPUT
my $FH;
if ($opts{i} =~ /\.gz$/){
    open ($FH, "gzip -dc $opts{i} |") or die "Can't open $opts{i} via gzip: $!\n";
}else{
    open ($FH, $opts{i} ) or die "Can't open $opts{i} for reading: $!\n";
}
my %hgmd_columns = ();
my $header = <$FH>;
chomp (my @head = split("\t", $header)); 
my $n = 0;
%hgmd_columns = map { lc($_) => $n++ } map { (my $trim = $_) =~ s/^#+//; $trim } @head;
my @hgmd_fields = 
    (
    "hgmd id",
    "disease",
    "variant class",
    "gene symbol",
    "chromosome",
    "coordinate start",
    "coordinate end",
    "strand",
    "hgvs",
);

foreach my $req (@hgmd_fields){
    if (not exists $hgmd_columns{$req}){
        die "Required column ($req) not found in HGMD file $opts{i}\nFound the following columns:\n" 
            .join("\n", @head) . "\n";
    }
}

#OPEN AND CHECK FASTA REFERENCE SEQUENCE

my $seqs_have_chr = 0;#check for fasta naming convention
my $seqs_are_mixed = 0;#check for fasta naming convention
my $fai = Bio::DB::Sam::Fai->load($opts{f});#should create index if it doesn't exist
#get hash of chromosomes to lengths from index and check naming convention
my %chroms = readFastaIndex($opts{f}); 

#SET UP OUTPUT

my $OUT = \*STDOUT;#main output defaults to STDOUT
my $UNCONV;#optional filehandle for writing unconverted variants
my $tmp; #filename of tempfile if we're sorting output
initializeOutput();

#START READING AND CONVERTING HGMD RECORDS

while (my $line = <$FH>){
    chomp($line);
    my @fields = split("\t", $line); 
    if (my $vcf_record = convertToVcf(\@fields) ){
        print $OUT "$vcf_record\n";
    }elsif($UNCONV){
        print $UNCONV "$line\n";
    }
}
close $FH;
close $OUT;
if ($opts{s}){
   sortOutput();   
}

###########################################################
sub sortOutput{
    my %args = (vcf => $tmp); 
    if ($opts{o}){
        $args{output} = $opts{o};
    }
    VcfReader::sortVcf( %args );
}
###########################################################
sub initializeOutput{
    if ($opts{s}){
        ($OUT, $tmp) = tempfile();
    }elsif ($opts{o}){
        open ($OUT, ">", $opts{o}) or die "Can't open $opts{o} for writing: $!\n";
    }
    print $OUT "##fileformat=VCFv4.2\n";
    foreach my $f (@hgmd_fields){
        (my $tag = $f) =~ s/\s/_/g; 
        print $OUT "##INFO=<ID=$tag,Number=1,Type=String,Description=\"$f field from HGMD mart download\",Source=\"$opts{i}\">\n"
    }
    print $OUT '#' . join("\t",  ( qw / CHROM POS ID REF ALT QUAL FILTER INFO / ) ) . "\n";

    if ($opts{u}){
        open ($UNCONV, ">", $opts{u}) or die "Can't open $opts{u} for writing: $!\n";
        print $UNCONV $header;
    }
}
###########################################################
sub convertToVcf{
    my $hgmd = shift;
    my @out = ();
    my @details = ();
    return if $hgmd->[ $hgmd_columns{"chromosome"} ] eq 'null';
    return if $hgmd->[ $hgmd_columns{"coordinate start"} ] eq 'null';
    my ($ref, $alt) = parseHgmdAlleles($hgmd);
    return if (not $ref and not $alt);
    my $chr = convertChrom($hgmd->[ $hgmd_columns{"chromosome"} ]);
    push @out, $chr; 
    my $pos = $hgmd->[ $hgmd_columns{"coordinate start"} ];
    if (length($alt) > length($ref)){#insertion
        #@out, $fields[ $hgmd_columns{"coordinate start"} ] -1 ;
        ($pos, $ref, $alt) = convertInsertion
            (
                $chr, 
                $ref, 
                $alt, 
                $hgmd->[ $hgmd_columns{"coordinate start"} ],
                $hgmd->[ $hgmd_columns{"coordinate end"} ],
            ); 
    }elsif (length($alt) < length($ref)){#deletion
        #push @out, $fields[ $hgmd_columns{"coordinate end"} ] ;
        ($pos, $ref, $alt) = convertDeletion
            (
                $chr, 
                $ref, 
                $alt, 
                $hgmd->[ $hgmd_columns{"coordinate start"} ],
                $hgmd->[ $hgmd_columns{"coordinate end"} ],
            ); 
    }
    push @out, $pos, $hgmd->[ $hgmd_columns{"hgmd id"} ], $ref, $alt, ".", "."; #blank QUAL and FILTER fields
    foreach my $d (@hgmd_fields){ 
        (my $tag = $d) =~ s/\s/_/g;
        (my $value = $hgmd->[ $hgmd_columns{$d} ]) =~ s/\s/_/g;
        push @details, "$tag=$value";
    }
    
    push @out, join(";", @details);
    return join("\t", @out) ;
}

###########################################################
sub convertInsertion{
    my ($chr, $ref, $alt, $start, $end) = @_;
    checkRegion ($chr, $start, $end);
    if (length($ref) > 0){
        #shouldn't be any ambiguity for a deletion-insertion (I think?)
        return ($start, $ref, $alt); 
    }
    my $s_coor = $start - (10 * length($alt));#initially look 10x the length of insertion backwards
    $s_coor = $s_coor > 0 ? $s_coor : 1;
    # $start should be the position preceding insertion if a plain old insertion/duplication
    my $seq = fetchFasta($chr, $s_coor, $start);
    my $subseq = substr($seq, length($seq) - length($alt)); 
    #for an insertion we should only need to move it backward if it's a duplication
    if ($subseq ne $alt){#i.e. not a duplication 
        $ref = substr($seq, -1);
        return ($start , $ref, "$ref$alt"); 
    }
    #dealing with a duplication now
    my $index = scanBack($seq, $alt); 
    if ($index < length($alt)){
        if ($s_coor > 1){
            return convertInsertion
                (
                    $chr, 
                    $ref, 
                    $alt, 
                    $s_coor - (10 * length($alt)),
                    $s_coor + $index,
                );
        }
    }
    my $pos = ($s_coor + $index) -1; #subtract one for our extra ref base
    my $first_base = substr($seq, $index-1, 1);
    return ($pos, "$first_base", "$first_base$alt");
}
###########################################################
sub convertDeletion{
    my ($chr, $ref, $alt, $start, $end) = @_;
    checkRegion ($chr, $start, $end);
    my $s_coor = $start - (10 * length($ref));#initially look 10x the length of deletion backwards
    $s_coor = $s_coor > 0 ? $s_coor : 1;
    my $seq = fetchFasta($chr, $s_coor, $end);
    my $subseq = substr($seq, length($seq) - length($ref), length($ref)); 
    if ($subseq ne $ref){ 
        die "Reference sequence doesn't match reference allele from variant at $chr:$start-$end\n"; 
    }
    my $index = scanBack($seq, $ref); 
    if ($index < length($ref)){
        if ($s_coor > 1){
            return convertDeletion
                (
                    $chr, 
                    $ref, 
                    $alt, 
                    $s_coor - (10 * length($ref)),
                    $s_coor + $index,
                ); 
        }
    }
    my $pos = ($s_coor + $index) -1; #subtract one for our extra ref base
    my $first_base = substr($seq, $index-1, 1);
    return ($pos, "$first_base$ref", "$first_base$alt");
}

###########################################################
sub checkRegion{
    my ($chr, $start, $end) = @_;
    if ($end < $start){
        die "ERROR retrieving fasta sequence for region '$chr:$start-$end' - start coordinate is greater than end coordinate!.\n";
    }
    if ($end > $chroms{$chr}){
        die "ERROR retrieving fasta sequence for region '$chr:$start-$end' - end coordinate is greater than chromosome length!.\n";
    }
}

###########################################################
sub getRepeatLength{
#finds shortest repeating unit of string
#if string is not repetitive returns length of string
    my $string = shift;
    my $l = length($string);
    #get first half of divisors (probable needless optimization, could just do this in one pass)
    my @div = grep{ $l % $_ == 0 } 1 .. sqrt($l);
    #get upper half of divisors 
    push @div, map {$l == $_*$_ ? () : $l/$_} reverse @div;
    for (my $i = 0; $i < @div; $i++){
        my $s = substr($string, 0, $div[$i]);
        my $j = $l/$div[$i];
        if ( $string eq ($s x $j)){
            return $div[$i];
        }
    }
    return $l;
}
###########################################################
sub scanBack{
    my ($seq, $allele) = @_;
    my $seq_length = length($seq);
    my $subpos = length($seq) - length($allele);#current matching index
    my $last_match = $subpos;
    my $r = getRepeatLength($allele); #get smallest repeating unit of $allele
    for (my $i = 1; $i * $r <= $seq_length; $i++){
    #move back by the length of the smallest repeating unit and check for match
        my $idx = $subpos - $r;
        last if $idx < 0;
        my $subseq = substr($seq, $idx, length($allele));
        if ($subseq ne $allele){
            return $subpos;
        }
        $subpos = $idx;
    }#got all the way back to beginning of seq always matching our REF allele
    return $subpos;
}

    
###########################################################
sub convertChrom{
    my $chrom = shift; 
    if ($seqs_have_chr){
        $chrom = "chr$chrom" if $chrom !~ /^chr/;
    }elsif(not $seqs_are_mixed){#do nothing if seqs are mixed style
        $chrom =~ s/^chr//;#otherwise remove prepending chr
    }
    return $chrom;
}
###########################################################
sub readFastaIndex{
    my $f = shift;
    my $faidx = "$f.fai"; 
    my %seqs = ();
    open (my $FAIDX, $faidx) or die "Can't open fata index ($faidx) for reading: $!\n";
    while (my $line = <$FAIDX>){
        my @rec = split("\t", $line); 
        $seqs{$rec[0]} = $rec[1];#contig name to contig length
        if ($rec[0] =~ /^chr/){
            $seqs_have_chr = 1;
        }elsif($seqs_have_chr){
            $seqs_are_mixed = 1;
        }
    }
    if ($seqs_are_mixed){
        print STDERR <<EOT
WARNING: Mixed chromosome styles detected. Some with prepended "chr", some without.
WARNING: Chromosome styles will not be converted for fasta retrieval. This may lead to errors.
EOT
;
    }
    return %seqs; 
}
###########################################################
sub fetchFasta{
    my ($chrom, $start, $stop) = @_;
    my $seq = $fai->fetch("$chrom:$start-$stop");
    if (length($seq) < (1 + $stop - $start)){
        die "ERROR retrieving fasta sequence for region '$chrom:$start-$stop'.\n";
    }
    return $seq;
}

###########################################################
sub revcomp{
    my $dna = shift;
    $dna = reverse($dna);
    $dna =~ tr/ACGTacgt/TGCAtgca/;
    return $dna;
}
###########################################################
sub parseHgmdAlleles{
    my $hgmd = shift;#array ref to split line from HGMD file
    my $hgvs = $hgmd->[ $hgmd_columns{hgvs} ] ;
    my $strand = $hgmd->[ $hgmd_columns{strand} ] ;
    my $ref = "";
    my $alt = "";
    if ($hgvs =~ /c\.\d+([\+-]\d+)*([ATGCN])>([ATGCN])$/i){ #SNV
        $ref = $2;
        $alt = $3;
    }elsif ($hgvs =~ /c\.\d+([\+-]\d+)*(_\d+([\+-]\d+)*)*del([ATGCN]+)$/i){
        $ref = $4;
    }elsif ($hgvs =~ /c\.\d+([\+-]\d+)*(_\d+([\+-]\d+)*)*(ins|dup)([ATGCN]+)$/i){
        $alt = $5;
    }elsif ($hgvs =~ /c\.\d+([\+-]\d+)*(_\d+([\+-]\d+)*)*del([ATGCN]+)ins([ATGCN]+)$/i){
        $ref = $4;
        $alt = $5;
    }
    if ($strand eq '-'){
        $ref = revcomp($ref);
        $alt = revcomp($alt);
    }   
    return ($ref, $alt);
}

#########################################################
sub usage{
    my $msg = shift;
    print STDERR "ERROR: $msg\n" if $msg;
    print <<EOT
    
    usage: $0 -i <hgmd_variants.txt> -f <reference fasta file> [options]

    Options:

    -i, --input FILE
        Input file of HGMD variants as obtained from HGMD professional MART. Must be tab separated values. Required.
    -f, --fasta FILE
        Fasta reference file for genome. Required.
    -o, --ouptut FILE
        VCF output file. Optional. Default is STDOUT.
    -u, --unconverted FILE
        Output file for variants that cannot be converted from HGMD. Optional.
    -s, --sort
        Use this flag to sort output in coordinate order.
    -?, -h, --help
        Show this help message.

EOT
;
    exit 1 if $msg;
    exit;
}
