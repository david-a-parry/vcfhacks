#!/usr/bin/env perl 
use strict;
use warnings;
use POSIX qw/strftime/;
use Getopt::Long;
use Bio::DB::Sam;
use Data::Dumper;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use FindBin qw($RealBin);
use lib "$RealBin/../lib/dapPerlGenomicLib";
use VcfReader 0.3;

my %opts = (b => 1000000);
GetOptions(
    \%opts,
    'i|input=s',        #vcf input
    'o|output=s',       #optional output file
    'x|mismatches=s',   #optional mismatches output file
    'f|fasta=s',        #genome fasta
    'b|seq_buffer=i',   #min no. bases to retrieve at a time
    'progress',
    'h|?|help',
) or die ("Syntax error.\n"); 

usage() if $opts{h};
usage("-i/--input is required" ) if (not $opts{i});
usage("-f/--fasta is required" ) if (not $opts{f});

my %seq = (); #buffer of 1Mb sequence from fasta

my $fai = Bio::DB::Sam::Fai->load($opts{f});#should create index if it doesn't exist
#get hash of chromosomes to lengths from index and check naming convention
my %contigs = readFastaIndex($opts{f}); 
 
my $OUT;
if ( $opts{o} ) {
    open( $OUT, ">$opts{o}" )
      || die "Can't open $opts{o} for writing: $!";
}
else {
    $OUT = *STDOUT;
}
my $MISMATCHES;
if ( $opts{x} ) {
    open( $MISMATCHES, ">$opts{x}" )
      || die "Can't open $opts{x} for writing: $!";
}

my $time = strftime( "%H:%M:%S", localtime );
print STDERR "[$time] INFO - Initializing input VCF...\n";
my ($header, $first_var, $VCF)  = VcfReader::getHeaderAndFirstVariant($opts{i});
die "Header not ok for input ($opts{i}) "
    if not VcfReader::checkHeader( header => $header );
my $total_vcf = 0;
my $next_update = 0;
my $prog_total = -1;
my $progressbar;
if ( defined $opts{progress} and not -p $opts{i} and  $opts{i} ne '-') {
    $total_vcf = VcfReader::countVariants( $opts{i} );
    $time = strftime( "%H:%M:%S", localtime );
    print STDERR "\n[$time] INFO - $opts{i} has $total_vcf variants.\n";
}elsif(defined $opts{progress}){
    $time = strftime( "%H:%M:%S", localtime );
    print STDERR "\n[$time] INFO - Input is from STDIN or pipe - will report progress per 10000 variants.\n";
}
if ( defined $opts{progress}){
    my $name = "Filtering";
    $progressbar = Term::ProgressBar->new(
        { 
          name => $name,
          count => $total_vcf,
          ETA => "linear" 
        } 
    );
}
my $n = 0;
print $OUT join("\n", @$header) . "\n";
print $MISMATCHES join("\n", @$header) . "\n" if $MISMATCHES;

processLine($first_var);
VAR: while (my $line = <$VCF> ) {
    processLine($line);
}
close $VCF;
close $OUT;
close $MISMATCHES if $MISMATCHES;

##################################################
sub processLine{
    my $line = shift;
    return if $line =~ /^#/;
    $n++;
    checkProgress(1);#always check progress here
    chomp $line;
    #our VcfReader methods should be more efficient on pre-split lines
    my @split_line = split( "\t", $line );
    my ($ref, $chrom, $pos) = VcfReader::getMultipleVariantFields
    (
        \@split_line, 
        'REF', 
        'CHROM',
        'POS',
    );
    my $span = VcfReader::getSpan(\@split_line);
    
    if (not %seq or varNotInRegion($chrom, $pos, $span)){
        %seq = getSeqForVar($chrom, $pos, $span);
    }
    if (refMatches($chrom, $pos, $ref)){
        print $OUT "$line\n";
    }elsif($MISMATCHES){
        print $MISMATCHES "$line\n";
    }

}

##################################################
sub refMatches{
    my ($chrom, $pos, $ref) = @_;
    return 0 if not %seq;
    my $refseq = substr($seq{dna}, $pos - $seq{start}, length($ref));
    if (uc($refseq) eq uc($ref)){
        return 1;
    }
    return 0;
}
##################################################
sub getSeqForVar{
    my ($chrom, $start, $end) = @_;
    if ( 1 + $end - $start < $opts{b}){#get at a minimum $opts{b} seq at a time
        $end = $start + $opts{b} - 1;
    }
    $end = $end <= $contigs{$chrom} ? $end : $contigs{$chrom};
    my $dna = $fai->fetch("$chrom:$start-$end");
    if (not $dna){
        warn "Could not retrieve DNA for $chrom:$start-$end\n";
        return ();
    }
    return (dna => $dna, chrom => $chrom, start => $start, end => $end);
}

##################################################
sub varNotInRegion{
    my ($chrom, $start, $end) = @_;
    return 1 if ($chrom ne $seq{chrom});
    return 1 if ($start < $seq{start});
    return 1 if ($end > $seq{end});
    return 0;
}

##################################################
sub checkProgress{
    return if not $progressbar;
    my $do_count_check = shift;
    if ($prog_total > 0){
        $next_update = $progressbar->update($n) if $n >= $next_update;
    }elsif($do_count_check){#input from STDIN/pipe
        if (not $n % 10000) {
            my $time = strftime( "%H:%M:%S", localtime );
            $progressbar->message( "[INFO - $time] $n variants read" );
        }
    }
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
    }
    return %seqs; 
}

###########################################################
sub usage{
    my $msg = shift;

    print STDERR "\n$msg\n" if $msg;

    print <<EOT

removeNonMatchingRefs.pl

Removes variants from a VCF where the REF allele does not match with the 
provided FASTA file.

The reason I wrote this script was in order to process the GRCh38 versions of 
the 1000 genomes VCFs as found here: 

ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/

There are a significant fraction of variants present in these VCFs 
for which the REF allele has changed between version 37 and 38 of the human 
genome, yet these changes are not represented in the VCF. These alleles prevent 
merging of the 1000 genomes data with true GRCh38/hg38 VCFs. This script removes
these problem variants (which can optionally be outputted to a separate file) in
order to be able to use the output with other data. 

In theory this can be used for other genomes as well should you have datasets 
with similar problems. This is NOT a means for converting a VCF between 
different genome assemblies; it is only a means of removing problem variants 
from a VCF.

Usage: $0 -i input.vcf -f GRCh38.fa [-o output.vcf -m mismatches.vcf] [options]


Options: 
        -i,--input
            Input vcf
        
        -o,--output
            Optional output file. Default = STDOUT
        
        -m,--mismatches
            Optional output file for variants that do not match the FASTA file.

        -f,--fasta
            Reference FASTA file.

        -b,--seq_buffer
            Min no. bases to retrieve at a time from the fasta file. A higher
            value here may decrease runtime at the expense of using more RAM.
            Default = 1000000.
    
        -p,--progress
            Show a progress bar.
    
        -h,?,--help
            Show this message and exit

EOT
;
    exit 1 if $msg;
    exit;
}
