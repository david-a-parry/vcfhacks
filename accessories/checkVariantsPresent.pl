#!/usr/bin/env perl
use warnings;
use strict;
use Pod::Usage;
use Data::Dumper;
use Bio::DB::Sam;
use POSIX qw/strftime/;
use FindBin qw($RealBin);
use Getopt::Long;
use lib "$RealBin/../lib/dapPerlGenomicLib";
use VcfReader 0.3;

my %opts = (gq => 0);
GetOptions(
    \%opts,
    "output|o=s",
    "input|i=s",
    "t|table=s",
    "f|fasta=s",
    "v|vcf=s",
    "gq=i",
    "help|h|?",
    "manual",
) or pod2usage( -exitval => 2, -message => "Syntax error" );

pod2usage( -verbose => 2, -exitval => 0  ) if $opts{manual};

pod2usage( -verbose => 1, -exitval => 0  ) if $opts{help};

pod2usage
( 
    -exitval => 2, 
    -message => "Missing required argument", 
) if not $opts{input} or (not $opts{t} and not $opts{v}) ;

pod2usage
( 
    -exitval => 2, 
    -message => "--fasta is required when using -t/--table option",
) if $opts{t} and not $opts{f};

my $time = strftime( "%H:%M:%S", localtime );
print STDERR "[$time] INFO - Initializing input VCF...\n";
my @header = VcfReader::getHeader($opts{input});
die "ERROR: Invalid VCF header in $opts{input}\n" 
  if not VcfReader::checkHeader(header => \@header);
my %sample_to_col = VcfReader::getSamples
(
    header => \@header,
    get_columns => 1,
);
my %search_args = VcfReader::getSearchArguments( $opts{input});
$time = strftime( "%H:%M:%S", localtime );
print STDERR "[$time] INFO - Finished initializing input VCF\n";
my $fai;
if ($opts{t}){
    $fai = Bio::DB::Sam::Fai->load($opts{f});#should create index if it doesn't exist
    parseTable();
}

if ($opts{v}){
    parseVcf();
}

##################################################
sub parseVcf{
    my ($head, $first_var, $FH) = VcfReader::getHeaderAndFirstVariant($opts{v});
    die "Header not ok for input ($opts{v}) "
        if not VcfReader::checkHeader( head => $head );
    my %v_sample_to_col = VcfReader::getSamples
    (
        header => $head,
        get_columns => 1,
    );
    processVcfLine($first_var);
    while (my $l = <$FH>){
        processVcfLine($l, \%v_sample_to_col);
    }
    close $FH;
}

##################################################
sub processVcfLine{
    my $line = shift;
    my $v_samp_to_col = shift;
    chomp $line;
    my @split = split("\t", $line); 
    my $chrom = VcfReader::getVariantField(\@split, 'CHROM');
    my $start = VcfReader::getVariantField(\@split, 'POS');
    my $filt_span = VcfReader::getSpan(\@split);
    my @hits = VcfReader::searchByRegion
    (
        %search_args,
        chrom => $chrom,
        start => $start,
        end   => $filt_span,
    );
    my ($f_sample_vars, @alleles) = getSamplesWithVariant(\@split, $v_samp_to_col);
    my %f_min_vars = VcfReader::minimizeAlleles(\@split);
    #TODO process each genotype and look for matching genotypes in hits 
    #- process one allele at a time?
    foreach my $h (@hits){
        my @h_split = split("\t", $h);
        my %h_min_vars = VcfReader::minimizeAlleles(\@h_split);
        ...
    }
    
}

##################################################
sub getSamplesWithVariant{
    my $var = shift;
    my $v_samp_to_col = shift;
    
    my %samp_to_gt = VcfReader::getSampleCall
    (
        line => $var,
        minGQ => $opts{gq},
        all => 1,
        sample_to_columns => $v_samp_to_col
    );
    my %samp_to_ad = VcfReader::getSampleGenotypeField
    (
        field => "AD",
        all => 1,
        sample_to_columns => $v_samp_to_col
    );
    ...
}
    
##################################################
sub parseTable{
    open (my $TAB, $opts{t}) or die "Could not open $opts{t} for reading: $!\n";
    while (my $line = <$TAB>){
        chomp $line;
        my ($chrom, $start, $end, $ref, $alt, $sample, $genotype) = split("\t", $line);
        if ($genotype =~ /^het/i){
            $genotype = 'het';
        }elsif ($genotype =~ /^hom/i){
            $genotype = 'hom';
        }else{
            die "Don't understand genotype '$genotype' for line:\n$line\n";
        }
        my $allele_match = 0;
        my $gt_match     = 0;
        my $var = tabToVcfAlleles($ref, $alt, $chrom, $start);
        if (exists ($sample_to_col{$sample}) ){
            my @hits = VcfReader::searchForPosition
            (
                %search_args,
                chrom => $chrom,
                pos   => $var->{POS},
            );
            foreach my $h (@hits){
                my @split = split("\t", $h);
                my %min_vars = VcfReader::minimizeAlleles(\@split);
                if ( my $match = getAlleleMatch
                    ( 
                        $var,
                        \%min_vars, 
                    )
                ) {
                    my $gt = VcfReader::getSampleCall
                    (
                        line => \@split, 
                        sample => $sample,
                        sample_to_columns => \%sample_to_col,
                        
                    );
                    my @g = split(/[\/\|]/, $gt);
                    if (grep {$_ eq $match} @g){ 
                        $allele_match = 1;
                    }
                    if ($allele_match and $genotype eq 'hom'){
                        $gt_match = ($g[0] eq $match and $g[1] eq $match);
                    }elsif ($allele_match and $genotype eq 'het'){
                        $gt_match = (
                            ($g[0] eq $match and $g[1] ne $match) or
                            ($g[1] eq $match and $g[0] ne $match) 
                        );
                    }
                    $gt_match ||= 0;
                }
            }
        }
        print "$line\t$allele_match\t$gt_match\n";
    }
    close $TAB;
}


##################################################
sub getAlleleMatch{
    my $var = shift;
    my $multivar = shift;
    foreach my $k (keys %{$multivar}){
        next if $var->{CHROM} ne $multivar->{$k}->{CHROM};
        next if $var->{POS} ne $multivar->{$k}->{POS};
        next if $var->{REF} ne $multivar->{$k}->{REF};
        next if $var->{ALT} ne $multivar->{$k}->{ALT};
        return $k;
    }
}

##################################################
sub tabToVcfAlleles{
    my ($ref, $alt, $chrom, $pos) = @_;
    if ($ref eq '-'){
        $ref = $fai->fetch("$chrom:$pos-$pos");
        $alt = "$ref$alt";
    }elsif($alt eq '-'){
        $pos--;
        $alt = $fai->fetch("$chrom:$pos-$pos");
        $ref = "$alt$ref";
    }else{
         ($pos, $ref, $alt) = VcfReader::reduceRefAlt($pos, $ref, $alt);
    }
    my %var = 
    (
        CHROM => $chrom,
        POS   => $pos,
        REF   => $ref,
        ALT   => $alt,
    );
    return \%var;
}   
        

