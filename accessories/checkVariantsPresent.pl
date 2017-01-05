#!/usr/bin/env perl
use warnings;
use strict;
use Pod::Usage;
use Data::Dumper;
use Bio::DB::Sam;
use POSIX qw/strftime/;
use FindBin qw($RealBin);
use List::Util qw /sum/;
use Getopt::Long;
use lib "$RealBin/../lib/dapPerlGenomicLib";
use VcfReader 0.3;

my %opts = 
(
    gq => 0,
    d  => 0,
    b  => 0,
);
GetOptions(
    \%opts,
    "o|output=s",
    "s|summary=s",
    "input|i=s",
    "t|table=s",
    "f|fasta=s",
    "v|vcf=s",
    "gq=i",
    "d|dp=i",
    "b|ab=f",
    "help|h|?",
    "manual",
) or pod2usage( -exitval => 2, -message => "Syntax error" );

pod2usage( -verbose => 2, -exitval => 0  ) if $opts{manual};

pod2usage( -verbose => 1, -exitval => 0  ) if $opts{help};

pod2usage
( 
    -exitval => 2, 
    -message => "-i/--input and either -t/--table or -v/--vcf arguments are required", 
) if not $opts{input} or (not $opts{t} and not $opts{v}) ;

pod2usage
( 
    -exitval => 2, 
    -message => "-o/--output and -s/summary arguments are required", 
) if not $opts{o} or not $opts{s} ;

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

open (my $OUT, ">", $opts{o}) or die "Could not open $opts{o} for reading: $!\n";
open (my $SUM, ">", $opts{s}) or die "Could not open $opts{s} for reading: $!\n";
#print header
print $OUT join
(
    "\t", 
    qw/
        CHROM
        POS
        REF
        ALT
        SAMPLE
        GT1
        GT2
        MATCH
        MISMATCH_TYPE
    /
) . "\n"; 

my %samp_summary = 
    map { $_ => 
        {
            TP => 0,
            TN => 0,
            FP => 0,
            FN => 0,
            M  => 0,
        } 
    } keys %sample_to_col;
   
$time = strftime( "%H:%M:%S", localtime );
print STDERR "[$time] INFO - starting to process variants...\n";
if ($opts{t}){
    $fai = Bio::DB::Sam::Fai->load($opts{f});#should create index if it doesn't exist
    parseTable();
}

if ($opts{v}){
    parseVcf();
}
close $OUT;
outputSampleSummary();
close $SUM;
##################################################
sub outputSampleSummary{
    print $SUM join
    (
        "\t",
        qw(
            SAMPLE 
            FRAC_CONCORDANT
            FRAC_DISCORDANT
            TRUE_POS/NON-TRUE_NEG
            CONCORDANT
            DISCORDANT
            TRUE_POS
            TRUE_NEG
            FALSE_POS
            FALSE_NEG
            MISMATCH
        )
    ) . "\n";
    foreach my $s (sort keys %samp_summary){
        my $conc = $samp_summary{$s}->{TP} + $samp_summary{$s}->{TN}; 
        my $disc = $samp_summary{$s}->{FP} + $samp_summary{$s}->{FN} + $samp_summary{$s}->{M}; 
        my $per_conc = 0;
        my $per_disc = 0;
        my $tp_over_exp = 'N/A';
        my $total = $conc + $disc;
        if ($total > 0){
            $per_disc = $disc/$total;
            $per_conc = $conc/$total;
        }
        my $total_non_ref = $samp_summary{$s}->{FN} + 
                            $samp_summary{$s}->{FP} +
                            $samp_summary{$s}->{TP};
        if ($total_non_ref > 0){
            $tp_over_exp =  $samp_summary{$s}->{TP}/$total_non_ref;
                            
        }
        print $SUM join
        (
            "\t",
            $s,
            $per_conc,
            $per_disc,
            $tp_over_exp,
            $conc,
            $disc,
            $samp_summary{$s}->{TP},
            $samp_summary{$s}->{TN},
            $samp_summary{$s}->{FP},
            $samp_summary{$s}->{FN},
            $samp_summary{$s}->{M},
        ) . "\n";
    }
}

##################################################
sub parseVcf{
    my ($head, $first_var, $FH) = VcfReader::getHeaderAndFirstVariant($opts{v});
    die "Header not ok for input ($opts{v}) "
        if not VcfReader::checkHeader( header => $head );
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
    my @alleles = VcfReader::readAlleles(line => \@split);
    return if @alleles > 2;#do biallelic variants only
    my $filt_span = VcfReader::getSpan(\@split);
    my @hits = VcfReader::searchByRegion
    (
        %search_args,
        chrom => $chrom,
        start => $start,
        end   => $filt_span,
    );
    my %f_sample_vars = getSampleGenotypes(\@split, $v_samp_to_col);
    my %f_min_vars = VcfReader::minimizeAlleles(\@split);
    my %comp_gts = ();
HIT: foreach my $h (@hits){
        my @h_split = split("\t", $h);
        my @h_alleles = VcfReader::readAlleles(line => \@h_split);
        my %h_min_vars = VcfReader::minimizeAlleles(\@h_split);
        if (my $m = getAlleleMatch($f_min_vars{1}, \%h_min_vars) ){
            my %h_sample_vars = getSampleGenotypes(\@h_split, \%sample_to_col);
            foreach my $s ( sort keys %sample_to_col ) { 
                next if not $f_sample_vars{$s}; #sample not in vcf
                next if $f_sample_vars{$s} !~ /\d[\/\|]\d/;#no call for this sample
                next if $h_sample_vars{$s} !~ /\d[\/\|]\d/;#no call for this sample
                my @f_al = split(/[\/\|]/, $f_sample_vars{$s});
                my @h_al = split(/[\/\|]/, $h_sample_vars{$s});
                $comp_gts{$s}->{i} = "$h_alleles[$h_al[0]]/$h_alleles[$h_al[1]]";
                $comp_gts{$s}->{v} = "$alleles[$f_al[0]]/$alleles[$f_al[1]]";
            }
            last HIT;#only expect one matching variant
        }
    }
    foreach my $s (sort keys %comp_gts){
        my $match = 0;
        my $mismatch = 0;
        my $ref = $alleles[0];
        if ($comp_gts{$s}->{i} eq $comp_gts{$s}->{v}){
            #matching genotypes
            $match = 1;
            if ($comp_gts{$s}->{i} =~ /$ref\Q\E[\/\|]$ref/){
                $samp_summary{$s}->{TN}++;
            }else{
                $samp_summary{$s}->{TP}++;
            }
        }else{
            if ($comp_gts{$s}->{v} =~ /$ref\Q\E[\/\|]$ref/){
                #false positive
                $mismatch = 'FP';
                $samp_summary{$s}->{FP}++;
            }elsif ($comp_gts{$s}->{i} =~ /$ref\Q\E[\/\|]$ref/){
                #false negative
                $mismatch = 'FN';
                $samp_summary{$s}->{FN}++;
            }else{
                #mismatch
                $mismatch = 'M';
                $samp_summary{$s}->{M}++;
            }
        }
        print $OUT join
        (
            "\t", 
            $chrom,
            $start,
            $alleles[0],
            $alleles[1],
            $s,
            $comp_gts{$s}->{v},
            $comp_gts{$s}->{i},
            $match,
            $mismatch,
        ) . "\n";
    }
}

##################################################
sub getSampleGenotypes{
    #returns hash of samples to genotypes
    my $var = shift;
    my $v_samp_to_col = shift;
    
    my %samp_to_gt = VcfReader::getSampleCall
    (
        line => $var,
        minGQ => $opts{gq},
        all => 1,
        sample_to_columns => $v_samp_to_col
    );
    return %samp_to_gt if not $opts{d} and not $opts{b};

    my %samp_to_ad = VcfReader::getSampleGenotypeField
    (
        line => $var,
        field => "AD",
        all => 1,
        sample_to_columns => $v_samp_to_col
    );
    foreach my $s (keys %samp_to_gt){
        my $dp = 0;
        my @ad = ();
        next if $samp_to_gt{$s} eq './.';
        if ($samp_to_ad{$s} and $samp_to_ad{$s} ne '.'){
            @ad = split(",", $samp_to_ad{$s});
            $dp = sum(@ad);
        }
        if ($opts{d} and $dp < $opts{d}){
            $samp_to_gt{$s} = "./.";
            next;
        }
        if ($opts{b}){
            my @new_alleles = (); 
            foreach my $allele (split(/[\/\|]/, $samp_to_gt{$s}) ){
                my $ab = 0;
                if ($dp > 0){
                    $ab = $ad[$allele]/$dp;
                }
                if ($opts{b} > $ab){
                    push @new_alleles, '.';
                }else{
                    push @new_alleles, $allele;
                }
            }
            $samp_to_gt{$s} = join("/", @new_alleles);
        }
    }
    return %samp_to_gt;
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
        my $mismatch     = 0;
        my $i_gt         = '';
        my $t_gt         = "$ref/$alt";
        if ($genotype eq 'hom'){
            $t_gt = "$alt/$alt";
        }
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
                        minGQ => $opts{gq},
                        
                    );
                    if ($opts{d} or $opts{b}){
                        my $allele_depth = VcfReader::getSampleGenotypeField
                        (
                            line => \@split, 
                            field => "AD",
                            sample => $sample,
                            sample_to_columns => \%sample_to_col
                        );
                        my @ad = ();
                        my $dp = 0;
                        if (defined $allele_depth and $allele_depth ne '.'){
                            @ad = split(",", $allele_depth);
                            $dp = sum(@ad);
                        }
                        if ($opts{d} > $dp){
                            $gt = './.';
                        }
                        if ($opts{b}){
                            if ($opts{d} > 0){
                                my $ab = $ad[$match]/$dp;
                                if ($opts{b} > $ab){
                                    $gt = './.';
                                }
                            }else{
                                $gt = './.';
                            }
                        }
                    }
                    my @g = split(/[\/\|]/, $gt);
                    if ($gt eq './.'){
                        $i_gt = $gt;
                    }else{
                        my @alleles = VcfReader::readAlleles(line => \@split);
                        $i_gt = "$alleles[$g[0]]/$alleles[$g[1]]";
                    }
                    if (grep {$_ eq $match} @g){
                        $allele_match = 1;
                    }
                    if ($allele_match and $genotype eq 'hom'){
                        $gt_match = ($g[0] eq $match and $g[1] eq $match);
                    }elsif ($allele_match and $genotype eq 'het'){
                        $gt_match = 
                        (
                            ($g[0] eq $match and $g[1] ne $match) or
                            ($g[1] eq $match and $g[0] ne $match) 
                        );
                    }
                    $gt_match ||= 0;
                    if ($gt_match){
                        $samp_summary{$sample}->{TP}++;
                    }else{
                        if ($gt =~ /^0[\/\|]0$/){
                            $mismatch = 'FN';
                            $samp_summary{$sample}->{FN}++;
                        }else{
                            $mismatch = 'M';
                            $samp_summary{$sample}->{M}++;
                        }
                    }
                }
            }
        }
        print $OUT join
        (
            "\t", 
            $chrom,
            $start,
            $ref,
            $alt,
            $sample,
            $t_gt,
            $i_gt,
            $gt_match,
            $mismatch,
        ) . "\n";


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
        if ($var->{ALT} ne $multivar->{$k}->{ALT}){
            next if 
            (
                $var->{ALT} ne '<NON_REF>' 
                and  
                $multivar->{$k}->{ALT} ne '<NON_REF>'
            );
        }
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
        

