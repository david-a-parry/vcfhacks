#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use POSIX qw/strftime/;
use Data::Dumper;
use Statistics::R;
use List::Util qw ( first max min ) ;
use File::Temp qw/ tempfile /;
use FindBin qw($RealBin);
use lib "$RealBin/lib";
use lib "$RealBin/lib/dapPerlGenomicLib";
use VcfReader 0.3;
use VcfhacksUtils;
use ParsePedfile;

my $R = Statistics::R->new() ;
$R->start_sharedR ;
my $progressbar;
my @classes = (); 
my %default_classes = (); 
my @add_classes = ();
my @damaging = ();
my @biotypes = ();
my @custom_af = ();
my @eval_filters = (); 

my %opts = 
(
    classes       => \@classes,
    add_classes   => \@add_classes,
    d             => \@damaging,
    biotypes      => \@biotypes,
    j             => \@custom_af,
    eval_filters  => \@eval_filters,
);

GetOptions(
    \%opts,
    "add_classes=s{,}" ,
    "a|af=f" ,
    "biotypes=s{,}",
    "b|progress" ,
    "canonical_only" ,
    "classes=s{,}" ,
    "consensus_splice_site",
    "c|cadd_filter=f",
    "d|damaging=s{,}" ,
    "eval_filters=s{,}",
    "f|ped_file=s",
    "g|gq=i",
    "h|?|help" ,
    "i|input=s" ,
    "j|custom_af=s{,}",
    "k|keep_any_damaging" ,
    "manual" ,
    "m|mode=s" ,
    "no_biotype_filtering",
    "no_indels",
    "o|output=s" ,
    "pass_filters" ,
    "pl=f", 
    "skip_unpredicted" ,
    "s|var_summaries=s",
    "v|var_quality=i",
    "u|max_sample_allele_frequency=f",
) or pod2usage(-message => "Syntax error", -exitval => 2);
pod2usage(-verbose => 2, -exitval => 0) if ($opts{manual});
pod2usage(-verbose => 1, -exitval => 0) if ($opts{h});

if (not $opts{i}){
    pod2usage
    (
        -message => "Syntax error - an input file (-i/--input) must be supplied.\n", 
        -exitval => 2
    );
}
if (not $opts{f}){
    pod2usage
    (
        -message => "Syntax error - a ped/fam file (-f/--ped_file) must be supplied.\n", 
        -exitval => 2
    );
}

if (defined $opts{af} && ($opts{af} < 0 or $opts{af} > 1.0)){
    pod2usage
    (
        -message =>  "-a/--allele_frequency option requires a value between 0.00 and 1.00 to ".
                     "filter on allele frequency.\n", 
        -exitval => 2
    );
}

if (defined $opts{u} ){
    if ($opts{u} < 0 or $opts{u} > 1.0){
        pod2usage
        (
            -message =>  "-u/--max_sample_allele_frequency option requires a ". 
                         "value between 0.00 and 1.00 to filter on allele ".
                         "frequency.\n", 
            -exitval => 2
        );
    }
}

#GQs are >= 0
$opts{v} = defined $opts{v} ? $opts{v} : 0;
pod2usage(
    -message => "SYNTAX ERROR: variant quality scores (-v/--var_quality)".
                " must be 0 or greater.\n",
    -exitval => 2
) if ( $opts{v} < 0 );

$opts{g} = defined $opts{g} ? $opts{g} : 20;
pod2usage(
    -message => "SYNTAX ERROR: Genotype quality (-g/--gq) ".
                "scores must be 0 or greater.\n",
    -exitval => 2
) if $opts{g} < 0;

if (defined $opts{pl}){
    pod2usage(
        -message => "SYNTAX ERROR: Genotype quality (--pl) ".
                    "scores must be 0 or greater.\n",
        -exitval => 2
    ) if $opts{pl} < 0;
}

#CADD phred score filter is >= 0 (0 means no CADD filtering)
if (defined $opts{c}){
    pod2usage(
        -message => "SYNTAX ERROR: -c/--cadd_filter cannot be a negative value.\n",
        -exitval => 2
    ) if $opts{c} < 0;
}
    
#get and check header

my @header = VcfReader::getHeader($opts{i});
die "ERROR: Invalid VCF header in $opts{i}\n" 
  if not VcfReader::checkHeader(header => \@header);
my %sample_to_col = VcfReader::getSamples
(
    header => \@header,
    get_columns => 1,
);

#read ped file and get case and control samples
my @samples  = ();
my @cases    = ();
my @controls = ();
getCasesAndControls();

#get available INFO fields from header
my %info_fields = VcfReader::getInfoFields(header => \@header);

#check mode

if ($opts{m}){
    if ( lc($opts{m}) ne 'vep' and lc($opts{m}) ne 'snpeff' ){
        die <<EOT
SYNTAX ERROR: Unrecognised value for --mode: '$opts{m}'. 

Valid values are 'vep' or 'snpeff'. 

EOT
;
    }
    $opts{m} = lc($opts{m}); 
}

#check VEP/SNPEFF header and get annotation order 
# - will set $opts{m} if not already defined
my %csq_header = getAndCheckCsqHeader();
  # hash of functional annotations from SnpEff/VEP to their annotation index
    
my $feature_id  = $opts{m} eq 'vep' ? "feature"     : "feature_id";
my $symbol_id   = $opts{m} eq 'vep' ? "symbol"      : "gene_name";
my $gene_id     = $opts{m} eq 'vep' ? "gene"        : "gene_id";
my $hgvsc       = $opts{m} eq 'vep' ? "hgvsc"       : "hgvs.c";
my $hgvsp       = $opts{m} eq 'vep' ? "hgvsp"       : "hgvs.p";
my $biotype     = $opts{m} eq 'vep' ? "biotype"     : "transcript_biotype";
my $consequence = $opts{m} eq 'vep' ? "consequence" : "annotation";

#set default consequence fields to retrieve 
my @csq_fields = getCsqFields();#consequence fields to retrieve from VCF

#and check variant consequence classes are acceptable
my %class_filters = map { $_ => undef } getAndCheckClasses();

#check in silico prediction classes/scores are acceptable
my %in_silico_filters = getAndCheckInSilicoPred();
  #hash of prediction program names and values to filter
my @eval_exp  = getAndCheckEvalFilters();

#and check biotype classes are acceptable
my %biotype_keep = map { $_ => undef } getAndCheckBiotypes();


# get available allele frequency annotations if filtering on AF
# NOTE: we do not use VEP annotations as they do not necessarily match the variant allele
my %af_info_fields = (); #check available INFO fields for AF annotations
if ( defined $opts{a} ) {
    %af_info_fields = getAfAnnotations(\%info_fields);   
}

#if filtering on CADD score check we have a CaddPhredScore header
if ( $opts{c} ){
    if (not exists $info_fields{CaddPhredScore}){
        informUser("WARNING: No 'CaddPhredScore' INFO field found in header ".
         "- your input probably does not contain annotations for filtering on" . 
         " CADD score.\n");
    }
}

#if using -t/--target_genes argument, check file

my ($OUT, $VAR) = openOutput();

#write header and program run parameters
writeOptionsToHeader();

my %contigs = (); 
my %counts  = ();   #counts for each transcript
my %trans_info = ();#store gene id and symbol id for each transcript we encounter
#keep track of contigs we've seen so we can check input is (roughly) sorted 
#... (only relevant if using --find_shared_genes option)
#progress bar vars
my $next_update      = 0;
my $total_vars       = 0; 
my $var_count        = 0;
my $functional_count = 0;
if ($opts{b}) {
    $progressbar = VcfhacksUtils::getProgressBar(input => $opts{i});
}

processByLine();

#################################################
sub getCasesAndControls{
    my $ped = ParsePedfile->new( file => $opts{f});
    #ParsePedfile should warn and croak if duplicate samples found
    push @cases, $ped->getAllAffecteds();       
    push @controls, $ped->getAllUnaffecteds();
    @samples = (@cases, @controls); 
    foreach my $s (@samples){
        if (not exists $sample_to_col{$s}){
            die "ERROR: Sample '$s' from $opts{f} does not exist in VCF!\n";
        }
    }
    die "ERROR: No cases identified from $opts{f}\n" if not @cases;
    die "ERROR: No controls identified from $opts{f}\n" if not @controls;
}
    
#################################################
sub processByLine{
    my $VCF = VcfReader::openVcf($opts{i}); 
    #read line
    while (my $line = <$VCF>){
        processLine($line);
        updateProgressBar();  
    }

    close $VCF;
    outputCounts();
    updateProgressBar();  
    close $OUT;

    if ($opts{b} and not $progressbar){
        print STDERR "\r$var_count variants processed...\n";
    }
    informUser
    (
        "$functional_count variants matching criteria found from $var_count total".
        " variants.\n" . 
        keys(%trans_info) . " transcripts analyzed.\n"
    );
}

#################################################
sub processLine{
    my $line = shift;
    return if $line =~ /^#/;#skip header
    $var_count++;#for progress bar
    chomp $line;
    my @split = split("\t", $line, 9);#do not need GT fields in the first 
    my ($chrom, $qual, $filter) = VcfReader::getMultipleVariantFields
    (
        \@split, 
        'CHROM', 
        'QUAL',
        'FILTER',
    );
    #skip whole line if not passing filter or variant quality settings
    if ($opts{pass_filters}){
        return if $filter ne 'PASS';
    }
    
    if ($opts{v}){
        return if $qual < $opts{v};
    }
    #get alleles to arrays of 'functional' consequences
    if (isNewChromosome($chrom)){
        #if new chrom check genes and clear collected data
        outputCounts();
    }
    $functional_count += countFunctionalCsq(\@split);
}

#################################################
sub isNewChromosome{
    my $chrom = shift;
    if (exists $contigs{$chrom} ){
        #require chroms to be ordered together
        if ($contigs{$chrom} != scalar(keys %contigs) - 1){
            die <<EOT
Encountered a variant for contig $chrom with at least one variant for a
different contig inbetween. Your VCF input must be sorted such that all contigs
are kept together when using the -f/--find_shared_genes or --gene_counts
options. Please sort your input and try again.
EOT
            ;
        }
        return 0;
    }
    #if new chrom record its contig index and return true
    $contigs{$chrom} = scalar(keys%contigs);
    return $contigs{$chrom} > 0;#return false if this is the first
                                #chrom we've encountered    
}

#################################################
sub countFunctionalCsq{
    my $split = shift;
    my %allele_to_csq = (benign => {}, damaging => {}); 
                            #keys are allele numbers, values = hash of 
                            #'functional' gene/#transcript IDs to 
                            #classification and counts
    my %benign_alleles = (); #override any other annotations for these alleles
    #get Allele frequency annotations if they exist and we're filtering on AF
    my %af_info_values = (); 
    if (%af_info_fields){ 
        %af_info_values = getAfInfoValues($split);
    }
    
    #get alleles and consequences for each allele
    my @alleles = VcfReader::readAlleles(line => $split);
    my @csq = getConsequences($split);
    my @alts_to_vep_allele = ();#get VEP style alleles if needed
    if ($opts{m} eq 'vep'){
        @alts_to_vep_allele = VcfReader::altsToVepAllele
        (
            ref => $alleles[0],
            alts => [ @alleles[1..$#alleles] ], 
        );
    }
    my @cadd_scores = (); 
    if ($opts{c}){
        my $cadd = VcfReader::getVariantInfoField
        (
            $split,
            'CaddPhredScore',
        );
        if ($cadd){
            @cadd_scores = split(",", $cadd);
        }
    }
    #assess each ALT allele (REF is index 0)
ALT: for (my $i = 1; $i < @alleles; $i++){
        next ALT if $alleles[$i] eq '*';
        if ($opts{no_indels}){
            next if length($alleles[0]) != length($alleles[$i]);
        }
        #if CADD filtering skip if less than user-specified cadd filter
        #if no CADD score for allele then we won't skip unless using --
        if (@cadd_scores){
            my $score = $cadd_scores[$i - 1];
            if ($score ne '.'){
                $benign_alleles{$i}++ if $score < $opts{c};
            }elsif ($opts{skip_unpredicted}){
                $benign_alleles{$i}++;
            }
        }
        
        #filter if AF annotations > MAF
        if (%af_info_fields){ #check for annotateSnps.pl etc. frequencies
           $benign_alleles{$i}++ if ( alleleAboveMaf($i - 1, \%af_info_values) );
        }

        #get all consequences for current allele
        my @a_csq = ();
        if ($opts{m} eq 'vep'){
            @a_csq = grep { $_->{allele} eq $alts_to_vep_allele[$i-1] } @csq;
        }else{
            @a_csq = grep { $_->{allele} eq $alleles[$i] } @csq;
        }
    
 
        #benign if VEP/SNPEFF annotation is not the correct class 
        # (or not predicted damaging if using that option)
        # and allele is not flagged as pathogenic in ClinVar
CSQ:    foreach my $annot (@a_csq){
            #intergenic variants have no feature associated with them - skip
            next CSQ if $annot->{$consequence} eq "intergenic_variant";
            #skip unwanted biotypes
            next CSQ if (not exists $biotype_keep{$annot->{$biotype}}) ;
            #skip non-canonical transcripts if --canonical_only selected
            if ($opts{canonical_only}) {
                next CSQ if ( not $annot->{canonical} );
            }
            my $t = $annot->{$feature_id};# get gene and symbol IDs for all 
                                          # qualifying transcripts
            if (not $trans_info{$t}){
                $trans_info{$t}->{gene} = $annot->{$gene_id};
                $trans_info{$t}->{symbol} = $annot->{$symbol_id};
            }
            if (not $benign_alleles{$i} and consequenceMatchesClass($annot, $split)){
                push @{ $allele_to_csq{damaging}->{$i} }, $annot;
            }else{
                push @{ $allele_to_csq{benign}->{$i} }, $annot;
            }
        }#each CSQ
    }#each ALT
    return 0 if not keys %{$allele_to_csq{damaging}} and 
          not keys %{$allele_to_csq{benign}};#no variant in valid transcript

    splitRemainingFields($split);
    # get allele frequencies in  @samples if using -u/--max_control_allele_frequency
    my %s_allele_counts = ();
    my $allele_count;
    if ($opts{u}){
        %s_allele_counts = VcfReader::countAlleles
        (
            minGQ   => $opts{g},
            line    => $split,
            samples => \@samples,
            sample_to_columns => \%sample_to_col,
        );
        map { $allele_count += $s_allele_counts{$_} } keys %s_allele_counts;
        # filter if frequency in @controls is greater or equal to 
        # -u/--max_sample_allele_frequency
        if ($opts{u} and $allele_count > 0){
            foreach my $i (keys %{$allele_to_csq{damaging}}){
                my $freq = $s_allele_counts{$i}/$allele_count;
                if ($freq >= $opts{u}){
                    push @{$allele_to_csq{benign}->{$i}}, @{$allele_to_csq{damaging}->{$i}};
                    delete $allele_to_csq{damaging}->{$i} 
                }
            }
        }
    }
    my %samp_to_gt = VcfReader::getSampleCall
    (
          line              => $split, 
          sample_to_columns => \%sample_to_col,
          minGQ             => $opts{g},
          multiple          => \@samples,
    );
    foreach my $k (keys %allele_to_csq){
        foreach my $j (keys %{$allele_to_csq{$k}}){
            foreach my $annot (@{ $allele_to_csq{$k}->{$j} }){
                my $t = $annot->{$feature_id};
                #add sample IDs as entries to this variant count class
                # e.g. $counts{ENST00001234567}->{benign}->{sample_1} = undef
                # we will count no. keys for each transcript at end of chromosome
                my %samps_with_allele = map {$_ => 1 } addSamplesWithAllele
                (
                    \%samp_to_gt, 
                    $j,
                    $split,
                );
                map { $counts{$t}->{$k}->{$_} = undef } keys %samps_with_allele;
                #output annotation details and sample counts if --q 
                if ($VAR and $k eq 'damaging'){
                    my @cases_with_allele = grep { $samps_with_allele{$_} } @cases;
                    my @hom_cases = grep { $samp_to_gt{$_} =~ /^$j[\|\/]$j$/ } @cases_with_allele;
                    my @het_cases = grep { $samp_to_gt{$_} !~ /^$j[\|\/]$j$/ } @cases_with_allele;
                    my @conts_with_allele = grep { $samps_with_allele{$_} } @controls;
                    my @hom_conts = grep { $samp_to_gt{$_} =~ /^$j[\|\/]$j$/ } @conts_with_allele;
                    my @het_conts = grep { $samp_to_gt{$_} !~ /^$j[\|\/]$j$/ } @conts_with_allele;
                    print $VAR join
                    (
                        "\t", 
                        @$split[0..1],
                        $alleles[0],
                        $alleles[$j],
                        $annot->{$symbol_id},
                        $annot->{$gene_id} ,
                        $t,
                        $annot->{$hgvsc} || '.',
                        $annot->{$hgvsp} || '.',
                        scalar(@cases_with_allele),
                        scalar(@conts_with_allele),
                        join(",", @cases_with_allele) || '.',
                        join(",", @conts_with_allele) || '.',
                        scalar(@het_cases),
                        scalar(@het_conts),
                        join(",", @het_cases) || '.',
                        join(",", @het_conts) || '.',
                        scalar(@hom_cases),
                        scalar(@hom_conts),
                        join(",", @hom_cases) || '.' ,
                        join(",", @hom_conts) || '.' ,
                    ) . "\n";
                }
            }
        }
    }
    return keys %{$allele_to_csq{damaging} } ; 
    #returns a count on no. 'damaging' alleles found in line
}

#################################################
sub outputCounts{
    #outputs counts for each transcript in %counts
    #and clears %counts hash to save memory
    foreach my $t (keys %counts){
        my %case_counts = (damaging => 0, benign => 0);
        my %cont_counts = (damaging => 0, benign => 0);
        my %fisher = (
            damaging => { 'p.value' => "NA", 'conf.int' => "NA", estimate => "NA"},
            benign   => { 'p.value' => "NA", 'conf.int' => "NA", estimate => "NA"},
        );
        foreach my $d (qw /damaging benign/){
            if (exists $counts{$t}->{$d} ){
                map { $case_counts{$d}++ if exists $counts{$t}->{$d}->{$_} } @cases;
                map { $cont_counts{$d}++ if exists $counts{$t}->{$d}->{$_} } @controls;
            }
            my $d_counts  = join
            (   ",", 
                $case_counts{$d},
                @cases - $case_counts{$d},
                $cont_counts{$d},
                @controls - $cont_counts{$d},
            );
            my $d_matrix = "matrix(c($d_counts), nrow = 2, ncol = 2)";
            $R->send
            (
                "fresult <- fisher.test($d_matrix, alternative='greater')"
            );
            foreach my $k (keys %{$fisher{$d}}){
                $R->send("cat(fresult\$$k)"); 
                $fisher{$d}->{$k} = $R->read();
            }
            my @ci = split(/\s+/, $fisher{$d}->{'conf.int'});
            $fisher{$d}->{l95} = defined $ci[0] ? $ci[0] : "NA";
            $fisher{$d}->{u95} = defined $ci[1] ? $ci[1] : "NA";
        }
        
        print $OUT join
        (
            "\t", 
            $trans_info{$t}->{symbol},
            $trans_info{$t}->{gene},
            $t,
            $case_counts{damaging},
            @cases - $case_counts{damaging},
            $cont_counts{damaging},
            @controls - $cont_counts{damaging},
            $case_counts{benign},
            @cases - $case_counts{benign},
            $cont_counts{benign},
            @controls - $cont_counts{benign},
            $fisher{damaging}->{'p.value'},
            $fisher{damaging}->{'estimate'},
            $fisher{damaging}->{l95},
            $fisher{damaging}->{u95},
            $fisher{benign}->{'p.value'},
            $fisher{benign}->{'estimate'},
            $fisher{benign}->{l95},
            $fisher{benign}->{u95},
        ) . "\n";
    }
    %counts = ();
}

#################################################
sub addSamplesWithAllele{
    my ($gts, $allele, $line) = @_; 
    my @samp = ();
    foreach my $s ( keys %{$gts} ){
        if ($gts->{$s} =~ /^($allele[\/\|]\d+|\d+[\/\|]$allele)$/){
            if (checkAllelePl($allele, $s, $line)){
                push @samp, $s;
            }
        }
    }
    return @samp;
}

#################################################
sub checkAllelePl{
    #returns 0 if PL for any genotype that does not include $allele
    # is lower than $opts{pl}. Otherwise returns 1.
    return 1 if not $opts{pl} ;
    my ($allele, $sample, $line) = @_;
    my $pl = VcfReader::getSampleGenotypeField
    (
        line => $line, 
        field => "PL",
        sample => $sample, 
        sample_to_columns => \%sample_to_col,
    );
    my @pls = split(",", $pl);
    my @idxs = VcfReader::calculateOtherAlleleGindexes($line, $allele); 
    foreach my $i (@idxs){
        if ($i > $#pls){
            my ($chrom, $pos) = VcfReader::getMultipleVariantFields
            (
                $line,
                'CHROM', 
                'POS'
            );
            informUser
            ( 
                "WARNING: Not enough PL scores ($pl) for allele " .
                "'$i', sample $sample at $chrom:$pos. Your VCF may be malformed.\n"
            );
            return 1;
        }
        return 0 if $pls[$i] < $opts{pl};
    }
    return 1;
}

#################################################
sub checkGtPl{
    #returns 0 if PL for any genotype other than those in @$gts
    # is lower than $opts{pl}. Otherwise returns 1.
    return 1 if not $opts{pl} ;
    my ($gts, $sample, $line) = @_;
    my $pl = VcfReader::getSampleGenotypeField
    (
        line => $line, 
        field => "PL",
        sample => $sample, 
        sample_to_columns => \%sample_to_col,
    );
    my @pls = split(",", $pl); 
    my @idxs = ();
    foreach my $gt (@$gts){
        push @idxs,  VcfReader::calculateGenotypeGindex($gt);
    }
    foreach my $idx (@idxs){
        if ($idx > $#pls){
            my ($chrom, $pos) = VcfReader::getMultipleVariantFields
            (
                $line,
                'CHROM', 
                'POS'
            );
            informUser
            ( 
                "WARNING: Not enough PL scores ($pl) for genotypes '".
                join(", ", @$gts) . "', sample $sample at $chrom:$pos.".
                " Your VCF may be malformed.\n"
            );
            return 1;
        }
    }

    for (my $i = 0; $i < @pls; $i ++){
        next if grep {$_ == $i} @idxs;
        return 0 if $pls[$i] < $opts{pl};
    }
    return 1;
}

#################################################
sub getAndCheckCsqHeader{
    my %csq_head = ();
    if (not $opts{m}){
        eval { 
            %csq_head = VcfReader::readVepHeader
            (
                header => \@header
            ); 
        } ;
        if (not $@){
            $opts{m} = 'vep';
            informUser("Found VEP header - using VEP 'CSQ' annotations.\n");
        }else{
            informUser("No VEP header found. Trying SnpEff...\n");
            eval { 
                %csq_head = VcfReader::readSnpEffHeader
                (
                    header => \@header
                ); 
            } ;
            if (not $@){
                $opts{m} = 'snpeff';
                informUser("Found SnpEff header - using SnpEff 'ANN' annotations.\n");
            }else{
                die "ERROR: Could not find VEP or SnpEff headers in input. ".
                    "Please annotate your input with either program and try again.\n";
            }
        }
    }else{
        if ($opts{m} eq 'vep'){
            %csq_head = VcfReader::readVepHeader
                (
                    header => \@header
                ); 
        }else{
            %csq_head = VcfReader::readSnpEffHeader
            (
                header => \@header
            ); 
        }
    }
    return %csq_head;
}

#################################################
sub getCsqFields{
   my @fields = ();
    if ($opts{m} eq 'vep'){
       @fields =  
        qw(
            allele
            gene
            feature
            feature_type
            consequence
            symbol
            biotype
            hgvsc
            hgvsp
        );
        if ($opts{canonical_only}){
            push @fields, "canonical";
        }
        if ($opts{consensus_splice_site}){
            push @fields, "splice_consensus";
        }
    }else{
        @fields = 
        qw(
            allele
            annotation
            annotation_impact
            gene_name
            gene_id
            feature_type
            feature_id
            transcript_biotype
            hgvs.c
            hgvs.p
        );
        if ($opts{canonical_only}){
            informUser
            ( 
                "WARNING: --canonical_only option is ignored when ".
                "working with SnpEff annotations.\n"
            );
            delete $opts{canonical_only};
        }
    }
    foreach my $f (@fields){
        if (not exists $csq_header{$f}){
            if ($f eq 'splice_consensus'){
                die "Could not find 'splice_consensus' field in header " . 
                  "- please ensure you use the SpliceConsensus plugin when ".
                  "running VEP.\n";
            }
            die "Could not find '$f' field in $opts{m} consequence header " .
              "- please ensure you have annotated your file including the appropriate ".
              "fields.\n";
        }
    }
    return @fields;
}

#################################################
sub getAndCheckClasses{
    my %all_classes = ();
    if ($opts{m} eq 'vep'){
        %all_classes =  VcfhacksUtils::readVepClassesFile();
    }else{
        %all_classes =  VcfhacksUtils::readSnpEffClassesFile();
    }
    grep { $all_classes{$_} eq 'default' } keys %all_classes;
    %default_classes = 
      map  { $_ => undef } 
      grep { $all_classes{$_} eq 'default' } 
      keys %all_classes;

    if (not @classes){
        @classes = keys %default_classes;
    }
    push @classes, @add_classes if (@add_classes);
    if ($opts{m} eq 'vep' and $opts{consensus_splice_site}){
        push @classes, "splice_region_variant" ;
    }
    @classes = map { lc($_) } @classes; 
    @classes = VcfhacksUtils::removeDups(@classes);
    foreach my $class (@classes) {
        die "Error - variant class '$class' not recognised.\n"
          if not exists $all_classes{lc($class)} ;
    }
    return @classes;
}

#################################################
sub getAndCheckBiotypes{
    return if $opts{no_biotype_filtering}; 
    my %all_biotypes = VcfhacksUtils::readBiotypesFile();
    if (not @biotypes){
        @biotypes = ('protein_coding');
    }
    @biotypes = map { lc($_) } @biotypes; 
    @biotypes = VcfhacksUtils::removeDups(@biotypes);
    foreach my $biotyp (@biotypes) {
        die "Error - biotype '$biotyp' not recognised.\n"
          if not exists $all_biotypes{lc($biotyp)} ;
    }
    return @biotypes;
}
    
#################################################
sub getAndCheckInSilicoPred{
    return if not @damaging;
    my %filters = VcfhacksUtils::getAndCheckInSilicoPred($opts{m}, \@damaging);
    if ($opts{m} eq 'vep'){#VEP prediction results will be in CSQ field
        foreach my $d (@damaging){
            if ($d ne 'all' and not exists $csq_header{$d}){
                informUser
                (
                    "WARNING: No '$d' field found in CSQ header of VCF. ".
                    "No in silico filtering will be performed for $d.\n"
                ); 
            }
        }
        %filters = map  { $_ => $filters{$_} } 
                   grep { exists $csq_header{$_} } 
                   keys %filters;
        push @csq_fields, keys %filters;
    }#SnpEff predictions will be added via SnpSift
    return %filters;
}

#################################################
sub getAndCheckEvalFilters{
    return if not @eval_filters;
    my @filters = (); 
    my %csq_add = ();
FLT: foreach my $s (@eval_filters){
        my %f =  VcfhacksUtils::getEvalFilter($s);
        foreach my $fld (@{$f{field}}){
            (my $fb = $fld) =~ s/^\(//;#we may have used ( to specify precedence
            if (not exists $csq_header{lc($fb)}){
                informUser
                (
                    "WARNING: No '$fb' field found in CSQ header of ".
                    "VCF. Cannot use --eval_filter expression '$s' ".
                    "for filtering.\n"
                ); 
                next FLT;
            }
            $csq_add{lc($fb)}++;
        }
        push @filters, \%f;
    }
    push @csq_fields, keys %csq_add;
    return @filters;
}

#################################################
sub splitRemainingFields{
    my $s = shift;#array ref to modify in place
    my $l = pop(@$s);
    push @$s, split("\t", $l);
}

#################################################
sub openOutput{
    my $OUT_FH;
    my $VAR_FH;
    if ($opts{o}) {
        open( $OUT_FH, ">", $opts{o} ) || die "Can't open $opts{o} for writing: $!\n";
    }
    else {
        $OUT_FH = \*STDOUT;
    }
    if ($opts{s}){
        open( $VAR_FH, ">", $opts{s} ) || die "Can't open $opts{s} for writing: $!\n";
    }
    return ($OUT_FH, $VAR_FH);
}

#################################################
sub writeOptionsToHeader{
    #print meta header lines
    #add header line detailing program options
    my $head_opts = VcfhacksUtils::getOptsVcfHeader(%opts) . "\n"; 
    print $OUT $head_opts;
    print $OUT join
    (
        "\t",
        qw/ 
            Symbol
            Gene
            Transcript
            CasesWithDamagingVariant
            CasesWithoutDamagingVariant
            ControlsWithDamagingVariant
            ControlsWithoutDamagingVariant
            CasesWithBenignVariant
            CasesWithoutBenignVariant
            ControlsWithBenignVariant
            ControlsWithoutBenignVariant
            Damaging_P_Value
            Damaging_OR
            Damaging_l95
            Damaging_u95
            Benign_P_Value
            Benign_OR
            Benign_l95
            Benign_u95
        /
    ) , "\n";
    if ($VAR){
        print $VAR $head_opts ;
        print $VAR join
        (
            "\t", 
            qw /
                CHROM
                POS
                REF
                ALT
                SYMBOL
                GENE
                TRANSCRIPT
                HGVSC
                HGVSP
                N_CASES
                N_CONTROLS
                CASE_IDS
                CONTROL_IDS
                N_HET_CASES
                N_HET_CONTROLS
                HET_CASE_IDS
                HET_CONTROL_IDS
                N_HOM_CASES
                N_HOM_CONTROLS
                HOM_CASE_IDS
                HOM_CONTROL_IDS
            /
        ) . "\n";
    }
}

#################################################
sub informUser{
    my $msg = shift;
    my $time = strftime( "%H:%M:%S", localtime );
    if ($progressbar){
        $progressbar->message( "[INFO - $time] $msg" );
    }else{
        print STDERR "[INFO - $time] $msg";
    }
}
#################################################
sub updateProgressBar{
    if ($progressbar) {
        if ($var_count >= $next_update){
            $next_update = $progressbar->update( $var_count )
        }
    }elsif($opts{b}){
        VcfhacksUtils::simpleProgress($var_count);
    }
}

#################################################
sub haveVariant{
    return 1 if not @samples;
    my $gts = shift;
    foreach my $k (keys %$gts){
        if ($gts->{$k} =~ /(\d+)[\/\|](\d+)/){
            return 1 if ( $1 > 0 or $2 > 0 );
        }
    }
    return 0;
}

#################################################
sub getAfAnnotations{
    my $info_fields = shift;
    my %af_found = ();
    my @af_fields =  qw ( 
        AS_CAF
        AS_G5A
        AS_G5
        AS_COMMON
        EVS_EA_AF
        EVS_AA_AF
        EVS_ALL_AF
    );
    foreach my $c (@custom_af){
        if (not exists $info_fields->{$c}){
            informUser( "WARNING: User specified custom allele frequency ".
                "(-j/--custom_af) field '$c' not found in VCF header. This ".
                "field may not exist in your VCF.\n");
        }
    }
    push @af_fields, @custom_af;
    foreach my $key (keys %$info_fields){
        my $warning = "WARNING: Found expected frequency annotation ($key) in ". 
          "INFO fields, but 'Number' field is $info_fields->{$key}->{Number}, ".
          "expected 'A'. Ignoring this field.\n";
        my $info = "Found allele frequency annotation: $key. ".
          "This will be used for filtering on allele frequency.\n";
        if (grep { $key eq $_ } @af_fields){
            if ($info_fields->{$key}->{Number} ne 'A'){
                informUser($warning);
            }else{
                informUser($info);
                $af_found{$key} = $info_fields->{$key};
            }
        }else{
            if ($key =~ /^FVOV_AF_\S+$/){
                if ($info_fields->{$key}->{Number} ne 'A'){
                    informUser($warning);
                }else{
                    informUser($info);
                    $af_found{$key} = $info_fields->{$key};
                }
            }
        }
    }
    return %af_found;
}

#################################################
sub alleleAboveMaf{
    my $i = shift; #1-based index of alt allele to assess
    my $af_values = shift; #hash ref of INFO fields to their values
    foreach my $k (keys %{$af_values}){
        next if not $af_values->{$k};
        next if $af_values->{$k} eq '.';
        my @split = split(",", $af_values->{$k}); 
        next if $split[$i] eq '.';
        if ($k eq "AS_G5" or $k eq "AS_G5A"){
            if ($opts{a} <= 0.05 and $split[$i] > 0){
                return 1;
            }
        }elsif ($k eq "AS_COMMON"){
            if ($opts{a} <= 0.01 and $split[$i] > 0){
                return 1;
            }
        }else{#should be a standard allele freq now
            if ($af_info_fields{$k}->{Type} eq 'Float' or $k eq 'AS_CAF'){
                return 1 if $split[$i] >= $opts{a};
            }else{
                inform_user("WARNING: Don't know how to parse INFO field: $k.\n");
            }
        }
    }
    return 0;
}
#################################################
sub getAfInfoValues{
    my $l = shift;
    my %values = ();
    foreach my $k (keys %af_info_fields){
        $values{$k} = VcfReader::getVariantInfoField($l, $k);
    }
    return %values;
}
#################################################
sub getConsequences{
    my $line = shift; 
    if ($opts{m} eq 'vep'){
        return VcfReader::getVepFields
        ( 
            line        => $line,
            field       => \@csq_fields,
            vep_header  => \%csq_header,
        );
    }else{
        return VcfReader::getSnpEffFields
        ( 
            line          => $line,
            field         => \@csq_fields,
            snpeff_header => \%csq_header,
        );
    }
}

#################################################
sub consequenceMatchesClass{
    my $annot = shift;
    my $l = shift;
    my $use_default_classes = shift;
    if ($opts{m} eq 'vep'){
        return consequenceMatchesVepClass($annot, $use_default_classes);
    }else{
        return consequenceMatchesSnpEffClass($annot, $l, $use_default_classes);
    }   
}

#################################################
sub consequenceMatchesVepClass{
    my $annot = shift;
    my $use_default = shift;
    #intergenic variants have no feature associated with them - skip
    return 0 if $annot->{consequence} eq "intergenic_variant";
    #skip unwanted biotypes
    return 0 if (not exists $biotype_keep{$annot->{biotype}}) ;
    #skip non-canonical transcripts if --canonical_only selected
    if ($opts{canonical_only}) {
        return 0 if ( not $annot->{canonical} );
    }
    
    my @anno_csq = split( /\&/, $annot->{consequence} );
    #skip NMD transcripts
    return 0 if ( grep { /NMD_transcript_variant/i } @anno_csq );

    #eval filters trump annotation class
    foreach my $evf (@eval_exp){
        return 1 if VcfhacksUtils::evalFilter($evf, $annot);
    }
    #score filters trump annotation class
        
ANNO: foreach my $ac (@anno_csq){
        $ac = lc($ac);#we've already converted %class_filters to all lowercase
        if ($use_default){
            return exists $default_classes{$ac} ;
        }
        if ( exists $class_filters{$ac} ){
            if ($ac eq 'missense_variant' and %in_silico_filters){
                return 1 if damagingMissenseVep($annot); 
            }elsif ( lc $ac eq 'splice_region_variant' 
                     and $opts{consensus_splice_site}){
                my $consensus = $annot->{splice_consensus};
                next if not $consensus;
                if ( $consensus !~ /SPLICE_CONSENSUS\S+/i ) {
                    inform_user(
"WARNING: SPLICE_CONSENSUS annotation '$consensus' is not " .
"recognised as an annotation from the SpliceConsensus VEP plugin.\n");
                }else{
                    return 1;
                }
            }else{
                return 1;
            }
        }
    }
    return 0;#no annotation matching %class_filters
}

#################################################
sub consequenceMatchesSnpEffClass{
    my $annot = shift;
    my $l = shift;
    my $use_default = shift;
    #skip variants with undef features (intergenic variants)
    return if not defined $annot->{feature_id};
    #skip unwanted biotypes
    return 0 if not exists $biotype_keep{lc $annot->{transcript_biotype} };

    #eval filters trump annotation class
    foreach my $evf (@eval_exp){
        return 1 if VcfhacksUtils::evalFilter($evf, $annot);
    }

    my @anno_csq = split( /\&/, $annot->{annotation} );
ANNO: foreach my $ac (@anno_csq){
        $ac = lc($ac);#we've already converted %class_filters to all lowercase
        if ($use_default){
            return exists $default_classes{$ac} ;
        }
        if ( exists $class_filters{$ac} ){
            if ($ac eq 'missense_variant' and %in_silico_filters){
                return 1 if damagingMissenseSnpEff($l); 
            }else{
                return 1;
            }
        }
    }
    return 0;#no annotation matching %class_filters
}

#################################################
sub damagingMissenseVep{
    #returns 1 if variant is damaging and should be kept
    #returns 0  if variant is benign and should be filtered
    my $anno = shift; 
    my %filter_matched = ();
PROG: foreach my $k ( sort keys %in_silico_filters) {
        my $score = $anno->{ lc $k };
        if ( not $score or $score =~ /^unknown/i ){ 
        #don't filter if score is not available for this variant 
        #unless skip_unpredicted is in effect
            $filter_matched{$k}++ unless $opts{skip_unpredicted};
            next;
        }
SCORE:  foreach my $f ( @{ $in_silico_filters{$k} } ) {
            my $do_filter = 0;
            if ($f =~ /^(polyphen|sift|condel)$/){
                $do_filter = isDamagingVepInsilico($k, $f, $score);
            }else{
                $do_filter = isDamagingDbnsfpInsilico($k, $f, $score);
            }
            if ($do_filter){
                return 1 if $opts{keep_any_damaging};
                $filter_matched{$k}++;
            }
        }
    }
    foreach my $k ( keys %in_silico_filters ) {
        #filter if any of sift/condel/polyphen haven't matched our deleterious settings
        return 0 if not exists $filter_matched{$k};
    }
    return 1;
}

#################################################
sub isDamagingDbnsfpInsilico{
    my ($tool, $f, $score) = @_;
    my @pred =  split("&", $score); #multiple scores are separated by '&'
    if ( $f =~ /^\d+(\.\d+)*$/ ) {#if we want to filter on score
        if ($tool =~ /^(fathmm_score|provean_score|sift_score)$/){
            #lower is more damaging for these tools
            my $min = min(@pred); 
            return $min <= $f;
        }else{
            #higher = more damaging
            my $max = max(@pred); 
            return $max >= $f;
        }
    }else{#if we want to filter on prediction
        return grep { lc($_) eq lc($f) } @pred;
    }
}

#################################################
sub isDamagingVepInsilico{
    my ($tool, $f, $score) = @_;
    if ( $f =~ /^\d(\.\d+)*$/ ) {#if we want to filter on score
        my $prob;
        if ( $score =~ /\((\d(\.\d+)*)\)/ ) {#get score
            $prob = $1;
        }else {
            return 0;
        }
        if ( lc $tool eq 'polyphen'){
            return $prob >= $f;
            #higher is more damaging for polyphen 
        }else{
            return $prob <= $f;
            #lower is more damaging for sift and condel 
        }
    }else{#if we want to filter on prediction
        $score =~ s/\(.*\)//;#get prediction
        return  lc $f eq lc $score ; #damaging if match
    }
}



#################################################
sub damagingMissenseSnpEff{
    #returns 1 if variant is damaging and should be kept
    #returns 0  if variant is benign and should be filtered
    my $l = shift; 
    my %filter_matched = ();
PROG: foreach my $k ( sort keys %in_silico_filters) {
        my $pred = VcfReader::getVariantInfoField($l, $k);
        if (not defined $pred){
            $filter_matched{$k}++ unless $opts{skip_unpredicted};
            next;
        }
        my @scores = split(",", $pred); 
         #snpsift does not annotate properly per allele...
         #... so as a fudge we have to just count any matching value
        foreach my $score (@scores){
            next if ( not $score or $score eq '.' );
            foreach my $f ( @{ $in_silico_filters{$k} } ) {
                if ( lc $f eq lc $score ) {    #damaging
                    return 1 if $opts{k};
                    $filter_matched{$k}++;
                    next PROG;
                }
            }
        }
    }
    foreach my $k ( keys %in_silico_filters ) {
        #filter if any of sift/condel/polyphen haven't matched our deleterious settings
        return 0 if not exists $filter_matched{$k};
    }
    return 1;
}


#################################################

=head1 NAME

geneBurdenTest.pl - identify damaging variants in genes and perform a gene burden test in cases vs controls 

=head1 SYNOPSIS

    geneBurdenTest.pl -i <variants.vcf> -f <cases_and_controls.ped> [options]
    geneBurdenTest.pl --help (show help message)
    geneBurdenTest.pl --manual (show manual page)

=cut 

=head1 ARGUMENTS

=over 8 

=item B<-i    --input>

VCF file annotated with Ensembl's variant_effect_predictor.pl (VEP) script or SnpEff.

=item B<-f    --ped_file>

PED file giving sample IDs of cases and controls. 

PED format - from the PLINK documentation:

       The PED file is a white-space (space or tab) delimited file: the first six columns are mandatory:

            Family ID
            Individual ID
            Paternal ID
            Maternal ID
            Sex (1=male; 2=female; other=unknown)
            Phenotype

       Affection status, by default, should be coded:

           -9 missing
            0 missing
            1 unaffected
            2 affected

In this manner, only samples present in this file and with an affectation status of '2' will be used as 'cases'. Only samples present in this file and with an affectation status of '1' will be used as 'controls'.


=item B<-o    --output>

File to print output table (optional). Will print to STDOUT by default.

=item B<-s    --var_summaries>

Optional file to print details of each 'damaging' variant used in the gene burden analyses.

=item B<-m    --mode>

This program will attempt to detect the format of your input automatically by looking for VEP or SnpEff annotations, but you may specify either 'vep' or 'snpeff' with this option to select the mode employed by the script if you have a file with both VEP and SnpEff annotations. By default, if both annotations are present and the program is run without this option, VEP annotations will be used. NOTE: Only SnpEff annotations using the more recent 'ANN' style annotations rather than the older 'EFF' style are recognised by this program.

=item B<--classes>

One or more mutation classes to retrieve. By default only variants labelled with one of the following classes will count as 'functional' variants:

B<VEP:>

        TFBS_ablation
        TFBS_amplification
        frameshift_variant
        inframe_deletion
        inframe_insertion
        initiator_codon_variant
        missense_variant
        protein_altering_variant
        regulatory_region_ablation
        regulatory_region_amplification
        splice_acceptor_variant
        splice_donor_variant
        stop_gained
        stop_lost
        transcript_ablation
        transcript_amplification

B<SnpEff:>

        chromosome
        coding_sequence_variant
        inframe_insertion
        disruptive_inframe_insertion
        inframe_deletion
        disruptive_inframe_deletion
        exon_loss_variant
        frameshift_variant
        missense_variant
        rare_amino_acid_variant
        splice_acceptor_variant
        splice_donor_variant
        stop_lost
        5_prime_UTR_premature_start_codon_gain_variant
        start_lost
        stop_gained
        exon_loss

Available classes that can be chosen instead of (or in addition to - see below) these classes can be found in the data/vep_classes.tsv and data/snpeff_classes.tsv files respectively. 

=item B<--add_classes>

Specify one or more classes, separated by spaces, to add to the default mutation classes used for finding functional variants.

=item B<--consensus_splice_site>

Use this flag in order to keep splice_region_variant classes only if they are in a splice consensus region as defined by the SpliceConsensus VEP plugin. You do not need to specify 'splice_region_variant' using --classes or --add_classes options when using this flag. You B<MUST> have used the SpliceConsensus plugin when running the VEP for this option to work correctly. This option is only used when running on VEP annotations. 

=item B<--canonical_only>

Only consider canonical transcripts (VEP annotations only). 

=item B<-c    --cadd_filter>

If you have annotated CADD phred scores for alleles using rankOnCaddScore.pl you may use this option to specify a CADD phred score threshold for variants. Any alleles that have a CADD phred score below the value specified here will be filtered. Variants without a CADD phred score will not be filtered unless using the --skip_unpredicted option. You should have run rankOnCaddScore.pl with the --do_not_sort option to maintain chromosome order of your variants if you use this option with the -f/--find_shared_genes option or else you will need to sort your VCF before running getFunctionalVariants.pl.

=item B<-d    --damaging>

Specify in silico prediction scores to filter on. 

If running on VEP annotations you may either use PolyPhen, SIFT and/or Condel scores given by VEP or annotations from dbNSFP added using the dbNSFP VEP plugin. If running on SnpEff, annotations scores provided by SnpSift's 'dbnsfp' mode will be used.

NOTE: when using SnpEff annotations, prediction program names are case-sensitive.

=over 12

=item B<VEP mode:> 

Specify SIFT, PolyPhen or Condel labels or scores to filter on. Add the names of the programs you want to use, separated by spaces, after the --damaging option. By default SIFT will keep variants labelled as 'deleterious', Polyphen will keep variants labelled as 'possibly_damaging' or 'probably_damaging' and  Condel will keep variants labelled as 'deleterious'.

If you want to filter on custom values specify values after each program name in the like so: 

    'polyphen=probably_damaging' 

Seperate multiple values with commas - e.g. 

    'polyphen=probably_damaging,possibly_damaging,unknown' 

You may specify scores between 0 and 1 to filter on scores rather than labels - e.g. 

    'sift=0.3' 

For polyphen, variants with scores lower than this score are considered benign and filtered, for SIFT and Condel higher scores are considered benign.

B<Valid labels for SIFT:> deleterious, tolerated

B<Valid labels for Polyphen:> probably_damaging, possibly_damaging, benign, unknown

B<Valid labels for Condel:> deleterious, neutral

To use default values for all three programs you may use 'all' (i.e. '--damaging all') BUT PLEASE NOTE: if you have added dbNSFP annotations to your input VCF these will also be included (see below).

=item B<SnpEff or VEP dbNSFP mode:> 

Your input must have been annotated using SnpSift's dbnsfp option (if using SnpEff annotations) or using the dbNSFP VEP plugin (if using VEP annotations). Recognised annotations are:

    fathmm-MKL_coding_pred
    FATHMM_pred
    LRT_pred
    MetaLR_pred
    MetaSVM_pred
    MutationAssessor_pred
    MutationTaster_pred
    PROVEAN_pred
    Polyphen2_HDIV_pred
    Polyphen2_HVAR_pred
    SIFT_pred

You may instead choose to use the 'score' or dbNSFP 'rankscore' annotations for these tools (e.g. 'fathmm-MKL_coding_score'). 

You may specify 'all' to use the default values for all programs. Choosing custom values is also performed in the same way as for VEP annotations above (e.g. dbNSFP_MutationAssessor_pred=H,M,L). 

=back

For more details of the available and default settings for these programs please see the files 'data/snpeff_insilico_pred.tsv' or 'data/vep_insilico_pred.tsv'.

The default behaviour is to only keep variants predicted as damaging by ALL programs specified, although if the value is not available for one or more programs than that program will be ignored for filtering purposes. See the next two options for alternative behaviours.

=item B<-k    --keep_any_damaging>

If using multiple programs for filters for --damaging argument use this flag to keep variants predicted to be damaging according to ANY of these programs.

=item B<--skip_unpredicted>

Skip alleles that do not have a score from one or more programs specified by the -d/--damaging argument or if using the -c/--cadd_filter option and an allele has no CADD phred score. The --keep_any_damaging argument will override this behaviour for -d/--damaging predictions if any of the available predictions are considered damaging.

=item B<--biotypes>

By default only features/transcripts with 'protein_coding' biotype are 
considered.

You may override this filter by specifying  list of biotypes here (separated 
with whitespace). 

The 'data/biotypes.tsv' file contains a list of valid biotypes.

=item B<--no_biotype_filtering>

Use this flag to consider consequences affecting ALL biotypes.

=item B<-a    --af>

Use a value between 0.00 and 1.00 to specify allele frequencey filtering for annotations from annotateSnps.pl, filterOnEvsMaf.pl or filterVcfOnVcf.pl if you've previously run these programs to annotate your VCF. If an allele frequency is available for an allele it will be filtered if equal to or greater than the value specfied here. 

Note: allele frequencies added by VEP are not used for filtering as they check the allele frequency at the site, not of the specific alleles in your variant call.

=item B<-j    --custom_af>

If using the --af/--allele_frequency option and your data contains allele frequency fields from sources not recognised by this program, you may give the name of these allele frequency INFO fields here and they will be used for filtering in addition to the default fields. Note that these annotations must contain an annotation per ALT allele (i.e. the INFO field header must specify 'Number=A') to work properly and the allele frequency should be expressed as a number between 0.00 and 1.00 in order to be compatible with the default allele frequency fields recognised by this program.

=item B<-u  --max_sample_allele_frequency>

Use this option to specify an allele frequency (between 0.00 and 1.00) for filtering alleles in your VCF. Alleles present at this frequency or higher in your case and control samples (combined) will be filtered. 

=item B<-v    --var_qual>

Minimum variant Phred-like quality score to consider. Variants with a QUAL field lower than this value will be filtered. Default is 20.

=item B<-g    --gq>

Minimum genotype qualities to consider. Only applicable when using the -s/--samples option. Any genotype call below this threshold will be considered a no call. Default is 20

=item B<--pl>

Minimum 0-based phred-scale genotype likelihood (see the VCF spec for details) for alternative genotypes. Only applicable when using the -s/--samples option. When considering a given allele, if the sample has a PL below this value for a genotype not including this allele, the sample will not be considered a carrier of that allele. Default - not used.

=item B<--pass_filters>

Only consider variants with a PASS filter field. If the FILTER field for variant is not PASS the variant will be skipped.

=item B<--eval_filters>

Use this option to create custom filters for KEEPING variants on the basis of values in the VEP or SnpEff consequence fields. 

Expressions must take the format of 'field name' 'comparator' 'value to compare' separated by white space. Multiple expressions can be used together along with the logical operators 'and', 'or' or 'xor'. The value for 'field name' will be used to extract the value for the given field from the VEP/SnpEff consequence INFO field. The resulting expression is evaluated using perl's built-in 'eval' function.

For example:

    --eval_filters "LoF eq 'HC'" 
    #keeps any variant with a LoF annotation of 'HC'

    --eval_filters "(ada_score >= 0.6 and rf_score >= 0.6) or maxentscan_diff > 5"
    #keeps variants with ada_score and rf_scores of 0.6 or higher or with 
    #maxenstscan diff of 5 or higher

=item B<-b    --progress>

Show a progress bar while working.

=item B<-h    --help>

Show the program's help message.

=item B<--manual>

Show the program's manual page.

=back

=cut

=head1 EXAMPLES

 Perform a case/control test using default options:

    geneBurdenText.pl -i input.vcf \
    -o out.tsv \
    -f cases_controls.ped

Exclude missense variants predicted as 'Tolerated' by FATHMM and with a minor allele
frequency of 0.1 \% or higher in dbSNP or ExAC. The former requires that you annotated
your VCF with FATHMM annotations using dbNSFP when running VEP or SnpSift. The latter
requires that you annotated frequencies using annotateSnps.pl and/or [for ExAC]
filterVcfOnVcf.pl

    geneBurdenText.pl -i input.vcf \
    -o out.tsv \
    -f cases_controls.ped \
    -d FATHMM_pred \
    --af 0.001

As above but specifying a subset of qualifying consequence classes:

 geneBurdenText.pl -i input.vcf \
 -o out.tsv \
 -f cases_controls.ped \
 -d FATHMM_pred \
 --af 0.001 \
 --classes frameshift_variant missense_variant splice_acceptor_variant splice_donor_variant stop_lost stop_gained

=cut

=head1 DESCRIPTION

This program reads functional consequence annotations in a VCF and performs a case/control burden test for qualifying 'damaging' variants in transcripts using a one-tailed Fisher's exact test. A p-value for every qualifying transcript is provided in the output as well as counts for cases/controls with 'damaging' and 'benign' variants.

=cut

=head1 AUTHOR

David A. Parry

University of Edinburgh

=head1 COPYRIGHT AND LICENSE

Copyright 2017, David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut



