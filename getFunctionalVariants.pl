#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use POSIX qw/strftime/;
use Data::Dumper;
use List::Util qw ( first max min ) ;
use File::Temp qw/ tempfile /;
use FindBin qw($RealBin);
use lib "$RealBin/lib";
use VcfReader;
use VcfhacksUtils;

my $progressbar;
my @samples = ();
my @classes = (); 
my @add_classes = ();
my @damaging = ();
my @biotypes = ();
my @custom_af = ();
my @score_filters = (); 

my %opts = 
(
    classes          => \@classes,
    add_classes      => \@add_classes,
    d                => \@damaging,
    s                => \@samples,
    biotype_filters  => \@biotypes,
    j                => \@custom_af,
    score_filters    => \@score_filters,
);

GetOptions(
    \%opts,
    "add_classes=s{,}" ,
    "a|af=f" ,
    "biotype_filters=s{,}",
    "b|progress" ,
    "canonical_only" ,
    "classes=s{,}" ,
    'clinvar=s',
    "consensus_splice_site",
    "c|cadd_filter=f",
    "d|damaging=s{,}" ,
    "e|equal_genotypes",
    "f|find_shared_genes:s", 
    "g|gq=i",
    "gene_counts=s",
    "h|?|help" ,
    "i|input=s" ,
    "j|custom_af=s{,}",
    "k|keep_any_damaging" ,
    "maf=f" ,
    "manual" ,
    "m|mode=s" ,
    "no_biotype_filtering",
    "n|num_matching=i",
    "o|output=s" ,
    "pass" ,
    "pl=f", 
    "s|samples=s{,}",
    "score_filters=s{,}",
    "skip_unpredicted" ,
    "t|target_genes=s",
    "v|var_quality=i",
    "u|max_sample_allele_frequency=f",
) or pod2usage(-message => "Syntax error", -exitval => 2);
pod2usage(-verbose => 2) if ($opts{manual});
pod2usage(-verbose => 1) if ($opts{h});

if (not $opts{i}){
    pod2usage
    (
        -message => "Syntax error - an input file must be supplied.\n", 
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
    
    if (not @samples){
        informUser
        (
            "WARNING: No 'samples specified - all samples will be used for ".
            "-u/--max_sample_allele_frequency calculation\n"
        );
#        pod2usage
#        (
#            -message =>  "-u/--max_sample_allele_frequency option requires ". 
#                         "samples to be specified with the -s/--samples ".
#                         "argument.\n", 
#            -exitval => 2
#        );
    }
}

if (defined $opts{f}){
    informUser
    (
        "WARNING: No 'samples specified - all samples will be used for ".
        "--find_shared_genes option\n"
    ) if not @samples;
#    pod2usage   
#    (
#        -message => "--find_shared_genes option requires at least one " . 
#                    "sample to be specified with the -s/--samples argument.\n", 
#        -exitval => 2
#    ) if not @samples;
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

#check --gene_counts arguments
my %gene_count_opts = (); 
if ($opts{gene_counts}){
    checkGeneCountArgs();
}

if ( ( not @samples and ($opts{f} or $opts{u} or $opts{gene_counts}) )
       or ( grep { /^all$/i } @samples )
    ){
    # if not samples specified but using --find_shared_genes, --gene_counts 
    # or --max_sample_allele_frequency option
    # or if 'all' is listed in --samples option 
    # add all samples
    if (not %sample_to_col ){
        if ( $gene_count_opts{count_mode} ne 'allele_counts'
        ){
            die "No samples found in VCF header!\n";
        }
    }
    @samples = ();#need to remove 'all' if nothing else
    push @samples, keys %sample_to_col;
}
@samples = VcfhacksUtils::removeDups(@samples);

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
    
my $feature_id = $opts{m} eq 'vep' ? "feature" : "feature_id";
my $symbol_id  = $opts{m} eq 'vep' ? "symbol"  : "gene_name";
my $gene_id    = $opts{m} eq 'vep' ? "gene"    : "gene_id";

#set default consequence fields to retrieve 
my @csq_fields = getCsqFields();#consequence fields to retrieve from VCF

#and check variant consequence classes are acceptable
my %class_filters = map { $_ => undef } getAndCheckClasses();

#check in silico prediction classes/scores are acceptable
my %in_silico_filters = getAndCheckInSilicoPred();
  #hash of prediction program names and values to filter
my @score_exp = getAndCheckScoreFilters();

#and check biotype classes are acceptable
my %biotype_filters = map { $_ => undef } getAndCheckBiotypes();


# get available allele frequency annotations if filtering on AF
# NOTE: we do not use VEP annotations as they do not necessarily match the variant allele
my %af_info_fields = (); #check available INFO fields for AF annotations
if ( defined $opts{a} ) {
    %af_info_fields = getAfAnnotations(\%info_fields);   
}

#by default keep variant if there's an associated pathogenic ClinVar field
my $keep_clinvar = checkClinVarInfo();

#if filtering on CADD score check we have a CaddPhredScore header
if ( $opts{c} ){
    if (not exists $info_fields{CaddPhredScore}){
        informUser("WARNING: No 'CaddPhredScore' INFO field found in header ".
         "- your input probably does not contain annotations for filtering on" . 
         " CADD score.\n");
    }
}

#if using -t/--target_genes argument, check file
my %targets = (); 
my %sargs   = ();
my $current_target = ''; 
my ($TMP, $tmpout);#have to write to temp file, sort and dedup if using targets
#keys are Ensembl Gene IDs, values are hashes of chrom, start, end coordinates
if ($opts{t}){
    if (not $opts{gene_counts}){
        die "-t/--target_genes option can only be used when also using "
            ."--gene_counts option\n";
    }
    readTargetGenesFile();
    %sargs = VcfReader::getSearchArguments($opts{i});
    ($TMP, $tmpout)  = tempfile(UNLINK => 1);
}

#Get filehandle for output (STDOUT or file) and optionally gene list 
#(STDERR or file) or gene counts (file)
my ($OUT, $LIST, $GENECOUNTS) = openOutput();

#write header and program run parameters
writeOptionsToHeader();

my %contigs = (); 
 #keep track of contigs we've seen so we can check input is (roughly) sorted 
 #... (only relevant if using --find_shared_genes option)
my %transcript_sample_vars = ();
 #stores variant IDs meeting our criteria per (key) transcript 
my %vcf_lines = ();
 #hash of variant IDs to matching VCF line
my %transcript_to_symbol = ();
 #transcript ID to gene symbol
my %genes_to_transcripts = ();
 #gene symbols => transcripts with variants
my %gene_burden = (); 
 #if using --gene_counts this has an entry for each gene with a functional 
 #variant. Entries are hashes where keys are sample names
my %allele_counts = ();
 #if using --gene_counts in allele_counts mode
my %allele_number = ();
 #if using --gene_counts in allele_counts mode

#progress bar vars
my $next_update      = 0;
my $total_vars       = 0; 
my $var_count        = 0;
my $functional_count = 0;
if ($opts{b}) {
    $progressbar = getProgressBar();
}

if ($opts{t}){
    processTargets();
}else{
    processByLine();
}

#################################################
sub processTargets{
    foreach my $t (keys %targets){
        $current_target = $t;
        my @lines = VcfReader::searchByRegion
        (
            %sargs,
            chrom => $targets{$t}->{chrom},
            start => $targets{$t}->{start},
            end   => $targets{$t}->{end},
        );
        processLine($_) for @lines;
        $var_count++;#for progress bar
        updateProgressBar();
        outputGeneCounts();
    }
    outputGeneCounts();
    close $TMP; 
    my %contigs = VcfReader::getContigOrder($opts{i});#will already be indexed
    informUser("Sorting and outputting variants...\n");
    VcfReader::sortVcf
    (
        vcf => $tmpout,
        output => $OUT,
        contig_order => \%contigs,
    );
    informUser("Done\n");
}

#################################################
sub processByLine{
    my $VCF = VcfReader::openVcf($opts{i}); 
    #read line
    while (my $line = <$VCF>){
        processLine($line);
        $var_count++;#for progress bar
        updateProgressBar();  
    }

    if (defined $opts{f}){
        $functional_count += checkMatchingVariants();
    }elsif($opts{gene_counts}){
        outputGeneCounts();
    }
    close $VCF;

    updateProgressBar();  
    outputGeneList();
    close $OUT;

    informUser
    (
        "$functional_count variants matching criteria found from $var_count total".
        " variants.\n"
    );
}

#################################################
sub processLine{
    my $line = shift;
    return if $line =~ /^#/;#skip header
    chomp $line;
    my @split = split("\t", $line); 
    my ($chrom, $qual, $filter) = VcfReader::getMultipleVariantFields
    (
        \@split, 
        'CHROM', 
        'QUAL',
        'FILTER',
    );
    #skip whole line if not passing filter or variant quality settings
    if ($opts{pass}){
        return if $filter ne 'PASS';
    }
    
    if ($opts{v}){
        return if $qual < $opts{v};
    }
    #get alleles to arrays of 'functional' consequences
    my %allele_to_csq = getFunctionalCsq(\@split);
    if (defined $opts{f}){
    #if using matching variants go through each allele and each consequence and
    #put transcript IDs into a hash to a hash of sample names to variant IDs
        if (isNewChromosome($chrom)){
            #if new chrom check genes and clear collected data
            $functional_count += checkMatchingVariants();
        }
        recordMatchedSamples
        (
            \@split, 
            \%allele_to_csq,
        ); 
    }elsif( $opts{gene_counts}){
    #if using gene_counts 
        if (isNewChromosome($chrom)){
            outputGeneCounts();
        }
        recordGeneCounts
        (
            \@split, 
            \%allele_to_csq,
        ); 
    }else{
    #if not doing gene counts or matching variants simply print if we've got a
    #functional allele
        if (%allele_to_csq){
            my $line_ref = VcfReader::addVariantInfoField
            (
                line  => \@split,
                id    => "getFunctionalVariantsMatchedAlleles",
                value => join(",",  (sort {$a <=> $b} keys %allele_to_csq)),
            );
            print $OUT join("\t", @$line_ref) ."\n";
            $functional_count++;
        }
    }
}

#################################################
sub recordGeneCounts{
    my $split = shift;
    my $alleles = shift;
    my @matched_samples = (); 
    return if not %{$alleles};
    my $out_line_ref;
    if ($gene_count_opts{count_mode} eq 'allele_counts'){
        my @allele_counts = (); 
        foreach my $i (sort {$a <=> $b} keys %{$alleles}){
            my ($sc, $sn) = getSampleCounts
            (
                $split, 
                $i,
            );
            push @allele_counts, $sc;
            foreach my $annot (@{$alleles->{$i}}){
                my $tr = $annot->{$feature_id};
                my $symbol = $annot->{$symbol_id};
                my $gene   = $annot->{$gene_id};
                $gene   ||= $tr;
                $symbol ||= $gene;
                if ($current_target){
                    next if $gene ne $current_target;
                }
                $transcript_to_symbol{$gene} = $symbol;
                $gene_burden{$gene}->{SC} += $sc;
                $gene_burden{$gene}->{SN} ||= 0 ;
                $gene_burden{$gene}->{SN} = $gene_burden{$gene}->{SN} > $sn ?
                                            $gene_burden{$gene}->{SN} : $sn;
            }
            
        }
        $out_line_ref = VcfReader::addVariantInfoField
        (
            line  => $split,
            id    => "getFunctionalVariantsAlleleCounts",
            value => join(",", @allele_counts),
        );
    }else{
        my @biallelics = ();#array ref per allele of samples w/ findBiallelic tags
        if ($gene_count_opts{model} eq 'recessive'){
        #if wanting to use --gene_counts to find recessives
            my @hom = split
            (
                ",", 
                VcfReader::getVariantInfoField($split, "findBiallelicSamplesHom") 
            );
            my @het = split
            (
                ",", 
                VcfReader::getVariantInfoField($split, "findBiallelicSamplesHet") 
            );
            #entries are per allele, with multiple samples separated by '|' 
            for (my $i = 0; $i < @hom; $i++){
                my @s = grep { $_ ne '.' } split(/\|/, $hom[$i]); 
                push @s, grep { $_ ne '.' } split(/\|/, $het[$i]);
                push @biallelics, \@s;
            }
        }
        foreach my $i (sort {$a <=> $b} keys %{$alleles}){
            my @samp_with_allele = ();
            if ($gene_count_opts{model} eq 'recessive'){
                if (@{$biallelics[$i - 1]}){
                    push @samp_with_allele, @{$biallelics[$i - 1]};
                }
            }else{#getting gene counts for specified inheritance model
                my %samp_to_gt = VcfReader::getSampleCall
                (
                  line              => $split, 
                  sample_to_columns => \%sample_to_col,
                  minGQ             => $opts{g},
                  multiple          => \@samples,
                );        
                @samp_with_allele = addSamplesForGeneCounts
                (
                        \%samp_to_gt, 
                        $i,
                        $split,
                );
            }
            if (@samp_with_allele){ 
                push @matched_samples, join("|", @samp_with_allele); 
                #for INFO field of output line
                foreach my $annot (@{$alleles->{$i}}){
                    my $tr     = $annot->{$feature_id};
                    my $symbol = $annot->{$symbol_id};
                    my $gene   = $annot->{$gene_id};
                    $gene   ||= $tr;
                    $symbol ||= $gene;
                    if ($current_target){
                        next if $gene ne $current_target;
                    }
                    foreach my $sample (@samp_with_allele){
                        $gene_burden{$gene}->{$sample} = undef;
                        $transcript_to_symbol{$gene}   = $symbol;
                        # keys for each %{gene_burden{$gene}} will be undef for
                        # counting burden (# samples with qualifying variant)
                    }
                }
            }else{
                push @matched_samples, '.';
            }
        }
        $out_line_ref = VcfReader::addVariantInfoField
        (
            line  => $split,
            id    => "getFunctionalVariantsMatchedSamples",
            value => join(",", @matched_samples),
        );
    }
    $out_line_ref = VcfReader::addVariantInfoField
    (
        line  => $out_line_ref,
        id    => "getFunctionalVariantsMatchedAlleles",
        value => join(",",  (sort {$a <=> $b} keys %{$alleles})),
    );
    my $FHOUT = $OUT;
    if ($opts{t}){
        $FHOUT = $TMP;
    }
    print $FHOUT join("\t", @$out_line_ref ) . "\n";
    $functional_count++;
}

#################################################
sub recordMatchedSamples{
# this returns the samples with a matching allele/genotype for adding to INFO
# and also adds information to %matched_transcript hash to record samples with
# functional variants in a given transcript
    my $split = shift;
    my $alleles = shift;
    return if not %{$alleles};
    my @matched_samples = (); 
    my %samp_to_gt = VcfReader::getSampleCall
    (
          line              => $split, 
          sample_to_columns => \%sample_to_col,
          minGQ             => $opts{g},
          multiple          => \@samples,
    );
    my ($chrom, $pos, $ref, $alt) = VcfReader::getMultipleVariantFields
    (
        $split, 
        'CHROM', 
        'POS', 
        'REF',
        'ALT',
    );
    my $var_id = "$chrom:$pos-$ref,$alt";
    my $store_this_line = 0;
    #skip if not identical GTs in all affecteds and identical variants required
    if ($opts{e}){
        my @matched_gts, identicalGenotypes(\%samp_to_gt, $opts{n});
        return if not @matched_gts;
        my $i_match ;
        foreach my $i (sort {$a <=> $b} keys %{$alleles}){
            foreach my $matched_gt( @matched_gts ){
                #if require equal genotypes skip if this is not a matched allele
                my @m = split("/", $matched_gt); 
                if (first { $_ eq $i } @m ){#first is ok cos we never touch allele 0
                    $i_match++; 
                    last; 
                }
            }
            if ($i_match){
                if (my @samp_with_allele = addSamplesWithGt
                    (
                        \%samp_to_gt, 
                        \@matched_gts,
                        $split,
                    )
                ){
                    push @matched_samples, join("|", @samp_with_allele);
                    recordTranscriptSamplesAndVariants
                    (
                        \@{$alleles->{$i}},
                        \@samp_with_allele,
                        $var_id,
                    );
                    $store_this_line++;
                }else{
                    push @matched_samples, '.';
                }
            }else{
                push @matched_samples, '.';
            }
        }
    }else{
        foreach my $i (sort {$a <=> $b} keys %{$alleles}){
            if (my @samp_with_allele = addSamplesWithAllele
                (
                    \%samp_to_gt, 
                    $i,
                    $var_id,
                ) 
            ){
                push @matched_samples, join("|", @samp_with_allele);
                recordTranscriptSamplesAndVariants
                (
                    \@{$alleles->{$i}},
                    \@samp_with_allele,
                    $split,
                );
                $store_this_line++;
            }else{
                push @matched_samples, '.';
            }
        }
    }
    if ($store_this_line){
        my $line_ref = VcfReader::addVariantInfoField
        (
            line  => $split,
            id    => "getFunctionalVariantsMatchedSamples",
            value => join(",", @matched_samples),
        );
        $line_ref = VcfReader::addVariantInfoField
        (
            line  => $split,
            id    => "getFunctionalVariantsMatchedAlleles",
            value => join(",",  (sort {$a <=> $b} keys %{$alleles})),
        );
        $vcf_lines{$var_id} = join("\t", @$line_ref);
    }
}

#################################################
sub recordTranscriptSamplesAndVariants{
    my $annotations = shift;
    my $samples = shift;
    my $var_id = shift;
    foreach my $annot (@$annotations){
        my $tr = $annot->{$feature_id};
        $transcript_to_symbol{$tr} = $annot->{$symbol_id};
        foreach my $sample (@$samples){
            push @{ $transcript_sample_vars{$tr}->{ $sample } }, $var_id;
        }
    }
}

#################################################
sub isNewChromosome{
    my $chrom = shift;
    return 0 if $opts{t};
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
sub getFunctionalCsq{
    my $split = shift;
    my %allele_to_csq = (); #keys are allele numbers, values = array of
                            #'functional' csq
    #skip if no variant allele in affecteds 
    my %samp_to_gt = VcfReader::getSampleCall
    (
          line              => $split, 
          sample_to_columns => \%sample_to_col,
          minGQ             => $opts{g},
          multiple          => \@samples,
    );
    return if not haveVariant(\%samp_to_gt);
    
    #skip if not identical GTs in all affecteds and identical variants required
    my @matched_gts = (); 
    if ($opts{e}){
        push @matched_gts, identicalGenotypes(\%samp_to_gt, $opts{n});
        return if not @matched_gts;
    }
    
    #get Allele frequency annotations if they exist and we're filtering on AF
    my %af_info_values = (); 
    if (%af_info_fields){ 
        %af_info_values = getAfInfoValues($split);
    }
    
    # get allele frequencies in  @samples if using -u/--max_sample_allele_frequency
    my %s_allele_counts = ();
    my $allele_count;
    if ($opts{u}){
        %s_allele_counts = VcfReader::countAlleles
        (
            minGQ => $opts{g},
            line  => $split,
        );
        map { $allele_count += $s_allele_counts{$_} } keys %s_allele_counts;
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
        #if CADD filtering skip if less than user-specified cadd filter
        #if no CADD score for allele then we won't skip unless using --
        if (@cadd_scores){
            my $score = $cadd_scores[$i - 1];
            if ($score ne '.'){
                next ALT if $score < $opts{c};
            }elsif ($opts{skip_unpredicted}){
                next ALT;
            }
        }
        
        #filter if AF annotations > MAF
        if (%af_info_fields){ #check for annotateSnps.pl etc. frequencies
           next ALT if ( alleleAboveMaf($i - 1, \%af_info_values) );
        }

        # filter if frequency in @samples is greater or equal to 
        # -u/--max_sample_allele_frequency
        if ($opts{u} and $allele_count > 0){
            #my $freq = $s_allele_counts{$i}/(@samples*2);
            my $freq = $s_allele_counts{$i}/$allele_count;
            next ALT if $freq >= $opts{u}; 
        }

        #get all consequences for current allele
        my @a_csq = ();
        if ($opts{m} eq 'vep'){
            @a_csq = grep { $_->{allele} eq $alts_to_vep_allele[$i-1] } @csq;
        }else{
            @a_csq = grep { $_->{allele} eq $alleles[$i] } @csq;
        }
        
        #0 or 1 for each allele if pathogenic - i.e. we keep all annotations
        # if allele is flagged as pathogenic in ClinVar
        # ClinVar annotations come via annotateSnps.pl using a ClinVar file
        my @clinvar_path = keepClinvar($split, scalar(@alleles)-1);
        #skip if VEP/SNPEFF annotation is not the correct class 
        # (or not predicted damaging if using that option)
        # and allele is not flagged as pathogenic in ClinVar
CSQ:    foreach my $annot (@a_csq){
            if ($clinvar_path[$i-1] or consequenceMatchesClass($annot, $split)){
                push @{$allele_to_csq{$i}}, $annot;
            }
        }
    }#each ALT
    return %allele_to_csq;
}


#################################################
sub outputGeneCounts{
    foreach my $g (keys %gene_burden){
        my $count = 0;
        my $without = 0;
        if ($gene_count_opts{count_mode} eq 'allele_counts'){
            $count = $gene_burden{$g}->{SC};
            $without = $gene_burden{$g}->{SN};
        }else{
            $count = keys %{$gene_burden{$g}};
            $without = @samples - $count;
        }
        
        print $GENECOUNTS join
        (
            "\t",
            $g,
            $transcript_to_symbol{$g},
            $count,
            $without,
        ) . "\n";
    }
    %gene_burden = ();
}

#################################################
sub getSampleCounts{
    #return no. samples with allele and total number samples
    my $line = shift;
    my $allele = shift;
    my $field = 'AC';
    my $count = 0;
    my $number = 0;
    my $chrom  =  VcfReader::getVariantField($line, 'CHROM');
    if (exists $info_fields{AC_Het} and exists $info_fields{AC_Hom}){
        my $hemi = 0;
        if ($chrom =~ /^(chr)X|Y/ and exists $info_fields{AC_hemi}){
            my $hemis = VcfReader::getVariantInfoField($line, "AC_Hemi");
            $hemi = (split ",", $hemis)[$allele - 1];
        }
        my $hets = VcfReader::getVariantInfoField($line, "AC_Het");
        my $homs = VcfReader::getVariantInfoField($line, "AC_Hom");
        my $het  = (split ",", $hets)[$allele - 1];
        my $hom  = (split ",", $homs)[$allele - 1];
        if ($gene_count_opts{model} eq 'dominant' and $hom > 0){
            $count = 0;#discard if there's a homozygote for a dominant variant(?)
        }else{
        #no reliable way to find comp. hets from ExAC allele counts
            $count = $het + $hemi + $hom;
        }
    }elsif (exists $info_fields{AC}){
        my $ac = VcfReader::getVariantInfoField($line, 'AC');
        $count =  ((split ",", $ac)[$allele - 1])/2;#assume diploidy
    }else{
        die "AC INFO field is required for using --gene_counts in ".
            "'allele_counts' mode\n";
    }
    if ($gene_count_opts{sample_count}){
        if ($chrom =~ /^(chr)Y/ and exists $gene_count_opts{males}){
            $number = $gene_count_opts{males};
        }else{
            $number = $gene_count_opts{sample_count};
        }
    }else{
        my $an_field; 
        if (exists $info_fields{AN_Adj}){
            $an_field = 'AN_Adj';
        }elsif (exists $info_fields{AN}){
            $an_field = 'AN';
        }else{
            die "AN INFO field is required for using --gene_counts in ".
                "'allele_counts' mode\n";
        }
        my $an = VcfReader::getVariantInfoField($line, $an_field);  
        $number = int($an/2) + ($an % 2);#assume diploidy
    }
    return ($count, $number);
}

#################################################
sub addSamplesForGeneCounts{
#will already have added samples if using 
#recessive gene count model
    my ($gts, $allele, $line) = @_; 
    my @samp = ();
    foreach my $s ( keys %{$gts} ){
        if ($gene_count_opts{model} eq 'dominant'){
            if ($gts->{$s} =~ /^$allele[\/\|]$allele$/){
            #found a homozygote - discard all(?)
                if (checkGtPl(["$allele/$allele"], $s, $line )){
                    return;
                }
            }
        }
        if ($gts->{$s} =~ /^($allele[\/\|]\d+|\d+[\/\|]$allele)$/){
            if (checkAllelePl($allele, $s, $line)){
                push @samp, $s;
            }
        }
    }
    return @samp;
}

#################################################
sub addSamplesWithGt{
    my ($gts, $matched, $line) = @_; 
    my @samp = ();
    foreach my $s ( keys %{$gts} ){
        next if ($gts->{$s} =~ /^\.[\/\|]\.$/ );
        (my $gt = $gts->{$s}) =~ s/[\|]/\//;
          #don't let allele divider affect whether a genotype matches
        if ($opts{pl}){
            if (checkGtPl($matched, $s, $line )){
                push @samp, $s;
            }
        }else{
            if ( first { $_ eq $gt } @$matched ){
                push @samp, $s ;
            }
        }
    }
    return @samp;
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
sub checkClinVarInfo{
    if ($opts{clinvar} and lc($opts{clinvar}) eq 'disable'){
        return 0;
    }
    my $iclvar = 0;
    foreach my $f ( qw / ClinVarPathogenic AS_CLNSIG / ) { 
        if (exists $info_fields{$f} and $info_fields{$f}->{Number} eq 'A'){
            informUser
            (
                "Identified '$f' INFO field - alleles marked as pathogenic ".
                "will be considered as having a 'functional' consequence ".
                "regardless of functional consequence. To change this ".
                "behaviour run with the option '--clinvar disable'\n"
            );
            $iclvar++;
        }elsif(exists $info_fields{$f}){
            informUser
            (
                "WARNING: Ignoring INFO field '$f' because 'Number' is given ".
                "as $info_fields{$f}->{Number} rather than the expected 'A'.\n"
            );
            delete $info_fields{$f};#remove so we don't try using this annotation later
        }
    }
    if ($iclvar){
        if (not $opts{clinvar}){
            return 'all';
        }elsif (lc($opts{clinvar}) eq 'all'){
            return 'all';
        }elsif (lc($opts{clinvar}) eq 'no_conflicted'){
            if (not exists $info_fields{ClinVarConflicted} and 
                not exists $info_fields{AS_CLNSIG}
                #can glean from CLNSIG field if more than one type of assertion made
            ){
                informUser
                ( 
                    "WARNING: no ClinVarConflicted or CLNSIG INFO field found in header.".
                    " All variants with ClinVarPathogenic annotations will be".
                    " kept.\n"
                );
                return 'all';
            }
            return 'no_conflicted';
        }else{
            die "Unrecognised value '$opts{clinvar}' passed to --clinvar option.\n";
        }
    }else{
        if ($opts{clinvar}){
            informUser
            (
                "WARNING: neither ClinVarPathogenic nor AS_CLNSIG INFO fields".
                " were found in header. ClinVarPathogenic variants will not be".
                " identified.\n"
            );
        }
        return 0;
    }
}

#################################################
sub keepClinvar{
#returns an array with a value for each allele
#values are 1 if pathogenic and and 0 if not
    my $v = shift;
    my $n_alts = shift;
    if (not $keep_clinvar){
        return map { 0 } 0..$n_alts 
    }
    my @c_path = ();
    my @d_path = ();
    my @c_conf = ();
    my @d_conf = ();
    if (exists $info_fields{ClinVarPathogenic}){
        my $path = VcfReader::getVariantInfoField
        (
            $v,
            'ClinVarPathogenic',
        );
        if (defined $path){ 
            @c_path = split(",", $path);
            if ($keep_clinvar eq 'no_conflicted'){
                @c_conf = split(",", VcfReader::getVariantInfoField
                    (
                        $v,
                        'ClinVarConflicted',
                    )
                );
                @c_path = map {$c_path[$_] == 1 and $c_conf[$_] == 0} 0..$#c_path;
            }
        }
    }
    if (exists $info_fields{AS_CLNSIG}){
        my $path =  VcfReader::getVariantInfoField
        (
            $v,
            'AS_CLNSIG',
        );
        if (defined $path){ 
            my @csig = split(",", $path);
            foreach my $c (@csig){
                my @path = split(/\|/, $c);
                if (grep {$_ eq '4' or $_ eq '5'} @path){
                    push @d_path, 1;
                    if (grep {$_ eq '2' or $_ eq '3'} @path){
                        push @d_conf, 1;
                    }else{
                        push @d_conf, 0;
                    }
                }else{
                    push @d_path, 0;
                    push @d_conf, 0;
                }
            }
            if ($keep_clinvar eq 'no_conflicted'){
                @d_path = map {$d_path[$_] == 1 and  $d_conf[$_] == 0} 0..$#d_path;
            }
        }
    }
    if (@c_path){
        if (@d_path){
            if ($keep_clinvar eq 'no_conflicted'){
                #an annotation of 0 from the ClinVarPathogenic might mean 
                #no info rather than designated as benign, so only consider as
                #conflicted if ClinVarConflicted annotation is 1 or $d_path[$_] == 0
                #in this case, if we have annotations from VCF in @d_path,
                #about what is in @c_path
                return map { $d_path[$_] == 1  and $c_conf[$_] == 0} 0..$#d_path;
            }else{  
                #don't care about conflicted annotations
                #keep if either source is flagged as pathogenic
                return map {$d_path[$_] == 1 or $c_path[$_] == 1} 0..$#d_path;
            }
        }else{
            return @c_path;
        }
    }elsif(@d_path){
        return @d_path;
    }else{
        return map { 0 } 0..$n_alts 
    }
}
 
#################################################
sub outputGeneList{
    return if not $LIST;
    foreach my $k (sort keys %genes_to_transcripts){
        print $LIST "$k:\t" . join
        (
            ",", 
            @{ $genes_to_transcripts{$k} }
        ) . "\n";
    }
    close $LIST;
}
#################################################
sub checkMatchingVariants{
    my $count = 0; 
    my $min_matching = $opts{n} ? $opts{n} : scalar @samples;
    my %lines_to_print = (); 
    foreach my $tr (keys %transcript_sample_vars){
        if (scalar (keys %{ $transcript_sample_vars{$tr} }) > $min_matching ){
            my %new_lines = map { $vcf_lines{$_} => undef } 
                            map { @{$transcript_sample_vars{$tr}->{$_}} } 
                            keys %{ $transcript_sample_vars{$tr} } ;
            %lines_to_print =  (%new_lines, %lines_to_print); 
            my $symbol = $transcript_to_symbol{$tr};
            push @{$genes_to_transcripts{$symbol}}, $tr;
        }
    }
    #convert lines to split array refs for sorting
    my @lines_out =  map { [ split("\t", $_ )] } keys %lines_to_print;
    @lines_out = VcfReader::sortByPos( \@lines_out ); 
    foreach my $l (@lines_out){
        print $OUT join("\t", @$l) . "\n";
        $count++;
    }
    %transcript_sample_vars = ();
    %vcf_lines = ();
    %transcript_to_symbol = ();
    return $count;
}

#################################################
sub readTargetGenesFile{
    return if not $opts{t}; 
    open (my $TARGETS, "<", $opts{t}) or die "Could not open -t/--target_genes"
                                             . " file '$opts{t}': $!\n";
    my @cols = 
    (
        "Ensembl Gene ID",
        "Chromosome Name",
        "Gene Start (bp)",
        "Gene End (bp)",
    );
    my $header = <$TARGETS>; 
    my @split = split("\t", $header);
    my %h = ();
    {
        no warnings 'uninitialized';
        foreach my $c (@cols){
            my $i = 0;
            $i++ until uc($split[$i]) eq uc($c) or $i > $#split;
            if ($i > $#split){
                die <<EOT
Required column '$c' not found in -t/--target_genes file '$opts{t}'. Appropriate
files can be obtained by outputting your gene lists from Ensembl's biomart with
the attributes "Ensembl Gene ID", "Chromosome Name", "Gene Start (bp)", "Gene
End (bp).
EOT
                ;  
            }else{
                $h{$c} = $i;
            }
        }
         
    }   
    while (my $l = <$TARGETS>){
        my @s = split("\t", $l);
        my $id = $s[$h{"Ensembl Gene ID"}];
        $targets{$id}->{chrom} =  $s[$h{"Chromosome Name"}];
        $targets{$id}->{start} =  $s[$h{"Gene Start (bp)"}];
        $targets{$id}->{end}   =  $s[$h{"Gene End (bp)"}];
    }
}

#################################################
sub checkGeneCountArgs{
    informUser
    (
        "WARNING: No 'samples specified - all samples will be used for ".
        "--gene_counts option\n"
    ) if not @samples;
    my @keys = qw / file model count_mode sample_count /;
    my @split = split(',', $opts{gene_counts}); 
    @gene_count_opts{@keys} = @split;
    $gene_count_opts{model} ||= 'dominant';
    $gene_count_opts{count_mode} ||= 'genotypes'; 
    if ( $gene_count_opts{model} ne 'dominant' and
         $gene_count_opts{model} ne 'recessive' and 
         $gene_count_opts{model} ne 'both' 
    ){
        die "ERROR: Unrecognised model '$gene_count_opts{model}' for --gene_counts argument\n";
    }

    if ($gene_count_opts{model} eq 'recessive'){
        foreach my $f (qw /findBiallelicSamplesHom findBiallelicSamplesHet/){
            if (not exists $info_fields{$f}){
                die "ERROR: Could not find '$f' INFO field entry in VCF header".
                    " In order to use recessive mode with the --gene_counts ".
                    "please run findBiallelic.pl on your input first.\n";
            }
        }
    }
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
        );
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
        );
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
    if (not @classes){
        @classes = grep { $all_classes{$_} eq 'default' } keys %all_classes;
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
        @biotypes = grep { $all_biotypes{$_} eq 'filter' } keys %all_biotypes;
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
sub getAndCheckScoreFilters{
    return if not @score_filters;
    my @filters = (); 
    my %csq_add = ();
FLT: foreach my $s (@score_filters){
        my %f =  VcfhacksUtils::getScoreFilter($s);
        foreach my $fld (@{$f{field}}){
            (my $fb = $fld) =~ s/^\(//;#we may have used ( to specify precedence
            if (not exists $csq_header{lc($fb)}){
                informUser
                (
                    "WARNING: No '$fb' field found in CSQ header of ".
                    "VCF. Cannot use --score_filter expression '$s' ".
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
sub openOutput{
    my $OUT_FH;
    my $LIST_FH;
    my $COUNT_FH;
    if ($opts{o}) {
        open( $OUT_FH, ">", $opts{o} ) || die "Can't open $opts{o} for writing: $!\n";
    }
    else {
        $OUT_FH = \*STDOUT;
    }
    if (defined $opts{f}) {
        if ( $opts{f} eq '' ){
          #user specified --list option but provided no argument
            $LIST_FH = \*STDERR;
        }else{
            open( $LIST_FH, ">", $opts{f} )
              or die "Can't open $opts{f} for writing: $!\n";
        }
    }
    if ($opts{gene_counts}){
        open( $COUNT_FH, ">", $gene_count_opts{file} )
          or die "Can't open $gene_count_opts{file} for writing: $!\n";
    }
    return ($OUT_FH, $LIST_FH, $COUNT_FH);        
}

#################################################
sub writeOptionsToHeader{
    #print meta header lines
    my $FH = $OUT;
    if ($TMP){
        $FH = $TMP;
    }
    print $FH join("\n", grep { /^##/ } @header) . "\n" ;
    #add header line detailing program options
    print $FH VcfhacksUtils::getOptsVcfHeader(%opts) . "\n"; 
    print $FH "$header[-1]\n";
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
    }
}

#################################################
sub getProgressBar{
    my $msg = <<EOT
WARNING: -b/--progress option requires installation of the perl module 'Term::ProgressBar'.
WARNING: 'Term::ProgressBar' module was not found. No progress will be shown.
EOT
;
    eval "use Term::ProgressBar; 1 " or informUser($msg) and return;
    
    if ($opts{t}){
        my $ts = keys %targets;
        informUser("Processing $ts gene targets\n");
        return Term::ProgressBar->new
        (
            {
                name  => "Filtering",
                count => $ts, 
                ETA   => "linear",
            }
        );
    }
    if ( $opts{i} eq "-" ) {
        informUser("Can't use -b/--progress option when input is from STDIN\n");
        return;
    }else{
        informUser("Counting variants in input for progress monitoring.\n"); 
        $total_vars = VcfReader::countVariants($opts{i});
        informUser("$opts{i} has $total_vars variants.\n");
        return Term::ProgressBar->new
        (
            {
                name  => "Filtering",
                count => $total_vars,
                ETA   => "linear",
            }
        );
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
sub identicalGenotypes{
    my $gts = shift;
    my $min_matching = shift;
    $min_matching = defined $min_matching ? $min_matching : scalar keys %$gts;
    my @matched = (); 
    my %matches = (); #keys are genotypes, values are arrays of samples with GT
    foreach my $s (keys %$gts){
        if ( $gts->{$s} =~ /^\.[\/\|]\.$/ ){    
            next;
        }elsif ( $gts->{$s} =~ /^0[\/\|]0$/ ){
            next;
        }
        (my $current_geno = $gts->{$s}) =~ s/[\|]/\//;
          #don't let allele divider affect whether a genotype matches or not
        push @{ $matches{$current_geno} } , $s;
        if ($min_matching == scalar keys %$gts){
        #if we require all to match, bail out if there's more than one genotype
            return 0 if keys %matches > 1;
        }
    }
    foreach my $k (keys %matches){
        my $m = 0;
        foreach my $s (@{$matches{$k}}){
            $m++;
        }
        push @matched, $k if $m >= $min_matching;
    }
    return @matched;
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
    if ($opts{m} eq 'vep'){
        return consequenceMatchesVepClass($annot);
    }else{
        return consequenceMatchesSnpEffClass($annot, $l);
    }   
}

#################################################
sub consequenceMatchesVepClass{
    my $annot = shift;
    #intergenic variants have no feature associated with them - skip
    return 0 if $annot->{consequence} eq "intergenic_variant";
    #skip unwanted biotypes
    return 0 if (exists $biotype_filters{$annot->{biotype}}) ;
    #skip non-canonical transcripts if --canonical_only selected
    if ($opts{canonical_only}) {
        return 0 if ( not $annot->{canonical} );
    }
    
    my @anno_csq = split( /\&/, $annot->{consequence} );
    #skip NMD transcripts
    return 0 if ( grep { /NMD_transcript_variant/i } @anno_csq );

    #score filters trump annotation class
    foreach my $scf (@score_exp){
        return 1 if VcfhacksUtils::scoreFilter($scf, $annot);
    }

ANNO: foreach my $ac (@anno_csq){
        $ac = lc($ac);#we've already converted %class_filters to all lowercase
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
    #skip variants with undef features (intergenic variants)
    return if not defined $annot->{feature_id};
    #skip unwanted biotypes
    return 0 if exists $biotype_filters{lc $annot->{transcript_biotype} };

    #score filters trump annotation class
    foreach my $scf (@score_exp){
        return 1 if VcfhacksUtils::scoreFilter($scf, $annot);
    }
    
    my @anno_csq = split( /\&/, $annot->{annotation} );
ANNO: foreach my $ac (@anno_csq){
        $ac = lc($ac);#we've already converted %class_filters to all lowercase
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

getFunctionalVariants.pl - retrieve variants according to their functional annotation

=head1 SYNOPSIS

    getFunctionalVariants.pl -i <variants.vcf> [options]
    getFunctionalVariants.pl --help (show help message)
    getFunctionalVariants.pl --manual (show manual page)

=cut 

=head1 ARGUMENTS

=over 8 

=item B<-i    --input>

VCF file annotated with Ensembl's variant_effect_predictor.pl (VEP) script or SnpEff.

=item B<-o    --output>

File to print output (optional). Will print to STDOUT by default.

=item B<-m    --mode>

This program will attempt to detect the format of your input automatically by looking for VEP or SnpEff annotations, but you may specify either 'vep' or 'snpeff' with this option to select the mode employed by the script if you have a file with both VEP and SnpEff annotations. By default, if both annotations are present and the program is run without this option, VEP annotations will be used. NOTE: Only SnpEff annotations using the more recent 'ANN' style annotations rather than the older 'EFF' style are recognised by this program.

=item B<-s    --samples>

Only keep variants present in one of these samples.  You may specify 'all' (without quotes) instead of a sample name to select all samples in the VCF.  

=item B<-f    --find_shared_genes>

Use this switch to only output variants that make up 'functional' variants in the same genes for the samples specified by the -s/--samples argument. If -s/--samples option is not specified this will act as if ALL samples in the VCF were specified.  This will also return a list of genes containing 'functional' variants in these samples. If a filename is passed to this option the list will be printed to file, otherwise the list will be printed to STDERR.

=item B<-n    --num_matching>

If s/--samples or -f/--find_shared_genes arguments are specified use this option to specify the minimum number of samples with variants in the same gene before outputting variants (and genes). 

=item B<-e    --equal_genotype>

If -s/--samples argument is specified use this flag if you only want to keep variants with identical genotypes in each sample (or a minimum number of samples as specified by the -n/--num_matching option).

=item B<--gene_counts>

Give a filename for outputting the counts of gene IDs vs number of samples with qualifying variants (e.g. for input to a burden test). Optionally the user may also add the 'mode' to use for the gene counts, separated by a comma after the filename (e.g. --gene_counts gene_counts.txt,recessive). Valid modes are 'dominant' (only samples with heterozygous variants counted), 'recessive' (requires annotations from findBiallelic.pl to identify samples with compound het or homozygous variants) or 'both' (count a sample regardless of whether it is heterozygous or homozygous). If using the 'recessive' mode you will need to annotate your variants with findBiallelic.pl first and use consistent settings for functional/in silico filters between both programs.

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

Specify in silico prediction scores to filter on. If running on VEP annotations PolyPhen, SIFT and Condel scores given by VEP will be used. If running on SnpEff annotations scores given by running SnpSift's 'dbnsfp' mode will be used.  

B<VEP mode:> Specify SIFT, PolyPhen or Condel labels or scores to filter on. Add the names of the programs you want to use, separated by spaces, after the --damaging option. By default SIFT will keep variants labelled as 'deleterious', Polyphen will keep variants labelled as 'possibly_damaging' or 'probably_damaging' and  Condel will keep variants labelled as 'deleterious'.

If you want to filter on custom values specify values after each program name in the like so: 'polyphen=probably_damaging'. Seperate multiple values with commas - e.g. 'polyphen=probably_damaging,possibly_damaging,unknown'. You may specify scores between 0 and 1 to filter on scores rather than labels - e.g. 'sift=0.3'. For polyphen, variants with scores lower than this score are considered benign and filtered, for SIFT and Condel higher scores are considered benign.

Valid labels for SIFT: deleterious, tolerated

Valid labels for Polyphen: probably_damaging, possibly_damaging, benign, unknown

Valid labels for Condel : deleterious, neutral

To use default values for all three programs use 'all' (i.e. '--damaging all').

The default behaviour is to only keep variants predicted as damaging by ALL programs specified, although if the value is not available for one or more programs than that program will be ignored for filtering purposes.

B<SnpEff mode:> Specify one of the following annotations provided by dbNSFP (your input must have been annotated using SnpSift's dbnsfp mode): 

    dbNSFP_LRT_pred
    dbNSFP_MutationAssessor_pred
    dbNSFP_MutationTaster_pred
    dbNSFP_Polyphen2_HVAR_pred
    dbNSFP_SIFT_pred

You may omit the 'dbNSFP_' from the beginning if you wish. As with the VEP scores, you may specify 'all' to use the default values for all programs. Choosing custom values is also performed in the same way as for VEP annotations (e.g. dbNSFP_MutationAssessor_pred=H,M,L). 

The default scores considered 'damaging' are the following:

    dbNSFP_LRT_pred=D
    dbNSFP_MutationAssessor_pred=H,M
    dbNSFP_MutationTaster_pred=A,D
    dbNSFP_Polyphen2_HVAR_pred=P,D
    dbNSFP_SIFT_pred=D

Other valid scores, not used by default:

    dbNSFP_LRT_pred=N
    dbNSFP_MutationAssessor_pred=L,N
    dbNSFP_MutationTaster_pred=N,P
    dbNSFP_Polyphen2_HVAR_pred=B
    dbNSFP_SIFT_pred=T

The default behaviour is to only keep variants predicted as damaging by ALL programs specified, although if the value is not available for one or more programs than that program will be ignored for filtering purposes.

B<WARNING:> At the time of writing, SnpSift (latest version - SnpSift version 4.1l, build 2015-10-03) does not properly annotate the prediction values per allele. Multiple scores for the same allele are separated by commas, as are values per allele, therefore it is not possible to determine if a prediction is for a particular allele. For this reason, getFunctionalVariants.pl will use the highest present prediction value found for a variant in order to decide whether or not to filter an allele when using SnpEff/SnpSift annotations.

=item B<-k    --keep_any_damaging>

If using multiple programs for filters for --damaging argument use this flag to keep variants predicted to be damaging according to ANY of these programs.

=item B<--skip_unpredicted>

Skip alleles that do not have a score from one or more programs specified by the -d/--damaging argument or if using the -c/--cadd_filter option and an allele has no CADD phred score. The --keep_any_damaging argument will override this behaviour for -d/--damaging predictions if any of the available predictions are considered damaging.

=item B<--biotype_filters>

By default features/transcripts with the following biotypes are ignored:

    IG_C_pseudogene
    IG_J_pseudogene
    IG_V_pseudogene
    nonsense_mediated_decay
    non_stop_decay
    processed_pseudogene
    pseudogene
    retained_intron
    transcribed_processed_pseudogene
    transcribed_unprocessed_pseudogene
    TR_J_pseudogene
    TR_V_pseudogene
    unitary_pseudogene
    unprocessed_pseudogene

You may override this filter by specifying biotypes to filter with this option or prevent biotype filtering using the --no_biotype_filtering option (see below). 

Filtering will not occur if any non-filtered biotype contains a 'functional' variant. For example, with the default settings, if a splice_donor_variant only affects a 'pseudogene' transcript then this variant will be filtered, but if it also affects a 'protein_coding' transcript then it will not be filtered. 

The 'data/biotypes.tsv' file contains a list of valid biotypes and the default behaviour of this program (i.e. 'keep' or 'filter'). 

=item B<--no_biotype_filtering>

Use this flag to consider consequences affecting ALL biotypes.

=item B<-a    --allele_frequency>

Use a value between 0.00 and 1.00 to specify allele frequencey filtering for annotations from annotateSnps.pl, filterOnEvsMaf.pl or filterVcfOnVcf.pl if you've previously run these programs to annotate your VCF. If an allele frequency is available for an allele it will be filtered if equal to or greater than the value specfied here. 

Note: allele frequencies added by VEP are not used for filtering as they check the allele frequency at the site, not of the specific alleles in your variant call.

=item B<-j    --custom_af>

If using the --af/--allele_frequency option and your data contains allele frequency fields from sources not recognised by this program, you may give the name of these allele frequency INFO fields here and they will be used for filtering in addition to the default fields. Note that these annotations must contain an annotation per ALT allele (i.e. the INFO field header must specify 'Number=A') to work properly and the allele frequency should be expressed as a number between 0.00 and 1.00 in order to be compatible with the default allele frequency fields recognised by this program.

=item B<-u  --max_sample_allele_frequency>

Use this option to specify an allele frequency (between 0.00 and 1.00) for filtering alleles in your VCF. Alleles present at this frequency or higher in your samples of interest will be filtered. If -s/--samples argument is specified, only these samples will be used for calculating the allele frequency, otherwise all samples in your VCF will be used.

=item B<-v    --var_qual>

Minimum variant Phred-like quality score to consider. Variants with a QUAL field lower than this value will be filtered. Default is 20.

=item B<-g    --gq>

Minimum genotype qualities to consider. Only applicable when using the -s/--samples option. Any genotype call below this threshold will be considered a no call. Default is 20

=item B<--pl>

Minimum 0-based phred-scale genotype likelihood (see the VCF spec for details) for alternative genotypes. Only applicable when using the -s/--samples option. When considering a given allele, if the sample has a PL below this value for a genotype not including this allele, the sample will not be considered a carrier of that allele. Default - not used.

=item B<--pass_filters>

Only consider variants with a PASS filter field. If the FILTER field for variant is not PASS the variant will be skipped.

=item B<-b    --progress>

Show a progress bar while working.

=item B<-h    --help>

Show the program's help message.

=item B<--manual>

Show the program's manual page.

=back

=cut

=head1 EXAMPLES

 getFunctionalVariants.pl -i input.vcf -o out.vcf 
 #output variants with a default 'functional' consequence (annotated by VEP or SnpEff)

 getFunctionalVariants.pl -i input.vcf -o out.vcf -s sample1 sample2 -f shared_genes.txt
 #as above but only for variants where sample1 or sample2 contain a variant allele with a 'functional' consequence

 getFunctionalVariants.pl -s all -i input.vcf -o out.vcf -n 2 -af 0.001  -f shared_genes.txt 
 #


=cut

=head1 DESCRIPTION

In its simplest form this program will print specific variant classes from a VCF file annotated with either Ensembl's variant_effect_predictor.pl program or SnpEff and filter out other variant classes. Input must be a VCF annotated by the variant_effect_predictor.pl program using the '--vcf' option or a VCF annotated with SnpEff using the (now default) ANN style annotations.

As well as outputting variants on the basis of functional annotation, this program can identify genes with functional variants in specific samples using the --samples and --find_shared genes options.

=cut

=head1 AUTHOR

David A. Parry

University of Edinburgh

=head1 COPYRIGHT AND LICENSE

Copyright 2015, David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut


