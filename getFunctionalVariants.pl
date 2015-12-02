#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use POSIX qw/strftime/;
use Data::Dumper;
use List::Util qw ( first ) ;
use FindBin;
use lib "$FindBin::Bin/lib";
use ParsePedfile;
use VcfReader;
use VcfhacksUtils;

my $progressbar;
my @samples = ();
my @classes = (); 
my @add_classes = ();
my @damaging = ();
my @biotypes = ();
my @gene_lists = ();
my @custom_af = ();

my %opts = 
(
    classes             => \@classes,
    add_classes         => \@add_classes,
    d                   => \@damaging,
    list                => \@gene_lists,
    s                   => \@samples,
    biotype_filters     => \@biotypes,
    j                   => \@custom_af,
);

GetOptions(
    \%opts,
    "add_classes=s{,}" ,
    "a|af=f" ,
    "biotype_filters=s{,}",
    "b|progress" ,
    "canonical_only" ,
    "classes=s{,}" ,
    "consensus_splice_site",
    "c|cadd_filter=f",
    "d|damaging=s{,}" ,
    "e|equal_genotypes",
    "f|find_shared_genes:s", 
    "g|gq=i",
    "h|?|help" ,
    "i|input=s" ,
    "j|custom_af=s{,}",
    "k|keep_any_damaging" ,
    "maf=f" ,
    "manual" ,
    "m|mode" ,
    "no_biotype_filtering",
    "n|num_matching=i",
    "o|output=s" ,
    "pass" ,
    "pl=f",
    "s|samples=s{,}",
    "skip_unpredicted" ,
    "v|var_quality=i",
    "u|max_sample_allele_frequency=f",
    #"l|list=s{,}",
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
        -message =>  "-a/--af option requires a value between 0.00 and 1.00 to ".
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
        pod2usage
        (
            -message =>  "-u/--max_sample_allele_frequency option requires ". 
                         "samples to be specified with the -s/--samples ".
                         "argument.\n", 
            -exitval => 2
        );
    }
}

if (defined $opts{f}){
    pod2usage   
    (
        -message => "--find_shared_genes option requires at least one " . 
                    "sample to be specified with the -s/--samples argument.\n", 
        -exitval => 2
    ) if not @samples;
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
}

#check VEP/SNPEFF header and get annotation order 
# - will set $opts{m} if not already defined
my %csq_header = getAndCheckCsqHeader();
  # hash of functional annotations from SnpEff/VEP to their annotation index
    
my $feature_id = $opts{m} eq 'vep' ? "feature" : "feature_id";
my $symbol_id = $opts{m} eq 'vep' ? "symbol" : "gene_name";

#set default consequence fields to retrieve 
my @csq_fields = getCsqFields();#consequence fields to retrieve from VCF

#and check variant consequence classes are acceptable
my %class_filters = map { $_ => undef } getAndCheckClasses();

#check in silico prediction classes/scores are acceptable
my %in_silico_filters = getAndCheckInSilicoPred();
  #hash of prediction program names and values to filter

#and check biotype classes are acceptable
my %biotype_filters = map { $_ => undef } getAndCheckBiotypes();


#if 'all' is listed in --samples option add all samples
if (grep { /^all$/i } @samples){
    @samples = ();
    push @samples, keys %sample_to_col;
}
@samples = VcfhacksUtils::removeDups(@samples);

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

######PRE-PROCESSING######

#Get filehandle for output (STDOUT or file) and optionally gene list (STDERR or file)
my ($OUT, $LIST) = openOutput();

#write header and program run parameters
writeOptionsToHeader();

#start progress bar 
my $next_update = 0;
my $total_vars  = 0; 
my $var_count   = 0;
if ($opts{b}) {
    if ( $opts{i} eq "-" ) {
        informUser("Can't use -b/--progress option when input is from STDIN\n");
    }else{
        my $msg = <<EOT
WARNING: -b/--progress option requires installation of the perl module 'Term::ProgressBar'.
WARNING: 'Term::ProgressBar' module was not found. No progress will be shown.
EOT
;
        eval "use Term::ProgressBar; 1 " or informUser($msg);
        if (not $@){
            informUser("Counting variants in input for progress monitoring.\n"); 
            $total_vars = VcfReader::countVariants($opts{i});
            informUser("$opts{i} has $total_vars variants.\n");
            $progressbar = Term::ProgressBar->new(
                {
                    name  => "Filtering",
                    count => $total_vars,
                    ETA   => "linear",
                }
            );
        }
    }
}

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


my $VCF = VcfReader::openVcf($opts{i}); 
#read line
LINE: while (my $line = <$VCF>){
    next if $line =~ /^#/;#skip header
    chomp $line;
    updateProgressBar();  
    $var_count++;#for progress bar
    my @split = split("\t", $line); 
    my ($chrom, $pos, $qual, $filter) = VcfReader::getMultipleVariantFields
    (
        \@split, 
        'CHROM', 
        'POS', 
        'QUAL',
        'FILTER',
    );
    if (defined $opts{f}){
        #check if new chromosome - 
        if (exists $contigs{$chrom}){
            #require chroms to be ordered together
            if ($contigs{$chrom} != scalar(keys %contigs) - 1){
                die
    "Encountered a variant for contig $chrom with at least one variant for a ".
    "different contig inbetween. Your VCF input must be sorted such that all contigs".
    " are kept together when using the -f/--find_shared_genes option. Please sort ".
    "your input and try again.\n";
            }
        }else{
            #if new chrom check biallelics and clear collected data
            $contigs{$chrom} = scalar(keys%contigs);
            checkMatchingVariants() if $contigs{$chrom} > 0;
        }
    }

    #skip if FILTER != PASS and PASS required
    if ($opts{pass}){
        next LINE if $filter ne 'PASS';
    }
    
    if ($opts{v}){
        next LINE if $qual < $opts{v};
    }
    #skip if no variant allele in affecteds 
    my %samp_to_gt = VcfReader::getSampleCall
    (
          line              => \@split, 
          sample_to_columns => \%sample_to_col,
          minGQ             => $opts{g},
          multiple          => \@samples,
    );
    next LINE if not haveVariant(\%samp_to_gt);
    
    #skip if not identical GTs in all affecteds and identical variants required
    my @matched_gts = (); 
    if ($opts{e}){
        push @matched_gts, identicalGenotypes(\%samp_to_gt, $opts{n});
        next LINE if not @matched_gts;
    }

    #get Allele frequency annotations if they exist and we're filtering on AF
    my %af_info_values = (); 
    if (%af_info_fields){ 
        %af_info_values = getAfInfoValues(\@split);
    }
    
    # get allele frequencies in  @samples is using -u/--max_sample_allele_frequency
    my %s_allele_counts = ();
    my $allele_count;
    if ($opts{u}){
        %s_allele_counts = VcfReader::countAlleles
        (
            minGQ => $opts{g},
            line  => \@split,
        );
        map { $allele_count += $s_allele_counts{$_} } keys %s_allele_counts;
    }
    #get alleles and consequences for each allele
    my @alleles = VcfReader::readAlleles(line => \@split);
    my @csq = getConsequences(\@split);
    my @alts_to_vep_allele = ();#get VEP style alleles if needed
    if ($opts{m} eq 'vep'){
        @alts_to_vep_allele = VcfReader::altsToVepAllele
        (
            ref => $alleles[0],
            alts => [ @alleles[1..$#alleles] ], 
        );
    }
    #assess each ALT allele (REF is index 0)
    my @cadd_scores = (); 
    if ($opts{c}){
        my $cadd = VcfReader::getVariantInfoField
        (
            \@split,
            'CaddPhredScore',
        );
        if ($cadd){
            @cadd_scores = split(",", $cadd);
        }
    }
    my @matched_samples = (); 
    my %matched_transcript = (); 
ALT: for (my $i = 1; $i < @alleles; $i++){
        push @matched_samples, '.'; #add placeholder for matched INFO field
        next ALT if $alleles[$i] eq '*';
        if (@matched_gts ){
            my $i_match ;
            foreach my $matched_gt( @matched_gts ){
                #if require equal genotypes skip if this is not a matched allele
                my @m = split("/", $matched_gt); 
                if (first { $_ eq $i } @m ){#first is ok cos we never touch allele 0
                    $i_match++; 
                    last; 
                }
            }
            next ALT if not $i_match;
        }
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
        
        #skip if VEP/SNPEFF annotation is not the correct class 
        # (or not predicted damaging if using that option)
        my @samp_with_allele = ();
CSQ:    foreach my $annot (@a_csq){
            if (consequenceMatchesClass($annot, \@split)){
                if (defined $opts{f}){# if looking for matching genes 
                    if (not @samp_with_allele){
                    #only need to add samples once per allele
                        if (@matched_gts){
                            addSamplesWithGt
                            (
                                \@samp_with_allele, 
                                \%samp_to_gt, 
                                \@matched_gts,
                                \@split,
                            );
                        }else{
                            addSamplesWithAllele
                            (
                                \@samp_with_allele, 
                                \%samp_to_gt, 
                                $i,
                                \@split,
                            );
                        }
                        if (@samp_with_allele){
                            pop @matched_samples; #remove '.' placeholder
                            push @matched_samples, join("|", @samp_with_allele);
                        }
                    }
                    if (@samp_with_allele){
                        push @{ $matched_transcript{ $annot->{$feature_id} } }, 
                                                             @samp_with_allele;
                        $transcript_to_symbol{$annot->{$feature_id}} = 
                                                          $annot->{$symbol_id};
                    }
                }else{#if not checking matching genes can print and bail out now
                    print $OUT "$line\n";
                    next LINE; 
                }
            }
        }#end of CSQ
    }#end of ALT
    if (defined $opts{f} and @matched_samples){
        my $line_ref = VcfReader::addVariantInfoField
        (
            line  => \@split,
            id    => "getFunctionalVariantsMatch",
            value => join(",", @matched_samples),
        );
        my $out_line = join("\t", @$line_ref);
        my $var_id = "$chrom:$pos-" . join("/", @alleles ); 
        $vcf_lines{$var_id} = $out_line;
        foreach my $tr (keys %matched_transcript){
            foreach my $sample ( @{ $matched_transcript{$tr} } ){
                push @{ $transcript_sample_vars{$tr}->{ $sample } }, $var_id;
                                       # transcript -> sample -> [var ids]
            }
        }
    }
}#end of LINE
if (defined $opts{f}){
    checkMatchingVariants();
}
close $VCF;

updateProgressBar();  
outputGeneList();
close $OUT;

#################################################
sub addSamplesWithGt{
    my ($samp_ar, $gts, $matched, $line) = @_; 
    foreach my $s ( keys %{$gts} ){
        next if ($gts->{$s} =~ /^\.[\/\|]\.$/ );
        (my $gt = $gts->{$s}) =~ s/[\|]/\//;
          #don't let allele divider affect whether a genotype matches
        if ($opts{pl}){
            if (checkGtPl(\$matched, $s, $line )){
                push @$samp_ar, $s;
            }
        }else{
            if ( first { $_ eq $gt } @$matched ){
                push @$samp_ar, $s ;
            }
        }
    }
}

#################################################
sub addSamplesWithAllele{
    my ($samp_ar, $gts, $allele, $line) = @_; 
    foreach my $s ( keys %{$gts} ){
        if ($gts->{$s} =~ /^($allele[\/\|]\d+|\d+[\/\|]$allele)$/){
            if (checkAllelePl($allele, $s, $line)){
                push @$samp_ar, $s;
            }
        }
    }
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
    }
    %transcript_sample_vars = ();
    %vcf_lines = ();
    %transcript_to_symbol = ();
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
        push @csq_fields, keys %filters;
    }#SnpEff predictions will be added via SnpSift
    return %filters;
}

#################################################
sub openOutput{
    my $OUT_FH;
    my $LIST_FH;
    if ($opts{o}) {
        open( $OUT_FH, ">$opts{o}" ) || die "Can't open $opts{o} for writing: $!\n";
    }
    else {
        $OUT_FH = \*STDOUT;
    }
    if (defined $opts{f}) {
        if ( $opts{f} eq '' ){
          #user specified --list option but provided no argument
            $LIST_FH = \*STDERR;
        }else{
            open( $LIST_FH, ">$opts{f}" )
              or die "Can't open $opts{f} for writing: $!\n";
        }
    }
    return ($OUT_FH, $LIST_FH);        
}

#################################################
sub writeOptionsToHeader{
    #print meta header lines
    print $OUT join("\n", grep { /^##/ } @header) . "\n" ;
    #add header line detailing program options
    print $OUT VcfhacksUtils::getOptsVcfHeader(%opts) . "\n"; 
    print $OUT "$header[-1]\n";
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
SCORE: foreach my $f ( @{ $in_silico_filters{$k} } ) {
            if ( $f =~ /^\d(\.\d+)*$/ ) {
                my $prob;
                if ( $score =~ /\((\d(\.\d+)*)\)/ ) {
                    $prob = $1;
                }else {
                    next SCORE
                      ; #if score not available for this feature ignore and move on
                }
                if ( lc $k eq 'polyphen'){
                    if ( $prob >= $f ){
                    #higher is more damaging for polyphen - damaging
                        return 1 if $opts{k};
                        $filter_matched{$k}++;
                        next PROG;
                    }
                }elsif( $prob <= $f ){
                    #lower is more damaging for sift and condel - damaging
                    return 1 if $opts{k};
                    $filter_matched{$k}++;
                    next PROG;
                }
            }else{
                $score =~ s/\(.*\)//;
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

If -s/--samples argument is specified use this switch to only output variants that make up 'functional' variants in the same genes for the given samples. This will also return a list of genes containing 'functional' variants in these samples. If a filename is passed to this option the list will be printed to file, otherwise the list will be printed to STDERR.

=item B<-n    --num_matching>

If -s/--samples and -f/--find_shared_genes arguments are specified use this option to specify the minimum number of samples with variants in the same gene before outputting variants (and genes). 

=item B<-e    --equal_genotype>

If -s/--samples argument is specified use this flag if you only want to keep variants with identical genotypes in each sample (or a minimum number of samples as specified by the -n/--num_matching option).

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
        splice_region_variant
        stop_lost
        5_prime_UTR_premature_start_codon_gain_variant
        start_lost
        stop_gained
        exon_loss

Available classes that can be chosen instead of (or in addition to - see below) these classes can be found in the data/vep_classes.tsv and data/snpeff_classes.tsv files respectively. 

=item B<-a    --add_classes>

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

=item B<--af    --allele_frequency>

Use a value between 0.00 and 1.00 to specify allele frequencey filtering for annotations from annotateSnps.pl, filterOnEvsMaf.pl or filterVcfOnVcf.pl if you've previously run these programs to annotate your VCF. If an allele frequency is available for an allele it will be filtered if equal to or greater than the value specfied here. 

Note: allele frequencies added by VEP are not used for filtering as they check the allele frequency at the site, not of the specific alleles in your variant call.

=item B<-j    --custom_af>

If using the --af/--allele_frequency option and your data contains allele frequency fields from sources not recognised by this program, you may give the name of these allele frequency INFO fields here and they will be used for filtering in addition to the default fields. Note that these annotations must contain an annotation per ALT allele (i.e. the INFO field header must specify 'Number=A') to work properly and the allele frequency should be expressed as a number between 0.00 and 1.00 in order to be compatible with the default allele frequency fields recognised by this program.

=item B<-u  --max_sample_allele_frequency>

If -s/--samples argument is specified, use this option to specify an allele frequency (between 0.00 and 1.00) for filtering alleles. Alleles present at this frequency or higher in your samples of interest will be filtered.

=item B<-v    --var_qual>

Minimum variant Phred-like quality score to consider. Variants with a QUAL field lower than this value will be filtered. Default is 20.

=item B<-g    --gq>

Minimum genotype qualities to consider. Only applicable when using the -s/--samples option. Any genotype call below this threshold will be considered a no call. Default is 20

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

    TODO!

=cut

=head1 AUTHOR

David A. Parry

University of Edinburgh

=head1 COPYRIGHT AND LICENSE

Copyright 2015, David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut


