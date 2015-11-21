#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Term::ProgressBar;
use Data::Dumper;
use List::Util qw ( first ) ;
use POSIX qw/strftime/;
use FindBin;
use lib "$FindBin::Bin/lib";
use ParsePedfile;
use VcfReader;

my $data_dir = "$FindBin::Bin/data";
my $progressbar;
my @samples = ();
my @reject = ();
my @reject_except = ();
my @classes = (); 
my @add_classes = ();
my @damaging = ();
my @biotypes = ();

my %opts = (
    s               => \@samples,
    r               => \@reject,
    x               => \@reject_except,
    classes         => \@classes,
    add_classes     => \@add_classes,
    d               => \@damaging,
    biotype_filters => \@biotypes,
);


GetOptions(
    \%opts,
    'i|input=s' ,
    'o|output=s',
    'l|list:s', 
    's|samples=s{,}',  
    'r|reject=s{,}',    
    'x|reject_all_except:s{,}',
    'f|family=s', # ped file
    'm|mode=s',
    'classes=s{,}',
    'add_classes=s{,}',
    'd|damaging=s{,}',
    'canonical_only',
    'biotype_filters=s{,}',
    'no_biotype_filtering',
    'pass_filters',
    'keep_any_damaging',
    'skip_unpredicted_missense',
    'a|af=f',
    'x_linked=i', #1 = look for x-linked recessive only, 2= look for x-linked recessive as well
    'check_all_samples',
    'e|equal_genotypes',
    'q|quality=i',
    'w|aff_quality=i',
    'u|un_quality=i',
    'n|num_matching=i',
    'y|num_matching_per_family=i',
    't|ignore_carrier_status',
    'b|progress',
    'z|homozygous_only',
    'consensus_splice_site',
    'lenient',
    'h|?|help',
    'manual'
) or pod2usage( -message => "SYNTAX ERROR", exitval => 2 );

pod2usage( -verbose => 2 ) if $opts{manual};
pod2usage( -verbose => 1 ) if $opts{h};

######SANITY CHECKS######
pod2usage( -message => "SYNTAX ERROR: input is required.", exitval => 2 )
  if not $opts{i};
pod2usage(
    -message =>
"SYNTAX ERROR: please specify samples to analyze using --samples (-s), --check_all_samples or --family (-f) arguments.",
    exitval => 2
  )
  if not @samples
  and not $opts{check_all_samples}
  and not $opts{f};

#af is >0 and <= 0.5
pod2usage(
    -message =>
"SYNTAX ERROR: --af option requires a value between 0.00 and 0.50 to filter on global minor allele frequency.\n",
    -exitval => 2
) if ( defined $opts{a} && ( $opts{a} < 0 or $opts{a} > 0.5 ) );

#GQs are >= 0
$opts{q} = defined $opts{q} ? $opts{q} : 20;
pod2usage(
    -message => "SYNTAX ERROR: Genotype quality scores must be 0 or greater.\n",
    -exitval => 2
) if ( $opts{q} < 0 );

if ( defined $opts{w} ) {
    pod2usage(
        -message => "SYNTAX ERROR: Genotype quality scores must be 0 or greater.\n",
        -exitval => 2
    ) if $opts{w} < 0;
}else {
    $opts{w} = $opts{q};
}

if ( defined $opts{u} ) {
    pod2usage(
        -message => "SYNTAX ERROR: Genotype quality scores must be 0 or greater.\n",
        -exitval => 2
    ) if $opts{u} < 0;
}
else {
    $opts{u} = $opts{q};
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


#check VEP/SNPEFF header and get annotation order - will set $opts{m} if not already defined
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

#check ped, if provided is ok and get specified families
my ($ped_obj, @fams) = checkPedAndFamilies();

#add all samples if --check_all_samples option is selected

if ($opts{check_all_samples}) {
    push @samples, keys %sample_to_col;
}

# check affected and unaffected samples from ped and specified on commandline

addSamplesFromPed();
#store ped samples in a hash so we can quickly look up whether we should be checking family IDs
my %ped_samples = ();
if ($opts{f}){
    %ped_samples = map { $_ => undef } getPedSamples();
}
checkSamples(); 

my $min_biallelic = checkMinMatching();

# get available allele frequency annotations if filtering on AF
# NOTE: we do not use VEP annotations as they do not necessarily match the variant allele

my %af_info_fields = (); #check available INFO fields for AF annotations
if ( defined $opts{a} ) {
    my %info_fields = VcfReader::getInfoFields(header => \@header);
    %af_info_fields = getAfAnnotations(\%info_fields);   
}
######PRE-PROCESSING######

#Get filehandle for output (STDOUT or file) and optionally gene list (STDERR or file)
my ($OUT, $LIST ) = openOutput();

#write header and program run parameters
writeOptionsToHeader();

#start progress bar 
my $next_update = 0;
my $total_vars  = 0; 
my $var_count   = 0;
if ($opts{b}) {
    if ( $opts{i} eq "-" ) {
        informUser("Can't use --progress option when input is from STDIN\n");
    }else{
        my $msg = <<EOT
WARNING: -p/--progress option requires installation of the perl module 'Term::ProgressBar'.
WARNING: 'Term::ProgressBar' module was not found. No progress will be shown.
EOT
;
        eval "use Term::ProgressBar; 1 " or informUser($msg);
        if (not $@){
            informUser("Counting variants in input for progress monitoring.\n"); 
            $total_vars = VcfReader::countVariants($opts{i});
            informUser("INFO: $opts{i} has $total_vars variants.\n");
            $progressbar = Term::ProgressBar->new(
                {
                    name  => "Biallelic",
                    count => $total_vars,
                    ETA   => "linear",
                }
            );
        }
    }
}

######READ VARIANTS AND PROCESS######

my %contigs = (); 
 #keep track of contigs we've seen so we can check input is (roughly) sorted
my %transcript_vars = ();
 #stores variant IDs meeting our criteria per (key) transcript 
my %sample_vars = (); 
 #hash of variant IDs to sample_id to allele counts
my %vcf_lines = ();
 #hash of variant IDs to matching VCF line
my %transcript_to_symbol = ();
 #transcript ID to gene symbol
my %gene_listing = ();
 #gene symbols => transcripts with biallelic variants
open (my $VCF, $opts{i}) or die "Cannot open $opts{i} for processing variants: $!\n";

#read line
LINE: while (my $line = <$VCF>){
    next if $line =~ /^#/;#skip header
    chomp $line;
    $var_count++;#for progress bar
    my @split = split("\t", $line); 
    my ($chrom, $pos, $filter) = VcfReader::getMultipleVariantFields
    (
        \@split, 
        'CHROM', 
        'POS', 
        'FILTER',
    );
    #check if new chromosome - 
    if (exists $contigs{$chrom}){
        #require chroms to be ordered together
        if ($contigs{$chrom} != scalar(keys %contigs) - 1){
            die
"Encountered a variant for contig $chrom with at least one variant for a different contig inbetween. Your VCF input must be sorted such that all contigs are kept together - please sort your input and try again.\n";
        }
    }else{
        #if new chrom check biallelics and clear collected data
        $contigs{$chrom} = scalar(keys%contigs);
        checkBiallelic() if $contigs{$chrom} > 0;
    }
    
    #skip if FILTER != PASS and PASS required
    if ($opts{pass}){
        next LINE if $filter ne 'PASS';
    }
    
    #check that this is an autosome (or X if we're interested in X-linked)
    next LINE if not checkChromosome($chrom); 
    
    #skip if no variant allele in affecteds 
    my %samp_to_gt = VcfReader::getSampleCall
    (
          line              => \@split, 
          sample_to_columns => \%sample_to_col,
          minGQ             => $opts{w},
          multiple          => \@samples,
    );
    next LINE if not haveVariant(\%samp_to_gt);
    
    #skip if not identical GTs in all affecteds and identical variants required
    if ($opts{e}){
        next LINE if not identicalGenotypes(\%samp_to_gt, $opts{n});
    }

    #collect genotypes from @reject samples so we can skip homozygous alleles
    #and add them to our per transcript hash of alleles
    my %reject_to_gt = VcfReader::getSampleCall
    (
          line              => \@split, 
          sample_to_columns => \%sample_to_col,
          minGQ             => $opts{u},
          multiple          => \@reject,
    );
    
    #get Allele frequency annotations if they exist and we're filtering on AF
    my %af_info_values = (); 
    if (%af_info_fields){ 
        %af_info_values = getAfInfoValues(\@split);
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
ALT: for (my $i = 1; $i < @alleles; $i++){
        next ALT if $alleles[$i] eq '*';
        #skip if homozygous in any unaffected
        next ALT if 
        (  
            grep { /^$i[\/\|]$i$/ } 
            map { $reject_to_gt{$_} }
            keys %reject_to_gt
        );
    
        #filter if AF annotations > MAF
        if (%af_info_fields){ #check for annotateSnps.pl etc. frequencies
           next ALT if ( alleleAboveMaf($i - 1, \%af_info_values) );
        }
    
        #get all consequences for current allele
        my @a_csq = ();
        if ($opts{m} eq 'vep'){
            @a_csq = grep { $_->{allele} eq $alts_to_vep_allele[$i-1] } @csq;
        }else{
            @a_csq = grep { $_->{allele} eq $alleles[$i] } @csq;
        }
        
        #skip if VEP/SNPEFF annotation is not the correct class (or not predicted damaging if using that option)
        foreach my $annot (@a_csq){
            if (consequenceMatchesClass($annot)){
                #create variant ID for this allele
                my $var_id = "$chrom:$pos-$i";
                #add variant ID to transcript variants
                push @{$transcript_vars{$annot->{$feature_id}}}, $var_id;
                #get allele counts per sample for this variant ID if we haven't already
                if (not exists $sample_vars{$var_id}){
                    $sample_vars{$var_id} = getSampleAlleleCounts
                    (
                        {%samp_to_gt, %reject_to_gt}, 
                        $i,
                    );
                }
                #store VCF line for this variant ID for printing out later
                if (not exists $vcf_lines{$var_id}){
                    $vcf_lines{$var_id} = \@split;
                }
                if (not exists $transcript_to_symbol{$annot->{$feature_id}}){
                    $transcript_to_symbol{$annot->{$feature_id}} = $annot->{$symbol_id};
                }
            }
        }
    }
    updateProgressBar();  
}
close $VCF;

checkBiallelic();
close $OUT;

outputGeneList();


#################################################
sub getSampleAlleleCounts{
    my ($gts, $i) = @_;
    my %counts = ();
    foreach my $s (keys %$gts){
        if ( $gts->{$s} =~ /^$i[\/\|]$i$/ ) {
            $counts{$s} = 2;    #homozygous for this alt allele
        }elsif ( $gts->{$s} =~ /^([\d+\.][\/\|]$i|$i[\/\|][\d+\.])$/ ) {
            $counts{$s} = 1;    #het for alt allele
        }elsif ( $gts->{$s} =~ /^\.[\/\|]\.$/ ) {
            $counts{$s} = -1;    #no call
        }else {
            $counts{$s} = 0;    #does not carry alt allele
        }
    }
    return \%counts;
}

#################################################
sub checkBiallelic{
#For each transcript (i.e. each key of our hash of transcript IDs to arrays of var IDs)
    my %lines_to_write = (); #collect vcf lines to write out as keys...
                            #...with values being a hash of hom and het samples like so:
                #$lines_to_write{vcfline}->{alt_num}->{2}->{sample} = n  [hom samples] 
                #$lines_to_write{vcfline}->{alt_num}->{1}->{sample} = n  [het samples] 
    foreach my $k (keys %transcript_vars) { 
        if (parseAlleles($transcript_vars{$k}, \%lines_to_write)){
            push @{ $gene_listing{$transcript_to_symbol{$k}} }, $k if $LIST;
        }
    }
    #get VCF lines with biallelic variants
    my @lines_out = ();
    foreach my $v (keys %lines_to_write){
        #add findBiallelicSamplesHom and findBiallelicSamplesHet INFO fields
        my @het_string = ();
        my @hom_string = ();
        my @sv = split("\t", $v);
        my @alleles = VcfReader::readAlleles(line => \@sv); 
        for (my $alt = 1; $alt < @alleles; $alt++){ 
            if (exists $lines_to_write{$v}->{$alt}->{1}){
                push @het_string, join("|", keys %{$lines_to_write{$v}->{$alt}->{1}});
            }else{
                push @het_string, ".";
            }
            if (exists $lines_to_write{$v}->{$alt}->{2}){
                push @hom_string, join("|", keys %{$lines_to_write{$v}->{$alt}->{2}});
            }else{
                push @hom_string, ".";
            }
        }
        
        my $line_ref = VcfReader::addVariantInfoField
        (
            line  => \@sv,
            id    => "findBiallelicSamplesHom",
            value => join(",", @hom_string),
        );
        $line_ref = VcfReader::addVariantInfoField
        (
            line  => \@sv,
            id    => "findBiallelicSamplesHet",
            value => join(",", @het_string),
        );
        push @lines_out, $line_ref;
    }
    @lines_out = VcfReader::sortByPos(\@lines_out); 
    foreach my $l (@lines_out){
        print $OUT join("\t", @$l) . "\n";
    }
    #clear collected data before moving on to next chromosome
    %transcript_vars = ();
    %sample_vars = (); 
    %vcf_lines = ();
    %transcript_to_symbol = ();
}

#################################################
sub parseAlleles{
    my $vars = shift;
    my $output_hash = shift;
    my %incompatible = (); #key is genotype can't be a pathogenic compound het combo
    my %biallelic = (); #key is sample, value is array of (putative) biallelic genotypes
    #foreach allele
VAR: for (my $i = 0; $i < @$vars; $i ++){ 
        #check unaffected genotypes...
        foreach my $r (@reject){
            #... and skip VAR if homozygous
            next VAR if $sample_vars{$vars->[$i]}->{$r} == 2;
            if ($sample_vars{$vars->[$i]}->{$r} == 1){#reject sample is het...
            #... find other alleles in same sample - these can't make up biallelic combos
                for (my $j = $i + 1; $j < @$vars; $j++){
                    if ($sample_vars{$vars->[$j]}->{$r} == 1){
                        $incompatible{"$vars->[$i]/$vars->[$j]"}++;
                    }
                }
            }
        }
        @{ $incompatible{$vars->[$i]} } = removeDups(@{ $incompatible{$vars->[$i]} });
        #...get combinations of biallelic alleles for each affected
SAMPLE: foreach my $s (@samples){
            next SAMPLE if $sample_vars{$vars->[$i]}->{$s} < 1;
            if ($sample_vars{$vars->[$i]}->{$s} == 2){
                push @{ $biallelic{$s} } , "$vars->[$i]/$vars->[$i]";
            }
            next SAMPLE if $opts{z};#move on if -z/--homozygous only
            if ($sample_vars{$vars->[$i]}->{$s} >= 1){
                for (my $j = $i + 1; $j < @$vars; $j++){
                    if (not exists $incompatible{"$vars->[$i]/$vars->[$j]"}){
                    #...check phase info for any potential compound hets (if available)
                        if (! allelesInCis($vars->[$i], $vars->[$j], $s) ){
                            push @{ $biallelic{$s} } , "$vars->[$i]/$vars->[$j]";
                        }
                    }
                }
            }
        }
    }
    
    #%biallelic now has a key for every sample with a potential biallelic genotype
    #Test no. samples with biallelic genotypes against min. no. required matching 
    # samples or if not specified all samples
    #bail out straight away if no. samples less than required;
    return if keys %biallelic < $min_biallelic;
    if ($ped_obj){ 
        #check segregation of putative biallelic combinations 
        #...for each family and remove those that do not segregate/not
        #...present in all/min affecteds
        checkSegregation(\%biallelic);
    } 
    my $matches = 0;
    my %counted_fam = ();
    foreach my $s (keys %biallelic){
        #count each family only once if we have a ped file
        if (exists $ped_samples{$s}){
            my $f_id = $ped_obj->getFamilyId($s);
            $matches++ if not exists $counted_fam{$f_id};
            $counted_fam{$f_id}++;
        }else{
            $matches++;
        }
        last if $matches >= $min_biallelic;
    }
    return if $matches < $min_biallelic;

    # get lines for each allele of biallelic genotypes
    # and add INFO field for samples with biallelic variants
    foreach my $s (keys %biallelic){
        foreach my $gt ( @{ $biallelic{$s} } ){
            my @alleles = split(/\//, $gt);
            foreach my $al (@alleles){
                #add information to output hash reference 
                my $alt = (split "-", $al)[1];#ALT allele number from VCF line
                my $line = join("\t", @{$vcf_lines{$al}}); 
                $output_hash->{$line}->{$alt}->{$sample_vars{$al}->{$s}}->{$s}++;
                #$output_hash->{vcfline}->{alt_num}->{2}->{sample} = n  [hom samples] 
                #$output_hash->{vcfline}->{alt_num}->{1}->{sample} = n  [het samples] 
            }
        }
    }
    return 1;#transcript has biallelic variants
}

#################################################
sub checkSegregation{
    my $bial = shift;#hash of samples to biallelic genotypes
    #each genotype represented as "$chrom:$pos1-$i/$chrom:$pos2-$j";

    #for each biallelic genotype check all/min affecteds in each family carry genotype
    foreach my $f ( getPedFamilies() ) {
        my @checked_gts = ();#genotypes that pass seg test for each family go here
        my @f_aff = grep { exists $sample_to_col{$_} } $ped_obj->getAffectedsFromFamily($f);
        my @f_unaff = grep { exists $sample_to_col{$_} } $ped_obj->getUnaffectedsFromFamily($f);
        my %gt_counts = ();
        foreach my $s (@f_aff){
            foreach my $gt (@{ $bial->{$s} }){
                push @{$gt_counts{$gt}}, $s;
            }
        }
        my $min_matches = $opts{y} ? $opts{y} : scalar(@f_aff);
GENOTYPE: foreach my $gt (keys %gt_counts){
            next if @{$gt_counts{$gt}} < $min_matches; 
            my @alleles = split(/\//, $gt);
            die "ERROR: Got " . scalar(@alleles) . " alleles in genotype counts!" 
              if @alleles != 2; 
            #we should have already rejected any combinations carried by any unaffected
            #...found in the VCF, so now we just need to check obligate carriers
            my ($chrom) = (split /\:/, $alleles[0]); 
            if (not $opts{y} and not $opts{t} and $chrom !~ /^(chr)*X/i){
                foreach my $aff ( $ped_obj->getAffectedsFromFamily($f) ){
                    my @par =
                      grep { exists $sample_to_col{$_} } $ped_obj->getParents($aff);
                # if no min_matching_per_family and no ignore_carrier_status
                #...assume we are not entertaining any phenocopies and therefore
                #...all obligate carriers must carry mutation
                    foreach my $un (@par) {
                        if ($sample_vars{$alleles[0]}->{$un} == 0 and
                            $sample_vars{$alleles[1]}->{$un} == 0){
                            #parent does not carry either allele
                            next GENOTYPE;
                        }
                    }
                } 
            }#if we haven't skipped then genotype should be considered valid
            push @checked_gts, $gt;    
        }#end of each genotype
        foreach my $s (@f_aff){
            delete  $bial->{$s};    
            foreach my $gt (@checked_gts){
                if ( grep {$s eq $_ } @{$gt_counts{$gt}}){
                    push @{$bial->{$s}}, $gt;
                }
            }
        }
    }#end of each family
}
#################################################
sub allelesInCis{
#returns 0 if alleles are in trans or if phase unknown
#returns 1 if alleles in cis
#alleles to check ($var1 and $var2) are in format "$chrom:$pos-$i";
    my ($var1, $var2, $sample) = @_;
    #get corresponding split VCF line for allele
    my ($line1, $line2) = ($vcf_lines{$var1}, $vcf_lines{$var2});
    my %format1 = VcfReader::getVariantFormatFields($line1);
    my %format2 = VcfReader::getVariantFormatFields($line2);
    my %phase1 = ();
    my %phase2 = ();
    foreach my $f (qw / PID PGT / ){
        return 0 if not exists $format1{$f};
        return 0 if not exists $format2{$f};
        $phase1{$f} = VcfReader::getSampleGenotypeField
        (
           line => $line1, 
           field => $f, 
           sample => $sample, 
           sample_to_columns => \%sample_to_col
        );
        $phase2{$f} = VcfReader::getSampleGenotypeField
        (
           line => $line2, 
           field => $f, 
           sample => $sample, 
           sample_to_columns => \%sample_to_col
        );
    }
    if ($phase1{PID} ne $phase2{PID}){
        return 0;
    }
    my $al1 = (split "-", $var1)[1];#ALT allele number from VCF line
    my $al2 = (split "-", $var2)[1];
    return 0 if $phase1{PGT} eq '.' or $phase2{PGT} eq '.';
    my @pgt1 = split(/\|/, $phase1{PGT});
    my @pgt2 = split(/\|/, $phase2{PGT});
    my $p1 = first { $pgt1[$_] eq $al1 } 0 .. $#pgt1;
    my $p2 = first { $pgt2[$_] eq $al2 } 0 .. $#pgt2;
    return 0 if not defined $p1 or not defined $p2;
    return 0 if $p1 != $p2;
    return 1;
}

#################################################
sub readClassesFile{
    my $classes_file = "$data_dir/$opts{m}_classes.tsv";
    open (my $CLASS, $classes_file) or die 
"Could not open $opts{m} effect classes file '$classes_file': $!\n";
    my %classes = (); 
    while (my $line = <$CLASS>){
        next if $line =~ /^#/;
        $line =~ s/[\r\n]//g; 
        next if not $line;
        my @split = split("\t", $line);
        die "Not enough fields in classes file line: $line\n" if @split < 2;
        $classes{lc($split[0])} = lc($split[1]);
    }
    close $CLASS;
    return %classes;
}

#################################################
sub readBiotypesFile{
    my $btypes_file = "$data_dir/biotypes.tsv";
    open (my $BT, $btypes_file) or die 
"Could not open biotypes file '$btypes_file': $!\n";
    my %biotypes = (); 
    while (my $line = <$BT>){
        next if $line =~ /^#/;
        $line =~ s/[\r\n]//g; 
        next if not $line;
        my @split = split("\t", $line);
        die "Not enough fields in biotypes file line: $line\n" if @split < 2;
        $biotypes{lc($split[0])} = lc($split[1]);
    }
    close $BT;
    return %biotypes;
}

#################################################
sub getAndCheckClasses{
    my %all_classes = readClassesFile();
    if (not @classes){
        @classes = grep { $all_classes{$_} eq 'default' } keys %all_classes;
    }
    push @classes, @add_classes if (@add_classes);
    if ($opts{m} eq 'vep' and $opts{consensus_splice_site}){
        push @classes, "splice_region_variant" ;
    }
    @classes = map { lc($_) } @classes; 
    @classes = removeDups(@classes);
    foreach my $class (@classes) {
        die "Error - variant class '$class' not recognised.\n"
          if not exists $all_classes{lc($class)} ;
    }
    return @classes;
}

#################################################
sub getAndCheckBiotypes{
    return if $opts{no_biotype_filtering}; 
    my %all_biotypes = readBiotypesFile();
    if (not @biotypes){
        @biotypes = grep { $all_biotypes{$_} eq 'filter' } keys %all_biotypes;
    }
    @biotypes = map { lc($_) } @biotypes; 
    @biotypes = removeDups(@biotypes);
    foreach my $biotyp (@biotypes) {
        die "Error - biotype '$biotyp' not recognised.\n"
          if not exists $all_biotypes{lc($biotyp)} ;
    }
    return @biotypes;
}
    
#################################################
sub getAndCheckInSilicoPred{
    return if not @damaging;
    my %pred = readInSilicoFile(); 
    my %filters = (); 
    foreach my $d (@damaging) {
        my ( $prog, $label ) = split( "=", $d );
        if ($prog eq 'all'){
            foreach my $k (keys %pred){
                foreach my $j ( keys %{$pred{$k}} ){
                    push @{ $filters{$k} }, $j if $pred{$k}->{$j} eq 'default';
                }
            }
        }elsif ( exists $pred{$prog} ){
            if (not $label){
                foreach my $j ( keys %{$pred{$prog}} ){
                    push @{ $filters{$prog} }, $j if $pred{$prog}->{$j} eq 'default';
                }
            }else{
                if ($label =~ /^\d+(\.\d+)*$/){#add scores
                    push @{ $filters{$prog} } , $label;
                }else{
                    if (not exists $pred{$prog}->{lc($label)}){
                        die <<EOT;
ERROR: Unrecognised filter ('$label') for $prog passed to --damaging argument.
See --help/--manual for more info.
EOT
;
                    }
                    push @{ $filters{$prog} } , lc($label);
                }
            }
        }else{
            die <<EOT
ERROR: Unrecognised value ($d) passed to --damaging argument. 
See --help/--manual for more info.
EOT
;
        }
    }
    if ($opts{m} eq 'vep'){#VEP prediction results will be in CSQ field
        push @csq_fields, keys %filters;
    }#SnpEff predictions will be added via SnpSift
    return %filters;
}
    
#################################################
sub readInSilicoFile{
    my $insilico_file = "$data_dir/$opts{m}_insilico_pred.tsv";
    open (my $INS, $insilico_file) or die 
"Could not open $opts{m} effect classes file '$insilico_file': $!\n";
    my %pred = (); 
    while (my $line = <$INS>){
        next if $line =~ /^#/;
        $line =~ s/[\r\n]//g; 
        next if not $line;
        my @split = split("\t", $line);
        die "Not enough fields in classes file line: $line\n" if @split < 3;
        $pred{lc($split[0])}->{lc($split[1])} = lc($split[2]) ;
    }
    close $INS;
    return %pred;
}

#################################################
sub checkPedAndFamilies{
    return if not $opts{f}; 
    my $ped;
    my @fams = (); #any comma separated fields after ped filename = families to analyze
    @fams = split(",", $opts{f}); 
    my $pedfile = shift @fams; 
    $ped = ParsePedfile->new( file => $pedfile );
    if (@fams){#if families specified check they exist in PED
        my %all_fam = map { $_ => undef } $ped->getAllFamilies();
        foreach my $f (@fams){
            die "Error - User specified family \"$f\" does not exist in PED file!\n"
              if not exists $all_fam{$f}; 
        }
    }
    return ($ped, @fams); 
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
                die <<EOT 
ERROR: Could not find VEP or SnpEff headers in input. Please annotate your input with either program and try again.
EOT
;
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
sub getPedFamilies{
    if (@fams){
        return @fams;
    }else{
        return $ped_obj->getAllFamilies();
    }
}

#################################################
sub getPedSamples{
    my @samples = ();
    if (@fams){
        foreach my $f (@fams){
            push @samples, $ped_obj->getSamplesFromFamily($f);
        }
    }else{
        push @samples, $ped_obj->getAllSamples();       
    }
    return @samples;
}

#################################################
sub getPedAffecteds{
    my @affecteds = ();
    if (@fams){
        foreach my $f (@fams){
            push @affecteds, $ped_obj->getAffectedsFromFamily($f);
        }
    }else{
        push @affecteds, $ped_obj->getAllAffecteds();       
    }
    return @affecteds;
}

#################################################
sub getPedUnaffecteds{
    my @unaffecteds = ();
    if (@fams){
        foreach my $f (@fams){
            push @unaffecteds, $ped_obj->getUnaffectedsFromFamily($f);
        }
    }else{
        push @unaffecteds, $ped_obj->getAllUnaffecteds();       
    }
    return @unaffecteds;
}

#################################################
sub addSamplesFromPed{
    return if not $ped_obj;
    my @aff     = ();
    my @un      = ();
    my @not_aff = ();
    my @not_un  = ();
    foreach my $s ( getPedAffecteds() ) {
        if ( exists $sample_to_col{$s} ) {
            push @aff, $s;
        }
        else {
            push @not_aff, $s;
        }
    }
    foreach my $s ( getPedUnaffecteds() ) {
        if ( exists $sample_to_col{$s} ) {
            push @un, $s;
        }
        else {
            push @not_un, $s;
        }
    }
    informUser( "Found "
      . scalar(@aff)
      . " affected samples from pedigree in VCF.\n");
    informUser( scalar(@not_aff)
      . " affected samples from pedigree were not in VCF.\n");
    informUser( "Found "
      . scalar(@un)
      . " unaffected samples from pedigree in VCF.\n");
    informUser( scalar(@not_un)
      . " unaffected samples from pedigree were not in VCF.\n");
    if ($opts{y}) {
        if ( not $opts{n} ) {
            informUser(
"WARNING: --num_matching_per_family (-y) set without setting --num_matching (-n). "
              . "Setting --num_matching to $opts{y}.\n");
            $opts{n} = $opts{y};
        }
        foreach my $f ( getPedFamilies() ) {
            my %affected_ped =
              map { $_ => undef } $ped_obj->getAffectedsFromFamily($f);
            my @affected = grep { exists $affected_ped{$_} } @aff;
            die "Number of affected found in family $f ("
              . scalar(@affected)
              . ") is less than value for "
              . "--num_matching_per_family ($opts{y})\n"
              if $opts{y} > @affected;
        }
    }
    push @samples, @aff;
    push @reject,  @un;
   
}

#################################################
sub checkSamples{
    #check @samples and @reject and @reject_except exist in file
    my @not_found = ();
    foreach my $s ( @samples, @reject, @reject_except ) {
        #reject_except may contain a single empty entry
        next if length($s) == 0;
        if ( not exists $sample_to_col{$s} ) {
            push @not_found, $s;
        }
    }
    if (@not_found) {
        die "Could not find the following samples in VCF:\n"
          . join( "\n", @not_found ) . "\n";
    }
    # at least one affected sample in @sample and 
    die "No affected samples found in VCF\n" if not @samples;

    #parse our --reject_except/-x samples
    if (@reject_except) {
        my @all = keys %sample_to_col; 
        push @reject_except, @samples;
        my %subtract = map { $_ => undef } @reject_except;
        @all = grep { !exists $subtract{$_} } @all;
        push @reject, @all;
    }

    #remove any duplicate samples in @samples or @reject
    @reject  = removeDups(@reject); 
    @samples = removeDups(@samples);
    #make sure no samples appear in both @samples and @reject
    my %dup = map { $_ => undef } @samples;
    foreach my $s (@samples) {
        my @doubles = grep { exists( $dup{$_} ) } @reject;
        die "Same sample(s) specified as both affected and unaffected:\n"
          . join( "\n", @doubles ) . "\n"
          if @doubles;
    }
}

#################################################
sub checkMinMatching{
    #'n|num_matching=i',
    #'y|num_matching_per_family=i',
    
    if (not defined $opts{n} and not defined $opts{y}){
        return calculateMaxMatching();
    }
    if (defined $opts{n}){
        if ($opts{n} < 1){
            die "ERROR: -n/--num_matching argument must be greater than 0!\n";
        }
    }
    if (defined $opts{y}){
        if ($opts{y} < 1){
            die "ERROR: -y/--num_matching_per_family argument must be greater than 0!\n";
        }
    }
    
    if ($opts{f} and not $opts{y}) {
        informUser("NOTE: --num_matching (-n) argument should be used in "
          . "conjunction with --num_matching_per_family (-y) when using a PED file "
          . "if you want to allow for missing variants within the same family.\n");
    }

    #min matching settings are not greater than no. affecteds
    if ($opts{n} > @samples){
        die
"ERROR: -n/--num_matching value [$opts{n}] is greater than number of --samples identified ("
      . scalar(@samples) . ")\n";
    }
    #min matching settings are not greater than no. families if using ped
    if ($opts{f}) {
        my $aff_count =  calculateMaxMatching(); 
        if ($opts{n} > $aff_count){
            die
"ERROR: -n/--num_matching value ($opts{n}) is greater than the number of families "
      . "with affected members identified in ". $ped_obj->{get_file} . " and --samples identified ($aff_count)\n";
        }
    }
    return $opts{n};
}

#################################################
sub calculateMaxMatching{
    #only count one sample per family for $min_matching_samples
    if ($opts{f}){
        my $aff_count = 0;
        foreach my $f ( getPedFamilies() ) {
            $aff_count++ if ( $ped_obj->getAffectedsFromFamily($f) );
        }
        foreach my $s (@samples) {
            $aff_count++ if not grep { $_ eq $s } getPedSamples();
        }
        return $aff_count;
    }else{
        return scalar @samples;
    }
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
    if (defined $opts{l}) {
        if ( $opts{l} eq '' ){
          #user specified --list option but provided no argument
            $LIST_FH = \*STDERR;
        }else{
            open( $LIST_FH, ">$opts{l}" )
              or die "Can't open $opts{l} for writing: $!\n";
        }
    }
    return ($OUT_FH, $LIST_FH);        
}

#################################################
sub writeOptionsToHeader{
    my @opt_string = ();
    #get options into an array
    foreach my $k ( sort keys %opts ) {
        if ( not ref $opts{$k} ) {
            push @opt_string, "$k=$opts{$k}";
        }
        elsif ( ref $opts{$k} eq 'SCALAR' ) {
            if ( defined ${ $opts{$k} } ) {
                push @opt_string, "$k=${$opts{$k}}";
            }
            else {
                push @opt_string, "$k=undef";
            }
        }
        elsif ( ref $opts{$k} eq 'ARRAY' ) {
            if ( @{ $opts{$k} } ) {
                push @opt_string, "$k=" . join( ",", @{ $opts{$k} } );
            }
            else {
                push @opt_string, "$k=undef";
            }
        }
    }
    print $OUT join("\n", grep { /^##/ } @header);
    print $OUT "##findBiallelicVep.pl\"";
    print $OUT join( " ", @opt_string ) . "\"\n";
    
    #add INFO fields 
    print $OUT
    "##INFO=<ID=findBiallelicSamplesHom,Number=A,Type=String,Description=\"For each allele a list of samples that were found to meet findBiallelic.pl's criteria for constituting a homozygous variant.\">\n";
    print $OUT
    "##INFO=<ID=findBiallelicSamplesHet,Number=A,Type=String,Description=\"For each allele a list of samples that were found to meet findBiallelic.pl's criteria for constituting a potential compound heterozygous variant.\">\n";
    print $OUT "$header[-1]\n";
}

#################################################
sub checkChromosome{
    my $chrom = shift;
    if ( not defined $opts{x_linked} ){
        return 0 if not is_autosome($chrom);
    }elsif ($opts{x_linked} == 0){
        return 0 if not is_autosome($chrom);
    }elsif ($opts{x_linked} == 1){
        return 0 if $chrom !~ /^(chr)*X$/i;
    }elsif ($opts{x_linked} == 2){
        if (not is_autosome($chrom)){
            return 0 if $chrom !~ /^(chr)*X$/i;
        }
    }
    return 1;
}

#################################################
sub haveVariant{
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
    return 1 if $min_matching < 2;#if only one sample needs to match then no need to check
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
        my %family_matches = (); 
        my $m = 0;
        foreach my $s (@{$matches{$k}}){
            if ($opts{f}){
                if (exists $ped_samples{$s}){
                    my $f_id = $ped_obj->getFamilyId($s);
                    $family_matches{$f_id}++;
                }else{
                    $m++;
                }
            }else{
                $m++;
            }
        }
        foreach my $f (keys %family_matches){#only count families once 
            if ($opts{y}){#min matching per family
                $m++ if ($family_matches{$f} >= $opts{y});
            }else{#otherwise need all affecteds from family in VCF to match
                my @aff, $ped_obj->getAffectedsFromFamily($f);
                @aff = grep { exists $gts->{$_} } @aff;
                $m++ if $family_matches{$f} >= scalar @aff; 
            }
        }
        return 1 if $m >= $min_matching - 1;
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
    foreach my $key (keys %$info_fields){
        my $warning = <<EOT
[WARNING] Found expected frequency annotation ($key) in INFO fields, but 'Number' field is $info_fields->{$key}->{Number}, expected 'A'. Ignoring this field.
EOT
;
        my $info = <<EOT
[INFO] Found allele frequency annotation: $key. This will be used for filtering on allele frequency.
EOT
;
        if (grep { $key eq $_ } @af_fields){
            if ($info_fields->{$key}->{Number} ne 'A'){
                print STDERR $warning;
            }else{
                print STDERR $info;
                $af_found{$key} = $info_fields->{$key};
            }
        }else{
            if ($key =~ /^FVOV_AF_\S+$/){
                if ($info_fields->{$key}->{Number} ne 'A'){
                    print STDERR $warning;
                }else{
                    print STDERR $info;
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
    if ($opts{m}){
        return consequenceMatchesVepClass($annot);
    }else{
        return consequenceMatchesSnpEffClass($annot);
    }   
}

#################################################
sub consequenceMatchesVepClass{
    my $annot = shift;
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
    #skip unwanted biotypes
    return 0 if exists $biotype_filters{lc $annot->{transcript_biotype} };
    my @anno_csq = split( /\&/, $annot->{annotation} );
ANNO: foreach my $ac (@anno_csq){
        $ac = lc($ac);#we've already converted %class_filters to all lowercase
        if ( exists $class_filters{$ac} ){
            if ($ac eq 'missense_variant' and %in_silico_filters){
                return 1 if damagingMissenseSnpEff($ac); 
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
        #unless skip_unpredicted_missense is in effect
            $filter_matched{$k}++ unless $opts{skip_unpredicted_missense};
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
                        return 1 if $opts{keep_any_damaging};
                        $filter_matched{$k}++;
                        next PROG;
                    }
                }elsif( $prob <= $f ){
                    #lower is more damaging for sift and condel - damaging
                    return 1 if $opts{keep_any_damaging};
                    $filter_matched{$k}++;
                    next PROG;
                }
            }else{
                $score =~ s/\(.*\)//;
                if ( lc $f eq lc $score ) {    #damaging
                    return 1 if $opts{keep_any_damaging};
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
#    my $anno = shift; 
#    my %filter_matched = ();
    #TODO!
    ...
}

#################################################
sub outputGeneList{
    return if not $LIST;
    foreach my $k (sort keys %gene_listing){
        print $LIST "$k:\t" . join
        (
            ",", 
            @{ $gene_listing{$k} }
        ) . "\n";
    }
    close $LIST;
}
#################################################
sub is_autosome {
    my ($chrom) = @_;
    if ( $chrom =~ /^(chr)*[XYM]/i ) {
        return 0;
    }
    return 1;
}
#################################################
sub removeDups{
    my %seen = ();
    return grep { ! $seen{$_}++ } @_;
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

=head1 NAME

findBiallelic.pl - identify variants that make up potential biallelic variation of a gene

=head1 SYNOPSIS

    findBiallelic.pl -i <variants.vcf> -s <sample1> <sample2> [options]
    findBiallelic.pl --help (show help message)
    findBiallelic.pl --manual (show manual page)

=cut 

=head1 ARGUMENTS

=over 8 

=item B<-i    --input>

VCF file annotated with Ensembl's variant_effect_predictor.pl script or SnpEff.

=item B<-o    --output>

File to print output (optional). Will print to STDOUT by default.

=item B<-l    --list>

File to print a list of genes containing biallelic variants to. If you use this argument without specifying a value the list will be printed to STDERR;

=item B<-s    --samples>

One or more samples to identify biallelic genes from.  When more than one sample is given only genes with biallelic variants in ALL samples will be returned unless --num_matching option is used.

=item B<-r    --reject>

ID of samples to exclude variants from. Biallelic variants identified in these samples will be used to filter those found in samples supplied via the --samples argument.

=item B<-x    --reject_all_except>

Reject variants present in all samples except these. If used without an argument all samples in VCF that are not specified by --samples argument will be used to reject variants. If one or more samples are given as argument to this option then all samples in VCF that are not specified by --samples argument or this argument will be used to reject variants.

=item B<--x_linked>

By default this program only looks for potential biallelic variants on autosomes. This option takes a value between 0 and 2 to specify how to handle the X chromsome:

    1 -  look for biallelic/hemizygous variants on the X chromosome only 
    2 -  look for biallelic/hemizygous variants on the X chromosome in addition to autosomal biallelic variants. 
    0 -  look for biallelic variants on autosomes only (default)

=item B<-f    --family>

A PED file (format described below) containing information about samples in the VCF. In this way you can specify one or more families to allow the program to analyze biallelic variation that segregates correctly among affected and unaffected members. This assumes that any given family will have the same monogenic cause of recessive disease (i.e. this will not find phenocopies segregating in a family). One advantage of using a PED file is that the program can identify obligate carriers of a recessive disease and filter variants appropriately.  Can be used instead of or in conjunction with --samples (-s), --reject (-r) and --reject_all_except (-x) arguments. 

Not all samples in a given PED file need be present in the VCF. For example, you may specify an affected child not present in a VCF to indicate that an unaffected sample that IS present in the VCF is an obligate carrier. 

If your VCF contains more than one family you can specify one or more of these families by adding a comma after your ped file followed by a comma separated list of family IDs. For example:

    --family mypedfile.ped,Smith,Jones

...where 'mypedfile.ped' is the name of your ped file, while 'Smith' and 'Jones' are the family IDs of the families to analyze. To specify only family 'Smith' you would use the following argument: 

    --family mypedfile.ped,Smith

You may specify as many families as you like or none to analyze all families found in the PED file. Note, that it is assumed that you expect all families to have a mutation in the SAME GENE unless using the --num_matching argument. 

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

This script will ignore any lines in a PED file starting with '#' to allow users to include comments or headers.

=item B<--classes>

One or more mutation classes to retrieve. By default only variants labelled with one of the following classes will count towards biallelic variants:

        transcript_ablation
        splice_donor_variant
        splice_acceptor_variant
        stop_gained
        frameshift_variant
        stop_lost
        initiator_codon_variant
        inframe_insertion
        inframe_deletion
        missense_variant
        protein_altering_variant
        transcript_amplification
        TFBS_ablation
        TFBS_amplification
        regulatory_region_ablation
        regulatory_region_amplification

The user can specify one or more of the following classes instead: 

        transcript_ablation
        splice_donor_variant
        splice_acceptor_variant
        stop_gained
        frameshift_variant
        stop_lost
        initiator_codon_variant
        inframe_insertion
        inframe_deletion
        missense_variant
        protein_altering_variant
        transcript_amplification
        splice_region_variant
        incomplete_terminal_codon_variant
        synonymous_variant
        stop_retained_variant
        coding_sequence_variant
        mature_miRNA_variant
        5_prime_UTR_variant
        3_prime_UTR_variant
        intron_variant
        NMD_transcript_variant
        non_coding_exon_variant
        nc_transcript_variant
        upstream_gene_variant
        downstream_gene_variant
        TFBS_ablation
        TFBS_amplification
        TF_binding_site_variant
        regulatory_region_variant
        regulatory_region_ablation
        regulatory_region_amplification
        feature_elongation
        feature_truncation
        intergenic_variant


=item B<-a    --add_classes>

Specify one or more classes, separated by spaces, to add to the default mutation classes used for finding biallelic variants.

=item B<--consensus_splice_site>

Use this flag in order to keep splice_region_variant classes only if they are in a splice consensus region as defined by the SpliceConsensus plugin. You do not need to specify 'splice_region_variant' using --classes or --add_classes options when using this flag. You B<MUST> have used the SpliceConsensus plugin when running the VEP for this option to work correctly.

=item B<--canonical_only>

Only consider canonical transcripts.

=item B<-d    --damaging>

Specify SIFT, PolyPhen or Condel labels or scores to filter on. Add the names of the programs you want to use, separated by spaces, after the --damaging option. By default SIFT will keep variants labelled as 'deleterious', Polyphen will keep variants labelled as 'possibly_damaging' or 'probably_damaging' and  Condel will keep variants labelled as 'deleterious'.

If you want to filter on custom values specify values after each program name in the like so: 'polyphen=probably_damaging'. Seperate multiple values with commas - e.g. 'polyphen=probably_damaging,possibly_damaging,unknown'. You may specify scores between 0 and 1 to filter on scores rather than labels - e.g. 'sift=0.3'. For polyphen, variants with scores lower than this score are considered benign and filtered, for SIFT and Condel higher scores are considered benign.

Valid labels for SIFT: deleterious, tolerated

Valid labels for Polyphen: probably_damaging, possibly_damaging, benign, unknown

Valid labels for Condel : deleterious, neutral


To use default values for all three programs use 'all' (i.e. '--damaging all').

The default behaviour is to only keep variants predicted as damaging by ALL programs specified, although if the value is not available for one or more programs than that program will be ignored for filtering purposes.


=item B<-k    --keep_any_damaging>

If using multiple programs for filters for --damaging argument use this flag to keep variants predicted to be damaging according to ANY of these programs.

=item B<--skip_unpredicted_missense>

Skip variants that do not have a score from one or more programs specified by the --damaging argument. The --keep_any_damaging argument will override this behaviour if any of the available predictions are considered damaging.

=item B<-g    --gmaf>

Use a value between 0.00 and 0.50 to specify global minor allele frequencey filtering. If GMAF is available for variant it will be filtered if equal to or greater than the value specfied here.

=item B<--maf>

Like gmaf but filter on any population specific minor allele frequency annotated by the VEP as well as the GMAF.

=item B<-q    --quality>

Minimum genotype qualities to consider. This applies to samples specified by both --sample and --reject. Anything below this threshold will be considered a no call. Default is 20.

=item B<-w    --aff_quality>

Minimum genotype qualities to consider for affected samples only (i.e. samples specified by --sample argument or affected samples from a given pedigree). Anything below this threshold will be considered a no call. Default is 20.

=item B<-u    --un_quality>

Minimum genotype qualities to consider for unaffected samples only (i.e. samples specified by --reject argument or unaffected samples from a given pedigree). Anything below this threshold will be considered a no call. Default is 20.

=item B<-e    --equal_genotypes>

Use this flag if you only want to consider genotypes that are identical in each sample to count towards biallelic variation. Potentially useful if looking at several related individuals segregating the same disease and not using a PED file to specify their relationships (or if for some reason you think several families will have the same causative mutation).

=item B<--check_all_samples>

Check all samples in VCF. Assumes all samples are affected.

=item B<-n    --num_matching>

By default only transcripts that are found to contain potentially pathogenic variation in ALL affected samples are considered. Use this option to specify a minimum number of matching samples if you want to consider the possibility that one or more samples may not have a causative mutation in the same gene or to allow for no calls/calls falling below --quality threshold. If using a ped file, multiple samples from the same family will only count as a single match. For example, if you have two families and another sample specified using --samples you may use "--num_matching 2" to identify transcripts with matching variants in either both families or one family and the extra sample.  You may also use the --num_matching_per_family (-y) argument to allow for missing variants within families.

=item B<-y    --num_matching_per_family>

As above but for affected members per family. By default, even when using --num_matching, biallelic variants from a given family (if a PED file was used) are only passed on for consideration if all affected samples from that family present in the VCF contain matching biallelic variants. This argument sets the minimum number of affected samples per family for a biallelic combination of variants to be considered. 

=item B<-t    --ignore_carrier_status>

When using a PED file to specify family members, use this flag to prevent checking of carrier status of unaffected parents of affected children. Without this flag biallelic combinations will be rejected if an unaffected parent does not contain either variant (and the parent's genotype call quality is equal to or greater than the minimum genotype quality - see the --quality argument above). For example, you may want to use this option if you want to consider the possibility of a de novo occurence of one mutant allele or an undetected large deletion in one parent which may make an affected child's genotype look homozygous while the true genotype is compound heterozygosity for the mutation in the VCF alongside a large (undetected) deletion.

=item B<--pass_filters>

Only consider variants with a PASS filter field. If the FILTER field for variant is not PASS the variant will be skipped.

=item B<-b    --progress>

Show a progress bar while working.

=item B<-z    --homozygous_only>

Only consider homozygous variants, ignore potential compound heterozygotes (i.e. if autozygosity is assumed). 

=item B<-h    ---help>

Show the program's help message.

=item B<--manual>

Show the program's manual page.

=back

=cut

=head1 EXAMPLES

    findBiallelic.pl -i <variants.vcf> -s <sample1> <sample2> -r <sample3> <sample4>  -o output.vcf -l genelist.txt
    #find genes with biallelic variants in two unrelated samples but not in two unaffected samples. 

    findBiallelic.pl -i <variants.vcf> -s <sample1> <sample2> -r <sample3> <sample4> -d polyphen --maf 0.01 -o output.vcf -l genelist.txt
    #as above but only consider missense variants predicted damaging by polyphen and with a minor allele frequency less than 1%. 

    findBiallelic.pl -i <variants.vcf> -s <sample1> <sample2> -e -o output.vcf -l genelist.txt
    #find genes with biallelic variants in two related samples where you expect them to share the same causative variant.

    findBiallelic.pl -i <variants.vcf> -f families.ped -o output.vcf -l genelist.txt -q 30
    #use families.ped file to describe affected and unaffected individuals, only consider calls with genotype quality of 30 or higher
    
    findBiallelic.pl -i <variants.vcf> -f families.ped -o output.vcf -l genelist.txt -w 10 -z 30
    #use a low gneotype quality threshold (10) to identify variants in affected samples but use a higher threshold (30) to identify genotypes in unaffecteds to reject

    findBiallelic.pl -i <variants.vcf> -f families.ped,family1,family2 -o output.vcf 
    #specify specific families to analyze from ped file rather than all families

=cut

=head1 DESCRIPTION

This program reads VCF files annotated with Ensembl's Variant Effect Predictor or SnpEff and identifies transcripts with potential biallelic variation matching the various options specified above for the purpose of identifying potential recessively inherited pathogenic variants.  When more than one sample is specified using the --samples (-s) argument transcripts are identified that contain (not necessarily identical) potential biallelic variation in all samples. If multiple samples are specified in a PED file passed to the script with the --family (-f) argument, the program will attempt to find identical biallelic combinations within families and transcripts that contain potential biallelic variation in all affected samples from different families.

Genes are considered to contain potential biallelic variation if they either contain homozygous variants or two or more heterozygous variants. Phase can not be determined for variants so variants in cis may be erroneously considered to be potential biallelic variation.  Using variant data from unaffected parents with the --reject (-r) option or by specifying a PED file with the --family (-f) option  can help get around this issue.  Any samples specified using the --reject option will be used to remove biallelic variant combinations present - specifically, genotype combinations identified in affected samples (--samples) but also present in samples specified using the --reject argument will be removed from output. In this manner, if you have data from unaffected parents you should be able to avoid the problem of false positives from variants in cis as well as removing any shared genuine biallelic but non-pathogenic variation. However, this program does not require parental samples and can attempt to infer phase from a single parent if only one is available if the relationship is specified in a PED file passed to the script with the --family (-f) argument.

While related samples will be most useful in filtering using the --reject argument, data from any unaffected sample can be used to remove apparently non-pathogenic biallelic variation. Furthermore, while unrelated affected individuals can be used to identify shared genes containing apparent biallelic variation (when you believe the disorder to be caused by variation in the same gene), if using several related affected individuals you may use the --equal_genotypes flag to tell the program to only look for variants that are shared among all affected individuals AND potentially biallelic.


=cut

=head1 AUTHOR

David A. Parry

University of Leeds

=head1 COPYRIGHT AND LICENSE

Copyright 2012, 2013, 2014, 2015, David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut

