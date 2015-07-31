#!/usr/bin/perl
#David Parry January 2014
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Term::ProgressBar;
use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/lib";
use ParseVCF;
use SortGenomicCoordinates;
use ParsePedfile;
use VcfReader;

my $genotype_quality = 20;
my $aff_genotype_quality
  ;    #will convert to $genotype_quality value if not specified
my $unaff_genotype_quality
  ;    #will convert to $genotype_quality value if not specified
my $vcf;
my @samples;
my @reject;
my @reject_except;
my @head;
my @classes;
my $reject;
my $identical_genotypes;
my $allow_missing;
my @add;
my $list_genes = 0
  ; #will change to '' if user specifies --list but provides no argument, in which case print to STDERR
my $out;
my @damaging = ()
  ; #missense prediction programs (SIFT, polyphen, condel) to use for filtering missense variants
my $keep_any_damaging
  ; #default is that all scores from programs specified in @damaging must be 'damaging'
my $filter_unpredicted
  ; #ignore missense variants that don't have a score for programs used in @damaging -NOT IMPLEMENTED YET!
my $canonical_only;    #canonical transcripts only
my $pass_filters;      #only keep variants with PASS filter field.
my $gmaf;    #filter GMAF if present if equal to or above this value (0.0 - 0.5)
my $any_maf
  ;    #filter MAF if present if equal to or above this value in any population
my $help;
my $man;
my $progress;
my $check_all_samples;
my $homozygous_only;
my $splice_consensus =
  0;    #use this flag to check for SpliceConsensus VEP plugin annotations
my $pedigree;
my $min_matching_samples;
my $min_matching_per_family;
my $ignore_carrier_status;
my $x_linked = 0;
my %opts = (
    'input'                     => \$vcf,
    'output'                    => \$out,
    'list'                      => \$list_genes,
    'samples'                   => \@samples,
    'reject'                    => \@reject,
    'reject_all_except'         => \@reject_except,
    'family'                    => \$pedigree,
    'classes'                   => \@classes,
    'canonical_only'            => \$canonical_only,
    'pass_filters'              => \$pass_filters,
    'damaging'                  => \@damaging,
    'keep_any_damaging'         => \$keep_any_damaging,
    'skip_unpredicted_missense' => \$filter_unpredicted,
    'gmaf'                      => \$gmaf,
    'maf'                       => \$any_maf,
    'check_all_samples'         => \$check_all_samples,
    'equal_genotypes'           => \$identical_genotypes,
    'quality'                   => \$genotype_quality,
    'aff_quality'               => \$aff_genotype_quality,
    'un_quality'                => \$unaff_genotype_quality,
    'x_linked'                  => \$x_linked,
    #'allow_missing_genotypes' => \$allow_missing,
    'num_matching'             => \$min_matching_samples,
    'num_matching_per_family'  => \$min_matching_per_family,
    'ignore_carrier_status'    => \$ignore_carrier_status,
    'progress'                 => \$progress,
    'add_classes'              => \@add,
    'homozygous_only'          => \$homozygous_only,
    'consensus_splice_site'    => \$splice_consensus,
    'help'                     => \$help,
    'manual'                   => \$man
);
GetOptions(
    \%opts,
    'i|input=s' => \$vcf,
    'output=s',
    'list:s',
    's|samples=s{,}'           => \@samples,
    'r|reject=s{,}'            => \@reject,
    'x|reject_all_except:s{,}' => \@reject_except,
    'family=s',
    'classes=s{,}',
    'canonical_only',
    'pass_filters',
    'damaging=s{,}',
    'keep_any_damaging',
    'skip_unpredicted_missense',
    'gmaf=f',
    'maf=f',
    'x_linked=i', #1 = look for x-linked recessive only, 2= look for x-linked recessive as well
    'check_all_samples',
    'equal_genotypes',
    'quality=i',
    'w|aff_quality=i' => \$aff_genotype_quality,
    'u|un_quality=i'  => \$unaff_genotype_quality,

    #           'allow_missing_genotypes',
    'n|num_matching=i'            => \$min_matching_samples,
    'y|num_matching_per_family=i' => \$min_matching_per_family,
    't|ignore_carrier_status'     => \$ignore_carrier_status,
    'b|progress'                  => \$progress,
    'a|add_classes=s{,}'          => \@add,
    'z|homozygous_only'           => \$homozygous_only,
    'consensus_splice_site',
    'h|?|help' => \$help,
    'manual'
) or pod2usage( -message => "Syntax error", exitval => 2 );

pod2usage( -verbose => 2 ) if $man;
pod2usage( -verbose => 1 ) if $help;
pod2usage( -message => "Syntax error: input is required.", exitval => 2 )
  if not $vcf;
pod2usage(
    -message =>
"Syntax error: please specify samples to analyze using --samples (-s), --check_all_samples or --family (-f) arguments.",
    exitval => 2
  )
  if not @samples
  and not $check_all_samples
  and not $pedigree;
pod2usage(
    -message =>
"--gmaf option requires a value between 0.00 and 0.50 to filter on global minor allele frequency.\n",
    -exitval => 2
) if ( defined $gmaf && ( $gmaf < 0 or $gmaf > 0.5 ) );
pod2usage(
    -message => "Genotype quality scores must be 0 or greater.\n",
    -exitval => 2
) if ( $genotype_quality < 0 );

if ( defined $aff_genotype_quality ) {
    pod2usage(
        -message => "Genotype quality scores must be 0 or greater.\n",
        -exitval => 2
    ) if $aff_genotype_quality < 0;
}
else {
    $aff_genotype_quality = $genotype_quality;
}
if ( defined $unaff_genotype_quality ) {
    pod2usage(
        -message => "Genotype quality scores must be 0 or greater.\n",
        -exitval => 2
    ) if $unaff_genotype_quality < 0;
}
else {
    $unaff_genotype_quality = $genotype_quality;
}

#QUICKLY CHECK PEDIGREE BEFORE DEALING WITH POTENTIALLY HUGE VCF
my $ped;
if ($pedigree) {
    $ped = ParsePedfile->new( file => $pedigree );

#die "Pedigree file must contain at least one affected member!\n" if not $ped->getAllAffecteds();
}

#SORT OUT VEP DETAILS
my @valid = qw (transcript_ablation
  splice_donor_variant
  splice_acceptor_variant
  stop_gained
  frameshift_variant
  stop_lost
  initiator_codon_variant
  inframe_insertion
  inframe_deletion
  missense_variant
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
  intergenic_variant);

if ( not @classes ) {

#@classes = qw(missense nonsense stoploss deletion insertion splicing splice_consensus);
    @classes = qw (transcript_ablation
      splice_donor_variant
      splice_acceptor_variant
      stop_gained
      frameshift_variant
      stop_lost
      initiator_codon_variant
      inframe_insertion
      inframe_deletion
      missense_variant
      transcript_amplification
      TFBS_ablation
      TFBS_amplification
      regulatory_region_ablation
      regulatory_region_amplification);
}
push( @classes, @add ) if (@add);
push @classes, "splice_region_variant" if $splice_consensus;
foreach my $class (@classes) {
    die "Error - variant class '$class' not recognised.\n"
      if not grep { /$class/i } @valid;
}
my @csq_fields = qw (allele gene feature feature_type consequence hgnc)
  ;    #default fields to retrieve from CSQ INFO field
my %allelic_genes
  ; #record transcript id as key to anon hash with sample as key, and number of mutations as value
my %geneline;        #store lines here
my %id_to_symbol;    #hash of transcript id to gene symbol
my %listing = ()
  ; #hash of gene symbols with values being arrays of transcript IDs - persistent

my %damage_filters = ();  #hash of prediction program names and values to filter

if (@damaging) {
    my %valid_damaging = (
        sift => [ "deleterious", "tolerated" ],
        polyphen =>
          [ "probably_damaging", "possibly_damaging", "benign", "unknown" ],
        condel => [ "deleterious", "neutral" ]
    );
    my %default_damaging = (
        sift     => [ "deleterious", ],
        polyphen => [ "probably_damaging", "possibly_damaging", ],
        condel   => [ "deleterious", ]
    );
    foreach my $d (@damaging) {
        my ( $prog, $label ) = split( "=", $d );
        if ( exists $valid_damaging{ lc $prog } ) {
            push @csq_fields, lc($prog);
            no warnings 'uninitialized';
            my @filters =
              add_damaging_filters( lc $prog, lc $label, \%valid_damaging,
                \%default_damaging );
            push @{ $damage_filters{ lc $prog } }, @filters;
        }
        elsif ( lc($prog) eq 'all' ) {
            %damage_filters = %default_damaging;
            push @csq_fields, qw(sift polyphen condel);
        }
        else {
            die
"Unrecognised value ($d) passed to --damaging argument. See --help for more info.\n";
        }
    }
}
else {
    if ( $keep_any_damaging and $filter_unpredicted ) {
        die
"--keep_any_damaging and --unpredicted_missense arguments can only be used when --damaging argument is in effect.\n";
    }
    elsif ($keep_any_damaging) {
        die
"--keep_any_damaging argument can only be used when --damaging argument is in effect.\n";
    }
    elsif ($filter_unpredicted) {
        die
"--unpredicted_missense argument can only be used when --damaging argument is in effect.\n";
    }
}

if ($canonical_only) {
    push @csq_fields, 'canonical';
}

if ( defined $gmaf or defined $any_maf ) {
    push @csq_fields, 'gmaf';
}
if ($splice_consensus) {
    push @csq_fields, 'splice_consensus';
}

#PARSED VEP ARGUMENTS
#INITIALIZE VCF

print STDERR "Initializing VCF input ($vcf)...\n";
my $vcf_obj = ParseVCF->new( file => $vcf );

#CHECK VEP ARGUMENTS AGAINST VCF
my $vep_header     = $vcf_obj->readVepHeader();
my @available_mafs = ();
if ( defined $any_maf ) {
    foreach my $key ( keys %{$vep_header} ) {
        if ( $key =~ /\w_MAF$/ ) {
            push @available_mafs, $key;
            push @csq_fields,     $key;
        }
    }
}
if ( @csq_fields > 1 ) {
    my %seen = ();
    @csq_fields = grep { !$seen{$_}++ } @csq_fields;
}
my $replace_hgnc = 0;
foreach my $c (@csq_fields) {
    if ( not exists $vep_header->{$c} ) {
        if ( $c eq 'hgnc' ) {
            if ( not exists $vep_header->{'symbol'} )
            {    #they've awkwardly replaced hgnc with symbol in v73 and above
                die
"Couldn't find 'hgnc' or 'symbol' VEP field in header - please ensure your VCF is annotated with "
                  . "Ensembl's variant effect precictor specifying the appropriate annotations.\n";
            }
            else {
                $replace_hgnc++;
            }
        }
        else {
            die
"Couldn't find '$c' VEP field in header - please ensure your VCF is annotated with "
              . "Ensembl's variant effect precictor specifying the appropriate annotations.\n";
        }
    }
}

if ($replace_hgnc) {
    @csq_fields = grep { !/^hgnc$/ } @csq_fields;
    push @csq_fields, 'symbol';
}

#DONE CHECKING VEP INFO
#CHECK OUR SAMPLES/PEDIGREE INFO
if ($check_all_samples) {
    push @samples, $vcf_obj->getSampleNames();
}
if ($ped) {
    my @aff     = ();
    my @un      = ();
    my @not_aff = ();
    my @not_un  = ();
    foreach my $s ( $ped->getAllAffecteds() ) {
        if ( $vcf_obj->checkSampleInVcf($s) ) {
            push @aff, $s;
        }
        else {
            push @not_aff, $s;
        }
    }
    foreach my $s ( $ped->getAllUnaffecteds() ) {
        if ( $vcf_obj->checkSampleInVcf($s) ) {
            push @un, $s;
        }
        else {
            push @not_un, $s;
        }
    }
    print STDERR "Found "
      . scalar(@aff)
      . " affected samples from pedigree in VCF.\n";
    print STDERR scalar(@not_aff)
      . " affected samples from pedigree were not in VCF.\n";
    print STDERR "Found "
      . scalar(@un)
      . " unaffected samples from pedigree in VCF.\n";
    print STDERR scalar(@not_un)
      . " unaffected samples from pedigree were not in VCF.\n";
    if ($min_matching_per_family) {
        if ( not $min_matching_samples ) {
            print STDERR
"WARNING: --num_matching_per_family (-y) set without setting --num_matching (-n). "
              . "Setting --num_matching to $min_matching_per_family.\n";
            $min_matching_samples = $min_matching_per_family;
        }
        foreach my $f ( $ped->getAllFamilies ) {
            my %affected_ped =
              map { $_ => undef } $ped->getAffectedsFromFamily($f);
            my @affected = grep { exists $affected_ped{$_} } @aff;
            die "Number of affected found in family $f ("
              . scalar(@affected)
              . ") is less than value for "
              . "--num_matching_per_family ($min_matching_per_family)\n"
              if $min_matching_per_family > @affected;
        }
    }
    push @samples, @aff;
    push @reject,  @un;
}

#check @samples and @reject exist in file
my @not_found = ();
foreach my $s ( @samples, @reject ) {
    if ( not $vcf_obj->checkSampleInVcf($s) ) {
        push @not_found, $s;
    }
    if (@not_found) {
        die "Could not find the following samples in VCF:\n"
          . join( "\n", @not_found ) . "\n";
    }
}
die "No affected samples found in VCF\n" if not @samples;

if (@reject_except) {
    my @all = $vcf_obj->getSampleNames();
    push @reject_except, @samples;
    my %subtract = map { $_ => undef } @reject_except;
    @all = grep { !exists $subtract{$_} } @all;
    push @reject, @all;
}

#remove any duplicate samples in @samples or @reject
my %seen = ();
@reject  = grep { !$seen{$_}++ } @reject;
%seen    = ();
@samples = grep { !$seen{$_}++ } @samples;
%seen    = ();

#make sure no samples appear in both @samples and @reject
my %dup = map { $_ => undef } @samples;
foreach my $s (@samples) {
    my @doubles = grep { exists( $dup{$_} ) } @reject;
    die "Same sample(s) specified as both affected and unaffected:\n"
      . join( "\n", @doubles ) . "\n"
      if @doubles;
}

if ( defined $min_matching_samples ) {
    die "ERROR: --num_matching (-n) argument must be greater than 0.\n"
      if $min_matching_samples < 1;
    die
"ERROR: --num_matching (-n) value [$min_matching_samples] is greater than number of --samples identified ("
      . scalar(@samples) . ")\n"
      if $min_matching_samples > @samples;
    if ( $min_matching_samples < @samples ) {
        $allow_missing++;
    }
    if ($pedigree) {
        print STDERR "NOTE: --num_matching (-n) argument should be used in "
          . "conjunction with --num_matching_per_family (-y) when using a PED file "
          . "if you want to allow for missing variants within the same family.\n"
          if not $min_matching_per_family;
        my $aff_count =
          0;    #only count one sample per family for $min_matching_samples
        foreach my $f ( $ped->getAllFamilies() ) {
            $aff_count++ if ( $ped->getAffectedsFromFamily($f) );
        }
        foreach my $s (@samples) {
            $aff_count++ if not grep { $_ eq $s } $ped->getAllSamples();
        }
        die
"ERROR: --num_matching (-n) value [$min_matching_samples] is greater than the number of families "
          . "with affected members identified in $pedigree and --samples identified ($aff_count)\n"
          if $min_matching_samples > @samples;
    }
}

#DONE CHECKING SAMPLE INFO
#SET UP OUR OUTPUT FILES

my $OUT;
if ($out) {
    open( $OUT, ">$out" ) || die "Can't open $out for writing: $!\n";
}
else {
    $OUT = \*STDOUT;
}
my $LIST;
if ($list_genes) {
    open( $LIST, ">$list_genes" )
      || die "Can't open $list_genes for writing: $!\n";
}
elsif ( $list_genes eq '' )
{    #user specified --list option but provided no argument
    $LIST       = \*STDERR;
    $list_genes = 1;
}

print $OUT $vcf_obj->getHeader(0) . "##findBiallelicVep.pl\"";
my @opt_string = ();
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
print $OUT join( " ", @opt_string ) . "\"\n";
print $OUT
"##INFO=<ID=findBiallelicVepSamplesHom,Number=A,Type=String,Description=\"For each allele a list of samples that were found to meet findBiallelicVep.pl's criteria for constituting a homozygous variant.\">\n";
print $OUT
"##INFO=<ID=findBiallelicVepSamplesHet,Number=A,Type=String,Description=\"For each allele a list of samples that were found to meet findBiallelicVep.pl's criteria for constituting a potential compound heterozygous variant.\">\n";
print $OUT $vcf_obj->getHeader(1);

#ALLOW FOR CUSTOM VCF BY LOGGING CHROM AND POS HEADER COL #s FOR SORTING OF VCF LINES
my $chrom_header = $vcf_obj->getColumnNumber("CHROM");
my $pos_header   = $vcf_obj->getColumnNumber("POS");
my $prev_chrom   = 0;
my $line_count   = 0;

#SET UP PROGRESSBAR
my $progressbar;
if ($progress) {
    if ( $vcf eq "-" ) {
        print STDERR "Can't use --progress option when input is from STDIN\n";
        $progress = 0;
    }
    else {
        $progressbar = Term::ProgressBar->new(
            {
                name  => "Biallelic",
                count => $vcf_obj->countLines("variants"),
                ETA   => "linear",
            }
        );
    }
}
my $next_update = 0;

my %contigs = (); 
#START PROCESSING OUR VARIANTS
LINE: while ( my $line = $vcf_obj->readLine ) {
    $line_count++;
    if ($progress) {
        $next_update = $progressbar->update($line_count)
          if $line_count >= $next_update;
    }
    my $chrom = $vcf_obj->getVariantField("CHROM");
    if (exists $contigs{$chrom}){
        if ($contigs{$chrom} != scalar(keys %contigs) - 1){
            die
"Encountered a variant for contig $chrom with at least one variant for a different contig inbetween. Your VCF input must be sorted such that all contigs are kept together - please sort your input and try again.\n";
        }
    }else{
        $contigs{$chrom} = scalar(keys%contigs);
    }
    if ( $prev_chrom && $chrom ne $prev_chrom ) {

        #print biallelic mutations for previous chromosome here...
        my %vcf_lines = ();
        foreach my $gene ( keys %allelic_genes ) {
            my $found = 0;

            #CHECK SEGREGATION IF PEDIGREE
            if ($pedigree) {
                 $found =  check_segregation( \%vcf_lines,  $ped, \%{ $allelic_genes{$gene} }, $chrom);
            }
            else {
                $found =
                  check_all_samples_biallelic( \%vcf_lines, \%{ $allelic_genes{$gene} });
            }
            if ($found) {
                push( @{ $listing{ $id_to_symbol{$gene} } }, $gene )
                  if $list_genes;
            }
        }
        my @lines = ();
        foreach my $k (keys %vcf_lines){#keys are vcf lines
            my @split = split("\t", $k); 
            my $converted_line = $k;
            foreach my $l (keys %{$vcf_lines{$k}}){#keys are new INFO fields
                my $new_split = VcfReader::addVariantInfoField(
                    line  => \@split,
                    id    => $l,
                    value => $vcf_lines{$k}->{$l},
                );
                $converted_line = join("\t", @$new_split);
            }
            push @lines, $converted_line;
        }
        my $sort = sort_vcf_lines( \@lines, $chrom_header, $pos_header );
        print $OUT join( "\n", @$sort ) . "\n" if @$sort;
        %allelic_genes = ();
        %geneline      = ();
        %id_to_symbol  = ();
    }

    $prev_chrom = $chrom;
    if ($pass_filters) {
        next if $vcf_obj->getVariantField("FILTER") ne 'PASS';
    }
    if ($x_linked == 0){
        next if not is_autosome($chrom);
    }elsif ($x_linked == 1){
        next if $chrom !~ /^(chr)*X$/i;
    }elsif ($x_linked == 2){
        if (not is_autosome($chrom)){
            next if $chrom !~ /^(chr)*X$/i;
        }
    }
    my $have_variant = 0;
    foreach my $sample (@samples) {
        if (
            $vcf_obj->getSampleCall(
                sample => $sample,
                minGQ  => $aff_genotype_quality
            ) =~ /(\d+)[\/\|](\d+)/
          )
        {
            $have_variant++ if ( $1 > 0 or $2 > 0 );
        }
        last if $have_variant;
    }
    next LINE if not $have_variant;

    if ($identical_genotypes) {
        unless ( defined $min_matching_samples
            && $min_matching_samples < @samples )
        {
            next LINE if not identical_genotypes(@samples);
        }
    }

    #check for identical genotypes within family if using a ped file
    if (    $pedigree
        and not $min_matching_per_family
        and not $min_matching_samples )
    {
        foreach my $fam ( $ped->getAllFamilies() ) {
            next LINE
              if not identical_genotypes( $ped->getAffectedsFromFamily($fam) );
        }
    }

    my @csq = $vcf_obj->getVepFields( \@csq_fields )
      ; #returns array of hashes e.g. $csq[0]->{Gene} = 'ENSG000012345' ; $csq[0]->{Consequence} = 'missense_variant'
    die
"No consequence field found for line:\n$line\nPlease annotated your VCF file with ensembl's variant effect precictor before running this program.\n"
      if not @csq;
  CSQ: foreach my $annot (@csq) {
        my @anno_csq = split( /\&/, $annot->{consequence} );

        #skip NMD transcripts
        if ( grep { /NMD_transcript_variant/i } @anno_csq ) {
            next CSQ;
        }
        my $matches_class = 0;
        next if ( $annot->{consequence} eq 'intergenic_variant' );
        if ($canonical_only) {
            next if ( not $annot->{canonical} );
        }
        if ( defined $gmaf ) {
            if ( $annot->{gmaf} ) {
                if ( $annot->{gmaf} =~ /\w+:(\d\.\d+)/ ) {
                    next if $1 >= $gmaf;
                }
            }
        }
        if ( defined $any_maf ) {
            foreach my $some_maf (@available_mafs) {
                if ( $annot->{$some_maf} ) {
                    if ( $annot->{$some_maf} =~ /\w+:(\d\.\d+)/ ) {
                        next if $1 >= $any_maf;
                    }
                }
            }
        }

      CLASS: foreach my $class (@classes) {
          ANNO: foreach my $ac (@anno_csq) {
                if ( lc $ac eq lc $class ) {
                    if ( lc $class eq 'missense_variant' and %damage_filters ) {
                        next ANNO
                          if (
                            filter_missense(
                                $annot,             \%damage_filters,
                                $keep_any_damaging, $filter_unpredicted
                            )
                          );
                    }
                    elsif ( lc $class eq 'splice_region_variant'
                        and $splice_consensus )
                    {
                        my $consensus = $annot->{splice_consensus};
                        next if not $consensus;
                        if ( $consensus !~ /SPLICE_CONSENSUS\S+/i ) {
                            print STDERR
"WARNING - SPLICE_CONSENSUS annotation $consensus is"
                              . " not recognised as an annotation from the SpliceConsensus VEP plugin.\n";
                            next;
                        }
                    }
                    $matches_class++;
                    last CLASS;
                }
            }
        }
        if ($matches_class) {
            my %var_hash =
              create_var_hash( $annot, $vcf_obj, \@samples, \@reject );
            foreach my $k ( keys %var_hash ) {

                #$allelic_genes{$subannot[2]}->{$k} =  $var_hash{$k};
                $allelic_genes{ $annot->{feature} }->{$k} = $var_hash{$k};
            }

            # creates a  structure like:
            # $hash{transcript}->{chr:pos/allele}->{sample} = count
            # and $hash{transcript}->{chr:pos/allele}->{mutation} = $annotation
            # and $hash{transcript}->{chr:pos/allele}->{vcf_line} = $line
            # containing info for all relevant @classes
            if ( $annot->{symbol} ) {
                $id_to_symbol{ $annot->{feature} } = $annot->{symbol};
            }
            elsif ( $annot->{hgnc} ) {
                $id_to_symbol{ $annot->{feature} } = $annot->{hgnc};
            }
            elsif ( $annot->{gene} ) {
                $id_to_symbol{ $annot->{feature} } = $annot->{gene};
            }
            else {
                $id_to_symbol{ $annot->{feature} } = $annot->{feature};
            }
        }
    }
}

my %vcf_lines = ();
foreach my $gene ( keys %allelic_genes ) {
    my $found = 0;

    #CHECK SEGREGATION IF PEDIGREE
    if ($pedigree) {
        $found = check_segregation( \%vcf_lines, $ped, \%{ $allelic_genes{$gene} }, $prev_chrom) ;
    }
    else {
        $found = check_all_samples_biallelic(\%vcf_lines,  \%{ $allelic_genes{$gene} }) ;
    }
    if ($found) {
        push( @{ $listing{ $id_to_symbol{$gene} } }, $gene ) if $list_genes;
    }
}
my @lines = (); 
foreach my $k (keys %vcf_lines){#keys are vcf lines
    my @split = split("\t", $k); 
    my $converted_line = $k;
    foreach my $l (keys %{$vcf_lines{$k}}){#keys are new INFO fields
        my $new_split = VcfReader::addVariantInfoField(
            line  => \@split,
            id    => $l,
            value => $vcf_lines{$k}->{$l},
        );
        $converted_line = join("\t", @$new_split);
    }
    push @lines, $converted_line;
}
my $sort = sort_vcf_lines( \@lines, $chrom_header, $pos_header );
print $OUT join( "\n", @$sort ) . "\n" if @$sort;

if ($progressbar) {
    $progressbar->update( $vcf_obj->countLines("variants") )
      if $vcf_obj->countLines("variants") >= $next_update;
}

if ($list_genes) {
    my $list = sort_gene_listing( \%listing );
    print $LIST join( "\n", @$list ) . "\n";
}

###########
sub is_autosome {
    my ($chrom) = @_;
    if ( $chrom =~ /^(chr)*[XYM]/i ) {
        return 0;
    }
    return 1;

    #if ($chrom =~ /^(chr)*(\d+|GL\d+|Un_)/i){
    #    return 1;
    #}
    #return 0;
}

###########
sub identical_genotypes {
    my (@samp) = @_;
    my %gts = $vcf_obj->getSampleCall(
        multiple => \@samp,
        minGQ    => $aff_genotype_quality
    );
    my %no_calls;
    for ( my $i = 0 ; $i < $#samp ; $i++ ) {
        if ( $allow_missing and $gts{ $samp[$i] } =~ /^\.[\/\|]\.$/ )
        {    #if we're allowing missing values then skip no calls
            $no_calls{ $samp[$i] }++;
            next;
        }
        elsif ( $gts{ $samp[$i] } =~ /^\.[\/\|]\.$/ )
        {    #otherwise a no call means we should go to the next line
            return 0;
        }
        for ( my $j = $i + 1 ; $j <= $#samp ; $j++ ) {
            if ( $allow_missing and $gts{ $samp[$j] } =~ /^\.[\/\|]\.$/ ) {
                $no_calls{ $samp[$j] }++;
                next;
            }
            elsif ( $gts{ $samp[$i] } =~ /^\.[\/\|]\.$/ ) {
                return 0;
            }
            return 0 if $gts{ $samp[$j] } ne $gts{ $samp[$i] };
        }
    }

    #so, if we're here all @samp are identical (or no calls if $allow_missing)
    return 0
      if keys %no_calls == @samp
      ; #even if we $allow_missing we don't want to print a variant if none of our @samphave a call
    return 1;
}

###########
sub filter_missense {

#returns 1 if missense should be filtered
#otherwise returns 0;
# uses $keep_any setting to return 0 if any of these predictions match, otherwise all available
# scores must be deleterious/damaging
# if $filter_missing is used a variant will be filtered if no score is available (overriden by $keep_any setting)
    my ( $anno, $filter_hash, $keep_any, $filter_missing ) = @_;

#my %default_damaging = (sift => ["deleterious", ],  polyphen => ["probably_damaging", "possibly_damaging",], condel => ["deleterious", ]);
    my %filter_matched = ();
  PROG: foreach my $k ( sort keys %$filter_hash ) {
        my $score = $anno->{ lc $k };
        if ( not $score or $score =~ /^unknown/i )
        { #don't filter if score is not available for this variant unless $filter_missing is in effect
            $filter_matched{$k}++ unless $filter_missing;
            next;
        }
      SCORE: foreach my $f ( @{ $filter_hash->{$k} } ) {
            if ( $f =~ /^\d(\.\d+)*$/ ) {
                my $prob;
                if ( $score =~ /^(\d(\.\d+)*)/ ) {
                    $prob = $1;
                }
                else {
                    next SCORE
                      ; #if score not available for this feature ignore and move on
                }
                if ( lc $k eq 'polyphen' ) {
                    if ( $prob >= $f )
                    {    #higher is more damaging for polyphen - damaging
                        return 0 if $keep_any;
                        $filter_matched{$k}++;
                        next PROG;
                    }
                    else {    #benign
                    }
                }
                else {
                    if ( $prob <= $f )
                    {    #lower is more damaging for sift and condel - damaging
                        return 0 if $keep_any;
                        $filter_matched{$k}++;
                        next PROG;
                    }
                    else {    #benign
                    }
                }
            }
            else {
                $score =~ s/\(.*\)//;
                if ( lc $f eq lc $score ) {    #damaging
                    return 0 if $keep_any;
                    $filter_matched{$k}++;
                    next PROG;
                }
            }
        }

    }
    foreach my $k ( sort keys %$filter_hash ) {

 #filter if any of sift/condel/polyphen haven't matched our deleterious settings
        return 1 if not exists $filter_matched{$k};
    }
    return 0;
}
###########
sub add_damaging_filters {
    my ( $prog, $label, $valid_hash, $default_hash ) = @_;
    if ($label) {
        if ( $label =~ /^\d(\.\d+)*$/ ) {
            die
"Numeric values for condel, polyphen or sift filtering must be between 0 and 1.\n"
              if ( $label < 0 or $label > 1 );
            return $label;
        }
        else {
            my @lb = split( ",", $label );
            foreach my $l (@lb) {
                die
"Invalid filter parameter '$l' used with --damaging argument.  See --help for valid arguments.\n"
                  if not grep { /^$l$/ } @{ $valid_hash->{$prog} };
            }
            return @lb;
        }

    }
    else {    #give default values for $prog
        return @{ $default_hash->{$prog} };
    }
}

###########
sub get_alleles_to_reject {

  #requires array ref to sample IDs and hash ref to gene count hash
  #returns two hash refs - see %reject_genotypes and %incompatible_alleles below
    my ( $rej, $gene_counts ) = @_;
    my @keys             = ( keys %{$gene_counts} );  #keys are "chr:pos-allele"
    my %reject_genotypes = ()
      ; #collect all genotype combinations present in @$reject samples - these cannot be pathogenic
    my %incompatible_alleles = ()
      ; #key is allele, value is an array of alleles each of which can't be pathogenic if key allele is pathogenic and vice versa
    for ( my $i = 0 ; $i < @keys ; $i++ ) {
        foreach my $r (@$rej) {
            if ( $gene_counts->{ $keys[$i] }->{$r} >= 1 ) {
                $reject_genotypes{"$keys[$i]/$keys[$i]"}++
                  if $gene_counts->{ $keys[$i] }->{$r} >=
                  2;    #homozygous therefore biallelic
                for ( my $j = $i + 1 ; $j < @keys ; $j++ )
                { #check other alleles to see if there are any compund het combinations
                    if ( $gene_counts->{ $keys[$j] }->{$r} >= 1 ) {
                        $reject_genotypes{"$keys[$i]/$keys[$j]"}++;
                        push @{ $incompatible_alleles{ $keys[$i] } }, $keys[$j];
                        push @{ $incompatible_alleles{ $keys[$j] } }, $keys[$i];
                    }
                }
            }
        }
    }
    foreach my $k ( keys %incompatible_alleles ) {
        my %seen = ();
        @{ $incompatible_alleles{$k} } =
          grep { !$seen{$_}++ } @{ $incompatible_alleles{$k} };
    }
    return ( \%reject_genotypes, \%incompatible_alleles );
}

###########
sub get_biallelic {

#arguments are array ref to samples, hash ref to reject genotypes, hash ref to incompatible genotypes and gene_counts hash ref
#if $count_per_family is true only count multiple members of family once
#returns hash of samples to arrays of potential biallelic genotypes
    my ( $aff, $reject_geno, $incompatible, $min_matches, $gene_counts,
        $count_per_family )
      = @_;
    $min_matches = @$aff if ( not $min_matches );
    my @keys = ( keys %{$gene_counts} );    #keys are "chr:pos-allele"
    my %possible_biallelic_genotypes =
      ();   #keys are samples, values are arrays of possible biallelic genotypes
    my %genotypes = ()
      ; #keys are all possible biallelic genotypes present in %possible_biallelic_genotypes
    my %biallelic = ()
      ; #same as above, for returning final list of valid biallelic combinations
    for ( my $i = 0 ; $i < @keys ; $i++ ) {
        next
          if $reject_geno->{"$keys[$i]/$keys[$i]"}
          ;    #allele $i can't be pathogenic if homozygous in a @reject sample
        foreach my $s (@$aff) {
            if ( $gene_counts->{ $keys[$i] }->{$s} >= 1 ) {
                if ( $gene_counts->{ $keys[$i] }->{$s} >= 2 )
                {    #homozygous therefore biallelic
                    push @{ $possible_biallelic_genotypes{$s} },
                      "$keys[$i]/$keys[$i]";
                    $genotypes{"$keys[$i]/$keys[$i]"}++;
                }
                if ( not $homozygous_only )
                {    #don't consider hets if --homozygous_only flag is in effect
                    for ( my $j = $i + 1 ; $j < @keys ; $j++ )
                    { #check other alleles to see if there are any compund het combinations
                        if ( $gene_counts->{ $keys[$j] }->{$s} >= 1 ) {
                            push @{ $possible_biallelic_genotypes{$s} },
                              "$keys[$i]/$keys[$j]"
                              if not $reject_geno->{"$keys[$i]/$keys[$j]"};
                            $genotypes{"$keys[$i]/$keys[$j]"}++;
                        }
                    }
                }
            }
        }
    }

#so now our $aff only have genotypes not present in %{$reject_geno}
#however, between all our samples we could be using pairs of alleles from %{$incompatible}
# i.e. if two alleles are present in an unaffected (@reject) sample at most one allele can be pathogenic
# so we can't use both in different affected (@$aff) samples

    foreach my $gt ( keys %genotypes ) {
        my %incomp_gt =
          ();    # keep the alleles incompatible with current $gt here
        my %compatible_genotypes = ()
          ; #samples are keys, values are genotypes that pass our test against incompatible alleles
        foreach my $allele ( split( /\//, $gt ) ) {
            if ( exists $incompatible->{$allele} ) {
                %incomp_gt = map { $_ => undef } @{ $incompatible->{$allele} };
            }
        }
        foreach my $s (@$aff) {
            foreach my $s_gt ( @{ $possible_biallelic_genotypes{$s} } ) {
                if ($identical_genotypes) {
                    next if $s_gt ne $gt;
                }
                my @s_alleles = ( split( "\/", $s_gt ) );
                if (    not exists $incomp_gt{ $s_alleles[0] }
                    and not exists $incomp_gt{ $s_alleles[1] } )
                {
                    push @{ $compatible_genotypes{$s} }, $s_gt;
                }
            }
        }

     #Test no. samples with compatible genotypes
     #against min. no. required matching samples or if not specified all samples
        if ( keys %compatible_genotypes >= $min_matches ) {
            my $match_fam = 0;
            if ($count_per_family) {
                my %counted_fam = ();
                foreach my $s ( keys %compatible_genotypes ) {
                    if ( grep { $s eq $_ } $ped->getAllSamples() ) {
                        my $fam_id = $ped->getFamilyId($s);
                        $match_fam++ if not $counted_fam{$fam_id};
                        $counted_fam{$fam_id}++;
                    }
                    else {
                        $match_fam++;
                    }
                }
            }
            if ( not $count_per_family or $match_fam >= $min_matches ) {
                foreach my $s ( keys %compatible_genotypes ) {
                    foreach my $sgt ( @{ $compatible_genotypes{$s} } ) {
                        push @{ $biallelic{$s} }, $sgt;
                    }
                }
            }
        }
    }

    foreach my $k ( keys %biallelic ) {

        #remove duplicate genotypes for each sample
        my %seen = ();
        @{ $biallelic{$k} } = grep { !$seen{$_}++ } @{ $biallelic{$k} };
    }
    return %biallelic;
}

###########
sub check_segregation {

    #checks whether variants segregate appropriately per family
    #can check between multiple families using check_all_samples_biallelic
    #returns hash from $gene_counts with only correctly segregating alleles
    my ( $vcf_lines, $p, $gene_counts, $chrom ) = @_;
    my @keys      = ( keys %{$gene_counts} );    #keys are "chr:pos-allele"
    my %seg_count = ()
      ; #return version of %{$gene_counts} containing only alleles that appear to segregate correctly

  #Iterate over all @reject to store %reject_genotypes and %incompatible_alleles
    my ( $reject_genotypes, $incompatible_alleles ) =
      get_alleles_to_reject( \@reject, $gene_counts );

    #Iterate overl all @samples to store possible biallelic genotypes
    my %biallelic_candidates =
      get_biallelic( \@samples, $reject_genotypes, $incompatible_alleles,
        $min_matching_samples, $gene_counts );

#Then for each family identify common biallelic genotypes but throw away anything if neither allele present in an obligate carrier

    foreach my $f ( $ped->getAllFamilies ) {
        my %affected_ped = map { $_ => undef } $ped->getAffectedsFromFamily($f);
        my @affected = grep { exists $affected_ped{$_} } @samples;

#get intersection of all @affected biallelic genotypes (we're assuming there's no phenocopy so we're looking for identical genotypes)
        my %intersect = ();
        foreach my $s (@affected) {
            foreach my $bi ( @{ $biallelic_candidates{$s} } ) {
                $intersect{$bi}++;
            }
        }
        if ($min_matching_per_family) {
            foreach my $k ( keys %intersect ) {
                delete $intersect{$k}
                  if $intersect{$k} < $min_matching_per_family;
            }
        }
        else {    #by default require all affecteds to match
            foreach my $k ( keys %intersect ) {
                delete $intersect{$k} if $intersect{$k} < @affected;
            }
        }

#we've already checked that these allele combinations are not present in unaffected members including parents.
        if ( not $ignore_carrier_status and $chrom !~ /(chr)*X/i) {

#check that parents (if available) contain one allele of a compound het (admittedly this does not allow for hemizygous variants in case of a deletion in one allele)
            my %unaffected =
              map { $_ => undef } $ped->getUnaffectedsFromFamily($f);
          KEY: foreach my $key ( keys %intersect ) {
                my @al = split( /\//, $key );
                foreach my $aff ( keys %affected_ped )
                {    #check parents of any affected even if affected not in VCF
                     #IF MIN MATCHING PER FAMILY IS IN USE MAKE SURE WE ONLY CHECK AGAINST SAMPLES WITH ALLELES
                    if ($min_matching_per_family) {
                        next
                          if not grep { $key eq $_ }
                          @{ $biallelic_candidates{$aff} };
                    }
                    my @par =
                      grep { exists $unaffected{$_} } $ped->getParents($aff);
                    foreach my $u (@par) {

                        #check each parent carries at least one of the alleles
                        if (    $gene_counts->{ $al[0] }->{$u} == 0
                            and $gene_counts->{ $al[1] }->{$u} == 0 )
                        { #0 means called as not having allele, -1 means no call
                            delete $intersect{$key};
                            next KEY;
                        }
                    }
                    if ( @par == 2 )
                    { #if we have both parents make sure there's one of each allele per parent if compound het
                        if (   $gene_counts->{ $al[0] }->{ $par[0] } == 0
                            && $gene_counts->{ $al[0] }->{ $par[1] } == 0 )
                        {
                            delete $intersect{$key};
                            next KEY;
                        }
                        elsif ($gene_counts->{ $al[1] }->{ $par[0] } == 0
                            && $gene_counts->{ $al[1] }->{ $par[1] } == 0 )
                        {
                            delete $intersect{$key};
                            next KEY;
                        }
                    }
                }
            }
        }

        #now %intersect only has viable genotypes for this family
        #Put these genotypes into %seg_counts
        foreach my $k ( keys %intersect ) {
            my @al = split( /\//, $k );

            #add viable alleles to %seg_count
            $seg_count{ $al[0] } = $gene_counts->{ $al[0] };
            $seg_count{ $al[1] } = $gene_counts->{ $al[1] };
        }
    }

    #deal with any samples not in pedigree
    my %all_ped = map { $_ => undef } $ped->getAllSamples();
    my @other_affected = grep { !exists $all_ped{$_} } @samples;
    foreach my $aff (@other_affected) {
        foreach my $bi ( @{ $biallelic_candidates{$aff} } ) {
            my @al = split( /\//, $bi );

            #add viable alleles to %seg_count
            $seg_count{ $al[0] } = $gene_counts->{ $al[0] };
            $seg_count{ $al[1] } = $gene_counts->{ $al[1] };
        }
    }

#Find any common variation by running check_all_samples_biallelic using \%seg_counts
    if ($min_matching_samples) {
        return check_all_samples_biallelic( $vcf_lines, \%seg_count, 1 );
    }
    else {
        return check_all_samples_biallelic( $vcf_lines, \%seg_count );
    }
}

###########
sub check_all_samples_biallelic {
    my ( $vcf_lines, $gene_counts, $count_per_family ) = @_;
    my $found = 0; #return 1 if we've found a biallelic combination for this gene
#we need to go through every possible combination of biallelic alleles
#(represented as chr:pos-allele) to compare between @samples and against @$reject
    my @keys = ( keys %{$gene_counts} )
      ; #keys are "chr:pos-allele" values are hashes of samples to allele counts
    my %possible_biallelic_genotypes =
      ();   #keys are samples, values are arrays of possible biallelic genotypes
     #first check @$reject alleles and collect non-pathogenic genotpes in %reject_genotypes - the assumption here is that
     # the disease alleles are rare so not likely to be in cis in a @reject sample while in trans in an affected sample therefore we can reject
     # assuming presence of two alleles in a reject sample means either they are in cis or are harmless when in trans
     #also note alleles that can't BOTH be pathogenic storing them in %incompatible alleles
    my ( $reject_genotypes, $incompatible_alleles ) =
      get_alleles_to_reject( \@reject, $gene_counts );
    my %biallelic_candidates =
      get_biallelic( \@samples, $reject_genotypes, $incompatible_alleles,
        $min_matching_samples, $gene_counts, $count_per_family );
    $found = keys %biallelic_candidates;
#%biallelic_candidates - keys are samples, values are arrays of biallelic genotypes
    foreach my $s ( keys %biallelic_candidates ) {
        foreach my $gt ( @{ $biallelic_candidates{$s} } ) {
            foreach my $allele ( split( /\//, $gt ) ) {
                my @split_line =
                  split( "\t", $gene_counts->{$allele}->{vcf_line} );
                my $field = "findBiallelicVepSamplesHom";
                if ( $gene_counts->{$allele}->{$s} == 1 ) {
                    $field = "findBiallelicVepSamplesHet";
                }
                my %allele_to_sample = ();
                if ( my $binfo = $vcf_lines->{$gene_counts->{$allele}->{vcf_line}}->{$field} )
                {
                    my @perAllele = split( /\,/, $binfo );
                    for ( my $i = 0 ; $i < @perAllele ; $i++ ) {
                        push @{ $allele_to_sample{ $i + 1 } },
                          grep { $_ ne '.' } split( /\|/, $perAllele[$i] );
                    }
                }
                my $c = ( split /-/, $allele )[1];
                if ( exists $allele_to_sample{$c} ) {
                    if ( not grep { $_ eq $s } @{ $allele_to_sample{$c} } ) {
                        push @{ $allele_to_sample{$c} }, $s;
                    }
                }
                else {
                    push @{ $allele_to_sample{$c} }, $s;
                }

                #write INFO field like so:
                #findBiallelicSamplesHom=sample1|sample2,sample5|sample6|sample7
                #findBiallelicSamplesHet=sample3,sample4
                @{ $allele_to_sample{$c} } = sort @{ $allele_to_sample{$c} };
                my @alts = VcfReader::readAlleles(
                    line        => \@split_line,
                    alt_alleles => 1
                );
                my @info_add = ();
                for ( my $i = 0 ; $i < @alts ; $i++ ) {
                    if ( not exists $allele_to_sample{ $i + 1 } ) {
                        push @info_add, '.';
                    }
                    else {
                        if ( not @{ $allele_to_sample{ $i + 1 } } ) {
                            push @info_add, '.';
                        }
                        else {
                            push @info_add,
                              join( "|", @{ $allele_to_sample{ $i + 1 } } );
                        }
                    }
                }
                $vcf_lines->{$gene_counts->{$allele}->{vcf_line}}->{$field} =
                    join( ",", @info_add );
            }
        }
    }
    return $found;
}
###########
sub check_keys_are_true {
    my ( $key_array, $hash ) = @_;
    foreach my $k (@$key_array) {
        return 0 if not $hash->{$k};
    }
    return 1;
}

###########
sub create_var_hash {
    my ( $annotation, $vcf_obj, $aff, $un ) = @_;
    my %var_hash;
    my $coord = $vcf_obj->getVariantField("CHROM") . ":"
      . $vcf_obj->getVariantField("POS");

    #we should check sample alleles against subannot alleles
    my @alts = $vcf_obj->readAlleles( alt_alleles => 1 );
    my $ref = $vcf_obj->getVariantField("REF");

    #(Allele Gene Feature Feature_type Consequence HGNC);
    my $i = 0;   #count alleles as 1-based list to correspond to GT field in VCF
    my @v_all = $vcf_obj->altsToVepAllele( ref => $ref, alt => \@alts );
    my %vep_alleles = ();
    @vep_alleles{@alts} = @v_all; 
    foreach my $alt (@alts) {
        $i++;
        if ( uc($vep_alleles{$alt}) eq uc( $annotation->{allele} ) ) {
            $var_hash{"$coord-$i"}->{mutation} = $annotation;
            $var_hash{"$coord-$i"}->{vcf_line} = $vcf_obj->get_currentLine;
            foreach my $s (@$aff){
                addSampleToVarHash
                (
                    $s, 
                    $coord, 
                    $i,
                    $aff_genotype_quality,
                    \%var_hash,
                    $vcf_obj,
                );
            }
            foreach my $s (@$un) {
                addSampleToVarHash
                (
                    $s, 
                    $coord, 
                    $i,
                    $unaff_genotype_quality,
                    \%var_hash,
                    $vcf_obj,
                );
            }
        }
    }
    return %var_hash;
}

###########
sub addSampleToVarHash{
    my ($sample, $coord, $allele, $gq, $hash, $vcf_obj) = @_;
    my $gt = $vcf_obj->getSampleCall(
        sample => $sample,
        minGQ  => $gq
    );
    if ( $gt =~ /^$allele[\/\|]$allele$/ ) {
        $$hash{"$coord-$allele"}->{$sample} =
            2;    #homozygous for this alt allele
    }
    elsif ( $gt =~ /^[\d+\.][\/\|]$allele$/ or $gt =~ /^$allele[\/\|][\d+\.]$/ ) {
        $$hash{"$coord-$allele"}->{$sample} = 1;    #het for alt allele
    }
    elsif ( $gt =~ /^\.[\/\|]\.$/ ) {
        $$hash{"$coord-$allele"}->{$sample} = -1;    #no call
    }
    else {
        $$hash{"$coord-$allele"}->{$sample} = 0;    #does not carry alt allele
    }
}

###########
sub sort_gene_listing {
    my ($gene_list) = @_;

#$gene_list is a ref to hash with keys =  GeneSymbol and values = array of transcript IDs
    my @sorted_list = ();
    foreach my $k ( sort keys %$gene_list ) {
        push @sorted_list, join( ":", $k, sort @{ $gene_list->{$k} } );
    }
    return \@sorted_list;
}
###########
sub sort_vcf_lines {
    my ( $v_lines, $chrom_col, $pos_col ) = @_;

    #remove duplicates
    my %seen = ();
    @$v_lines = grep { !$seen{$_}++ } @$v_lines;

    #sort in coordinate order
    my $sort_obj = SortGenomicCoordinates->new(
        array     => $v_lines,
        type      => "custom",
        col       => $chrom_col,
        start_col => $pos_col,
        stop_col  => $pos_col
    );
    $sort_obj->order();
    return $sort_obj->get_ordered;
}
###########
#=item B<--allow_missing>
#When multiple --samples are being analysed use this flag to stop the script rejecting variants that are missing (as in no genotype call) from other samples.

=head1 NAME

findBiallelicVep.pl - identify variants that make up potential biallelic variation of a gene

=head1 SYNOPSIS

    findBiallelicVep.pl -i <variants.vcf> -s <sample1> <sample2> [options]
    findBiallelicVep.pl --help (show help message)
    findBiallelicVep.pl --manual (show manual page)

=cut 

=head1 ARGUMENTS

=over 8 

=item B<-i    --input>

VCF file annotated with Ensembl's variant_effect_predictor.pl script.

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

    findBiallelicVep.pl -i <variants.vcf> -s <sample1> <sample2> -r <sample3> <sample4>  -o output.vcf -l genelist.txt
    #find genes with biallelic variants in two unrelated samples but not in two unaffected samples. 

    findBiallelicVep.pl -i <variants.vcf> -s <sample1> <sample2> -r <sample3> <sample4> -d polyphen --maf 0.01 -o output.vcf -l genelist.txt
    #as above but only consider missense variants predicted damaging by polyphen and with a minor allele frequency less than 1%. 

    findBiallelicVep.pl -i <variants.vcf> -s <sample1> <sample2> -e -o output.vcf -l genelist.txt
    #find genes with biallelic variants in two related samples where you expect them to share the same causative variant.

    findBiallelicVep.pl -i <variants.vcf> -f families.ped -o output.vcf -l genelist.txt -q 30
    #use families.ped file to describe affected and unaffected individuals, only consider calls with genotype quality of 30 or higher
    
    findBiallelicVep.pl -i <variants.vcf> -f families.ped -o output.vcf -l genelist.txt -w 10 -z 30
    #use a low gneotype quality threshold (10) to identify variants in affected samples but use a higher threshold (30) to identify genotypes in unaffecteds to reject

=cut

=head1 DESCRIPTION

This program reads VCF files annotated with Ensembl's Variant Effect Predictor and identifies transcripts with potential biallelic variation matching the various options specified above for the purpose of identifying potential recessively inherited pathogenic variants.  When more than one sample is specified using the --samples (-s) argument transcripts are identified that contain (not necessarily identical) potential biallelic variation in all samples. If multiple samples are specified in a PED file passed to the script with the --family (-f) argument, the program will attempt to find identical biallelic combinations within families and transcripts that contain potential biallelic variation in all affected samples from different families.

Genes are considered to contain potential biallelic variation if they either contain homozygous variants or two or more heterozygous variants. Phase can not be determined for variants so variants in cis may be erroneously considered to be potential biallelic variation.  Using variant data from unaffected parents with the --reject (-r) option or by specifying a PED file with the --family (-f) option  can help get around this issue.  Any samples specified using the --reject option will be used to remove biallelic variant combinations present - specifically, genotype combinations identified in affected samples (--samples) but also present in samples specified using the --reject argument will be removed from output. In this manner, if you have data from unaffected parents you should be able to avoid the problem of false positives from variants in cis as well as removing any shared genuine biallelic but non-pathogenic variation. However, this program does not require parental samples and can attempt to infer phase from a single parent if only one is available if the relationship is specified in a PED file passed to the script with the --family (-f) argument.

While related samples will be most useful in filtering using the --reject argument, data from any unaffected sample can be used to remove apparently non-pathogenic biallelic variation. Furthermore, while unrelated affected individuals can be used to identify shared genes containing apparent biallelic variation (when you believe the disorder to be caused by variation in the same gene), if using several related affected individuals you may use the --equal_genotypes flag to tell the program to only look for variants that are shared among all affected individuals AND potentially biallelic.

Note that this program only considers autosomal recessive disease, it ignores X, Y and mitochondrial chromosomes.



=cut

=head1 AUTHOR

David A. Parry

University of Leeds

=head1 COPYRIGHT AND LICENSE

Copyright 2014, 2013, 2012  David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut

