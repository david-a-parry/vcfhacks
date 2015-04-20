#!/usr/bin/perl
#
use warnings;
use strict; 
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use Term::ProgressBar;
use FindBin;
use lib "$FindBin::Bin/lib";
use ParseVCF;

my $help;
my $manual;
my $infile;
my $outfile;
my @classes;
my $no_head;
my $progress;
my $pass;
my $total_lines;
my @add;
my $canonical_only;
my $gmaf;#filter GMAF if present if equal to or above this value (0.0 - 0.5)
my $any_maf;#filter MAF if present if equal to or above this value in any population
my @damaging;
my $keep_any_damaging;
my $filter_unpredicted;
my @gene_lists; #array of files containing Ensembl Gene IDs to filter
my @samples;
my $in_every_sample;
my $matching_genes = 0;#will change to '' if user specifies --find_shared_genes but provides no argument, in which case print to STDERR
my $splice_consensus = 0; #use this flag to check for SpliceConsensus VEP plugin annotations
my %opts = ("help" => \$help,
            "manual" => \$manual,
            "pass" => \$pass,
            "progress" => \$progress,
            "input" => \$infile,
            "output" => \$outfile,
            "classes" => \@classes,
            "additional_classes" => \@add,
            "remove_headers" => \$no_head,
            "canonical_only" => \$canonical_only,
            "maf" => \$any_maf,
            "gmaf" => \$gmaf,
            "damaging" => \@damaging,
            "keep_any_damaging" => \$keep_any_damaging,
            "unpredicted_missense" => \$filter_unpredicted,
            "list" => \@gene_lists,
            "samples" => \@samples,
            "each_sample" => \$in_every_sample,
            "find_shared_genes" => \$matching_genes,
            "consensus_splice_site" => \$splice_consensus,
             );
GetOptions(\%opts,
            "help" ,
            "manual" ,
            "pass" ,
            "progress" ,
            "input=s" ,
            "output=s" ,
            "classes=s{,}" ,
            "additional_classes=s{,}" ,
            "remove_headers" ,
            "canonical_only" ,
            "maf=f" ,
            "gmaf=f" ,
            "damaging=s{,}" ,
            "keep_any_damaging" ,
            "unpredicted_missense" ,
            "list=s{,}",
            "samples=s{,}",
            "each_sample",
            "find_shared_genes:s", 
            "consensus_splice_site",
            ) or pod2usage(-message => "Syntax error", -exitval => 2);
pod2usage(-verbose => 2) if ($manual);
pod2usage(-verbose => 1) if ($help);
pod2usage(-message => 
          "Syntax error - an input file must be supplied.\n", 
          -exitval => 2) if (not $infile);
pod2usage(-message => 
        "--gmaf option requires a value between 0.00 and 0.50 to filter on global minor allele frequency.\n", 
        -exitval => 2) if (defined $gmaf && ($gmaf < 0 or $gmaf > 0.5));
if ($matching_genes){
    pod2usage(-message => 
    "--find_shared_genes option requires at least one sample to be specified with the --samples argument.\n", 
    -exitval => 2) if not @samples;
}elsif ($matching_genes eq ''){
    pod2usage(-message => 
    "--find_shared_genes option requires at least one sample to be specified with the --samples argument.\n", 
    -exitval => 2) if not @samples;
}

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

if (not @classes){

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
push (@classes, @add) if (@add);
push @classes, "splice_region_variant" if $splice_consensus;
foreach my $class (@classes){
    die "Error - variant class '$class' not recognised.\n" if not grep {/^$class$/i} @valid;
}


my $OUT;#output filehandle
if ($outfile){ 
    open ($OUT, ">$outfile") || die "Can't open file $outfile for writing.\n";
}else{
    $OUT = \*STDOUT;
}
my $LIST; #output filehandle for matching genes
if ($matching_genes){
    open ($LIST, ">$matching_genes") || die "Can't open $matching_genes for writing: $!\n";
}elsif($matching_genes eq ''){#user specified --list option but provided no argument
    $LIST = \*STDERR;
    $matching_genes = 1;
}
my @csq_fields = qw (allele gene feature feature_type consequence hgnc);#default fields to retrieve from CSQ INFO field

my %damage_filters = (); #hash of prediction program names and values to filter

if (@damaging){
    my %valid_damaging = (sift => ["deleterious", "tolerated"],  polyphen => ["probably_damaging", "possibly_damaging", "benign", "unknown"], condel => ["deleterious", "neutral"]);
    my %default_damaging = (sift => ["deleterious", ],  polyphen => ["probably_damaging", "possibly_damaging",], condel => ["deleterious", ]);
    foreach my $d (@damaging){
        my ($prog, $label) = split("=", $d);
        if (exists $valid_damaging{lc$prog}){
            push @csq_fields, lc($prog);
            no warnings 'uninitialized';
            my @filters = add_damaging_filters(lc$prog, lc$label, \%valid_damaging, \%default_damaging);
            push @{$damage_filters{lc$prog}}, @filters;
        }elsif (lc($prog) eq 'all'){
            %damage_filters = %default_damaging;
            push @csq_fields, qw(sift polyphen condel);
        }else{
            die "Unrecognised value ($d) passed to --damaging argument. See --help for more info.\n";
        }
    }
}else{
    if ($keep_any_damaging and $filter_unpredicted){
        die "--keep_any_damaging and --unpredicted_missense arguments can only be used when --damaging argument is in effect.\n";
    }elsif ($keep_any_damaging){
        die "--keep_any_damaging argument can only be used when --damaging argument is in effect.\n";
    }elsif ($filter_unpredicted){
        die "--unpredicted_missense argument can only be used when --damaging argument is in effect.\n";
    }
}

if ($canonical_only){
        push @csq_fields, 'canonical';
}
if (defined $gmaf or defined $any_maf){
    push @csq_fields, 'gmaf';
}
if ($splice_consensus){
    push @csq_fields, 'splice_consensus';
}
my @genes_to_filter;
if (@gene_lists){
    foreach my $gene_file (@gene_lists){
        open (my $GENE, $gene_file) or die "Can't open gene list file '$gene_file': $!\n";
        chomp (my @head = split("\t", <$GENE>));
        no warnings 'uninitialized';
        my $id_field = 0;
        $id_field++ until $head[$id_field] eq 'Ensembl Gene ID' or $id_field > $#head;
        if ($id_field > $#head){
            die "Can't find 'Ensembl Gene ID' field in header of gene list file '$gene_file'.\n";
        }
        while (chomp (my $line = <$GENE>)){
            my @split = split("\t", $line);
            push @genes_to_filter, $split[$id_field];
        }
    }
    my %seen = ();
    @genes_to_filter = grep { ! $seen{$_}++ } @genes_to_filter;
    @genes_to_filter = sort @genes_to_filter;
}
my $file = 0;
my $got = 0;
my $count = 0 ;
my %genes_per_sample = (); # will take on the format of $genes_per_sample{sample} = anonymous hash with keys = transcript id values = gene symbol
                            # e.g. $genes_per_sample{sample}->{ENST0012345} = {AGENE}
my @available_mafs = ();
my $vcf_obj = ParseVCF->new(file=> $infile);
my $vep_header = $vcf_obj->readVepHeader();
die "No 'consequence' field identified in header for file $infile - " .
"please annotate with Ensembl's variant_effect_precictor.pl script.\n" if (not exists $vep_header->{consequence});
die "No GMAF field in header for file $infile - please annotate GMAF with Ensembl's variant_effect_precictor.pl script.\n" 
    if (defined $gmaf and not exists $vep_header->{gmaf});
if (defined $any_maf){
    foreach my $key (keys %{$vep_header}){
        if ($key =~ /\w_MAF$/){
            push @available_mafs, $key;
            push @csq_fields, $key;
        }
    }
}

if (@csq_fields > 1){
    my %seen = ();
    @csq_fields = grep { ! $seen{$_}++ } @csq_fields;
}

foreach my $c (@csq_fields){
    if (not exists $vep_header->{$c}){
        die "Couldn't find '$c' VEP field in header - please ensure your VCF is annotated with " .
        "Ensembl's variant effect precictor specifying the appropriate annotations.\n";
    }
}

my $progressbar;
my $line_count = 0;
my $next_update = 0;
if ($progress){
    if ($infile eq "-"){
        print STDERR "Can't use --progress option when input is from STDIN\n";
        $progress = 0;
    }else{
        $progressbar = Term::ProgressBar->new({name => "Retrieving:", count => $vcf_obj->countLines("variants"), ETA => "linear", });
    }
}
unless ($no_head){
    print $OUT  $vcf_obj->getHeader(0) ."##getFunctionalVariantsVep.pl=\"";
    my @opt_string = ();
    foreach my $k (sort keys %opts){
        if (not ref $opts{$k}){
            push @opt_string, "$k=$opts{$k}";
        }elsif (ref $opts{$k} eq 'SCALAR'){
            if (defined ${$opts{$k}}){
                push @opt_string, "$k=${$opts{$k}}";
            }else{
                push @opt_string, "$k=undef";
            }
        }elsif (ref $opts{$k} eq 'ARRAY'){
            if (@{$opts{$k}}){
                push @opt_string, "$k=" .join(",", @{$opts{$k}});
            }else{
                push @opt_string, "$k=undef";
            }
        }
    }
    print $OUT join(" ", @opt_string) . "\"\n" .  $vcf_obj->getHeader(1);
}
LINE: while (my $line = $vcf_obj->readLine){
    $count++;
    $line_count++;
    my $printed_line = 0;
    if ($progress){
        $next_update = $progressbar->update($line_count) if $line_count >= $next_update;
    }
    if ($pass){
        next LINE if $vcf_obj->getVariantField("FILTER") ne 'PASS';
    }
    #get vep consequence fields
    my @csq = $vcf_obj->getVepFields(\@csq_fields); #returns array of hashes e.g. $csq[0]->{Gene} = 'ENSG000012345' ; $csq[0]->{Consequence} = 'missense_variant'
    die "No consequence field found for line:\n$line\nPlease annotated your VCF file with ensembl's variant effect precictor before running this program.\n" if not @csq;
    #assign each ALT genotype to a VEP allele (VEP cuts off the first nt of insertion or deletion alt alleles)
    my %vep_alleles = ();#hash of VEP alleles to VCF's ALT alleles
    my @alts = $vcf_obj->readAlleles(alt_alleles => 1);
    my $ref = $vcf_obj->getVariantField("REF");
    my @v_all = $vcf_obj->altsToVepAllele( ref => $ref, alt => \@alts); 
    @vep_alleles{ @v_all } = @alts;#$vep{$allele} now corresponds to GT value 

    #START FILTERING on CSQ fields
    #check whether canonical transcript
ANNOT: foreach my $annot (@csq){
    my @anno_csq = split(/\&/, $annot->{consequence});
    #skip NMD transcripts
    if (grep {/NMD_transcript_variant/i} @anno_csq){
        next ANNOT;
    }
    #check against gene filter list
    if (@genes_to_filter){#this will be empty if --list argument wasn't specified
        my $i = binsearch($annot->{gene}, \@genes_to_filter);
        if ($i > -1){
            next ANNOT;
        }
    }
    if($canonical_only){
       next if (not $annot->{canonical});
    }
    #check minor allele frequency
    if (defined $gmaf){
        if ($annot->{gmaf}){
            if ($annot->{gmaf} =~ /\w+:(\d\.\d+)/){
                next if $1 >= $gmaf;
            }
        }
    }
    if (defined $any_maf){
        foreach my $some_maf (@available_mafs){
            if ($annot->{$some_maf}){
                if ($annot->{$some_maf} =~ /\w+:(\d\.\d+)/){
                    next if $1 >= $any_maf;
                }
            }
        }
    }   
                
    #check whether annotation is present in @samples if @samples specified
    my @found_in_sample = ();       
    foreach my $sample (@samples){
        my @sample_alleles = $vcf_obj->getSampleActualGenotypes(sample => $sample, return_alleles_only => 1);
        if ($in_every_sample){
            next ANNOT if not grep {/^$vep_alleles{$annot->{allele}}$/} @sample_alleles;
                push @found_in_sample, $sample;
            }else{
                if (grep {/^$vep_alleles{$annot->{allele}}$/} @sample_alleles){
                    push @found_in_sample, $sample;
                }
            }
        }
        if (@samples){
            next ANNOT if not @found_in_sample; 
        }

        #check vep consequence class
CLASS:  foreach my $class (@classes){
            foreach my $ac (@anno_csq){
                if (lc$ac eq lc$class){
                    if (lc$class eq 'missense_variant' and %damage_filters){
                        next if (filter_missense($annot, \%damage_filters, $keep_any_damaging, $filter_unpredicted));
                    }elsif(lc$class eq 'splice_region_variant' and $splice_consensus){
                        my $consensus = $annot->{splice_consensus};
                        next if not $consensus;
                        if ($consensus !~/SPLICE_CONSENSUS\S+/i){
                            print STDERR "WARNING - SPLICE_CONSENSUS annotation $consensus is".
                            " not recognised as an annotation from the SpliceConsensus VEP plugin.\n";
                            next;
                        }
                    }
                    print $OUT "$line\n" if not $printed_line;
                    $printed_line++;
                    $got++;
                    #need to cycle through each annotation if checking for matching genes between samples
                    #otherwise we can skip to next line
                    if ($matching_genes){
                        my $symbol ;
                        if ($annot->{hgnc}){
                            $symbol = $annot->{hgnc};
                        }elsif($annot->{gene}){
                            $symbol = $annot->{gene};
                        }else{
                            $symbol = $annot->{feature};
                        }
                        foreach my $sample (@found_in_sample){
                            $genes_per_sample{$sample}->{$annot->{feature}} = $symbol;
                        }
                    }else{
                        next LINE;
                    }
                }
            }
        }
    }
}

close $OUT;
print STDERR "$got variants of $count identified.\n";
if ($matching_genes){
    print STDERR "Identifying matching genes between samples...\n";
    my %matches = ();#symbol = key, value = array of transcript ids
    my $transcript_count = 0;
TRANSCRIPT: foreach my $transcript ( keys %{$genes_per_sample{$samples[0]}}){
        #check all transcripts from first sample in all other samples
        for (my $i = 1; $i < @samples; $i++){
            if (not $genes_per_sample{$samples[$i]}->{$transcript}){
                next TRANSCRIPT;
            }
        }
        #if we haven't hit 'next TRANSCRIPT' gene must be shared in all samples
        $transcript_count++;
        push @{$matches{$genes_per_sample{$samples[0]}->{$transcript}}}, $transcript;
    }
    if (keys %matches){
        foreach my $gene (sort keys %matches){
            my %seen = ();
            @{$matches{$gene}} = grep {! $seen{$_}++ } @{$matches{$gene}};#remove duplicates
            print $LIST "$gene:" . join(",", @{$matches{$gene}}) . "\n";
        }
        print STDERR "$transcript_count matching transcripts found between samples.\n";
        print STDERR scalar (keys %matches) . " matching genes found between samples.\n";
    }else{
        print STDERR "0 matching transcripts found between samples.\n";
    }
}

###################
###########
sub filter_missense{
#returns 1 if missense should be filtered
#otherwise returns 0;
# uses $keep_any setting to return 0 if any of these predictions match, otherwise all available
# scores must be deleterious/damaging 
# if $filter_missing is used a variant will be filtered if no score is available (overriden by $keep_any setting)
    my ($anno, $filter_hash, $keep_any, $filter_missing) = @_;
#my %default_damaging = (sift => ["deleterious", ],  polyphen => ["probably_damaging", "possibly_damaging",], condel => ["deleterious", ]);
    my %filter_matched = ();
PROG: foreach my $k (sort keys %$filter_hash){
        my $score = $anno->{lc$k};
        if (not $score or $score =~ /^unknown/i){ #don't filter if score is not available for this variant unless $filter_missing is in effect
            $filter_matched{$k}++ unless $filter_missing;
            next;
        }
SCORE:  foreach my $f (@{$filter_hash->{$k}}){
            if ($f =~ /^\d(\.\d+)*$/){
                my $prob;
                if ($score =~ /^(\d(\.\d+)*)/){
                    $prob = $1;
                }else{
                    next SCORE;#if score not available for this feature ignore and move on 
                }
                if (lc$k eq 'polyphen'){
                    if ($prob >= $f){#higher is more damaging for polyphen - damaging
                        return 0 if $keep_any;
                        $filter_matched{$k}++;
                        next PROG;
                    }else{#benign
                    }
                }else{
                    if ($prob <= $f){#lower is more damaging for sift and condel - damaging
                        return 0 if $keep_any;
                        $filter_matched{$k}++;
                        next PROG;
                    }else{#benign
                    }
                }
            }else{
                $score =~ s/\(.*\)//;
                if (lc$f eq lc$score){#damaging
                    return 0 if $keep_any;
                    $filter_matched{$k}++;
                    next PROG;
                }
        }
        }

    }
    foreach my $k (sort keys %$filter_hash){
        #filter if any of sift/condel/polyphen haven't matched our deleterious settings
        return 1 if not exists $filter_matched{$k};
    }
    return 0;
}

###########
sub add_damaging_filters{
    my ($prog, $label, $valid_hash, $default_hash) = @_;
    if ($label){
        if ($label =~ /^\d(\.\d+)*$/){
            die "Numeric values for condel, polyphen or sift filtering must be between 0 and 1.\n" if ( $label < 0 or $label > 1);
            return $label;
        }else{
            my @lb = split(",", $label);
            foreach my $l (@lb){
                die "Invalid filter parameter '$l' used with --damaging argument.  See --help for valid arguments.\n" if not grep{/^$l$/} @{$valid_hash->{$prog}};
            }
            return @lb;
        }

    }else{#give default values for $prog
        return @{$default_hash->{$prog}};
    }
}



###########
sub binsearch{
    my ($x, $ar) = @_;
    my $u = $#{$ar};
    my $l = 0;
    while ($l <= $u){
        my $i = int(($l + $u)/2);
        if ($x lt $ar->[$i]){
            $u = $i -1;
        }elsif ($x gt $ar->[$i]){
            $l = $i +1;
        }else{
            return $i;
        }
    }
    return -1;
}
=head1 NAME

getFunctionalVariantsVep.pl  -  retrieve specific variant classes from an annotated vcf file.

=head1 SYNOPSIS

=over

=item getFunctionalVariantsVep.pl --input <annotated vcf file> [options]

=item getFunctionalVariantsVep.pl --help for help message

=item getFunctionalVariantsVep.pl --manual for manual page

=back
=cut

=head1 ARGUMENTS

=over 8

=item B<-i    --input>

VCF file with functional annotations from variant_effect_predictor.pl script.  

=item B<-o    --output>

File to write output. Default is STDOUT.

=item B<--canonical_only>

Only look at consequences for canonical transcripts.

=item B<--classes>

If you wish to specify a custom set of variant classes to retrieve enter them here separated by spaces.

Default classes are:
          
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

Specify one or more classes, separated by spaces, to add to the default mutation classes used.
 
=item B<--consensus_splice_site>

Use this flag in order to keep splice_region_variant classes only if they are in a splice consensus region as defined by the SpliceConsensus plugin. You do not need to specify 'splice_region_variant' using --classes or --add_classes options when using this flag. You B<MUST> have used the SpliceConsensus plugin when running the VEP for this option to work correctly.

=item B<-g    --gmaf>

Use a value between 0.00 and 0.50 to specify global minor allele frequencey filtering. If GMAF is available for variant it will be filtered if equal to or greater than the value specfied here.

=item B<--maf>

Like gmaf but filter on any population specific minor allele frequency annotated by the VEP as well as the GMAF.

=item B<-d    --damaging>

Specify SIFT, PolyPhen or Condel labels or scores to filter on. Add the names of the programs you want to use, separated by spaces, after the --damaging option. By default SIFT will keep variants labelled as 'deleterious', Polyphen will keep variants labelled as 'possibly_damaging' or 'probably_damaging' and  Condel will keep variants labelled as 'deleterious'.

If you want to filter on custom values specify values after each program name in the like so: 'polyphen=probably_damaging'. Separate multiple values with commas - e.g. 'polyphen=probably_damaging,possibly_damaging,unknown'. You may specify scores between 0 and 1 to filter on scores rather than labels - e.g. 'sift=0.3'. For polyphen, variants with scores lower than this score are considered benign and filtered, for SIFT and Condel higher scores are considered benign.

=over

=item Valid labels for SIFT: deleterious, tolerated

=item Valid labels for Polyphen: probably_damaging, possibly_damaging, benign, unknown

=item Valid labels for Condel : deleterious, neutral


=back

To use default values for all three programs use 'all' (i.e. '--damaging all').

The default behaviour is to only keep variants predicted as damaging by ALL programs specified, although if the value is not available for one or more programs than that program will be ignored for filtering purposes.  For best results/maximum flexibility it is recommended that you annotate both 'scores' and 'predictions' when generating your input VCF with variant_effect_predictor.pl.


=item B<-k    --keep_any_damaging>

If using multiple programs with the --damaging argument use this flag to keep variants predicted to be damaging according to ANY of these programs.

=item B<-u    --unpredicted_missense>

Skip variants that do not have a score from one or more programs specified by the --damaging argument. The --keep_any_damaging argument will override this behaviour if any of the available predictions are considered damaging.

=item B<-l    --list>

File containing a list of genes that should be filtered (e.g. LOF-tolerant genes). The program expects a tab delimited list with a header line specifying 'Ensembl Gene ID' for the column containing gene identifiers. 

=item B<-r    --remove_headers>

Use this flag to prevent printing of headers to output. 

=item B<--pass>

Keep only variants passing filters (i.e. FILTER field is PASS).

=item B<-s    --samples>

Only keep variants present in at least one of these samples.  Can change behaviour to keep only variants present in ALL of these samples using the --each_sample switch.

=item B<-e    --each_sample>

When --samples arguments are specified this switch will cause the program to only return variants present in ALL of the samples specified by the --samples argument.

=item B<-f    --find_shared_genes>

If --samples argument is specified this switch will return a list of genes shared by each sample according to the filtering criteria.  If a filename is given the list will be printed to file, otherwise the list will be printed to STDERR.

=item B<--progress>

Display progress bar. 

=item B<-h    --help>

Display help message.

=item B<--manual>

Show manual page.

=back
=cut

=head1 EXAMPLES

=over

=item getFunctionalVariantsVep.pl -i var.vcf 

(filter on default set of consequences, print output to STDOUT.)


=item getFunctionalVariantsVep.pl -i var.vcf --classes stop_gained stop_lost

(filter on 'stop_gained' and 'stop_lost' classes only.)


=item getFunctionalVariantsVep.pl -i var.vcf --maf 0.01 

(filter on default set of consequences, filter anything with a MAF greater than or equal to 0.01 [1 %].)


=item getFunctionalVariantsVep.pl -i var.vcf --maf 0.01 --damaging all 

(As above but only keep missense variation predicted damaging by all annotation programs.)


=item getFunctionalVariantsVep.pl -i var.vcf --maf 0.01 --damaging all --list gene_list.txt 

(As above but also filter genes in gene_list.txt file.)


=item getFunctionalVariantsVep.pl -i var.vcf --maf 0.01 --damaging all --list gene_list.txt -o output.vcf

(As above but specifying output file)

=back 

=head1 DESCRIPTION

In its simplest form this program will print specific variant classes from a VCF file annotated with Ensembl's variant_effect_predictor.pl program and filter out others. Input must be a VCF annotated by the variant_effect_predictor.pl program using the '--vcf' option, annotations are read from the INFO field of the VCF. By default the following classes are kept:
          
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

You can add extra classes using the --add option or specify a totally different set of classes using the --classes option. Options are also available to additionally filter on Polyphen/SIFT/Condel prediction scores, GMAF, canonical transcripts or gene lists as detailed above.  

=head1 AUTHOR

David A. Parry

University of Leeds

=head1 COPYRIGHT AND LICENSE

Copyright 2013  David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

