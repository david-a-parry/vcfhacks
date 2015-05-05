#!/usr/bin/perl

use warnings;
use strict; 
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use Scalar::Util qw(looks_like_number);
use Term::ProgressBar;
use FindBin;
use lib "$FindBin::Bin/lib";
use ParseVCF;
use DbnsfpVcfFilter;

my $help;
my $manual;
my $infile;
my $outfile;
my @classes;#Effect classes to keep
my @add;#additional classes if adding to default @classes
my @effect_impact; #Effect_Impact classes to keep (High, Moderate, Low, Modifier)
my $no_head;
my $progress;
my $pass;
my $total_lines;
my @filter_dbnsfp;#array of strings in format "[dbNSFP field][comparator][value] e.g. "1000Gp1_EUR_AF<=0.01" or "genename=MUC4"
my @match_dbnsfp; #as above
my $multi_value_operator = 'or';#if a dbNSFP field has more than one value by default we use or to match any
my @damaging; #array of missense prediction program names to create dbNSFP expressions for 
my $gmaf;#filter 1000 genome allele frequency if present if equal to or above this value (0.0 - 1.00)
my $any_maf;#filter on Alelle Frequency if present if equal to or above this value in any population
my $keep_any_damaging;
my $filter_unpredicted;
my @gene_lists; #array of files containing Ensembl Gene IDs to filter
my @samples;
my $in_every_sample;
my $matching_genes = 0;#will change to '' if user specifies --find_shared_genes but provides no argument, in which case print to STDERR
my $gene_id_mode;
my %opts = ("help" => \$help, 
            "manual" => \$manual,
            "pass" => \$pass,
            "progress" => \$progress,
            "input" => \$infile,
            "output" => \$outfile,
            "gene_id_mode" => \$gene_id_mode,
            "classes" => \@classes,
            "additional_classes" => \@add,
            "effect_impact" => \@effect_impact,
            "remove_headers" => \$no_head,
            "allele_frequency" => \$any_maf,
            "1000_genomes_allele_frequency" => \$gmaf,
            'filter_dbnsfp' => \@filter_dbnsfp,
            'match_dbnsfp' => \@match_dbnsfp,
            'multi_value_operator' => \$multi_value_operator, 
            'damaging' => \@damaging,
            'keep_any_damaging' => \$keep_any_damaging,
            'list' => \@gene_lists,
            'samples' => \@samples,
            'each_sample' => \$in_every_sample,
            'find_shared_genes' => \$matching_genes );

GetOptions(\%opts, 
            "help" ,
            "manual" ,
            "pass" ,
            "progress" ,
            "input=s" ,
            "output=s" ,
            "classes=s{,}" ,
            "additional_classes=s{,}" ,
            "effect_impact=s{,}",
            "remove_headers" ,
            "allele_frequency=f" ,
            "1000_genomes_allele_frequency=f" ,
            'f|filter_dbnsfp=s{,}' =>  \@filter_dbnsfp,
            'match_dbnsfp=s{,}' ,
            'multi_value_operator=s',
            'damaging=s{,}' ,
            'keep_any_damaging' ,
            'unpredicted_missense' ,
            'list=s{,}',
            'samples=s{,}',
            'each_sample',
            'x|find_shared_genes:s' => \$matching_genes,
            ) or pod2usage(-message => "Syntax error",  -exitval => 2);

pod2usage(-verbose => 2) if ($manual);
pod2usage(-verbose => 1) if ($help);
pod2usage(-message => 
          "Syntax error - an input file must be supplied.\n", 
          -exitval => 2) if (not $infile);
pod2usage(-message => 
        "--maf option requires a value between 0.00 and 1.o0 to filter on global minor allele frequency.\n", 
        -exitval => 2) if (defined $any_maf && ($any_maf < 0 or $any_maf > 1.0));
pod2usage(-message => 
        "--gmaf option requires a value between 0.00 and 1.o0 to filter on global minor allele frequency.\n", 
        -exitval => 2) if (defined $gmaf && ($gmaf < 0 or $gmaf > 1.0));
if ($matching_genes){
    pod2usage(-message => 
    "--find_shared_genes option requires at least one sample to be specified with the --samples argument.\n", 
    -exitval => 2) if not @samples;
}elsif ($matching_genes eq ''){
    pod2usage(-message => 
    "--find_shared_genes option requires at least one sample to be specified with the --samples argument.\n", 
    -exitval => 2) if not @samples;
}
my @csq_fields = qw (
                Effect
                Effect_Impact
                Functional_Class
                Gene_Name
                Transcript_ID
                Transcript_BioType
                Genotype_Number);#default fields to retrieve from EFF INFO field

if (@effect_impact){
    checkEffectImpacts(@effect_impact);
    #if --effect_impact is used we don't use @classes
    #unless they've been specified by user
    if (@classes or @add){
	    push @classes, @add;
        checkAndParseClasses(\@classes);
    }

}else{
    checkAndParseClasses(\@classes, \@add);
}
my @genes_to_filter = checkGeneList(@gene_lists) if @gene_lists;

my $got = 0;
my $count = 0 ;
my %genes_per_sample = (); # will take on the format of $genes_per_sample{sample} = anonymous hash with keys = transcript id values = gene symbol
                            # e.g. $genes_per_sample{sample}->{ENST0012345} = {AGENE}
my @available_mafs = ();
print STDERR "Initializing input vcf...\n";
my $vcf_obj = ParseVCF->new(file=> $infile);
my $meta_header = $vcf_obj->getHeader(0);
my $snp_eff_header = $vcf_obj->readSnpEffHeader();
die "No 'Effect' field identified in header for file $infile - " .
"please annotate with SnpEff.\n" if (not exists $snp_eff_header->{Effect});
if($gmaf){
    push @filter_dbnsfp, "dbNSFP_1000Gp1_AF<=$gmaf";
}
my @filters = checkAndParseDbnsfpExpressions(\@filter_dbnsfp, $meta_header);
my @match = checkAndParseDbnsfpExpressions(\@match_dbnsfp, $meta_header);
if ($any_maf){
    push @filters, getFrequencyFilters( $any_maf, ">=",);
}

my @missense_filters = ();
if (@damaging){
    if ($keep_any_damaging){
        if (grep {/^all$/i} @damaging){
            push @missense_filters, getAnyDamagingFilters();
        }elsif (grep {/^any$/i} @damaging){
            push @missense_filters, getAnyDamagingFilters();
        }else{
            push @missense_filters, getAnyDamagingFilters(undef, \@damaging);
        }
    }else{
        if (grep {/^all$/i} @damaging){
            push @missense_filters, getAllDamagingFilters();
        }elsif (grep {/^any$/i} @damaging){
            push @missense_filters, getAnyDamagingFilters();
        }else{
            push @missense_filters, getAllDamagingFilters(undef, \@damaging);
        }
    }
}

#open output files
my $OUT;#output filehandle, defaults to STDOUT
if ($outfile){ 
    open ($OUT, ">$outfile") || die "Can't open file $outfile for writing.\n";
}else{
    $OUT = \*STDOUT;
}
my $LIST; #output filehandle for matching genes, defaults to STDERR
if ($matching_genes){
    open ($LIST, ">$matching_genes") || die "Can't open $matching_genes for writing: $!\n";
}elsif($matching_genes eq ''){#user specified --list option but provided no argument
    $LIST = \*STDERR;
    $matching_genes = 1;
}


my $progressbar;
my $line_count = 0;
my $next_update = 0;
my $total_variants = 0;
if ($progress){
    if ($infile eq "-"){
        print STDERR "Can't use --progress option when input is from STDIN\n";
        $progress = 0;
    }else{
        $total_variants = $vcf_obj->countLines("variants");
        $progressbar = Term::ProgressBar->new({
                name => "Retrieving:", 
                count => $total_variants, 
                ETA => "linear", 
            });
    }
}
unless ($no_head){
    print $OUT  $meta_header ."##getFunctionalVariantsSnpEff.pl=\"";
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

    #get SnpEff consequence fields
    my @csq = $vcf_obj->getSnpEffFields(\@csq_fields); #returns array of hashes e.g. $csq[0]->{Gene_Name} = 'ABCA1' ; $csq[0]->{EFFECT} = 'DOWNSTREAM'
    die "No consequence field found for line:\n$line\nPlease annotated your VCF file with SnpEff before running this program.\n" if not @csq;

    #START FILTERING on CSQ fields
ANNOT: foreach my $annot (@csq){
        #skip pseudogenes or NMD transcripts
        if ($annot->{Transcript_BioType} =~ /nonsense_mediated_decay|pseudogene$/i){
            next ANNOT;
        }
        #check against gene filter list
        if (@genes_to_filter){#this will be empty if --list argument wasn't specified
            my $i = binsearch($annot->{Gene_Name}, \@genes_to_filter);
                if ($i > -1){
                    next ANNOT;
            }
        }
                
        #check whether genotype is present in @samples if @samples specified
        my @found_in_sample = ();       
        foreach my $sample (@samples){
            my @sample_alleles = $vcf_obj->getSampleCall(sample => $sample, return_alleles_only => 1);
            if ($in_every_sample){
                next ANNOT if not grep {/^$annot->{Genotype_Number}$/} @sample_alleles;
                push @found_in_sample, $sample;
            }else{
                if (grep {/^$annot->{Genotype_Number}$/} @sample_alleles){
                    push @found_in_sample, $sample;
                }
            }
        }
        if (@samples){
            next ANNOT if not @found_in_sample; 
        }

        #do dbNSFP filtering
        if (@filters){
            if (filterDbnsfpExpressions(\@filters, $line, $multi_value_operator)){
                next LINE;
            }
        }
        if (@match){
            if (not filterDbnsfpExpressions(\@match, $line, $multi_value_operator )){
                next LINE;
            }
        }
        if ($annot->{Effect} eq 'NON_SYNONYMOUS_CODING'){
            if (@missense_filters){
                if (filterDbnsfpExpressions(\@missense_filters, $line, $multi_value_operator)){
                    next ANNOT;
                }
            }
        }
        #finally, check SnpEff consequence 
CLASS:  foreach my $class (@classes){
            if (lc$annot->{Effect} eq lc$class){
                print $OUT "$line\n" if not $printed_line;
                $printed_line++;
                $got++;
                #need to cycle through each annotation if checking for matching genes between samples
                #otherwise we can skip to next line
                if ($matching_genes){
                    my $symbol;
                    if ($annot->{Gene_Name}){
                        $symbol = $annot->{Gene_Name};
                    }else{
                        $symbol = $annot->{Transcript_ID};
                    }
                    foreach my $sample (@found_in_sample){
                        $genes_per_sample{$sample}->{$annot->{Transcript_ID}} = $symbol;
                    }
                }else{
                    next LINE;
                }
            }
        }
EFFECT_IMPACT:  foreach my $impact (@effect_impact){
            if (lc$annot->{Effect_Impact} eq lc$impact){
                print $OUT "$line\n" if not $printed_line;
                $printed_line++;
                $got++;
                #need to cycle through each annotation if checking for matching genes between samples
                #otherwise we can skip to next line
                if ($matching_genes){
                    my $symbol;
                    if ($annot->{Gene_Name}){
                        $symbol = $annot->{Gene_Name};
                    }else{
                        $symbol = $annot->{Transcript_ID};
                    }
                    foreach my $sample (@found_in_sample){
                        $genes_per_sample{$sample}->{$annot->{Transcript_ID}} = $symbol;
                    }
                }else{
                    next LINE;
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


##################################################
sub checkGeneList{
    my @genes_to_filter = ();
    foreach my $gene_file (@_){
        open (my $GENE, $gene_file) or die "Can't open gene list file '$gene_file': $!\n";
        chomp (my @head = split("\t", <$GENE>));
        no warnings 'uninitialized';
        my $id_field = 0;

        #TO DO(?)
        #CHECK TO SEE WHETHER INPUT USES GENE SYMBOLS OR ENSEMBL IDs?

        #$id_field++ until $head[$id_field] eq 'Ensembl Gene ID' or $id_field > $#head;
        if ($gene_id_mode){
            $id_field++ until $head[$id_field] eq 'Ensembl Gene ID' or $id_field > $#head;
        }else{
            $id_field++ until $head[$id_field] eq 'Associated Gene Name' or $id_field > $#head;
        }
        if ($id_field > $#head){
            die "Can't find 'Ensembl Gene ID' field in header of gene list file '$gene_file'.\n" 
                if $gene_id_mode;
            die "Can't find 'Associated Gene Name' field in header of gene list file '$gene_file'.\n";
        }
        while (chomp(my $line = <$GENE>)){
            my @split = split("\t", $line);
            push @genes_to_filter, $split[$id_field];
        }
    }
    my %seen = ();
    @genes_to_filter = grep { ! $seen{$_}++ } @genes_to_filter;
    @genes_to_filter = sort @genes_to_filter;
}
##################################################
sub checkEffectImpacts{
    my @valid = qw (High Moderate Low Modifier);
    foreach my $impact (@_){
        die "Error - Effect_Impact '$impact' not recognised.\n" if not grep {/$impact/i} @valid;
    }
}
    
##################################################
sub checkAndParseClasses{
    my ($classes, $additional) = @_;

    my @valid = qw (INTERGENIC
                UPSTREAM
                UTR_5_PRIME
                UTR_5_DELETED
                START_GAINED
                SPLICE_SITE_ACCEPTOR
                SPLICE_SITE_DONOR
                START_LOST
                SYNONYMOUS_START
                CDS
                GENE
                TRANSCRIPT
                EXON
                EXON_DELETED
                NON_SYNONYMOUS_CODING
                SYNONYMOUS_CODING
                FRAME_SHIFT
                CODON_CHANGE
                CODON_INSERTION
                CODON_CHANGE_PLUS_CODON_INSERTION
                CODON_DELETION
                CODON_CHANGE_PLUS_CODON_DELETION
                STOP_GAINED
                SYNONYMOUS_STOP
                STOP_LOST
                INTRON               
                UTR_3_PRIME
                UTR_3_DELETED
                DOWNSTREAM
                INTRON_CONSERVED
                INTERGENIC_CONSERVED
                INTRAGENIC
                RARE_AMINO_ACID
                NON_SYNONYMOUS_START);

    if (not @$classes){
        @$classes = qw (
                UTR_5_DELETED
                SPLICE_SITE_ACCEPTOR
                SPLICE_SITE_DONOR
                START_LOST
                EXON_DELETED
                NON_SYNONYMOUS_CODING
                FRAME_SHIFT
                CODON_INSERTION
                CODON_CHANGE_PLUS_CODON_INSERTION
                CODON_DELETION
                CODON_CHANGE_PLUS_CODON_DELETION
                STOP_GAINED
                STOP_LOST
                UTR_3_DELETED
                RARE_AMINO_ACID
                );
    }
    push (@classes, @$additional) if (@$additional);
    foreach my $class (@$classes){
        die "Error - variant class '$class' not recognised.\n" if not grep {/^$class$/i} @valid;
    }
}


##################################################
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

##################################################
=head1 NAME

getFunctionalVariantsSnpEff.pl  -  retrieve specific variant classes from a SnpEff annotated vcf file.

=head1 SYNOPSIS

=over

=item getFunctionalVariantsSnpEff.pl --input <annotated vcf file> [options]

=item getFunctionalVariantsSnpEff.pl --help for help message

=item getFunctionalVariantsSnpEff.pl --manual for manual page

=back
=cut

=head1 ARGUMENTS

=over 8

=item B<-i    --input>

VCF file with functional annotations from SnpEff.jar.  

=item B<-o    --output>

File to write output. Default is STDOUT.

=item B<--classes>

If you wish to specify a custom set of variant classes (SnpEff "Effects") to retrieve enter them here separated by spaces.

Default classes are:

                UTR_5_DELETED
                SPLICE_SITE_ACCEPTOR
                SPLICE_SITE_DONOR
                START_LOST
                EXON_DELETED
                NON_SYNONYMOUS_CODING
                FRAME_SHIFT
                CODON_INSERTION
                CODON_CHANGE_PLUS_CODON_INSERTION
                CODON_DELETION
                CODON_CHANGE_PLUS_CODON_DELETION
                STOP_GAINED
                STOP_LOST
                UTR_3_DELETED
                RARE_AMINO_ACID


The user can specify one or more of the following classes instead: 

                INTERGENIC
                UPSTREAM
                UTR_5_PRIME
                UTR_5_DELETED
                START_GAINED
                SPLICE_SITE_ACCEPTOR
                SPLICE_SITE_DONOR
                START_LOST
                SYNONYMOUS_START
                CDS
                GENE
                TRANSCRIPT
                EXON
                EXON_DELETED
                NON_SYNONYMOUS_CODING
                SYNONYMOUS_CODING
                FRAME_SHIFT
                CODON_CHANGE
                CODON_INSERTION
                CODON_CHANGE_PLUS_CODON_INSERTION
                CODON_DELETION
                CODON_CHANGE_PLUS_CODON_DELETION
                STOP_GAINED
                SYNONYMOUS_STOP
                STOP_LOST
                INTRON
                UTR_3_PRIME
                UTR_3_DELETED
                DOWNSTREAM
                INTRON_CONSERVED
                INTERGENIC_CONSERVED
                INTRAGENIC
                RARE_AMINO_ACID
                NON_SYNONYMOUS_START

For details see http://snpeff.sourceforge.net/SnpEff_manual.html#eff
                
=item B<--additional_classes>

Specify one or more additional classes, separated by spaces, to add to the default mutation classes used.

=item B<--effect_impact>

Specify one or more "Effect_Impact" values to only keep variants with a SnpEff "Effect_Impact" matching one of these values. Valid values are:

                High
                Moderate
                Low
                Modifier

For details see http://snpeff.sourceforge.net/SnpEff_manual.html#eff

As suggested by the SnpEff documentation, this is essentially a shortcut for specifying multiple Effect classes with the --classes option. If this option is used no "Effect" values will be used for filtering unless specifically specified using the --classes or --additional_classes options.

=item B<-1    --1000_genomes_allele_frequency>

Use a value between 0.00 and 1.00 to specify allele frequencey filtering from values in the 1000 genomes phase 1 global allele frequency. If the allele frequency is available for a variant allele it will be filtered if equal to or greater than the value specfied here.

=item B<--allele_frequency>

As above but for any of the following allele frequencies:

=over

1000 genomes phase 1 Global 

1000 genomes phase 1 African descendent

1000 genomes phase 1 European descendent

1000 genomes phase 1 American descendent

1000 genomes phase 1 Asian descendent

NHLBI ESP6500 African American

NHLBI ESP6500 European American

=back

=item B<-f    --filter_dbnsfp>

Filtering expressions for dbNSFP annotations added using SnpSift.jar (http://snpeff.sourceforge.net/SnpSift.html#dbNSFP). Each expression must formatted as shown in the following examples:

=over

"dbNSFP_GERP++_RSE<lt>2"

"dbNSFP_SIFT_pred=T or dbNSFP_Polyphen2_HVAR_pred=B"

"dbNSFP_ESP6500_EA_AFE<gt>=0.01 and dbNSFP_ESP6500_AA_AFE<lt>=0.02"

=back

Each expression consists of a dbNSFP annotation (e.g. "dbNSFP_GERP++_RS"), a comparator (e.g. "<") and a value (e.g. "2"). Multiple expressions can be joined with 'and', 'or' or 'xor' operators if desired. If joined by 'and' BOTH expressions must be true to filter. If joined by 'or' if either expression (or both) are true the variant will be filtered. If joined by 'xor' if either expression (but not both) are true the variant will be filtered. Valid comparators are:

=over

= (or ==): equal to

! (or !=): not equal to

<: less than

<=: less than or equal to

>: more than

>=: more than or equal to

=back 

In the first example any line with a GERP Rejected Substitution score less than 2 will be filtered. In the second example any variant that SIFT predicts as 'tolerated' (i.e. the dbNSFP_SIFT_pred annotation is 'T') or predicted by Polyphen to be benign (i.e. the dbNSFP_Polyphen2_HVAR_pred annotation is 'B') will be filtered. In the third example a variant is filtered if the NHLBI ESP6500 European American allele frequency is greater than or equal to 0.01 (1 %) and the NHLBI ESP6500 African American allele frequency if less than or equal to 0.02 (2 %). 

If a variant doesn't have an annotation for a given expression it will not be filtered (you can use the --match_dbnsfp argument to only print lines matching an expression if you prefer). By default, if a variant has multiple values for a given annotation it will be filtered if any of the annotations match the given expression.  You can change this behaviour using the --multi_value_operator option.

You can specify as many expressions as you wish. B<Make sure to enclose each expression in quotes or you will run into trouble using E<gt> symbols and the like.> Descriptions of the various dbNSFP annotations can be found at https://sites.google.com/site/jpopgen/dbNSFP

See the DESCRIPTION section to see shortened names available for dbNSFP fields.

=item B<--match_dbnsfp>

As above except lines will be kept if they match this expression and filtered if they do not. For example:

"dbNSFP_Interpro_domain=GPCR"

...will only return variants in genes encoding GPCR Interpro domains.

=item B<--multi_value_operator>

By default, when using --filter_dbnsfp or --match_dbnsfp expressions when an annotation has multiple values, the expression will be considered to be true if ANY of those values match.  You may specify 'and', 'or' or 'xor' after the --multi_value_operator to specify how multiple values should be treated instead ('or' is the default, 'and' means all values must match and 'xor' means only one value should match).   

=item B<-d    --damaging>

Filter missense variants unless predicted to be damaging by specified prediction programs.  Available programs are LRT, MutationAssessor, MutationTaster, Polyphen2 (using the HVAR prediction) and SIFT. Specify one or more of these programs, 'all' for all of them or 'any' for any of them (equivelent to using 'all' with the --keep_any_damaging option). 

Using option without the --keep_any_damaging argument and without passing 'any' as an argument essentially acts as a shortcut for using the --filter_dbnsfp expressions:

=over

"dbNSFP_LRT_pred=N"

"dbNSFP_MutationAssessor_pred=L or dbNSFP_MutationAssessor_pred=N"

"dbNSFP_MutationTaster_pred=N or dbNSFP_MutationTaster_pred=P"

"dbNSFP_Polyphen2_HVAR_pred=B"

"dbNSFP_SIFT_pred=T"


=back

The default behaviour is to only REMOVE variants predicted as damaging by ALL programs specified, but this behaviour can be changed using the --keep_any_damaging argument.

=item B<-k    --keep_any_damaging>

If using multiple programs with the --damaging argument use this flag to keep variants predicted to be damaging according to ANY of these programs. Missenses without predictions by any of the programs will still be kept.

=item B<-l    --list>

File containing a list of genes that should be filtered (e.g. LOF-tolerant genes). The program expects a tab delimited list with a header line specifying 'Associated Gene Name' for the column containing gene symbols unless --gene_id flag is used, in which case 'Ensembl Gene ID' column is required. 

=item B<--g    --gene_id_mode>

Use this flag if you have used SnpEff with -geneId flag to annotate IDs rather than gene names and you are using the --list option.

=item B<-r    --remove_headers>

Use this flag to prevent printing of headers to output. 

=item B<--pass>

Keep only variants passing filters (i.e. FILTER field is PASS).

=item B<-s    --samples>

Only keep variants present in at least one of these samples.  Can change behaviour to keep only variants present in ALL of these samples using the --each_sample switch.

=item B<--each_sample>

When --samples arguments are specified this switch will cause the program to only return variants present in ALL of the samples specified by the --samples argument.

=item B<-x    --find_shared_genes>

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

=item getFunctionalVariantsSnpEff.pl -i var.vcf 

(filter on default set of consequences, print output to STDOUT)


=item getFunctionalVariantsSnpEff.pl -i var.vcf --classes STOP_LOST STOP_GAINED

(keep variants with 'STOP_GAINED' or 'STOP_LOST' classes only)


=item getFunctionalVariantsSnpEff.pl -i var.vcf --1000_genomes_allele_frequency 0.01 

(filter on default set of consequences, filter anything with a 1000 genomes pilot allele frequency greater than or equal to 0.01)


=item getFunctionalVariantsSnpEff.pl -i var.vcf --1000_genomes_allele_frequency 0.01 --damaging all 

(As above but only keep missense variation predicted damaging by all annotation programs)


=item getFunctionalVariantsSnpEff.pl -i var.vcf --1000_genomes_allele_frequency 0.01 --damaging all --list gene_list.txt 

(As above but also filter genes in gene_list.txt file)


=item getFunctionalVariantsSnpEff.pl -i var.vcf --1000_genomes_allele_frequency 0.01 --damaging all --list gene_list.txt -o output.vcf

(As above but specifying output file)

=item getFunctionalVariantsSnpEff.pl -i var.vcf --allele_frequency 0.01 --damaging all --list gene_list.txt -o output.vcf

(Filtering if any given population specific or global allele frequency is greater than or equal to 0.01)

=item getFunctionalVariantsSnpEff.pl -i var.vcf o output.vcf --filter_dbnsfp "dbNSFP_GERP++_RSE<lt>2"

(Filter anything with a GERP RS score below 2)

=item getFunctionalVariantsSnpEff.pl -i var.vcf -o output.vcf --filter_dbnsfp "dbNSFP_GERP++_RSE<lt>2" --match_dbnsfp "dbNSFP_Interpro_domain=GPCR"

(As above but only printing variants in a gene containing a GPCR Interpro domain)

=back 

=head1 DESCRIPTION

In its simplest form this program will print specific variant classes (SnpEff "Effect" annotations) from a VCF file annotated with Ensembl's SnpEff.jar program and filter out others. Input must be a VCF annotated by the variant_effect_predictor.pl program using the '--vcf' option, annotations are read from the INFO field of the VCF. By default the following SnpEff "Effect" classes are kept:
                
                UTR_5_DELETED
                START_GAINED
                SPLICE_SITE_ACCEPTOR
                SPLICE_SITE_DONOR
                START_LOST
                EXON_DELETED
                NON_SYNONYMOUS_CODING
                FRAME_SHIFT
                CODON_INSERTION
                CODON_CHANGE_PLUS_CODON_INSERTION
                CODON_DELETION
                CODON_CHANGE_PLUS_CODON_DELETION
                STOP_GAINED
                STOP_LOST
                UTR_3_DELETED

You can add extra classes using the --add option or specify a totally different set of classes using the --classes option. Alternatively you can use SnpEff "Effect_Impact" shortcuts with the --effect_impact  option to filter on these values.

Options are also available to additionally filter on dbNSFP fields as detailed above.  You may use the following shorter names for the given dbNSFP fields in your --filter_dbnsfp or --match_dbnsfp expressions:

                1000Gp1_AC: dbNSFP_1000Gp1_AC
                1000Gp1_AF: dbNSFP_1000Gp1_AF
                1000Gp1_AFR_AC: dbNSFP_1000Gp1_AFR_AC
                1000Gp1_AFR_AF: dbNSFP_1000Gp1_AFR_AF
                1000Gp1_AMR_AC: dbNSFP_1000Gp1_AMR_AC
                1000Gp1_AMR_AF: dbNSFP_1000Gp1_AMR_AF
                1000Gp1_ASN_AC: dbNSFP_1000Gp1_ASN_AC
                1000Gp1_ASN_AF: dbNSFP_1000Gp1_ASN_AF
                1000Gp1_EUR_AC: dbNSFP_1000Gp1_EUR_AC
                1000Gp1_EUR_AF: dbNSFP_1000Gp1_EUR_AF
                29way_logOdds: dbNSFP_29way_logOdds
                29way_pi: dbNSFP_29way_pi
                AA: dbNSFP_ESP6500_AA_AF
                AA_AF: dbNSFP_ESP6500_AA_AF
                AC: dbNSFP_1000Gp1_AC
                AF: dbNSFP_1000Gp1_AF
                AFR: dbNSFP_1000Gp1_AFR_AC
                AFR: dbNSFP_1000Gp1_AFR_AF
                AFR_AC: dbNSFP_1000Gp1_AFR_AC
                AFR_AF: dbNSFP_1000Gp1_AFR_AF
                AMR: dbNSFP_1000Gp1_AMR_AC
                AMR: dbNSFP_1000Gp1_AMR_AF
                AMR_AC: dbNSFP_1000Gp1_AMR_AC
                AMR_AF: dbNSFP_1000Gp1_AMR_AF
                ASN: dbNSFP_1000Gp1_ASN_AC
                ASN: dbNSFP_1000Gp1_ASN_AF
                ASN_AC: dbNSFP_1000Gp1_ASN_AC
                ASN_AF: dbNSFP_1000Gp1_ASN_AF
                EA: dbNSFP_ESP6500_EA_AF
                EA_AF: dbNSFP_ESP6500_EA_AF
                ESP6500_AA_AF: dbNSFP_ESP6500_AA_AF
                ESP6500_EA_AF: dbNSFP_ESP6500_EA_AF
                EUR: dbNSFP_1000Gp1_EUR_AC
                EUR: dbNSFP_1000Gp1_EUR_AF
                EUR_AC: dbNSFP_1000Gp1_EUR_AC
                EUR_AF: dbNSFP_1000Gp1_EUR_AF
                Ensembl_geneid: dbNSFP_Ensembl_geneid
                Ensembl_transcriptid: dbNSFP_Ensembl_transcriptid
                FATHMM: dbNSFP_FATHMM_pred
                FATHMM_pred: dbNSFP_FATHMM_pred
                FATHMM_score: dbNSFP_FATHMM_score
                FATHMM_score_converted: dbNSFP_FATHMM_score_converted
                GERP: dbNSFP_GERP++_RS
                GERP++: dbNSFP_GERP++_RS
                GERP++_NR: dbNSFP_GERP++_NR
                GERP++_RS: dbNSFP_GERP++_RS
                Interpro_domain: dbNSFP_Interpro_domain
                LRT: dbNSFP_LRT_pred
                LRT_Omega: dbNSFP_LRT_Omega
                LRT_pred: dbNSFP_LRT_pred
                LRT_score: dbNSFP_LRT_score
                LRT_score_converted: dbNSFP_LRT_score_converted
                MutationAssessor: dbNSFP_MutationAssessor_pred
                MutationAssessor_pred: dbNSFP_MutationAssessor_pred
                MutationAssessor_score: dbNSFP_MutationAssessor_score
                MutationAssessor_score_converted: dbNSFP_MutationAssessor_score_converted
                MutationTaster: dbNSFP_MutationTaster_pred
                MutationTaster_pred: dbNSFP_MutationTaster_pred
                MutationTaster_score: dbNSFP_MutationTaster_score
                MutationTaster_score_converted: dbNSFP_MutationTaster_score_converted
                Polyphen2: dbNSFP_Polyphen2_HVAR_pred
                Polyphen2_HDIV_pred: dbNSFP_Polyphen2_HDIV_pred
                Polyphen2_HDIV_score: dbNSFP_Polyphen2_HDIV_score
                Polyphen2_HVAR: dbNSFP_Polyphen2_HVAR_pred
                Polyphen2_HVAR_pred: dbNSFP_Polyphen2_HVAR_pred
                Polyphen2_HVAR_score: dbNSFP_Polyphen2_HVAR_score
                Polyphen2_pred: dbNSFP_Polyphen2_HVAR_pred
                Polyphen2_score: dbNSFP_Polyphen2_HVAR_score
                SIFT: dbNSFP_SIFT_pred
                SIFT_pred: dbNSFP_SIFT_pred
                SIFT_score: dbNSFP_SIFT_score
                SIFT_score_converted: dbNSFP_SIFT_score_converted
                UniSNP_ids: dbNSFP_UniSNP_ids
                Uniprot_aapos: dbNSFP_Uniprot_aapos
                Uniprot_acc: dbNSFP_Uniprot_acc
                Uniprot_id: dbNSFP_Uniprot_id
                aapos: dbNSFP_aapos
                cds_strand: dbNSFP_cds_strand
                codonpos: dbNSFP_codonpos
                gene: dbNSFP_Ensembl_geneid
                geneid: dbNSFP_Ensembl_geneid
                genename: dbNSFP_genename
                name: dbNSFP_genename
                phyloP: dbNSFP_phyloP
                refcodon: dbNSFP_refcodon
                strand: dbNSFP_cds_strand
                transcript: dbNSFP_Ensembl_transcriptid
                transcriptid: dbNSFP_Ensembl_transcriptid


=head1 AUTHOR

David A. Parry

University of Leeds

=head1 COPYRIGHT AND LICENSE

Copyright 2014  David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

