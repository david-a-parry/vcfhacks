#!/usr/bin/perl 

use warnings;
use strict;
use Config;
use Getopt::Long qw(:config no_ignore_case);
use POSIX qw/strftime/;
use Pod::Usage;
use Data::Dumper;
use Net::FTP;
use Cwd;
use IO::File;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use File::Path qw(make_path remove_tree);
use File::Temp qw/ tempfile tempdir /;
use File::Copy;
use File::Basename;
use File::stat;
use Cwd;
use Term::ProgressBar;

#use Bio::SeqIO::entrezgene_interactants;
use FindBin;
use lib "$FindBin::Bin/lib";
use ParseVCF;

my $genedir;
my $help;
my $man;
my $out;
my $vcf;
my $prep;
my $downdb;
my $repair;
my @classes;
my @add;
my $progress; 
my $functional;
my @annotations;
my $snpEff_mode;
GetOptions(
    "input=s"                   => \$vcf,
    "directory=s"               => \$genedir,
    "output=s"                  => \$out,
    "snpEff"                    => \$snpEff_mode,
    "functional"                => \$functional,
    "classes=s{,}"              => \@classes,
    "add=s{,}"                  => \@add,
    "gene_annotations=s{,}"     => \@annotations,
    "progress"                  => \$progress, 
    "PREPARE"                   => \$prep,
    "DOWNLOAD_NEW"              => \$downdb,
    "REPAIR"                    => \$repair,
    "help"                      => \$help,
    "manual"                    => \$man
) or pod2usage( -message => "Syntax error", -exitval => 2 );
pod2usage( -verbose => 2, ) if $man;
pod2usage( -verbose => 1, ) if $help;
pod2usage( -exitval => 2, -message => "--input is required", )
  if not $vcf
  and not $downdb
  and not $prep
  and not $repair;

my ( $script, $script_dir ) = fileparse($0);
$script_dir = getcwd() if $script_dir eq "./" or $script_dir eq ".\\";
$genedir = "$script_dir/gene_ref_files" if ( not $genedir );
$genedir =~ s/\/$//;
my $OUT;
if ($out) {
    open( $OUT, ">$out" ) or die "Can't open $out for writing: $!\n";
}
else {
    $OUT = \*STDOUT;
}
my %file_lengths = ();
my %database     = (
    ensemblToEntrez => {
        localfile => "$genedir/ensemblToEntrez",
        col       => 0,
        delimiter => "\t",
        url       => undef,
        dir       => undef,
        file      => "ensemblToEntrez" 
    },    #create this one on the fly from human_summary

#    gene2refseq => {
#        localfile => "$genedir/gene2refseq",
#        col       => 3,
#        delimiter => "\t",
#        url       => "ftp.ncbi.nlm.nih.gov",
#        dir       => "gene/DATA",
#        file      => "gene2refseq.gz"
#    },#this file is big and not used by anything except maybe by a custom made SnpEff database
#    #might be worth including at a later date but will need to implement some processing
#    to reduce the filesize (filter on taxonomy and throwaway entries without nucleotide IDs)

    orthologs => {
        localfile => "$genedir/HMD_Human5.rpt",
        col       => 1,
        delimiter => "\t",
        url       => "ftp.informatics.jax.org",
        dir       => "pub/reports/archive/orthology_reports.old/",
        file      => "HMD_Human5.rpt"
    },
    human_summary => {
        localfile => "$genedir/Homo_sapiens_ncbi_gene_all_summaries.txt",
        col       => 0,
        delimiter => "\t",
        url       => "ftp.ncbi.nlm.nih.gov",
        dir       => "gene/DATA/ASN_BINARY/Mammalia",
        file      => "Homo_sapiens.ags.gz"
    },
    mouse_phenotype_acc => {
        localfile => "$genedir/MGI_PhenotypicAllele.rpt",
        col       => 6,
        delimiter => "\t",
        url       => "ftp.informatics.jax.org",
        dir       => "pub/reports/",
        file      => "MGI_PhenotypicAllele.rpt"
    },
    mouse_phenotype_desc => {
        localfile => "$genedir/VOC_MammalianPhenotype.rpt",
        col       => 0,
        delimiter => "\t",
        url       => "ftp.informatics.jax.org",
        dir       => "pub/reports/",
        file      => "VOC_MammalianPhenotype.rpt"
    },
    mgi2entrez => {
        localfile => "$genedir/MGI_EntrezGene.rpt",
        col       => 8,
        delimiter => "\t",
        url       => "ftp.informatics.jax.org",
        dir       => "pub/reports/",
        file      => "MGI_EntrezGene.rpt"
    },
    mim_morbid => {
        localfile => "$genedir/morbidmap",
        col       => 2,
        delimiter => "\\|",
        url       => "ftp.omim.org",
        dir       => "OMIM",
        file      => "morbidmap"
    },
);

#NEED TO FIND ALTERNATIVE TO "HMD_Human5.rpt" - deprecated in current MGI scheme

my @ref_files = ();
foreach my $k ( keys %database ) {
    push @ref_files, $database{$k};
}
my @missing = ();
foreach (@ref_files) {
    push( @missing, $_ ) if not -e $_->{localfile};
}
if ($downdb) {
    prepare_database( \%database );
    exit if not $vcf;
}
elsif ($repair) {
    prepare_database( \%database, "repair" );
    exit if not $vcf;
}
else {
    my @missing_files = ();
    foreach my $m (@missing) {
        push @missing_files, $m->{localfile};
    }
    display_error_and_exit(
        "Missing database files - rerun with the --DOWNLOAD_NEW option to correct this.",
        "Can't find following files:\n" . join( "\n", @missing_files ) . "\n"
    ) if @missing;
}
if ($prep) {
    prepare_files( \@ref_files, \%database );
    exit;
}

@missing = ();

CHECK_INDEXES: foreach my $file (@ref_files) {
    if ( not -e $file->{localfile} . ".idx" ) {

        #print STDERR "Can't find index for $file\n";
        push( @missing, $file );
    }
}

if (@missing) {
    my @files_to_index = ();
    foreach (@missing) {
        push @files_to_index, $_->{localfile};
    }
    display_error_and_continue(
        "Missing Indexes",
        "Need to make indexes for files:\n"
          . join( "\n", @files_to_index ) . "\n"
    );
    prepare_files( \@missing, \%database );
}

my $next_update  = 0;
my $pre_progress = 0;
my @all_annotations = qw(ENSGENE_ID ENTREZ_ID SYMBOL GO_ID GO_DESCRIPTION GENERIFS SUMMARY OMIM MGI_PHENOTYPE);
if (@annotations){
    foreach my $anno (@annotations){
        if (not grep {/uc($_) eq uc($anno)/} @annotations){
            die "Unrecognised --gene_annotation \"$anno\" specified\n";
        }
    }
}else{
    @annotations = @all_annotations;
}
    
print STDERR "Initialising input VCF...\n";
my $vcf_obj = ParseVCF->new( file => $vcf );
if ($snpEff_mode){
    my $snp_eff_header = $vcf_obj->readSnpEffHeader();
    die "No 'Effect' field identified in header for file $vcf - " .
    "please annotate with SnpEff.\n" if (not exists $snp_eff_header->{Effect});
    checkAndParseSnpEffClasses(\@classes, \@add);
}else{
    my $vep_header = $vcf_obj->readVepHeader();
    die "No 'GENE' field identified in header for file $vcf - "
    . "please annotate with Ensembl's variant_effect_precictor.pl script.\n"
    if ( not exists $vep_header->{gene} );
    if ( @classes or @add ) {
        $functional = 1;
    }
    checkAndParseVepClasses();
}

my $pre_progressbar;
if ($progress){
    $pre_progressbar = Term::ProgressBar->new(
        {
            name  => "Parsing database",
            count => scalar(keys(%database)) ,
            ETA   => "linear",
        }
    );
}
foreach my $k ( keys %database ) {
    open( my $IN, $database{$k}->{localfile} )
      or die "Can't open reference file $database{$k}->{localfile}: $!\n";
    $database{$k}->{fh} = $IN;
}
foreach my $k ( keys %database ) {
    open( $database{$k}->{idx}, $database{$k}->{localfile} . ".idx" )
      or die
      "Can't open reference index file $database{$k}->{localfile}.idx: $!\n";
    binmode $database{$k}->{idx};
    $database{$k}->{length} =
      get_file_length_from_index( $database{$k}->{idx} );
    close $database{$k}->{idx};
    open( $database{$k}->{idx}, $database{$k}->{localfile} . ".idx" )
      or die
      "Can't open reference index file $database{$k}->{localfile}.idx: $!\n";
    binmode $database{$k}->{idx};
    $pre_progress++;
    if ($progress){
        $next_update = $pre_progressbar->update($pre_progress)
            if $pre_progress >= $next_update;
    }
}

my $progressbar; 
if ($progress){
    $progressbar = Term::ProgressBar->new(
    {
        name  => "Annotating",
        count => $vcf_obj->get_totalLines,
        ETA   => "linear",
    });
}
my $vcf_line = 0;
$next_update = 0;
print $OUT join( "", @{ $vcf_obj->get_metaHeader() } );
print $OUT
    "##INFO=<ID=GeneAnno,Number=.,Type=String,Description=\"Collected Entrez/MGI gene annotations from for VEP/SnpEff annotated human genes. ". 
    "Multiple values per annotation are separated using two colons ('::'), spaces are replaced with underscores, commas are replaced with the ` ".
    "symbol, and semi-colons are replaced with the ^ symbol so that regexes can be used to extract the original text programmatically. ".
    "Format: " . join("|", @annotations) ."\">\n";

my $header = $vcf_obj->getHeader(1);
print $OUT $header;
my $gene_field = "gene";#if using VEP annotations we get Gene ID from this key
if ($snpEff_mode){
    $gene_field = "Gene_Name";#if using SnpEff annotations we get Gene ID from this key
}
LINE: while ( my $line = $vcf_obj->readLine ) {
    $vcf_line++;
    my %annot = ();
    my @entrez_ids = ();
    my @csq = ();
    if ($snpEff_mode){
        @csq = $vcf_obj->getSnpEffFields([ "Effect", "Gene_Name", "Transcript_BioType" ] ); #returns array of hashes e.g. $csq[0]->{Gene_Name} = 'ABCA1' ; $csq[0]->{EFFECT} = 'DOWNSTREAM'
        die
"No consequence field found for line:\n$line\nPlease annotated your VCF file with SnpEff before running this program with the --snpEff_mode option.\n"
        if not @csq;
    }else{
        @csq = $vcf_obj->getVepFields( [ "Gene", "Consequence" ] );
        die
"No consequence field found for line:\n$line\nPlease annotated your VCF file with ensembl's variant effect precictor before running this program.\n"
        if not @csq;
    }
    foreach my $c (@csq) {
        if ($functional) {
            next if ( not check_consequence( $c ) );
        }
        if ($c->{$gene_field} =~ /^ENS/){
            push @entrez_ids, search_ensembl_id($c); 
        }elsif($c->{$gene_field} =~ /^\d+$/){
            push @entrez_ids, $c->{$gene_field} ; 
        }elsif($c->{$gene_field} =~ /^[NX][MR]_\d+$/){
            #not implemented unless we decide to include gene2refseq in database
           # push @entrez_ids, search_refseq_id($c); 
        }
    }

    #USE ENTREZ IDs TO SEARCH DATABSE FILES
    foreach my $entrez_id ( @entrez_ids ) {
        my $i = binSearchLineWithIndex(
            $entrez_id,
            $database{human_summary}->{fh},
            $database{human_summary}->{idx},
            $database{human_summary}->{length},
            $database{human_summary}->{col}
        );
        if ( $i > 0 ) {
            chomp(
                my @humsum_line = split(
                    "\t",
                    line_with_index(
                        $database{human_summary}->{fh},
                        $database{human_summary}->{idx},
                        $i
                    )
                )
            );
            $annot{$entrez_id}->{symbol}  = $humsum_line[1];
            $annot{$entrez_id}->{summary} = $humsum_line[2];
            push @{ $annot{$entrez_id}->{ensgene_id} }, split(/\|/, $humsum_line[3]);
            $annot{$entrez_id}->{go_id}   = $humsum_line[4];
            $annot{$entrez_id}->{go_description} = $humsum_line[5];
            $annot{$entrez_id}->{generifs}     = $humsum_line[6];

            #$annot{$entrez_id}->{interactants} = $humsum_line[8];
            my $mim_accession = $humsum_line[7];
            if ( $mim_accession ne "-" ) {
                my $j = binSearchLineWithIndex(
                    $mim_accession,
                    $database{mim_morbid}->{fh},
                    $database{mim_morbid}->{idx},
                    $database{mim_morbid}->{length},
                    $database{mim_morbid}->{col},
                    $database{mim_morbid}->{delimiter}
                );
                if ( $j > 0 ) {
                  MIM_UP:
                    for ( my $line_no = $j ; $line_no > 0 ; $line_no-- ) {
                        my @mim_line = split(
                            /\|/,
                            line_with_index(
                                $database{mim_morbid}->{fh},
                                $database{mim_morbid}->{idx},
                                $line_no
                            )
                        );
                        if ( $mim_line[2] eq $mim_accession ) {
                            push( @{ $annot{$entrez_id}->{omim} },
                                $mim_line[0] );
                        }
                        else {
                            last MIM_UP;
                        }
                    }
                  MIM_DOWN:
                    for (
                        my $line_no = $j + 1 ;
                        $line_no <= $database{mim_morbid}->{length} ;
                        $line_no++
                      )
                    {
                        my @mim_line = split(
                            /\|/,
                            line_with_index(
                                $database{mim_morbid}->{fh},
                                $database{mim_morbid}->{idx},
                                $line_no
                            )
                        );
                        if ( $mim_line[2] eq $mim_accession ) {
                            push( @{ $annot{$entrez_id}->{omim} },
                                $mim_line[0] );
                        }
                        else {
                            last MIM_DOWN;
                        }
                    }
                }
                else {
                    push( @{ $annot{$entrez_id}->{omim} }, "" );
                }
            }
            else {
                push( @{ $annot{$entrez_id}->{omim} }, "" );
            }
        }else {
            $annot{$entrez_id}->{symbol}  = "";
            $annot{$entrez_id}->{summary} = "";
            $annot{$entrez_id}->{go_id}   = "";
            $annot{$entrez_id}->{go_description} = "";
            $annot{$entrez_id}->{generifs}     = "";

            #$annot{$entrez_id}->{interactants} = "-";
            push( @{ $annot{$entrez_id}->{omim} }, "" );
        }
        my $ortholog = getMouseOrtholog(
            $entrez_id,
            $database{orthologs}->{fh},
            $database{orthologs}->{idx},
            $database{orthologs}->{length}
        );
        if ($ortholog) {
            push(
                my @mgi,
                get_MGI_phenotype(
                    $ortholog,
                    $database{mouse_phenotype_acc}->{fh},
                    $database{mouse_phenotype_acc}->{idx},
                    $database{mouse_phenotype_acc}->{length},
                    $database{mouse_phenotype_desc}->{fh},
                    $database{mouse_phenotype_desc}->{idx},
                    $database{mouse_phenotype_desc}->{length},
                    $database{mgi2entrez}->{fh},
                    $database{mgi2entrez}->{idx},
                    $database{mgi2entrez}->{length}
                )
            );
            @mgi
              ? push @{ $annot{$entrez_id}->{mgi_phenotype} }, @mgi
              : push @{ $annot{$entrez_id}->{mgi_phenotype} }, "";
        }
        else {
            push( @{ $annot{$entrez_id}->{mgi_phenotype} }, "" );
        }
    }
    if ( not keys %annot) {
        my $ens_annot = "|" x @annotations;
        my $new_line = $vcf_obj->addVariantInfoField('GeneAnno', $ens_annot);
        #print $OUT "-\t-\t-\t-\t-\t-\t-\t-\t$line\n";
        print $OUT "$new_line\n";
        next LINE;
    }
    my @gene_anno = ();
    foreach my $entrez_id ( keys %annot ) {
        my @single_anno = ();
        foreach my $annot (@annotations){
            if (lc$annot eq 'entrez_id'){
                push @single_anno, $entrez_id;    
            }elsif (exists $annot{$entrez_id}->{lc($annot)}){
                if (ref ($annot{$entrez_id}->{lc($annot)}) eq 'ARRAY'){
                    remove_duplicates( \@{ $annot{$entrez_id}->{lc($annot)} } );
                    my @conv = ();
                    foreach my $g (@{$annot{$entrez_id}->{lc($annot)}}){
                        push @conv, convert_text($g);
                    }
                    push @single_anno, join("::", @conv);
                }else{
                    my $converted = convert_text($annot{$entrez_id}->{lc($annot)});
                    push @single_anno, $converted;
                }
            }else{
                push @single_anno, '';
            }
        }
        push @gene_anno, join("|", @single_anno);
    }
    my $ens_annot = join(",", @gene_anno);
    my $new_line = $vcf_obj->addVariantInfoField('GeneAnno', $ens_annot);
    print $OUT "$new_line\n";
    if ($progress){
        $next_update = $progressbar->update($vcf_line) if $vcf_line >= $next_update;
    }
}
if ($progress){
    $progressbar->update( $vcf_obj->get_totalLines )
        if $vcf_obj->get_totalLines >= $next_update;
}
for my $k ( keys %database ) {
    close $database{$k}->{fh};
    close $database{$k}->{idx};
}
print STDERR "\nAnnotation finished.\n";
$vcf_obj -> DESTROY();

########################################
sub search_refseq_id{
    my ($csq) = @_;
    my @ids = ();
    (my $nm_short = $csq->{$gene_field}) =~ s/\.\d+$//;
    my $i = binSearchRefseqWithIndex(
        $nm_short,
        $database{gene2refseq}->{fh},
        $database{gene2refseq}->{idx},
        $database{gene2refseq}->{length},
        $database{gene2refseq}->{col}
    );
    return if ( $i < 1 );
    for ( my $j = $i - 1 ; $j > 0 ; $j-- ) {
        my @ens_line = split(
            "\t",
            line_with_index(
                $database{gene2refseq}->{fh},
                $database{gene2refseq}->{idx},
                $j
            )
        );
        chomp @ens_line;
        (my $ref_short = $ens_line[ $database{gene2refseq}->{col} ]) =~ s/\.\d\+$//;
        if ( $ref_short eq $nm_short )
        {
            push @ids, $ens_line[1] ; 
        }
        else {
            last;
        }
    }
    my @ens_line = split(
        "\t",
        line_with_index(
            $database{ensemblToEntrez}->{fh},
            $database{ensemblToEntrez}->{idx},
            $i
        )
    );
    chomp @ens_line;
    push @ids, $ens_line[1] ; 
    for (
        my $j = $i + 1 ;
        $j <= $database{ensemblToEntrez}->{length} ;
        $j++
        )
    {
        my @ens_line = split(
            "\t",
            line_with_index(
                $database{ensemblToEntrez}->{fh},
                $database{ensemblToEntrez}->{idx},
                $j
            )
        );
        chomp @ens_line;
        (my $ref_short = $ens_line[ $database{gene2refseq}->{col} ]) =~ s/\.\d\+$//;
        if ( $ref_short eq $nm_short )
        {
            push @ids, $ens_line[1] ; 
        }
        else {
            last;
        }
    }
    return @ids;
}


########################################
sub search_ensembl_id{
    my ($csq) = @_;
    my @ids = (); 
    my $i = binSearchLineWithIndex(
        $csq->{$gene_field},
        $database{ensemblToEntrez}->{fh},
        $database{ensemblToEntrez}->{idx},
        $database{ensemblToEntrez}->{length},
        $database{ensemblToEntrez}->{col}
    );
    return if ( $i < 1 );
    for ( my $j = $i - 1 ; $j > 0 ; $j-- ) {
        my @ens_line = split(
            "\t",
            line_with_index(
                $database{ensemblToEntrez}->{fh},
                $database{ensemblToEntrez}->{idx},
                $j
            )
        );
        chomp @ens_line;
        if ( $ens_line[ $database{ensemblToEntrez}->{col} ] eq $csq->{$gene_field} )
        {
            push @ids, $ens_line[1];
        }
        else {
            last;
        }
    }
    my @ens_line = split(
        "\t",
        line_with_index(
            $database{ensemblToEntrez}->{fh},
            $database{ensemblToEntrez}->{idx},
            $i
        )
    );
    chomp @ens_line;
    push @ids, $ens_line[1];
    for (
        my $j = $i + 1 ;
        $j <= $database{ensemblToEntrez}->{length} ;
        $j++
        )
    {
        my @ens_line = split(
            "\t",
            line_with_index(
                $database{ensemblToEntrez}->{fh},
                $database{ensemblToEntrez}->{idx},
                $j
            )
        );
        chomp @ens_line;
        if ( $ens_line[ $database{ensemblToEntrez}->{col} ] eq $csq->{$gene_field} )
        {
            push @ids, $ens_line[1];
        }
        else {
            last;
        }
    }
    return @ids;
}

########################################
sub checkAndParseSnpEffClasses{

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

    if (not @classes){
        @classes = qw (
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
    push (@classes, @add) if (@add);
    foreach my $class (@classes){
        die "Error - variant class '$class' not recognised.\n" if not grep {/^$class$/i} @valid;
    }
}
########################################
sub checkAndParseVepClasses{
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
        @classes = qw (transcript_ablation
        splice_donor_variant
        splice_acceptor_variant
        splice_region_variant
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
    foreach my $class (@classes) {
        die "Error - variant class '$class' not recognised.\n"
        if not grep { /$class/i } @valid;
    }
}
########################################
sub convert_text{
    #replace characters not compatible with VCF
    #format for annotations
    my ($string) = @_;
    $string =~ s/\|/::/g;
    $string =~ tr/;, /^`_/;
    return $string;
}
########################################
sub check_consequence {
    my ( $annot ) = @_;
    if ($snpEff_mode){
        return check_snpeff_consequence($annot);
    }else{
        return check_vep_consequence($annot);
    }
}
########################################
sub check_snpeff_consequence {
    my ( $annot ) = @_;
    if ($annot->{Transcript_BioType} =~ /nonsense_mediated_decay|pseudogene$/i){
        return 0;
    }
CLASS: foreach my $class (@classes) {
        if (lc$annot->{Effect} eq lc$class){
            return 1;
        }
    }
    return 0;
}
########################################
sub check_vep_consequence {

    #return 1 if matches class and isn't a NMD_transcript_variant
    my ( $annot ) = @_;
    my @anno_csq = split( /\&/, $annot->{consequence} );
    if ( grep { /NMD_transcript_variant/i } @anno_csq ) {
        return 0;
    }
CLASS: foreach my $class (@classes) {
        foreach my $ac (@anno_csq) {
            if ( lc $ac eq lc $class ) {
                return 1;
            }
        }
    }
    #no match
    return 0;
}
########################################
sub prepare_files {
    my ( $file_array_ref, $database_ref ) = @_;
    my $increment    = 100 / @$file_array_ref;
    my $prep_percent = 0;
    my $next_update  = 0;
    my $prep_bar     = Term::ProgressBar->new(
        { name => "Preparing Database", count => 100, ETA => "linear", } );
    foreach my $file (@$file_array_ref) {
        sort_and_index_gene_files(
            $file->{localfile}, $file->{col}, \$prep_percent,
            $increment,         \$prep_bar,   \$next_update,
            $file->{delimiter}
        );
    }
}

########################################
sub prepare_database {
    my ( $database_ref, $mode ) = @_;
    $mode = "replace" if not $mode;
    my $dir = getcwd();
    my @files;
    my $oldest = 0;
    if ( $mode eq "replace" ) {
        foreach my $k ( keys %{$database_ref} ) {
            push @files, $database_ref->{$k};
        }
    }
    elsif ( $mode eq "repair" ) {
        foreach my $k ( keys %$database_ref ) {
            if ( not -e $database_ref->{$k}->{localfile} ) {
                push @files, $database_ref->{$k};
            }
            else {
                my $age = -M $database_ref->{$k}->{localfile};
                $oldest < $age ? $oldest = $age : ();
            }
        }
        $oldest = int $oldest;
        if ( @files < 1 ) {
            if ($oldest) {
                display_error_and_continue(
                    "No files missing in database.",
"Oldest file is $oldest days old. Use the \"Update Database\" option to update."
                );
            }
            else {
                display_error_and_continue( "Database OK",
                    "No files missing in database." );
            }
            return 0;
        }
    }
    else {
        die "Unrecognised mode ($mode) for prepare_database subroutine\n";
    }
    my $db_percent  = 0;
    my $next_update = 0;
    #my $progressbar = Term::ProgressBar->new(
    #    { name => "Prep Database", count => 100, ETA => "linear", } );
    if ( not -e $genedir ) {
        unless (mkdir $genedir){
            display_error_and_exit(
                    "Permissions Error",
"Can't create directory $genedir for database files - please check permissions"
            );
        }
    }
    foreach my $file (@files) {
        if (defined $file->{file} && $file->{file} eq "Homo_sapiens.ags.gz"){
            eval "use Bio::SeqIO::entrezgene; 1" 
              or die "The Bio::SeqIO::entrezgene module must be installed in order ".
              "to extract NCBI gene summaries for database creation. ".
              "Please install bioperl and try again\n";
        }
    }
    my $t = @files;
    if ( grep { $_->{file} eq "ensemblToEntrez" } @files ){
        $t--;
        if ( !grep { $_->{file} eq 'Homo_sapiens.ags.gz' } @files ){
            $t++;
            push @files, $database_ref->{human_summary};
        }
    }
    my $n = 0;
    foreach my $file (@files) {
        my ( $file_name, $file_dir ) = fileparse( $file->{localfile} );
        if ( $file_name eq 'ensemblToEntrez' )
        {    #we create this one when processing the human summaries
            next;
        }
        $n++;
        my $time = strftime( "%H:%M:%S", localtime );
        print STDERR "[$time] Processing $file_name, file $n of $t...\n";
        my $increment = 100 / @files;
        $increment /= 10;
        chdir $genedir or display_error_and_exit(
            "Directory Error",
            "Can't move to directory $genedir",
        );
        my $file_exists = 0;

        if ( -e $file ) {
            move( $file, "$file.bkup" ) or display_error_and_exit(
                "File Error",
                "Error creating file backup of $file->{localfile}",
                "Check permissions and/or disk space."
            );
            $file_exists++;
        }
        $time = strftime( "%H:%M:%S", localtime );
        print STDERR "[$time] Connecting to $file->{url}...\n";
        my $ftpobj =
          Net::FTP->new( $file->{url} )
          or restore_file( $file_exists, $file )
          && display_error_and_exit( "Can't connect to site $file->{url}",
            "Could not download $file->{url}/$file->{dir}/$file->{file}" );
        $ftpobj->login( "anonymous", "" )
          or restore_file( $file_exists, $file )
          && display_error_and_exit( "Can't login to $file->{url}",
            "Could not download $file->{url}/$file->{dir}/$file->{file}" );
        $ftpobj->cwd( $file->{dir} )
          or restore_file( $file_exists, $file ) && display_error_and_exit(
            "Can't locate directory $file->{dir} at $file->{url}",
            "Could not download $file->{url}/$file->{dir}/$file->{file}"
          );
        $ftpobj->binary();
        $time = strftime( "%H:%M:%S", localtime );
        print STDERR "[$time] Downloading $file->{file}...\n";
        $ftpobj->get( $file->{file} )
          or restore_file( $file_exists, $file )
          && display_error_and_exit( "Download error",
            "Could not download $file->{url}/$file->{dir}/$file->{file}" );
        $db_percent += $increment;
        #$next_update = $progressbar->update($db_percent)
          #if $db_percent >= $next_update;
        $increment *= 9;

        if ( $file->{file} =~ /\.gz$/ ) {
            $time = strftime( "%H:%M:%S", localtime );
            print STDERR "[$time] Decompressing $file->{file}...\n";
            $increment /= 3;
            ( my $output = $file->{file} ) =~ s/\.gz$//i;
            if ( $file->{file} =~ /\.ags/ )
            {    #make sure we use binmode for .ags files
                 #for some reason gzipped ASN1 file seems to give premature EOF even in binmode, so use MultiStream option
                my $z = new IO::Uncompress::Gunzip $file->{file},
                  MultiStream => 1
                  or restore_file( $file_exists, $file )
                  && display_error_and_exit( "Error decompressing file",
                    "Error decompressing $file->{file}" );
                open( my $ZOUT, ">$output" )
                  or restore_file( $file_exists, $file )
                  && display_error_and_exit( "Error decompressing file",
                    "Error decompressing $file->{file}" );
                binmode $ZOUT;
                $z->binmode;
                my $buffer;
                while ( $z->read($buffer) ) {
                    print $ZOUT $buffer;
                }
                $db_percent += $increment;
                #$next_update = $progressbar->update($db_percent)
                  #if $db_percent >= $next_update;

            }
            else {
                gunzip( $file->{file} => $output )
                  or restore_file( $file_exists, $file )
                  && display_error_and_exit( "Error decompressing file",
                    "Error decompressing $file->{file}" );
                $db_percent += $increment;
                #$next_update = $progressbar->update($db_percent)
                  #if $db_percent >= $next_update;
            }

            #$file_name = $output unless $file->{file} =~ /\.ags/;
            #$next_update = $progressbar->update($db_percent)
              #if $db_percent >= $next_update;
            $increment *= 2;
        }
        if ( $file->{file} =~ /\.ags/ ) {
            $increment /= 2;

            #use gene2xml script to extract summaries...
            my $gene2xml = "./gene2xml";
            if ( not -e $gene2xml ) {
                $time = strftime( "%H:%M:%S", localtime );
                print STDERR "[$time] Retrieving gene2xml executable...\n";
                download_gene2xml("./");
            }
            ( my $decomp_file = $file->{file} ) =~ s/\.gz$//;

#my $command = "\"$gene2xml\" -i \"$decomp_file\" -b -o \"$xml_out\"";
#my $exit_status = `$command`;
#display_error_and_continue("Error processing gene2xml command", "Exit status of $exit_status from $command") if $exit_status;
            $db_percent += $increment / 2;
            #$next_update = $progressbar->update($db_percent)
              #if $db_percent >= $next_update;

#unlink $file->{file} or display_error_and_exit( "Can't delete xml output ($file->{file})", "Check permissions - it is safe to manually delete this file now");
            my ( $enstoEntrez_file_name, $file_dir ) = fileparse( $database_ref->{ensemblToEntrez}->{localfile} );
            extract_ncbi_summaries( $gene2xml, $decomp_file, "$file_name.tmp",
                $enstoEntrez_file_name .".tmp" );
                move( "$file_name.tmp", $file_name ) or display_error_and_exit(
                "File Error",
                "Error creating file $file_name",
                "Check permissions and/or disk space."
                );
                move( $enstoEntrez_file_name .".tmp", $enstoEntrez_file_name) or display_error_and_exit(
                "File Error",
                "Error creating file $database_ref->{ensemblToEntrez}->{localfile}",
                "Check permissions and/or disk space."
                );
            $time = strftime( "%H:%M:%S", localtime );
            print STDERR "[$time] Sorting and indexing ensemblToEntrez file...\n";
            sort_and_index_gene_files(
                $enstoEntrez_file_name,
                $database_ref->{ensemblToEntrez}->{col},
                \$db_percent,
                $increment,
                $progressbar,
                $next_update,
                $database_ref->{ensemblToEntrez}->{delimiter}
            );
            $db_percent += $increment / 2;
            #$next_update = $progressbar->update($db_percent)
              #if $db_percent >= $next_update;
            unlink $decomp_file or display_error_and_continue(
                "Can't delete decompressed ags file ($decomp_file)",
"Check permissions - it is safe to manually delete this file now"
            );

#unlink $xml_out or display_error_and_exit( "Can't delete xml output ($xml_out)", "Check permissions - it is safe to manually delete this file now");
        }
        if ( $file->{file} =~ /HMD_Human5\.rpt/ ) {
            $increment /= 3;
            move( $file->{file}, "$file->{file}.bak" )
              or display_error_and_exit(
                "File Error",
                "Error creating file backup of $file->{file}",
                "Check permissions and/or disk space."
              );
            open( my $HMD_DOWN, "$file->{file}.bak" )
              or display_error_and_exit( "File Read Error",
                "Can't open $file->{file}.bak to read" );
            open( my $HMD_MOD, ">$file->{file}" )
              or display_error_and_exit( "Write Error",
                "Can't open $file->{file} to write." );
            while ( my $line = <$HMD_DOWN> ) {
                next if $line !~ /^[0-9MXYU]/;
                print $HMD_MOD $line;
            }
            $db_percent += $increment;
            #$next_update = $progressbar->update($db_percent)
              #if $db_percent >= $next_update;
            close $HMD_MOD;
            close $HMD_DOWN;
            unlink "$file->{file}.bak" or display_error_and_continue(
                "Can't delete backup file $file->{file}.bak",
"Check permissions - it is safe to manually delete this file now"
            );
            $increment *= 2;
        }
        chdir $dir;
        $time = strftime( "%H:%M:%S", localtime );
        print STDERR "[$time] Sorting and indexing $file_name...\n";
        sort_and_index_gene_files( "$file_dir/$file_name", $file->{col},
            \$db_percent, $increment, $progressbar, $next_update,
            $file->{delimiter} );
        if (-e "$file.bkup"){
            unlink "$file.bkup"
                or display_error_and_continue(
                "Can't delete backup file \"$file.bkup\"",
                "Check permissions - it is safe to manually delete this file now" );
        }
    }
    #$progressbar->update(100);
    my $time = strftime( "%H:%M:%S", localtime );
    print STDERR "[$time] Database update finished.\n";
}

#########################################
sub sort_and_index_gene_files {
    my (
        $file,        $sort_column, $prog_ref, $increment,
        $progressbar, $next_update, $delimiter
    ) = @_;    #column no. is 0 based
    open( my $FILE, $file )
      or die "Can't open $file for sorting and indexing!\n";
    my ( $file_short, $dir ) = fileparse($file);
    my @lines;
    $delimiter = "\t" if not $delimiter;
    my $line_count = 0;
    $line_count += tr/\n/\n/ while sysread( $FILE, $_, 2**16 );
    $line_count ||= 1;

    #print STDERR "$file has $line_count lines\n";/
    close $FILE;
    my $incr_per_line = ( $increment / 6 ) / $line_count;

    #print STDERR "increment per line is $incr_per_line\n";
    open( $FILE, $file ) or die "Can't open $file for sorting and indexing!\n";
    my $check_taxon = 0;
    $check_taxon = 1
      if ( $file =~ /gene2go/ or $file =~ /generifs/ )
      ;    #for these files we'll ignore non-mouse and non-human genes
    my $line_counter = 0;
    my $prev_counter = 0;
    while (<$FILE>) {
        $$prog_ref += $incr_per_line;
        #$next_update = $progressbar->update($$prog_ref)
          #if $$prog_ref >= $next_update;
        next if /^#/;
        if ($check_taxon) {
            my $tax_id = ( split "\t" )[0];
            next if ( $tax_id != 9606 and $tax_id != 10090 );
        }
        chomp;
        if ( ( split /$delimiter/ )[$sort_column] )
        {    #some MGI files have empty columns
            push( @lines, [ $_, ( split /$delimiter/ )[$sort_column] ] );
        }
    }

    #$$prog_ref += $increment/6;
    #print $bar "$$prog_ref";
    close $FILE;
    $incr_per_line = ( $increment / 2 ) / @lines;
    @lines = sort { $a->[1] cmp $b->[1] } @lines;
    $$prog_ref += $increment / 2;
    #$next_update = $progressbar->update($$prog_ref)
      #if $$prog_ref >= $next_update;
    my ( $tmp, $TEMP );
    ( $TEMP, $tmp ) = tempfile( "$dir/tmp_dharmaXXXX", UNLINK => 1 )
      or die "Can't create temporary sort file\n";
    $incr_per_line = ( $increment / 6 ) / @lines;
    $line_counter  = 0;
    $prev_counter  = 0;

    foreach my $line (@lines) {
        $$prog_ref += $incr_per_line;
        #$next_update = $progressbar->update($$prog_ref)
          #if $$prog_ref >= $next_update;
        $prev_counter = int($line_counter);
        print $TEMP "$line->[0]\n";
    }
    close $TEMP;
    move( $tmp, $file );

    #print STDERR "$file replaced with sorted version.\n";
    my $indexfile = $file . ".idx";
    open( my $INDEX, "+>$indexfile" )
      or die "can't open $indexfile for writing index ";
    binmode $INDEX;
    open( my $NEWFILE, "$file" ) or die "Can't open $file for reading ";
    build_index( $NEWFILE, $INDEX );
    $$prog_ref += $increment / 6;
    #$next_update = $progressbar->update($$prog_ref)
      #if $$prog_ref >= $next_update;

}

#####################################
sub build_index {

    # usage: build_index(*DATA_HANDLE, *INDEX_HANDLE)
    my $data_file  = shift;
    my $index_file = shift;
    my $offset     = 0;
    while (<$data_file>) {
        print $index_file pack( "N", $offset );
        $offset = tell($data_file);
    }
}

#####################################
sub line_with_index {

    # usage: line_with_index(*DATA_HANDLE, *INDEX_HANDLE, $LINE_NUMBER)
    # returns line or undef if LINE_NUMBER was out of range
    my $data_file   = shift;
    my $index_file  = shift;
    my $line_number = shift;
    my $size;        # size of an index entry
    my $i_offset;    # offset into the index of the entry
    my $entry;       # index entry
    my $d_offset;    # offset into the data file
    $size = length( pack( "N", 0 ) );
    $i_offset = $size * ( $line_number - 1 );
    seek( $index_file, $i_offset, 0 ) or return;
    read( $index_file, $entry, $size );
    $d_offset = unpack( "N", $entry );
    seek( $data_file, $d_offset, 0 );
    return scalar(<$data_file>);
}
######################################################################
sub binSearchLineWithIndex {
    my ( $x, $FILE, $FILE_INDEX, $u, $column, $delimiter ) = @_;
    my $l = 1;
    $delimiter = "\t" if not $delimiter;
    while ( $l <= $u ) {
        my $i = int( ( $u + $l ) / 2 );
        my $line = line_with_index( $FILE, $FILE_INDEX, $i );
        if ( $line =~ /^#/ ) {
            $l = $i + 1;
        }
        else {
            my @line = split( /$delimiter/, $line );
            if ( $x lt $line[$column] ) {
                $u = $i - 1;
            }
            elsif ( $x gt $line[$column] ) {
                $l = $i + 1;
            }
            else {
                return $i;
            }
        }
    }
    return 0;
}
######################################################################
sub binSearchRefseqWithIndex {
    my ( $x, $FILE, $FILE_INDEX, $u, $column, $delimiter ) = @_;
    my $l = 1;
    $delimiter = "\t" if not $delimiter;
    while ( $l <= $u ) {
        my $i = int( ( $u + $l ) / 2 );
        my $line = line_with_index( $FILE, $FILE_INDEX, $i );
        if ( $line =~ /^#/ ) {
            $l = $i + 1;
        }
        else {
            my @line = split( /$delimiter/, $line );
			(my $refshort = $line[$column]) =~ s/\.\d+$//;
            if ( $x lt $refshort ) {
                $u = $i - 1;
            }
            elsif ( $x gt $refshort ) {
                $l = $i + 1;
            }
            else {
                return $i;
            }
        }
    }
    return 0;
}

######################################################################
sub extract_ncbi_summaries {

    #produces file of gene ID plus Summary plus Interactants
    my ( $gene2xml, $agsfile, $sum_out, $ensToEntrezout ) = @_;
    eval "use Bio::SeqIO::entrezgene; 1" 
        or die "The Bio::SeqIO::entrezgene module must be installed in order ".
        "to extract NCBI gene summaries. Please install bioperl and try again\n";
    my $time = strftime( "%H:%M:%S", localtime );
    my $cmd = "$gene2xml -i $agsfile -b -x > $agsfile.asn.1";
    print STDERR "[$time] Converting NCBI binary gene data to text from $agsfile using the following command:\n$cmd\n";
    system($cmd); 
    if ($? == -1) {
        die "failed to execute conversion: $!\n";
    }elsif ($?){
        die "Conversion command failed: $?\n";
    }
    my $io = Bio::SeqIO->new(
        -format => 'entrezgene',
        -file   => "$agsfile.asn.1"
    );
    open( my $ENSOUT, ">$ensToEntrezout" )
      or die "Can't open $ensToEntrezout for writing\n";

    open( my $SUMOUT, ">$sum_out" ) or die "Can't open $sum_out for writing\n";

    $time = strftime( "%H:%M:%S", localtime );
    print STDERR "[$time] Extracting NCBI gene data from $agsfile - this may take some time...\n";
    while ( my ( $gene, $genestructure, $uncaptured ) = $io->next_seq ) {
    
        my $chrom;
        my @go_ids      = ();
        my @go_descs    = ();
        my @rifs        = ();
        my @mim         = ();
        my @ensembl     = ();
        my $annotations = $gene->annotation;
        my @annotations = $annotations->get_Annotations();
        
        print $SUMOUT "$gene->{primary_seq}->{accession_number}\t";
        exists $gene->{primary_seq}->{display_id}
          ? print $SUMOUT "$gene->{primary_seq}->{display_id}\t"
          : print $SUMOUT "-\t";
        exists $gene->{primary_seq}->{desc}
          ? print $SUMOUT "$gene->{primary_seq}->{desc}\t"
          : print $SUMOUT "-\t";
        
        foreach my $anno (@annotations) {
            if ( $anno->tagname eq 'dblink' ) {
                if (   ( exists $anno->{database} )
                    && ( $anno->{database} eq 'generif' ) )
                {
                    push( @rifs, $anno->{comment}->{text} )
                      if defined $anno->{comment}->{text};
                }
            }
            elsif ( $anno->tagname eq 'MIM' ) {
                push( @mim, $anno->{value} );
            }
            elsif ( $anno->tagname eq 'chromosome' ) {
                $chrom = $anno->{value};
            }
            elsif ( $anno->tagname eq 'Ensembl' ) {
                push @ensembl, $anno->{value};
                print $ENSOUT
                  "$anno->{value}\t$gene->{primary_seq}->{accession_number}\n";
            }
        }
        @annotations = $annotations->get_Annotations('OntologyTerm');
        foreach my $anno (@annotations) {
            next
              if ( ( defined $anno->term->authority )
                && ( $anno->term->authority eq 'STS marker' ) )
              ;    # No STS markers
            push( @go_ids,   $anno->identifier );
            push( @go_descs, $anno->name );
        }
         
        #my (@seq_member) = $genestructure->get_members;
        @ensembl
          ? print $SUMOUT join( "|", @ensembl ) . "\t"
          : print $SUMOUT "-\t";
        @go_ids
          ? print $SUMOUT join( "|", @go_ids ) . "\t"
          : print $SUMOUT "-\t";
        @go_descs
          ? print $SUMOUT join( "|", @go_descs ) . "\t"
          : print $SUMOUT "-\t";
        @rifs ? print $SUMOUT join( "|", @rifs ) . "\t" : print $SUMOUT "-\t";
        @mim  ? print $SUMOUT join( "|", @mim ) . "\t"  : print $SUMOUT "-\t";
   
# my @interactants;
# foreach my $k (keys %{$inters}){
#         push (@interactants, "$inters->{$k}($k)");
# }
# @interactants ? print $SUMOUT join("|", @interactants) ."\t"  : print $SUMOUT "-\t";
        print $SUMOUT "\n";
    }
    
    close $SUMOUT;
    close $ENSOUT;
    unlink "$agsfile.asn.1" or warn "Error deleting $agsfile.asn.1 file - you may want to delete this manually.\n";
}

#####################################
sub display_error_and_exit {
    my ( $text, $informative_text, $button ) = @_;

#we require the main text at a minimum.  Follow this with informative text and $button label if desired.
    $informative_text = "" if not $informative_text;
    print STDERR ">$text\n>$informative_text\n";
    exit 2;
}
#####################################
sub display_error_and_continue {
    my ( $text, $informative_text ) = @_;

#we require the main text at a minimum.  Follow this with informative text if desired.
    $informative_text = "" if not $informative_text;
    print STDERR ">$text\n>$informative_text\n";
}

#########################################
sub getMouseOrtholog {
    my ( $gene, $fh, $index, $u ) = @_;
    my $i = binSearchLineWithIndex( $gene, $fh, $index, $u, 1 );
    if ( $i > 0 ) {
        return ( split "\t", line_with_index( $fh, $index, $i ) )[6];
    }
}

#########################################
sub get_MGI_phenotype {
    my (
        $gene_id,               $PHENO,
        $PHENO_INDEX,           $pheno_line_count,
        $PHENO_DESC,            $PHENO_DESC_INDEX,
        $pheno_desc_line_count, $MGI_ENTREZ,
        $MGI_ENTREZ_INDEX,      $mgi_entrez_line_count
    ) = @_;
    my @desc = ();
    my $i    = binSearchLineWithIndex( $gene_id, $MGI_ENTREZ, $MGI_ENTREZ_INDEX,
        $mgi_entrez_line_count, 8 );
    if ( $i > 0 ) {
        my $accession =
          ( split "\t", line_with_index( $MGI_ENTREZ, $MGI_ENTREZ_INDEX, $i ) )
          [0];
        my $j = binSearchLineWithIndex( $accession, $PHENO, $PHENO_INDEX,
            $pheno_line_count, 6 );
        if ( $j > 0 ) {
            my @phenos = split( ",",
                ( split "\t", line_with_index( $PHENO, $PHENO_INDEX, $j ) )[10]
            );
            if (@phenos) {
                foreach my $ph (@phenos) {
                    my $k = binSearchLineWithIndex( $ph, $PHENO_DESC,
                        $PHENO_DESC_INDEX, $pheno_desc_line_count, 0 );
                    if ( $k > 0 ) {
                        push(
                            @desc,
                            (
                                split "\t",
                                line_with_index(
                                    $PHENO_DESC, $PHENO_DESC_INDEX, $k
                                )
                            )[1]
                        );
                    }
                    else {
                        display_error_and_continue(
                            "Missing phenotype description",
"Can't find matching phenotype description for phenotype ID $ph\nTry updating your database."
                        );
                    }
                }
            }
        }
    }
    return @desc;
}

######################################################################
sub restore_file {
    my ( $exists, $file ) = @_;
    if ($exists) {
        move( "$file.bkup", $file );
    }
}

######################################################################
sub get_file_length_from_index {
    my $index_file  = shift;
    my $size        = length( pack( "N", 0 ) );
    my $index_size  = ( -s $index_file );
    my $file_length = $index_size / $size;
    return $file_length;
}
######################################################################
sub remove_duplicates {

    #remove dups from arrays passed as array references
    foreach my $array_ref (@_) {
        die "Value passed to remove_duplicates not an array reference!\n"
          if ref $array_ref ne 'ARRAY';
        my %seen = ();
        @$array_ref = grep { !$seen{$_}++ } @$array_ref;
    }
}

######################################################################
sub download_gene2xml {
    my ($destination) = @_;
    my $pwd = getcwd();
    chdir $destination;
    my $site = "ftp.ncbi.nlm.nih.gov";
    my $dir  = "toolbox/ncbi_tools/converters/by_program/gene2xml";
    my $ftpobj =
      Net::FTP->new($site) or display_error_and_exit( "Can't connect to $site ",
        "Could not download gene2xml" );
    $ftpobj->login( "anonymous", "" )
      or display_error_and_exit( "Login to $site failed",
        "Login failed when attempting to download gene2xml" );
    $ftpobj->cwd($dir) or display_error_and_exit(
        "Couldn't change directory to $dir at $site",
        "cwd failed when attempting to download gene2xml"
    );
    $ftpobj->binary();
    my $prog = "";

    if ( $^O eq 'darwin' ) {
        $prog = "mac.gene2xml.gz";
    }
    elsif ( $^O eq 'linux' ) {
        if ( $Config{archname} =~ 'x86_64' ) {
            $prog = "linux64.gene2xml.gz";
        }
        else {
            $prog = "linux.gene2xml.gz";
        }
    }
    else {
        display_error_and_exit(
            "Unsupported OS: $^O",
"Failed to retrieve gene2xml because your OS is not supported.  If you have a binary place it in your database folder and rerun the database update"
        );
    }
    $ftpobj->get($prog)
      or display_error_and_exit( "Failed to retrieve $prog from $site", "$!" );
    my $output = "gene2xml";
    gunzip( $prog => $output )
      or display_error_and_exit( "Extract of $prog failed", "$!" );
    chmod 0755, $output
      or display_error_and_exit( "Could not make $output executable", "$!" );
    chdir $pwd;
}


=head1 NAME

ensemblGeneAnnotator.pl - add gene annotations to a VCF file annotated by Ensembl's variant_effect_predictor.pl

=head1 SYNOPSIS

        ensemblGeneAnnotator.pl -i [VEP annotated vcf file] -d [directory containing reference files] [options]
        ensemblGeneAnnotator.pl -h (display help message)
        ensemblGeneAnnotator.pl -m (display manual page)

=cut

=head1 ARGUMENTS

=over 8

=item B<-i    --input>

Input VCF file annotated with variant_effect_predictor.pl or snpEff. Required.

=item B<-o    --output>

Output file name.

=item B<-d    --directory>

Directory containing reference files. Will look for reference files in the same directory as this program if not supplied.

=item B<-s    --snpEff>

Use SnpEff annotations instead of VEP annotations to identify genes. You must have used the -geneId flag when annotating with snpEff for this script to work.

=item B<-f    --functional>

Use this flag to only annotate standard 'functional' variant classes (transcript_ablation, splice_donor_variant, splice_acceptor_variant, splice_region_variant, stop_gained, frameshift_variant, stop_lost, initiator_codon_variant, inframe_insertion, inframe_deletion, missense_variant, transcript_amplification, TFBS_ablation, TFBS_amplification, regulatory_region_ablation, regulatory_region_amplification). Prevents annotation of gene information for genes/transcripts that overlap a variant but are not affected in a way defined by one of these variant classes.

=item B<-c    --classes>

Use this to specify VEP/SnpEff variant classes to annotate. Gene information will only be annotated for genes affected by one of these variant classes. Overrides --functional option.

=item B<-a    --additional_classes>

Use this to specify additional VEP variant classes to annotate as well as those used by --functional.

=item B<-g    --gene_annotations>

List of gene annotations to include in output. By default all of the following classes are included:

    ENSGENE_ID
    ENTREZ_ID
    SYMBOL
    GO_ID
    GO_DESCRIPTION
    GENERIFS
    SUMMARY
    OMIM
    MGI_PHENOTYPE

Specify one or more of these to limit the annotations in your output to these classes only.

=item B<-p    --progress>

Show a progress bar.

=item B<-P    --PREPARE>

Sort and index reference files (e.g. if you've downloaded them manually).

=item B<-D    --DOWNLOAD_NEW>

Download and prepare reference files. Will download to directory specified by --directory argument.

=item B<-R    --REPAIR>

Download missing reference files only (e.g. after an interrupted --DOWNLOAD_NEW command). Will download to directory specified by --directory argument.

=item B<-h    --help>

Show help message.

=item B<-m    --manual>

Show manual page.


=back 

=cut

=head1 EXAMPLES 

        ensemblGeneAnnotator.pl -i vep_annotated.vcf -d ~/ensGeneAnnotatorDatabase -o vep_and_ensGene_annotated.vcf
        #annotate vep_annotated.vcf file using database folder located at ~/ensGeneAnnotatorDatabase
        
        ensemblGeneAnnotator.pl -i vep_annotated.vcf -d ~/ensGeneAnnotatorDatabase -o vep_and_ensGene_annotated.vcf --functional
        #as above but only considering genes which have at least one 'functional' variant consequence
        
        ensemblGeneAnnotator.pl -d ~/newEnsGeneAnnotatorDatabase -D
        #download files and create database in ~/newEnsGeneAnnotatorDatabase folder

        ensemblGeneAnnotator.pl -d ~/newEnsGeneAnnotatorDatabase -R
        #replace missing files in ~/newEnsGeneAnnotatorDatabase 


=head1 DESCRIPTION

This script reads VCF lines annotated with Ensembl's variant_effect_predictor.pl, identifies the corresponding human Entrez Gene ID for each ensembl gene and annotates information from Gene RIFS, Gene Ontology, NCBI summaries, OMIM and MGI phenotypes to each variants INFO field. In order to conform to VCF format, text in annotations has spaces replaced with underscores, semi-colons replaced with the ^ symbol and commas replaced with the ` symbol. Multiple values for annotations are separated with two colons ("::"). This is designed to allow for conversion of fields to readable text with a simple regex when attempting to report output in a human readable format.

This program must create a local database of gene annotations to search, which is a time consuming process but only needs to be done occasionally, when you wish to update the stored annotations. Alternatively, a prebuilt database may be downloaded from http://sourceforge.net/projects/vcfhacks/files/ensAnnotatorDatabase/ 

=cut

=head1 AUTHOR

David A. Parry
University of Leeds


=head1 COPYRIGHT AND LICENSE

Copyright 2013, 2014  David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut

