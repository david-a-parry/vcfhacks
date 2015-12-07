#!/usr/bin/env perl 

use warnings;
use strict;
use Config;
use LWP::Simple;
use Getopt::Long qw(:config no_ignore_case);
use POSIX qw/strftime/;
use Pod::Usage;
use Data::Dumper;
use Net::FTP;
use Cwd;
use IO::File;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use IO::Uncompress::Unzip qw(unzip $UnzipError) ;
use File::Path qw(make_path remove_tree);
use File::Temp qw/ tempfile tempdir /;
use File::Copy;
use File::Basename;
use File::stat;
use Cwd;
use Term::ProgressBar;

use FindBin qw($RealBin);
use lib "$FindBin::Bin/lib";
use lib "$FindBin::Bin/lib/Bioperl";
use lib "$FindBin::Bin/lib/BioASN1EntrezGene/lib";
use VcfhacksUtils;
use VcfReader; 

my @add_classes;
my @annotations;
my @biotypes;
my @classes;
my $rest_server = "http://grch37.rest.ensembl.org" ;
my %entid_cache = (); #store results of ENTREZ ID searches here
my %opts = (
    a   => \@add_classes,
    b   => \@biotypes,
    c   => \@classes,
    g   => \@annotations,
); 

GetOptions(
    \%opts,
    "D|DOWNLOAD_NEW",
    "P|PREPARE",
    "R|REPAIR",
    "a|add=s{,}",
    "b|biotype_filters=s{,}",
    "c|classes=s{,}",
    "d|directory=s",
    "f|functional",
    "g|gene_annotations=s{,}",
    "h|?|help",
    "i|input=s",
    "l|list_mode",
    "manual", 
    "m|mode=s",
    "n|no_rest_queries",
    "no_biotype_filtering",
    "o|output=s",
    "p|progress", 
    "r|rest_server=s{,}" => \$rest_server,
) or pod2usage( -message => "Syntax error", -exitval => 2 );
pod2usage( -verbose => 2, ) if $opts{manual};
pod2usage( -verbose => 1, ) if $opts{h};
pod2usage( -exitval => 2, -message => "--input is required", )
  if not $opts{i}
  and not $opts{D}
  and not $opts{P}
  and not $opts{R};

my $http;
if (not $opts{n}){
    if (eval "use HTTP::Tiny; 1"){
        if (eval "use JSON; 1"){
            $http = HTTP::Tiny->new();
        }else{
            informUser
            ( 
                "[WARNING] JSON module is required for using rest ".
                "queries to find missing ENTREZ IDs. If you want to ".
                "use this feature please install this module using CPAN.\n"
            );  
            $opts{n} = 1;
        }
    }else{
        informUser
        (
            "[WARNING] HTTP::Tiny is required for using rest queries to ".
            "find missing ENTREZ IDs. If you want to use this feature please ".
            "install this module using CPAN.\n"
        );  
        $opts{n} = 1;
    }
}

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



my ( $script, $script_dir ) = fileparse($0);
$opts{d} = "$RealBin/ensAnnotatorDb" if ( not $opts{d} );
$opts{d} =~ s/\/$//;
my $OUT;
if ($opts{o}) {
    open( $OUT, ">$opts{o}" ) or die "Can't open $opts{o} for writing: $!\n";
}
else {
    $OUT = \*STDOUT;
}
my %file_lengths = ();
my %database     = (
    ensemblToEntrez => {
        localfile => "$opts{d}/ensemblToEntrez",
        col       => 0,
        delimiter => "\t",
        url       => undef,
        dir       => undef,
        file      => "ensemblToEntrez" 
    },    #create this one on the fly from human_summary

#    gene2refseq => {
#        localfile => "$opts{d}/gene2refseq",
#        col       => 3,
#        delimiter => "\t",
#        url       => "ftp.ncbi.nlm.nih.gov",
#        dir       => "gene/DATA",
#        file      => "gene2refseq.gz"
#    },#this file is big and not used by anything except maybe by a custom made SnpEff database
#    #might be worth including at a later date but will need to implement some processing
#    to reduce the filesize (filter on taxonomy and throwaway entries without nucleotide IDs)

    human_summary => {
        localfile => "$opts{d}/Homo_sapiens_ncbi_gene_all_summaries.txt",
        col       => 0,
        delimiter => "\t",
        url       => "ftp.ncbi.nlm.nih.gov",
        dir       => "gene/DATA/ASN_BINARY/Mammalia",
        file      => "Homo_sapiens.ags.gz"
    },
    mouse_phenotype_desc => {
        localfile => "$opts{d}/VOC_MammalianPhenotype.rpt",
        col       => 0,
        delimiter => "\t",
        url       => "ftp.informatics.jax.org",
        dir       => "pub/reports/",
        file      => "VOC_MammalianPhenotype.rpt"
    },
    mgi2entrez => {
        localfile => "$opts{d}/HMD_HumanPhenotype.rpt",
        col       => 1,
        delimiter => "\t",
        url       => "ftp.informatics.jax.org",
        dir       => "pub/reports/",
        file      => "HMD_HumanPhenotype.rpt"
    },
    mim_morbid => {
        localfile => "$opts{d}/morbidmap",
        col       => 2,
        delimiter => "\\|",
        url       => "ftp.omim.org",
        dir       => "OMIM",
        file      => "morbidmap"
    },
    biogrid    => {
        localfile => "$opts{d}/BIOGRID-ALL-3.4.128.tab2.txt",
        col       => 1,
        delimiter => "\t",
        url       => "http://thebiogrid.org",
        dir       => "downloads/archives/Release%20Archive/BIOGRID-3.4.128/",
        file      => "BIOGRID-ALL-3.4.128.tab2.zip"
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
if ($opts{D}) {
    prepare_database( \%database );
    exit if not $opts{i};
}
elsif ($opts{R}) {
    prepare_database( \%database, "repair" );
    exit if not $opts{i};
}
else {
    my @missing_files = ();
    foreach my $m (@missing) {
        push @missing_files, $m->{localfile};
    }
    display_error_and_exit
    (
        "Missing database files - rerun with the --DOWNLOAD_NEW option ".
        "to correct this. Can't find following files:\n" . 
        join( "\n", @missing_files ) . "\n"
    ) if @missing;
}
if ($opts{P}) {
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
my @all_annotations = 
qw(
    ENSGENE_ID
    ENTREZ_ID
    SYMBOL
    GO_ID
    GO_DESCRIPTION
    GENERIFS
    BIOGRID_INTERACTANTS
    SUMMARY
    OMIM
    MGI_PHENOTYPE
);
if (@annotations){
    foreach my $anno (@annotations){
        if (not grep {/uc($_) eq uc($anno)/} @annotations){
            die "Unrecognised --gene_annotation \"$anno\" specified\n";
        }
    }
}else{
    @annotations = @all_annotations;
}

my $pre_progressbar;
if ($opts{p}){
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
    if ($opts{p}){
        $next_update = $pre_progressbar->update($pre_progress)
            if $pre_progress >= $next_update;
    }
}

my $progressbar; 
my %symbol_to_id = ();
if ($opts{l}){
    parseAsList(); 
}else{
    parseAsVcf();
}


######################################################
sub getEntrezIdFromSymbol{
    my $s = shift;
    if (not %symbol_to_id){
        getSymbolsToIds();
    }
    if (exists $symbol_to_id{$s}){
        return $symbol_to_id{$s};
    }
    return;
}

######################################################
sub getSymbolsToIds{
    open (my $FH, $database{human_summary}->{localfile}) 
      or die "Cannot open database file for reading: $!\n";
    while (my $line = <$FH>){
        chomp $line;
        my @spl = split("\t", $line);
        $symbol_to_id{$spl[1]} = $spl[0];
    }
    close $FH;
}


######################################################
sub parseAsList{
    open (my $LIST, $opts{i}) or die "Could not open input ($opts{i}): $!\n";
    print $OUT join(",", @annotations) . "\n";
    while (my $line = <$LIST>){
        next if $line =~ /^#/;
        chomp $line;
        next if not $line;
        my @entrez_ids = (); 
        my $id = (split "\t", $line)[0];
        my $entrez_id;
        if ($id =~ /^\d+$/){
            push @entrez_ids, $id;
        }elsif ($id =~ /^ENSG\d+(\.\d+)*$/){
            if ($1){
                $id =~ s/\.\d+$//; 
            }
            push @entrez_ids, search_ensembl_id($id); 
        }else{
            push @entrez_ids, getEntrezIdFromSymbol($id);
        }
        if (not @entrez_ids){
            print $OUT join(",", ( "$id-NOT_FOUND")  x @annotations)  . "\n";
            next;
        }
        foreach my $entrez_id (@entrez_ids){
            my @single_anno = ();
            my $gene_annot = searchWithEntrezId($entrez_id);
            foreach my $annot (@annotations){
                if (lc$annot eq 'entrez_id'){
                    push @single_anno, $entrez_id; 
                }elsif (exists $gene_annot->{lc($annot)}){
                    if (ref ($gene_annot->{lc($annot)}) eq 'ARRAY'){
                        @{$gene_annot->{lc($annot)}} = VcfhacksUtils::removeDups
                        ( 
                            @{ $gene_annot->{lc($annot)} } 
                        );
                        my @conv = ();
                        foreach my $g (@{$gene_annot->{lc($annot)}}){
                            push @conv, $g;
                        }
                        
                        my $joined = join("|", @conv);
                        if ($joined =~ /\s/){
                            push @single_anno, "\"$joined\"";
                        }else{
                            push @single_anno, $joined;
                        }
                    }else{
                        if ($gene_annot->{lc($annot)} =~ /\s/){
                            push @single_anno, "\"$gene_annot->{lc($annot)}\"";
                        }else{
                            push @single_anno, "$gene_annot->{lc($annot)}";
                        }
                    }
                }else{
                    push @single_anno, '';
                }
            }
            print $OUT join(",", @single_anno) . "\n";
        }
    }
    close $LIST;
}

######################################################
sub parseAsVcf{
    print STDERR "Initialising input VCF...\n";
    my @header = VcfReader::getHeader($opts{i});
    die "ERROR: Invalid VCF header in $opts{i}\n" 
      if not VcfReader::checkHeader(header => \@header);
    
    my %csq_header = getAndCheckCsqHeader(\@header);
    my @csq_fields = getCsqFields(\%csq_header);#consequence fields to retrieve from VCF
    my %class_filters = map { $_ => undef } getAndCheckClasses();
    my %biotype_filters = map { $_ => undef } getAndCheckBiotypes();
    
    my $total_vars  = 0; 
    my $vcf_line = 0;
    $next_update = 0;
    if ($opts{p}){
        informUser("Counting variants in input for progress monitoring.\n"); 
        $total_vars = VcfReader::countVariants($opts{i});
        informUser("$opts{i} has $total_vars variants.\n");
        $progressbar = Term::ProgressBar->new(
        {
            name  => "Annotating",
            count => $total_vars,
            ETA   => "linear",
        });
    }
    my %gene_inf = 
    (
        ID          => "GeneAnno",
        Number      => ".",
        Type        => "String",
        Description => "Collected Entrez/MGI/BioGrid gene annotations from for".
                       " VEP/SnpEff annotated human genes. Multiple values per".
                       " annotation are separated using two colons ('::'), " . 
                       "spaces are replaced with underscores, commas are " .
                       "replaced with the ` symbol, and semi-colons are " . 
                       "replaced with the ^ symbol so that regexes can be used".
                       " to extract the original text programmatically. ".
                       "Format: " . join("|", @annotations) 
    );

    print $OUT join("\n", grep { /^##/ } @header) . "\n";
    #add header line detailing program options
    print $OUT VcfhacksUtils::getInfoHeader(%gene_inf) . "\n"; 
    print $OUT VcfhacksUtils::getOptsVcfHeader(%opts) . "\n"; 
    print $OUT "$header[-1]\n";

    my $gene_field = "gene";#if using VEP annotations we get Gene ID from this key
    if ($opts{m} eq 'snpeff'){
        $gene_field = "gene_id";#if using SnpEff annotations we get Gene ID from this key
    }
    
    my $VCF = VcfReader::openVcf($opts{i}); 
LINE: while ( my $line = <$VCF> ) {
        next if $line =~ /^#/;#skip header
        chomp $line;
        $vcf_line++;
        my @split = split("\t", $line); 
        my @entrez_ids = ();
        my @csq = getConsequences(\@split, \@csq_fields, \%csq_header);
    
        my %gene_annot = ();
        die "No consequence field found for line:\n$line\nPlease annotate ".
            " your VCF file with ensembl's variant effect precictor or SnpEff ".
            "before running this program.\n" if not @csq;
        foreach my $c (@csq) {
            if ($opts{f}) {
                next if ( not check_consequence
                    ( 
                        $c, 
                        \%biotype_filters, 
                        \%class_filters 
                    ) 
                );
            }
            if ($c->{$gene_field} =~ /^ENS/){
                my $ensid = $c->{$gene_field};
                if ($ensid =~ /ENS[A-Z]*\d+-ENS/){
                    my @e = split("-", $ensid); 
                    foreach my $id (@e){
                        push @entrez_ids, search_ensembl_id($id);
                    }
                }else{ 
                    push @entrez_ids, search_ensembl_id($ensid); 
                }
            }elsif($c->{$gene_field} =~ /^\d+$/){
                if ($c->{$gene_field} =~ /\d+-\d/){
                    push @entrez_ids, split("-", $c->{$gene_field}); 
                }else{ 
                    push @entrez_ids, $c->{$gene_field} ; 
                }
            }elsif($c->{$gene_field} =~ /^[NX][MR]_\d+$/){
                #not implemented unless we decide to include gene2refseq in database
               # push @entrez_ids, search_refseq_id($c, $gene_field); 
            }
        }
        foreach my $id (@entrez_ids){
            $gene_annot{$id} = searchWithEntrezId($id);
        }

        if ( not keys %gene_annot) {
            my $ens_annot = "|" x @annotations;
            my $new_line = VcfReader::addVariantInfoField
            (
                line  => \@split, 
                id    => 'GeneAnno', 
                value => $ens_annot
            );
            print $OUT join("\t", @$new_line) . "\n";
            next LINE;
        }
        my @gene_anno = ();
        foreach my $entrez_id ( keys %gene_annot ) {
            my @single_anno = ();
            foreach my $annot (@annotations){
                if (lc$annot eq 'entrez_id'){
                    push @single_anno, $entrez_id;    
                }elsif (exists $gene_annot{$entrez_id}->{lc($annot)}){
                    if (ref ($gene_annot{$entrez_id}->{lc($annot)}) eq 'ARRAY'){
                        
                        @{$gene_annot{$entrez_id}->{lc($annot)}} =
                          VcfhacksUtils::removeDups
                          ( 
                             @{ $gene_annot{$entrez_id}->{lc($annot)} }
                          );
                        my @conv = ();
                        foreach my $g (@{$gene_annot{$entrez_id}->{lc($annot)}}){
                            push @conv, convert_text($g);
                        }
                        push @single_anno, join("::", @conv);
                    }else{
                        my $converted = convert_text($gene_annot{$entrez_id}->{lc($annot)});
                        push @single_anno, $converted;
                    }
                }else{
                    push @single_anno, '';
                }
            }
            push @gene_anno, join("|", @single_anno);
        }
        my $ens_annot = join(",", @gene_anno);
        my $new_line = VcfReader::addVariantInfoField
        (
            line  => \@split, 
            id    => 'GeneAnno', 
            value => $ens_annot
        );
        print $OUT join("\t", @$new_line) . "\n";
        if ($opts{p}){
            $next_update = $progressbar->update($vcf_line) 
              if $vcf_line >= $next_update;
        }
    }
    close $VCF;
    close $OUT;
    if ($opts{p}){
        $progressbar->update( $total_vars )
            if $total_vars >= $next_update;
    }
    for my $k ( keys %database ) {
        close $database{$k}->{fh};
        close $database{$k}->{idx};
    }
    print STDERR "\nAnnotation finished.\n";

}

#################################################
sub getCsqFields{
    my $csq_header = shift;
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
    }else{
        @fields = 
        qw(
            allele
            annotation
            gene_name
            gene_id
            feature_type
            feature_id
            transcript_biotype
        );
    }
    foreach my $f (@fields){
        if (not exists $csq_header->{$f}){
            die "Could not find '$f' field in $opts{m} consequence header " .
              "- please ensure you have annotated your file including the appropriate ".
              "fields.\n";
        }
    }
    return @fields;
}

#################################################
sub getConsequences{
    my ($line, $csq_fields, $csq_head) = @_;
    if ($opts{m} eq 'vep'){
        return VcfReader::getVepFields
        ( 
            line        => $line,
            field       => $csq_fields,
            vep_header  => $csq_head,
        );
    }else{
        return VcfReader::getSnpEffFields
        ( 
            line          => $line,
            field         => $csq_fields,
            snpeff_header => $csq_head,
        );
    }
}

#################################################


sub getAndCheckCsqHeader{
    my $head = shift;
    my %csq_head = ();
    if (not $opts{m}){
        eval { 
            %csq_head = VcfReader::readVepHeader
            (
                header => $head
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
                    header => $head
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
                    header => $head
                ); 
        }else{
            %csq_head = VcfReader::readSnpEffHeader
            (
                header => $head
            ); 
        }
    }
    return %csq_head;
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
 
######################################################
sub searchWithEntrezId{ 
    #USE ENTREZ IDs TO SEARCH DATABSE FILES
    my $entrez_id = shift;
    my %an = ();
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
        $an{symbol}  = $humsum_line[1];
        $an{summary} = $humsum_line[2];
        push @{ $an{ensgene_id} }, split(/\|/, $humsum_line[3]);
        $an{go_id}   = $humsum_line[4];
        $an{go_description} = $humsum_line[5];
        $an{generifs}     = $humsum_line[6];

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
                        push( @{ $an{omim} },
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
                        push( @{ $an{omim} },
                            $mim_line[0] );
                    }
                    else {
                        last MIM_DOWN;
                    }
                }
            }
            else {
                push( @{ $an{omim} }, "" );
            }
        }
        else {
            push( @{ $an{omim} }, "" );
        }
    }else {
        $an{symbol}  = "";
        $an{summary} = "";
        $an{go_id}   = "";
        $an{go_description} = "";
        $an{generifs}     = "";

        push( @{ $an{omim} }, "" );
    }
    my @interactants = get_interactants($entrez_id);
    @interactants 
      ? push @{ $an{biogrid_interactants} }, @interactants,
      : push @{ $an{biogrid_interactants} }, "";
    my @mgi = get_MGI_phenotype($entrez_id) ;
    @mgi
      ? push @{ $an{mgi_phenotype} }, @mgi
      : push @{ $an{mgi_phenotype} }, "";
    return \%an;
}
######################################################
sub informUser{
    my $msg = shift;
    my $time = strftime( "%H:%M:%S", localtime );
    if ($progressbar){
        $progressbar->message( "[INFO - $time] $msg" );
    }else{
        print STDERR "[INFO - $time] $msg";
    }
}
######################################################
sub ensRestQuery{
    my $url = shift;
    my $response = $http->get($url, {
          headers => { 'Content-type' => 'application/json' }
    });
    if (not $response->{success} ){
        informUser("[WARNING] Ensembl REST query ('$url') failed!\n");
        return ;
    }
    if(length $response->{content}) {
      return decode_json($response->{content});
    }
    informUser("[WARNING] No content for Ensembl REST query ('$url')!\n");
}
########################################
sub search_refseq_id{
    my ($csq, $gene_field) = @_;
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
sub query_rest_ensid{
    my $ensid = shift;
    my @ids = ();
    my $rest_url = "$rest_server/xrefs/id/$ensid?external_db=ENTREZGENE";
    informUser("[INFO] Querying Ensembl REST server for Gene ID $ensid...\n");
    my $xref_array = ensRestQuery($rest_url);
    if (ref $xref_array eq 'ARRAY'){
        foreach my $xref (@$xref_array){
            if (exists $xref->{primary_id}){
                push @ids, $xref->{primary_id};
            }
        }
    }
    
    $entid_cache{$ensid} = \@ids;
    return @ids;
}

########################################
sub search_ensembl_id{
    my ($ensid) = @_;
    if (exists $entid_cache{$ensid}){
        if (ref $entid_cache{$ensid} eq 'ARRAY' ){ 
            return @{ $entid_cache{$ensid} } ;
        }
        return;
    }
    my @ids = (); 
    my $i = binSearchLineWithIndex(
        $ensid,
        $database{ensemblToEntrez}->{fh},
        $database{ensemblToEntrez}->{idx},
        $database{ensemblToEntrez}->{length},
        $database{ensemblToEntrez}->{col}
    );
    if ( $i < 1 ){
        return query_rest_ensid($ensid) if (not $opts{n});
        return;
    }
            
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
        if ( $ens_line[ $database{ensemblToEntrez}->{col} ] eq $ensid )
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
        if ( $ens_line[ $database{ensemblToEntrez}->{col} ] eq $ensid )
        {
            push @ids, $ens_line[1];
        }
        else {
            last;
        }
    }
    $entid_cache{$ensid} = \@ids;
    return @ids;
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
    if ($opts{m} eq 'snpeff'){
        return check_snpeff_consequence(@_);
    }else{
        return check_vep_consequence(@_);
    }
}
########################################
sub check_snpeff_consequence {
    my ($annot, $biotype_filters, $class_filters) = @_;
    #skip variants with undef features (intergenic variants)
    return 0 if not defined $annot->{feature_id};
    #skip unwanted biotypes
    return 0 if exists $biotype_filters->{lc $annot->{transcript_biotype} };
    my @anno_csq = split( /\&/, $annot->{annotation} );
ANNO: foreach my $ac (@anno_csq){
        $ac = lc($ac);#we've already converted %class_filters to all lowercase
        if ( exists $class_filters->{$ac} ){
            return 1;
        }
    }
    return 0;#no annotation matching %class_filters
}

########################################
sub check_vep_consequence {
    my ($annot, $biotype_filters, $class_filters) = @_;
    #intergenic variants have no feature associated with them - skip
    return 0 if $annot->{consequence} eq "intergenic_variant";
    #skip unwanted biotypes
    return 0 if (exists $biotype_filters->{$annot->{biotype}}) ;
    #skip non-canonical transcripts if --canonical_only selected
    if ($opts{canonical_only}) {
        return 0 if ( not $annot->{canonical} );
    }
    
    my @anno_csq = split( /\&/, $annot->{consequence} );
    #skip NMD transcripts
    return 0 if ( grep { /NMD_transcript_variant/i } @anno_csq );

ANNO: foreach my $ac (@anno_csq){
        $ac = lc($ac);#we've already converted %class_filters to all lowercase
        if ( exists $class_filters->{$ac} ){
            return 1;
        }
    }
    return 0;#no annotation matching %class_filters
}
########################################
sub prepare_files {
    my ( $file_array_ref, $database_ref ) = @_;
    foreach my $file (@$file_array_ref) {
        sort_and_index_gene_files(
            $file->{localfile}, 
            $file->{col}, 
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
"Oldest file is $oldest days old. Use the -D/--DOWNLOAD_NEW option to update."
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
    #my $progressbar = Term::ProgressBar->new(
    #    { name => "Prep Database", count => 100, ETA => "linear", } );
    if ( not -e $opts{d} ) {
        unless (mkdir $opts{d}){
            display_error_and_exit(
                    "Permissions Error",
"Can't create directory $opts{d} for database files - please check permissions"
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
        chdir $opts{d} or display_error_and_exit(
            "Directory Error",
            "Can't move to directory $opts{d}",
        );
        my $file_exists = 0;

        if ( -e $file->{localfile} ) {
            move( $file->{localfile} , "$file->{localfile}.bkup" ) or display_error_and_exit(
                "File Error",
                "Error creating file backup of $file->{localfile}",
                "Check permissions and/or disk space."
            );
            $file_exists++;
        }
        if ($file->{url} eq "http://thebiogrid.org"){
            chdir $dir;
            downloadBiogrid($file, $file_exists);
            next; 
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
        if ( $file->{file} =~ /\.gz$/ ) {
            $time = strftime( "%H:%M:%S", localtime );
            print STDERR "[$time] Decompressing $file->{file}...\n";
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
            }
            else {
                gunzip( $file->{file} => $output )
                  or restore_file( $file_exists, $file )
                  && display_error_and_exit( "Error decompressing file",
                    "Error decompressing $file->{file}" );
            }

        }
        if ( $file->{file} =~ /\.ags/ ) {
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
#unlink $file->{file} or display_error_and_exit( "Can't delete xml output ($file->{file})", "Check permissions - it is safe to manually delete this file now");
            my ( $enstoEntrez_file_name, $file_dir ) 
               = fileparse( $database_ref->{ensemblToEntrez}->{localfile} );
            extract_ncbi_summaries( $gene2xml, $decomp_file, "$file_name.tmp",
                $enstoEntrez_file_name .".tmp" );
                move( "$file_name.tmp", $file_name ) or display_error_and_exit(
                "File Error",
                "Error creating file $file_name",
                "Check permissions and/or disk space."
                );
                move( $enstoEntrez_file_name .".tmp", $enstoEntrez_file_name) 
                  or display_error_and_exit
                (
                    "File Error",
                    "Error creating file $database_ref->{ensemblToEntrez}->{localfile}",
                    "Check permissions and/or disk space."
                );
            $time = strftime( "%H:%M:%S", localtime );
            print STDERR "[$time] Sorting and indexing ensemblToEntrez file...\n";
            sort_and_index_gene_files(
                $enstoEntrez_file_name,
                $database_ref->{ensemblToEntrez}->{col},
                $database_ref->{ensemblToEntrez}->{delimiter}
            );
            unlink $decomp_file or display_error_and_continue(
                "Can't delete decompressed ags file ($decomp_file)",
"Check permissions - it is safe to manually delete this file now"
            );

#unlink $xml_out or display_error_and_exit( "Can't delete xml output ($xml_out)", "Check permissions - it is safe to manually delete this file now");
        }
        chdir $dir;
        $time = strftime( "%H:%M:%S", localtime );
        print STDERR "[$time] Sorting and indexing $file_name...\n";
        sort_and_index_gene_files( "$file_dir/$file_name", $file->{col},
            $file->{delimiter} );
        if (-e "$file->{localfile}.bkup"){
            unlink "$file->{localfile}.bkup"
                or display_error_and_continue(
                "Can't delete backup file \"$file->{localfile}.bkup\"",
                "Check permissions - it is safe to manually delete this file now" );
        }
    }
    #$progressbar->update(100);
    my $time = strftime( "%H:%M:%S", localtime );
    print STDERR "[$time] Database update finished.\n";
}

#########################################
sub downloadBiogrid{
    my $file = shift;
    my $exists = shift;
    my $url = "$file->{url}/$file->{dir}/$file->{file}";
    my $dl  = "$opts{d}/$file->{file}";
    my $time = strftime( "%H:%M:%S", localtime );
    print STDERR "[$time] Downloading $url...\n";
    getstore($url, $dl )
          or restore_file( $exists, $file )
          && display_error_and_exit( 
            "Download error",
            "Error downloading $url!\n"
    );
    $time = strftime( "%H:%M:%S", localtime );
    print STDERR "[$time] Decompressing $dl...\n";
    unzip $dl => $file->{localfile} or die "Unzip failed: $UnzipError\n";   
    unlink($dl) 
                or display_error_and_continue(
                "Can't delete backup file \"$dl\"",
                "Check permissions - it is safe to manually delete this file now" );
    print STDERR "[$time] Sorting and indexing $file->{localfile}...\n";
    sort_and_index_gene_files(
        $file->{localfile},
        $file->{col},
        $file->{delimiter},
    );
    if (-e "$file->{localfile}.bkup"){
        unlink "$file->{localfile}.bkup"
            or display_error_and_continue(
            "Can't delete backup file \"$file->{localfile}.bkup\"",
            "Check permissions - it is safe to manually delete this file now" );
    }
}

#########################################
sub sort_and_index_gene_files {
    my ($file,$sort_column, $delimiter) = @_;    #column no. is 0 based
    open( my $FILE, $file )
      or die "Can't open $file for sorting and indexing!\n";
    my ( $file_short, $dir ) = fileparse($file);
    my @lines;
    $delimiter = "\t" if not $delimiter;
    my $check_taxon = 0;
    $check_taxon = 1
      if ( $file =~ /gene2go/ or $file =~ /generifs/ )
      ;    #for these files we'll ignore non-mouse and non-human genes
    while (<$FILE>) {
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
    close $FILE;
    @lines = sort { $a->[1] cmp $b->[1] } @lines;
    my ( $tmp, $TEMP );
    ( $TEMP, $tmp ) = tempfile( "$dir/tmp_dharmaXXXX", UNLINK => 1 )
      or die "Can't create temporary sort file\n";
    foreach my $line (@lines) {
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


#####################################
sub look_forward_and_back{
    my ( $x, $FILE, $FILE_INDEX, $u, $column, $start, $delimiter ) = @_;
    $delimiter = "\t" if not $delimiter;
    my @lines = (); 
    push @lines, line_with_index( $FILE, $FILE_INDEX, $start );
    for (my $i = $start - 1; $i > 0; $i--){#look back
        my $line = line_with_index( $FILE, $FILE_INDEX, $i );
        my @split = split(/$delimiter/, $line); 
        if ($split[$column] eq $x){
            push @lines, $line;
        }else{
            last;
        }
    }
    for (my $i = $start + 1; $i < $u; $i++){#look forward 
        my $line = line_with_index( $FILE, $FILE_INDEX, $i );
        my @split = split(/$delimiter/, $line); 
        if ($split[$column] eq $x){
            push @lines, $line;
        }else{
            last;
        }
    }
    return @lines;
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
            (my $value = $line[$column] ) =~ s/\s+//g; #some files have extra whitespace...
            if ( $x lt $value ) {
                $u = $i - 1;
            }
            elsif ( $x gt $value ) {
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
    print STDERR "[$time] Converting NCBI binary gene data to text from ". 
                 "$agsfile using the following command:\n$cmd\n";
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
    print STDERR "[$time] Extracting NCBI gene data from $agsfile - ". 
                 "this may take some time...\n";
    while ( my ( $gene, $genestructure, $uncaptured ) = $io->next_seq ) {
    
        my $chrom;
        my @go_ids      = ();
        my @go_descs    = ();
        my @rifs        = ();
        my @mim         = ();
        my @ensembl     = ();
        my $annotation_ref = $gene->annotation;
        my @annotations = $annotation_ref->get_Annotations();
        
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
        @annotations = $annotation_ref->get_Annotations('OntologyTerm');
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
   
        print $SUMOUT "\n";
    }
    
    close $SUMOUT;
    close $ENSOUT;
    unlink "$agsfile.asn.1" or warn "Error deleting $agsfile.asn.1 file - ".
           "you may want to delete this manually.\n";
}

#####################################
sub display_error_and_exit {
    my ( $text, $informative_text, ) = @_;

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
sub get_interactants{
    my $gene_id = shift;
    my @inter = ();
    my @symbols = ();
    my $i = binSearchLineWithIndex
      ( 
        $gene_id, 
        $database{biogrid}->{fh}, 
        $database{biogrid}->{idx},
        $database{biogrid}->{length},
        1 
      );
    if ( $i > 0 ) {
        my @hits = look_forward_and_back
          ( 
            $gene_id, 
            $database{biogrid}->{fh},
            $database{biogrid}->{idx},
            $database{biogrid}->{length},
            1,
            $i,
        );
        foreach my $hit (@hits){
            push @inter, (split "\t", $hit)[2];
        }
    }
    my %seen = ();
    @inter = grep { !$seen{$_}++ }  @inter;
    foreach my $inter (@inter){
        my $i = binSearchLineWithIndex(
            $inter,
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
            push @symbols, $humsum_line[1];
        }
        
    }
    @symbols = sort @symbols;
    return @symbols;
}
     
#########################################
sub get_MGI_phenotype {
    my $gene_id = shift;
    my @desc = ();
    my $i = binSearchLineWithIndex
      ( 
        $gene_id, 
        $database{mgi2entrez}->{fh}, 
        $database{mgi2entrez}->{idx},
        $database{mgi2entrez}->{length},
        1 
      );
    if ( $i > 0 ) {
        my @phenos = ();
        my @hits = look_forward_and_back
          ( 
            $gene_id, 
            $database{mgi2entrez}->{fh},
            $database{mgi2entrez}->{idx},
            $database{mgi2entrez}->{length},
            1,
            $i,
        );
        foreach my $hit (@hits){
            push @phenos, split( /\s+/,
              ( split "\t", $hit )[5]
            );
            @phenos = grep { !/^\s*$/ } @phenos;#some MGI files have extra whitespace
        }
        foreach my $ph (@phenos) {
            my $k = binSearchLineWithIndex
              ( 
                $ph, 
                $database{mouse_phenotype_desc}->{fh},
                $database{mouse_phenotype_desc}->{idx},
                $database{mouse_phenotype_desc}->{length},
                0,
              );

            if ( $k > 0 ) {
                push(
                    @desc,
                    (
                        split "\t",
                        line_with_index(
                            $database{mouse_phenotype_desc}->{fh},
                            $database{mouse_phenotype_desc}->{idx},
                            $k,
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
        display_error_and_exit
        (
            "Unsupported OS: $^O",
            "Failed to retrieve gene2xml because your OS is not supported. ".
            "If you have a binary place it in your database folder and rerun ".
            "the database update"
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

Input VCF file annotated with variant_effect_predictor.pl or snpEff. Alternatively, if using --list_mode, a list of gene symbols or Enztrez Gene IDs. Required.

=item B<-o    --output>

Output file name.

=item B<-d    --directory>

Directory containing reference files. Will look for a folder called 'ensAnnotatorDb' in the same directory as this program if not supplied. This directory will be created and populated if using the --DOWNLOAD_NEW option.

=item B<-f    --functional>

Use this flag to only annotate standard 'functional' variant classes.  Prevents annotation of gene information for genes/transcripts that overlap a variant but are not affected in a way defined by one of these variant classes. The default classes are:

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

Available classes that can be chosen instead of (or in addition to - see below) these classes can be found in the data/vep_classes.tsv and data/snpeff_classes.tsv files respectively using the -c/--classes. 

=item B<-c    --classes>

Use this to specify VEP/SnpEff variant classes to annotate. Gene information will only be annotated for genes affected by one of these variant classes. Overrides --functional option.

=item B<-a    --additional_classes>

Use this to specify additional VEP variant classes to annotate as well as those used by --functional.

=item B<-b    --biotype_filters>

When using the -f/--functional or -c/--classes options, features/transcripts with the following biotypes are ignored by default:

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

The 'data/biotypes.tsv' file contains a list of valid biotypes and the default behaviour of this program (i.e. 'keep' or 'filter'). 

=item B<--no_biotype_filtering>

Use this flag to consider consequences affecting ALL biotypes when using -f/--functional or -c/--classes options.

=item B<-g    --gene_annotations>

List of gene annotations to include in output. By default all of the following classes are included:

    ENSGENE_ID
    ENTREZ_ID
    SYMBOL
    GO_ID
    GO_DESCRIPTION
    GENERIFS
    BIOGRID_INTERACTANTS 
    SUMMARY
    OMIM
    MGI_PHENOTYPE

Specify one or more of these to limit the annotations in your output to these classes only.

=item B<-l    --list_mode>

Use this flag to indicate your input is not VCF but a list of gene symbols or Entrez IDs to search and annotate. If multiple columns are supplied only the first will be read. Lines beginning with a '#' will be ignored. Output will by in .csv format if using this flag. 

=item B<-r    --rest_server>

URL for Ensmbl REST server for gene ID queries. Default is "http://grch37.rest.ensembl.org". If using GRCh38 you may which to change to  "http://rest.ensembl.org". REST queries are only made when a gene ID can't be found in the local database. An internet connection is required.

=item B<-n    --no_rest_queries>

Use this flag to disable REST queries (e.g. if you have no internet connection). 

=item B<-m    --mode>

This program will attempt to detect the format of your input automatically by looking for VEP or SnpEff annotations, but you may specify either 'vep' or 'snpeff' with this option to select the mode employed by the script if you have a file with both VEP and SnpEff annotations. By default, if both annotations are present and the program is run without this option, VEP annotations will be used. NOTE: Only SnpEff annotations using the more recent 'ANN' style annotations rather than the older 'EFF' style are recognised by this program.

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

=item B<--manual>

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

=head1 COPYRIGHT AND LICENSE

Copyright 2013, 2014, 2015  David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut

