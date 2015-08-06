#!/usr/bin/env perl 
#David A. Parry, June 2015

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Excel::Writer::XLSX;
use Excel::Writer::XLSX::Utility;
use List::Util qw(first sum);
use Pod::Usage;
use File::Basename;
use FindBin;
use Tabix; 
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use lib "$FindBin::Bin/lib";
use VcfReader;
use TextToExcel;

my @allele_balance = ();
my $min_gq = 0; 
my %opts = (b => \@allele_balance, g => $min_gq);
GetOptions(
    \%opts,
    'i|input=s',        #vcf input
    'o|output=s',       #xlsx output
    'm|hgmd=s',         #vcf of HGMD variations converted with hgmd_to_vcf.pl
    't|transcripts=s',  #optional tsv list of transcript IDs in order of priority
    'c|clinvar=s',      #optional ClinVar VCF to add ClinVar CLINSIG annotations
    'b|allele_balance=f{,}', #min and optional max alt allele ratio per sample call
    'd|depth=i',         #optional min depth for sample call
    'a|allele_cutoff=f', #remove if allele is present in this proportion or more calls
    'g|gq=f',            #min GQ quality for calls
    'h|?|help',
) or pod2usage(-exitval => 2, -message => "Syntax error.\n"); 

pod2usage( -verbose => 1 ) if $opts{h};
pod2usage( -exitval => 2, -message => "-i/--input is required" ) if (not $opts{i});
pod2usage( -exitval => 2, -message => "-m/--hgmd is required" ) if (not $opts{m});

#open VCF, get samples and get VEP annotations

my @vhead              = VcfReader::getHeader($opts{i});
die "VCF header not OK for $opts{i}\n" if not VcfReader::checkHeader(header => \@vhead);
my %samples_to_columns = VcfReader::getSamples(header => \@vhead, get_columns => 1);
my %vep_fields         = VcfReader::readVepHeader(header => \@vhead);
my $FH                 = VcfReader::_openFileHandle($opts{i}); 

#check HGMD VCF and retrieve search arguments 

my %search_args = getSearchArgs();

# check ClinVar TSV file (https://github.com/macarthur-lab/clinvar) 
# if supplied and retrieve search arguments 

my %clinvar_sargs = getClinVarSearchArgs();
my %clnsig_codes = getClnSigCodes();

#set consequence ranks;
my %so_ranks = ();
setConsequenceRanks();

#open and check transcripts list

my %refseq_ranks = ();
setTranscriptsRanks();

#setup our hash of headers

my %headers = getHeaders();

#setup output
#Excel related global variables
my %sheets = ();
my $header_formatting;    
my $std_formatting;     
my $url_format;
my $xl_obj;
setupOutput();

#set up counts for variant indexes of each sheet
my %variant_counts = map { $_ => 0 } qw / HGMD LOF DamagingMissense BenignMissense Other / ;

#store variants per sample in this hash for writing sample sheet
my %sample_vars = (); 

#start reading variants
while (my $l = <$FH>){
    next if $l =~ /^#/;
    assessAndWriteVariant($l);
}
close $FH;
print STDERR "Done writing variant sheets - writing sample summary...\n"; 
writeSampleSummary();
print STDERR "Done.\n";
if ($xl_obj){
    $xl_obj->DESTROY();
}


###########################################################
sub writeSampleSummary{
    my $sample_sheet = $xl_obj->get_worksheets()->[$sheets{SampleSummary}];
    my $row = 1;
    foreach my $s (
            sort {$samples_to_columns{$a} <=> $samples_to_columns{$b}} 
            keys %samples_to_columns
    ){
        my $col = 0;
        foreach my $field ( @{$headers{sample}} ){
            if ($field eq 'Sample'){
                $sample_sheet->write($row, $col, $s)
            }elsif (exists $sample_vars{$s}->{$field} ){
                $sample_sheet->write($row, $col, scalar @{$sample_vars{$s}->{$field}});
                $sample_sheet->write_comment($row, $col,  join("\n", @{$sample_vars{$s}->{$field}}), visible => 0);
            }else{
                $sample_sheet->write($row, $col, 0);
            }
            $col++;
        }
        $row++;
    }
}

###########################################################
sub getClinVarSearchArgs{
    return if not $opts{c};
    #use bgzip compressed version of clinvar file from https://github.com/macarthur-lab/clinvar
    my ($bgz, $clinVarCols) = checkClinvarFile($opts{c});
    my $index = "$bgz.tbi";
    my $iterator = Tabix->new(-data =>  $bgz, -index => $index) ;
    my %sargs = ( tabix_iterator => $iterator, columns => $clinVarCols ); 
    return %sargs;
}
###########################################################
sub checkClinvarFile{
#returns columns hash of col name to col number
    my $cv = shift;
    if ($cv !~ /\.(b)*gz$/){
        print STDERR "ClinVar file ($cv) is missing a .gz extension - " .
            "attempting to compress with bgzip.\n";
        return compressClinVarTsv($cv); 
#    }elsif( not -e "$cv.tbi"){
#        print STDERR "No index found for ClinVar file ($cv). Attempting to index with tabix...\n";
#        return indexClinVar($cv);
    }else{
        my %col = getClinVarColumns($cv);
        return $cv, \%col;
    }
}
###########################################################
sub compressClinVarTsv{
    my $cv = shift;
    if (not `which bgzip`){
        die "Could not find bgzip executable - please install bgzip and ensure ".
            "it is in your PATH or compress and index $cv manually.\n";
    }
    open (my $CVAR, $cv) or die "Can't open $cv for reading header: $!\n";
    my @header = split("\t", <$CVAR>); 
    my %columns = getClinVarColumns($cv);
    print STDERR "Compressing $cv with bgzip...\n";
    system("bgzip $cv"); 
    checkExit($?);
    $cv = "$cv.gz";
#    indexClinVar($cv, \%columns); 
    return ($cv, \%columns);
}
###########################################################
sub indexClinVar{
    my $cv = shift;
    my %columns = ();
    if (my $col = shift){
        %columns = %$col;
    }else{
        %columns = getClinVarColumns($cv);
    }
    system("tabix -s $columns{chrom} -b $columns{pos} -e $columns{pos} $cv"); 
    checkExit($?);
    return %columns;
}
###########################################################
sub checkExit{
    my $exit = shift;
    return if not $exit;
    if ($exit == -1) {
        print "failed to execute: $!\n";
    }elsif ($exit & 127) {
        printf "command died with signal %d, %s coredump\n",
        ($exit & 127),  ($exit & 128) ? 'with' : 'without';
    }elsif ($exit) {
        printf "command exited with value %d\n", $exit >> 8;
    }
    die "Error executing command. Exiting.\n";
}
    
###########################################################
sub getClinVarColumns{
    my $cv = shift;
    my $CVAR; 
    if ($cv =~ /\.(b)*gz$/){
        $CVAR = new IO::Uncompress::Gunzip $cv, MultiStream => 1 
          or die "IO::Uncompress::Gunzip failed while opening $cv for reading: \n$GunzipError";
    }else{
        open ($CVAR, $cv) or die "Can't open $cv for reading header: $!\n";
    }
    my @header = split("\t", <$CVAR>); 
    my $n = 0;
    my %columns = map { lc($_) => $n++} @header;
    for my $c ( qw / 
                    chrom
                    pos
                    ref
                    alt
                    mut
                    clinical_significance
                    pathogenic
                    all_traits
                    all_pmids
                / ) { 
        if (not exists $columns{$c}){
            die "Required column ($c) not found in ClinVar file ($cv) header.\n";
        }
    }
    return %columns;
}
###########################################################
sub getSearchArgs{
    die "Header not ok for input ($opts{m}) "
      if not VcfReader::checkHeader( vcf => $opts{m} );
    return VcfReader::getSearchArguments($opts{m});
}
    
###########################################################
sub assessAndWriteVariant{
#identify transcript to report (make sure variant is functional and choose highest ranked transcript)
#write variant details to HGMD Excel sheet
#add var info to sample summary sheet and link to variant in main worksheet;
    my $line = shift;
    chomp ($line); 
    my @split = split("\t", $line); 
    #see if variant overlaps any variants in HGMD file
    my %min= VcfReader::minimizeAlleles(\@split);
    foreach my $al (sort {$a<=>$b} keys %min){
        my @hgmd_matches = searchHgmd($min{$al});
        my @clinvar_matches = searchClinVar($min{$al});
        writeToSheet(\@split, $min{$al}, \@hgmd_matches, \@clinvar_matches);
    }
}



###########################################################
sub writeToSheet{
    my ($l, $var, $matches, $clinvar_matches) = @_;
    #split line, minimized variant hash, HGMD matching lines, ClinVar matching lines
    
    #filter line if allele is too common
    if ($opts{a}){
        my %allele_counts = VcfReader::countAlleles
            (
                  line => $l,
                  minGQ => $min_gq,
            );
        my $total_allele_count = 0;
        foreach my $k (keys %allele_counts){
            $total_allele_count += $allele_counts{$k}; 
        }
        return if $allele_counts{$var->{ALT_INDEX}}/$total_allele_count >= $opts{a};
    }
    my @row = (); #values to write out to excel sheet
    my @split_cells = (); #values to write in split cells spanning row
    my %hgmd = ();#keys are HGMD fields, values are array refs of values
    #collect HGMD database annotations if found
    if (@$matches){
        my @hgmd_fields = 
        ( qw /
                hgmd_id
                disease
                variant_class
                gene_symbol
                hgvs
            /
        );
        foreach my $h (@$matches){
            my @match = split("\t", $h);
            foreach my $f (@hgmd_fields){
                push @{$hgmd{$f}}, 
                        VcfReader::getVariantInfoField(\@match, $f);
            }
        }
        foreach my $f (@hgmd_fields){
            push @row, join(",", @{$hgmd{$f}});
        }
    }
    #collect ClinVar database annotations if found
    my ($clnSig, $clnTraits, $cvarPathognic, $cvarConflicted) 
                         = getClinSig($clinvar_matches, $var);
    push @row, $clnSig, $clnTraits; 
    if ($cvarPathognic){
        push @row, $cvarConflicted;
    }
    
    #add standard VCF fields
    foreach my $f ( qw / CHROM ORIGINAL_POS ORIGINAL_REF ORIGINAL_ALT / ){
        push @row, $var->{$f};
    }
    #deal with VEP annotations
    my @get_vep =  qw /Symbol Consequence Allele feature canonical hgvsc hgvsp polyphen sift GERP++_RS/ ; 
    my @vep_csq = VcfReader::getVepFields
    (
        line       => $l,
        vep_header => \%vep_fields,
        field      => \@get_vep,
    );
    my $ref      = VcfReader::getVariantField($l, 'REF');
    my @alts     = split(",", VcfReader::getVariantField($l, 'ALT'));
    my @vep_alts = VcfReader::altsToVepAllele
    (
        ref  => $ref,
        alts => \@alts,
    );
    my %alt_to_vep = ();
    @alt_to_vep{@alts} = @vep_alts;
    my @csq_to_rank = ();
    foreach my $csq (@vep_csq){ 
        #collect info for each transcript before selecting which to use in report
        #skip consequences for different alleles if variant site is multiallelic
        next if $csq->{allele} ne $alt_to_vep{ $var->{ORIGINAL_ALT} };
        #if we've provided a list of RefSeq ranks skip if this gene symbol is not listed
        if (keys %refseq_ranks){
            my $sym = uc($csq->{symbol});
            next if not exists $refseq_ranks{$sym};
        }
        push @csq_to_rank , $csq;
    }
    my $csq_to_report = rankTranscriptsAndConsequences(\@csq_to_rank); 
    
    my @vep_fields = (qw / symbol feature consequence hgvsc hgvsp GERP++_RS/ );
    my $most_damaging_csq = getMostDamagingConsequence($csq_to_report);
    if (@$matches or 
        $most_damaging_csq eq 'missense_variant' or  
        $most_damaging_csq =~  /^inframe_(inser|dele)tion$/
    ){
        push @vep_fields, qw /polyphen sift/;
    }
    
    foreach my $f (@vep_fields){
        push @row, $csq_to_report->{$f};
    }
    my @cadd = split(",", VcfReader::getVariantInfoField($l, "CaddPhredScore"));
    my $allele_cadd = $cadd[ ($var->{ALT_INDEX} -1) ];
    push @row, $allele_cadd;
    #choose sheet to write to
    my $s_name;
    my %lofs = map {$_ => undef} qw /
        transcript_ablation  
        splice_acceptor_variant
        splice_donor_variant
        stop_gained
        frameshift_variant
        stop_lost
        start_lost
        transcript_amplification
    /;
    if (@$matches){
        $s_name = "HGMD";
    }elsif($cvarPathognic){
        $s_name = "ClinVarPathogenic";
    }elsif (exists $lofs{$most_damaging_csq}){
        $s_name = "LOF";
    }elsif ($most_damaging_csq eq 'missense_variant'){
        #check if damaging or benign... should we use Polyphen/SIFT as well/instead?
        if ($allele_cadd  >= 10 and 
            $csq_to_report->{polyphen} =~ /damaging/i and 
            $csq_to_report->{sift} =~ /deleterious/i
        ){
            $s_name = "DamagingMissense";
        }else{
            $s_name = "BenignMissense";
        }
    }elsif($most_damaging_csq =~  /^inframe_(inser|dele)tion$/){
        if ($allele_cadd >= 10){
            $s_name = "DamagingMissense";
        }else{
            $s_name = "BenignMissense";
        }
    }else{
        $s_name = "Other";
    }
    $variant_counts{$s_name}++;
    unshift @row, $variant_counts{$s_name};
    
    #add sample info to new array of array refs to be written 
    # in split cells alongside global variant fields
    my %samp_gts = VcfReader::getSampleActualGenotypes
        (
              line => $l, 
              all => 1,
              sample_to_columns => \%samples_to_columns,
              minGQ => $min_gq,
        );
    my %samp_gqs = VcfReader::getSampleGenotypeField
        (
              line => $l, 
              field => "GQ",
              all => 1,
              sample_to_columns => \%samples_to_columns,
              minGQ => $min_gq,
        );
    my %samp_ads = VcfReader::getSampleGenotypeField
         (
              line => $l,
              field => "AD",
              all => 1,
              sample_to_columns => \%samples_to_columns,
              minGQ => $min_gq,
        );
    
    foreach my $s (
            sort {$samples_to_columns{$a} <=> $samples_to_columns{$b}} 
            keys %samples_to_columns
    ){
        my @alts = split(/[\/\|]/, $samp_gts{$s});
        if (grep { $_ eq $var->{ORIGINAL_ALT} } @alts ){ 
            my @ads = split(",", $samp_ads{$s}); 
            my $depth = sum(@ads);
            if ($opts{d}){
                next if $opts{d} > $depth;
            }   
            
            my $ab = 0;
            if ( $depth > 0){
                $ab = $ads[$var->{ALT_INDEX}]/$depth;
            }
            if (@{$opts{b}}){
                next if $ab < $opts{b}->[0];
                if (scalar @{$opts{b}} > 1){
                    next if $ab > $opts{b}->[1];
                }
            }
            push @split_cells, [$s, $samp_gts{$s}, $samp_ads{$s}, $ab, $samp_gqs{$s}];
            my $var_class = $s_name; 
            if ($s_name eq 'HGMD'){
                if (grep {$_ eq 'DM'} @{$hgmd{variant_class}}){
                    $var_class = "HGMD_DM";
                }else{
                    $var_class = "HGMD_other";
                }
            }
            {
                no warnings 'uninitialized';
                push @{$sample_vars{$s}->{$var_class}}, join("|", @row);
            }
        }
    }
    $xl_obj->writeLine
    (
        line       => \@row, 
        worksheet  => $sheets{$s_name},
        succeeding => \@split_cells,
    );
}
###########################################################
#kept for legacy in case need to switch back to ClinVar VCF
sub getClinVarCodeVcf{
    my ($clinvars, $var) = @_;
    #array ref of clinvar VCF lines and a minimized variant hash from input
    my @annots = ();
    foreach my $cl (@$clinvars){
        my @match = split("\t", $cl);
        my @alts = split(",", VcfReader::getVariantField(\@match, "ALT"));
        my @clnsigs = split(",", VcfReader::getVariantInfoField(\@match, "CLNSIG"));
        my @codes = ();
        if (@alts == @clnsigs){#CLNSIG codes are a bit sketchy, sometimes one per allele sometimes not
            if (my $i = alleleMatches($var, $cl)){
                push @codes, split(/\|/, $clnsigs[$i-1]);
            }
        }else{
            foreach my $c (@clnsigs){
                push @codes, split(/\|/, $c);
            }
        }
        foreach my $code (@codes){
            push @annots, "$code ($clnsig_codes{$code})";
        }
    }
    return join(", ", @annots); 
}

###########################################################
sub getClinSig{
    my ($clinvars, $var) = @_;
    #array ref of clinvar lines and a minimized variant hash from input
    my @clnsig = ();
    my @trait  = ();
    my $isPathogenic = 0;
    foreach my $cl (@$clinvars){
        my @match = split("\t", $cl);
        push @clnsig, $match[$clinvar_sargs{columns}->{clinical_significance}] ;
        push @trait, $match[$clinvar_sargs{columns}->{all_traits}] ;
        $isPathogenic +=  $match[$clinvar_sargs{columns}->{pathogenic}] ;
    }
    return join(", ", @clnsig), join(", ", @trait), $isPathogenic; 
}
    
###########################################################
#kept for legacy in case need to switch back to ClinVar VCF
sub searchClinVarVcf{
    return if not $opts{c};
    my $var = shift;#$var is a ref to a single entrey from minimized alleles hash
    #simplify alleles and check if there's a match in HGMD file
    my @matches = ();
    #below is a hack because ClinVar file should be tsv.gz not VCF
    my @hits = VcfReader::searchByRegion(
        %clinvar_sargs,
        chrom => $var->{CHROM},
        start => $var->{POS},
        end   => $var->{POS} + length($var->{REF}) - 1,
    );
    foreach my $h (@hits){
        if (alleleMatchesClinVar($var, $h)){
            push @matches, $h;
        }
    }
    return @matches;
}

###########################################################

sub searchClinVar{
    return if not $opts{c};
    my $var = shift;#$var is a ref to a single entrey from minimized alleles hash
    #simplify alleles and check if there's a match in HGMD file
    my @matches = ();
    #below is a hack because ClinVar file should be tsv.gz not VCF
    my @hits = VcfReader::searchByRegion(
        %clinvar_sargs,
        chrom => $var->{CHROM},
        start => $var->{POS},
        end   => $var->{POS} + length($var->{REF}) - 1,
    );
    foreach my $h (@hits){
        if (alleleMatchesClinVar($var, $h)){
            push @matches, $h;
        }
    }
    return @matches;
}

###########################################################
sub searchHgmd{
    my $var = shift;#$var is a ref to a single entrey from minimized alleles hash
    #simplify alleles and check if there's a match in HGMD file
    my @matches = ();
    my @hits = VcfReader::searchByRegion(
        %search_args,
        chrom => $var->{CHROM},
        start => $var->{POS},
        end   => $var->{POS} + length($var->{REF}) - 1,
    );
    foreach my $h (@hits){
        if (alleleMatches($var, $h)){
            push @matches, $h;
        }
    }
    return @matches;
}


##########################################################
sub alleleMatches{
    my ($var, $line) = @_;
    #$var is a ref to a single entrey from minimized alleles hash
    #$line is a VCF entry
    #returns allele code for matching alt
    my @split = split ("\t", $line);
    my $chrom = VcfReader::getVariantField(\@split, 'CHROM');
    next if $chrom ne $var->{CHROM}; 
    my %l_min= VcfReader::minimizeAlleles(\@split);
    foreach my $k (keys %l_min){
        next if $l_min{$k}->{POS} != $var->{POS};
        next if $l_min{$k}->{REF} ne $var->{REF};
        next if $l_min{$k}->{ALT} ne $var->{ALT};
        return $k;
    }
    return 0;
}
##########################################################
sub alleleMatchesClinVar{
    my ($var, $line) = @_;
    #$var is a ref to a single entrey from minimized alleles hash
    #$line is a clinvar entry
    #returns 1 if it matches
    my @split = split ("\t", $line);
    my $chrom = $split[$clinvar_sargs{columns}->{chrom}] ;
    next if $chrom ne $var->{CHROM}; 
    my $pos = $split[$clinvar_sargs{columns}->{pos}] ;
    #alleles should already be minimized
    my $ref = $split[$clinvar_sargs{columns}->{ref}] ;
    my $alt = $split[$clinvar_sargs{columns}->{alt}] ;
    return 0 if $pos != $var->{POS};
    return 0 if $ref ne $var->{REF};
    return 0 if $alt ne $var->{ALT};
    return 1;
}
###########################################################
sub setupOutput{
    my @suffixes = (".vcf", ".txt");
    if (not $opts{o}){
        my ($out, $dir, $extension) = fileparse($opts{i}, @suffixes);
        $opts{o} = "$dir/$out.xlsx";
    }else{
        $opts{o} .= ".xlsx" if $opts{o} !~ /\.xlsx$/;
    }
    $xl_obj = TextToExcel->new( file=> $opts{o}, name => "HGMD");
    $sheets{HGMD} = 0;
    $header_formatting = $xl_obj->createFormat(bold => 1);
    $url_format = $xl_obj->createFormat(
        color     => 'blue',
        underline => 1,
    );  
    writeHeader($sheets{HGMD}, $headers{HGMD});

    $sheets{ClinVarPathogenic} = addSheet("ClinVarPathogenic", $headers{clinvar});

    $sheets{LOF} = addSheet("LOF", $headers{LOF});
    
    $sheets{DamagingMissense} = addSheet("DamagingMissense", $headers{missense});

    $sheets{BenignMissense} = addSheet("BenignMissense", $headers{missense});

    $sheets{Other} = addSheet("Other", $headers{LOF});
    
    $sheets{SampleSummary} = addSheet("Sample Summary", $headers{sample});
}

###########################################################
sub getHeaders{
    my %h = ();
    @{$h{HGMD}} =  ( 
        qw / 
            index
            Hgmd_ID
            variant_class
            Disease
            HGMD_Symbol
            HGVS
            ClinVarSig
			ClinVarTrait
            Chrom
            Pos
            Ref
            Alt
            Symbol
            Feature
            Consequence 
            HGVSc 
            HGVSp 
            GERP
            Polyphen
            SIFT
            CADD
            Sample
            GT
            AD
            AB
            GQ
     /);
    @{$h{LOF}} =  ( 
        qw / 
            index
            ClinVarSig
			ClinVarTrait
            Chrom
            Pos
            Ref
            Alt
            Symbol
            Feature 
            Consequence 
            HGVSc 
            HGVSp 
            GERP
            CADD
            Sample
            GT
            AD
            AB
            GQ
     /);
    
    @{$h{clinvar}} =  ( 
        qw / 
            index
            ClinVarSig
			ClinVarTrait
            Conflicted
            Chrom
            Pos
            Ref
            Alt
            Symbol
            Feature 
            Consequence 
            HGVSc 
            HGVSp 
            GERP
            Polyphen
            SIFT
            CADD
            Sample
            GT
            AD
            AB
            GQ
     /);



    @{$h{missense}} =  ( 
        qw / 
            index
            ClinVarSig
			ClinVarTrait
            Chrom
            Pos
            Ref
            Alt
            Symbol
            Feature 
            Consequence 
            HGVSc 
            HGVSp 
            GERP
            Polyphen
            SIFT
            CADD
            Sample
            GT
            AD
            AB
            GQ
     /);

    @{$h{sample}} =  (  
        qw / 
            Sample
            HGMD_DM
            ClinVarPathogenic
            HGMD_other
            LOF
            DamagingMissense
            BenignMissense
            Other
            /
    );
    return %h;
}
###########################################################
sub writeHeader{
    my $sheet = shift; 
    my $header = shift;
    $xl_obj->writeLine(line => $header, format => $header_formatting, worksheet => $sheet);
}

###########################################################
sub addSheet{
    my ($name, $header) = @_; 
    my $sheet = $xl_obj->addWorksheet($name);
    writeHeader($sheet, $header);
    return $sheet;
}

###########################################################
sub rankTranscriptsAndConsequences{
    my $csq_array = shift;#ref to array of VEP consequences 
    @$csq_array = rankConsequences(@$csq_array); 
    my $most_damaging = $csq_array-> [0] -> {consequence} ;
    @$csq_array = rankTranscripts(@$csq_array); 
    return first { $_ -> {consequence} eq $most_damaging } @$csq_array;
}
###########################################################
sub rankTranscripts{
    my @vars = @_;#array of hashes of VEP consequences
    return sort { getTranscriptsRanks( uc($a -> {feature}) ) <=> getTranscriptsRanks( uc($b -> {feature}) ) } @vars;
}
###########################################################
sub getTranscriptsRanks{
    my $symbol = shift;
    my $transcript = shift; 
    return -1 if not exists $refseq_ranks{$symbol};
    my $idx = first { $refseq_ranks{$symbol}->[$_] eq $transcript } 0 ..  $#{$refseq_ranks{$symbol}}; #also returns undef if transcript not found
    return -1 if not defined $idx;
    return $idx;
}
###########################################################
sub setTranscriptsRanks{
    return if ( not $opts{t} );
    open (my $TR, $opts{t}) or die "Can't open --transcripts file ($opts{t}) for reading: $!\n";
    my $header = <$TR>;
    chomp (my @head = split("\t", $header)); 
    my $n = 0;
    my %tr_columns = map { lc($_) => $n++ } map { (my $trim = $_) =~ s/^#+//; $trim } @head;
    foreach my $req ( qw / symbol refseq / ){
        if (not exists  $tr_columns{$req}){
            die "Could not find required column '$req' in header of --transcripts file $opts{t}. Found the following columns:\n" . join("\n", @head) . "\n";
        }
    }

    while (my $line = <$TR>){
        my @split = split("\t", $line); 
        my $symbol = uc($split[ $tr_columns{symbol} ]);
        my $transcript = uc ($split[ $tr_columns{refseq} ]); 
        push @{$refseq_ranks{$symbol} }, $transcript; 
    }
}


###########################################################
sub rankConsequences{
    my @vars = @_;#array of hashes of VEP consequences
    return sort { 
        $so_ranks{getMostDamagingConsequence($a)} <=> 
        $so_ranks{getMostDamagingConsequence($b)} 
    } @vars;
    
}

###########################################################
sub getMostDamagingConsequence{
    my $csq = shift;#hash ref to VEP consequences for single transcript/feature
    my @s_csq = split("&", $csq->{consequence} );
    @s_csq = sort { $so_ranks{$a} <=> $so_ranks{$b} } @s_csq;
    return $s_csq[0];
}
    
###########################################################
sub setConsequenceRanks{
    my @so_terms = qw /
        transcript_ablation  
        splice_acceptor_variant
        splice_donor_variant
        stop_gained
        frameshift_variant
        stop_lost
        start_lost
        transcript_amplification
        inframe_insertion
        inframe_deletion
        missense_variant
        protein_altering_variant
        splice_region_variant
        incomplete_terminal_codon_variant
        stop_retained_variant
        synonymous_variant
        coding_sequence_variant
        mature_miRNA_variant
        5_prime_UTR_variant
        3_prime_UTR_variant
        non_coding_transcript_exon_variant
        intron_variant
        NMD_transcript_variant
        non_coding_transcript_variant
        upstream_gene_variant
        downstream_gene_variant
        TFBS_ablation
        TFBS_amplification
        TF_binding_site_variant
        regulatory_region_ablation
        regulatory_region_amplification
        feature_elongation
        regulatory_region_variant
        feature_truncation
        intergenic_variant 
    /;
    my $n = 0;
    %so_ranks = map { $_ => $n++ } @so_terms; 
}

###########################################################
sub getClnSigCodes{
    return (
        0 => "Uncertain significance",
        1 => "not provided",
        2 => "Benign",
        3 => "Likely benign",
        4 => "Likely pathogenic",
        5 => "Pathogenic",
        6 => "drug response",
        7 => "histocompatibility",
        255 => "other"
    );
}


