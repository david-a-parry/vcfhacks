#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Term::ProgressBar;
use Data::Dumper;
use POSIX qw/strftime/;
use FindBin;
use lib "$FindBin::Bin/lib";
use ParsePedfile;
use VcfReader;
use FindBiallelic;

my $data_dir = "$FindBin::Bin/data";
my $progressbar;
my @samples = ();
my @reject = ();
my @reject_except = ();
my @classes = (); 
my @add_classes = ();
my @damaging = ();

my %opts = (
    s           => \@samples,
    r           => \@reject,
    x           => \@reject_except,
    classes     => \@classes,
    add_classes => \@add_classes,
    d           => \@damaging,
);


GetOptions(
    \%opts,
    'i|input=s' ,
    'o|output=s',
    'l|list:s', 
    's|samples=s{,}'           => \@samples,
    'r|reject=s{,}'            => \@reject,
    'x|reject_all_except:s{,}' => \@reject_except,
    'f|family=s', # ped file
    'm|mode=s',
    'classes=s{,}'             => \@classes,
    'add_classes=s{,}'         => \@add_classes,
    'd|damaging=s{,}',         => \@damaging,
    'canonical_only',
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

#check VEP/SNPEFF header and get annotation order
my %csq_header = getAndCheckCsqHeader();
  # hash of functional annotations from SnpEff/VEP to their annotation index
    
#set default consequence fields to retrieve 
my @csq_fields = getCsqFields();#consequence fields to retrieve from VCF

#and check variant consequence classes are acceptable
getAndCheckClasses();

#check in silico prediction classes/scores are acceptable
my %in_silico_filters = getAndCheckInSilicoPred();
  #hash of prediction program names and values to filter

#check ped, if provided is ok and get specified families
my ($ped_obj, @fams) = checkPedAndFamilies();

#check affected and unaffected samples
addSamplesFromPed();



#at least one affected sample


#samples not duplicated between affected and unaffected

#min matching settings are not greater than no. affecteds

#min matching settings are not greater than no. families if using ped

######PRE-PROCESSING######

#Open output

#write meta header

#add program run parameters

#add INFO fields 

#start progress bar 

######READ VARIANTS AND PROCESS######

#read line

#check if new chromosome - if so go to checking biallelics and clear collected data

#skip if not PASS and PASS required

#skip if no variant allele in affecteds 

#skip if not identical GTs in all affecteds and identical variants required

#skip if homozygous in any unaffected

#filter if AF annotations > MAF

#skip if VEP/SNPEFF annotation is not the correct class (or not predicted damaging if using that option)

#create variant hash if passed all of the above 

#add variant hash to hash of transcript ID to arrays of variant hashes 

######CHECK BIALLELIC######

#For each transcript (i.e. each key of our hash of transcript IDs to arrays of var hashes)...

#...get combinations of biallelic alleles for each affected

#...get combinations of biallelic alleles for each unaffected

#...remove biallelic alleles present in unaffecteds from putative biallelic alleles in affecteds

#...if using ped files check segregation of putative biallelic combinations for each family and remove those that do not segregate properly

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
    my %seen = ();
    @classes = grep { ! $seen{$_}++ } @classes;
    foreach my $class (@classes) {
        die "Error - variant class '$class' not recognised.\n"
          if not exists $all_classes{lc($class)} ;
    }
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
        die "Not enough fields in classes file line: $line\n" if @split < 2;
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
    if ($opts{m} eq 'vep'){
        return 
        qw(
            allele
            gene
            feature
            feature_type
            consequence
            hgnc
        );
    }else{
        return 
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
#################################################
#################################################
#################################################
#################################################
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
