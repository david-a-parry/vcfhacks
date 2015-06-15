#!/usr/bin/perl
#David Parry June 2015
use strict;
use warnings;
use FileHandle;
use Getopt::Long;
use Pod::Usage;
use Term::ProgressBar;
use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/lib";
use ParsePedfile;
use VcfReader;

my %opts = ();
GetOptions(
    \%opts,
    'i|input=s', 
    'o|output=s',
    'f|family=s', #.ped file
    'p|pass_filters',
    'q|quality=i',
    'g|genotype_quality=i',
    'b|progress',
    'h|?|help' ,
    'm|manual'
) or pod2usage( -message => "Syntax error", exitval => 2 );

pod2usage( -verbose => 2 ) if $opts{m};
pod2usage( -verbose => 1 ) if $opts{h};
pod2usage( -message => "Syntax error: --input is required.", exitval => 2 )
  if not $opts{i};
pod2usage( -message => "Syntax error: a pedigree file must be passed to the --family argument.", exitval => 2 )
  if not $opts{f};
$opts{q} = 0 if not defined $opts{q};
$opts{g} = 0 if not defined $opts{g};
my %children = ();
print STDERR "Parsing pedigree file...\n";
my $ped = ParsePedfile->new( file => $opts{f});
foreach my $sample ($ped->getAllSamples()){
    my @parents = $ped->getParents($sample);
    if (@parents == 2){
        push @{$children{$sample}}, @parents;
    }
}
if (not keys %children){
    die "ERROR: At least one parent-child relationship is required in pedigree file!\n";
}
print STDERR "Found ". scalar(keys %children) . " samples with two parents in ped file.\n";

my %samples_to_columns = VcfReader::getSamples(vcf => $opts{i}, get_columns => 1);
foreach my $c (sort keys %children){
    print STDERR "Checking child $c inheritance from ". join(" and ",  @{$children{$c}}) . ".\n";
    foreach my $s ($c, @{$children{$c}}){
        if (not exists $samples_to_columns{$s}){
            die "Sample $s does not exist in VCF. Samples are: " .join(", ", sort keys %samples_to_columns) . "\n";
        }
    }
} 
my %output = ();
foreach my $c (sort keys %children){
    my $fid = $ped->getFamilyId($c);
    foreach my $file (qw / mendelian non-mendelian / ){
        $output{$c}->{$file}->{file} = "$fid-$c-$file.txt";
        if ($opts{o}){
            $output{$c}->{$file}->{file} = $opts{o} ."_$output{$c}->{$file}->{file}";
        }
        my $fh = FileHandle->new("> $output{$c}->{$file}->{file}") or die 
                "Can't open ouptut file $output{$c}->{$file}->{file} for writing.\n";
        $output{$c}->{$file}->{fh} = $fh;
    }
}

writeHeaders();
my %counts = ();
foreach my $c (sort keys %children){
    $counts{$c}->{mendelian} = 0;
    $counts{$c}->{"non-mendelian"} = 0;
    $counts{$c}->{informative} = 0;
    $counts{$c}->{uninformative} = 0;
}
my $progressbar;
my $next_update = 0;
my $total_vcf = 0;
if ( defined $opts{b} ) {
    print STDERR "Counting variants and setting up progressbar...\n";
    $total_vcf = VcfReader::countVariants( $opts{i} );
    print STDERR "$opts{i} has $total_vcf variants. \n";
    $progressbar = Term::ProgressBar->new(
        { name => "Analyzing", 
          count => ($total_vcf), 
          ETA => "linear" 
        } 
    );
}
my $n = 0 ;
my $FH = VcfReader::_openFileHandle($opts{i}); 
while (my $line = <$FH>){
    next if $line =~ /^#/;
    $n++;
    if ($progressbar) {
        $next_update = $progressbar->update($n) if $n >= $next_update;
    }
    chomp($line); 
    my @split_line = split("\t", $line);
    my $chr = VcfReader::getVariantField(\@split_line, "CHROM");
    my $pos = VcfReader::getVariantField(\@split_line, "POS");
    my $qual = VcfReader::getVariantField(\@split_line, "QUAL");
    my $ref = VcfReader::getVariantField(\@split_line, "REF");
    my $alt = VcfReader::getVariantField(\@split_line, "ALT");
    next if (length($ref) > 1);#SNVs only
    next if (length($alt) > 1);#biallelic SNVs only
    next if ($qual < $opts{q}); 
    foreach my $c (sort keys %children){
        my $inh = checkInheritance(\@split_line, $c, $chr);
        #returns 0 if uninformative, 1 if mendelian inheritance
        # and -1 if non-mendelian inheritance
        if (not $inh){#uninformative
            $counts{$c}->{uninformative}++;
        }else{
            $counts{$c}->{informative}++;
            if ($inh > 0){
                $counts{$c}->{mendelian}++;
                writeVar(\@split_line, $c, "mendelian");
            }else{
                $counts{$c}->{"non-mendelian"}++;
                writeVar(\@split_line, $c, "non-mendelian");
            }
        }
    }
}
my $summary_file = "summary.txt";
if ($opts{o}){
    $summary_file = $opts{o} ."_$summary_file";
}
open (my $SUM, "> $summary_file") or die "Can't open $summary_file for writing.\n";
print $SUM join("\t", 
qw /
    Family 
    Sample 
    Mendelian 
    Non-Mendelian
    Informative
    Uninformative
/) . "\n";

foreach my $c (sort {$ped->getFamilyId($a) cmp $ped->getFamilyId($b)} keys %children){
    my $fid = $ped->getFamilyId($c);
    my @fields = ();
    push @fields, 
    (
        $fid, 
        $c, 
        $counts{$c}->{mendelian},
        $counts{$c}->{"non-mendelian"},
        $counts{$c}->{informative},
        $counts{$c}->{uninformative},
    );
    print $SUM join("\t", @fields) . "\n"; 
}
print STDERR "Done.\n"; 

###########################################################
sub writeVar{
    my ($var, $c, $file) = @_;
    my $f = $ped->getFather($c);
    my $m = $ped->getMother($c);
    foreach my $field ( qw /CHROM POS REF ALT QUAL/){
        $output{$c}->{$file}->{fh} -> print
            (
                VcfReader::getVariantField
                    (
                        $var, 
                        $field,
                    ) 
                ."\t"
            );
    }
    my @gqs = ();
    my @gts = ();
    foreach my $s ($c, $f, $m){ 
        push @gqs, VcfReader::getSampleGenotypeField
            (
                line => $var, 
                field => "GQ",
                sample => $s,
                sample_to_columns => \%samples_to_columns
            );
        push @gts, VcfReader::getSampleCall
            (
                line => $var, 
                sample => $s, 
                minGQ => $opts{g},
                sample_to_columns => \%samples_to_columns
            );
    }
    $output{$c}->{$file}->{fh} -> print
        (
            join("\t", @gqs, @gts)."\n"
        );
}
###########################################################
sub checkInheritance{
#returns 0 if uninformative, 1 if mendelian inheritance
    my ($var, $c, $chrom) = @_;
    my $f = $ped->getFather($c);
    my $m = $ped->getMother($c);
    my $gt = VcfReader::getSampleCall
     (
         line => $var, 
         sample => $c, 
         minGQ => $opts{g},
         sample_to_columns => \%samples_to_columns
     );
    return 0 if $gt =~ /\.[\/\|]\./;
    my $f_gt;
    if ($f){
        $f_gt = VcfReader::getSampleCall
         (
             line => $var, 
             sample => $f, 
             minGQ => $opts{g},
             sample_to_columns => \%samples_to_columns
         );
        return 0 if $f_gt =~ /\.[\/\|]\./;
    }
    my $m_gt;
    if ($m){
        $m_gt = VcfReader::getSampleCall
         (
             line => $var, 
             sample => $m, 
             minGQ => $opts{g},
             sample_to_columns => \%samples_to_columns
         );
        return 0 if $m_gt =~ /\.[\/\|]\./;
    }
    if ($gt eq $m_gt and $gt eq $f_gt){
        return 0;
    }
    $chrom =~ s/^chr//;
    if ($chrom eq 'Y'){
        if ($ped->getSex($c) == 1){#male
            if ($f_gt =~ /(\d)[\/\|]\g1/ and $gt =~ /(\d)[\/\|]\g1/){
            #only makes sense if homo(hemi)zygous on father's Y and child
                return 1 if $f_gt eq $gt;
                return -1;
            }else{
                return 0;
            }
        }else{
            return 0;
        }
    }elsif($chrom eq 'X' and $ped->getSex($c) == 1){#male X chrom
        if ($gt !~ /(\d)[\/\|]\g1/){
        #only makes sense if homo(hemi)zygous on male child's X
            return 0;
        }
        my %mat_alleles = map {$_ => undef} split(/[\/\|]/, $m_gt);
        my $c_al = (split /[\/\|]/, $m_gt)[0];
        foreach my $m_al (keys %mat_alleles){
            if ($c_al eq $m_al){
                return 1;
            }
        }
        return -1;
    }else{
        my %mat_alleles = map {$_ => undef} split(/[\/\|]/, $m_gt);
        my %pat_alleles = map {$_ => undef} split(/[\/\|]/, $f_gt);
        my @mat = keys %mat_alleles;
        my @pat = keys %pat_alleles;
        for (my $i = 0; $i < @pat; $i++){
            for (my $j = 0; $j < @mat; $j++){
                if ($gt =~ /^\Q$pat[$i]\E[\/\|]$mat[$j]$/){
                    return 1;
                }        
                if ($gt =~ /^\Q$mat[$j]\E[\/\|]$pat[$i]$/){
                    return 1;
                }
            }
        }
        return -1;
    }
}
###########################################################
sub writeHeaders{
    foreach my $c (sort keys %children){
        foreach my $file (qw / mendelian non-mendelian / ){
            $output{$c}->{$file}->{fh} -> print(
                    join("\t", qw /Chr Pos Ref Alt Qual/)); 
            $output{$c}->{$file}->{fh} -> print("\t$c-GQ");
            if ($ped->getFather($c)){
                $output{$c}->{$file}->{fh} -> print("\tPat-GQ");
            }
            if ($ped->getMother($c)){
                $output{$c}->{$file}->{fh} -> print("\tMat-GQ");
            }
            $output{$c}->{$file}->{fh} -> print("\t$c-GT");
            if ($ped->getFather($c)){
                $output{$c}->{$file}->{fh} -> print("\tPat-GT");
            }
            if ($ped->getMother($c)){
                $output{$c}->{$file}->{fh} -> print("\tMat-GT");
            }
            $output{$c}->{$file}->{fh} -> print("\n");
        }
    }
}
        

    
