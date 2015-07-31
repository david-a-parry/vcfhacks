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

##########################################################
       
=head1 NAME

checkInheritance.pl - output simple stats about variants and expected inheritance patterns

=head1 SYNOPSIS

    checkInheritance.pl -i <variants.vcf> -f family.ped [options]
    checkInheritance.pl --help (show help message)
    checkInheritance.pl --manual (show manual page)

=cut 

=head1 ARGUMENTS

=over 8 

=item B<-i    --input>

VCF file input.

=item B<-o    --output>

Prefix for output files. By default the following output files are produced for each child identified in the ped file:

    <family_id>-<child_id>-mendelian.txt
    <family_id>-<child_id>-non-mendelian.txt
    
A single summary file (summary.txt) is also produced. If a prefix is specified here, it will be prepended to each of these output files.

=item B<-f    --family>

A PED file (format described below) containing information about samples in the VCF. At least one parent-child relationship must be included. All parent child relationships identified in the file will be examined to assess expected inheritance patterns.

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


=item B<-p    --pass_filters>

Only count variants with a value of PASS in its FILTER field.

=item B<-q    --quality>

Minimum variant quality score to consider. Variants with QUAL scores below this value will be skipped.

=item B<-g    --genotype_quality>

Minimum genotype quality score to consider. Variants will be skipped if any of the parent-child trios/duos have a genotype quality score below this value.

=item B<-b    --progress>

Show a progress bar while running.

=item B<-h    ---help>

Show the program's help message.

=item B<--manual>

Show the program's manual page.

=back

=cut

=head1 DESCRIPTION
 
This program reads a VCF file and associated PED file in order to assess the amount of variants following expected inheritance patterns between parents and children. For each parent-child relationship identified the following files are produced:
  
    <family_id>-<child_id>-mendelian.txt
    <family_id>-<child_id>-non-mendelian.txt
    
These files each produce a tab delimited file with the following columns:
    
    Chr Pos Ref Alt Qual Child-GQ Mat-GQ Pat-GQ Child-GT Mat-GT Pat-GT

Each row gives details for a variant following either expected Mendelian inheritance patterns or not following expected inheritance patterns for the respective files. This may be useful to identify sample mixups in conjunction with tools such as the relatedness tools available from vcftools (https://vcftools.github.io/examples.html).

A summary file (summary.txt) is also produced which summarises for each child the number of Mendelian inherited variants, the number of non-mendelian inherited variants, number of informative variants and number of uninformative variants.

 
=cut

=head1 AUTHOR

David A. Parry

University of Leeds

=head1 COPYRIGHT AND LICENSE

Copyright 2015  David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut

