#!/usr/bin/perl
use warnings;
use strict; 
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use Term::ProgressBar;
use FindBin;
use lib "$FindBin::Bin/lib";
use ParseVCF;

my $input;
my $out;
my $find_only;
my $help;
my $manual;

GetOptions('help' => \$help, 'manual' => \$manual, 'input=s{,}' => \$input, 'output=s' => \$out, 'find' => \$find_only,
) or pod2usage(-exitval => 2, -message => "Syntax error") ;

pod2usage (-verbose => 2) if $manual;
pod2usage (-verbose => 1) if $help;
pod2usage(-exitval => 2, -message => "--input argument is required") if not $input;
if ($find_only){
    print STDERR "Identifying and outputting unbalanced multiallelic variants only. Not splitting.\n";
}
my $OUT;
if ($out){
   open ($OUT,">$out") or die "Can't open $out for writing: $!\n"; 
}else{
    $OUT = \*STDOUT;
}

my $vcf_obj = ParseVCF->new(file=>$input);
print $OUT  $vcf_obj->getHeader(0) ."##split_multiallelic_variants.pl=\"--input $input\"\n";
print $OUT "##INFO=<ID=ParseVCF_split_variant-sample_fields_may_be_inaccurate,Number=0,Type=Flag,Description=\"Variant has been split from being a multiallelic variant containing both an SNV and an Indel into separate lines by ParseVCF.pm. Care should be taken when interpreting variant fields.\">\n";
print $OUT $vcf_obj->getHeader(1);
LINE:    while (my $line = $vcf_obj->readLine){
    my $ref = $vcf_obj->getVariantField("REF");
    my @alts = $vcf_obj->readAlleles(alt_alleles => 1);
    my $same_length = 0;
    my $diff_length = 0;
    foreach my $alt (@alts){
        if (length($alt) == length($ref)){
            $same_length++;
        }else{
            $diff_length++;
        }
    }
    if ($same_length and $diff_length){#SNV/MNV and Indel at same site
        if ($find_only){
            print $OUT "$line\n";
        }else{
            my @splits = $vcf_obj->splitMultiAllelicVariants();
            print $OUT join("\n", @splits) ."\n";
        }
    }else{
        print $OUT "$line\n" unless ($find_only);
    }
}



=head1 NAME

splitMultiAllelicVariants.pl - split indels and SNVs/MNVs present on the same line in a VCF

=head1 SYNOPSIS

        splitMultiAllelicVariants.pl -i [vcf file] [options]
        splitMultiAllelicVariants.pl -h (display help message)
        splitMultiAllelicVariants.pl -m (display manual page)

=cut

=head1 ARGUMENTS

=over 8

=item B<-i    --input>

VCF file input.

=item B<-o    --output>

Output file. Optional - default is STDOUT.

=item B<-f    --find>

Use this flag to find and print unbalanced variants that would normally be split. Only these (unedited) variants will be printed.

=item B<-h    --help>

Display help message.

=item B<-m    --manual>

Show manual page.

=back 

=cut


=head1 DESCRIPTION

This program will identify lines of a VCF containing both SNV/MNVs and Indel variants and write these onto different lines as individual variants.  For example:

1   2587458 .   C    CTA,A 104.49   PASS

will be split into one SNV and one insertion variant.

1   10845114    rs111606813 GTC TCC,G   390.10  PASS

will be split into one MNV and one deletion. 

The motivation for this was due to a bug in the variant_effect_predictor.pl (present at v73 and before) that mishandled SNVs on the same line as 'unbalanced' variants, bearing in mind the GATK HaplotypeCaller does not offer separate SNP and INDEL genotyping nor does the GATK offer a way to split such variants from each other when present on the same line. Please note that this script is a BIG HACK and that certain more complicated VCF fields will not be accurate following use of this script, but sample genotypes should be as you would expect. It's intended only for use before using VEP or to rescue variants that GATKs SelectVariants tool may ignore. Hopefully the GATK people will come up with a better fix for this in time.

(splitMultiallelicVariants is a bit of a misnomer, a better suggestion would be appreciated!).


=cut

=head1 AUTHOR

David A. Parry

=head1 COPYRIGHT AND LICENSE

Copyright 2013, 2014  David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut

