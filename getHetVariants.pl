#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/lib";
use ParseVCF;
use Pod::Usage;
use Getopt::Long;
my $vcf;
my $out;
my @samples = ();
my $min_gq  = 0;
my @ab = ();
my $hom;
my $help;
my $man;
GetOptions(
    "input=s"      => \$vcf,
    "output=s"     => \$out,
    "reverse"      => \$hom,
    "samples=s{,}" => \@samples,
    "gq=f"         => \$min_gq,
    "ab=f{,}"      => \@ab,
    "help|?"       => \$help,
    "manual"       => \$man
) or pod2usage( -exitval => 2, -message => "Syntax error" );
pod2usage( -verbose => 2 ) if $man;
pod2usage( -verbose => 1 ) if $help;
pod2usage( -exitval => 2, -message => "Syntax error" ) if not $vcf;

foreach my $bal (@ab){
    die "--ab argument must be between 0 and 1\n" if $bal > 1 or $bal < 0;
}

my $OUT;

if ($out) {
    open( $OUT, ">$out" ) || die "Can't open $out for writing: $!\n";
}
else {
    $OUT = \*STDOUT;
}

my $obj = ParseVCF->new( file => $vcf );

if ( not @samples ) {
    @samples = $obj->getSampleNames;
}

print $OUT join( "", @{ $obj->get_metaHeader } );
print $OUT $obj->get_header . "\n";
while ( my $line = $obj->readLine ) {
    if ($hom) {
        my $hom =
          $obj->sampleIsHomozygous( multiple => \@samples, minGQ => $min_gq );
        print $OUT "$line\n" if $hom;
    }
    else {
        my $het;
        if (@ab){
            $het = $obj->sampleIsHeterozygous(allele_balance => \@ab, multiple => \@samples, minGQ => $min_gq );
        }else{
            $het = $obj->sampleIsHeterozygous(multiple => \@samples, minGQ => $min_gq );
        }
        print $OUT "$line\n" if $het;
    }
}

=head1 NAME

getHetVariants.pl - print heterozygous variants from VCF file


=head1 SYNOPSIS

getHetVariants.pl -i [vcf_file] [options]

=head1 ARGUMENTS

=over 8

=item B<-i    --input>

Input VCF file.

=item B<-o    --output>

Output filename.

=item B<-s    --samples>

One or more samples to check variants for. Default is to check all samples in vcf.

=item B<-g    --gq>

Minimum genotype quality to specify for het variants.  Anything less will be counted as a no call.

=item B<-r    --reverse>

Use this flag to print homozygous variants rather than heterozygous variants.

=item B<-h    --help>

Show help message.

=item B<-m    --manual>

Show manual page.

=back

=cut

=head1 DESCRIPTION

Prints only heterozygous variants from a VCF file.  Will only print if all variants or all variants for samples selected with --samples argument are heterozygous. More sophisticated implementations may follow.


=head1 AUTHOR

David A. Parry

University of Leeds

=head1 COPYRIGHT AND LICENSE

AUTHOR

       David A. Parry

       University of Leeds


COPYRIGHT AND LICENSE

       Copyright 2013  David A. Parry


This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut
