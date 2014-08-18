#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Term::ProgressBar;
use FindBin;
use lib "$FindBin::Bin";
use VcfReader;
my %opts = ();
GetOptions(\%opts,
        "output=s",
        "input=s",
        "help",
        "manual",
        ) or pod2usage(-exitval => 2, -message => "Syntax error") ;
pod2usage (-verbose => 2) if $opts{manual};
pod2usage (-verbose => 1) if $opts{help};
pod2usage(-exitval => 2, -message => "Syntax error") if not $opts{input}; 

if ($opts{output}){
    VcfReader::sortVcf(vcf => $opts{input}, output => $opts{output});
}else{
    VcfReader::sortVcf(vcf => $opts{input});
}


=head1 NAME

sortVcf.pl - Sort a VCF file in coordinate order.

=head1 SYNOPSIS

        sortVcf.pl -i [vcf file] -o [output]
        sortVcf.pl -h (display help message)
        sortVcf.pl -m (display manual page)

=cut

=head1 ARGUMENTS

=over 8

=item B<-i    --input>

VCF file input.

=item B<-o    --output>

Output file. Optional - default is STDOUT.

=item B<-h    --help>

Display help message.

=item B<-m    --manual>

Show manual page.

=back 

=cut


=head1 DESCRIPTION

This program will attempt to identify the contig order from a VCF file header and sort all variants in contig and coordinate order. If contigs are not present in your VCF header (i.e. lines beginning '##contig=<ID=') variants will be sorted numerically, followed by 'X', 'Y', 'MT' and then ascibetically

Although VCFs usually are generated already in coordinate sorted order, occasionally programs such as the VEP output the odd line out of order so this program is designed as a quick fix.

=cut

=head1 AUTHOR

David A. Parry

University of Leeds

=head1 COPYRIGHT AND LICENSE

Copyright 2013  David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut

