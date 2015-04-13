#!/usr/bin/perl

=head1 NAME

getSampleNames - print the sample IDs in one or more VCF files


=head1 SYNOPSIS
 
getSampleNames -i [vars.vcf] [options]

=cut

=head1 ARGUMENTS

=over 8

=item B<-i    --input>

One or more VCF files to get sample names from. Sample names and counts will be given for each file separately.

=item B<-s    --sort>

Sort sample names in alphabetical order (normally samples are given in the order they appear in the VCF).

=item B<-c    --count_only>

Output the number of samples only, do not output sample names.

=item B<-n    --new_lines>

Print each sample name on a new line instead of separating them with commas and spaces.

=item B<-h    --help>

Show this script's help information.

=item B<-m    --manual>

Show this script's manual page.

=back

=head1 DESCRIPTION

Reports the number of samples and sample names for given VCF files.

=cut

=head1 AUTHOR

David A. Parry

University of Leeds


=head1 COPYRIGHT AND LICENSE

Copyright 2014 David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin/lib";
use ParseVCF;

my @vcfs = ();
my $delimiter = ", ";
my $use_newlines;
my $count_only;
my $sort;
my $help;
my $manual;
GetOptions(
    'input=s{,}' => \@vcfs,
    'count_only' => \$count_only,
    'new_lines' => \$use_newlines,
    'sort' => \$sort,
    'help' => \$help,
    'manual' => \$manual,
) or pod2usage (-message => 'Syntax error.', -exitval => 2);
pod2usage( -verbose => 2 ) if ($manual);
pod2usage( -verbose => 1 ) if ($help);
pod2usage (-message => '--input argument is required. For help run the --help or --manual option', -exitval => 2) if not @vcfs;
$delimiter = "\n" if $use_newlines;
foreach my $vcf (@vcfs){ 
    my $obj = ParseVCF->new(file => $vcf, noLineCount => 1);
    my @samples = ();
    eval{@samples = $obj->getSampleNames();};
    print "$vcf has " . scalar @samples . " sample";
    print "s" if @samples != 1;
    if (not $count_only){
        my $spacer = "\t";
        $spacer = "\n" if $use_newlines;
        @samples = sort @samples if $sort;
        print ":$spacer" . join($delimiter, @samples)  if @samples;
    }
    print "\n";
}
