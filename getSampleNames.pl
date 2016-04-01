#!/usr/bin/env perl

=head1 NAME

getSampleNames.pl - print the sample IDs in one or more VCF files


=head1 SYNOPSIS
 
getSampleNames.pl -i [vars.vcf] [options]

=cut

=head1 ARGUMENTS

=over 8

=item B<-i    --input>

One or more VCF files to get sample names from. Sample names and counts will be given for each file separately.

=item B<-l    --list_only>

Output the sample names only, ommiting the number of samples.

=item B<-c    --count_only>

Output the number of samples only, do not output sample names.

=item B<-d    --delimiter>

By default sample names are separated by commas and spaces. Specify a different delimiter with this option.

=item B<-n    --new_lines>

Print each sample name on a new line instead of separating them with commas and spaces (i.e. like specifying a new line to the --delimiter option)

=item B<-s    --sort>

Sort sample names in alphabetical order (by default samples are given in the order they appear in the VCF, which is usually alphabetical anyway).

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


=head1 COPYRIGHT AND LICENSE

Copyright 2014,2015 David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($RealBin);
use lib "$RealBin/lib";
use ParseVCF;

my @vcfs = ();
my %opts = 
(
    input       => \@vcfs,
    delimiter   => ", ",
);
GetOptions(
    \%opts,
    'input=s{,}' => \@vcfs,
    'delimiter=s',
    'count_only',
    'list_only',
    'new_lines',
    'sort',
    'help',
    'manual',
) or pod2usage (-message => 'Syntax error.', -exitval => 2);
pod2usage( -verbose => 2 ) if ($opts{manual});
pod2usage( -verbose => 1 ) if ($opts{help});
pod2usage (-message => '--input argument is required. For help run the --help or --manual option', -exitval => 2) if not @vcfs;
if ($opts{count_only} and $opts{list_only}){
    pod2usage 
    (
        -message => 
        '--list_only and --count_only are mutually exculsive. For help run the --help or --manual option', -exitval => 2
    );
}
$opts{delimiter} = "\n" if $opts{new_lines};
foreach my $vcf (@vcfs){ 
    my $obj = ParseVCF->new(file => $vcf, noLineCount => 1);
    my @samples = ();
    eval{@samples = $obj->getSampleNames();};
    unless ($opts{list_only}){
        print "$vcf has " . scalar @samples . " sample";
        print "s" if @samples != 1;
        my $spacer = "\t";
        $spacer = "\n" if $opts{new_lines};
        print ":$spacer" unless $opts{count_only};
    }
    if (not $opts{count_only}){
        @samples = sort @samples if $opts{sort};
        print join($opts{delimiter}, @samples)  if @samples;
    }
    print "\n";
}
