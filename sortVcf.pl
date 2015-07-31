#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Term::ProgressBar;
use FindBin;
use lib "$FindBin::Bin/lib";
use VcfReader;
my %opts = ();
GetOptions(\%opts,
        "output=s",
        "input=s",
        "dict=s",
        "help",
        "manual",
        ) or pod2usage(-exitval => 2, -message => "Syntax error") ;
pod2usage (-verbose => 2) if $opts{manual};
pod2usage (-verbose => 1) if $opts{help};
pod2usage(-exitval => 2, -message => "Syntax error") if not $opts{input}; 

my %sort_args = (vcf => $opts{input});
if ($opts{output}){
    $sort_args{output} = $opts{output};
}
if($opts{dict}){
    my ($contigs, $id_head) = readDictFile($opts{dict});
    $sort_args{contig_order} = $contigs;
    $sort_args{dict} = $id_head;
}

VcfReader::sortVcf(%sort_args);

##################################################
sub readDictFile{
    my ($dict) = @_;
    open (my $DICT, $dict) or die "Can't open DICT file $dict for reading: $!\n";
    my @seqs = ();
    my @dict = ();
    while (<$DICT>){
        chomp;
        if (/^\@SQ\s/){
            if (/\sSN:(\S+)/){
                my $seq = $1;
                my $length;
                if (/\sLN:(\d+)/){
                    $length = $1;
                }
                push @seqs, $seq;
                my $head_line = "##contig=<ID=$seq";
                $head_line .= ",length=$length" if $length;
                $head_line .= ">";
                push @dict, $head_line;
            }else{
                die "Couldn't find sequence name (SN) field for \@SQ line:\n$_\n";
            }
        }
    }
    close $DICT;
    die "No contigs found in $dict!\n" if not @seqs;
    print STDERR "Found " . scalar(@seqs) . " contigs in $dict.\n";
    my $n = 0;
    my %chroms = map {$_ => $n++} @seqs;
    return (\%chroms, \@dict);
}
##################################################


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

=item B<-d    --dict>

An optional sequence dictionary (.dict) file to specify the sort order of contigs. Sequences should be specified in this file by whitespace delimited lines beginning @SQ and containing at a minimum 'SN:[sequence name]'.

=item B<-h    --help>

Display help message.

=item B<-m    --manual>

Show manual page.

=back 

=cut


=head1 DESCRIPTION

This program will attempt to identify the contig order from a VCF file header and sort all variants in contig and coordinate order. If contigs are not present in your VCF header (i.e. lines beginning '##contig=<ID=') and a sequence dictionary is not provided variants will be sorted numerically, followed by 'X', 'Y', 'MT' and then ascibetically.

Although VCFs usually are generated already in coordinate sorted order, occasionally programs such as the VEP output the odd line out of order so this program is designed as a quick fix.

=cut

=head1 AUTHOR

David A. Parry

=head1 COPYRIGHT AND LICENSE

Copyright 2013, 2014, 2015  David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut

