#!/usr/bin/perl
#David Parry University of Leeds April 2011

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.If not, see <http://www.gnu.org/licenses/>.

=head1 NAME

getVariantsByLocation.pl - print variants from a vcf file that lie within a list of genomic regions


=head1 SYNOPSIS
 
getVariantsByLocation.pl -i [vars.vcf] -b [regions.bed] -o [output.vcf]

=cut

=head1 ARGUMENTS
   

=over 8

=item B<-i    --input>

VCF format file containing list of variants to filter (required). 

=item B<-b    --bed>

Bed file containing regions to filter on (required unless using --regions argument). Can also be a file containing a list of intervals, one interval per line.

=item B<-r    --regions>

Specify one or more regions to filter on seperated with spaces in format "chr1:1-2000". Can be used in conjunction with or instead of --bed argument.

=item B<-o    --output>

File to print output (optional). Will print to STDOUT by default.

=item B<-n    --no_header>

Use this flag to prevent outputting the VCF header.

=item B<-q    --quiet>

Use this flag to supress printing of information to STDERR.

=item B<-s    --silent>

Use this flag to supress warnings and information to STDERR.

=item B<-h    --help>

Show this script's help information.

=item B<-m    --manual>

Show this script's manual page.

=back

=cut

=head1 EXAMPLES

getVariantsByLocation.pl -i [vars.vcf] -b [regions.bed] 

getVariantsByLocation.pl -i [vars.vcf] -r chr1:2000000-50000000 -o [filtered.vcf]

   
=cut

=head1 DESCRIPTION

Takes a bed file or a list of regions and outputs variants from a vcf file that lie in these regions.


=cut

=head1 AUTHOR

David A. Parry

University of Leeds


=head1 COPYRIGHT AND LICENSE

Copyright 2014 David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut

use warnings;
use strict;
use Getopt::Long qw(:config no_ignore_case);
use Data::Dumper;
use Pod::Usage;
use POSIX qw/strftime/;
use FindBin;
use lib "$FindBin::Bin";
use SortGenomicCoordinates;
use VcfReader;

my $help;
my $manual;
my $vcf;
my $outfile;
my @bedfile;
my @reg;
my $offset;
my $exclude;
my $total_lines = 0;
my $no_header;
my $quiet;
my $silent;

my %opts = (
    regions   => \@reg,
    no_header => \$no_header,
    output    => \$outfile,
    input     => \$vcf,
    bed       => \@bedfile,
    quiet     => \$quiet,
    silent    => \$silent,
    help      => \$help,
    manual    => \$manual
);
GetOptions(
    'exclude'      => \$exclude,
    'regions=s{,}' => \@reg,
    'output=s'     => \$outfile,
    'no_header'    => \$no_header,
    'input=s'      => \$vcf,
    'bed=s{,}'     => \@bedfile,
    'quiet'        => \$quiet,
    'silent'       => \$silent,
    'help'         => \$help,
    'manual'       => \$manual
) or pod2usage( -message => "Syntax error.", -exitval => 2 );
pod2usage( -verbose => 2 ) if ($manual);
pod2usage( -verbose => 1 ) if ($help);
pod2usage( -message => "Syntax error.", -exitval => 2 )
  if ( not $vcf or ( not @bedfile and not @reg ) );
$quiet++ if $silent;
my $time = strftime( "%H:%M:%S", localtime );
print STDERR "Time started: $time\n" unless $quiet;
$offset = 0 if not $offset;
my $total_variants;
print STDERR "[$time] Initializing input VCF... " unless $quiet;
die "Header not ok for input ($vcf) "
  if not VcfReader::checkHeader( vcf => $vcf );
my %sargs = VcfReader::getSearchArguments($vcf);
my %contigs = ();
if (exists $sargs{contig_order}){ 
    foreach my $k (keys %{$sargs{contig_order}}){
        if (ref $sargs{contig_order}->{$k} eq 'HASH' && exists $sargs{contig_order}->{$k}->{order}){
	    $contigs{$k} = $sargs{contig_order}->{$k}->{order};
        }
    }
}else{
    %contigs = VcfReader::getContigOrder($vcf);
}

$time = strftime( "%H:%M:%S", localtime );
print STDERR "\n[$time] Finished initializing input VCF\n" unless $quiet;
$time = strftime( "%H:%M:%S", localtime );
print STDERR "[$time] Preparing regions... " unless $quiet;
my @regions = ();
if (@bedfile) {
    foreach my $bedfile (@bedfile) {
        open( BED, $bedfile ) || die "can't open $bedfile: $!";
        my @temp = ();
        while ( my $temp = <BED> ) {
            chomp $temp;
            next if not $temp;
            $temp =~ s/[:-]/\t/g;
            $temp =~ s/\,//g;
            push( @temp, $temp );
        }

        #    my @temp = (<BED>) =~ s/[:-]/\t/g;
        close BED;

        my $invalid ;
        foreach my $t (@temp){
            if ($t !~ /\S+\t\d+\t\d+\s/ and $t !~ /\S+\t\d+\t\d+$/){
                $invalid++;
            }else{
                push( @regions, $t);
            }
        }
        warn "$invalid invalid region(s) found in $bedfile\n" 
            if $invalid and not $silent;
    }
}
if (@reg) {
    foreach my $reg (@reg) {
        $reg =~ s/\,//g;
        if ($reg =~ /^\S+:\d+-\d+$/){
            (my $r = $reg) =~ s/[\:\-]/\t/g;#keep original intact for option string
            push( @regions, $r);
        }elsif($reg =~ /^\S+:(\d+)$/){
            my $r = "$reg-$1";
            $r =~ s/[\:\-]/\t/g;
            push @regions, $r; 
        }else{
            die "invalid region specied: $reg";
        }
    }
}

my $reg_obj;
if (%contigs){
    my %invalid_contigs = ();
    my @indices_to_remove = ();
    for (my $i = 0; $i <  @regions; $i++){
        my $c = (split "\t", $regions[$i])[0];
        if (not exists $contigs{$c}){
            $invalid_contigs{$c}++;
            unshift @indices_to_remove, $i;
        }
    }
    if (@indices_to_remove){
        my $inv_contigs = keys %invalid_contigs;
        warn "$inv_contigs contigs in regions not found in VCF (". scalar(@indices_to_remove)." regions).\n" unless $silent;
        foreach my $i (@indices_to_remove){
            splice(@regions, $i, 1);
        }
    }
    die "No valid regions to parse.\n" if not @regions;
    $reg_obj = SortGenomicCoordinates->new( 
         array => \@regions, 
         type => 'bed', 
         col => 1, 
         contig_order => \%contigs
    );
}else{
    die "No valid regions to parse.\n" if not @regions;
     $reg_obj = SortGenomicCoordinates->new( 
         array => \@regions, 
         type => 'bed', 
         col => 1, 
    );
}

my $region_ref = $reg_obj->prep();
$time = strftime( "%H:%M:%S", localtime );
print STDERR "\n[$time] Finished preparing regions.\n" unless $quiet;
my $OUT;
if ($outfile) {
    open( $OUT, ">$outfile" ) || die "Can't open $outfile for writing: $!";
}
else {
    $OUT = *STDOUT;
}
my $lines        = 0;
my $next_update  = 0;
my $printed      = 0;
my $meta_head    = VcfReader::getMetaHeader($vcf);
if (not $no_header){
    print $OUT "$meta_head\n";
    print $OUT "##getVariantsByLocation.pl=\"";
    my @opt_string = ();

    foreach my $k ( sort keys %opts ) {
        if ( not ref $opts{$k} ) {
            push @opt_string, "$k=$opts{$k}";
        }
        elsif ( ref $opts{$k} eq 'SCALAR' ) {
            if ( defined ${ $opts{$k} } ) {
                push @opt_string, "$k=${$opts{$k}}";
            }
            else {
                push @opt_string, "$k=undef";
            }
        }
        elsif ( ref $opts{$k} eq 'ARRAY' ) {
            if ( @{ $opts{$k} } ) {
                push @opt_string, "$k=" . join( ",", @{ $opts{$k} } );
            }
            else {
                push @opt_string, "$k=undef";
            }
        }
    }
    my $col_header = VcfReader::getColumnHeader($vcf);
    print $OUT join( " ", @opt_string ) . "\"\n" . $col_header . "\n";
}
my %potential_overlaps = ();
my $prev_chrom = '';
foreach my $r (@$region_ref) {
    foreach my $k (keys %potential_overlaps){
        if ($potential_overlaps{$k} < $r->{start} 
            or $prev_chrom ne $r->{chrom}){
            delete $potential_overlaps{$k};
        }
    }
    my @hits = VcfReader::searchByRegion(
        %sargs,
        chrom => $r->{chrom},
        start => $r->{start},
        end   => $r->{end},
    );

    #one variant may overlap two regions 
    # - need to ensure we only print these variants once...
    #  - store deletion lines in %potential_overlaps 
    #until current region has overtaken the end of that deletion

    foreach my $h (@hits) {
        #print join( "\t", @$h ) . "\n";
        my @split = split("\t", $h);
        my $span = VcfReader::getSpan(\@split);
        if (not exists $potential_overlaps{$h}){
            print $OUT "$h\n";
            $printed++;
            if ($span > $r->{end}){
                $potential_overlaps{$h} = $span;
            }
        }
    }
    $prev_chrom = $r->{chrom};
}

$reg_obj->DESTROY();
$time = strftime( "%H:%M:%S", localtime );
print STDERR "Time finished: $time\n" unless $quiet;
print STDERR "$printed variants retained.\n" unless $quiet;
if (exists $sargs{file_handle}){
    close $sargs{file_handle};
}
close $OUT;
