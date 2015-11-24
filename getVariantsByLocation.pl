#!/usr/bin/env perl

=head1 NAME

getVariantsByLocation.pl - print variants from a vcf file that lie within a list of genomic regions


=head1 SYNOPSIS
 
getVariantsByLocation.pl -i [vars.vcf] -b [regions.bed] -o [output.vcf]
getVariantsByLocation.pl -i [vars.vcf] -r [regions] -o [output.vcf]
getVariantsByLocation.pl -i [vars.vcf] -v [other.vcf] -o [output.vcf]

=cut

=head1 ARGUMENTS
   

=over 8

=item B<-i    --input>

VCF format file containing list of variants to filter (required). 

=item B<-b    --bed>

Bed file containing regions to filter on (required unless using --regions argument). Can also be a file containing a list of intervals, one interval per line.

=item B<-r    --regions>

Specify one or more regions to filter on seperated with spaces in format "X:1-2000" or (for single coordinates) "X:2000". Can be used in conjunction with or instead of --bed argument.

=item B<-v    --vcf_filter> 

A VCF file to use to find variants overlapping variants from your input. By default, variants from your input file will be printed if their coordinates overlap with any variants in these VCF files and variants are not checked to see whether they share the same alleles. To only print variants with matching alleles use the --matching argument. 

Note, that the VCFs used in this argument do not have to be coordinate sorted, but that it is possible that you will output duplicated lines from your input if they are not. Your output will be in the same order as this VCF.

This option can not be used in conjunction with --regions or --bed arguments.

=item B<-e    --matching> 

When using the --vcf_filter argument, use this flag to only output variants with matching alleles. If variants are multiallelic they will be printed if any of the alleles match.

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
 #get variants from vars.vcf that lie within regions in bedfile regions.bed
 
 getVariantsByLocation.pl -i [vars.vcf] -r 1:2000000-50000000 -o [filtered.vcf]
 #get variants from vars.vcf that lie within the region 1:2000000-50000000 and output to file filtered.vcf

 getVariantsByLocation.pl -i [vars.vcf] -b [regions.bed] -r 1:2000000-50000000 -o [filtered.vcf]
 #get variants from vars.vcf that lie within regions in bedfile regions.bed or the region 1:2000000-50000000 and output to file filtered.vcf

 getVariantsByLocation.pl -i [vars.vcf] -v [other.vcf]
 #get variants from vars.vcf that overlap any variant in other.vcf 
   
 getVariantsByLocation.pl -i [vars.vcf] -v [other.vcf] --matching
 #get variants from vars.vcf that match a variant in other.vcf
 
=cut

=head1 DESCRIPTION

Takes a bed file and/or a list of regions/coordinates and outputs variants from a VCF file that lie within/overlap these regions. The VCF input must be sorted in coordinate order but can be either uncompressed (a .vridx index file will be created if it does not exist) or bgzip compressed (a tabix index will be created if it does not exist). Use with bgzip compressed VCFs requires Tabix.pm to be installed and tabix index creation requires the tabix executable to be installed and in your PATH. 


=cut

=head1 AUTHOR

David A. Parry


=head1 COPYRIGHT AND LICENSE

Copyright 2014,2015 David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut

use warnings;
use strict;
use Getopt::Long;
use Data::Dumper;
use Pod::Usage;
use POSIX qw/strftime/;
use FindBin;
use lib "$FindBin::Bin/lib";
use SortGenomicCoordinates;
use VcfReader;
use IdParser;

my @bedfile;
my @reg;

my %opts = 
(
    r => \@reg,
    b => \@bedfile,
);

GetOptions
(
    'e|exclude',
    'r|regions=s{,}',
    'o|output=s',
    'n|no_header',
    'i|input=s',
    'b|bed=s{,}',
    'v|vcf_filter=s',
    'e|matching',
    'q|quiet',
    's|silent',
    'h|?|help',
    'm|manual',
) or pod2usage( -message => "Syntax error.", -exitval => 2 );

pod2usage( -verbose => 2 ) if ($opts{m});
pod2usage( -verbose => 1 ) if ($opts{h});

pod2usage( -message => "Syntax error.", -exitval => 2 )
  if ( not $opts{i} or ( not @bedfile and not $opts{l} and not @reg and not $opts{v}) );

pod2usage( -message => "ERROR: --vcf_filter option can not be used in conjunction with --bed or --region arguments.", -exitval => 2 )
  if ( (@bedfile or @reg) and $opts{v} );

if ($opts{e} and not $opts{v}){
    print STDERR "WARNING: redundant use of --matching argument without --vcf_filter argument.\n" unless $opts{s};
}
$opts{q}++ if $opts{s};
my $time = strftime( "%H:%M:%S", localtime );
print STDERR "Time started: $time\n" unless $opts{q};
my $total_variants;
print STDERR "[$time] Initializing input VCF... " unless $opts{q};
die "Header not ok for input ($opts{i}) "
  if not VcfReader::checkHeader( vcf => $opts{i} );
my %sargs = VcfReader::getSearchArguments($opts{i});
my %contigs = ();
if (exists $sargs{contig_order}){ 
    foreach my $k (keys %{$sargs{contig_order}}){
        if (ref $sargs{contig_order}->{$k} eq 'HASH' && exists $sargs{contig_order}->{$k}->{order}){
	    $contigs{$k} = $sargs{contig_order}->{$k}->{order};
        }
    }
}else{
    %contigs = VcfReader::getContigOrder($opts{i});
}

$time = strftime( "%H:%M:%S", localtime );
print STDERR "\n[$time] Finished initializing input VCF\n" unless $opts{q};
my $printed      = 0;
my $OUT;
my %potential_overlaps = ();
my $prev_chrom = '';

if ($opts{v}){
    process_vcf_filter();
}else{
    process_regions();
}

$time = strftime( "%H:%M:%S", localtime );
print STDERR "Time finished: $time\n" unless $opts{q};
print STDERR "$printed variants retained.\n" unless $opts{q};
if (exists $sargs{file_handle}){
    close $sargs{file_handle};
}
close $OUT;

sub process_regions{
    $time = strftime( "%H:%M:%S", localtime );
    my @regions = ();
    print STDERR "[$time] Preparing regions... " unless $opts{q};
    if (@bedfile) {
        foreach my $bedfile (@bedfile) {
            open( BED, $bedfile ) || die "can't open $bedfile: $!";
            my @temp = ();
            my $invalid ;
            while ( my $line = <BED> ) {
                chomp $line;
                next if not $line;
                (my $temp = $line) =~ s/\,//g;
                my ($chr, $start, $stop); 
                if ( $temp =~ /\S+:(\d+)(-\d+)*/){
                    ($chr, $start, $stop) = split(/[\:\-]/, $temp); 
                    $stop ||= $start ;
                }elsif ( $temp =~ /\S+\t\d+\t\d+\s/ or $temp =~ /\S+\t\d+\t\d+$/){
                    ($chr, $start, $stop) = split(/\t/, $temp); 
                    $start++;#BED should be 0-based 
                    $stop--;
                }else{
                    $invalid++;
                    next;
                }
                if ($start > $stop){
                    warn "Invalid region found \"$line\" in bedfile - start is greater than end.\n";
                    $invalid++;
                    next;
                }
                push @regions, "$chr:$start-$stop";
            }
            close BED;

            warn "$invalid invalid region(s) found in $bedfile\n" 
                if $invalid and not $opts{s};
        }
    }
    if (@reg) {
        foreach my $reg (@reg) {
            $reg =~ s/\,//g;
            if ($reg =~ /^\S+:(\d+)-(\d+)$/){
                if ($1 > $2){
                    die "Invalid region \"$reg\" - start is greater than end.\n";
                }
                push( @regions, $reg);
            }elsif($reg =~ /^\S+:(\d+)$/){
                my $r = "$reg-$1";
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
            my $c = (split /[:-]/, $regions[$i])[0];
            if (not exists $contigs{$c}){
                $invalid_contigs{$c}++;
                unshift @indices_to_remove, $i;
            }
        }
        if (@indices_to_remove){
            my $inv_contigs = keys %invalid_contigs;
            warn "$inv_contigs contigs in regions not found in VCF (". scalar(@indices_to_remove)." regions).\n" unless $opts{s};
            foreach my $i (@indices_to_remove){
                splice(@regions, $i, 1);
            }
        }
        die "No valid regions to parse.\n" if not @regions;
        $reg_obj = SortGenomicCoordinates->new( 
            array => \@regions, 
            type => 'regions', 
            col => 1, 
            contig_order => \%contigs
        );
    }else{
        die "No valid regions to parse.\n" if not @regions;
        $reg_obj = SortGenomicCoordinates->new( 
            array => \@regions, 
            type => 'regions', 
            col => 1, 
        );
    }

    my $region_ref = $reg_obj->prep();
    $time = strftime( "%H:%M:%S", localtime );
    print STDERR "\n[$time] Finished preparing regions.\n" unless $opts{q};
    open_output();
    write_header();
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
            start => $r->{start}+1,
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
}

sub process_vcf_filter{
    my $VCF = VcfReader::_openFileHandle( $opts{v} );
    open_output();
    write_header();
    my $prev_start = 0;
    while ( my $line = <$VCF> ) {
        next if $line =~ /^#/;
        my @filter_split = split("\t", $line);
        my $chrom = VcfReader::getVariantField(\@filter_split, 'CHROM');
        my $start = VcfReader::getVariantField(\@filter_split, 'POS');
        my $filt_span = VcfReader::getSpan(\@filter_split);
        foreach my $k (keys %potential_overlaps){
            if ($potential_overlaps{$k} < $start 
                or $prev_chrom ne $chrom){
                delete $potential_overlaps{$k};
            }
        }
        my @hits = VcfReader::searchByRegion(
            %sargs,
            chrom => $chrom,
            start => $start,
            end   => $filt_span,
        );
        foreach my $h (@hits) {
            my @split = split("\t", $h);
            my $span = VcfReader::getSpan(\@split);
            if (not exists $potential_overlaps{$h}){
                if ($opts{e}){
                    if (VcfReader::variantsHaveMatchingAlleles(\@split, \@filter_split)){
                        print $OUT "$h\n";
                        $printed++;
                        $potential_overlaps{$h} = $span;
                    }
                }else{
                    print $OUT "$h\n";
                    $printed++;
                    $potential_overlaps{$h} = $span;
                }
            }
        }
        $prev_chrom = $chrom;
    }
}
sub getRegionFromGene{
    my $l = shift;
    my @s = split(/\s+/, $l); 
    
    my $parser = new IdParser();
    
    
}
sub open_output{
    if ($opts{o}) {
        open( $OUT, ">$opts{o}" ) || die "Can't open $opts{o} for writing: $!";
    }
    else {
        $OUT = *STDOUT;
    }
}

sub write_header{
    if ($opts{n}){
        return ;
    }
    my $meta_head    = VcfReader::getMetaHeader($opts{i});
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
    my $col_header = VcfReader::getColumnHeader($opts{i});
    print $OUT join( " ", @opt_string ) . "\"\n" . $col_header . "\n";
}
