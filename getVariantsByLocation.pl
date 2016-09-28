#!/usr/bin/env perl

=head1 NAME

getVariantsByLocation.pl - print variants from a vcf file that lie within a list of genomic regions


=head1 SYNOPSIS
 
        getVariantsByLocation.pl -i vars.vcf -r chr1:100000-200000 [options]
        getVariantsByLocation.pl -h (display help message)
        getVariantsByLocation.pl -m (display manual page)

=cut

=head1 ARGUMENTS
   

=over 8

=item B<-i    --input>

VCF format file containing list of variants to filter (required). 

=item B<-b    --bed>

Bed file containing regions to filter on (required unless using --regions argument). Can also be a file containing a list of intervals, one interval per line.

=item B<-r    --regions>

Specify one or more regions to filter on seperated with spaces in format "X:1-2000" or (for single coordinates) "X:2000". Can be used in conjunction with or instead of --bed argument.

=item B<-g    --gene_ids>

One or more gene symbols or IDs to use to identify regions. This method uses Ensembl's REST server to identify coordinates of genes, by default using the GRCh37 reference.

=item B<-l    --gene_list>

A list of gene symbols or IDs to use to identify regions. 

=item B<-s    --species>

If you want to use a species other than human to identify matching genes, specify the species name here.

=item B<-d    --defaultRestServer>

Use this option if you want to use the default Ensembl REST server for gene queries (if you want to use GRCh38 rather than GRCh37, for example).

=item B<-v    --vcf_filter> 

A VCF file to use to find variants overlapping variants from your input. By default, variants from your input file will be printed if their coordinates overlap with any variants in these VCF files and variants are not checked to see whether they share the same alleles. To only print variants with matching alleles use the --matching argument. 

Note, that the VCFs used in this argument do not have to be coordinate sorted, but that it is possible that you will output duplicated lines from your input if they are not. Your output will be in the same order as this VCF.

This option can not be used in conjunction with --regions, --bed, --gene_list or --gene_ids  arguments.

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

 getVariantsByLocation.pl -i [vars.vcf] -g SMG1 -o [filtered.vcf]
 #get variants from vars.vcf that lie within the SMG1 gene (coordinates from GRCh37) and output to file filtered.vcf
 
 getVariantsByLocation.pl -i [vars.vcf] -g SMG1 -d -o [filtered.vcf]
 #get variants from vars.vcf that lie within the SMG1 gene using the current Ensembl genome build and output to file filtered.vcf
 
 getVariantsByLocation.pl -i [vars.vcf] -v [other.vcf]
 #get variants from vars.vcf that overlap any variant in other.vcf 
   
 getVariantsByLocation.pl -i [vars.vcf] -v [other.vcf] --matching
 #get variants from vars.vcf that match a variant in other.vcf
 
=cut

=head1 DESCRIPTION

Outputs variants from a VCF that overlap given regions. Regions can be specified on the commandline or in a BED file. Alternatively, another VCF file can be supplied in order to output overlapping or matching variants or gene IDs can be supplied to output variants that lie withing those genes.
 
The VCF input must be sorted in coordinate order but can be either uncompressed (a .vridx index file will be created if it does not exist) or bgzip compressed (a tabix index will be created if it does not exist). Use with bgzip compressed VCFs requires the Bio::DB::HTS::Tabix module to be installed and tabix index creation requires the tabix executable to be installed and in your PATH.


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
use FindBin qw($RealBin);
use lib "$RealBin/lib";
use SortGenomicCoordinates;
use VcfReader;
use IdParser;
use EnsemblRestQuery;

my @bedfile;
my @reg;
my @gene_ids;
my $restQuery ;
my $id_parser = new IdParser();

my %opts = 
(
    r => \@reg,
    b => \@bedfile,
    g => \@gene_ids,
);

GetOptions
(
    \%opts,
   # 'exclude', not implmented
    'r|regions=s{,}',
    'o|output=s',
    'n|no_header',
    'i|input=s',
    'b|bed=s{,}',
    'v|vcf_filter=s',
    'l|gene_list=s',
    'g|gene_ids=s{,}',
    'e|matching',
    'q|quiet',
    's|silent',
    'species=s',
    'd|defaultRestServer',
    'h|?|help',
    'm|manual',
) or pod2usage( -message => "Syntax error.", -exitval => 2 );

pod2usage( -verbose => 2 ) if ($opts{m});
pod2usage( -verbose => 1 ) if ($opts{h});

pod2usage
(
     -message => "ERROR: -i/--input argument is required.\n", 
     -exitval => 2 
) if ( not $opts{i}); 

pod2usage
(
     -message => "ERROR: A regions argument must be supplied. ".
                 "See -r/--regions, -b/--bed, -g/--gene_ids, ".
                 "-l/--gene_list and -v/--vcf_filter arguments in ". 
                 "the  help or manual documentation of this program.\n", 
     -exitval => 2 
) if ( not @bedfile and not @reg and not @gene_ids and not $opts{v}  and not $opts{l} ) ;

pod2usage( -message => "ERROR: --vcf_filter option can not be used in conjunction with --bed, --region, --gene_ids or --gene_list arguments.\n", -exitval => 2 )
  if ( (@bedfile or @reg or @gene_ids) and $opts{v} );

if ($opts{e} and not $opts{v}){
    print STDERR "WARNING: redundant use of --matching argument without --vcf_filter argument.\n" unless $opts{s};
}

$opts{q}++ if $opts{s};
$opts{species} ||= 'human'; 

if (lc($opts{species}) !~  /^human|homo sapiens|h_sapiens|homo|enshs|hsap|9606|homsap|hsapiens$/){
    $opts{d}++;#always use default server if not querying human
}

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

#########################################################
sub process_regions{
    $time = strftime( "%H:%M:%S", localtime );
    my @regions = ();
    if ($opts{l} or @gene_ids){
    #make EnsemblRestQuery modules non-essential if not using genes
        eval "use EnsemblRestQuery; 1 " or die "Missing required module " . 
          "for --gene_list based queries. Please install the module indicated below " .
          "and try again:\n\n$@\n"; 
        $restQuery = new EnsemblRestQuery();
        $restQuery->useGRCh37Server() unless $opts{d};
    }
    if ($opts{l}){
        open (my $GENES, $opts{l}) or die "Could not open --gene_list '$opts{l}' for reading: $!\n";
        while (my $line = <$GENES>){
            my @s = split(/\s+/, $line); 
            push @gene_ids, $s[0];
        }
    }
    my %seen = ();
    @gene_ids = grep { ! $seen{$_}++} @gene_ids;
    foreach my $g (@gene_ids){
        my $region = get_region_from_gene($g);
        push @regions, $region if defined $region;
    }
    
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
        my %converted_contigs = ();
        my @indices_to_remove = ();
        for (my $i = 0; $i <  @regions; $i++){
            my $c = (split /[:-]/, $regions[$i])[0];
            if (not exists $contigs{$c}){
                if (exists $contigs{"chr$c"}){
                    $converted_contigs{$c} = "chr$c";
                    $regions[$i] =~ s/^$c/chr$c/;
                }else{
                    (my $cc = $c) =~ s/^chr//;
                    if (exists $contigs{$cc}){
                        $converted_contigs{$c} = $cc;
                        $regions[$i] =~ s/^$c/$cc/;
                    }else{
                        $invalid_contigs{$c}++;
                        unshift @indices_to_remove, $i;
                    }
                }
            }
        }
        if (@indices_to_remove){
            my $inv_contigs = keys %invalid_contigs;
            warn "$inv_contigs contigs in regions not found in VCF ("
                . scalar(@indices_to_remove)." regions).\n" unless $opts{s};
            foreach my $i (@indices_to_remove){
                splice(@regions, $i, 1);
            }
        }
        if (my $conv = keys %converted_contigs){
            warn "$conv contigs in regions were converted assuming only ". 
                "difference between VCF contigs is prepended 'chr'.\n".
                "The Following contigs were altered:\n" 
                .join
                (  "\n", 
                    map { 
                        sprintf
                        (
                                "%38s => %s", 
                                $_ , 
                                $converted_contigs{$_}
                        ) 
                    } sort keys %converted_contigs
                ) . "\n";
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
}

#########################################################
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

#########################################################
sub geneFromEnst{
    my $id = shift;
    if (not $opts{q}){
        print STDERR "Identifying parent gene from Ensembl transcript $id...\n";
    }
    return $restQuery->getParent($id, 1);
}

#########################################################
sub geneFromEnsp{
    my $id = shift;
    if (not $opts{q}){
        print STDERR "Identifying parent gene from Ensembl protein $id...\n";
    }
    my $par = $restQuery->getParent($id);
    if ($par){
        if (exists $par->{id}){
            return geneFromEnst($par->{id});
        }
    }
}

#########################################################
sub get_region_from_gene{
    my $id = shift;
    $id_parser->parseId($id);
    my $gene_hash; 
    my @lookups = ();
    if (not $opts{q}){
        print STDERR "Interpretting ID \"$id\" as of type \"" . 
          $id_parser->get_identifierType() . "\"...\n";
    }
    if ($id_parser->get_isEnsemblId()){
        if ( $id_parser->get_isTranscript() ){
            $gene_hash = geneFromEnst($id);
        }elsif( $id_parser->get_isProtein() ) {
            $gene_hash = geneFromEnsp($id);
        }else{
            $gene_hash = $restQuery->lookUpEnsId($id, 1);
        }
    }elsif($id_parser->get_isTranscript()  or $id_parser->get_isProtein() ) {
        if (not $opts{q}){
            print STDERR "Identifying Ensembl gene via transcript cross-reference...\n";
        }
        my $transcript = $restQuery->getTranscriptViaXreg($id, $opts{species});
        if ($transcript and ref $transcript eq 'ARRAY'){
            if (@$transcript > 1){
                print STDERR "WARNING: Multiple transcripts identified by ".
                  "cross-reference search for $id - picking the first.\n";
            }
            my $tr = $transcript->[0];
            if (exists $tr->{id}){
                $gene_hash = geneFromEnst($tr->{id});
            }
        }else{
            if (not $opts{s}){
                print STDERR "WARNING: No transcript identified for ID \"$id\"\n";
            }
        }
    }else{
        if (not $opts{q}){
            print STDERR "Identifying Ensembl gene via gene cross-reference...\n";
        }
        my $gene = $restQuery->getGeneViaXreg($id, $opts{species});
        if (ref $gene eq 'ARRAY'){
            foreach my $ge (@$gene){
                if ($ge->{id}){
                    my $ge_hash = $restQuery->lookUpEnsId($ge->{id}, 1);
                    if (uc($ge_hash->{display_name}) eq uc($id)){
                    #if gene symbol matches then we use this entry
                        $gene_hash = $ge_hash;
                        last;
                    }else{
                        push @lookups, $ge_hash;
                    }
                }
            }
            if (not $gene_hash){
                if (@lookups == 1){
                    $gene_hash = $lookups[0];
                }
            }
        }
    }
    if (not $gene_hash){
        if (not $opts{s}){
            print STDERR "WARNING: Could not identify gene for ID \"$id\"\n";
            if (@lookups){
                my $idstring = join("\n", map { $_->{display_name} } @lookups );
                print STDERR "Identified the following non-matching display names:\n".
                             "$idstring\n";
            }
            return;
        }
    }
    my $chrom = $gene_hash->{seq_region_name};
    my $start = $gene_hash->{start};
    my $end = $gene_hash->{end};
    my $r = "$chrom:$start-$end";
    if (not $opts{q}){
        print STDERR "For gene ID '$id', found gene " . 
          $gene_hash->{display_name} . "/" . $gene_hash->{id} .
          " with coordintes $r ($gene_hash->{assembly_name})\n";
    }
    return $r;
}

#########################################################
sub open_output{
    if ($opts{o}) {
        open( $OUT, ">$opts{o}" ) || die "Can't open $opts{o} for writing: $!";
    }
    else {
        $OUT = *STDOUT;
    }
}

#########################################################
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
