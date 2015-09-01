=head1 NAME

ClinVarReader.pm - reads and parses the ClinVar file curated by the MacArthur lab (https://github.com/macarthur-lab/clinvar).

=head1 SYNOPSIS

    use ClinVarReader;

    my ($bgzipped, $col_hash) = ClinVarReader::checkClinVarFile('clinvar.tsv');
    #checks input file and creates a bgzip compressed and tabix indexed version if not already compressed and indexed

    bgzip clinvar.tsv && tabix -S 1 -s 1 -b 2 -e 2 clinvar.tsv.gz
    ClinVarReader::checkClinVarFile('clinvar.tsv.gz');
    #(alternatively compress and index manually before checking)

    my %search_args = ClinVarReader::getClinVarSearchArg('clinvar.tsv.gz');
    my @matches = ClinVarReader::searchClinVar
    (
        %search_args,
        chrom => 1, 
        start => 5922867, 
        end   => 6052533,
    );
    #(get all ClinVar variants between specified genomic coordinates)

    my @matches = ClinVarReader::searchForMatchingVariant
    (
        %search_args,   
        chrom => 1, 
        pos   => 5947454, 
        ref   => 'G',
        alt   => 'A', 
    );
    #(search for matching variant in ClinVar)

    my %col = ClinVarReader::getClinVarColumns('clinvar.gz');
    #get column names and numbers in a hash
    
    my @pathogenic = ClinVarReader::getPathogenic(\@matches, \%col);
    #get pathogenic variants from search results

    my @nonpathogenic = ClinVarReader::getNonPathogenic(\@matches, %col);
    #get non-pathogenic variants from search results

    my @conflicted = ClinVarReader::getConflicted(\@pathogenic, %col);
    #get conflicted variants from search results

    my @unconflicted = ClinVarReader::getUnconflicted(\@pathognic, %col);
    #get unconflicted variants from search results

    foreach my $m (@matches){
        push @traits, ClinVarReader::getColumnValues($m, 'all_traits', %col);
    }


=head1 DESCRIPTION

=head2 Overview

This module provides simple methods for retrieving information from the ClinVar file curated by the MacArthur lab (https://github.com/macarthur-lab/clinvar).

=over 8 

=cut

package ClinVarReader;

#parses and searches clinvar file from https://github.com/macarthur-lab/clinvar

use strict;
use warnings;
use Carp;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use Tabix;

=item B<checkClinVarFile>

Checks input file and croaks if it is not a valid input file. If the input is not bgzip compressed and tabix indexed bgzip compression and indexing will be performed (depends on having bgzip and tabix installed). 

Returns the bgzip compressed filename and a reference to a hash of column names to column nhumbers.

    my ($cvar, $col_hash) = ClinVarParser::checkClinVarFile('clinvar.tsv');

=cut 

sub checkClinVarFile{
#returns columns hash of col name to col number
    my $cv = shift;
    if ($cv !~ /\.(b)*gz$/){
        print STDERR "ClinVar file ($cv) is missing a .gz extension - " .
            "attempting to compress with bgzip.\n";
        return _compressClinVarTsv($cv); 
    }elsif( not -e "$cv.tbi"){
        print STDERR "No index found for ClinVar file ($cv). Attempting to index with tabix...\n";
        return indexClinVar($cv);
    }elsif( -M "$cv.tbi" > -M $cv){
        print STDERR "Tabix index is older than ClinVar file ($cv). Attempting to re-index with tabix...\n";
        return indexClinVar($cv);
    }else{
        my %col = getClinVarColumns($cv);
        return $cv, \%col;
    }
}


sub _compressClinVarTsv{
    my $cv = shift;
    if (not `which bgzip`){
        die "Could not find bgzip executable - please install bgzip and ensure ".
            "it is in your PATH or compress and index $cv manually.\n";
    }
    open (my $CVAR, $cv) or die "Can't open $cv for reading header: $!\n";
    my @header = split("\t", <$CVAR>); 
    my %columns = getClinVarColumns($cv);
    print STDERR "Compressing $cv with bgzip...\n";
    system("bgzip -c $cv > $cv.gz"); 
    _checkExit($?);
    $cv = "$cv.gz";
    print STDERR "Attempting to index with tabix...\n";
    indexClinVar($cv, \%columns); 
    return ($cv, \%columns);
}


=item B<indexClinVar>

Takes a bgzip compressed clinvar.tsv.gz file and indexes with tabix automatically determining chromosome and position columns. Requires tabix to be installed and in your PATH. Returns a hash of column names to column numbers.

  my %col = ClinVarReader::indexClinVar('clinvar.gz');

=cut

sub indexClinVar{
    my $cv = shift;
    my %columns = ();
    if (my $col = shift){
        %columns = %$col;
    }else{
        %columns = getClinVarColumns($cv);
    }
    my $seqcol = $columns{chrom} + 1;
    my $poscol = $columns{pos}   + 1;
    system("tabix -S 1 -s $seqcol -b $poscol -e $poscol $cv"); 
    _checkExit($?);
    return %columns;
}
   
=item B<getClinVarColumns>

Reads a clinvar.tsv file and produces a hash of column names to column numbers (0-based).

  my %col = ClinVarReader::getClinVarColumns('clinvar.gz');


=cut

sub getClinVarColumns{
    my $cv = shift;
    my $CVAR; 
    if ($cv =~ /\.(b)*gz$/){
        $CVAR = new IO::Uncompress::Gunzip $cv, MultiStream => 1 
          or die "IO::Uncompress::Gunzip failed while opening $cv for reading: \n$GunzipError";
    }else{
        open ($CVAR, $cv) or die "Can't open $cv for reading header: $!\n";
    }
    chomp (my @header = split("\t", <$CVAR>)); 
    my $n = 0;
    my %columns = map { lc($_) => $n++} @header;
    for my $c ( qw / 
                    chrom
                    pos
                    ref
                    alt
                    mut
                    clinical_significance
                    pathogenic
                    conflicted
                    all_traits
                    all_pmids
                / ) { 
        if (not exists $columns{$c}){
            die "Required column ($c) not found in ClinVar file ($cv) header.\n";
        }
    }
    return %columns;
}

=item B<getClinVarSearchArgs>

Reads a clinvar.tsv file and produces a hash of arguments to be passed to the search methods below. Specifically it gives a hash with the keys 'tabix_iterator' pointing to a tabix iterator and 'col_hash' pointing to an anonymous hash of column names to numbers (0-based).

    my %sargs = getClinVarSearchArgs('clinvar.tsv');

=cut

sub getClinVarSearchArgs{
    my $file = shift;
    #use bgzip compressed version of clinvar file from https://github.com/macarthur-lab/clinvar
    my ($bgz, $clinVarCols) = checkClinVarFile($file);
    my $index = "$bgz.tbi";
    my $iterator = Tabix->new(-data =>  $bgz, -index => $index) ;
    my %sargs = ( tabix_iterator => $iterator, col_hash => $clinVarCols ); 
    return %sargs;
}
=item B<searchClinVar>

Returns all lines representing variants between given genomic coordinates.

Arguments:

=over 16

=item tabix_iterator

A tabix iterator for the clinvar file. Included in the search arguments hash produced by the 'getClinVarSearchArgs' method.

=item chrom

Chromosome to search. Required.

=item start

Starting coordinate to search. Required. 

=item end

End coordinate of region to search. Optional - defaults to start coordinate.

=back


    my %sargs = ClinVarReader::getClinVarSearchArgs('clinvar.tsv');
    my @matches = ClinVarReader::searchClinVar
    (
        %search_args,
        chrom => 1, 
        start => 5922867, 
        end   => 6052533, 
    );

=cut 

################################################################
sub searchClinVar {
    my %args = @_;
    croak "tabix_iterator argument is required for searchClinVar method "
      if not exists $args{tabix_iterator};
    croak "chrom argument is required for searchClinVar method "
      if not exists $args{chrom};
    croak "start argument is required for searchClinVar method "
      if not exists $args{start};
    if ( not exists $args{end} ) {
        $args{end} = $args{start};
    }

    my $iter = $args{tabix_iterator}
      ->query( $args{chrom}, $args{start} - 1, $args{end} );
    return
      if not defined $iter->{_}
      ;    #$iter->{_} will be undef if our chromosome isn't in the vcf file
    my @matches = ();
    while ( my $m = $args{tabix_iterator}->read($iter) ) {
        push @matches, $m;
    }
    return @matches if defined wantarray;
    carp "searchClinVar called in void context ";
}


=item B<searchForMatchingVariant>

Returns all lines representing ClinVar entries matching a given variant.

Arguments:

=over 16

=item tabix_iterator

A tabix iterator for the clinvar file. Included in the search arguments hash produced by the 'getClinVarSearchArgs' method.

=item col_hash

A reference to a hash of column names to column numbers. Included in the search arguments hash produced by the 'getClinVarSearchArgs' method.

=item chrom

Chromosome to search. Required.

=item pos

Chromosome position to search. Required. 

=item ref

Reference allele to match. Required.

=item alt

Alt (variant) allele to match. Required.

=back

    my %sargs = ClinVarReader::getClinVarSearchArgs('clinvar.tsv');
    my @matches = ClinVarReader::searchForMatchingVariant
    (
        %sargs,
        chrom => 1, 
        pos   => 5947454, 
        ref   => 'G',
        alt   => 'A', 
    );

=cut

################################################################
sub searchForMatchingVariant {
    my %args = @_;

    #args are chrom, pos, ref and alt for variant
    croak "chrom argument is required for searchForMatchingVariant method "
      if not exists $args{chrom};
    croak "pos argument is required for searchForMatchingVariant method "
      if not exists $args{pos};
    croak "ref argument is required for searchForMatchingVariant method "
      if not exists $args{ref};
    croak "alt argument is required for searchForMatchingVariant method "
      if not exists $args{alt};
    croak "col_hash argument is required for searchForMatchingVariant method "
      if not exists $args{col_hash};
    croak "tabix_iterator argument is required for searchClinVar method "
      if not exists $args{tabix_iterator};
    my @matches = ();
    my @hits    = searchClinVar(
        chrom => $args{chrom},
        start => $args{pos},
        end   => $args{pos} + length( $args{ref} ) - 1,
        tabix_iterator => $args{tabix_iterator},
    );
    foreach my $h (@hits) {

        if (
             alleleMatchesClinVar(
                %args, cvar => $h,
            )
          )
        {
            push @matches, $h;
        }
    }
    return @matches;
}

##########################################################
sub alleleMatchesClinVar {
    my %args = @_;

    #args are chrom, pos, ref and alt for variant
    # and cvar is a single clinvar entry
    croak "chrom argument is required for searchForMatchingVariant method "
      if not exists $args{chrom};
    croak "pos argument is required for searchForMatchingVariant method "
      if not exists $args{pos};
    croak "ref argument is required for searchForMatchingVariant method "
      if not exists $args{ref};
    croak "alt argument is required for searchForMatchingVariant method "
      if not exists $args{alt};
    croak "col_hash argument is required for searchForMatchingVariant method "
      if not exists $args{col_hash};
    croak "cvar argument is required for searchForMatchingVariant method "
      if not exists $args{cvar};

    #returns 1 if it matches
    my @split = split( "\t", $args{cvar} );
    my $chrom = $split[ $args{col_hash}->{chrom} ];
    next if uc($chrom) ne uc( $args{chrom} );
    my $pos = $split[ $args{col_hash}->{pos} ];

    #alleles should already be minimized
    my $ref = $split[ $args{col_hash}->{ref} ];
    my $alt = $split[ $args{col_hash}->{alt} ];
    return 0 if $pos != $args{pos};
    return 0 if uc($ref) ne uc( $args{ref} );
    return 0 if uc($alt) ne uc( $args{alt} );
    return 1;
}
=item B<getPathogenic> 

Returns an array of lines with a pathogenic or likely pathogenic annotation as indicated by the 'pathogenic' column of the clinvar.tsv file. Takes a reference to an array of lines (as retrieved by searchClinVar or searchForMatchingVariant methods) as the first argument and a reference to a hash of columns names to numbers as retrieved from the 'getClinVarColumns' method as the second argument. 

    my %col = ClinVarReader::getClinVarColumns('clinvar.gz');
    my @pathogenic = ClinVarReader::getPathogenic(\@matches, \%col);
=cut

sub getPathogenic {
    my ( $lines, $col_hash ) = @_;
    my @path = ();
    foreach my $l (@$lines) {
        if ( getColumnValue( $l, 'pathogenic', $col_hash) ) {
            push @path, $l;
        }
    }
    return @path;
}

=item B<getNonPathogenic> 

Returns an array of lines without a pathogenic or likely pathogenic annotation as indicated by the 'pathogenic' column of the clinvar.tsv file. Takes a reference to an array of lines (as retrieved by searchClinVar or searchForMatchingVariant methods) as the first argument and a reference to a hash of columns names to numbers as retrieved from the 'getClinVarColumns' method as the second argument. 

    my %col = ClinVarReader::getClinVarColumns('clinvar.gz');
    my @nonpathogenic = ClinVarReader::getNonPathogenic(\@matches, \%col);

=cut
sub getNonPathogenic {
    my ( $lines, $col_hash ) = @_;
    my @path = ();
    foreach my $l (@$lines) {
        if ( ! getColumnValue( $l, 'pathogenic', $col_hash ) ) {
            push @path, $l;
        }
    }
    return @path;
}

=item B<getConflicted> 

Returns an array of lines indicated as 'conflicted' (that is have been described as both pathogenic and non-pathogenic previously) as indicated by the 'conflicted' column of the clinvar.tsv file. Takes a reference to an array of lines (as retrieved by searchClinVar or searchForMatchingVariant methods) as the first argument and a reference to a hash of columns names to numbers as retrieved from the 'getClinVarColumns' method as the second argument. 

    my %col = ClinVarReader::getClinVarColumns('clinvar.gz');
    my @conflicted = ClinVarReader::getConflicted(\@pathogenic, \%col);

=cut

sub getConflicted {
    my ( $lines, $col_hash ) = @_;
    my @path = ();
    foreach my $l (@$lines) {
        if ( getColumnValue( $l, 'conflicted', $col_hash ) ) {
            push @path, $l;
        }
    }
    return @path;
}

=item B<getUnconflicted> 

Returns an array of lines not indicated as 'conflicted' (that is have not been described as both pathogenic and non-pathogenic previously) as indicated by the 'conflicted' column of the clinvar.tsv file. Takes a reference to an array of lines (as retrieved by searchClinVar or searchForMatchingVariant methods) as the first argument and a reference to a hash of columns names to numbers as retrieved from the 'getClinVarColumns' method as the second argument. 

    my %col = ClinVarReader::getClinVarColumns('clinvar.gz');
    my @unconflicted = ClinVarReader::getUnconflicted(\@pathognic, \%col);

=cut

sub getUnconflicted {
    my ( $lines, $col_hash ) = @_;
    my @path = ();
    foreach my $l (@$lines) {
        if ( ! getColumnValue( $l, 'conflicted', $col_hash ) ) {
            push @path, $l;
        }
    }
    return @path;
}

=item B<getColumnValue>

For a given line returns the value for a specified column. The first argument is a line from the clinvar.tsv file, the second argument is the name of the column to return and the third is a reference to a hash of columns names to numbers as retrieved from the 'getClinVarColumns' method.

    my %col = ClinVarReader::getClinVarColumns('clinvar.gz');
    my $traits = ClinVarReader::getColumnValue($line, 'all_traits', \%col);

=cut

sub getColumnValue {
    my ( $line, $column, $col_hash ) = @_;
    if ( not exists $col_hash->{$column} ) {
        carp "Column $column does not exist in ClinVar file!\n";
        return;
    }
    return ( split "\t", $line )[ $col_hash->{$column} ];
}



sub _checkExit{
    my $exit = shift;
    return if not $exit;
    if ($exit == -1) {
        print "failed to execute: $!\n";
    }elsif ($exit & 127) {
        printf "command died with signal %d, %s coredump\n",
        ($exit & 127),  ($exit & 128) ? 'with' : 'without';
    }elsif ($exit) {
        printf "command exited with value %d\n", $exit >> 8;
    }
    die "Error executing command. Exiting.\n";
}
 

=back


=head1 AUTHOR

David A. Parry

=head1 COPYRIGHT AND LICENSE

Copyright 2015  David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.


=cut 
1;
