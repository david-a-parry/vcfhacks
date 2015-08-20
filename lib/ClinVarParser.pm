package ClinVarParser;

#parses and searches clinvar file from https://github.com/macarthur-lab/clinvar

use strict;
use warnings;
use Carp;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use Tabix;
our $AUTOLOAD;
{
    my $_count = 0;
    my %_attrs = (
        _file    => [ "", "read/required" ],
        _columns => [ "", "read" ],
    );

    sub _all_attrs {
        keys %_attrs;
    }

    sub _accessible {
        my ( $self, $attr, $mode ) = @_;
        $_attrs{$attr}[1] =~ /$mode/;
    }

    sub _attr_default {
        my ( $self, $attr ) = @_;
        $_attrs{$attr}[0];
    }

    sub get_count {
        $_count;
    }

    sub _incr_count {
        $_count++;
    }

    sub _decr_count {
        $_count--;
    }

}

sub DESTROY {
    my ($self) = @_;
    $self->_decr_count();
}

sub new {
    my ( $class, %args ) = @_;
    my $self = bless {}, $class;
    foreach my $attr ( $self->_all_attrs() ) {
        my ($arg) = ( $attr =~ /^_(.*)/ );
        if ( exists $args{$arg} ) {
            $self->{$attr} = $args{$arg};
        }
        elsif ( $self->_accessible( $attr, "required" ) ) {
            croak "$attr argument required";
        }
        else {
            $self->{$attr} = $self->_attr_default($attr);
        }
    }
    $self->_parseClinVarFile();
    $class->_incr_count();
    return $self;
}

###########################################################
sub getTabixIterator {
#use bgzip compressed version of clinvar file from https://github.com/macarthur-lab/clinvar
    my $self = shift;
    my $index = $self->{_file} .".tbi";
    return Tabix->new( -data => $self->{_file}, -index => $index );
}

################################################################
sub searchClinVar {
    my ( $self, %args ) = @_;
    croak "chrom argument is required for searchClinVar method "
      if not exists $args{chrom};
    croak "start argument is required for searchClinVar method "
      if not exists $args{start};
    if ( not exists $args{end} ) {
        $args{end} = $args{start};
    }

    my $iter = $self->{tabixIterator}
      ->query( $args{chrom}, $args{start} - 1, $args{end} );
    return
      if not defined $iter->{_}
      ;    #$iter->{_} will be undef if our chromosome isn't in the vcf file
    my @matches = ();
    while ( my $m = $self->{tabixIterator}->read($iter) ) {
        push @matches, $m;
    }
    return @matches if defined wantarray;
    carp "searchClinVar called in void context ";
}

################################################################
sub searchForMatchingVariant {
    my ( $self, %args ) = @_;

    #args are chrom, pos, ref and alt for variant
    croak "chrom argument is required for searchForMatchingVariant method "
      if not exists $args{chrom};
    croak "pos argument is required for searchForMatchingVariant method "
      if not exists $args{pos};
    croak "ref argument is required for searchForMatchingVariant method "
      if not exists $args{ref};
    croak "alt argument is required for searchForMatchingVariant method "
      if not exists $args{alt};
    my @matches = ();
    my @hits    = $self->searchClinVar(
        chrom => $args{chrom},
        start => $args{pos},
        end   => $args{pos} + length( $args{ref} ) - 1,
    );
    foreach my $h (@hits) {

        if (
            $self->alleleMatchesClinVar(
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
    my ( $self, %args ) = @_;

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
    croak "cvar argument is required for searchForMatchingVariant method "
      if not exists $args{cvar};

    #returns 1 if it matches
    my @split = split( "\t", $args{cvar} );
    my $chrom = $split[ $self->{_columns}->{chrom} ];
    next if uc($chrom) ne uc( $args{chrom} );
    my $pos = $split[ $self->{_columns}->{pos} ];

    #alleles should already be minimized
    my $ref = $split[ $self->{_columns}->{ref} ];
    my $alt = $split[ $self->{_columns}->{alt} ];
    return 0 if $pos != $args{pos};
    return 0 if uc($ref) ne uc( $args{ref} );
    return 0 if uc($alt) ne uc( $args{alt} );
    return 1;
}
######################################################
sub grepResults {
    my ( $self, %args ) = @_;
    croak "lines argument is required for grepResults method "
      if not exists $args{lines};

#croak "lines argument in grepResults method must be an array reference " if ref $args{lines} ne 'ARRAY';
    croak "column argument is required for grepResults method "
      if not exists $args{column};
    croak "term argument is required for grepResults method "
      if not exists $args{term};

    if ( not exists $self->{_columns}->{ $args{column} } ) {
        carp "Column $args{column} does not exist in ClinVar file!\n";
        return;
    }
    my @matches = ();
    foreach my $l ( @{ $args{lines} } ) {
        my $field = ( split "\t", $l )[ $self->{_columns}->{ $args{column} } ];
        if ( $args{exact} ) {
            if ( $field eq $args{term} ) {
                if ( $args{field_only} ) {
                    push @matches, $field;
                }
                else {
                    push @matches, $l;
                }
            }
        }
        elsif ( $field =~ /$args{term}/ ) {
            if ( $args{field_only} ) {
                push @matches, $field;
            }
            else {
                push @matches, $l;
            }
        }
    }
    return @matches;
}

######################################################
sub getPathogenic {
    my ( $self, $lines ) = @_;
    my @path = ();
    foreach my $l (@$lines) {
        if ( $self->getColumnValue( $l, 'pathogenic' ) ) {
            push @path, $l;
        }
    }
    return @path;
}

######################################################
sub getNonPathogenic {
    my ( $self, $lines ) = @_;
    my @path = ();
    foreach my $l (@$lines) {
        if ( !$self->getColumnValue( $l, 'pathogenic' ) ) {
            push @path, $l;
        }
    }
    return @path;
}

######################################################
sub getConflicted {
    my ( $self, $lines ) = @_;
    my @path = ();
    foreach my $l (@$lines) {
        if ( $self->getColumnValue( $l, 'conflicted' ) ) {
            push @path, $l;
        }
    }
    return @path;
}

######################################################
sub getUnconflicted {
    my ( $self, $lines ) = @_;
    my @path = ();
    foreach my $l (@$lines) {
        if ( !$self->getColumnValue( $l, 'conflicted' ) ) {
            push @path, $l;
        }
    }
    return @path;
}

######################################################
sub getColumnValue {
    my ( $self, $line, $column ) = @_;
    if ( not exists $self->{_columns}->{$column} ) {
        carp "Column $column does not exist in ClinVar file!\n";
        return;
    }
    return ( split "\t", $line )[ $self->{_columns}->{$column} ];
}

######################################################
sub _parseClinVarFile {
    my $self = shift;
    if ( $self->{_file} !~ /\.(b)*gz$/ ) {
        print STDERR
          "ClinVar file ($self->{_file}) is missing a .gz extension - "
          . "attempting to compress with bgzip.\n";
        $self->_compressClinVarTsv();
    }
    elsif ( not -e "$self->{_file}.tbi" ) {
        print STDERR
"No index found for ClinVar file ($self->{_file}). Attempting to index with tabix...\n";
        $self->_indexClinVar();
    }
    elsif ( -M "$self->{_file}.tbi" > -M $self->{_file} ) {
        print STDERR
"Tabix index is older than ClinVar file ($self->{_file}). Attempting to re-index with tabix...\n";
        $self->_indexClinVar();
    }
    else {
        $self->{_columns} = $self->_getClinVarColumns();
    }
    $self->{tabixIterator} = $self->getTabixIterator();
}
###########################################################
sub _compressClinVarTsv {
    my $self = shift;
    if ( not `which bgzip` ) {
        croak
          "Could not find bgzip executable - please install bgzip and ensure "
          . "it is in your PATH or compress and index $self->{_file} manually.\n";
    }
    open( my $CVAR, $self->{_file} )
      or croak "Can't open $self->{_file} for reading header: $!\n";
    my @header = split( "\t", <$CVAR> );
    $self->{_columns} = $self->_getClinVarColumns();
    print STDERR "Compressing $self->{_file} with bgzip...\n";
    system("bgzip -c $self->{_file} > $self->{_file}.gz");
    $self->_checkExit($?);
    $self->{_file} = "$self->{_file}.gz";
    $self->_indexClinVar();
}
###########################################################
sub _indexClinVar {
    my $self = shift;
    if ( not $self->{_columns} ) {
        $self->{_columns} = $self->_getClinVarColumns( $self->{_file} );
    }
    my $seqcol = $self->{_columns}->{chrom} + 1;
    my $poscol = $self->{_columns}->{pos} + 1;
    system("tabix -S 1 -s $seqcol -b $poscol -e $poscol $self->{_file}");
    $self->_checkExit($?);
}
###########################################################
sub _getClinVarColumns {
    my $self = shift;
    my $CVAR;
    if ( $self->{_file} =~ /\.(b)*gz$/ ) {
        $CVAR = new IO::Uncompress::Gunzip $self->{_file}, MultiStream => 1
          or croak
"IO::Uncompress::Gunzip failed while opening $self->{_file} for reading: \n$GunzipError";
    }
    else {
        open( $CVAR, $self->{_file} )
          or croak "Can't open $self->{_file} for reading header: $!\n";
    }
    chomp( my @header = split( "\t", <$CVAR> ) );
    my $n = 0;
    my %columns = map { lc($_) => $n++ } @header;
    for my $c (
        qw /
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
        /
      )
    {

        if ( not exists $columns{$c} ) {
            croak
"Required column ($c) not found in ClinVar file ($self->{_file}) header.\n";
        }
    }
    return \%columns;
}
###########################################################
sub _checkExit {
    my $self = shift;
    my $exit = shift;
    return if not $exit;
    if ( $exit == -1 ) {
        print "failed to execute: $!\n";
    }
    elsif ( $exit & 127 ) {
        printf "command died with signal %d, %s coredump\n",
          ( $exit & 127 ), ( $exit & 128 ) ? 'with' : 'without';
    }
    elsif ($exit) {
        printf "command exited with value %d\n", $exit >> 8;
    }
    croak "Error executing command. Exiting.\n";
}

###########################################################

1;

=head1 NAME

ClinVarParser.pm - reads and parses the ClinVar file curated by the MacArthur lab (https://github.com/macarthur-lab/clinvar).

=head1 SYNOPSIS

    use ClinVarParser;

    $cvar = ClinVarParser -> new(file => 'clinvar.tsv');
    #(initialise object with a clinvar file. ClinVarParser.pm will try to bgzip compress and tabix index this file)

    bgzip clinvar.tsv && tabix -S 1 -s 1 -b 2 -e 2 clinvar.tsv.gz
    $cvar = ClinVarParser -> new(file => 'clinvar.tsv.gz');
    #(initialise object with a bgzip compressed clinvar file)

    my @matches = $cvar->searchClinVar
    (
        chrom => 1, 
        start => 5922867, 
        end   => 6052533, 
    );
    #(get all ClinVar variants between specified genomic coordinates)

    my @matches = $cvar->searchForMatchingVariant
    (
        chrom => 1, 
        pos   => 5947454, 
        ref   => 'G',
        alt   => 'A', 
    );
    #(search for matching variant in ClinVar)

    my @pathogenic = $cvar->getPathogenic(\@matches);
    #get pathogenic variants from search results

    my @nonpathogenic = $cvar->getNonPathogenic(\@matches);
    #get non-pathogenic variants from search results

    my @conflicted = $cvar->getConflicted(\@pathogenic);
    #get conflicted variants from search results

    my @unconflicted = $cvar->getUnconflicted(\@pathognic);
    #get unconflicted variants from search results

    my @nphp_var = $cvar->grepResults
    (
        lines  => \@matches, 
        column => "all_traits",
        term   => "Nephronophthisis",
    );

    foreach my $m (@matches){
        push @traits, $cvar->getColumnValues($m, 'all_traits');
    }



=head1 DESCRIPTION

=head2 Overview

This module provides a simple object-oriented interface for retrieving information from the ClinVar file curated by the MacArthur lab (https://github.com/macarthur-lab/clinvar).


=head2 Constructor and initialization

To initialise call the 'new' method specifying a copy of the clinvar.tsv file downloaded from https://github.com/macarthur-lab/clinvar.

    $cvar = ClinVarParser -> new(file => 'clinvar.tsv');

A bgzip compressed and tabix indexed copy will be created if it does not already exist, but you must have both bgzip and tabix installed and in your PATH. Alternatively, you may use a precompressed and indexed version. To compress and index the file you can use the following commands:

    bgzip clinvar.tsv && tabix -S 1 -s 1 -b 2 -e 2 clinvar.tsv.gz

After initializing your ClinVarParser object you can search for matching variants or by region. However, you must have the Tabix.pm perl module installed.


=head2 Class and object methods


=head3 Accessors and mutators

The following features can be accessed using get_[feature] substituting [feature] for the feature of interest. 

=over 8

=item B<file>

filename/path of clinvar.tsv file used.

=item B<columns>

a hash of column names to the 0-based index of the column

=back

=head3 Methods

=over 8

=item B<searchClinVar>

Returns all lines representing variants between given genomic coordinates.


Arguments:

=over 16

=item chrom

Chromosome to search. Required.

=item start

Starting coordinate to search. Required. 

=item end

End coordinate of region to search. Optional - defaults to start coordinate.

=back


    my @matches = $cvar->searchClinVar
    (
        chrom => 1, 
        start => 5922867, 
        end   => 6052533, 
    );


=item B<searchForMatchingVariant>

Returns all lines representing ClinVar entries matching a given variant.

Arguments:

=over 16

=item chrom

Chromosome to search. Required.

=item pos

Chromosome position to search. Required. 

=item ref

Reference allele to match. Required.

=item alt

Alt (variant) allele to match. Required.

=back

    my @matches = $cvar->searchForMatchingVariant
    (
        chrom => 1, 
        pos   => 5947454, 
        ref   => 'G',
        alt   => 'A', 
    );

=item B<grepResults> 

Search specific fields within matching lines (as retrieved by searchClinVar or searchForMatchingVariant methods). Returns an array of all matching lines.

Arguments:

=over 16

=item lines

A reference to an array of lines from the ClinVar file (e.g. as retrieved by searchClinVar or searchForMatchingVariant methods). Required.

=item column

Name of the column to search. Required.

=item term

Search term. Required. By default if this string is matched in a line anywhere within the given column it is considered a match. This behaviour can be changed using the 'exact' or 'regex' arguments.

=item exact

Set this argument to any non-zero value to require that the fields searched must exactly match your search term.

=item field_only

Set this argument to any non-zero value to only require the columns searched of any matches rather than the whole lines.

=back

    my @nphp_var = $cvar->grepResults
    (
        lines  => \@matches, 
        column => "all_traits",
        term   => "Nephronophthisis",
    );

    my @prof = $cvar->grepResults
    (
        lines  => \@matches, 
        column => "review_status",
        term   => "reviewed by professional society",
        exact  => 1,
    );

    my @changed_gly = $cvar->grepResults
    (
        lines      => \@matches, 
        column     => "hgvs_p",
        term       => qr"^NP_\d+\.\d+:p.Gly\d+[^=\s]+$",
        field_only => 1,
    );


=item B<getPathogenic> 

Returns an array of lines with a pathogenic or likely pathogenic annotation as indicated by the 'pathogenic' column of the clinvar.tsv file. Takes a reference to an array of lines (as retrieved by searchClinVar or searchForMatchingVariant methods) as the only argument. 

    my @pathogenic = $cvar->getPathogenic(\@matches);

=item B<getNonPathogenic> 

Returns an array of lines without a pathogenic or likely pathogenic annotation as indicated by the 'pathogenic' column of the clinvar.tsv file. Takes a reference to an array of lines (as retrieved by searchClinVar or searchForMatchingVariant methods) as the only argument. 

    my @nonpathogenic = $cvar->getNonPathogenic(\@matches);

=item B<getConflicted> 

Returns an array of lines indicated as 'conflicted' (that is have been described as both pathogenic and non-pathogenic previously) as indicated by the 'conflicted' column of the clinvar.tsv file. Takes a reference to an array of lines (as retrieved by searchClinVar or searchForMatchingVariant methods) as the only argument. 

    my @conflicted = $cvar->getConflicted(\@pathogenic);

=item B<getUnconflicted> 

Returns an array of lines not indicated as 'conflicted' (that is have not been described as both pathogenic and non-pathogenic previously) as indicated by the 'conflicted' column of the clinvar.tsv file. Takes a reference to an array of lines (as retrieved by searchClinVar or searchForMatchingVariant methods) as the only argument. 

    my @unconflicted = $cvar->getUnconflicted(\@pathognic);

=item B<getColumnValue>

For a given line returns the value for a specified column. The first argument is a line from the clinvar.tsv file and the second argument is the name of the column to return.

    my $traits = $cvar->getColumnValue($line, 'all_traits');

=back


=head1 AUTHOR

David A. Parry

=head1 COPYRIGHT AND LICENSE

Copyright 2015  David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.




