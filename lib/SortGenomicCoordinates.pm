package SortGenomicCoordinates;

use strict;
use warnings;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use List::Util qw(max);
use Carp;
our $AUTOLOAD;
{
    my $_count = 0;
    my %_attrs = (
        _file         => [ "",      "read/write" ],
        _array        => [ "",      "read/write" ],
        _col          => [ 0,       "read/write" ],
        _start_col    => [ 1,       "read/write" ],
        _stop_col     => [ 2,       "read/write" ],
        _type         => [ "bed",   "read/write" ],
        _delimiter    => [ "\t",    "read/write" ],
        _contig_order => [ "",      "read/write" ],
        _ordered      => [ "",      "read/write" ],
        _merged       => [ "",      "read" ],
        _chrom        => [ "",      "read/write" ],
        _position     => [ "",      "read/write" ],
        _location     => [ "",      "read" ],
        _return_type  => [ "match", "read/write" ],
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
    $class->_incr_count();
    if ($args{type}){
        $self->set_type( %args); 
    }else{
        $self->set_type( type => $self->{_type} );
    }
    #    $self -> _strip_chr() if $self -> {_chrom};
    return $self;
}

sub clone {
    my ( $caller, %arg ) = @_;
    my $class = ref($caller);
    my $self = bless {}, $class;
    foreach my $attr ( $self->_all_attrs() ) {
        my ($arg) = ( $attr =~ /^_(.*)/ );
        if ( exists $arg{$arg} ) {
            $self->{$attr} = $arg{$arg};
        }
        else {
            $self->{$attr} = $caller->{$attr};
        }
    }
    $self->_incr_count();
    #$self -> _strip_chr() if $self -> {_chrom};
    return $self;
}

sub set_contig_order{
    my ($self, $contig_order) = @_;
    if (ref $contig_order ne 'HASH'){
        croak "Argument passed to set_contig_order function must be a HASH reference ";
    }
    $self->{_contig_order} = $contig_order;
}

sub set_type {
    my ( $self, %args ) = @_;
    croak "'type' argument is required for set_type method " if not defined $args{type};
    $self->{_type} = $args{type};
    $self->{_col} = 0;
    if ($args{type} eq 'bed'){
        $self->{_start_col} = 1;
        $self->{_stop_col} = 2;
    }elsif($args{type} eq 'gene'){
        $self->{_start_col} = 2;
        $self->{_stop_col} = 3;
    }elsif($args{type} eq 'regions'){
        $self->{_start_col} = 1;
        $self->{_stop_col} = 2;
        $self->{_delimiter} = '[:-]';
    }elsif($args{type} eq 'custom'){
        foreach my $r (qw(col start_col stop_col)){
            croak "'$r' argument is required when using set_type on a custom format "
                if not exists ($args{$r});
            $self->_check_col($args{$r});
        }
        $self->{_col} = $args{col};
        $self->{_start_col} = $args{start_col};
        $self->{_stop_col} = $args{stop_col};
        $self->{_delimiter} = $args{delimiter} if defined $args{delimiter};
    }else{
        croak "invalid type passed to set_type method ";
    }
}

sub _check_col {
    my ( $self, $col ) = @_;
    croak "column arguments must be integers" if $col !~ /^-*\d+$/;
}

sub _strip_chr {
    my ($self) = @_;
    $self->{_chrom} =~ s/^chr//;
}

#use autoload for standard ->get and ->set methods
sub AUTOLOAD {
    my ( $self, $val ) = @_;
    no strict 'refs';
    if ( $AUTOLOAD =~ /.*::get(_\w+)/ and $self->_accessible( $1, "read" ) ) {
        my $attr = $1;
        croak "No such attribute \"$attr\"" unless exists $self->{$attr};
        *{$AUTOLOAD} = sub { return $_[0]->{$attr} };
        return $self->{$attr};
    }
    elsif ( $AUTOLOAD =~ /.*::set(_\w+)/ and $self->_accessible( $1, "write" ) )
    {
        my $attr = $1;
        croak "No such attribute \"$attr\"" unless exists $self->{$attr};
        #$self -> _strip_chr() if $attr eq "chrom";
        *{$AUTOLOAD} = sub { $_[0]->{$attr} = $_[1]; return; };
        $self->{$attr} = $val;
        return;
    }
    else {
        croak "Method name \"$AUTOLOAD\" not available";
    }
}

#if a search produces a hit with -> locate method return the array of info for the relevant merged interval
#otherwise return undef if not hit
sub return_info {
    my ( $self, %args ) = @_;
    croak "A region must be located before info can be returned"
      if not defined $_[0]->{_location};
    if ( $_[0]->{_location} < 0 ) {
        return if defined(wantarray);
        carp "return_info method called in void context";
    }
    elsif ( $self->{_return_type} eq "all" ) {
        return @{ $_[0]->{_merged}->[ $_[0]->{_location} ]->{info} }
          if defined(wantarray);
        carp "return_info method called in void context"
          ;    #complain if not sent to either scalar or array
    }
    elsif ( $self->{_return_type} =~ /^(match|matches)$/i ) {
        my @ret = ();
        my $spl = $self->{_delimiter};
        foreach
          my $region ( @{ $_[0]->{_merged}->[ $_[0]->{_location} ]->{info} } )
        {
            my @split = ( split /$spl/, $region );
            if ( $split[ $self->{_col} ] eq $self->{_chrom}
                and $self->{_position} >
                $split[ $self->{_start_col} ]
                and $self->{_position} <=
                uc $split[ $self->{_stop_col} ] )
            {
                push( @ret, $region );
            }
        }
        croak
"Error, region was located but no matching regions for current chromosome and position values exist "
          if ( not @ret );
        return @ret
          if defined(wantarray);   
        carp "return_info method called in void context"
          ;    
    }
    else {
        croak "unrecognised return type given to $self -> return_info method ";
    }
}

# a simple binary search checking if position lies within an interval
# returns the location in our merged array of a hit or -1 if no hit
# the returned value is also stored in $self->{_location} for access
# by ->return_info method, so it is fine (and probably recommended)
# to call ->locate in void context
sub locate {
    my ( $self, %args ) = @_;
    croak
"You must order and merge your reference coordinates before performing a search"
      if not $self->{_merged};
    set_chrom( $self, $args{chrom} ) if $args{chrom};
    croak "You must set a chromosome to search for" if not $self->{_chrom};
    $self->{_chrom} =~ s/^chr//i;
    set_position( $self, $args{position} ) if $args{position};
    croak "You must set a position to search for"
      if not defined $self->{_position};
    $self->{_position} =~ s/\,//g;
    croak "Position argument must be an integer"
      if $self->{_position} !~ /^\d+$/;
    my $l = 0;
    my $u = @{ $self->{_merged} } - 1;

    while ( $l <= $u ) {
        my $i = int( ( $u + $l ) / 2 );
        if (
            ( $self->_chromIsLessThanChrom($self->{_chrom}, $self->{_merged}->[$i]->{chrom}) )
            or (    $self->{_chrom} eq $self->{_merged}->[$i]->{chrom}
                and $self->{_position} <= $self->{_merged}->[$i]->{start}
                and $self->{_position} < $self->{_merged}->[$i]->{end} )
          )
        {
            $u = $i - 1;
        }
        elsif (
            ( $self->_chromIsLessThanChrom($self->{_merged}->[$i]->{chrom}, $self->{_chrom}) )
            or (  $self->{_chrom} eq $self->{_merged}->[$i]->{chrom}
                and $self->{_position} > $self->{_merged}->[$i]->{start}
                and $self->{_position} > $self->{_merged}->[$i]->{end} )
          )
        {
            $l = $i + 1;
        }
        elsif ( $self->{_chrom} eq $self->{_merged}->[$i]->{chrom}
            and $self->{_position} >= $self->{_merged}->[$i]->{start}
            and $self->{_position} <= $self->{_merged}->[$i]->{end} )
        {
            $self->{_location} = $i;
            return $i;
        }else{
            croak "Logic error - please report this bug! ";
        }
    }
    $self->{_location} = -1;
    return -1;
}

sub _chromIsLessThanChrom{
    my ($self, $a_chrom, $b_chrom) = @_;
    if ($self->{_contig_order}){
        if (not exists $self->{_contig_order}->{$a_chrom}){
            carp "$a_chrom does not exist in contig_order passed to SortGenomicCoordinates ";
            return -1 if exists $self->{_contig_order}->{$b_chrom};
            return 1 if $a_chrom  lt $b_chrom ;
            return 0;
        }
        if (not exists $self->{_contig_order}->{$b_chrom}){
            carp "$b_chrom does not exist in contig_order passed to SortGenomicCoordinates ";
            return 1;
        }
        return 1 if $self->{_contig_order}->{$a_chrom} < $self->{_contig_order}->{$b_chrom} ;
    }else{
        return 1 if $a_chrom  lt $b_chrom ;
    }
    return 0;
}

#usually intervals will need sorting and merging so ->prep invokes both methods
sub prep {
    my ( $self, %args ) = @_;
    set_file( $self, $args{file} ) if $args{file};
    set_array( $self, $args{array} ) if $args{array};
    croak "an input file or array is required before running \"order\" method"
      if not $self->{_file} and not $self->{_array};
    set_type( $self, $args{type} ) if $args{type};
    set_col( $self, $args{col} ) if $args{col};
    $self->order();
    my $merged = $self->merge();
    return $merged if defined wantarray;
}

#for effective merging and searching intervals must be sorted
sub order {    #returns reference to a sorted array of lines in chromosome order
    my ( $self, %args ) = @_;
    $self->set_file( $args{file} ) if $args{file};
    $self->set_array($args{array} ) if $args{array};
    $self->set_type( $args{type} ) if $args{type};
    $self->set_col( $args{col} ) if $args{col};
    $self->set_contig_order( $args{contig_order} ) if $args{contig_order};
    croak "an input file or array is required before running \"order\" method"
      if not $self->{_file} and not $self->{_array};
    my $spl = $self->{_delimiter}; 
    my ($start, $stop ) = ( $self->{_start_col}, $self->{_stop_col} );
    my $FILE;
    if ( $self->{_file} ) {
        if ( $self->{_file} =~ /\.gz$/ ) {
            $FILE = new IO::Uncompress::Gunzip $self->{_file}
              || croak(
"IO::Uncompress::Gunzip failed while opening $self->{_file} for sorting:\n$GunzipError\n"
              );
        }
        else {
            open( $FILE, $self->{_file} ) || croak "Can't open file ",
              $self->{_file}, " for sorting";
        }
        @{ $self->{_array} } = grep { !/^#/ } <$FILE>;
        close $FILE;
    }
    @{ $self->{_array} } = sort {_by_coordinate($self, $a, $b)} @{ $self->{_array} };
    $self->{_ordered} = $self->{_array};
    return $self->{_ordered} if defined(wantarray);
}

sub _by_coordinate{
    my ($self, $a, $b) = @_;
    my @a_split = split(/$self->{_delimiter}/, $a);
    my @b_split = split(/$self->{_delimiter}/, $b);
    my $max = max($self->{_col}, $self->{_start_col}, $self->{_stop_col});
    if ($max >= @a_split){
        return 0 if $max >= @b_split;
        return -1;
    }elsif ($max >= @b_split){
        return 1;
    }
    my $a_chrom = $a_split[$self->{_col}];
    my $a_start = $a_split[$self->{_start_col}];
    my $a_end = $a_split[$self->{_stop_col}];
    my $b_chrom = $b_split[$self->{_col}];
    my $b_start = $b_split[$self->{_start_col}];
    my $b_end = $b_split[$self->{_stop_col}];
    if ($self->{_contig_order}){
        if (not exists $self->{_contig_order}->{$a_chrom}){
            carp "$a_chrom does not exist in contig_order passed to SortGenomicCoordinates ";
            if (exists $self->{_contig_order}->{$b_chrom}){
                return 1;
            }else{
                return (
                    $a_chrom  cmp $b_chrom || 
                    $a_start  <=> $b_start ||
                    $a_end    <=> $b_end 
                    );
            }
        }
        if (not exists $self->{_contig_order}->{$b_chrom}){
            carp "$b_chrom does not exist in contig_order passed to SortGenomicCoordinates ";
            return -1;
        }
        return (
            $self->{_contig_order}->{$a_chrom} <=> $self->{_contig_order}->{$b_chrom} ||
                                      $a_start <=> $b_start                           ||
                                        $a_end <=> $b_end 
            );
    }else{
        return (
            $a_chrom  cmp $b_chrom || 
            $a_start  <=> $b_start ||
            $a_end    <=> $b_end 
            );
    }
}

#for quick binary searching we need to merge all overlapping intervals from our sorted array
#we return an array of anonymous hashes each with keys 'chrom', 'start', 'stop' and 'info'
#original info for merged intervals is stored as an anonymous array accessed by key 'info'
sub merge {
    my ( $self, $merging ) =
      @_;    #input must be a ref to sorted array of bed regions
    if ($merging) {
        carp
"Warning, merged array provided independent of regions file \"$self -> {_file}\" "
          if $self->{_file};
        carp
"Warning, merged array overwriting existing sorted array \"$self -> {_ordered}\" "
          if $self->{_ordered};
    }
    else {
        croak
"An ordered array must be specified either by using \"->order\" method or by providing a reference to an array when calling \"->merge\" method "
          unless $self->{_ordered};
        $merging = $self->{_ordered};
    }
    ref $merging eq "ARRAY"
      || croak
      "Scalar passed to \"merge\" method must be reference to an array";
    my @merged = ();
    my %current_region_hash = ();
    my $prev_chrom;
    my $prev_start;
    foreach my $line (@$merging) {
        my @regions = split( /$self->{_delimiter}/, $line);
        my $cur_start = $regions[ $self->{_start_col} ];
        my $cur_end   = $regions[ $self->{_stop_col} ];
        my $cur_chrom = $regions[ $self->{_col} ];
        if (%current_region_hash){
            croak "Array passed to merge subroutine is not sorted properly "
            . "- try sorting with \"SortGenomicCoordinates -> order\" subroutine or ensure chromosomes are in"
            . " contig/ascibetical order and coordinates are in numerical order"
                if (not $self->_check_chrom_order($prev_chrom, $cur_chrom)) or ( $prev_chrom eq $cur_chrom and $cur_start < $prev_start);

            if (   ( $cur_chrom eq $current_region_hash{chrom} and $cur_start > $current_region_hash{end} )
                or ( $cur_chrom ne $current_region_hash{chrom} ) )
            {
                push( @merged, {%current_region_hash} ) ;
                %current_region_hash = (
                    chrom => $cur_chrom,
                    start => $cur_start,
                    end   => $cur_end,
                    info  => [$line],
                );
            }

            elsif ( uc $cur_chrom eq uc $current_region_hash{chrom}
                and $cur_start <= $current_region_hash{end}
                and $cur_end > $current_region_hash{end} )
            {    #if overlaps at 3' end
                $current_region_hash{end} = $cur_end;
                push( @{$current_region_hash{info}}, $line );
            }
            elsif ( uc $cur_chrom eq uc $current_region_hash{chrom}
                and $cur_start <= $current_region_hash{end}
                and $cur_end <= $current_region_hash{end} )
            {    #if contained within larger region
                push( @{$current_region_hash{info}}, $line );
            }
        }else{
            %current_region_hash = (
                    chrom => $cur_chrom,
                    start => $cur_start,
                    end   => $cur_end,
                    info  => [$line],
            );
        }
        $prev_start = $cur_start;
        $prev_chrom = $cur_chrom;
    }
    push( @merged, {%current_region_hash} ) ;
    $self->{_merged} = \@merged;
    return \@merged if defined(wantarray);
}

sub _check_chrom_order{
    my ($self, $first_chrom, $next_chrom) = @_;
    if ($self->{_contig_order}){
        if (not exists $self->{_contig_order}->{$first_chrom}){
            if (not exists $self->{_contig_order}->{$next_chrom}){
                return 0 if $first_chrom gt $next_chrom;
                return 1;
            }else{
                return 1;
            }
        }elsif (not exists $self->{_contig_order}->{$next_chrom}){
            return 1;
        }
        return 0 if $self->{_contig_order}->{$first_chrom} > $self->{_contig_order}->{$next_chrom};
        return 1;
    }else{
        return 0 if $first_chrom gt $next_chrom;
        return 1;
    }
}

1;

=head1 NAME

SortGenomicCoordinates.pm - sort, merge  and search genomic regions.

=head1 SYNOPSIS

$obj = SortGenomicCoordinates -> new( );
#(initialise object)

$obj = SortGenomicCoordinates -> new(file => "regions_file", type => "bed");
#(initialise object with file of regions, specifying format)

$obj -> prep;
#(sort and merge regions keeping reamining info from each region)

$obj -> locate(chrom => "chr1", position => 1000000);
#(find location of regions containing chr1:1000000  within sorted and merged array contained within $obj)

@found_info = $obj -> return_info;
#(return array of all lines (regions) found by locate function)

=head1 DESCRIPTION

=head2 Overview

This modeule can be used to read bed files, files of regions (chr1:100000-200000) or gene files from UCSC genome browser (e.g. refGene.txt) into memory, sort and merge the regions within in order to perform rapid searches for regions containing user-specified locations. Overlapping regions are merged and ordered to allow binary searching, matching regions or all regions contained within a merged region can be returned for further manipulation. Usually an object will be initialised and regions prepared once and locations searched multiple times on the same object. A reference to an array of regions can be provided instead of a filename if desired. Custom formats may also be specified (see below). 

=head2 Constructor and initialization

Minimal requirement for initialisation is to call the 'new' method without arguments. File and other arguments can be set in a separate step or a reference to a manually sorted array can be provided instead.

$obj = SortGenomicCoordinates -> new; 

Usually initialisation will be performed specifying arguments for file or array and region type at least:

$obj = SortGenomicCoordinates -> new(file => "regions_file", type => "gene");

$obj = SortGenomicCoordinates -> new(array => \@array, type => "bed");

Specifying 'type' is a shortcut to specify what each column represents. Acceptable values are 'bed', 'gene', 'region' or 'custom'. See below for details on formats. The default type is 'bed'. If type is set to "custom" the user MUST specify the chromosome, start coordinate and stop coordinate columns. 

$obj = SortGenomicCoordinates -> new(file => "regions.txt", type => "custom", col => 1, start_col => 3, stop_col => 4);


=head2 Class and object methods

=head3 Accessors and mutators

The following features can be accessed using get_[feature] or set using set_[feature], substituting [feature] for the feature of interest. 

=over 12

=item B<file>

filename containing regions of interest in bed, region or ucsc gene format.

=item B<array>

reference to an array containing regions of interest in bed, region or ucsc gene format.

=item B<type>

the format of the regions (i.e. bed, regions, gene or custom). Bed format expects columns separated by tabs with the chromosome as the first column, start coordinate as the second column and end cooridinate as the third column. For region formats regions are expected in the format 'chromosome:start-end'. Gene formats have chromosome as the first column, start coordinate as the third column and end cooridinate as the fourth column. For custom formats the chromosome column must be specified with the 'col' feature, and start and stop coordinates by the 'start_col' and 'stop_col' features. The column delimiter can also be specified using the 'delimiter' feature.

=item B<col>

0-based column containing the chr field of region. This defaults to 0 (i.e. the first column).

=item B<start_col>

0-based column containing the start coordinate of region. This is usually determined by the region format but must be manually specified when specifying regions of type 'custom'. Defaults to 1 (i.e. the second column).

=item B<stop_col>

0-based column containing the stop coordinate of region. This is usually determined by the region format but must be manually specified when specifying regions of type 'custom'. This can be the same as 'start_col' (e.g. for VCF files). Defaults to 2 (i.e. the third column).

=item B<delimiter>

The character(s) separating each column. Defaults to tab characters unless 'type' is set to 'region' in which case ':' and '-" characters are used to separate chromosome and coordinates.

=item B<ordered>

reference to sorted array. The user may wish to specify their own sorted array of regions using '$obj -> set_ordered(\@sorted_regions)' rather than setting the file and using '->order' or '->prep' methods. 

=item B<merged>

reference to merged array of hashes (keys are 'chrom', 'start', 'end' and 'info'). Read only.

=item B<chrom>

chromosome to locate, to be used in conjunction with B<position> for searching.

=item B<position>

position to locate in conjunction with B<chrom> argument.

=item B<location>

index in merged array of any hit from B<locate> method. If not match was found this value will be -1 so make sure any code checks this value is > -1 if using this value directly. Read only. 

=item B<return_type>

Setting to give to return info. Valid types are "all" or "match". Using "match" only returns regions from the merged array that match the location searched, "all" will return all regions from within the merged array. Default is "match";

=item B<contig_order>

A reference to a hash where keys are the chromosome names and values are the relative order of each chromosome. By default, chromosomes are sorted in ascibetical order, but by passing a hash reference to the contig_order argument chromosomes can be sorted in any order desired. Methods will croak if this argument is used and a chromosome not present in the contig_order hash reference is encountered.

=back

=head3 Methods

=over 12

=item B<prep>

Sort and merge regions to prepare for binary searching. This method simply invokes the B<order> and B<merge> methods. 'file', 'col' and 'type' arguments can be provided if not already specified as desired.

$obj -> prep( );

$obj -> prep(file => "regions.bed", type => "bed");


=item B<order>

Sort regions from file putting chromosomes in ascibetical order and positions in numeric order. Header lines (starting "#") are skipped. 'file', 'col' and 'type' arguments can be specified here too.

$obj -> order( );

$obj -> order(file => "regions.bed", type => "bed");

=item B<merge>

Merge overlapping regions of sorted array.  A reference to an ordered array can be specified here, default and recommended usage is to use the sorted array already created by 'order' method.

$obj -> merge( );

$obj -> merge(\@sorted_array);

=item B<locate>

Search merged and sorted array for a specific chromosome and position. 'chrom' and 'position' arguments can be specified here.

$obj -> locate( );

$obj -> locate(chrom => "chr1", position => 10000000);


=item B<return_info>

Return array of matching regions and accompanying info if a region was found following 'locate' method. Returns nothing if no hit was found. Calling with return type "all" argument returns all regions from the matching merged region regardless of whether each individual region matches, otherwise argument returns only matching regions.

@info = $obj -> return_info;

@all_info = $obj -> return_info(return_type => "all");

=back

=head1 AUTHOR

David A. Parry

=head1 COPYRIGHT AND LICENSE

Copyright 2011, 2012, 2014, 2015  David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut
