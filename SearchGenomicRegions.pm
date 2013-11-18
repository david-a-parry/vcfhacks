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

package SearchGenomicRegions;

use strict;
use warnings;
use Carp;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

our $AUTOLOAD;
use base ( "SortGenomicCoordinates" );
{
 	my $_count = 0;
	my %_attrs = (
		_file => ["", "read/write"], 
		_array => ["", "read/write"],
		_col => [0, "read/write"], 
		_type => ["bed", "read/write"],
		_chrom => ["", "read/write"],
		_start => ["", "read/write"],
		_stop => ["", "read/write"],
		_ordered => ["", "read/write"],
		_merged => ["", "read"],
		_search_file => ["", "read/write"],
		_search_list => ["", "read/write"],
		_list_type => ["bed", "read/write"],
		_return_type => ["match", "read/write"],
		_flanks => [0, "read/write"],
	);
	sub _all_attrs{
		keys %_attrs;
	}
	sub _accessible{
		my ($self, $attr, $mode) = @_;
		$_attrs{$attr}[1] =~ /$mode/
	}
	sub _attr_default{
		my ($self, $attr) = @_;
		$_attrs{$attr}[0];
	}
		
}
sub DESTROY{
	my ($self) = @_;
	$self -> _decr_count( );
}

#use autoload for standard ->get and ->set methods
sub AUTOLOAD{
	my ($self, $val) = @_;
	no strict 'refs';
	if ($AUTOLOAD =~ /.*::get(_\w+)/ and $self -> _accessible($1, "read")){
		my $attr = $1;
		croak "No such attribute \"$attr\"" unless exists $self->{$attr};
		*{$AUTOLOAD} = sub { return $_[0] -> {$attr} };
		return $self->{$attr};
	}elsif ($AUTOLOAD =~ /.*::set(_\w+)/ and $self -> _accessible($1, "write")){
		my $attr = $1;
		croak "No such attribute \"$attr\"" unless exists $self->{$attr};
		$self -> _check_col($val) if $attr eq "col";
		$self -> _check_type($val) if $attr eq "type";
		*{$AUTOLOAD} = sub { $_[0] -> {$attr} = $_[1]; return ; };
		$self -> {$attr} = $val;
		return
	}else{
		croak "Method name \"$AUTOLOAD\" not available";
	}
}
sub search_list{
	my ($self, %args) = @_;
	my $SEARCH; #filehandle
	croak "You must order and merge your reference coordinates before performing a search" if not  $self -> {_merged};
	set_search_file($self, $args{search_file}) if $args{search_file};
	if ($self-> {_search_file}){
		if ($self-> {_search_file} =~ /\.gz$/){
			$SEARCH = new IO::Uncompress::Gunzip $self-> {_search_file} || die("IO::Uncompress::Gunzip failed while opening $self->{_search_file} for searching:\n$GunzipError\n");
		}else{
			open ($SEARCH, $self->{_search_file}) || croak "Can't open search file \"$self->{_search_file}\" ";
		}
		my @search_array = grep {!/^#/} (<$SEARCH>);
		croak "No regions in search file" if @search_array < 1;
		$self -> {_search_list} = \@search_array;
	}
	set_return_type($self, $args{return_type}) if $args{return_type};
	set_search_list($self, $args{search_list}) if $args{search_list};
	croak "search_list or search_file argument must be set for search_list subroutine " if not $self->{_search_list};
	ref $self->{_search_list} eq "ARRAY" || die "_search_list argument must be a reference to an array ";
	my $spl = "\t";
	$spl = '[:-]' if $self->{_list_type} eq "regions";
	my ($start, $stop) = (1, 2); #fine for regions and bed files
	($start, $stop) = (2, 3) if $self -> {_list_type} eq "gene"; #gene tables have strand field between chr and start
	my $min = 3;
	$min = 4 if $self -> {_list_type} eq "gene";
	my @ret = ();
	foreach my $search (@{$self->{_search_list}}){
		my $off = 0;
		my @srch = split(/$spl/, $search);
		croak "Minimum of $min fields must be specified for list type $self -> {_list_type} " if @srch < $min;
		no warnings 'uninitialized';
		$off++ until $srch[$off] =~ /^chr[\dXYMUG][TLn\d_]*/i or $off > $#srch;
		croak ("no valid region found in line $search ") if $off > $#srch;
		set_chrom($self, $srch[$off]);
		set_start($self, $srch[$off + $start]);
		set_stop($self, $srch[$off + $stop]);
		push (@ret, $search) if ($self -> locate >= 0);
	}
	return @ret if defined(wantarray);
	carp "search_list method called in a void context";
}
sub search_against_list{
	my ($self, %args) = @_;
	croak "You must order and merge your reference coordinates before performing a search" if not  $self -> {_merged};
	set_search_file($self, $args{search_file}) if $args{search_file};
	my $SEARCH;
	if ($self-> {_search_file}){
		if ($self-> {_search_file} =~ /\.gz$/){
			$SEARCH = new IO::Uncompress::Gunzip $self-> {_search_file} || die("IO::Uncompress::Gunzip failed while opening $self->{_search_file} for searching:\n$GunzipError\n");
		}else{
			open ($SEARCH, $self->{_search_file}) || croak "Can't open search file \"$self->{_search_file}\" ";
		}
		my @search_array = grep {!/^#/} (<$SEARCH>);
		croak "No regions in search file" if @search_array < 1;
		$self -> {_search_list} = \@search_array;
	}
	set_return_type($self, $args{return_type}) if $args{return_type};
	set_search_list($self, $args{search_list}) if $args{search_list};
	croak "search_list or search_file argument must be set for search_list subroutine " if not $self->{_search_list};
	ref $self->{_search_list} eq "ARRAY" || die "_search_list argument must be a reference to an array ";
	my $spl = "\t";
	$spl = '[:-]' if $self->{_list_type} eq "regions";
	my ($start, $stop) = (1, 2); #fine for regions and bed files
	($start, $stop) = (2, 3) if $self -> {_list_type} eq "gene"; #gene tables have strand field between chr and start
	my $min = 3;
	$min = 4 if $self -> {_list_type} eq "gene";
	my @ret = ();
	foreach my $search (@{$self->{_search_list}}){
		my $off = 0;
		my @srch = split(/$spl/, $search);
		croak "Minimum of $min fields must be specified for list type $self -> {_list_type} " if @srch < $min;
		no warnings 'uninitialized';
		$off++ until $srch[$off] =~ /^chr[\dXYMUG][TLn\d_]*/i or $off > $#srch;
		croak ("no valid region found in line $search ") if $off > $#srch;
		set_chrom($self, $srch[$off]);
		set_start($self, $srch[$off + $start]);
		set_stop($self, $srch[$off + $stop]);
		$self -> locate;
		push (@ret, $self-> return_info);
	}
	my %seen = ();
	@ret = grep {! $seen{$_}++} @ret;
	return @ret if defined(wantarray);
	carp "search_against_list method called in a void context";
}

#if a search produces a hit with -> locate method return the array of info for the relevant merged interval
#otherwise return undef if not hit
sub return_info{
	my ($self, %args) = @_;
	set_return_type($self, $args{return_type}) if $args{return_type};
	croak "A region must be located before info can be returned" if not defined $self->{_location};
	if ($_[0]->{_location} < 0){
		return if defined(wantarray);
		carp "return_info method called in void context";
	}elsif ($self->{_return_type} =~ /^all$/i){
		return @{$_[0]->{_merged}->[$_[0]->{_location}]->{info}} if defined(wantarray);
		carp "return_info method called in void context"; #complain if not sent to either scalar or array
	}elsif ($self->{_return_type} =~ /^(match|matches)$/i){
		my @ret = ();
		my $spl = "\t";
		$spl = '[:-]' if $self->{_type} eq "regions";
		$self->{_chrom}  =~ s/^chr//i;
		my ($start, $stop) = (1, 2); #fine for regions and bed files
		($start, $stop) = (2, 3) if $self -> {_type} eq "gene"; #gene tables have strand field between chr and start
		foreach my $region (@{$_[0]->{_merged}->[$_[0]->{_location}]->{info}}){
			if ($self->{_col} == 0){ #with no valid column value we'll try and find the region field ourself
				my $off = 0;
				my @split = (split /$spl/, $region);
				no warnings 'uninitialized';
				$off++ until $split[$off] =~ /^chr[\dXYMUG][TLn\d_]*/i or $off > $#split;
				croak ("no valid region found in line $region ") if $off > $#split;
				$split[$off] =~ s/^chr//i;
				$self->{_chrom}  =~ s/^chr//i;
				if (uc$split[$off] eq uc$self->{_chrom}	and 
				($self->{_start} <=($split[$off+$stop] + $self->{_flanks}) and $self->{_stop} >= ($split[$off+$start] - $self->{_flanks}))){ 
					push (@ret, $region);
				}
			}else{
				my @split = (split /$spl/, $region);
				$split[$self->{_col}-1] =~ s/^chr//i;
                                $self->{_chrom}  =~ s/^chr//i;
				if (uc$split[$self->{_col} -1] eq uc$self->{_chrom}	and 
				(($self->{_start} <=($split[$self->{_col}-1+$stop] + $self->{_flanks})) and ($self->{_stop} >= ($split[$self->{_col}-1+$start] - $self->{_flanks}) ))){
					push (@ret, $region);
				}
			}
		}
		croak "Error, region was located but no matching regions for current chromosome and position values exist " if (not @ret);
		return @ret if defined(wantarray); #allow scalar or array for return value
		carp "return_info method called in void context"; #complain if not sent to either scalar or array
	}else{
		croak "unrecognised return type given to $self -> return_info method ";
	}
}

# a simple binary search testing whether a position lies within an interval from our merged array
# returns the location in our merged array of a hit or -1 if no hit
# the returned value is also stored in $self->{_location} for access 
# by ->return_info method, so it is fine (and probably recommended)
# to call ->locate in void context
sub locate{
	my ($self, %args) = @_;
	croak "You must order and merge your reference coordinates before performing a search" if not  $self -> {_merged};
	my ($search_chr, $search_start, $search_stop);
	set_chrom($self, $args{chrom}) if $args{chrom};
	croak "chrom argument is required for locate_region method " if not $self->{_chrom};
	$self->{_chrom}  =~ s/^chr//i;
	set_start($self, $args{start}) if $args{start};
	croak "start argument is required for locate_region method " if not defined $self->{_start};
	$self->{_start} =~ s/\,//g;
	croak "Start argument must be an integer" if $self->{_start} !~ /^\d+$/;
	set_stop($self, $args{stop}) if $args{stop};
	croak "stop argument is required for locate_region method " if not defined $self->{_stop};
	$self->{_stop} =~ s/\,//g;
	croak "Stop argument must be an integer" if $self->{_stop} !~ /^\d+$/;
	$self->{_flanks} =~ s/\,//g; 
	croak "Flanks argument must be an integer" if $self->{_flanks} !~ /^\d+$/;
	my $l = 0;
	my $u = @{$self->{_merged}} -1;
	while ($l <= $u){
		my $i = int(($u+$l)/2);
		if((uc$self->{_chrom} lt uc$self->{_merged}->[$i]->{chrom}) or 
		( uc$self->{_chrom} eq uc$self->{_merged}->[$i]->{chrom} and $self->{_stop}  <= ($self->{_merged}->[$i]->{start} - $self->{_flanks})  )){
			$u = $i - 1;
		}elsif ((uc$self->{_chrom} gt uc$self->{_merged}->[$i]->{chrom}) or 
		( uc$self->{_chrom} eq uc$self->{_merged}->[$i]->{chrom} and $self->{_start}  >= ($self->{_merged}->[$i]->{end} + $self->{_flanks})  )){
			$l = $i + 1;
		}else{
			$self->{_location} = $i;
			return $i;
		}
	}
	$self->{_location} = -1;
	return -1;
}

1;

=head1 NAME

SearchGenomicRegions.pm - sort, merge  and search genomic regions for overlapping genomic regions.

=head1 SYNOPSIS

$obj = SearchGenomicRegions -> new( );
#(initialise object)

$obj = SearchGenomicRegions -> new(file => "regions_file", type => "bed");
#(initialise object with file of regions, specifying format)

$obj -> prep;
#(sort and merge regions keeping remaining info from each region)

$obj -> locate(chrom => "chr1", start => 1000000, stop => 2000000);
#(find location of regions containing chr1:1000000  within sorted and merged array contained within $obj)

@found_info = $obj -> return_info;
#(return array of all lines (regions) found by locate function)

@matches = $obj-> search_list(search_file => "file.bed", list_type => "bed");
#specify and search a list of regions to compare search your regions of interest with. This method returns matching lines of the former list.

@matches = $obj-> search_list(search_against_file => "file.bed", list_type => "bed");
#specify and search a list of regions to compare search your regions of interest with. This method returns matching lines of the latter list (i.e. your regions of interest).

=head1 DESCRIPTION

=head2 Overview

This modeule inherits from SortGenomicCoordinates.pm to read bed files, files of regions (chr1:100000-200000) or gene files from UCSC genome browser (e.g. refGene.txt) into memory, sort and merge the regions within in order to perform rapid searches for regions overlapping with user-specified regions. The key diference between this module and SortGenomicCoordinates.pm is that instead of searching for regions that contain a discrete user-specified location this module is used to search for regions that overlap with each other. Lists of regions can be used to search either returning matches from those regions ($obj->search_list) or from the reference regions ($obj->search_against_list).


=head2 Constructor and initialization

Minimal requirement for initialisation is to call the 'new' method without arguments. File and other arguments can be set in a separate step or a reference to a manually sorted array can be provided instead.

$obj = SearchGenomicRegions -> new; 

Usually initialisation will be performed specifying arguments for file or array and region type at least:

$obj = SearchGenomicRegions -> new(file => "regions_file", type => "gene");
$obj = SearchGenomicRegions -> new(array => \@array, type => "gene");

The user may wish to specify which column of the input file contains the start (i.e. chromosome field) of the region field, perhaps to speed up the merging and sorting of regions or if the file contains more than one set of fields that could be interpreted as regions by the program. 

$obj = SearchGenomicRegions -> new(file => "regions_file", type => "gene", col => 3);

Objects can also be created using the 'clone' method, although this is seldom likely to be a good idea.

$another_obj = $obj -> clone(file => "another_regions_file", type => "regions", col => 0);


=head2 Class and object methods

DOCUMENTATION STILL UNDER CONSTRUCTION

=head3 Accessors and mutators

The following features can be accessed using get_[feature] or set using set_[feature], substituting [feature] for the feature of interest. 

=over 12

=item B<file>

filename containing regions of interest in bed, region or ucsc gene format.

=item B<array>

reference to an array containing regions of interest in bed, region or ucsc gene format.

=item B<type>

the format of the regions (i.e. bed, regions or gene)

=item B<col>

column containing the chr field of region (or entire region field if type is region). This defaults to 0 in which case the column will be determined for each region by searching for the first "chr[0-9XYMU][0-9Tn]*" field.

=item B<ordered>

reference to sorted array. The user may wish to specify their own sorted array of regions using '$obj -> set_ordered(\@sorted_regions)' rather than setting the file and using '->order' or '->prep' methods. 

=item B<merged>

reference to merged array of hashes (keys are 'chrom', 'start', 'end' and 'info'). Read only.

=item B<location>

location in merged array of any hit from B<locate> method. Read only.

=item B<chrom>

chromosome of region to search against.

=item B<start>

5' position of a region to search against. 

=item B<stop>

3' position of a region to search against. 

=item B<search_file>

file of regions to search with (i.e. each region will be used to search the merged array created previously). 

=item B<search_list>

array of regions to search with (to be used instead of 'search_file').

=item B<flanks>

number of nucleotides to add to either end of a region when searching regions

=back

=head3 Methods

=over 12

=item B<prep>

Sort and merge regions to prepare for binary searching. This method simply invokes the B<order> and B<merge> methods. 'file', 'col' and 'type' arguments can be provided if not already specified as desired.

$obj -> prep( );

$obj -> prep(file => "regions.bed", type => "bed", col => 1);


=item B<order>

Sort regions from file putting chromosomes in ascibetical order and positions in numeric order. Header lines (starting "#") are skipped. 'file', 'col' and 'type' arguments can be specified here too.

$obj -> order( );

$obj -> order(file => "regions.bed", type => "bed", col => 1);

=item B<merge>

Merge overlapping regions of sorted array.  A reference to an ordered array can be specified here, default and recommended usage is to use the sorted array already created by 'order' method.

$obj -> merge( );

$obj -> merge(\@sorted_array);

=item B<locate>

Search merged and sorted array for a specific chromosome start position and end position. Region arguments 'chrom', 'start' and 'stop' must be specified.

$obj -> locate_region(chrom => "chr1", start => 10000000, stop => 12000000);

=item B<return_info>

Return array of matching regions and accompanying info if a region was found following 'locate' method. Returns nothing if no hit was found. Calling with "all" argument returns all regions from the matching merged region regardless of whether each individual region matches, calling without arguments or with "matches" argument returns only matching regions.

@info = $obj -> return_info;

@all_info = $obj -> return_info("all");

=item B<search_list>

Use this method to return matching lines from search_list. The arguments 'search_file' or 'search_list' can be specified here. This method will error if 'search_file' or 'search_list' have not been specified either prior to or when calling this method. 

@matches = $obj-> search_list;

@matches = $obj-> search_list(search_file => "file.bed", list_type => "bed");

=item B<search_against_list>

Use this method to return lines from regions of interest (i.e. those in merged array) that overlap those in search_list. The arguments 'search_file' or 'search_list' can be specified here. This method will error if 'search_file' or 'search_list' have not been specified either prior to or when calling this method. 

@matches = $obj-> search_against_list;

@matches = $obj-> search_against_list(search_file => "file.bed", list_type => "bed");

=back

=head1 AUTHOR

David A. Parry
University of Leeds


=head1 COPYRIGHT AND LICENSE

Copyright 2011, 2012  David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
 implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU Gen
eral Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut


