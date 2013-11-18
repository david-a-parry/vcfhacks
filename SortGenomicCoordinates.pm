package SortGenomicCoordinates;

use strict;
use warnings;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

use Carp;
our $AUTOLOAD;
{
 	my $_count = 0;
	my %_attrs = (
		_file => ["", "read/write"], 
		_array => ["", "read/write"],
		_col => [0, "read/write"], 
		_type => ["bed", "read/write"],
		_chrom => ["", "read/write"],
		_position => ["", "read/write"],
		_ordered => ["", "read/write"],
		_merged => ["", "read"],
		_location => ["", "read"],
		_return_type => ["match", "read/write"],
		_start_col => ["", "read/write"], 
		_stop_col => ["", "read/write"], 
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
	sub get_count{
		$_count;
	}
	sub _incr_count{
		$_count++;
	}
	sub _decr_count{
		$_count--;
	}
		
}
sub DESTROY{
	my ($self) = @_;
	$self -> _decr_count( );
}

sub new {
	my ($class, %args) = @_;
	my $self = bless { }, $class;
	foreach my $attr ($self -> _all_attrs( ) ){
		my ($arg) = ($attr =~ /^_(.*)/);
		if (exists $args{$arg}){
			$self->{$attr} = $args{$arg};
		}elsif($self->_accessible($attr, "required")){
			croak "$attr argument required";
		}else{
			$self->{$attr} = $self->_attr_default($attr);
		}
	}
	$class -> _incr_count();
	$self -> set_type($self -> {_type});
	$self -> _check_col($self -> {_col});
#	$self -> _strip_chr() if $self -> {_chrom};
	return $self;
}


sub clone {
	my ($caller, %arg) = @_;
	my $class = ref($caller);
	my $self = bless { }, $class;
	foreach my $attr ($self ->_all_attrs ( ) ){
		my ($arg) = ($attr =~ /^_(.*)/);
		if (exists $arg{$arg}){
			$self -> {$attr} = $arg{$arg};
		}else{
			$self -> {$attr} = $caller->{$attr};
		}
	}
	$self -> _incr_count();
	$self -> _check_type($self -> {_type});
	$self -> _check_col($self -> {_col});
	#$self -> _strip_chr() if $self -> {_chrom};
	return $self;
}

sub set_type{
    my ($self, $type) = @_;
    $self->_check_type($type);
    $self->{_type} = $type;
    $self->{_col} = 1 if $type ne 'custom' and $self->{_col} == 0;#if not
}

sub _check_type{
	my @valid_types = qw/bed gene regions custom/;
	my ($self, $type) = @_;
	croak "invalid type" if not grep {/$type/} @valid_types;
}
sub _check_col{
	my ($self, $col) = @_;
	croak "col argument must be an integer" if $col !~ /^-*\d+$/;
}
sub _strip_chr{
	my ($self) = @_;
	$self->{_chrom} =~ s/^chr//;
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
		#$self -> _strip_chr() if $attr eq "chrom";
		*{$AUTOLOAD} = sub { $_[0] -> {$attr} = $_[1]; return ; };
		$self -> {$attr} = $val;
		return
	}else{
		croak "Method name \"$AUTOLOAD\" not available";
	}
}
sub check_custom{
	my ($self) = @_;
	foreach ($self->{_col}, $self->{_start_col}, $self->{_stop_col}){
		croak "col, custom_start and custom_stop arguments must be defined when using \"custom\" type mode" if not defined $_;
	}
	croak "col argument $self->{_col} must be an integer" if $self->{_col} !~ /^-*\d+$/;
	croak "custom_start argument ($self->{_start_col}) must be an integer" if $self->{_start_col} !~ /^-*\d+$/;
	croak "custom_stop argument ($self->{_stop_col}) must be an integer" if $self->{_stop_col} !~ /^-*\d+$/;
}
#if a search produces a hit with -> locate method return the array of info for the relevant merged interval
#otherwise return undef if not hit
sub return_info{
	my ($self, %args) = @_;
	croak "A region must be located before info can be returned" if not defined $_[0]->{_location};
	if ($_[0]->{_location} < 0){
		return if defined(wantarray);
		carp "return_info method called in void context";
	}elsif ($self->{_return_type} eq "all"){
		return @{$_[0]->{_merged}->[$_[0]->{_location}]->{info}} if defined(wantarray);
		carp "return_info method called in void context"; #complain if not sent to either scalar or array
	}elsif ($self->{_return_type} =~ /^(match|matches)$/i){
		my @ret = ();
		my $spl = "\t";
		$spl = '[:-]' if $self->{_type} eq "regions";
		my ($start, $stop) = (1, 2); #fine for regions and bed files
		($start, $stop) = (2, 3) if $self -> {_type} eq "gene"; #gene tables have strand field between chr and start
		if ($self -> {_type} eq "custom"){
			check_custom($self);
			($start, $stop) = ($self->{_start_col}, $self->{_stop_col});
		}
		foreach my $region (@{$_[0]->{_merged}->[$_[0]->{_location}]->{info}}){
			if ($self->{_col} == 0){ #with no valid column value we'll try and find the region field ourself
				my $off = 0;
				my @split = (split /$spl/, $region);
				no warnings 'uninitialized';
				$off++ until $split[$off] =~ /chr[\dXYMUG][TLn\d_]*/i or $off > $#split;
				croak ("no valid region found in line $region ") if $off > $#split;
				$split[$off] =~ s/^chr//i;
				$self->{_chrom}  =~ s/^chr//i;
				if (uc$split[$off] eq uc$self->{_chrom} and $self->{_position} > $split[$off + $start] and $self->{_position} <= uc$split[$off + $stop]){
					push (@ret, $region);
				}
			}else{
				my @split = (split /$spl/, $region);
				$split[$self->{_col}-1] =~ s/^chr//i;
				$self->{_chrom}  =~ s/^chr//i;
				if (uc$split[$self->{_col} -1] eq uc$self->{_chrom} and $self->{_position} > $split[$self->{_col}+$start -1] and $self->{_position} <= uc$split[$self->{_col}+$stop -1]){
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

# a simple binary search checking if position lies within an interval
# returns the location in our merged array of a hit or -1 if no hit
# the returned value is also stored in $self->{_location} for access 
# by ->return_info method, so it is fine (and probably recommended)
# to call ->locate in void context
sub locate{
	my ($self, %args) = @_;
	croak "You must order and merge your reference coordinates before performing a search" if not  $self -> {_merged};
	set_chrom($self, $args{chrom}) if $args{chrom};
	croak "You must set a chromosome to search for" if not  $self -> {_chrom};
	$self -> {_chrom} =~ s/^chr//i;
	set_position($self, $args{position}) if $args{position};
	croak "You must set a position to search for" if not  defined $self -> {_position};
	$self->{_position} =~ s/\,//g; 
	croak "Position argument must be an integer" if $self->{_position} !~ /^\d+$/;
	my $l = 0;
	my $u = @{$self->{_merged}} -1;
	while ($l <= $u){
		my $i = int(($u+$l)/2);
		if((uc$self->{_chrom} lt uc$self->{_merged}->[$i]->{chrom}) or 
		( uc$self->{_chrom} eq uc$self->{_merged}->[$i]->{chrom} and $self->{_position} <= $self->{_merged}->[$i]->{start} 
		and $self->{_position} < $self->{_merged}->[$i]->{end})){
			$u = $i - 1;
		}elsif ((uc$self->{_chrom} gt uc$self->{_merged}->[$i]->{chrom}) or 
		( uc$self->{_chrom} eq uc$self->{_merged}->[$i]->{chrom} and $self->{_position} > $self->{_merged}->[$i]->{start} 
		and $self->{_position} > $self->{_merged}->[$i]->{end})){
			$l = $i + 1;
		}elsif  (uc$self->{_chrom} eq uc$self->{_merged}->[$i]->{chrom} and 
		$self->{_position} >= $self->{_merged}->[$i]->{start} and $self->{_position} <= $self->{_merged}->[$i]->{end}){
			$self->{_location} = $i;
			return $i;
		}
	}
	$self->{_location} = -1;
	return -1;
}

#usually intervals will need sorting and merging so ->prep invokes both methods
sub prep{
	my ($self, %args) = @_;
	set_file($self, $args{file}) if $args{file};
	set_array($self, $args{array}) if $args{array};
	croak "an input file or array is required before running \"order\" method" if not  $self -> {_file} and not $self -> {_array};
	set_type($self, $args{type}) if $args{type};
	set_col($self, $args{col}) if $args{col};
	$self -> order();
	$self -> merge();
}

#for effective merging and searching intervals must be sorted
sub order { #returns reference to a sorted array of lines in chromosome order
	my ($self, %args) = @_;
	my $FILE;
	set_file($self, $args{file}) if $args{file};
	set_array($self, $args{array}) if $args{array};
	set_type($self, $args{type}) if $args{type};
	set_col($self, $args{col}) if $args{col};
	croak "an input file or array is required before running \"order\" method" if not  $self -> {_file} and not $self -> {_array};
	my $spl = "\t";
	$spl = '[:-]' if $self -> {_type} eq "regions";
	my ($start, $stop) = (1, 2); #fine for regions and bed files
	($start, $stop) = (2, 3) if $self -> {_type} eq "gene"; #gene tables have strand field between chr and start
	if ($self -> {_type} eq "custom"){
		check_custom($self);
		($start, $stop) = ($self->{_start_col}, $self->{_stop_col});
	}
	if ($self->{_file}){
		if ($self->{_file} =~ /\.gz$/){
			$FILE = new IO::Uncompress::Gunzip $self->{_file} || croak("IO::Uncompress::Gunzip failed while opening $self->{_file} for sorting:\n$GunzipError\n");
		}else{
			open ($FILE, $self -> {_file} ) || croak "Can't open file " , $self->{_file} ," for sorting";
		}
		#set split parameters depending on file type
		my @sorted = ();
		if ($self -> {_col} == 0){ #with no valid column value we'll try and find the region field ourself
			chomp(@sorted = map { $_ -> [0] } 
				sort{ uc$a -> [1] cmp uc$b -> [1]  || 
					$a -> [2] <=> $b -> [2]  || 
					$a -> [3] <=> $b -> [3]  }
				map { 	my $off = 0;
					#WE SHOULD ADD A REGEX TO CHECK THAT THE REGION IS CORRECT FORMAT REALLY
					chomp (my @split = (split /$spl/));
					no warnings 'uninitialized';
					$off++ until $split[$off] =~ /^chr[\dXYMUG][TLn\d_]*/i or $off > $#split;
					$split[$off] =~ s/^chr//i unless $off > $#split;
					$off <= $#split ? [$_, $split[$off], $split[$off+$start], $split[$off+$stop]] : carp ("no valid region found in line $_") #only keep lines with a valid chr field
				} grep {!/^#/}  (<$FILE>));
		}else{
			chomp(@sorted = map { $_ -> [0] } #otherwise use user supplied column (using subtraction of 1 to convert to 0-based)
				sort{ uc$a -> [1] cmp uc$b -> [1]  || 
					$a -> [2] <=> $b -> [2]  || 
					$a -> [3] <=> $b -> [3]  }
				map { 
					chomp (my @split = (split /$spl/));
					$split[$self -> {_col} -1] =~ s/^chr//i;
					[$_, $split[$self -> {_col} -1], $split[$self -> {_col}+$start -1], $split[$self -> {_col}+$stop -1]]
				} grep {!/^#/}  (<$FILE>));
		}
		close $FILE;
		$self -> {_ordered}  = \@sorted;
		return \@sorted if defined(wantarray);
	}else{
		if ($self -> {_col} == 0){ #with no valid column value we'll try and find the region field ourself
			chomp(@{$self->{_array}} = map { $_ -> [0] } 
				sort{ uc$a -> [1] cmp uc$b -> [1]  || 
					$a -> [2] <=> $b -> [2]  || 
					$a -> [3] <=> $b -> [3]  }
				map { 	my $off = 0;
					#WE SHOULD ADD A REGEX TO CHECK THAT THE REGION IS CORRECT FORMAT REALLY
					chomp (my @split = (split /$spl/));
					no warnings 'uninitialized';
					$off++ until $split[$off] =~ /^chr[\dXYMUG][TLn\d_]*/i or $off > $#split;
					$split[$off] =~ s/^chr//i unless $off > $#split;
					$off <= $#split ? [$_, $split[$off], $split[$off+$start], $split[$off+$stop]] : carp ("no valid region found in line $_") #only keep lines with a valid chr field
				} grep {!/^#/} @{$self->{_array}} );
		}else{
			chomp(@{$self->{_array}} = map { $_ -> [0] } #otherwise use user supplied column (using subtraction of 1 to convert to 0-based)
				sort{ uc$a -> [1] cmp uc$b -> [1]  || 
					$a -> [2] <=> $b -> [2]  || 
					$a -> [3] <=> $b -> [3]  }
				map { 
					chomp (my @split = (split /$spl/));
					$split[$self -> {_col} -1] =~ s/^chr//i;
					[$_, $split[$self -> {_col} -1], $split[$self -> {_col}+$start -1], $split[$self -> {_col}+$stop -1]]
				} grep {!/^#/} @{$self->{_array}} );
		}
		$self->{_ordered}  = $self->{_array};
		return $self->{_ordered} if defined(wantarray);
	}
}
#for quick binary searching we need to merge all overlapping intervals from our sorted array
#we return an array of anonymous hashes each with keys 'chrom', 'start', 'stop' and 'info'
#original info for merged intervals is stored as an anonymous array accessed by key 'info' 
sub merge {
	my ($self, $merging) = @_; #input must be a ref to sorted array of bed regions
	if ($merging){
		carp "Warning, merged array provided independent of regions file \"$self -> {_file}\" " if $self -> {_file};
		carp "Warning, merged array overwriting existing sorted array \"$self -> {_ordered}\" " if $self -> {_ordered};
	}else{
		croak "An ordered array must be specified either by using \"->order\" method or by providing a reference to an array when calling \"->merge\" method " unless $self -> {_ordered};
		$merging = $self -> {_ordered};
	}
	ref $merging eq "ARRAY" || croak "Scalar passed to \"merge\" method must be reference to an array";
        my @merged = ();
        my ($cur_chrom,  $cur_start,  $cur_end);
        my  ($prev_chrom, $prev_start, $prev_end, $count) = (0, 0, 0, 0);
	#set split parameters depending on file type
	my $spl = "\t";
	$spl = '[:-]' if $self -> {_type} eq "regions";
	my ($start, $stop) = (1, 2); #fine for regions and bed files
	($start, $stop) = (2, 3) if $self -> {_type} eq "gene"; #gene tables have strand field between chr and start
	if ($self -> {_type} eq "custom"){
		check_custom($self);
		($start, $stop) = ($self->{_start_col}, $self->{_stop_col});
	}
	my @info;
        foreach my $table(@$merging){
		my $off = 0;
                my @regions = split(/$spl/, $table);
		if ($self -> {_col} == 0){ #with no valid column value we'll try and find the region field ourself
			no warnings 'uninitialized';
			$off++ until $regions[$off] =~ /chr*[\dXYMUG][TnL\d_]*/i or $off > $#regions;
			croak ("no valid region found in line $_") if $off > $#regions;
		}else{
			$off = $self -> {_col} -1;
		}
                ($cur_chrom = $regions[$off]) =~ s/^chr//i;
                $cur_start = $regions[$off + $start];
                $cur_end = $regions[$off + $stop];
                croak "Array passed to merge subroutine is not sorted properly "
			."- try sorting with \"SortGenomicCoordinates -> order\" subroutine or ensure chromosomes are in" 
				." ascibetical order and coordinates are in numerical order" 
					if (uc$cur_chrom lt uc$prev_chrom) or (uc$cur_chrom eq uc$prev_chrom and $cur_start < $prev_start) and $count >0;
		if ($cur_start > $cur_end){
			carp "Start value is after end value for region $cur_chrom:$cur_start-$cur_end. Swapping start and end values. ";
                	$cur_end = $regions[$off + $start];
                	$cur_start = $regions[$off + $stop];
			
		}
		
		if ((uc$cur_chrom eq uc$prev_chrom and $cur_start > $prev_end) or (uc$cur_chrom ne uc$prev_chrom)){
                        my $hash = {chrom=>$prev_chrom, start=>$prev_start, end=>$prev_end, info=>[@info]} unless ($count < 1);
                        push (@merged, $hash) unless ($count < 1);
                        $prev_chrom = $cur_chrom;
                        $prev_start = $cur_start;
                        $prev_end = $cur_end;
			@info = ();
			push (@info, $table); 
                        $count++;
                }
                elsif (uc$cur_chrom eq uc$prev_chrom and $cur_start <= $prev_end and $cur_end > $prev_end){ #if overlaps at 3' end
                        $prev_end = $cur_end;
			push (@info, $table);
                }
                elsif (uc$cur_chrom eq uc$prev_chrom and $cur_start <= $prev_end and $cur_end <= $prev_end){ #if contained within larger region
			push (@info, $table);
		}
        }
        my $hash = {chrom=>$prev_chrom, start=>$prev_start, end=>$prev_end, info=>[@info]};
        push (@merged, $hash);
        $self -> {_merged} = \@merged;
	return \@merged if defined(wantarray);
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
$obj = SortGenomicCoordinates -> new(array => \@array, type => "gene");

The user may wish to specify which column of the input file contains the start (i.e. chromosome field) of the region field, perhaps to speed up the merging and sorting of regions or if the file contains more than one set of fields that could be interpreted as regions by the program. 

$obj = SortGenomicCoordinates -> new(file => "regions_file", type => "gene", col => 3);

If type is set to "custom" the user MUST specify the chromosome column and the relative positions of the start and stop columns for regions.

$obj = SortGenomicCoordinates -> new(file => "custom_regions_file", type => "custom", col => 3, start_col => -2, stop_col => -1);

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

the format of the regions (i.e. bed, regions, gene or custom)

=item B<col>

column containing the chr field of region (or entire region field if type is region). This defaults to 0 in which case the column will be determined for each region by searching for the first "chr[0-9XYMU][0-9Tn]*" field.

=item B<ordered>

reference to sorted array. The user may wish to specify their own sorted array of regions using '$obj -> set_ordered(\@sorted_regions)' rather than setting the file and using '->order' or '->prep' methods. 

=item B<merged>

reference to merged array of hashes (keys are 'chrom', 'start', 'end' and 'info'). Read only.

=item B<chrom>

chromosome to locate, to be used in conjunction with B<position> for searching.

=item B<position>

position to locate in conjunction with B<chrom> argument.

=item B<location>

location in merged array of any hit from B<locate> method. Read only.

=item B<return_type>

Setting to give to return info. Valid types are "all" or "match". Using "match" only returns regions from the merged array that match the location searched, "all" will return all regions from within the merged array. Default is "match";

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
University of Leeds


=head1 COPYRIGHT AND LICENSE

Copyright 2011, 2012  David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut
