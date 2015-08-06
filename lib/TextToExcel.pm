#    Copyright 2014, 2015  David A. Parry
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.


package TextToExcel;
use strict;
use warnings;
use Data::Dumper;
use Excel::Writer::XLSX;
use Excel::Writer::XLSX::Utility;
use Carp;
require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw();
our $AUTOLOAD;
{
    my $_count = 0;
    my %_attrs = (
        _file => ["", "read"],
        _rows => [ [0] , "read"],
        _columns => [ [0] , "read"],
        _worksheets => [ []  , "read/write"],
        _workbook => ["", "read"],
        _name => ["Sheet 1", "read"], # name of first worksheet
        _delimiter => ["\t", "read/write"],
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
        $self->{_workbook}->close() if defined $self->{_workbook};
        $self -> _decr_count( );
}

sub new {
    #requires either 'file' or 'workbook' to be passed to constructor

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
    if ($self->{_file}){
        $self->{_workbook}  = Excel::Writer::XLSX->new($self->{_file});
    }elsif(not $self->{_workbook}){
        croak "Either file or workbook items must be passed to constructor ";
    }
    $self->{_worksheets}->[0] = $self->{_workbook}->add_worksheet($self->{_name});
    $self->{_std_format} = $self->{_workbook}->add_format();
    return $self;
}

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
        *{$AUTOLOAD} = sub { $_[0] -> {$attr} = $_[1]; return ; };
        $self -> {$attr} = $val;
        return
    }else{
        croak "Method name \"$AUTOLOAD\" not available";
    }
}




sub _checkWriteArguments{
    my ($self, $method, %args) = @_;
    my @required = qw (
                    line
                    );
    my @valid = qw (
                    line
                    row
                    column
                    vcf_object
                    worksheet
                    format
                    preceding
                    succeeding
                    delimiter
                    ) ;
    #line is either a ref to an array (one value per cell)
    # or a string with each value separated by a delimiter (tab = default)
    #delimiter is a string delimiter used if if line or 
    #row is row to starting writing to on worksheet (0-based)
    #column is starting column to write to on  worksheet (0-based)
    #preceding is a ref to a 2D array of columns and rows to write before line, the line will span these rows if more than one
    #or preceding can be a ref to an array of strings to be split using delimiter
    #succeeding is a ref to a 2D array of columns and rows to write after VCF line, the line will span these rows if more than one
    #or succceding can be a ref to an array of strings to be split using delimiter
    #format must be a format value created using Excel::Writer::XLSX
    foreach my $arg (@required){
        croak "$arg is required by $method " if not exists $args{$arg};
    }
    foreach my $arg (keys %args){
        croak "Argument '$arg' passed to $method is not recognised " if not grep {$_ eq $arg} @valid;
    }
    if (exists $args{preceding}){
        if (not ref $args{preceding} eq 'ARRAY'){
            croak "value for 'preceding' passed to $method must be an ARRAY reference " ;
        }
    }

    if (exists $args{succeeding}){
        if (not ref $args{succeeding} eq 'ARRAY'){
            croak "value for 'preceding' passed to $method must be an ARRAY reference " ;
        }
    }
}

sub addWorksheet{
#adds worksheet and returns its index from @{$self->{_worksheets}}
    my ($self, $name) = @_;
    my $sheet = $self->{_workbook}->add_worksheet($name);
    push @{$self->{_worksheets}}, $sheet;
    push @{$self->{_rows}}, 0;
    push @{$self->{_columns}}, 0;
    return @{$self->{_worksheets}} - 1;
}


sub createFormat{
    #creates a format for use with workbook
    #returns a format that can be passed to the writeLine 
    #function using 
    my ($self, %args) = @_;
    my $format = $self->{_workbook}->add_format(%args);
    return $format;
}

sub writeLine{
    #returns next row number of worksheet
    my ($self, %args) = @_;
    $self->_checkWriteArguments("writeLine", %args);
    return $self->_writeLineToExcel(%args);
}

sub _writeLineToExcel{
    #returns the number of the row in the worksheet after writing this line
    #usually this will just be used by writeLine function 
    my ($self, %args) = @_;
    $self->_checkWriteArguments("_writeLineToExcel", %args);

    my $columns_before = 0;
    my $additional_rows = 0;
    my @fields = ();
    my $sheet_index = 0;
    if (defined $args{worksheet}){
        $sheet_index = $args{worksheet};
    }
    my $worksheet = $self->{_worksheets}->[$sheet_index];
    my $row = $self->{_rows}->[$sheet_index] ;
    my $column = $self->{_columns}->[$sheet_index] ;
    my $format = $self->{_std_format};
    my $delimiter = $self->{_delimiter};

    if (defined $args{row}){
        $row = $args{row};
    }
        
    if (defined $args{column}){
        $column = $args{column};
    }
        
    if (defined $args{format}){
        $format = $args{format};
    }
        
    if (defined $args{delimiter}){
        $delimiter = $args{delimiter};
    }

    if (ref $args{line} eq 'ARRAY'){
        @fields = @{$args{line}};
    }else{
        @fields = split(/$delimiter/, $args{line});
    }
    if (exists $args{preceding}){
        if (not ref $args{preceding} eq 'ARRAY'){
            croak "value for 'preceding' passed to _writeLineToExcel must be an ARRAY reference " ;
        }
        my $added = @{$args{preceding}} - 1; 
        $additional_rows = $added > $additional_rows ? $added : $additional_rows;
    }
    if (exists $args{succeeding}){
        if (not ref $args{succeeding} eq 'ARRAY'){
            croak "value for 'succeeding' passed to _writeLineToExcel must be an ARRAY reference " ;
        }
        my $added = @{$args{succeeding}} - 1; 
        $additional_rows = $added > $additional_rows ? $added : $additional_rows;
    }
    if (exists $args{preceding}){
        my $added = @{$args{preceding}} - 1; 

        my $rows_per_value;
        my $row_extra;
        if ($added < $additional_rows and @{$args{preceding}} > 0){
            $rows_per_value = int (($additional_rows + 1)/@{$args{preceding}});
            $row_extra = (($additional_rows + 1) % @{$args{preceding}});
        }

        my $temp_row = $row;
        for (my $i = 0; $i < @{$args{preceding}}; $i++){
            my $temp_column = $column;
            my @p_fields = ();
            if (ref $args{preceding}->[$i] eq 'ARRAY'){
                @p_fields = @{$args{preceding}->[$i]};
            }else{
                @p_fields = split(/$delimiter/, $args{preceding}->[$i]);
            }
            foreach my $p (@p_fields){
                if (not defined $p){
                    $p = '';
                }

                if ($rows_per_value){
                    my $type = 'string';
                    if ($p =~ /^\d+(\.\d+)*]$/){
                        $type = 'number';
                    }
                    my $top_cell = xl_rowcol_to_cell($temp_row, $temp_column);
                    my $bottom_cell = xl_rowcol_to_cell($temp_row + $rows_per_value - 1, $temp_column );
                    if ($i == $#{$args{preceding}}){
                        $bottom_cell = xl_rowcol_to_cell($temp_row + $rows_per_value + $row_extra - 1, $temp_column );
                    }
                    if ($bottom_cell eq $top_cell){
                        $worksheet->write($temp_row, $temp_column, $p, $format);
                    }else{
                        $worksheet->merge_range_type($type, "$top_cell:$bottom_cell", $p, $format);
                    }
                    $temp_column++;
                }else{
                    $worksheet->write($temp_row, $temp_column++, $p, $format);
                }
            }
            if ($rows_per_value){
                $temp_row += $rows_per_value;
                if ($i == $#{$args{preceding}}){
                    $temp_row += $row_extra;
                }
            }else{
                $temp_row++;
            }
            $columns_before = @p_fields > $columns_before ? @p_fields : $columns_before;
        }
    }
    if (exists $args{succeeding}){
        my $added = @{$args{succeeding}} - 1; 
        my $rows_per_value;
        my $row_extra;
        if ($added < $additional_rows and @{$args{succeeding}} > 0){
            $rows_per_value = int (($additional_rows + 1)/@{$args{succeeding}});
            $row_extra = (($additional_rows + 1) % @{$args{succeeding}});
        }

        my $temp_row = $row; 

        for (my $i = 0; $i < @{$args{succeeding}}; $i++){
            my $temp_column = $column + @fields + $columns_before;
            my @p_fields = ();
            if (ref $args{succeeding}->[$i] eq 'ARRAY'){
                @p_fields = @{$args{succeeding}->[$i]};
            }else{
                @p_fields = split(/$delimiter/, $args{succeeding}->[$i]);
            }
            
            foreach my $p  (@p_fields){
                if (not defined $p){
                    $p = ' ';
                }
                if ($rows_per_value){
                    my $type = 'string';
                    if ($p =~ /^\d+(\.\d+)*]$/){
                        $type = 'number';
                    }
                    my $top_cell = xl_rowcol_to_cell($temp_row, $temp_column);
                    my $bottom_cell = xl_rowcol_to_cell($temp_row + $rows_per_value - 1, $temp_column );
                    if ($i == $#{$args{succeeding}}){
                        $bottom_cell = xl_rowcol_to_cell($temp_row + $rows_per_value + $row_extra - 1, $temp_column );
                    }
                    if ($bottom_cell eq $top_cell){
                        $worksheet->write($temp_row, $temp_column, $p, $format);
                    }else{
                        $worksheet->merge_range_type($type, "$top_cell:$bottom_cell", $p, $format);
                    }
                    $temp_column++;
                }else{
                    $worksheet->write($temp_row, $temp_column++, $p, $format);
                }
            }
            if ($rows_per_value){
                $temp_row += $rows_per_value;
                if ($i == $#{$args{succeeding}}){
                    $temp_row += $row_extra;
                }
            }else{
                $temp_row++;
            }
        }
    }

    my $temp_column = $column + $columns_before;
    foreach my $f (@fields){
        if (not defined $f){
            $f = '';
        } 
        if ($additional_rows > 0){
            my $type = 'string';
            if ($f =~ /^\d+(\.\d+)*]$/){
                $type = 'number';
            }
            my $top_cell = xl_rowcol_to_cell($row, $temp_column);
            my $bottom_cell = xl_rowcol_to_cell($row + $additional_rows, $temp_column );
            $worksheet->merge_range_type($type, "$top_cell:$bottom_cell", $f, $format);
            $temp_column++;
        }else{
            $worksheet->write($row, $temp_column++, $f, $format);
        }
    }
    $self->{_rows}->[$sheet_index] = $row + $additional_rows + 1;
    $self->{_columns}->[$sheet_index]  = 0;
    return $row + $additional_rows + 1;
}

1;
