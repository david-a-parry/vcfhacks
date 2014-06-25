package TextToExcel;
use strict;
use warnings;
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
        _row => ["0", "read"],
        _column => ["0", "read"],
        _worksheet => ["", "read/write"],
        _workbook => ["", "read"],
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
        $self->{_workbook}->close;
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
    $self->{_worksheet} = $self->{_workbook}->add_worksheet();
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

sub createFormat{
    #creates a format for use with workbook
    #returns a format that can be passed to the writeLine 
    #function using 
    my ($self, %args) = @_;
    my $format = $self->{_workbook}->add_format(%args);
    return $format;
}

sub writeLine{
    my ($self, %args) = @_;
    $self->_checkWriteArguments("writeLine", %args);
    $self->{_row} = $self->_writeLineToExcel(%args);
}

sub _writeLineToExcel{
    #returns the number of the row in the worksheet after writing this line
    #usually this will just be used by writeLine function 
    my ($self, %args) = @_;
    $self->_checkWriteArguments("_writeLineToExcel", %args);

    my $columns_before = 0;
    my $additional_rows = 0;
    my @fields = ();
    my $worksheet = $self->{_worksheet};
    my $row = $self->{_row};
    my $column = $self->{_column};
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
        $delimiter = $args{worksheet};
    }
        
    if (defined $args{worksheet}){
        $worksheet = $args{worksheet};
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
        my $temp_row = $row;
        foreach my $p_row (@{$args{preceding}}){
            my $temp_column = $column;
            my @p_fields = ();
            if (ref $p_row eq 'ARRAY'){
                @p_fields = @$p_row;
            }else{
                @p_fields = split(/$delimiter/, $p_row);
            }
            foreach my $p (@p_fields){
                if (not defined $p){
                    $p = '';
                }
                $worksheet->write($temp_row, $temp_column++, $p, $format);
            }
            $columns_before = @p_fields > $columns_before ? @p_fields : $columns_before;
            $temp_row++;
            $added++;
        }
    }
    if (exists $args{succeeding}){
        if (not ref $args{succeeding} eq 'ARRAY'){
            croak "value for 'succeeding' passed to _writeLineToExcel must be an ARRAY reference " ;
        }
        my $added = @{$args{succeeding}} - 1; 
        $additional_rows = $added > $additional_rows ? $added : $additional_rows;
        my $temp_row = $row; 
        foreach my $p_row (@{$args{succeeding}}){
            my $temp_column = $column + @fields;
            my @p_fields = ();
            if (ref $p_row eq 'ARRAY'){
                @p_fields = @$p_row;
            }else{
                @p_fields = split(/$delimiter/, $p_row);
            }
            foreach my $p (@p_fields){
                $worksheet->write($temp_row, $temp_column++, $p, $format);
            }
            $columns_before = @p_fields > $columns_before ? @p_fields : $columns_before;
            $temp_row++;
        }
        $additional_rows = ($temp_row - $row) > $additional_rows ? $temp_row - $row : $additional_rows;
    }
    my $temp_column = $column + $columns_before;
    foreach my $f (@fields){
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
    return $row + $additional_rows + 1;
}

1;
