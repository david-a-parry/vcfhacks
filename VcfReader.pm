package VcfReader;
use strict;
use warnings;
use Carp;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use IO::Compress::Gzip qw(gzip $GzipError) ;
use Fcntl 'SEEK_SET';
use Data::Dumper;
require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();

my %vcf_fields = (
    CHROM   => 0, 
    POS     => 1,
    ID      => 2,
    REF     => 3,
    ALT     => 4, 
    QUAL    => 5, 
    FILTER  => 6, 
    INFO    => 7, 
    FORMAT  => 8,
    );

#index every 0-9999 bp of a chrom
my $REGION_SPANS = 10000;

#header utilities
sub getHeader{
    my $vcf = shift; 
    croak "getHeader method requires a file as an argument" if not $vcf;
    my $FH = _openFileHandle($vcf);
    my @header = ();
    while (my $vcf_line = scalar readline $FH){
        if ($vcf_line =~ /^#/){
            chomp $vcf_line;
            push(@header, $vcf_line);
        }else{
            croak "No header found for VCF file $vcf " if not @header;
            last;
        }
    }
    close $FH;
    croak "No header found for VCF file $vcf " if not @header;
    return @header if wantarray;
    return join("\n", @header) if defined wantarray;
    carp "getHeader called in void context ";
}

sub getMetaHeader{
    my $vcf = shift; 
    croak "printMetaHeader method requires a file as an argument" if not $vcf;
    my @header = grep {/^##/} getHeader($vcf);
    return @header if wantarray;
    return join("\n", @header) if defined wantarray;
    carp "getMetaHeader called in void context ";
}

sub getColumnHeader{
    my $vcf = shift; 
    croak "printColumnHeader method requires a file as an argument" if not $vcf;
    my @header = grep {/^#CHROM/} getHeader($vcf);
    if (@header < 1){
        croak "No column header found for $vcf ";
    }
    if (@header > 1){
        carp "Warning - more than 1 column header found for $vcf ";
    }
    return "$header[-1]";
}

sub printHeader{
    my $vcf = shift; 
    croak "printHeader method requires a file as an argument" if not $vcf;
    print join("\n", getHeader($vcf));
}

sub checkHeader{
    #returns 1 if header ok
    #returns 0 if not
    my (%args) = @_;
    my @header = ();
    if ($args{vcf}){
        @header = getHeader($args{vcf});
    }elsif($args{header}){
        if (ref $args{header} eq 'ARRAY'){
            @header = @{$args{header}};
        }else{
            @header = split("\n", $args{header});
        }
    }else{
        croak "readHeader method requires either 'vcf' or 'header arguments ";
    }
    if ($header[-1] !~ /^#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO(\tFORMAT)*/){
        return 0;
    }else{
        return 1;
    }
}

sub getHeaderColumns{
    my (%args) = @_;
    croak "Invalid header " if not checkHeader(%args);
    my @header = ();
    if ($args{vcf}){
        @header = getHeader($args{vcf});
    }elsif($args{header}){
        if (ref $args{header} eq 'ARRAY'){
            @header = @{$args{header}};
        }else{
            @header = split("\n", $args{header});
        }
    }else{
        croak "getHeaderColumns method requires either 'vcf' or 'header arguments ";
    }
    my @columns = split("\t", $header[-1]);
    return @columns;
}

sub getSamples{
#return either an array of sample names
#of a hash of sample_name => column number (if called with "get_columns => 1")
    my (%args) = @_;
    croak "Invalid header " if not checkHeader(%args);
    my @header = ();
    if ($args{vcf}){
        @header = getHeader($args{vcf});
    }elsif($args{header}){
        if (ref $args{header} eq 'ARRAY'){
            @header = @{$args{header}};
        }else{
            @header = split("\n", $args{header});
        }
    }else{
        croak "readHeader method requires either 'vcf' or 'header arguments ";
    }
    my @columns = split("\t", $header[-1]);
    return if @columns < 10;
    my @samples = @columns[9..$#columns];
    if (exists $args{get_columns} && $args{get_columns}){
        my $n = 9;
        my %samp = map {$_ => $n++} @samples;
        return %samp;
    }else{
        return @samples;
    }
}
            
sub getInfoFields{
    #return a hash of INFO IDs, to anon hashes of 
    #Number, Type and Description
    my (%args) = @_;
    my @header = ();
    if ($args{vcf}){
        @header = getHeader($args{vcf});
    }elsif($args{header}){
        if (ref $args{header} eq 'ARRAY'){
            @header = @{$args{header}};
        }else{
            @header = split("\n", $args{header});
        }
    }else{
        croak "getInfoFields method requires either 'vcf' or 'header arguments ";
    }
    my %info = ();
    foreach my $h (@header){
        if ($h =~ /^##INFO=<ID=(\w+),Number=([\.\w+]),Type=(\w+),Description=(.+)>$/){
            $info{$1} = 
                {
                Number => $2,
                Type => $3,
                Description => $4,
                };
        }
    }
    return %info;
}

sub getContigOrder{
    my ($vcf) = shift;
    croak "getContigOrder method requires a file as an argument" if not $vcf;
    my %contigs = ();
    my @meta = getMetaHeader($vcf);
    my @con = grep {/##contig=</} @meta;
    if (@con){
        my $n = 0;
        foreach my $c (@con){
            if ($c =~/ID=([^,>]+)[>,]/){
                $contigs{$1} = $n++;
            }else{
                carp "ERROR - failed to parse header contig line: $c ";
            }
        }
        return %contigs if %contigs;
    }
    print STDERR "Failed to retrive contigs from header - reading/creating index.\n";
    if ($vcf =~ /\.gz$/){
        eval "use Tabix; 1" 
            or croak "Tabix module is not installed and VCF file $vcf appears to be (b)gzip compressed.  ".
            "  Please install Tabix.pm in order to quickly extract contigs from bgzip compressed VCFs.\n";
        my $index = "$vcf.tbi";
        if (not -e $index){
            print STDERR "Indexing $vcf with tabix...\n";
            indexVcf($vcf);
            croak "Tabix indexing failed? $index does not exist " if (not -e $index);
        }
        my $t = getTabixIterator($vcf, $index);
        my $n = 0;
        %contigs = map {$_ => $n++} $t->getnames();
    }else{
        my %idx = readIndex($vcf);
        foreach my $k (keys %idx){
            if (ref $idx{$k} eq 'HASH' && exists $idx{$k}->{order}){
                $contigs{$k} = $idx{$k}->{order};
            }
        }
    }
    return %contigs;
}

sub getLineCount{
    my $vcf = shift; 
    croak "getLineCount method requires a file as an argument" if not $vcf;
    my $line_count = 0;
    my $FH = _openFileHandle($vcf);
    if ($vcf =~ /\.gz$/){
       # $line_count++ while (<$FH>);
        $line_count += tr/\n/\n/ while sysread($FH, $_, 2 ** 20);
    }else{
        my $index = "$vcf.vridx";
        if (-e $index){
            $line_count = getFileLengthFromIndex($vcf, $index); 
        }else{
            $line_count += tr/\n/\n/ while sysread($FH, $_, 2 ** 20);
        }
    }
    close $FH;
    return $line_count;
}

sub checkCoordinateSorted{
#return 1 if vcf is coordinate sorted (not checking chrom order)
    my ($vcf) = @_;
    my $FH = _openFileHandle($vcf);
    my %contigs = ();
    my $prev_pos = 0;
    while (my $line = <$FH>){
        next if $line =~ /^#/;
        my @split = split("\t", $line);
        if (exists $contigs{$split[0]}){
            if ($contigs{$split[0]}  != scalar(keys%contigs) -1 ){
                return 0; #not sorted - encountered contig twice with another inbetween
            }
            return 0 if $split[1] < $prev_pos;
        }else{
            $contigs{$split[0]}  = scalar(keys%contigs);
        }
        $prev_pos = $split[1];
    }
    return 1;
}

sub indexVcf{
    local $Data::Dumper::Terse = 1;
    local $Data::Dumper::Useqq = 1;
    my ($vcf) = @_;
    if ($vcf =~ /\.gz$/){
        chomp (my $tabix = `which tabix`);
        croak "Can't find tabix executable to index compressed VCF.  Please ensure it is ".
            "installed and in your PATH or index your VCF with tabix manually. "
            if not $tabix; 
        my $er = `$tabix -p vcf $vcf 2>&1`;
        croak "Tabix indexing failed, code $?: $er\n" if $?;
        return;
    }
    my $index = "$vcf.vridx";
    #open (my $INDEX, "+>$index") or croak "can't open $index for writing index: $! ";
    my $gz = new IO::Compress::Gzip $index
        or croak "IO::Compress::Gzip failed to write index: $GzipError\n";
    my $offset = 0;
    my $prev_offset = 0;
    my $prev_pos = 0;
    my %contigs = ();
    my $is_var = 0;
    my $n = 0;
    my $last_line_indexed = 0;
    
    my $FH = _openFileHandle($vcf);
    while (my $line = scalar readline $FH){
        $n++;
        $prev_offset = $offset;
        #print $INDEX pack($pack, $offset);
        $offset = tell $FH;
        next if $line =~ /^#/;
        $is_var++;
        if ($is_var == 1){
            $contigs{first_line} = $n;
            $contigs{first_offset} = $prev_offset;
        }
        my @s = split ("\t", $line);
        my $chrom = $s[$vcf_fields{CHROM}];
        my $pos = $s[$vcf_fields{POS}];
        my $span = $s[$vcf_fields{POS}] + length($s[$vcf_fields{REF}]) - 1;
        my $pos_rounddown = int($pos/$REGION_SPANS) * $REGION_SPANS;
        my $span_rounddown = int($span/$REGION_SPANS) * $REGION_SPANS;
        if (exists $contigs{$chrom}){
            if ($contigs{$chrom}->{order}  != scalar(keys%contigs) -1 or $pos < $prev_pos){
                $gz->close();
                unlink $index or carp "Couldn't delete partial index $index - please delete manually: $!" ;
                croak "Can't index VCF - $vcf not sorted properly ";
            }
            if (exists $contigs{$chrom}->{regions}->{$pos_rounddown}){
                my $merged = 0;
                foreach my $offs (@{$contigs{$chrom}->{regions}->{$pos_rounddown}}){
                    if ($offs->{line_end} +1  eq $n){
                        #if lines are contiguous create a single region
                        $offs->{offset_end} = $offset;
                        $offs->{line_end} = $n;
                        $offs->{pos_end} = $span if ($span > $offs->{pos_end});
                        $merged++;
                        last;
                    }
                }
                if (not $merged){
                    push @{$contigs{$chrom}->{regions}->{$pos_rounddown}},
                        {offset_start=> $prev_offset,  offset_end => $offset,
                        line_start => $n, line_end => $n, pos_start => $pos, pos_end => $span};
                }
            }else{
                push @{$contigs{$chrom}->{regions}->{$pos_rounddown}},
                    {offset_start=> $prev_offset,  offset_end => $offset,
                    line_start => $n, line_end => $n, pos_start => $pos, pos_end => $span};
            }
            if (exists $contigs{$chrom}->{regions}->{$span_rounddown}){
                my $merged = 0;
                foreach my $offs (@{$contigs{$chrom}->{regions}->{$span_rounddown}}){
                    if ($offs->{line_end} eq $n){
                        #already merged
                        $merged++;
                        last;
                    }elsif ($offs->{line_end} +1  eq $n){
                        #if lines are contiguous create a single region
                        $offs->{offset_end} = $offset;
                        $offs->{line_end} = $n;
                        $offs->{pos_end} = $span if ($span > $offs->{pos_end});
                        $merged++;
                        last;
                    }
                }
                if (not $merged){
                    push @{$contigs{$chrom}->{regions}->{$span_rounddown}},
                        {offset_start=> $prev_offset,  offset_end => $offset,
                        line_start => $n, line_end => $n, pos_start => $pos, pos_end => $span};
                }
            }else{
                push @{$contigs{$chrom}->{regions}->{$span_rounddown}},
                    {offset_start=> $prev_offset,  offset_end => $offset,
                    line_start => $n, line_end => $n, pos_start => $pos, pos_end => $span};
            }
        }else{
            $contigs{$chrom}->{order}  = scalar(keys%contigs);
            push @{$contigs{$chrom}->{regions}->{$pos_rounddown}},
                {offset_start=> $prev_offset,  offset_end => $offset,
                line_start => $n, line_end => $n, pos_start => $pos, pos_end => $span};
        }
        $prev_pos = $pos;
    }
    $contigs{last_line} = $n;
    $contigs{last_offset} = $prev_offset;
    close $FH;
    print $gz Dumper \%contigs;
    close $gz;
}

sub _openFileHandle{
    my $vcf = shift;
    croak "_openFileHandle method requires a file as an argument" if not $vcf;
    my $FH; 
    if ($vcf =~ /\.gz$/){
        $FH = new IO::Uncompress::Gunzip $vcf, MultiStream => 1 or croak "IO::Uncompress::Gunzip failed while opening $vcf for reading: \n$GunzipError";
    }else{
        open ($FH, $vcf) or croak "Failed to open $vcf for reading: $! ";
    }
    return $FH;
}


#variant utilitiea
sub getVariantField{
    my ($line, $field) = @_;
    $field =~ s/^#+//;
    croak "Invalid field ($field) passed to getVariantField method " 
        if not exists $vcf_fields{$field};
    my @split = ();
    if (ref $line eq 'ARRAY'){
        @split = @$line;
    }else{
        @split = split("\t", $line);
    }
    if ($vcf_fields{$field} > $#split){
        if ($field eq 'FORMAT'){
            #carp "No FORMAT field for line " . join("\t", @split) . " " ;
            return;
        }else{
            croak "Line has too few fields: " . join("\t", @split) . " " ;
        }
    }
    return $split[$vcf_fields{$field}];
}

sub getVariantInfoField{
    my ($line, $info_field) = @_;
    my @split = ();
    if (ref $line eq 'ARRAY'){
        @split = @$line;
    }else{
        @split = split("\t", $line);
    }
    my @info = split(';', getVariantField(\@split, "INFO"));
    foreach my $inf (@info){
        if ($inf =~ /^$info_field=(.+)/){
            return $1;
        }elsif ($inf eq $info_field){
            return 1;
        }
    }
    return;
}

sub readAlleles{
    my (%args) = @_;
    croak "line argument is required for readAlleles method " if not defined $args{line};
    my @alleles = split(",", getVariantField($args{line}, "ALT"));
    #if no variant at position then ALT will be "." and we want to remove it from our alleles array
    @alleles = grep {! /\./ } @alleles; 
    if (defined $args{alt_alleles} && $args{alt_alleles}){
        return @alleles if defined wantarray;
        carp "readAlleles method called in void context ";
    }
    unshift(@alleles, getVariantField($args{line}, "REF"));
    #now ref is at index 0 of @alleles, and the alts correspond to the number in the call field
    return @alleles if defined wantarray;
    carp "readAlleles method called in void context ";
}

sub isMultiAllelic{
    my ($line) = @_;
    my $alts = readAlleles(line => $line, alt_alleles => 1);
    return 1 if $alts > 1;
    return 0;
}

sub getFormatFields{
    my ($line) = @_;
    my $format = getVariantField($line, "FORMAT");
    if (not defined $format){
        return;
    }
    my @form = split(":", $format);
    my %format_fields = ();
    my $i = 0;
    foreach my $f (@form){
        $format_fields{$f} = $i++;
    }
    return %format_fields if wantarray;
    return keys %format_fields if defined wantarray;
    carp "getFormatFields method called in void context ";
}

sub getSampleVariant{
    my ($line, $column) = @_;
    croak "Invalid column ($column) for getSampleVariant - " .
        "samples are only present in columns 10 and onwards.\n"
        if $column < 9;
    my @split = ();
    if (ref $line eq 'ARRAY'){
        @split = @$line;
    }else{
        @split = split("\t", $line);
    }
    if ($column){
        return $split[$column];
    }else{#otherwise just return first sample
        croak "Line has too few fields: " . join("\t", @split) . " " if 9 > $#split;
        return $split[9];
    }
}

sub getSampleGenotypeField{
#returns scalar value for genotype field 
#unless multiple argument is used in which
#case a hash of sample=>value key/values is returned
#assumes identifier is column unless a hash of sample names to columns
#is provided as "sample_to_columns" argument, in which case identifier is assumed to be sample ID
#column can be used as an argument instead of sample to explicitly specify it as a column rather than id
    my (%args) = @_;
    croak "\"field\" argument must be passed to getSampleGenotypeField - e.g. getSampleGenotypeField(field=>\"GQ\") " if not defined $args{field};
    croak "\"line\" argument must be passed to getSampleGenotypeField " if not defined $args{line};
    carp "WARNING Both multiple and sample arguments supplied to getSampleGenotypeField method - only multiple argument will be used " if (defined $args{multiple} and defined $args{sample});
    if (defined $args{sample_to_columns}){
        if (ref $args{sample_to_columns} ne 'HASH'){
            croak "\"sample_to_columns\" argument passed to getSampleGenotypeField must be a hash reference ";
        }
    #}elsif(defined $args{sample} or $args{multiple}){
#        croak "\"multiple\" and \"sample\" arguments can only be used in conjunction with \"sample_to_columns\" option for getSampleGenotypeField method ";
    }
    my %var_format = getFormatFields($args{line});
    if (not defined $var_format{$args{field}}){
        carp "Field $args{field} not found for getSampleGenotypeField ";
        return;
    }
    my @split = ();
    if (ref $args{line} eq 'ARRAY'){
        @split = @{$args{line}};
    }else{
        @split = split("\t", $args{line});
    }
    #croak "Line has too few fields: " . join("\t", @split) . " " if $self->{_samples}->{$sample} > $#split;
    my $var; 
    if ($args{all}){
        my %values = ();
        if (defined $args{sample_to_columns}){
            foreach my $sample (keys %{$args{sample_to_columns}}){
                my $mvar = getSampleVariant(
                                    $args{line}, 
                                    $args{sample_to_columns}->{$sample},
                                    );
                my $value = (split ":", $mvar)[$var_format{$args{field}}];
                $values{$sample} =  $value;
            }
        }else{
            foreach my $col (9..$#split){
                my $mvar = getSampleVariant(
                                    $args{line}, 
                                    $col,
                                    );
                my $value = (split ":", $mvar)[$var_format{$args{field}}];
                $values{$col} =  $value;

            }
        }
        return %values if defined wantarray;
        carp "getSampleGenotypeField called in a void context ";
    }elsif($args{multiple}){
        croak "multiple argument must be an array reference " if ref $args{multiple} ne 'ARRAY';
        my %values = ();
        foreach my $sample (@{$args{multiple}}){
            if (defined $args{sample_to_columns}){
                my $mvar = getSampleVariant(
                                    $args{line}, 
                                    $args{sample_to_columns}->{$sample},
                                    );
                my $value = (split ":", $mvar)[$var_format{$args{field}}];
                $values{$sample} =  $value;
            }else{
                my $mvar = getSampleVariant(
                                    $args{line}, 
                                    $sample,
                                    );
                my $value = (split ":", $mvar)[$var_format{$args{field}}];
                $values{$sample} =  $value;
            }
        }
        return %values if defined wantarray;
        carp "getSampleGenotypeField called in a void context ";
    }elsif (defined $args{sample}){
        my $col = $args{sample};
        if (defined $args{sample_to_columns}){
            $col = $args{sample_to_columns}->{$args{sample}};
        }
        $var = getSampleVariant(
                            $args{line}, 
                            $args{sample_to_columns}->{$args{sample}},
                            );
    }elsif (defined $args{column}){
        $var = getSampleVariant(
                                $args{line}, 
                                $args{column},
                                );
    }else{#otherwise just look at first sample
        $var = getSampleVariant($args{line});
    }
    return if not defined $var and defined wantarray;
    my $value = (split ":", $var)[$var_format{$args{field}}];
    return $value if defined wantarray;
    carp "getSampleGenotypeField called in a void context ";
}


sub getSampleCall{
#returns scalar value for genotype called 
#unless multiple argument is used in which
#case a hash of sample=>call key/values is returned
#returns './.' for samples below $args{minGQ}
#use return_alleles_only => 1 to only return allele codes, not genotypes (always returns an array)
    my (%args) = @_;
    croak "\"field\" argument must be passed to getSampleCall - e.g. getSampleCall(field=>\"GQ\") " if not defined $args{field};
    croak "\"line\" argument must be passed to getSampleCall " if not defined $args{line};
    carp "WARNING Both multiple and sample arguments supplied to getSampleCall method - only multiple argument will be used " if (defined $args{multiple} and defined $args{sample});
    if (defined $args{sample_to_columns}){
        if (ref $args{sample_to_columns} ne 'HASH'){
            croak "\"sample_to_columns\" argument passed to getSampleCall must be a hash reference ";
        }
#    }elsif(defined $args{sample} or $args{multiple}){
#        croak "\"multiple\" and \"sample\" arguments can only be used in conjunction with \"sample_to_columns\" option for getSampleCall method ";
    }
    my @split = ();
    if (ref $args{line} eq 'ARRAY'){
        @split = @{$args{line}};
    }else{
        @split = split("\t", $args{line});
    }
    my $var; 
    my %calls = ();
    if ($args{all}){
        if (defined $args{sample_to_columns}){
            foreach my $sample (keys %{$args{sample_to_columns}}){
                if (exists $args{minGQ} and $args{minGQ} > 0){
                    $calls{$sample} = _getGenotype(
                                            $args{line},
                                            $args{sample_to_columns}->{$sample}, 
                                            $args{minGQ}
                                            );
                }else{
                    $calls{$sample} = _getGenotype(
                                                    $args{line},
                                                    $args{sample_to_columns}->{$sample}, 
                                                    );
                }
            }
        }else{
            foreach my $col (9..$#split){
                if (exists $args{minGQ} and $args{minGQ} > 0){
                    $calls{$col} = _getGenotype(
                                            $args{line},
                                            $col, 
                                            $args{minGQ}
                                            );
                }else{
                    $calls{$col} = _getGenotype(
                                                $args{line},
                                                $col, 
                                            );
                }
            }
        }
        if ($args{return_alleles_only}){
            my %rev_calls = reverse %calls;
            my %allele_codes = ();
            foreach my $g ( keys %rev_calls){
                my @al = split(/[\/\|]/, $g);
                foreach my $a (@al){
                    $allele_codes{$a}++;
                }
            }
                return keys %allele_codes if defined wantarray;
        }else{
            return %calls if defined wantarray;
        }
        carp "getSampleCall called in a void context ";
    }elsif($args{multiple}){
        croak "multiple argument must be an array reference " if ref $args{multiple} ne 'ARRAY';
        croak "Sample \"$sample\" does not exist in samples_to_columns hash passed to getSampleCall " 
            if not exists $args{sample_to_columns}->{$sample};
        if (defined $args{sample_to_columns}){
            foreach my $sample (@{$args{multiple}}){
                if (exists $args{minGQ} and $args{minGQ} > 0){
                    $calls{$sample} = _getGenotype(
                                                $args{line},
                                                $args{sample_to_columns}->{$sample}, 
                                                $args{minGQ}
                                                );
                }else{
                    $calls{$sample} = _getGenotype(
                                                    $args{line},
                                                    $args{sample_to_columns}->{$sample}, 
                                                    );
                }
            }
        }else{
            foreach my $col (@{$args{multiple}}){
                if ($col !~ /^\d+$/){
                    croak "columns passed to getSampleCall must be integers not \"$col\". ".
                    "To parse sample names you must pass a hash reference of sample IDs ".
                    "to columns as generated by the getSamples sub as the sample_to_columns argument ";
                    
                if (exists $args{minGQ} and $args{minGQ} > 0){
                    $calls{$col} = _getGenotype(
                                                $args{line},
                                                $col, 
                                                $args{minGQ}
                                                );
                }else{
                    $calls{$col} = _getGenotype(
                                                    $args{line},
                                                    $col, 
                                                    );
                }
            }
        }
        if ($args{return_alleles_only}){
            my %rev_calls = reverse %calls;
            my %allele_codes = ();
            foreach my $g ( keys %rev_calls){
                my @al = split(/[\/\|]/, $g);
                foreach my $a (@al){
                    $allele_codes{$a}++;
                }
            }
                return keys %allele_codes if defined wantarray;
        }else{
            return %calls if defined wantarray;
        }
        carp "getSampleCall called in a void context ";
    }else{
        my $call;
        my $col = 9;#default is to get first sample
        if (defined $args{sample}){
            if (defined $args{sample_to_columns}){
                if (not exists $args{sample_to_columns}->{$sample}){
                    croak "Sample \"$sample\" does not exist in samples_to_columns hash passed to getSampleCall " 
                }
                $col = $args{sample_to_columns}->{$args{sample}};
            }else{
                $col = $args{sample};
            }
        }elsif (defined $args{column}){
            $col = $args{column};
        }
        if (exists $args{minGQ} and $args{minGQ} > 0){
            $call = _getGenotype(
                        $args{line},
                        $col,
                        $args{minGQ}
                        );
        }else{
            $call = _getGenotype(
                        $args{line},
                        $col,
                        );
        }
        if ($args{return_alleles_only}){
            if ($call eq './.'){
                return '.' if defined wantarray;
            }
            my @al = split(/[\/\|]/, $call);
            return @al if defined wantarray;
        }else{
            return $call if defined wantarray;
        }
    }
    carp "getSampleCall called in a void context ";
}


sub getSampleActualGenotypes{
    my ($self, %args) = @_;
    croak "\"line\" argument must be passed to getSampleActualGenotypes " if not defined $args{line};
    carp "WARNING Both multiple and sample arguments supplied to getSampleActualGenotypes method - only multiple argument will be used " if (defined $args{multiple} and defined $args{sample});
    my @split = ();
    if (ref $args{line} eq 'ARRAY'){
        @split = @{$args{line}};
    }else{
        @split = split("\t", $args{line});
    }
    if (defined $args{sample_to_columns}){
        if (ref $args{sample_to_columns} ne 'HASH'){
            croak "\"sample_to_columns\" argument passed to getSampleActualGenotypes must be a hash reference ";
        }
    #}elsif(defined $args{sample} or $args{multiple}){
    #    croak "\"multiple\" and \"sample\" arguments can only be used in conjunction with \"sample_to_columns\" option for getSampleCall method ";
    }
    my @alleles = $self -> readAlleles(line => \@split);
    my %multiple = ();
    my $genotype;
    my @sample_alleles = ();
    my $var;
    if ($args{all}){
        if (defined $args{sample_to_columns}){
            foreach my $sample (keys %{$args{sample_to_columns}}){
                $var = getSampleVariant(\@split, $args{sample_to_columns}->{$sample});
                my $call = (split ":", $var)[0];
                if ($call =~ /(\d+)[\/\|](\d+)/){
                    if ($args{return_alleles_only}){
                        push (@sample_alleles, ($alleles[$1], $alleles[$2]));
                    }else{
                        $multiple{$sample} = "$alleles[$1]/$alleles[$2]";
                    }
                }else{
                    if (not $args{return_alleles_only}){
                        $multiple{$sample} = "-/-";
                    }
                }
            }
        }else{
            foreach my $col (9..$#split){
                $var = getSampleVariant(\@split, $col);
                my $call = (split ":", $var)[0];
                if ($call =~ /(\d+)[\/\|](\d+)/){
                    if ($args{return_alleles_only}){
                        push (@sample_alleles, ($alleles[$1], $alleles[$2]));
                    }else{
                        $multiple{$col} = "$alleles[$1]/$alleles[$2]";
                    }
                }else{
                    if (not $args{return_alleles_only}){
                        $multiple{$col} = "-/-";
                    }
                }
            }
        }
        if ($args{return_alleles_only}){
            my %seen = ();
            @sample_alleles = grep {!$seen{$_}++} @sample_alleles;#remove duplicates
            return @sample_alleles;
        }else{
            return %multiple;
        }
    }elsif($args{multiple}){
        croak "multiple argument must be an array reference " if ref $args{multiple} ne 'ARRAY';
        if (defined $args{sample_to_columns}){
            foreach my $sample (@{$args{multiple}}){
                $var = getSampleVariant(\@split, $args{sample_to_columns}->{$sample});
                my $call = (split ":", $var)[0];
                if ($call =~ /(\d+)[\/\|](\d+)/){
                    if ($args{return_alleles_only}){
                        push (@sample_alleles, ($alleles[$1], $alleles[$2]));
                    }else{
                        $multiple{$sample} = "$alleles[$1]/$alleles[$2]";
                    }
                }else{
                    if (not $args{return_alleles_only}){
                        $multiple{$sample} = "-/-";
                    }
                }
            }
        }else{
            foreach my $col (@{$args{multiple}}){
                $var = getSampleVariant(\@split, $col);
                my $call = (split ":", $var)[0];
                if ($call =~ /(\d+)[\/\|](\d+)/){
                    if ($args{return_alleles_only}){
                        push (@sample_alleles, ($alleles[$1], $alleles[$2]));
                    }else{
                        $multiple{$col} = "$alleles[$1]/$alleles[$2]";
                    }
                }else{
                    if (not $args{return_alleles_only}){
                        $multiple{$col} = "-/-";
                    }
                }
            }
        }
        if ($args{return_alleles_only}){
            my %seen = ();
            @sample_alleles = grep {!$seen{$_}++} @sample_alleles;#remove duplicates
            return @sample_alleles;
        }else{
            return %multiple;
        }
    }else{
        my $col = 9;
        if ($args{sample}){
            if (defined $args{sample_to_columns}){
                $col = $args{sample_to_columns}->{$args{sample}};
            }
        }elsif ($args{column}){
            $col = $args{column};
        }
        $var = getSampleVariant(\@split, $col);
        my $call = (split ":", $var)[0];
        if ($call =~ /(\d+)[\/\|](\d+)/){
            if ($args{return_alleles_only}){
                push (@sample_alleles, ($alleles[$1], $alleles[$2]));
                my %seen = ();
                @sample_alleles = grep {!$seen{$_}++} @sample_alleles;#remove duplicates
                return @sample_alleles;
                
            }else{
                $genotype = "$alleles[$1]/$alleles[$2]";
            }
        }else{
            $genotype = "-/-";
        }
        return if $args{return_alleles_only};
        return $genotype 
    }
}

sub _getGenotype{
    my ($line, $col, $min_gq) = @_;
    $col = 9 if not $col;
    $min_gq = 0 if not $min_gq;
    my $mvar = getSampleVariant( $line,
                                 $col,
                                 );
    if ($min_gq > 0){
        if (not defined getSampleGenotypeField(
                        line => $line, 
                        column=>$col, 
                        field=>'GQ')
        ){
            #no GQ field - return no call (the above generally means that the call is './.' anyway)
            return './.';
        }
        if (getSampleGenotypeField(
                        line => $line, 
                        column=>$col ,
                        field=>'GQ',
        ) < $min_gq){
            return './.';
        }
    }
    my $call = getSampleGenotypeField(
        line => $line, 
        column => $col, 
        field=>'GT'
    ); 
    return $call;
}




#returns genotype codes - e.g. 0/0, 0/1, 1/1
sub getAllPossibleGenotypeCodes{
    my ($line) = @_;
    my @alleles = readAlleles(line => $line);
    my @combinations = ();
    for (my $n = 0; $n < @alleles; $n++){
        for (my $m = 0; $m <= $n; $m++){
            push (@combinations, "$m/$n");
        }
    }
    return @combinations if defined wantarray;
    carp "getAllPossibleGenotypes called in void context ";
}

#returns actual allele genotypes - e.g. A/A, A/T, T/T
sub getAllPossibleGenotypes{
    my ($line) = @_;
    my @alleles = readAlleles(line => $line);
    my @combinations = ();
    for (my $n = 0; $n < @alleles; $n++){
        for (my $m = 0; $m <= $n; $m++){
            push (@combinations, "$alleles[$m]/$alleles[$n]");
        }
    }
    return @combinations if defined wantarray;
    carp "getAllPossibleGenotypes called in void context ";
}

sub replaceVariantField{
    my ($line, $field, $replacement) = @_;
    $field =~ s/^#+//;
    croak "Invalid field ($field) passed to getVariantField method " if not exists $vcf_fields{$field};
    my @split = ();
    if (ref $line eq 'ARRAY'){
        @split = @$line;
    }else{
        @split = split("\t", $line);
    }
    splice(@split, $vcf_fields{$field}, 1, $replacement);
    return join("\t", @split) if defined wantarray;
    carp "replaceVariantField called in void context ";
}

sub getTabixIterator{
    my ($vcf, $index) = @_;
    $index ||= "$vcf.tbi";
    return Tabix->new(-data =>  $vcf, -index => $index) ;
}

sub getSearchArguments{
#returns hash of arguments and values for passing to searchForPosition method
#return hash of vcf, file_handle, index_handle and contig_order 
#for uncompressed files
#or tabix_iterator for bgzip compressed files
    my ($vcf, $contig_index) = @_;
    if ($vcf =~ /\.gz/){ 
        eval "use Tabix; 1" 
            or croak "Tabix module is not installed and VCF file $vcf appears to be (b)gzip compressed.  ".
            "  Please install Tabix.pm in order to search bgzip compressed VCFs.\n";
        my $index = "$vcf.tbi";
        if (not -e $index){
            print STDERR "Indexing $vcf with tabix...";
            indexVcf($vcf);
            croak "Tabix indexing failed? $index does not exist " if (not -e $index);
            print STDERR " Done.\n";
        }
        return (tabix_iterator => getTabixIterator($vcf, $index));
    }else{
        my $FH = _openFileHandle($vcf);
        if ($contig_index){
            if (ref ($contig_index) ne 'HASH'){
                croak "second argument passed to getSearchArguments method must be a hash reference to a contig index ";
            }
        }else{
            my $index = "$vcf.vridx" ;
            if (not -e $index){
                print STDERR "$index does not exist - indexing $vcf...\n";
                indexVcf($vcf);
                croak "Indexing failed? $index does not exist " if (not -e $index);
                print STDERR " Done.\n";
            }
            my %contig_order  = readIndex($vcf);
            $contig_index = \%contig_order;
        }
        return (file_handle => $FH, contig_order => $contig_index);
    }
}

sub searchForPosition{
#if vcf argument is provided will use Tabix.pm (searchForPositionCompressed) or internal method (searchForPositionUncompressed)
#depending on file extension
#otherwise will use Tabix.pm if tabix_iterator argument is provided
#or internal method if file_handle argument is provided
#my @matching lines = $obj->searchForPosition(chrom => 1, pos = 2000000)
    my (%args) = @_;
    croak "chrom argument is required for searchForPosition method " if not exists $args{chrom};
    croak "pos argument is required for searchForPosition method " if not exists $args{pos};
    if (exists $args{vcf}){
        if ($args{vcf} =~ /\.gz$/){
            return searchForPositionCompressed(%args);
        }else{
            return searchForPositionUncompressed(%args);
        }
    }elsif(exists $args{tabix_iterator}){
        return searchForPositionCompressed(%args);
    }elsif(exists $args{file_handle}){
        croak "file_handle argument can only be used without vcf argument if contig_order is provided "
            if not $args{contig_order};
        return searchForPositionUncompressed(%args);
    }else{
        croak "vcf or tabix_iterator arguments are required for searchForPosition method " ;
    }
}
 
sub searchForPositionCompressed{
    my (%args) = @_;
    croak "vcf or tabix_iterator arguments are required for searchForPositionCompressed method " 
        if not exists $args{vcf} and not exists $args{tabix_iterator};
    eval "use Tabix; 1" 
        or croak "Tabix module is not installed and VCF file $args{vcf} appears to be (b)gzip compressed.  ".
        "  Please install Tabix.pm in order to search bgzip compressed VCFs.\n";
    
    my $tabixIterator; 
    if ($args{tabix_iterator}){
        $tabixIterator = $args{tabix_iterator};
    }else{
        my $index = defined $args{index} ? $args{index} : "$args{vcf}.tbi";
        if (not -e $index){
            print STDERR "Indexing $args{vcf} with tabix...";
            indexVcf($args{vcf});
            croak "Tabix indexing failed? $index does not exist " if (not -e $index);
            print STDERR " Done.\n";
        }
        $tabixIterator = Tabix->new(-data =>  $args{vcf}, -index => $index) ;
    }
    my $iter = $tabixIterator->query($args{chrom}, $args{pos} -1, $args{pos});
    return if not defined $iter->{_}; #$iter->{_} will be undef if our chromosome isn't in the vcf file
    my @matches = ();
    while (my $m =  $tabixIterator->read($iter)){
        push @matches, $m;
    } 
    return @matches if defined wantarray;
    carp "searchForPosition called in void context ";     
}           

sub searchForPositionUncompressed{
    my (%args) = @_;
    croak "vcf or file_handle arguments are required for searchForPositionCompressed method " 
        if not exists $args{vcf} and not exists $args{file_handle};
    my $contig_order;
    my $blocks;
    my $FH = exists $args{file_handle} ? $args{file_handle} : _openFileHandle($args{vcf});
    my $index;
    my $contig_index;
    if ($args{vcf}){
        $index = defined  $args{index} ?  $args{index} : "$args{vcf}.vridx" ;
    }
    if (exists $args{contig_order}){
        if (ref $args{contig_order} eq 'HASH'){
            $contig_order = $args{contig_order};
        }else{
            croak "contig_order argument passed to searchForPosition method must be a hash reference ";
        }
    }else{
        croak "contig_order argument is required to use searchForPositionUncompressed without vcf argument "
            if not exists $args{vcf};
    }
    if (not $contig_order){
        my %c  = readIndex($args{vcf});
        $contig_order = \%c;
        if (not %{$contig_order}){
            croak "Could not find any contigs in contig index $args{vcf}.vridx. Try deleting $args{vcf}.vridx and rerunning " ;
        }
    }

    my @matches = _searchVcf(
                            chrom           => $args{chrom}, 
                            pos             => $args{pos},
                            contig_order    => $contig_order,
                            fh              => $FH,
                            );
    return @matches if defined wantarray;
    carp "searchForPosition called in void context ";     
}

sub readIndex{
    my ($vcf) = @_;
    if ($vcf =~/\.gz/){
        #if compressed just create index if it doesn't exist and return
        my $index = "$vcf.tbi"; 
        if (not -e $index){
            print STDERR "$index does not exist - indexing $vcf...\n";
            indexVcf($vcf);
            croak "Indexing failed? $index does not exist " if (not -e $index);
            print STDERR " Done.\n";
        }
        return;
    }
    my $index = "$vcf.vridx"; 
    my %contigs = ();
    my $block_dump;
    if (not -e $index){
        print STDERR "$index does not exist - indexing $vcf...\n";
        indexVcf($vcf);
        croak "Indexing failed? $index does not exist " if (not -e $index);
        print STDERR " Done.\n";
    }
    my $z = new IO::Uncompress::Gunzip $index
        or die "gunzip failed to read index $index: $GunzipError\n";
    {
        local $/;
        $block_dump = <$z>;
    }
    close ($z);
    %contigs = %{ eval $block_dump };
    return %contigs;
}

sub _searchVcf{
    my (%args) = @_;
    croak "chrom argument is required for _searchVcf method " if not exists $args{chrom};
    croak "pos argument is required for _searchVcf method " if not exists $args{pos};
    croak "contig_order argument is required for _searchVcf method " if not exists $args{contig_order};
    #croak "index argument is required for _searchVcf method " if not exists $args{index};
    croak "fh argument is required for _searchVcf method " if not exists $args{fh};
    #my $total_lines = exists $args{length} ? $args{length} : get_file_length_from_index($args{fh}, $args{index}); 
    my @matches = ();
    my $total_lines = $args{contig_order}->{last_line};
    my $pos_rounddown = int($args{pos}/$REGION_SPANS) * $REGION_SPANS;
    return if not (exists $args{contig_order}->{$args{chrom}}->{regions}->{$pos_rounddown});
    foreach my $reg (@{$args{contig_order}->{$args{chrom}}->{regions}->{$pos_rounddown}}){
        next if $reg->{pos_start} > $args{pos};
        next if $reg->{pos_end} < $args{pos};
        my @lines = _readLinesByOffset($reg->{offset_start}, $reg->{offset_end}, $args{fh});
        foreach my $l (@lines){
            my @sp = split("\t", $l);
            my $l_pos =  $sp[$vcf_fields{POS}]; 
            last if $l_pos > $args{pos};
            if ($l_pos == $args{pos}){
                push @matches, $l;
                next;
            }
            my $span = $l_pos + length($sp[$vcf_fields{REF}]) -1;
            if ($l_pos <= $args{pos} and $span >= $args{pos}){
                push @matches, $l;
            }
        }
    }
    return @matches;
}

sub _readLinesByOffset{
    my ($start, $end, $fh) = @_;
    my $data = '';
    sysseek ($fh, $start, SEEK_SET);
    sysread($fh, $data, $end - $start, 0);
    return split("\n", $data);
}


sub getFileLengthFromIndex{
    my $vcf = shift;
    my %idx = readIndex($vcf);
    return $idx{last_line};
}

sub countVariants{
    my $vcf = shift;
    croak "countVariants method requires a file as an argument" if not $vcf;
    my @head = getHeader($vcf);
    my $all = getLineCount($vcf);
    return $all - @head;
}


sub sortVariants{
#sort a list of vcf lines
    my ($list, $contig_order) = @_;
    croak "sortVariants required an array reference as an argument " if (ref $list ne 'ARRAY');
    croak "second argument passed to sortVariants must be a hash reference of contigs and their orders " if (ref $contig_order ne 'HASH');
    my %contigs = ();
    my $add_ids = 0;
    my $i = 0;
    if ($contig_order){
        %contigs = %{$contig_order};
    }else{
        my %temp = ();
        foreach my $l (@$list){
            $temp{getVariantField($l, 'CHROM')}++;
        }
        my $n = 0;
        %contigs = map {$_ => $n++} sort byContigs(keys %temp);
    }
    my @sort = sort {
        $contigs{getVariantField($a, 'CHROM')} <=> $contigs{getVariantField($b, 'CHROM')} ||
        getVariantField($a, "POS") <=> getVariantField($b, "POS")
    } @$list;
    return @sort if defined wantarray;
    carp "sortVariants method called in void context ";
}

sub minimizeAlleles{
    #reduce alleles to their simplest representation
    #so that multiallelic variants can be represented 
    #in their most basic form
    #e.g. (from http://www.cureffi.org/2014/04/24/converting-genetic-variants-to-their-minimal-representation/)
    #1  1001    .   CTCC    CCC,C,CCCC 
    #becomes
    #1  1001   CTCC    CCC →   1001    CT  C
    #1001   CTCC    C   →   1001    CTCC    C
    #1001   CTCC    CCCC    →   1002    T   C
    #
    my ($line) = @_;
    my %min_alleles = ();#key is allele number, each entry is anon hash of CHROM, REF, POS, ALT
    my @al =  readAlleles(line => $line);
    for (my $i = 1; $i < @al; $i++){
        my ($pos, $ref, $alt) = reduceRefAlt(getVariantField($line, "POS"), $al[0], $al[$i]);
        $min_alleles{$i} = {
            CHROM           => getVariantField($line, "CHROM"),
            POS             => $pos,
            REF             => $ref,
            ALT             => $alt,
            ORIGINAL_POS    => getVariantField($line, "POS"),
            ORIGINAL_REF    => getVariantField($line, "REF"),
            ORIGINAL_ALT    => getVariantField($line, "ALT"),
        };
    }
    return %min_alleles if defined wantarray;
    carp "minimizeAlleles called in void context ";
}

sub reduceRefAlt{
    #reduce a single ref/alt pair to their simplest representation
    my ($pos, $ref, $alt) = @_;
    if (length($ref) > 1 and length($alt) > 1){
        #can only reduce if both REF and ALT are longer than 1
        my @r = split('', $ref);
        my @al = split('', $alt);
        while ($r[-1] eq $al[-1] and @r > 1 and @al > 1){
            #remove identical suffixes
            pop @r;
            pop @al;
        }
        while ($r[0] eq $al[0] and @r > 1 and @al > 1){
            #remove identical prefixes
            #increment position accordingly
            shift @r;
            shift @al;
            $pos++;
        }
        $ref = join('', @r);
        $alt = join('', @al);
    }
    return ($pos, $ref, $alt);
}


sub byContigs{
    $a =~ s/^chr//;
    $b =~ s/^chr//;
    if ($a =~ /^\d+$/){
        if ($b =~ /^\d+$/){
            return $a <=> $b;
        }else{
            return -1;
        }
    }elsif ($b =~ /^\d+$/){
        return 1;
    }elsif ($a =~ /^[XY]$/){
        if ($b =~ /^[XY]$/){
            return $a cmp $b;
        }else{
            return -1;
        }
    }elsif ($b =~ /^[XY]$/){
        return 1;
    }elsif ($a =~ /^MT*$/){
        return $b cmp $a;
    }elsif ($b =~ /^MT*$/){
        return 1;
    }else{
        return $a cmp $b;
    }
}


