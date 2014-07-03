package VcfReader;
use strict;
use warnings;
use Carp;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
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
    return @header;
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
    croak "getHeader method requires a file as an argument" if not $vcf;
    my $contig_index = "$vcf.vrdict";
    my %contigs = ();
    if (not -e $contig_index){
        print STDERR "indexing $vcf... ";
        indexVcf($vcf);
        print STDERR "Done.\n";
        croak "Indexing failed? $contig_index does not exist " if not -e $contig_index;
    }
    open (my $CONTIG, $contig_index) or croak "Can't open $contig_index for reading: $! ";
    while (my $line = <$CONTIG>){
        chomp $line;
        croak "Contig index $contig_index is corrupt - extraneous white space found - please delete and re-index" if $line =~ /\s/;
        if (exists $contigs{$line}){
            if ($contigs{$line}  != scalar(keys%contigs) -1 ){
                croak "Contig index $contig_index is corrupt - duplicate entries ";
            }
        }else{
            $contigs{$line}  = scalar(keys%contigs);
        }
    }
    close $CONTIG;
    return %contigs;
}

#file utilities

sub getLineCount{
    my $vcf = shift; 
    croak "getLineCount method requires a file as an argument" if not $vcf;
    my $line_count = 0;
    my $FH = _openFileHandle($vcf);
    if ($vcf =~ /\.gz$/){
        $line_count ++ while (<$FH>);
    }else{
        $line_count += tr/\n/\n/ while sysread($FH, $_, 2 ** 16);
    }
    close $FH;
    return $line_count;
}


sub indexVcf{
    my ($vcf) = @_;
    if ($vcf =~ /\.gz$/){
        chomp (my $tabix = `which tabix`);
        croak "Can't find tabix executable to index compressed VCF.  Please ensure it is ".
            "installed and in your PATH or index your VCF with tabix manually. "
            if not $tabix; 
        print STDERR "Indexing $vcf with tabix.\n";
        my $er = `$tabix -p vcf $vcf 2>&1`;
        croak "Tabix indexing failed, code $?: $er\n" if $?;
        return;
    }
    my $index = "$vcf.vridx";
    open (my $INDEX, "+>$index") or croak "can't open $index for writing index: $! ";
    my $offset = 0;
    my $contig_index = "$vcf.vrdict";
    open (my $CONTIGS, ">$contig_index") or croak "Can't open $contig_index file for writing contig order: $! ";
    my %contigs = ();
    my $pack = "N";
    $pack = "Q" if -s $vcf >= 4294967296;
    my $prev_pos = 0;
    my $FH = _openFileHandle($vcf);
    while (my $line = scalar readline $FH){
        print $INDEX pack($pack, $offset);
        $offset = tell $FH;
        next if $line =~ /^#/;
        my $chrom = (split "\t", $line)[$vcf_fields{CHROM}];
        my $pos = (split "\t", $line)[$vcf_fields{POS}];
        if (exists $contigs{$chrom}){
            if ($contigs{$chrom}  != scalar(keys%contigs) -1 or $pos < $prev_pos){
                close $INDEX;
                unlink $index or carp "Couldn't delete partial index $index - please delete manually: $!" ;
                croak "Can't index VCF - $vcf not sorted properly ";
            }
        }else{
            $contigs{$chrom}  = scalar(keys%contigs);
            print $CONTIGS "$chrom\n";
        }
        $prev_pos = $pos;
    }
    close $CONTIGS;
    close $INDEX;
    close $FH;
}

sub _lineWithIndex{
    my ($FH, $IDX, $line_number) =@_;
    my $size;               # size of an index entry
    my $i_offset;           # offset into the index of the entry
    my $entry;              # index entry
    my $d_offset;           # offset into the data file
    my $pack = "N";
    $pack = "Q" if -s $FH >= 4294967296;
    $size = length(pack($pack, 0));
    $i_offset = $size * ($line_number-1);
    seek($IDX, $i_offset, 0) or return;
    read($IDX, $entry, $size);
    $d_offset = unpack($pack, $entry);
    seek($FH, $d_offset, 0);
    return scalar readline $FH;
}

sub _openFileHandle{
    my $vcf = shift;
    croak "_openFileHandle method requires a file as an argument" if not $vcf;
    my $FH; 
    if ($vcf =~ /\.gz$/){
        $FH = new IO::Uncompress::Gunzip $vcf, MultiStream => 1 or croak "IO::Uncompress::Gunzip failed while opening $vcf for reading: \n$GunzipError";
    }else{
        open (my $FH, $vcf) or croak "Failed to open $vcf for reading: $! ";
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
        if (defined $args{sample_to_columns}){
            $var = getSampleVariant(
                                $args{line}, 
                                $args{sample_to_columns}->{$args{sample}},
                                );
        }else{
            $var = getSampleVariant(
                                $args{line}, 
                                $args{sample},
                                );
        }
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
    }elsif(defined $args{sample} or $args{multiple}){
        croak "\"multiple\" and \"sample\" arguments can only be used in conjunction with \"sample_to_columns\" option for getSampleCall method ";
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
    }elsif (defined $args{sample}){
        my $call;
        if (exists $args{minGQ} and $args{minGQ} > 0){
            $call = _getGenotype(
                        $args{line},
                        $args{sample_to_columns}->{$args{sample}}, 
                        $args{minGQ}
                        );
        }else{
            $call = _getGenotype(
                        $args{line},
                        $args{sample_to_columns}->{$args{sample}}, 
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
    }else{#otherwise just look at first sample
        my $call;
        if (exists $args{minGQ} and $args{minGQ} > 0){
            $call = _getGenotype(
                        $args{line},
                        9,
                        $args{minGQ}
                        );
        }else{
            $call = _getGenotype(
                        $args{line},
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

=cut

#returns genotype codes - e.g. 0/0, 0/1, 1/1
sub getAllPossibleGenotypeCodes{
    my ($self, $line) = @_;
    my @alleles = $self -> readAlleles($line);
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
    my ($self, $line) = @_;
    my @alleles = $self -> readAlleles($line);
    my @combinations = ();
    for (my $n = 0; $n < @alleles; $n++){
        for (my $m = 0; $m <= $n; $m++){
            push (@combinations, "$alleles[$m]/$alleles[$n]");
        }
    }
    return @combinations if defined wantarray;
    carp "getAllPossibleGenotypes called in void context ";
}

sub getSampleActualGenotypes{
    my ($self, %args) = @_;
    croak "Can't invoke getSampleActualGenotypes method when no samples/genotypes are present in VCF " if not defined $self->{_samples};
    croak "\"line\" argument must be passed to getSampleActualGenotypes " if not defined $args{line};
    carp "WARNING Both multiple and sample arguments supplied to getSampleActualGenotypes method - only multiple argument will be used " if (defined $args{multiple} and defined $args{sample});
    my @split = ();
    if (ref $args{line} eq 'ARRAY'){
        @split = @{$args{line}};
    }else{
        @split = split("\t", $args{line});
    }
    my @alleles = $self -> readAlleles(line => \@split);
    my %multiple = ();
    my $genotype;
    my @sample_alleles = ();
    my $var;
    if ($args{all}){
        foreach my $sample (keys %{$self->{_samples}}){
            $var = $self->getSampleVariant(\@split, $sample);
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
        if ($args{return_alleles_only}){
            my %seen = ();
            @sample_alleles = grep {!$seen{$_}++} @sample_alleles;#remove duplicates
            return @sample_alleles;
        }else{
            return %multiple;
        }
    }elsif($args{multiple}){
        croak "multiple argument must be an array reference " if ref $args{multiple} ne 'ARRAY';
        foreach my $sample (@{$args{multiple}}){
            $var = $self->getSampleVariant(\@split, $sample);
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
        if ($args{return_alleles_only}){
            my %seen = ();
            @sample_alleles = grep {!$seen{$_}++} @sample_alleles;#remove duplicates
            return @sample_alleles;
        }else{
            return %multiple;
        }
    }elsif ($args{sample}){
        $var = $self->getSampleVariant(\@split, $args{sample});
    }else{
        $var = $self->getSampleVariant(\@split);
    }
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

sub replaceVariantField{
    my ($self, $line, $field, $replacement) = @_;
    $field =~ s/^#+//;
    $self->{_currentLine} ||= $self->readLine; #get line if we haven't already
    my @valid_fields = $self->getValidFields;#qw(CHROM POS ID REF ALT QUAL FILTER INFO FORMAT);
    croak "Invalid field ($field) passed to getVariantField method " if not grep {$field eq $_} @valid_fields;
    my @split = ();
    if (ref $line eq 'ARRAY'){
        @split = @$line;
    }else{
        @split = split("\t", $line);
    }
    my @new_line = ();
    foreach my $f ( split ("\t", $self->{_header})){
        $f =~ s/^#+//;
        if ($f eq $field){
            push @new_line, $replacement;
        }else{
            push @new_line, $self->getCustomField(\@split, $f);
        }
    }
    return join("\t", @new_line) if defined wantarray;
    carp "replaceVariantField called in void context ";
}

sub searchForPosition{
#my @matching lines = $obj->searchForPosition(chrom => 1, pos = 2000000)
    my ($self, %args) = @_;
    croak "Can't use searchForPosition method when noLineCount option is set and input is not bgzip compressed " if $self->{_noLineCount} and not $self->{_file} =~ /\.gz$/;
    croak "Can't use searchForPosition method when reading from STDIN" if $self->{_inputIsStdin} ;
    carp "chrom argument is required for searchForPosition method " if not exists $args{chrom};
    carp "pos argument is required for searchForPosition method " if not exists $args{pos};
    return if not exists $args{pos} or not exists $args{chrom};
    #if compressed assume bgzip compression and try to read with tabix
    if ($self->{_file} =~ /\.gz$/){
        eval "use Tabix; 1" 
            or croak "Tabix module is not installed and VCF file $self->{_file} appears to be (b)gzip compressed.  ".
            "  Please install Tabix.pm in order to search compressed VCFs.\n";
        if (not -e "$self->{_index}"){
            $self->setTabix if not $self->{_tabix}; 
            print STDERR "Indexing $self->{_file} with tabix.\n";
            my $er = `$self->{_tabix} -p vcf $self->{_file} 2>&1`;
            croak "Tabix indexing failed, code $?: $er\n" if $?;
        }

#FIX THIS!
#DO WE NEED ONE TABIX ITERATOR PER POTENTIAL THREAD?!

        if (not $self->{_tabixIterator}){
            $self->{_tabixIterator} = Tabix->new(-data =>  $self->{_file}, -index => $self->{_index});
        }
        my $iter = $self->{_tabixIterator}->query($args{chrom}, $args{pos} -1, $args{pos});
        return if not defined $iter->{_}; #$iter->{_} will be undef if our chromosome isn't in the vcf file
        my @matches;
        while (my $result = $self->{_tabixIterator} ->read($iter)){ 
            push @matches, $result;
        }
        $self->{_searchMatches} = \@matches;
        return @matches if defined wantarray;
        return;
    }

    #otherwise use our own line indexing method
    if (not $self->{_indexHandle}){
        if (not -e $self->{_index}){
            $self->{_indexHandle} = $self->indexVcf;
        }else{
            open ($self->{_indexHandle}, $self->{_index}) or croak "Can't open index file $self->{_index} for reading: $! ";
        }
    }
    if (not $self->{_contigOrder}){
        $self->{_contigIndex} ||= $self->{_file}.".pvdict";
        if (not -e $self->{_contigIndex}){
            $self->{_indexHandle} = $self->indexVcf;
        }else{
            open (my $CONTIGS, $self->{_contigIndex}) or croak "Can't open contig index $self->{_contigIndex} for reading: $! ";
            my %contigs = ();
            while (my $line = <$CONTIGS>){
                chomp $line;
                croak "Contig index $self->{_contigIndex} is corrupt - extraneous white space found " if $line =~ /\s/;
                if (exists $contigs{$line}){
                    if ($contigs{$line}  != scalar(keys%contigs) -1 ){
                        croak "Contig index $self->{_contigIndex} is corrupt - duplicate entries ";
                    }
                }else{
                    $contigs{$line}  = scalar(keys%contigs);
                }
            }
            close $CONTIGS;
            croak "Could not find any contigs in contig index ($self->{_contigIndex}). Try deleting $self->{_contigIndex} and rerunning "  if (not %contigs); 
            $self->{_contigOrder} = \%contigs;
        }
    }
    my @matches = $self->_searchVcf(chrom => $args{chrom}, pos => $args{pos});
    $self->{_searchMatches} = \@matches;
    return @matches if defined wantarray;
}

sub _binSearch{
    my ($self, %args) = @_;
    my $u = $self->{_totalLines};
    my $l = 1;
    return 0 if not exists $self->{_contigOrder}->{$args{chrom}};
    my $i = 0;
    while ($l <= $u){
        $i = int(($l + $u)/2);
        my $line =  $self->_lineWithIndex($i);
        if ($line =~ /^#/){
            $l = $i+1;
            next;
        }
        my @split = split("\t", $line);
        if ($self->{_contigOrder}->{$args{chrom}}  < $self->{_contigOrder}->{$split[$self->{_fields}->{CHROM}]} or
            ($self->{_contigOrder}->{$args{chrom}} == $self->{_contigOrder}->{$split[$self->{_fields}->{CHROM}]} and $args{pos} < $split[$self->{_fields}->{POS}]) ){
                $u = $i -1;
        }elsif($self->{_contigOrder}->{$args{chrom}}  > $self->{_contigOrder}->{$split[$self->{_fields}->{CHROM}]} or
            ($self->{_contigOrder}->{$args{chrom}} == $self->{_contigOrder}->{$split[$self->{_fields}->{CHROM}]} and $args{pos} > $split[$self->{_fields}->{POS}]) ){
                $l = $i + 1;
        }else{
            return $i;
        }
    }
    #FIND OVERLAPPING VARIANTS (DELETIONS/MNVs) NOT NECESSARILY OF SAME COORDINATE
    for (my $j = $i - 1; $j > 0; $j--){#look at previous lines
        my $line =  $self->_lineWithIndex($j);
        last if ($line =~ /^#/);
        my @split = split("\t", $line);
        #skip if we're at the next chrom...
        next if 
            $self->{_contigOrder}->{$args{chrom}}  < $self->{_contigOrder}->{$split[$self->{_fields}->{CHROM}]};
        #...or downstream
        next if 
            $self->{_contigOrder}->{$args{chrom}} == $self->{_contigOrder}->{$split[$self->{_fields}->{CHROM}]}
	            and $args{pos} < $split[$self->{_fields}->{POS}];
        #we're done if we've got to the previous chromosome...
        last if 
            $self->{_contigOrder}->{$args{chrom}}  > $self->{_contigOrder}->{$split[$self->{_fields}->{CHROM}]};
        #and let's assume we're not going to have detected any del/mnv larger than 200 bp (is there a better way?)
        last if 
            $self->{_contigOrder}->{$args{chrom}} == $self->{_contigOrder}->{$split[$self->{_fields}->{CHROM}]} 
                and $args{pos} < ($split[$self->{_fields}->{POS}] - 200);
        my $ref_length = length($split[$self->{_fields}->{REF}]);
        if ($ref_length > 1 and ($ref_length + $split[$self->{_fields}->{POS}] - 1) >= $args{pos}){
            #we've found a deletion/mnv that overlaps our pos
            return $j;
        }
    }
    return 0;
}

sub _searchVcf{
    my ($self, %args) = @_;
    croak "chrom argument is required for _searchVcf method " if not exists $args{chrom};
    croak "pos argument is required for _searchVcf method " if not exists $args{pos};
    my $i = $self->_binSearch(chrom=>$args{chrom}, pos => $args{pos});
    if ($i){
        my @hits;
        for (my $j = $i - 1; $j > 0; $j--){#look at previous lines
            my $line = $self->_lineWithIndex($j);
            last if ($line =~ /^#/);
            my @split = split("\t", $line);
            if ($split[$self->{_fields}->{CHROM}] eq $args{chrom} and $split[$self->{_fields}->{POS}] == $args{pos}){
                push @hits, $line;
            }else{
                #FIND OVERLAPPING VARIANTS (DELETIONS/MNVs) NOT NECESSARILY OF SAME COORDINATE
                #we're done if we've got to the previous chromosome...
                last if 
                    $self->{_contigOrder}->{$args{chrom}}  > $self->{_contigOrder}->{$split[$self->{_fields}->{CHROM}]};
                #and let's assume we're not going to have detected any del/mnv larger than 200 bp (is there a better way?)
                last if 
                    $self->{_contigOrder}->{$args{chrom}} == $self->{_contigOrder}->{$split[$self->{_fields}->{CHROM}]} 
                        and $args{pos} < ($split[$self->{_fields}->{POS}] - 200);
                my $ref_length = length($split[$self->{_fields}->{REF}]);
                if ($ref_length > 1 and ($ref_length + $split[$self->{_fields}->{POS}] - 1) >= $args{pos}){
                    #we've found a deletion/mnv that overlaps our pos
                    push @hits, $line;
                }
            }
        }
        @hits = reverse(@hits);#put in original order
        push @hits, $self->_lineWithIndex($i);#add original hit
        for (my $j = $i + 1; $j <= $self->{_totalLines}; $j++){#look at next lines
            my $line = $self->_lineWithIndex($j);
             my @split = split("\t", $line);
            if ($split[$self->{_fields}->{CHROM}] eq $args{chrom} and $split[$self->{_fields}->{POS}] == $args{pos}){
                push @hits, $line;
            }else{
                last;
            }
        }
        return @hits;
    }else{
        return;
    }
}
=cut

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


