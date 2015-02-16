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

package ParseVCF;

use strict;
use warnings;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use Data::Dumper;
use Carp;
use Scalar::Util qw (looks_like_number);
use List::Util qw(sum);
our $AUTOLOAD;
{
        my $_count = 0;
        my %_attrs = (
        _file => ["", "read/required"],
        _header => ["", "read/write"],
        _metaHeader => [[], "read/write"],
        _variantCount => [0, "read"],
        _totalLines => [0, "read"],
        _inputIsStdin => [0, "read"],
       # _sampleOrder => [[], "read"],
        _currentLine => ["", "read"],
        _index => ["", "read"],
        _noHeaderCheck => [0, "read/write"],
        _noLineCount => [0, "read/write"],
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
        close $self->{_filehandle} if $self->{_filehandle};
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
    if ($self->{_file} eq '-'){
        $self->{_inputIsStdin} = 1;
        $self->{_filehandle} = $self -> _openFileHandle();
    }else{
        $self->{_filehandle} = $self -> _openFileHandle();
        $self->{_totalLines} = $self -> _getLineCount() unless $self->{_noLineCount};#get line count
        $self -> reopenFileHandle();#close and repopen file
        if ($self->{_file} =~ /\.gz$/){
            $self->{_index} = $self->{_file}.".tbi";
        }else{
            $self->{_index} = $self->{_file}.".pvidx";
        }
        $self->{_fileSize} = -s $self->{_file};
    }
    $self->{_validFields} = $self->_getValidFields();
    #if user has supplied header seperately
    if ($self->{_header}){
        my $success = $self->changeHeader(string=>$self->{_header});
        if (not $success){
            croak "Could not initialise with user supplied header ";
        }
    }else{
        $self -> _readHeader();
        $self -> _readHeaderInfo();
    }
    $class -> _incr_count();
    return $self;
}


sub _getValidFields{
    my ($self) = @_;
    my %valid_fields = map {$_ => undef } qw(CHROM POS ID REF ALT QUAL FILTER INFO FORMAT);
    return \%valid_fields;
}

sub _getLineCount{
    my ($self) = shift;
    my $line_count = 0;
    $line_count += tr/\n/\n/ while sysread($self->{_filehandle}, $_, 2 ** 16);
    return $line_count;
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

sub reopenFileHandle{
    my ($self) = @_;
    croak "Error - can't use reopenFileHandle method on input from STDIN!\n" if ($self->{_file} eq '-');
    close $self->{_filehandle};
    $self->{_filehandle} = $self -> _openFileHandle();
    return $self->{_filehandle} if defined wantarray;
}

sub _openFileHandle{
    my ($self) = @_;
    my $FH;
    croak "file argument is not defined! " if not defined $self->{_file};
    if ($self->{_file} =~ /\.gz$/){
        $FH = new IO::Uncompress::Gunzip $self->{_file}, MultiStream => 1 or croak "IO::Uncompress::Gunzip failed while opening $self->{_file} for reading: \n$GunzipError";
    }else{
        open ($FH, $self->{_file}) or croak "Can't open VCF file $self->{_file}: $! ";
    }
    return $FH;
}

sub _byContigs{
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

sub sortVcf{
#creates a new vcf sorted in coordinate order
#returns filename of new sorted file
#does not alter $self
    my ($self, %args) = @_;
    my %contigs = ();
    my $add_ids = 0;
    my $i = 0;
    if (not $self->{_contigOrder}){
        foreach my $head (@{$self->{_metaHeader}}){
            if ($head =~ /##contig=<ID=([^,>]+),\S+/){
                $contigs{$1} = $i++;
            }
        }
        if (not keys(%contigs)){
            print "WARNING - no contig lines found in header for $self->{_file}, will gather contigs and sort arbitrarily.\n";
            $add_ids++;
            $self->reopenFileHandle();#ensure we're at the start of the file
            my %chroms = ();
            while (my $line = scalar readline $self->{_filehandle}){
                next if $line =~ /^#/;
                $chroms{(split "\t", $line)[$self->{_fields}->{CHROM}]}++;
            }
            $i = 0;
            %contigs = map {$_ => $i++} sort _byContigs(keys%chroms);
        }
        croak "No contigs found in $self->{_file}, cannot sort\n" if not keys(%contigs);
        $self->{_contigOrder} = \%contigs;
    }else{
        %contigs = %{$self->{_contigOrder}};
    }
    my $SORTOUT;
    if (exists $args{output}){
        open ($SORTOUT,  ">$args{output}") or croak "Can't open file $args{output} for output of sortVcf: $! ";
        print STDERR "Sorting $self->{_file} to $args{output}.\n"; 
    }else{
        $SORTOUT = \*STDOUT;
    }
    eval "use Sort::External; 1" or carp "The Sort::External module was not found - will attempt to sort in memory. For huge files it is recommended to install Sort::External via CPAN.\n";
    $self->reopenFileHandle();#ensure we're at the start of the file
    if ($@){#no Sort::External , sort in memory
        print STDERR "Reading variants into memory...\n";
        my @sort = readline $self->{_filehandle};
        @sort = grep { ! /^#/ } @sort;
        print STDERR "Performing sort...\n";
        @sort = sort {$contigs{(split "\t", $a)[$self->{_fields}->{CHROM}]} <=> $contigs{(split "\t", $b)[$self->{_fields}->{CHROM}]} || 
                        (split "\t", $a)[$self->{_fields}->{POS}] <=> (split "\t", $b)[$self->{_fields}->{POS}]
                    } @sort;
        print STDERR "Printing output...";
        print $SORTOUT $self->getHeader();
        if ($add_ids){
            foreach my $c (sort {$contigs{$a} <=> $contigs{$b} } keys %contigs){
                print $SORTOUT "##contig=<ID=$c>\n";
            }
        }
        print $SORTOUT join('', @sort);
    }else{
        my $vcfsort = 
            sub { 
                $contigs{(split "\t", $Sort::External::a)[$self->{_fields}->{CHROM}]} <=> $contigs{(split "\t", $Sort::External::b)[$self->{_fields}->{CHROM}]} || 
                (split "\t", $Sort::External::a)[$self->{_fields}->{POS}] <=> (split "\t", $Sort::External::b)[$self->{_fields}->{POS}]
            };
        my $sortex = Sort::External->new(mem_threshold => 1024**2 * 24, 
                                        sortsub => $vcfsort);
        my @feeds = ();
        my $n = 0;
        while (my $line = scalar readline $self->{_filehandle}){
            $n++;
            if ($line =~ /^##/){
                print $SORTOUT $line;
                next;
            }elsif ($line =~ /^#/){
                if ($add_ids){
                    foreach my $c (sort {$contigs{$a} <=> $contigs{$b} } keys %contigs){
                        print $SORTOUT "##contig=<ID=$c>\n";
                    }
                }
                print $SORTOUT $line;
                next;
            }
            push @feeds, $line;
            if (@feeds > 99999){
                $sortex->feed(@feeds);
                @feeds = ();
                print STDERR "Fed $n lines";
                print STDERR " of " . ($self->{_totalLines})  if $self->{_totalLines};
                print STDERR " to sort...\n";
            }
        }
        $sortex->feed(@feeds) if @feeds;
        print STDERR "Fed $n lines";
        print STDERR " of " . ($self->{_totalLines})  if $self->{_totalLines};
        print STDERR ".\n";
        print STDERR "Sorting and writing output...";
        $sortex->finish; 
        $n = 0;
        while ( defined( $_ = $sortex->fetch ) ) {
            print $SORTOUT $_;
            $n++;
        }
    }
    if (exists $args{output}){
        close $SORTOUT or croak "Couldn't finish writing to sort output file $args{output}: $!\n" ;
    }
    print STDERR " Done.\n";
    $self->reopenFileHandle();#reset file
    return $args{output} if defined wantarray
}

sub sortVariants{
#sort a list of vcf lines
    my ($self, $list) = @_;
    croak "sortVariants required an array reference as an argument " if (ref $list ne 'ARRAY');
    my %contigs = ();
    my $add_ids = 0;
    my $i = 0;
    if (not $self->{_contigOrder}){
        foreach my $head (@{$self->{_metaHeader}}){
            if ($head =~ /##contig=<ID=([^,>]+),\S+/){
                $contigs{$1} = $i++;
            }
        }
        if (not keys(%contigs)){
            print "WARNING - no contig lines found in header for $self->{_file}, will gather contigs and sort arbitrarily.\n";
            $add_ids++;
            $self->reopenFileHandle();#ensure we're at the start of the file
            my %chroms = ();
            while (my $line = scalar readline $self->{_filehandle}){
                next if $line =~ /^#/;
                $chroms{(split "\t", $line)[$self->{_fields}->{CHROM}]}++;
            }
            $i = 0;
            %contigs = map {$_ => $i++} sort _byContigs(keys%chroms);
        }
        croak "No contigs found in $self->{_file}, cannot sort\n" if not keys(%contigs);
        $self->{_contigOrder} = \%contigs;
    }else{
        %contigs = %{$self->{_contigOrder}};
    }
    my @sort = sort {$contigs{(split "\t", $a)[$self->{_fields}->{CHROM}]} <=> $contigs{(split "\t", $b)[$self->{_fields}->{CHROM}]} || 
                        (split "\t", $a)[$self->{_fields}->{POS}] <=> (split "\t", $b)[$self->{_fields}->{POS}]
                    } @$list;
    return @sort if defined wantarray;
    carp "sortVariants method called in void context ";
}

sub altsToVepAllele{
    my ($self, %args) = @_;
    my @vep_alleles = ();
    my $ref ;
    if ($args{ref}){
        $ref = $args{ref};
    }else{
        $ref = $self->getVariantField('REF');
    }
    if (not $args{alt}){
        foreach my $alt ($self->readAlleles(alt_alleles => 1)){
            push @vep_alleles, _altToVep($alt, $ref);
        }
    }else{
        if (ref $args{alt} eq 'ARRAY'){
            foreach my $alt (@{$args{allele}}){
                push @vep_alleles, _altToVep($alt, $ref);
            }
        }else{
            return _altToVep($args{alt}, $ref) if defined wantarray;
        }
    }
    return @vep_alleles if defined wantarray;
    carp "altsToVepAllele called in void context. ";
}
                
sub _altToVep{
    my ($alt, $ref) = @_;
    if (length($alt) == length($ref)){
        return $alt;
    }elsif(length($alt) > length($ref)){#insertion - VEP trims first base
       return substr($alt, 1);
    }else{#deletion - VEP trims first base or gives '-' if ALT is only 1 nt long
        if (length($alt) > 1){
            return substr($alt, 1);
        }else{
            return '-';
        }
    }
}

sub readVepHeader{
    my ($self) = @_;
    if (not @{$self->{_metaHeader}}){
        croak "Method 'readVepHeader' requires meta header but no meta header lines found ";
    }
    my @info = grep{/^##INFO=<ID=CSQ/} (@{$self->{_metaHeader}});
    if (not @info){
        croak "Method 'readVepHeader' requires CSQ INFO field in meta header (e.g. '##INFO=<ID=CSQ,Number...') but no matching lines found ";
    }
    carp "Warning - multiple CSQ fields found, ignoring all but the most recent field " if @info > 1;
    my $csq_line = $info[-1] ;#assume last applied VEP consequences are what we are looking for 
    my @csq_fields = ();    
    if ($csq_line =~ /Format:\s(\S+\|\S+)">/){
        @csq_fields = split(/\|/, $1);
    }else{
        croak "Method 'readVepHeader' couldn't properly read the CSQ format from the corresponding INFO line: $csq_line ";
    }
    if (not @csq_fields){
        croak "Method 'readVepHeader' didn't find any VEP fields from the corresponding CSQ INFO line: $csq_line ";
    }
    for (my $i = 0; $i < @csq_fields; $i++){
        $self->{_vepCsqHeader}->{lc($csq_fields[$i])} = $i;
    }
    #from v73 of VEP 'hgnc' has been replaced with 'symbol'
    #the code below is aimed to maintain compatibility between versions
    if (grep {/^hgnc$/i} @csq_fields){
        if (not grep {/^symbol$/i} @csq_fields){
            $self->{_vepCsqHeader}->{symbol} = $self->{_vepCsqHeader}->{hgnc};
        }
    }elsif(grep {/^symbol$/i} @csq_fields){
        if (not grep {/^hgnc$/i} @csq_fields){
            $self->{_vepCsqHeader}->{hgnc} = $self->{_vepCsqHeader}->{symbol};
        }
    }
    $self->{_vepCsqArrayField} = \@csq_fields;
    return $self->{_vepCsqHeader} if defined wantarray;
}

sub vepFieldToString{
#take single vep hash and order it according to vep header - useful if splitting vep fields to sep lines
    my ($self, $csq) = @_;
    if (not $self->{_vepCsqArrayField}){
        $self->readVepHeader();
    }
    my @string = ();
    foreach my $field (@{$self->{_vepCsqArrayField}}){
        if (exists $csq->{lc$field}){
            $csq->{lc$field} ? push @string, $csq->{lc$field} : push @string, "";
        }
    }
    return join("|", @string);
}

sub getVepFields{
#if more than one VEP annotation has been performed this will only look at the field corresponding to the last 
#applied annotation
    #ALLELE FIELD IS REQUIRED
    #my @all_genes =  $obj->getVepFields("Gene")
    #my @all_gene_cons_pph = $obj->getVepFields([Gene, Consequence, PolyPhen]); # returns an array of hashes with the keys being the field 
    #my @everything = $obj->getVepFields("All");
    #    e.g. $all_gene_cons_pph[0] -> {Gene} = ENSG000012345

    my ($self, $fields, $no_warn) = @_;
    if (not $fields){
        croak "Method 'getVepFields' requires either an ARRAY reference or a scalar value as an argument";
    }elsif(ref $fields && ref $fields ne 'ARRAY'){
        croak "Method 'getVepFields' takes either an ARRAY reference or a scalar value as an argument, not a " . ref $fields ." reference ";
    }
    if (not $self->{_vepCsqHeader}){#get VEP CSQ format from header if we haven't already
        $self->readVepHeader();
    }
    my @info_field = split(';', $self->getVariantField("INFO"));
    my @csqs = grep{/^CSQ/} @info_field;
    if (not @csqs){
        carp "No CSQ field found in INFO field for line: $self->{_currentLine} " unless $no_warn;
        return;
    }
    my @vep = split(",", $csqs[-1]);
    my @return = ();
    if (not ref $fields){
        foreach my $v (@vep){
            my @v_part = split(/[\|]/, $v);
            if (lc$fields eq 'all'){
                my %consequence = ();
                foreach my $f (keys %{$self->{_vepCsqHeader}}){
                    my $value = $v_part[$self->{_vepCsqHeader}->{$f}];
                    $value =~ s/^CSQ=//i if $value;
                    $consequence{lc$f} = $value;
                }
                push @return, \%consequence;
            }else{
                if (exists $self->{_vepCsqHeader}->{lc($fields)}){
                    my $value = $v_part[$self->{_vepCsqHeader}->{lc($fields)}] ;
                    $value =~ s/^CSQ=//i;
                    push @return, $value;
                }else{
                    carp "$fields feature does not exist in CSQ field " unless $no_warn;
                }
            }
        }
    }elsif (ref $fields eq 'ARRAY'){
        my %warned = ();
        foreach my $v (@vep){
            my %consequence = ();
            my @v_part = split(/\|/, $v);
            foreach my $f (@$fields){
                if (exists $self->{_vepCsqHeader}->{lc($f)}){
                    my $value = $v_part[$self->{_vepCsqHeader}->{lc($f)}];
                    $value =~ s/^CSQ=//i if $value;
                    $consequence{lc$f} = $value;
                }else{
                    carp "$f feature does not exist in CSQ field " unless $warned{lc$f} or $no_warn;
                    $consequence{lc$f} = '';
                    $warned{lc$f}++;
                }
            }
            push @return, \%consequence;
        }
    }
    return @return if wantarray;
    return \@return if defined wantarray;
    carp "getVepFields called in void context ";
}

#implementation of readVepHeader but for snpEff annotations
sub readSnpEffHeader{
    my ($self) = @_;
    if (not @{$self->{_metaHeader}}){
        croak "Method 'readSnpEffHeader' requires meta header but no meta header lines found ";
    }
    my @info = grep{/^##INFO=<ID=EFF/} (@{$self->{_metaHeader}});
    if (not @info){
        croak "Method 'readSnpEffHeader' requires EFF INFO field in meta header (e.g. '##INFO=<ID=EFF,Number...') but no matching lines found ";
    }
    carp "Warning - multiple EFF fields found, ignoring all but the most recent field " if @info > 1;
    my $eff_line = $info[-1] ;#assume last applied VEP consequences are what we are looking for 
    my @eff_fields = ();  
    if ($eff_line =~ /Format:\s+'Effect\s+\(\s+(.*)\)'\s*">/){
        @eff_fields = split(/\s*\|\s*/, $1);
        unshift @eff_fields, 'Effect';
    }else{
        croak "Method 'readSnpEffHeader' couldn't properly read the EFF format from the corresponding INFO line: $eff_line ";
    }
    if (not @eff_fields){
        croak "Method 'readSnpEffHeader' didn't find any snpEff fields from the corresponding EFF INFO line: $eff_line ";
    }

    for (my $i = 0; $i < @eff_fields; $i++){
        $eff_fields[$i] =~ s/[\s\[\]]//g;
        $self->{_snpEffHeader}->{$eff_fields[$i]} = $i;
    }

    $self->{_snpEffArrayField} = \@eff_fields;
    return @eff_fields if wantarray;#return an array of fields in order
    return $self->{_snpEffHeader} if defined wantarray;
}

#implementation of getVepFields but for snpEff annotations
sub getSnpEffFields{
    my ($self, $fields) = @_;
    if (not $fields){
        croak "Method 'getSnpEffFields' requires either an ARRAY reference or a scalar value as an argument";
    }elsif(ref $fields && ref $fields ne 'ARRAY'){
        croak "Method 'getSnpEffFields' takes either an ARRAY reference or a scalar value as an argument, not a " . ref $fields ." reference ";
    }
    if (not $self->{_snpEffHeader}){#get EFF format from header if we haven't already
        $self->readSnpEffHeader();
    }
    my @info_field = split(';', $self->getVariantField("INFO"));
    my @csqs = grep{/^EFF/} @info_field;
    if (not @csqs){
        carp "No EFF field found in INFO field for line: $self->{_currentLine} ";
        return;
    }
    $csqs[-1]  =~ s/^EFF=//i ;
    my @eff = split(",", $csqs[-1]);
    my @return = ();
    if (not ref $fields){
        foreach my $e (@eff){
            if ($e =~ /^(\S+)\((\S+)\)$/){
                my @e_part = ($1, split(/\|/, $2));
                if (lc$fields eq 'all'){
                    my %consequence = ();
                    foreach my $f (keys %{$self->{_snpEffHeader}}){
                        my $value = $e_part[$self->{_snpEffHeader}->{$f}];
                        $consequence{$f} = $value;
                    }
                    push @return, \%consequence;
                }else{
                    if (exists $self->{_snpEffHeader}->{$fields}){
                        my $value = $e_part[$self->{_snpEffHeader}->{$fields}] ;
                        push @return, $value;
                    }else{
                        carp "$fields feature does not exist in EFF field ";
                    }
                }
            }else{
                carp "Unexpected format for SnpEff entry - ignored: $e\n";
            }
        }
    }elsif (ref $fields eq 'ARRAY'){
        my %warned = ();
        foreach my $e (@eff){
            if ($e =~ /^(\S+)\((\S+)\)$/){
                my @e_part = ($1, split(/\|/, $2));
                my %consequence = ();
                foreach my $f (@$fields){
                    if (exists $self->{_snpEffHeader}->{$f}){
                        my $value = $e_part[$self->{_snpEffHeader}->{$f}];
                        $consequence{$f} = $value;
                    }else{
                        carp "$f feature does not exist in EFF field " unless $warned{$f};
                        $consequence{$f} = '';
                        $warned{$f}++;
                    }
                }
                push @return, \%consequence;
            }else{
                carp "Unexpected format for SnpEff entry - ignored: $e\n";
            }
        }
    }
    return @return if wantarray;
    return \@return if defined wantarray;
    carp "getSnpEffFields called in void context ";
}

#this method is required for finding tabix to allow
#indexing of vcf if tabix index is not present
#region retrievals are done using Tabix.pm module
sub setTabix{
    my ($self, $tabix) = @_;
    if ($tabix){
        if (not -e $tabix){
            croak "$tabix passed to setTabix method does not exist!\n";
        }
        elsif (not -x $tabix){
            croak "$tabix passed to setTabix method is not executable!\n";
        }else{
            $self->{_tabix} = $tabix;
        }
    }else{
        chomp ($self->{_tabix} = `which tabix`);
        croak "Can't find tabix executable to index compressed VCF.  Please ensure it is ".
            "installed and in your PATH or index your VCF with tabix manually.\n" 
            if not $self->{_tabix};
    }
    return $self->{_tabix} if defined wantarray;
}


sub readPosition{
#first use:
#     if ($obj->searchForPosition(chrom => 1, pos = 2000000)){
#then: 
#        while($obj->readPosition){
#            <do something with each line/variant (_currentLine/_currentVar)>
#        }
#    }
    my ($self) = @_;
    return if not $self->{_searchMatches};
    return if not @{$self->{_searchMatches}};
    $self->{_currentLine} = shift @{$self->{_searchMatches}};
    chomp ($self->{_currentLine});
    @{$self->{_split_line}} = split("\t", $self->{_currentLine});
    $self->readVariant;
    return $self->{_currentLine} if defined wantarray;
}

sub searchForPosition{
#my @matching lines = $obj->searchForPosition(chrom => 1, pos = 2000000)
#OR
#     if ($obj->searchForPosition(chrom => 1, pos = 2000000)){
#        while($obj->readPosition){
#            <do something with each matching line (_currentLine/_currentVar)>
#        }
#    }
    my ($self, %args) = @_;
    croak "Can't use searchForPosition method when noLineCount option is set and input is not bgzip compressed " if $self->{_noLineCount} and not $self->{_file} =~ /\.gz$/;
    croak "Can't use searchForPosition method when reading from STDIN" if $self->{_inputIsStdin} ;
    carp "chrom argument is required for searchForPosition method " if not exists $args{chrom};
    carp "pos argument is required for searchForPosition method " if not exists $args{pos};
    return if not exists $args{pos} or not exists $args{chrom};
    #if compressed assume bgzip compression and try to read with tabix
    if ($self->{_file} =~ /\.gz$/){
        eval "use Tabix; 1" 
            or croak "Tabix module is not installed and VCF file $self->{_file} appears to be bgzip compressed.  ".
            "  Please install Tabix.pm in order to search compressed VCFs.\n";
        if (not -e "$self->{_index}"){
            $self->setTabix if not $self->{_tabix}; 
            print STDERR "Indexing $self->{_file} with tabix.\n";
            my $er = `$self->{_tabix} -p vcf $self->{_file} 2>&1`;
            croak "Tabix indexing failed, code $?: $er\n" if $?;
        }
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
            binmode $self->{_indexHandle};
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


sub indexVcf{
    my ($self) = @_;
    if ($self->{_file} =~ /\.gz$/){
        $self->setTabix if not $self->{_tabix};
        print STDERR "Indexing $self->{_file} with tabix.\n";
        my $er = `$self->{_tabix} -p vcf $self->{_file} 2>&1`;
        croak "Tabix indexing failed, code $?: $er\n" if $?;
        return;
    }
    open (my $INDEX, "+>$self->{_index}") or croak "can't open $self->{_index} for writing index: $!";
    $self->reopenFileHandle();
    my $offset = 0;
    $self->{_contigIndex} = $self->{_file}.".pvdict";
    open (my $CONTIGS, ">$self->{_contigIndex}") or croak "Can't open $self->{_contigIndex} file for writing contig order: $! ";
    my %contigs = ();
    my $pack = "N";
    print STDERR "Indexing $self->{_file} ($self->{_totalLines} lines)...";
    $pack = "Q" if $self->{_fileSize} >= 4294967296;
    my $prev_pos = 0;
    while (my $line = scalar readline $self->{_filehandle}){
        print $INDEX pack($pack, $offset);
        $offset = tell $self->{_filehandle};
        next if $line =~ /^#/;
        my $chrom = (split "\t", $line)[$self->{_fields}->{CHROM}];
        my $pos = (split "\t", $line)[$self->{_fields}->{POS}];
        if (exists $contigs{$chrom}){
            if ($contigs{$chrom}  != scalar(keys%contigs) -1 or $pos < $prev_pos){
                close $INDEX;
                unlink $self->{_index} or carp "Couldn't delete partial index $self->{_index} - please delete manually: $!" ;
                croak "Can't index VCF - $self->{_file} not sorted properly ";
            }
        }else{
            $contigs{$chrom}  = scalar(keys%contigs);
            print $CONTIGS "$chrom\n";
        }
        $prev_pos = $pos;
    }
    print STDERR " Done indexing.\n";
    $self->{_contigOrder} = \%contigs;
    $self->reopenFileHandle();
    close $CONTIGS;
    close $INDEX;
    open ($INDEX, $self->{_index}) or croak "can't open $self->{_index} for reading index: $!";
    binmode $INDEX;
    return $INDEX if defined wantarray;
}

sub _lineWithIndex{
    my ($self, $line_number) =@_;
    my $size;               # size of an index entry
    my $i_offset;           # offset into the index of the entry
    my $entry;              # index entry
    my $d_offset;           # offset into the data file
    my $pack = "N";
    $pack = "Q" if $self->{_fileSize} >= 4294967296;
    $size = length(pack($pack, 0));
    $i_offset = $size * ($line_number-1);
    seek($self->{_indexHandle}, $i_offset, 0) or return;
    read($self->{_indexHandle}, $entry, $size);
    $d_offset = unpack($pack, $entry);
    seek($self->{_filehandle}, $d_offset, 0);
    return scalar readline $self->{_filehandle};
}
sub _readHeader{
    my ($self) = @_;
    my @meta_header = ();
    my @header = ();
    my $pos = tell $self->{_filehandle};
    while (my $vcf_line = scalar readline $self->{_filehandle}){
        if ($vcf_line =~ /^##/){
            push(@meta_header, $vcf_line);
        }elsif ($vcf_line =~ /^#/){
            push(@header, $vcf_line);
        }else{
            croak "No header found for VCF file $self->{_file} " if not @header;
            if ($self->{_inputIsStdin} or $self->{_file} =~ /\.gz$/){
                $self->{_firstLine} = $vcf_line;
                last;
            }else{
                #move back so that readLine method won't skip first line
                seek ($self->{_filehandle},0,0);
                seek ($self->{_filehandle}, $pos, 0);
                last;
            }
        }
        $pos = tell $self->{_filehandle};
    }
    croak "No header found for VCF file $self->{_file} " if not @header;
    if (@header > 1){
        carp "WARNING - Multiple headers found for VCF file $self->{_file} ";
        #in case there are multiple header lines we'll assume the last one is the REAL header
        #and put the others in the _metaHeader
        push (@meta_header, @header[0..$#header-1]);
    }
    $self->{_metaHeader} = \@meta_header;
    #in case there are multiple header lines we'll assume the last one is the REAL header
    chomp ($self->{_header} = $header[-1]);
    return if not defined wantarray;
    return (@meta_header, $self->{_header}) if wantarray;
    return $self->{_header};
}

sub _readHeaderInfo{
    my ($self) = @_;
    croak "Header not defined!" if not defined $self->{_header};#this should never happen - should have exited in _getHeader called by new
    my @split = split ("\t", $self->{_header});
    #Essential fields are CHROM POS ID REF ALT QUAL FILTER INFO
    foreach my $f (qw (CHROM POS ID REF ALT QUAL FILTER INFO)){
        my $i = 0;
        no warnings 'uninitialized';
        $i++ until $split[$i] =~ /^#*$f$/ or $i > $#split;
        if (not $self->{_noHeaderCheck}){
            croak "Can't find required $f field in header:\n$self->{_header}\n" if $i > $#split;
        }else{
            next if $i > $#split;
        }
        $self->{_fields} ->{$f} = $i;
        #$self->{_fields} is a ref to an anonymous hash giving the column of 
        #each of the essential fields in the header. While these fields are 
        #usually the same this allows for custom fields in VCFs e.g. before CHROM
    }
    #FORMAT FIELD ONLY PRESENT IF GENOTYPE INFORMATION IS AVAILABLE
    my $i = 0; 
    no warnings 'uninitialized';
    $i++ until $split[$i] =~ /#*FORMAT$/ or $i > $#split;
    $self->{_fields}->{FORMAT} = $i unless $i > $#split;
    #sample fields must be after FORMAT if genotpyes are present - assume anything after FORMATis a sample
    if (defined ($self->{_fields}->{FORMAT})){
        croak "No Samples found in header:\n$self->{_header}\n" if $self->{_fields}->{FORMAT} >= $#split;
        for (my $i = 1 + $self->{_fields}->{FORMAT}; $i < @split; $i++){
            $self->{_samples}->{$split[$i]} = $i;
            push @{$self->{_sampleOrder}}, $split[$i];
        }
        return ($self->{_fields}, $self->{_samples}) if defined wantarray;
    }else{#otherwise there should be no samples/genotypes
        return ($self->{_fields}) if defined wantarray;
    }
        
}

sub checkCoordinateSorted{
#return 1 if vcf is coordinate sorted (not checking chrom order)
    my ($self) = @_;
    $self->reopenFileHandle();
    my %contigs = ();
    my $prev_pos = 0;
    while (my $line = scalar readline $self->{_filehandle}){
        next if $line =~ /^#/;
        my $chrom = (split "\t", $line)[$self->{_fields}->{CHROM}];
        my $pos = (split "\t", $line)[$self->{_fields}->{POS}];
        if (exists $contigs{$chrom}){
            if ($contigs{$chrom}  != scalar(keys%contigs) -1 ){
                return 0; #not sorted - encountered contig twice with another inbetween
            }
            return 0 if $pos < $prev_pos;
        }else{
            $contigs{$chrom}  = scalar(keys%contigs);
        }
        $prev_pos = $pos;
    }
    $self->reopenFileHandle();
    return 1;
}

sub getInfoFields{
    #return a hash of INFO IDs, to anon hashes of 
    #Number, Type and Description
    my ($self) = @_;
    return %{$self->{INFO_FIELDS}} if exists $self->{INFO_FIELDS};
    my %info = ();
    foreach my $h (@{$self->{_metaHeader}}){
        if ($h =~ /^##INFO=<ID=(\w+),Number=([\.\w+]),Type=(\w+),Description=(.+)>$/){
            $info{$1} = 
                {
                Number => $2,
                Type => $3,
                Description => $4,
                };
        }
    }
    $self->{INFO_FIELDS} = \%info;
    return %info;
}


sub getHeaderColumns{
    my ($self) = @_;
    my %columns = ();
    foreach my $k (keys %{$self->{_fields}}){
        $columns{$k} = $self->{_fields}->{$k};
    }
    foreach my $k (keys %{$self->{_samples}}){
        $columns{$k} = $self->{_samples}->{$k};
    }
    return %columns;
}


sub getSampleNames{
    my ($self) =@_;
    croak "No samples found in header " if not defined $self->{_samples};
    return @{$self->{_sampleOrder}} if defined wantarray;
    #return keys %{$self->{_samples}} if defined wantarray;
    carp "getSampleNames method called in void context ";    
}

sub changeHeader{
    my ($self, %args) = @_;
    my @header = ();
    my @meta_header = ();
    if ($args{file}){
        open (my $HEAD, $args{file}) or croak "Can't open header file $args{file}: $! ";
        while (my $vcf_line = <$HEAD>){
            if ($vcf_line =~ /^##/){
                push(@meta_header, $vcf_line);
            }elsif ($vcf_line =~ /^#/){
                push(@header, $vcf_line);
            }else{
               if (not @header){
                    carp "No header found for VCF file $args{file} ";
                    close $HEAD;
                    return 0;
                }
                last;
            }
        }
    }elsif ($args{array}){
        croak "array argument passed to changeHeader method must be a reference to an arrray " if ref $args{array} ne 'ARRAY';
        foreach my $vcf_line (@{$args{array}}){
            if ($vcf_line =~ /^##/){
                            push(@meta_header, $vcf_line);
                    }elsif ($vcf_line =~ /^#/){
                            push(@header, $vcf_line);
                    }else{
                if (not @header){
                    carp "No header found for array ";
                    return 0;
                }
                last;
            }
        }
    }elsif($args{string}){
        my @lines = split("\n", $args{string});
        foreach my $vcf_line (@lines){
            if ($vcf_line =~ /^##/){
                push(@meta_header, $vcf_line);
            }elsif ($vcf_line =~ /^#/){
                push(@header, $vcf_line);
            }else{
                if (not @header){
                    carp "No header found for string ";
                    return 0;
                }
                last;
            }
        }
    }
    if (not @header){
        carp "No header found in changeHeader method " ;
        return 0;
    }
    if (@header > 1){
        carp "WARNING - Multiple headers found in changeHeader method";
        #in case there are multiple header lines we'll assume the last one is the REAL header
        #and put the others in the _metaHeader
        push (@meta_header, @header[0..$#header-1]);
    }
    if (@meta_header){
        $self->{_oldHeaderCount} += @{$self->{_metaHeader}} + 1;
        $self->{_metaHeader} = \@meta_header;
    }
    #in case there are multiple header lines we'll assume the last one is the REAL header
    chomp ($self->{_header} = $header[-1]);
    $self -> _readHeaderInfo();
    return 1;
}


sub getHeader{
    my ($self, $head_only) = @_;
    if (defined $head_only){
        if ($head_only){
            return $self->{_header} ."\n";
        }else{
            return join("", @{$self->{_metaHeader}});
        }
    }else{
        return join("", @{$self->{_metaHeader}}) . $self->{_header} ."\n";
    }
}

sub readLine{
    my ($self) = @_;
    if (($self->{_inputIsStdin} or $self->{_file} =~ /\.gz$/) and $self->{_firstLine}){
        $self->{_currentLine} = $self->{_firstLine};
        $self->{_firstLine} = '';
    }else{
        $self->{_currentLine} = scalar readline $self->{_filehandle};
    }
    if (defined $self->{_currentLine}){
        if ($self->{_currentLine} =~ /^#/){
#we should only encounter headers if we've specified a diff header on initialisation
#we need to keep track of these lines tho in order for countLines method to work
            no warnings 'recursion';#we might recurse through >100 headers
            $self->{_oldHeaderCount}++;
            $self->readLine;
        }else{
            $self->{_variantCount}++;
            chomp ($self->{_currentLine});
            @{$self->{_split_line}} = split("\t", $self->{_currentLine});
            $self->readVariant;
            return $self->{_currentLine} if defined wantarray;
        }
    }else{
        return;
    }
}

sub countLines{
    my ($self, $type) = @_;
    if ($self->{_noLineCount}){
        carp "Can't call countLines method when noLineCount argument is set ";
        return;
    }
    if ($self->{_inputIsStdin}){
        carp "Can't call countLines method when input is STDIN ";
        return;
    }
    $type = "variants" if not $type;
    my @valid_types = qw/total header remaining variants/;
    croak "invalid type $type passed to countLines method " if not grep {$type eq $_} @valid_types;
    my $header_count = 0;
    if ($self->{_oldHeaderCount}){#if we've replaced the header already
        $header_count = $self->{_oldHeaderCount};
    }else{
        $header_count = @{$self->{_metaHeader}} + 1; 
    }
    return $header_count if $type eq "header";
    
    return $self->{_totalLines} if $type eq "total";
    return ($self->{_totalLines} - $header_count - $self->{_variantCount}) if $type eq "remaining";
    return ($self->{_totalLines} - $header_count)  if $type eq "variants";
}
    
#this is not a perfect way of doing things but is a fudged attempt to
#split multiallelic records that may cause a problem with VEP
#makes no attempt at modifying INFO fields, has a go at modifying GT, AD and PL
#fields appropriately
sub splitMultiAllelicVariants{
    my ($self) = @_;
    $self->{_currentLine} ||= $self->readLine; #get line if we haven't already
    my $ref = $self->{_currentVar}->{REF};
    my @alleles = split(",", $self->{_currentVar}->{ALT});
    return  $self->{_currentLine} if (@alleles < 2);#not multiallelic
    my @splitLines = ();
    $self->getFormatFields();
    my %ac = $self->countAlleles();
    my $an = 0;
    map {$an += $ac{$_}} keys %ac;
    my @possible_codes = $self->getAllPossibleGenotypeCodes();
    for (my $i = 0; $i < @alleles; $i++){
        my $call_code = $i + 1;
        my @line = ();
        foreach my $field (qw(CHROM POS ID REF)){
            push @line, $self->getVariantField($field);
        }
        push @line, $alleles[$i];#add current ALT allele
        foreach my $field (qw(QUAL FILTER INFO FORMAT)){
            if ($field eq 'INFO'){
                my @new_inf;
                my $inf = $self->getVariantField($field);
                my @split_inf = split(";", $inf);
                foreach my $f (@split_inf){
                    if ($f =~ /^AC=/){
                        if (exists $ac{$i}){
                            push @new_inf, "AC=$ac{$i}";
                        }
                    }elsif ($f =~ /^AF=/){
                        if (exists $ac{$i} and $an > 0){
                            push @new_inf, "AF=" .($ac{$i}/$an);
                        }
                    }else{
                        push @new_inf, $f;
                    }
                }
                push @new_inf, "ParseVCF_split_variant-sample_fields_may_be_inaccurate" ;
                push @line, join(";", @new_inf);    
            }else{
                push @line, $self->getVariantField($field);
            }
        }
        foreach my $sample ($self->getSampleNames){
            my @sample_string = ();
            my $call = $self->getSampleVariant($sample);
            if ($call =~ /^\.[\/\|]\.$/){#no call
               push @line, $call;
               next;
            }else{
                my $delimiter = '/';
                foreach my $format (sort { 
                   $self->{_currentVar}->{varFormat}->{$a} <=> $self->{_currentVar}->{varFormat}->{$b} 
                    } keys %{$self->{_currentVar}->{varFormat}}){
                    if ($format eq 'GT'){#convert genotype to only include REF and $alleles[$i]
                        my $geno = $self->getSampleGenotypeField(sample => $sample, field => $format);
                        my @gt = split(/[\/\|]/, $geno);
                        foreach my $g (@gt){
                            if ($g ne '.' && $g != $call_code){#not a no call and not current allele call code
                                $g = 0;
                            }elsif($g == $call_code){
                                $g = 1;
                            }
                        }
                        if ($geno =~ /\d+(\||\/)\d+/){
                            $delimiter = $1;
                        }
                        push @sample_string, join($delimiter, @gt);
                    }elsif ($format eq 'AD'){
                        my $allele_depth = $self->getSampleGenotypeField(sample => $sample, field => $format);
                        if ($allele_depth eq '.'){
                            push @sample_string, $allele_depth;
                        }else{
                            my @ad = split(",", $allele_depth);
                            push @sample_string, "$ad[0],$ad[$call_code]";
                        }
                    }elsif ($format eq 'PL'){
                    #there will be a PL for each possible genotype.
                    #this cycles through each possible genotype and if it contains the current ALT allele we add 
                    #it to our new PL array
                        my @pl = split(",", $self->getSampleGenotypeField(sample => $sample, field => $format));
                        my @new_pl = ();
                        for (my $j = 0; $j < @possible_codes; $j++){
                            if ($possible_codes[$j] =~ /^(0|$call_code)[\/\|](0|$call_code)$/){
                                push @new_pl, $pl[$j];
                            }
                        }
                        push @sample_string, join(",", @new_pl);
                    }else{
                        push @sample_string, $self->getSampleGenotypeField(sample => $sample, field => $format);
                    }
                }
                push @line, join(":", @sample_string);
            }
        }
        push @splitLines, join("\t", @line);
    }
    return @splitLines;
}

sub getVariantID{
    my ($self) = @_;
    $self->{_currentLine} ||= $self->readLine; #get line if we haven't already
    return $self->{_currentVar}->{ID};
}

sub getVariantInfoField{
    my ($self, $info_field) = @_;
    if (not defined $self->{INFO_FIELDS}){
        $self->getInfoFields();
    }
    #if (not exists $self->{INFO_FIELDS}->{$info_field}){
       #return;
       # carp "INFO field $info_field does not exist in VCF header.\n";
    #}
    my @info = split(';', $self->getVariantField("INFO"));
    if (exists $self->{INFO_FIELDS}->{$info_field} && 
        $self->{INFO_FIELDS}->{$info_field}->{Type} eq 'Flag'){
        return 1 if grep {/^$info_field$/} @info;
    }else{
        #check type here?
        #check Number here?
        foreach my $inf (@info){
            if ($inf =~ /^$info_field=(\S+)/){
                return $1;
            }
        }
    }
    return;
}



sub getVariantField{
    my ($self, $field) = @_;
    $field =~ s/^#+//;
    $self->{_currentLine} ||= $self->readLine; #get line if we haven't already
    croak "Invalid field ($field) passed to getVariantField method " if not exists $self->{_validFields}->{$field} ;
    if (exists $self->{_currentVar}->{$field}){
        return $self->{_currentVar}->{$field};
    }else{
        return;
    }
}


sub getCustomField{
    my ($self, $field) = @_;
    $self->{_currentLine} ||= $self->readLine; #get line if we haven't already
    if (not exists $self->{_fields}->{$field}){
        my @split = split ("\t", $self->{_header});
        my $i = 0;
        no warnings 'uninitialized';
        $i++ until $split[$i] =~ /^#*$field$/ or $i > $#split;
        return if $i > $#split;
        $self->{_fields}->{$field} = $i;
    }
    return $self->{_split_line}->[$self->{_fields}->{$field}];
}

sub replaceVariantField{
    my ($self, $field, $replacement) = @_;
    $field =~ s/^#+//;
    $self->{_currentLine} ||= $self->readLine; #get line if we haven't already
    croak "Invalid field ($field) passed to replaceVariantField method " if not exists $self->{_validFields}->{$field} ;
    $self->{_currentVar}->{$field} = $replacement;
    my @line = ();
    foreach my $f ( split ("\t", $self->{_header})){
        $f =~ s/^#+//;
        if ($f eq $field){
            push @line, $replacement;
        }else{
            push @line, $self->getCustomField($f);
        }
    }
    @{$self->{_split_line}} = @line;
    $self->{_currentLine} = join("\t", @line);
    return $self->{_currentLine} if defined wantarray;
}

sub addVariantInfoField{
    my ($self, $field, $value) = @_;
    croak "a field must be passed to addVariantInfoField method " if not defined $field;
    $self->{_currentLine} ||= $self->readLine; #get line if we haven't already
    if ($field =~ /[,;]/){
        croak "INFO field ($field) passed to addVariantInfoField method has invalid characters.\n";
    }
    if (defined $value && $value =~ /;/){
        croak "value ($value) passed to addVariantInfoField method has invalid characters.\n";
    }
    my @info = split(";", $self->getVariantField('INFO'));
    my @new_inf = ();
    foreach my $inf (@info){
        if ($inf ne $field and $inf !~ /^$field=/){
            push @new_inf, $inf;
        }
    }
    if (defined $value){
        push @new_inf, "$field=$value";
    }else{
        push @new_inf, $field; 
    }
    return $self->replaceVariantField('INFO', join(";", @new_inf));
}

sub getColumnNumber{
    my ($self, $field) = @_;
    if (not exists $self->{_fields}->{$field}){
                my @split = split ("\t", $self->{_header});
                my $i = 0;
        no warnings 'uninitialized';
                $i++ until $split[$i] =~ /^#*$field$/ or $i > $#split;
                return if $i > $#split;
        return $i;
        }else{
            return $self->{_fields}->{$field};
    }
}


sub readVariant{
#e.g. while ($obj->readLine){
#         my %var = $obj->readVariant;
#     }
#
    my ($self, %args) = @_;
    $self->{_currentLine} ||= $self->readLine; #get line if we haven't already
    if (defined $self->{_currentLine}){
        my %var = map {$_ => $self->{_split_line}->[$self->{_fields}->{$_}]} keys %{$self->{_fields}};#e.g. $var{CHROM} = 1, $var{POS} = 10000
        if (defined $self->{_samples}){
            %{$var{calls}} = map{ $_ => $self->{_split_line}->[$self->{_samples}->{$_}] }  keys %{$self->{_samples}};#e.g. $var{calls}->{sample1} = '0/1:5,6:...'
        }
        $self->{_currentVar} = \%var;
        return %var if wantarray;
    }else{
        return;
    }
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
    my ($self) = @_;
    my %min_alleles = ();#key is allele number, each entry is anon hash of CHROM, REF, POS, ALT
    my @al =  $self->readAlleles();
    for (my $i = 1; $i < @al; $i++){
        my ($pos, $ref, $alt) = _reduceRefAlt($self->getVariantField("POS"), $al[0], $al[$i]);
        $min_alleles{$i} = {
            CHROM => $self->getVariantField("CHROM"),
            POS => $pos,
            REF => $ref,
            ALT => $alt,
            ORIGINAL_POS => $self->getVariantField("POS"),
            ORIGINAL_REF => $self->getVariantField("REF"),
            ORIGINAL_ALT => $self->getVariantField("ALT"),
        };
    }
    $self->{_minimizedAlleles} = \%min_alleles;
    return %min_alleles if wantarray;
}

sub _reduceRefAlt{
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

sub readAlleles{
    my ($self, %args) = @_;
    $self->{_currentLine} ||= $self->readLine; #get line if we haven't already
    if (defined $self->{_currentLine}){
        my @alleles = split(",", $self->{_currentVar}->{ALT});
        #if no variant at position then ALT will be "." and we want to remove it from our alleles array
        @alleles = grep {! /\./ } @alleles; 
        if (defined $args{alt_alleles} && $args{alt_alleles}){
            return @alleles if defined wantarray;
            carp "readAlleles method called in void context ";
        }
        unshift(@alleles, $self->{_currentVar}->{REF});
        #now ref is at index 0 of @alleles, and the alts correspond to the number in the call field
        return @alleles if defined wantarray;
        carp "readAlleles method called in void context ";
    }else{
        return;
    }
}

sub isMultiAllelic{
    my ($self, %args) = @_;
    my $alts = $self->readAlleles(alt_alleles => 1);
    return 1 if $alts > 1;
    return 0;
}

sub getFormatFields{
    my ($self) = @_;
    if (not defined ($self->{_fields}->{FORMAT})){
        carp "No FORMAT field in VCF "; 
        return;
    }
    $self->{_currentLine} ||= $self->readLine; #get line if we haven't already
    my %format_fields = ();
    my @format = split(":", $self->{_currentVar}->{FORMAT});
    my $i = 0;
    foreach my $f (@format){
        $format_fields{$f} = $i++;
    }
    $self->{_currentVar}->{varFormat} = \%format_fields;
    return %format_fields if wantarray;
    return keys %format_fields if defined wantarray;
}

sub checkSampleInVcf{
    my ($self, $sample) = @_;
    return 1 if exists $self->{_samples}->{$sample};
    return 0;
}

sub getSampleVariant{
    my ($self, $sample) = @_;
    croak "Can't invoke getSampleVariant method when no samples/genotypes are present in VCF " if not defined $self->{_samples};
    $self->{_currentLine} ||= $self->readLine; #get line if we haven't already
    if (defined $self->{_currentVar}){
        if ($sample){
            croak "Sample $sample does not exist - samples in header are: ". join(", ",  keys %{$self->{_samples}}) ." "  if not exists $self->{_samples}->{$sample};
            return $self->{_currentVar} -> {calls} -> {$sample} if defined wantarray;
        }else{#otherwise just return first sample
            return $self->{_split_line}->[$self->{_fields}->{FORMAT}+1];    
        }
    }else{ 
        return;
    }
}

sub getSampleCall{
#returns scalar value for genotype called 
#unless multiple argument is used in which
#case a hash of sample=>call key/values is returned
#returns './.' for samples below $args{minGQ}
#use return_alleles_only => 1 to only return allele codes, not genotypes (always returns an array)
    my ($self, %args) = @_;
    croak "Can't invoke getSampleCall method when no samples/genotypes are present in VCF " if not defined $self->{_samples};
    carp "WARNING Both multiple and sample arguments supplied to getSampleCall method - only multiple argument will be used " if (defined $args{multiple} and defined $args{sample});
    $self->{_currentLine} ||= $self->readLine; #get line if we haven't already
    if (defined $self->{_currentVar}){
        my $var; 
        if ($args{all}){
            my %calls = ();
            foreach my $sample (keys %{$self->{_samples}}){
                if (exists $args{minGQ} and $args{minGQ} > 0){
                    if (not defined $self->getSampleGenotypeField(sample=>$sample, field=>'GQ')){
                        #no GQ field - return no call (the above generally means that the call is './.' anyway)
                        $calls{$sample} = './.';
                        next;
                    }
                    if ($self->getSampleGenotypeField(sample=>$sample, field=>'GQ') < $args{minGQ}){
                        $calls{$sample} = './.';
                        next;
                    }
                }
                my $call = $self->getSampleGenotypeField(sample=>$sample, field=>'GT'); 
                $calls{$sample} = $call;
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
            my %calls = ();
            foreach my $sample (@{$args{multiple}}){
                if (exists $args{minGQ} and $args{minGQ} > 0){
                    if (not defined $self->getSampleGenotypeField(sample=>$sample, field=>'GQ')){
                    #no GQ field - return no call (the above generally means that the call is './.' anyway)
                        $calls{$sample} = './.';
                        next;
                    }
                    if ($self->getSampleGenotypeField(sample=>$sample, field=>'GQ') < $args{minGQ}){
                        $calls{$sample} = './.';
                        next;
                    }
                }
                my $call = $self->getSampleGenotypeField(sample=>$sample, field=>'GT'); 
                $calls{$sample} = $call;
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
            if (exists $args{minGQ} and $args{minGQ} > 0){
                if (not defined $self->getSampleGenotypeField(sample=>$args{sample}, field=>'GQ') ){
                #the above generally means that the call is './.' anyway
                    if ($args{return_alleles_only}){
                        return '.';
                    }else{
                        return './.';
                    }
                }
                if ($self->getSampleGenotypeField(sample=>$args{sample}, field=>'GQ') < $args{minGQ}){
                    if ($args{return_alleles_only}){
                        return '.';
                    }else{
                        return './.';
                    }
                }
            }
            my $call = $self->getSampleGenotypeField(sample=>$args{sample}, field=>'GT'); 
            if ($args{return_alleles_only}){
                my @al = split(/[\/\|]/, $call);
                return @al if defined wantarray;
            }else{
                return $call if defined wantarray;
            }
        }else{#otherwise just look at first sample
            if (exists $args{minGQ} and $args{minGQ} > 0){
                #the above generally means that the call is './.' anyway
                if ($args{return_alleles_only}){
                    return '.' if defined wantarray;
                }else{
                    return './.' if defined wantarray;
                }
            }    
            if ($self->getSampleGenotypeField( field=>'GQ') < $args{minGQ}){
                if ($args{return_alleles_only}){
                    return '.' if defined wantarray;
                }else{
                    return './.' if defined wantarray;
                }
            }
            my $call = $self->getSampleGenotypeField(field=>'GT');     
            if ($args{return_alleles_only}){
                my @al = split(/[\/\|]/, $call);
                return @al if defined wantarray;
            }else{
                return $call if defined wantarray;
            }
        }
    }else{
        return if defined wantarray;
    }
    carp "getSampleCall called in a void context ";
}

sub getSampleGenotypeField{
#returns scalar value for genotype field 
#unless multiple argument is used in which
#case a hash of sample=>value key/values is returned
    my ($self, %args) = @_;
    croak "Can't invoke getSampleGenotypeField method when no samples/genotypes are present in VCF " if not defined $self->{_samples};
    carp "WARNING Both multiple and sample arguments supplied to getSampleGenotypeField method - only multiple argument will be used " if (defined $args{multiple} and defined $args{sample});
    croak "\"field\" argument must be passed to getSampleGenotypeField - e.g. getSampleGenotypeField(field=>\"GQ\") " if not defined $args{field};
    $self->getFormatFields();
    if (not defined $self->{_currentVar}->{varFormat}->{$args{field}}){
        carp "Field $args{field} not found for getSampleGenotypeField ";
        return;
    }
    $self->{_currentLine} ||= $self->readLine; #get line if we haven't already
    if (defined $self->{_currentVar}){
        my $var; 
        if ($args{all}){
            my %values = ();
            foreach my $sample (keys %{$self->{_samples}}){
                my $mvar = $self->getSampleVariant($sample);
                my $value = (split ":", $mvar)[$self->{_currentVar}->{varFormat}->{$args{field}}];
                $values{$sample} =  $value;
            }
                        return %values if defined wantarray;
                        carp "getSampleGenotypeField called in a void context ";
        }elsif($args{multiple}){
            croak "multiple argument must be an array reference " if ref $args{multiple} ne 'ARRAY';
            my %values = ();
            foreach my $sample (@{$args{multiple}}){
                my $mvar = $self->getSampleVariant($sample);
                my $value = (split ":", $mvar)[$self->{_currentVar}->{varFormat}->{$args{field}}];
                $values{$sample} =  $value;
            }
            return %values if defined wantarray;
            carp "getSampleGenotypeField called in a void context ";
        }elsif (defined $args{sample}){
            $var = $self->getSampleVariant($args{sample});
        }else{#otherwise just look at first sample
            $var = $self->getSampleVariant();
        }
        return if not defined $var and defined wantarray;
        my $value = (split ":", $var)[$self->{_currentVar}->{varFormat}->{$args{field}}];
        return $value if defined wantarray;
    }else{
        return if defined wantarray;
    }
    carp "getSampleGenotypeField called in a void context ";
}

#returns genotype codes - e.g. 0/0, 0/1, 1/1
sub getAllPossibleGenotypeCodes{
    my ($self, %args) = @_;
    my @alleles = $self -> readAlleles;
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
    my ($self, %args) = @_;
    my @alleles = $self -> readAlleles;
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
    carp "WARNING Both multiple and sample arguments supplied to getSampleActualGenotypes method - only multiple argument will be used " if (defined $args{multiple} and defined $args{sample});
    $self->{_currentLine} ||= $self->readLine; #get line if we haven't already
    my @alleles = $self -> readAlleles;
    my %multiple = ();
    my $genotype;
    if (defined $self->{_currentVar}){
        my @sample_alleles = ();
        my $var;
        if ($args{all}){
            foreach my $sample (keys %{$self->{_samples}}){
                $var = $self->getSampleVariant($sample);
                my $call = (split ":", $var)[0];
                if ($call =~ /(\d+)([\/\|])(\d+)/){
                    if ($args{return_alleles_only}){
                        push (@sample_alleles, ($alleles[$1], $alleles[$3]));
                    }else{
                        $multiple{$sample} = "$alleles[$1]$2$alleles[$3]";
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
                $var = $self->getSampleVariant($sample);
                my $call = (split ":", $var)[0];
                if ($call =~ /(\d+)([\/\|])(\d+)/){
                    if ($args{return_alleles_only}){
                        push (@sample_alleles, ($alleles[$1], $alleles[$3]));
                    }else{
                        $multiple{$sample} = "$alleles[$1]$2$alleles[$3]";
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
            $var = $self->getSampleVariant($args{sample});
        }else{
            $var = $self->getSampleVariant();
        }
        my $call = (split ":", $var)[0];
        if ($call =~ /(\d+)([\/\|])(\d+)/){
            if ($args{return_alleles_only}){
                push (@sample_alleles, ($alleles[$1], $alleles[$3]));
                my %seen = ();
                @sample_alleles = grep {!$seen{$_}++} @sample_alleles;#remove duplicates
                return @sample_alleles;
                
            }else{
                $genotype = "$alleles[$1]$2$alleles[$3]";
            }
        }else{
            $genotype = "-/-";
        }
        return if $args{return_alleles_only};
        return $genotype 
    }
}

sub sampleIsHomozygous{
# use minGQ => 30 to count any sample with less than GQ = 30 as a no call
# use variants => 1 to only include homozygous variants, not homozygous variants
#returns 1 if all samples are homozygous or a mixture of homozygous and no calls
#returns 0 if any sample is homozygous or all are no calls - we should improve this for multiple samples.
    my ($self, %args) = @_;
    croak "Can't invoke sampleIsHomozygous method when no samples/genotypes are present in VCF " if not defined $self->{_samples};
    $self->{_currentLine} ||= $self->readLine; #get line if we haven't already
    carp "WARNING Both multiple and sample arguments supplied to sampleIsHomozygous method - only multiple argument will be used " if (defined $args{multiple} and defined $args{sample});
    if (defined $self->{_currentVar}){
        my $var;
        if ($args{all}){
            my $no_call = 0;
            foreach my $sample (keys %{$self->{_samples}}){
                if ($args{minGQ}){
                    my $gq = $self->getSampleGenotypeField(sample=>$sample, field=>'GQ');
                    if (defined $gq && $gq < $args{minGQ}){
                        $no_call++;
                        next;
                    }
                }
                $var = $self->getSampleVariant($sample);
                my $call = (split ":", $var)[0];
                if ($call =~ /(\d+)[\/\|](\d+)/){
                    return 0 if $1 != $2;
                    if ($args{variants}){
                        return 0 if ($1 == 0 and $2 == 0);
                    }
                }else{
                    $no_call++;
                }
            }
            return 0 if $no_call == keys %{$self->{_samples}};#i.e. return 0 if all our samples are no calls
            return 1;
        }elsif(defined $args{multiple}){
            croak "multiple argument must be an array reference " if ref $args{multiple} ne 'ARRAY';
            my $no_call = 0;
            foreach my $sample (@{$args{multiple}}){
                if ($args{minGQ}){
                    my $gq = $self->getSampleGenotypeField(sample=>$sample, field=>'GQ');
                    if (defined $gq && $gq < $args{minGQ}){
                        $no_call++;
                        next;
                    }
                }
                $var = $self->getSampleVariant($sample);
                my $call = (split ":", $var)[0];
                if ($call =~ /(\d+)[\/\|](\d+)/){
                    return 0 if $1 != $2;
                    if ($args{variants}){
                        return 0 if ($1 == 0 and $2 == 0);
                    }
                }else{
                    $no_call++;
                }
            }
            return 0 if $no_call == @{$args{multiple}};#i.e. return 0 if all our samples are no calls
            return 1;
        }elsif (defined $args{sample}){
            if ($args{minGQ}){
                my $gq = $self->getSampleGenotypeField(sample=>$args{sample}, field=>'GQ');
                if (defined $gq && $gq < $args{minGQ}){
                    return 0;
                }
            }
            $var = $self->getSampleVariant($args{sample});
        }else{ 
            if ($args{minGQ}){
                my $gq = $self->getSampleGenotypeField(field=>'GQ');
                if (defined $gq && $gq < $args{minGQ}){
                    return 0;
                }
            }
            $var = $self->getSampleVariant();
        }
        my $call = (split ":", $var)[0];
        if ($call =~ /(\d+)[\/\|](\d+)/){
            if ($args{variants}){
                return 0 if ($1 == 0 and $2 == 0);
            }
            return 1 if $1 == $2;
            return 0;
        }
        return 0;
    }
}

sub sampleIsHeterozygous{
# use minGQ => 30 to count any sample with less than GQ = 30 as a no call
#returns 1 if all samples are heterozygous or a mix of hets and no calls
#returns 0 if any sample is homozygous or all are no calls - we should improve this for multiple samples.
    my ($self, %args) = @_;
    croak "Can't invoke sampleIsHeterozygous method when no samples/genotypes are present in VCF " if not defined $self->{_samples};
    $self->{_currentLine} ||= $self->readLine; #get line if we haven't already
    carp "WARNING Both multiple and sample arguments supplied to sampleIsHomozygous method - only multiple argument will be used " if (defined $args{multiple} and defined $args{sample});
    if (defined $self->{_currentVar}){
        my $var;
        if ($args{all}){
            my $no_call = 0;
            foreach my $sample (keys %{$self->{_samples}}){
                if ($args{minGQ}){
                    my $gq = $self->getSampleGenotypeField(sample=>$sample, field=>'GQ');
                    if (defined $gq && $gq < $args{minGQ}){
                        $no_call++;
                        next;
                    }
                }
                $var = $self->getSampleVariant($sample);
                my $call = (split ":", $var)[0];
                if ($call =~ /(\d+)[\/\|](\d+)/){
                    return 0 if $1 == $2;
                }else{
                    $no_call++;
                }
            }
            return 0 if $no_call == keys%{$self->{_samples}};#i.e. return 0 if all our samples are no calls
            return 1;#otherwise return 1
        }elsif(defined $args{multiple}){
            my $no_call = 0;
            croak "multiple argument must be an array reference " if ref $args{multiple} ne 'ARRAY';
            foreach my $sample (@{$args{multiple}}){
                if ($args{minGQ}){
                    my $gq = $self->getSampleGenotypeField(sample=>$sample, field=>'GQ');
                    if (defined $gq && $gq < $args{minGQ}){
                        $no_call++;
                        next;
                    }
                }
                $var = $self->getSampleVariant($sample);
                my $call = (split ":", $var)[0];
                if ($call =~ /(\d+)[\/\|](\d+)/){
                    return 0 if $1 == $2;
                }else{
                    $no_call++;
                }
            }
            return 0 if $no_call == @{$args{multiple}};#i.e. return 0 if all our samples are no calls
            return 1;#otherwise return 1
        }elsif (defined $args{sample}){
            if ($args{minGQ}){
                my $gq = $self->getSampleGenotypeField(sample=>$args{sample}, field=>'GQ');
                if (defined $gq && $gq < $args{minGQ}){
                    return 0;
                }
            }
            $var = $self->getSampleVariant($args{sample});
        }else{    
            if ($args{minGQ}){
                my $gq = $self->getSampleGenotypeField(field=>'GQ');
                if (defined $gq && $gq < $args{minGQ}){
                    return 0;
                }
            } 
            $var = $self->getSampleVariant();
        }
        my $call = (split ":", $var)[0];
        if ($call =~ /(\d+)[\/\|](\d+)/){
            return 1 if $1 != $2;
            return 0;
        }
        return 0;
    }
}
sub countAlleles{
    #%allele_counts = $obj->countGenotypes(); 
    #returns count for all samples and all genotypes
    #%allele_counts = $obj->countGenotypes(samples=>["sample1","sample2"])
    #$count = $obj->countGenotypes(samples=>["sample1","sample2"], genotypes => '0/1');
    my ($self, %args) = @_;
    if ($args{minGQ}){
        croak "Can't invoke countAlleles method with minGQ argument when no samples/genotypes are present in VCF " if not defined $self->{_samples};
        croak "minGQ argument must be a number " if not looks_like_number($args{minGQ});
    }
    my %counts; 
    $self->{_currentLine} ||= $self->readLine; #get line if we haven't already
    defined $self->{INFO_FIELDS} || $self->getInfoFields();
    if (not defined $args{samples}){
        if (not $args{minGQ} and exists $self->{INFO_FIELDS}->{AC} 
          and exists $self->{INFO_FIELDS}->{AN}){
            my @ac = split(",", $self->getVariantInfoField('AC'));
            unshift @ac,  $self->getVariantInfoField('AN') - sum(@ac);
            for (my $i = 0; $i < @ac; $i++){
                $counts{$i} = $ac[$i];
            }
            return %counts;
        }
    }
    my @alleles = $self->readAlleles();
    for (my $i = 0; $i < @alleles; $i++){
        $counts{$i} = 0;
    }
    if(defined $args{samples}){
        croak "samples argument must be an array reference " if ref $args{samples} ne 'ARRAY';
        foreach my $sample (@{$args{samples}}){
            my $call;
            if (defined $args{minGQ}){
                $call = $self->getSampleCall(sample=>$sample, minGQ => $args{minGQ});
            }else{
                $call= $self->getSampleCall(sample=>$sample);
            }
            my @ca = split(/[\/\|]/, $call);
            @ca = grep {$_ ne '.'} @ca;
            map {$counts{$_}++} @ca;
        }
    }else{
        foreach my $sample (keys %{$self->{_samples}}){
            my $call;
            if (defined $args{minGQ}){
                $call = $self->getSampleCall(sample=>$sample, minGQ => $args{minGQ});
            }else{
                $call= $self->getSampleCall(sample=>$sample);
            }
            my @ca = split(/[\/\|]/, $call);
            @ca = grep {$_ ne '.'} @ca;
            map {$counts{$_}++} @ca;
        }
    }
    return %counts;
}
sub countGenotypes{
    #%genotypes = $obj->countGenotypes(); 
    #returns count for all samples and all genotypes
    #%genotypes = $obj->countGenotypes(genotypes => ['0/1', '0/2']);
    #%genotypes = $obj->countGenotypes(samples=>["sample1","sample2"], genotypes => ['0/1', '0/2']);
    #$count = $obj->countGenotypes(samples=>["sample1","sample2"], genotypes => '0/1');
    my ($self, %args) = @_;
    croak "Can't invoke countGenotypes method when no samples/genotypes are present in VCF " if not defined $self->{_samples};
    if (defined $args{minGQ}){
        croak "minGQ argument must be a number " if not looks_like_number($args{minGQ});
    }
    $self->{_currentLine} ||= $self->readLine; #get line if we haven't already
    my %genotypes = ();
    if(defined $args{samples}){
        croak "samples argument must be an array reference " if ref $args{samples} ne 'ARRAY';
        foreach my $sample (@{$args{samples}}){
            if (defined $args{minGQ}){
                my $call = $self->getSampleCall(sample=>$sample, minGQ => $args{minGQ});
                $genotypes{$call}++;
            }else{
                my $call = $self->getSampleCall(sample=>$sample);
                $genotypes{$call}++;
            }
        }
    }else{
        foreach my $sample (keys %{$self->{_samples}}){
            if (defined $args{minGQ}){
                my $call = $self->getSampleCall(sample=>$sample, $args{minGQ});
                $genotypes{$call}++;
            }else{
                my $call = $self->getSampleCall(sample=>$sample);
                $genotypes{$call}++;
            }
        }
    }
    if (defined $args{genotypes}){
        my %user_gts = ();
        if (ref $args{genotypes} eq 'ARRAY'){
            foreach my $gt (@{$args{genotypes}}){
                if (exists $genotypes{$gt}){
                    $user_gts{$gt} =  $genotypes{$gt};
                }else{
                    $user_gts{$gt} = 0;
                }
            }
            return %user_gts;
        }else{
            if (exists $genotypes{$args{genotypes}}){
                return $genotypes{$args{genotypes}};
            }else{
                return 0;
            }
        }
    }
        
    return %genotypes;
}

sub countZygosity{
# use minGQ => 30 to count any sample with less than GQ = 30 as a no call
# use samples => [sam1, sam2, sam3,] to only count specific samples
# use return => het to count no. of het calls
# use return => hom to return no. of hom calls 
# use return => nocall to return no. of no calls 
# use return => both to return 2 element array of het calls and hom calls
# use return => all to return 3 element array of het calls, hom calls and no calls
    my ($self, %args) = @_;
    croak "Can't invoke countZygosity method when no samples/genotypes are present in VCF " if not defined $self->{_samples};
    if (not defined $args{return}){
        $args{return} = "all";
    }else{
        my @valid_types = qw(het hom nocall both all);
        croak "argument return value $args{return} not recognised in countZygosity method " if not grep { /$args{return}/ } @valid_types;
    }
    $self->{_currentLine} ||= $self->readLine; #get line if we haven't already
    my $hets = 0;
    my $homs = 0;
    my $no_calls = 0;
    if (defined $self->{_currentVar}){
        if(defined $args{samples}){
            croak "samples argument must be an array reference " if ref $args{samples} ne 'ARRAY';
            foreach my  $sample (@{$args{samples}}){
                my $var = $self->getSampleVariant($sample);
                if ($args{minGQ}){
                    my $gq = $self->getSampleGenotypeField(sample=>$sample, field=>'GQ');
                    if (defined $gq && $gq < $args{minGQ}){
                        $no_calls++;
                        next;
                    }
                } 
                my $call = (split ":", $var)[0];
                    if ($call =~ /(\d+)[\/\|](\d+)/){
                        $hets++ if $1 != $2;
                        $homs++ if $1 == $2;
                    }else{
                        $no_calls++;
                    }
                }
        }else{
            foreach my $sample (keys %{$self->{_samples}}){
                my $var = $self->getSampleVariant($sample);
                if ($args{minGQ}){
                    my $gq = $self->getSampleGenotypeField(sample=>$sample, field=>'GQ');
                    if (defined $gq && $gq < $args{minGQ}){
                        $no_calls++;
                        next;
                    }
                } 
                my $call = (split ":", $var)[0];
                if ($call =~ /(\d+)[\/\|](\d+)/){
                    $hets++ if $1 != $2;
                    $homs++ if $1 == $2;
                }else{
                    $no_calls++;
                }
            }
                
        }
    }
    if ($args{return} eq "all"){
        return ($hets, $homs, $no_calls);
        carp "countZygosity method called in void context ";
    }elsif ($args{return} eq "het"){
        return $hets;
        carp "countZygosity method called in void context ";
    }elsif ($args{return} eq "hom"){
        return $homs;
        carp "countZygosity method called in void context ";
    }elsif ($args{return} eq "nocall"){
        return $no_calls;
        carp "countZygosity method called in void context ";
    }elsif ($args{return} eq "both"){
        return ($hets, $homs);
        carp "countZygosity method called in void context ";
    }
}

1;

=head1 NAME

ParseVCF.pm - read standard and custom Variant Call Format (VCF) files and return information.

=head1 SYNOPSIS


DOCUMENTATION NEEDS SERIOUS ATTENTION
#TO DO - ADD filtering on GQ etc. for sample variants, test changeHeader method
#NEED DOCUMENTATION FOR getSampleGenotypeField, getSampleCall, getFormatFields, getAltAlleles, reopenFileHandle, searchForPosition, readPosition, getVepFields, checkSampleInVcf methods

$obj = ParseVCF -> new(file => $vcf);
#(initialise object with a VCF file)

my $meta_header_string = $obj -> getHeader;
#(get full header including meta information)

my $header_string = $obj ->  getHeader( 1 ) ;
#(get header line only)

my $line_count = $obj ->  countLines;
#(get length of VCF file)

while (my $var_line = $obj ->  readLine){
    #(do something with each variant)
}

my $rs_id = $obj ->  getVariantField('ID');
#(get value for a standard field in current variant)

my $custom = $obj ->  getCustomField('AdditionalColumn');
#(get value for a custom field in current variant)

my %var = $obj ->  readVariant;
#(get a hash of fields and values for current variant)

my @alleles = $obj ->  readAlleles;
#(get an array of all alleles present at variant site)

my $variant_field = $obj ->  getSampleVariant('SAMPLE1');
#(returns the variant field with genotype information etc.)

my @possible_genotypes = $obj ->  getAllPossibleGenotypes;
#(returns an array of all possible allele combinations at variant site)

my $genotpye = $obj ->  getSampleActualGenotypes(sample => 'SAMPLE1');
#(returns a string containing a samples called genotype - e.g. "A/G")

my %genotypes = $obj ->  getSampleActualGenotypes(multiple => ['SAMPLE1', 'SAMPLE2', 'SAMPLE3']);
#(returns a hash of each sample and its genotype - e.g. SAMPLE1 => "A/G")

my @observed_alleles = $obj ->  getSampleActualGenotypes(return_alleles_only => 1, multiple => ['SAMPLE1', 'SAMPLE2', 'SAMPLE3']);
#(returns an array of alleles found in specified samples)

my $is_hom = $obj ->  sampleIsHomozygous(sample => 'SAMPLE1');
#(returns 1 if sample is homozygous, 0 if not)

my $is_hom = $obj ->  sampleIsHomozygous(multiple => ['SAMPLE1', 'SAMPLE2', 'SAMPLE3']);
#(returns 1 all samples are homozygous, 0 if not)

my $is_het = $obj ->  sampleIsHeterozygous(sample => 'SAMPLE1');
#(returns 1 if sample is heterozygous, 0 if not)

my $is_het = $obj ->  sampleIsHeterozygous(multiple => ['SAMPLE1', 'SAMPLE2', 'SAMPLE3']);
#(returns 1 all samples are heterozygous, 0 if not)

my $het_calls = $obj ->  countZygosity(return => 'het');    
#(returns no. of het calls)

my $hom_calls = $obj ->  countZygosity(return => 'hom');    
#(returns no. of hom calls)

my ($hets_calls, $hom_calls)  = $obj ->  countZygosity(return => 'both');    
#(returns 2 element array of het calls and hom calls)

my $no_calls = $obj ->  countZygosity(return => 'nocall');    
#(returns no. of no calls)

my ($hets_calls, $hom_calls, $no_calls) = $obj ->  countZygosity(return => 'all');    
#(returns 3 element array of het calls, hom calls and no calls)

my $variants_processed = $obj -> get_variantCount;
#(returns number of variants read so far, including current variant)

my $success = $obj -> changeHeader(string=>"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tNEW_SAMPLE_NAME");
#(change header info for VCF.  Won't change the original file but can be printed to new file.  Header info will be read from the new header if successful. Meta header information will be replaced if new meta header lines are present)

my $success = $obj -> changeHeader(array=>\@new_header_lines);
#(as above but new header is read from array reference)

my $success = $obj -> changeHeader(file=>"new_header.txt");
#(as above but new header is read from file)


=head1 DESCRIPTION

=head2 Overview

This module can be used to read Variant Call Format (VCF) files and return useful information for variant processing. 

THIS MODULE AND DOCUMENTATION IS STILL UNDER CONSTRUCTION

=head2 Constructor and initialization

Minimal requirement for initialisation is to call the 'new' method specifying a VCF file to read.

$obj = ParseVCF -> new(file => $vcf);

The length of the file will be determined and header information read as part of initialisation.

If header is incorrect or missing you must correct this upon initialisation.

$obj = ParseVCF -> new(file => $vcf, header => $header);


=head2 Class and object methods

DOCUMENTATION STILL UNDER CONSTRUCTION

=head3 Accessors and mutators

The following features can be accessed using get_[feature] or set using set_[feature], substituting [feature] for the feature of interest. 

=over 12

=item B<file>

filename of VCF.

=item B<currentLine>

last read line of vcf. Read only.

=item B<variantCount>

Number of variants processed so far. Read only.

=item B<totalLines>

Total number of lines in VCF file. Read only.

=item B<noHeaderCheck>

If set to true then headers won't be checked for required fields. 

=back

=head3 Methods

=over 12

=back


=head1 AUTHOR

David A. Parry
University of Leeds


=head1 COPYRIGHT AND LICENSE

Copyright 2011, 2012  David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut
