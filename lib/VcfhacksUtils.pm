=head1 NAME

VcfhacksUtils.pm - shared convenience methods for vcfhacks scripts

=head1 VERSION

version 0.1

=head1 SYNOPSIS

 use VcfhacksUtils;
 
 ...


=cut

package VcfhacksUtils;
use strict;
use warnings;
use Carp;
use File::Basename; 
use POSIX qw/strftime/;
use FindBin qw($RealBin);

my $data_dir = "$RealBin/data";

=head1 FUNCTIONS


=head2 General Output Utilities

=over 12

=item B<getOptsVcfHeader>

Takes an options hash and returns a header line documenting the options passed to the program.

 my $optstring = VcfhacksUtils::getOptsVcfHeader(%opts); 

=cut

sub getOptsVcfHeader{
    my %opts = @_;
    my @opt_string = ();
    #get options into an array
    foreach my $k ( sort keys %opts ) {
        if ( not ref $opts{$k} ) {
            push @opt_string, "$k=$opts{$k}";
        }elsif ( ref $opts{$k} eq 'SCALAR' ) {
            if ( defined ${ $opts{$k} } ) {
                push @opt_string, "$k=${$opts{$k}}";
            }else {
                push @opt_string, "$k=undef";
            }
        }elsif ( ref $opts{$k} eq 'ARRAY' ) {
            if ( @{ $opts{$k} } ) {
                push @opt_string, "$k=" . join( ",", @{ $opts{$k} } );
            }else {
                push @opt_string, "$k=undef";
            }
        }
    }
    my $caller = fileparse($0);
    return "##$caller\"" . join( " ", @opt_string ) . "\"";
    
}


=item B<getInfoHeader>

Takes a hash containing the keys 'ID, Number, Type and Description' and outputs an appropriately formatted INFO header.

 my %info_field = 
 (
     ID          => "AnInfoField",
     Number      => "A",
     Type        => "String",
     Description => "A made up VCF INFO field.",
 );
 my $inf_string = VcfhacksUtils::getInfoHeader(%info_field); 
 
=cut

sub getInfoHeader{
    my %info = @_;
    foreach my $f ( qw / ID Number Type Description / ){
        if (not exists $info{$f}){
            carp "INFO field '$f' must be provided to getInfoHeader method!\n";
            return;
        }
    }
    return "##INFO=<ID=$info{ID},Number=$info{Number},".
      "Type=$info{Type},Description=\"$info{Description}\">";
}

=item B<getFormatHeader>

Takes a hash containing the keys 'ID, Number, Type and Description' and outputs an appropriately formatted FORMAT header.

 my %format_field = 
 (
     ID          => "AFormatField",
     Number      => "A",
     Type        => "String",
     Description => "A made up VCF INFO field.",
 );
 my $f_string = VcfhacksUtils::getFormatHeader(%format_field); 
 
=cut

sub getFormatHeader{
    my %info = @_;
    foreach my $f ( qw / ID Number Type Description / ){
        if (not exists $info{$f}){
            carp "INFO field '$f' must be provided to getFormatHeader method!\n";
            return;
        }
    }
    return "##FORMAT=<ID=$info{ID},Number=$info{Number},".
      "Type=$info{Type},Description=\"$info{Description}\">";
}

=item B<getFilterHeader>

Takes a hash containing the keys 'ID, and DESCRIPTION' and outputs an appropriately formatted FORMAT header.

 my %filter_field = 
 (
     ID          => "AFilterField",
     Description => "A made up VCF FILTER field.",
 );
 my $f_string = VcfhacksUtils::getFilterHeader(%filter_field); 
 
=cut

sub getFilterHeader{
    my %info = @_;
    foreach my $f ( qw / ID Description / ){
        if (not exists $info{$f}){
            carp "INFO field '$f' must be provided to getFilterHeader method!\n";
            return;
        }
    }
    return "##FILTER=<ID=$info{ID},Description=\"$info{Description}\">";
}


=back

=head2 VEP/SnpEff Data Utilities

=over 12

=item B<readVepClassesFile>

Return a hash of VEP classes with values being either 'default' or 'valid'. 

 my %all_classes = VcfhacksUtils::readVepClassesFile();

=cut

sub readVepClassesFile{
    my $classes_file = "$data_dir/vep_classes.tsv";
    return _readClassesFile($classes_file);
}

=item B<readSnpEffClassesFile>

Return a hash of SnpEff classes with values being either 'default' or 'valid'. 

 my %all_classes = VcfhacksUtils::readSnpEffClassesFile();

=cut

sub readSnpEffClassesFile{
    my $classes_file = "$data_dir/snpeff_classes.tsv";
    return _readClassesFile($classes_file);
}


sub _readClassesFile{
    my $classes_file = shift;
    open (my $CLASS, $classes_file) or die 
"Could not open effect classes file '$classes_file': $!\n";
    my %classes = (); 
    while (my $line = <$CLASS>){
        next if $line =~ /^#/;
        $line =~ s/[\r\n]//g; 
        next if not $line;
        my @split = split("\t", $line);
        die "Not enough fields in classes file line: $line\n" if @split < 2;
        $classes{lc($split[0])} = lc($split[1]);
    }
    close $CLASS;
    return %classes;
}

=item B<readBiotypesFile>

Return a hash of biotypes classes with values being either 'filter' or 'keep' indicating default behaviour of functional consequence filtering scripts. 

 my %all_biotypes = VcfhacksUtils::readBiotypesFile();

=cut

sub readBiotypesFile{
    my $btypes_file = "$data_dir/biotypes.tsv";
    open (my $BT, $btypes_file) or die 
"Could not open biotypes file '$btypes_file': $!\n";
    my %biotypes = (); 
    while (my $line = <$BT>){
        next if $line =~ /^#/;
        $line =~ s/[\r\n]//g; 
        next if not $line;
        my @split = split("\t", $line);
        die "Not enough fields in biotypes file line: $line\n" if @split < 2;
        $biotypes{lc($split[0])} = lc($split[1]);
    }
    close $BT;
    return %biotypes;
}

=item B<readVepInSilicoFile>

Return a hash of in silico prediction program names with values being anonymous hashes of prediction scores to either 'default' or 'valid' indicating default behaviour of functional consequence filtering scripts. 

 my %in_silico = VcfhacksUtils::readVepInSilicoFile();

=cut

sub readVepInSilicoFile{
    my $insilico_file = "$data_dir/vep_insilico_pred.tsv";
    return _readInSilicoFile($insilico_file, 1);
}


=item B<readSnpEffInSilicoFile>

Return a hash of dbNSFP in silico prediction program names with values being anonymous hashes of prediction scores to either 'default' or 'valid' indicating default behaviour of functional consequence filtering scripts. 

 my %in_silico = VcfhacksUtils::readSnpEffInSilicoFile();

=cut

sub readSnpEffInSilicoFile{
    my $insilico_file = "$data_dir/snpeff_insilico_pred.tsv";
    return _readInSilicoFile($insilico_file, 0);
}

sub _readInSilicoFile{
    my $insilico_file = shift;
    my $convert_case = shift;
    open (my $INS, $insilico_file) or die 
"Could not open in silico classes file '$insilico_file': $!\n";
    my %pred = (); 
    while (my $line = <$INS>){
        next if $line =~ /^#/;
        $line =~ s/[\r\n]//g; 
        next if not $line;
        my @split = split("\t", $line);
        die "Not enough fields in classes file line: $line\n" if @split < 3;
        if ($convert_case){
            $pred{lc($split[0])}->{lc($split[1])} = lc($split[2]) ;
        }else{
            $pred{$split[0]}->{$split[1]} = $split[2] ;
        }
    }
    close $INS;
    return %pred;
}


=item B<getAndCheckInSilicoPred>

Takes the mode (either 'vep' or 'snpeff') as its first argument and a reference to an array of consequence filters (e.g. "polyphen=probably_damaging" or "all") as its second argument. Returns a hash of prediction programs and values to filter on.

 my %filters = VcfhacksUtils::getAndCheckInSilicoPred('vep', \@damaging); 

=cut

sub getAndCheckInSilicoPred{
    my $mode = shift;
    my $damaging = shift;
    if (ref $damaging ne 'ARRAY'){
        croak "Second argument to 'getAndCheckInSilicoPred' method must be an".
            " array reference!\n";
    }
    my $conv_case = $mode eq 'vep' ? 1 : 0;
    my %pred = _readInSilicoFile("$data_dir/$mode"."_insilico_pred.tsv", $conv_case);
    my %filters = (); 
    foreach my $d (@$damaging) {
        if ($conv_case){
            $d = lc($d); 
        }
        my ( $prog, $label ) = split( "=", $d );
        if ($mode eq 'snpeff' and $prog !~ /^all$/i){
            $prog = "dbNSFP_$prog" if $prog !~ /^dbNSFP_/;
            $prog =~ s/\-/_/g; 
        }
        if ($prog =~ /^all$/i){
            foreach my $k (keys %pred){
                foreach my $j ( keys %{$pred{$k}} ){
                    push @{ $filters{$k} }, $j if $pred{$k}->{$j} eq 'default';
                }
            }
        }elsif ( exists $pred{$prog} ){
            if (not $label){
                foreach my $j ( keys %{$pred{$prog}} ){
                    push @{ $filters{$prog} }, $j if $pred{$prog}->{$j} =~ /(default|score)/;
                }
            }else{
                if ($label =~ /^\d+(\.\d+)*$/){#add scores
                    push @{ $filters{$prog} } , $label;
                }else{
                    if (not exists $pred{$prog}->{lc($label)}){
                        croak <<EOT
ERROR: Unrecognised filter ('$label') for $prog passed to --damaging argument.
See --help/--manual for more info.
EOT
;
                    }
                    push @{ $filters{$prog} } , lc($label);
                }
            }
        }else{
            croak <<EOT
ERROR: Unrecognised value ($d) passed to in silico predictions (--damaging) argument. 
See --help/--manual for more info.
EOT
;
        }
    }
    return %filters;
}

=item B<getEvalFilter>

Takes an expression to eval and creates a hash that can be used by the 'evalFilter' function of this module. Expressions must take the format of 'field name' 'comparator' 'value to compare' separated by white space. Multiple expressions can be used together along with the logical operators 'and', 'or' or 'xor'.

The 'field name' component is used to extract the value from the corresponding key in the values hash passed to the 'evalFilter' function. The main purpose for this is to be able to use an eval expression to evaluate values from VEP/SnpEff consequences retrieved using the 'getVepFields' method from VcfReader. 

 my %exp = VcfhacksUtils::getEvalFilter("ada_score > 0.6");
 
 my %exp = VcfhacksUtils::getEvalFilter("LoF eq 'HC'");
 
 my %exp = VcfhacksUtils::getEvalFilter
 (
     "(ada_score > 0.6 and rf_score > 0.6) or maxentscan_diff > 5"
 );

=cut

sub getEvalFilter{
    my $s = shift;
    my $open_brackets = () = $s =~ /\(/g;   
    my $close_brackets = () = $s =~ /\)/g;
    if ($open_brackets > $close_brackets){
        croak "ERROR: Unclosed brackets in expression '$s' passed to --eval_filter\n";
    }elsif ($open_brackets < $close_brackets){
        croak "ERROR: Trailing brackets in expression '$s' passed to --eval_filter\n";
    }
    my %exp = 
    (
        field => [], 
        value => [], 
        comparator => [], 
        operator => [], 
    );
    $s =~ s/^\s+//;#remove preceding whitespace
    $s =~ s/\s+$//;#remove trailing whitespace
    my @split = split(/\s+/, $s);#split on whitespace
    for(my $i = 0; $i<@split; $i++){
    #every fourth array element must be an operator (and/or/xor)
        if ($i and $i % 4 == 3){
            if ($split[$i] !~ /^(and|or|xor|\|\||\&\&)$/){#check operator
                my $pos = $i + 1;
                croak "ERROR: Expected operator (or/and/xor) at position $pos in --eval_filter argument '$s' but got '$split[$i]'\n";
            }
            push @{$exp{operator}}, $split[$i];
        }elsif($i % 4 == 2){ #third arg must be value
            push @{$exp{value}}, $split[$i];
        }elsif($i % 4 == 1){ #second arg must be comparator
            if ($split[$i] !~ /^(eq|ne|[!=]~|[<>][=]{0,1}|[=!]=)$/){#check comparator
                my $pos = $i + 1;
                croak "ERROR: Expected comparator at position $pos in --eval_filter argument '$s' but got '$split[$i]'\n";
            }
            push @{$exp{comparator}}, $split[$i];
        }else{#0th arg is consequence field
            push @{$exp{field}}, $split[$i];
        }
    }
    return %exp;
}

=item B<evalFilter>

Using a hash created using the getEvalFilter function (above), and a hash of field names to values, this function returns the result of evalutaing the resulting expressions using perl's eval function. The main purpose for this is to query VEP/SnpEff consequences retrieved using VcfReader (e.g. return true if a score is above a certain threshold or a string is matched).

 my %exp = VcfhacksUtils::getScoreFilter("ada_score > 0.6 or LoF eq 'HC'");
 my @vep_csq = VcfReader::getVepFields
 (
    line       => \@split_line,
    vep_header => \%vep_header,
    field      => 'all',
 );
 foreach my $csq (@vep_csq){
     if (VcfhacksUtils::scoreFilter(\%exp, $csq)){
         ...
     }
 }

=cut 


sub evalFilter{
    my ($exps, $vals) = @_;
    my @eval = (); 
    for (my $i = 0; $i < @{$exps->{field}}; $i++){
        (my $field = $exps->{field}->[$i]) =~ s/^\(//;#remove preceding bracket
        my $v = $vals->{$field};
        push @eval, "(" if $exps->{field}->[$i] =~ /^\(/;
        if ($exps->{comparator}->[$i] =~ /^([<>][=]{0,1}|[=!]=)$/){
            if ($v ne ''){
                push @eval, $v; 
                push @eval, "$exps->{comparator}->[$i]";
                push @eval, "$exps->{value}->[$i]";
            }else{
            #if empty convert to false and eval rest of expression
                push @eval, 0 ;
                push @eval, ")" if $exps->{value}->[$i] =~ /\)$/;
            }
        }else{
            push @eval, "'$v'"; 
            push @eval, "$exps->{comparator}->[$i]";
            push @eval, "$exps->{value}->[$i]";
        }
        if ($i < @{$exps->{operator}}){
            push @eval, $exps->{operator}->[$i];
        }
    }
    my $ev =eval join(" ", @eval);
    carp "$@\n" if $@;
    return $ev;
}


=back

=head2 Misc Utilities

=over 12

=item B<removeDups>

Remove duplicate entries from an array.

 @uniq = VcfhacksUtils::removeDups(@array); 

=cut

sub removeDups{
    my %seen = ();
    return grep { ! $seen{$_}++ } @_;
}

=item B<getProgressBar>

Given an input file return a progress bar and total number of variants if 
Term::Progressbar is installed and input file is not STDIN, pipe or huge. 
Otherwise returns 0.

Arguments

=over 16

=item input

Input VCF filename or filehandle. Required.

=item total

Number to use as total for progressbar - providing this argument stops any 
calculation from input and will return a progressbar regardless of input type 
or size.

=item name

Name for progressbar. Default = "processing".

=item factor

Multiply total count calculated from input by this number. Default = 1.

=back

 my ($bar, $total) = VcfhacksUtils::getProgressBar(input => $vcf, namne =>  "filtering");

=cut

sub getProgressBar{
    my %args = @_;
    croak "input or total argument is required for getProgressBar method " 
      if not defined $args{input} and not defined $args{total};
    $args{name}   ||= "Processing";
    $args{factor} ||= 1;
    my $msg = <<EOT
WARNING: 'Term::ProgressBar' module was not found. Progress will be displayed 
as the number of variants processed.
EOT
    ;
    eval "use Term::ProgressBar; 1 " or informUser($msg) and return;
    if ($args{total}){
        return (
            Term::ProgressBar->new
            (
                {
                    name  => $args{name},
                    count => $args{total} * $args{factor},
                    ETA   => "linear",
                }
            ), 
            $args{total}
            );
    }
    if ( $args{input} eq "-" or -p $args{input}) {
        informUser(
            "Progress will be shown as counts of variants processed ".
            "when input is from STDIN or pipe\n"
        );
    }elsif( ($args{input} =~ /\.(b)*gz$/ and -s $args{input} > 500*1024*1024 ) or 
            ($args{input} !~ /\.(b)*gz$/ and -s $args{input} > 4*1024*1024*1024 ) ) {
        informUser( "Input file size is large - progress will be shown as ".
                    "counts of variants processed.\n"
        );
    }else{
        informUser("Counting variants in input for progress monitoring.\n"); 
        my $total_vars = VcfReader::countVariants($args{input});
        informUser("$args{input} has $total_vars variants.\n");
        return (
            Term::ProgressBar->new
            (
                {
                    name  => $args{name},
                    count => $total_vars * $args{factor},
                    ETA   => "linear",
                }
            ),
            $total_vars
        );
    }
    return 0;
}

#################################################
sub informUser{
    my $msg = shift;
    my $time = strftime( "%H:%M:%S", localtime );
    print STDERR "[INFO - $time] $msg";
}

=item B<simpleProgress>

Print an updating progress count to STDERR
    
Requires the current count as first argument. Takes optional second and third 
arguments:

The second argument determines whether to increase the count only 
proportional to the total count (default, or if second argument evaluates to 
false) or to increment granularly (if second argument evaluates to true). 

The third argument is the message to output alongside the count (default = 
"variants processed...").

 while (<$FH>){
     ...
     VcfhacksUtils::simpleProgress(++$n, 0, "variants converted...");
 } 

=cut

sub simpleProgress{
    my ($count, $granular, $msg) = @_;
    $msg ||= "variants processed...";
    local $| = 0;
    my $mod = 1;
    if ($granular){
        $mod = 0;
    }elsif ($count >= 100000){#update at least every 100 counts
        $mod = 100;
    }elsif($count > 0){
        my $exp = int ( log($count)/log(10)); 
        $mod *= 10**($exp-2);
    }
    if ($mod  < 1 or not $count % $mod){
        print STDERR "\r$count $msg";
    }
}

=back

=cut

1;
