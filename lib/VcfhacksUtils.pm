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
    return _readInSilicoFile($insilico_file);
}


=item B<readSnpEffInSilicoFile>

Return a hash of dbNSFP in silico prediction program names with values being anonymous hashes of prediction scores to either 'default' or 'valid' indicating default behaviour of functional consequence filtering scripts. 

 my %in_silico = VcfhacksUtils::readSnpEffInSilicoFile();

=cut

sub readSnpEffInSilicoFile{
    my $insilico_file = "$data_dir/snpeff_insilico_pred.tsv";
    return _readInSilicoFile($insilico_file);
}

sub _readInSilicoFile{
    my $insilico_file = shift;
    open (my $INS, $insilico_file) or die 
"Could not open in silico classes file '$insilico_file': $!\n";
    my %pred = (); 
    while (my $line = <$INS>){
        next if $line =~ /^#/;
        $line =~ s/[\r\n]//g; 
        next if not $line;
        my @split = split("\t", $line);
        die "Not enough fields in classes file line: $line\n" if @split < 3;
        $pred{lc($split[0])}->{lc($split[1])} = lc($split[2]) ;
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
    my %pred = _readInSilicoFile("$data_dir/$mode"."_insilico_pred.tsv");
    my %filters = (); 
    foreach my $d (@$damaging) {
        $d = lc($d); 
        my ( $prog, $label ) = split( "=", $d );
        if ($mode eq 'snpeff' and $prog ne 'all'){
            $prog = "dbnsfp_$prog" if $prog !~ /^dbnsfp_/;
        }
        if ($prog eq 'all'){
            foreach my $k (keys %pred){
                foreach my $j ( keys %{$pred{$k}} ){
                    push @{ $filters{$k} }, $j if $pred{$k}->{$j} eq 'default';
                }
            }
        }elsif ( exists $pred{$prog} ){
            if (not $label){
                foreach my $j ( keys %{$pred{$prog}} ){
                    push @{ $filters{$prog} }, $j if $pred{$prog}->{$j} eq 'default';
                }
            }else{
                if ($label =~ /^\d+(\.\d+)*$/){#add scores
                    push @{ $filters{$prog} } , $label;
                }else{
                    if (not exists $pred{$prog}->{lc($label)}){
                        die <<EOT;
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

=back

=cut

1;
