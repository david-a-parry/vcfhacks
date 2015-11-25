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
    return "##$caller\"" . join( " ", @opt_string ) . "\"\n";
    
}


=item B<getInfoHeader>

Takes a hash containing the keys 'ID, NUMBER, TYPE and DESCRIPTION' and outputs an appropriately formatted INFO header.

 my %info_field = 
 (
     ID          => "AnInfoField",
     NUMBER      => "A",
     TYPE        => "String",
     DESCRIPTION => "A made up VCF INFO field.";
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
    return "##INFO=<ID=ID=$info{ID},Number=$info{Number},".
      "Type=$info{Type},Description=\"$info{Description}\">";
}


1;
