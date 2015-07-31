package ParsePedfile;

use strict;
use warnings;
use Carp;
our $AUTOLOAD;
{
        my $_count = 0;
        my %_attrs = (
        _file => ["", "read/required"],
#        _families => ["", "read"],
#        _samples => ["", "read"],
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
    $self->_convertPedToHash($self->{_file});
    $class -> _incr_count();
    return $self;
}

sub _convertPedToHash{
#%{$self->{samples}} hash has Sample IDs as keys to anon hashes of ped info
#%{$self->{families}} hash has Family IDs as keys to anon arrays of Sample IDs
    my ($self, $ped_file) = @_;
    my %ped = ();
    my %samples = ();
    open (my $IN, $ped_file) or croak "Can't open $ped_file for reading.";
    my $l = 0; 
    while (my $line = <$IN>){   
        $l++;
        chomp $line;
        next if not $line;
        #next if $line =~ /^Family ID\s+Individual ID\s+Paternal ID\s+Maternal ID\s+Sex\s+Phenotype/i;
        next if $line =~ /^#/;#allow comments/headers
        my @f = split (/\s+/, $line);
        croak "Incorrect number of fields in ped file $ped_file line $l: found " .(scalar @f). " but require.6 " 
            if @f != 6;
        croak "Invalid value '$f[5]' for phenotype field at line $l of $ped_file." 
            if ($f[5] !~ /^[012]$/ and $f[5] ne '-9');
        croak "ERROR - Sample $f[1] encountered twice in pedigree!" if exists  $self->{samples}->{$f[1]};
        $self->{samples}->{$f[1]}->{family} = $f[0];
        $self->{samples}->{$f[1]}->{father} = $f[2];
        $self->{samples}->{$f[1]}->{mother} = $f[3];
        $self->{samples}->{$f[1]}->{sex} = $f[4];
        $self->{samples}->{$f[1]}->{phenotype} = $f[5];
        push @{$self->{families}->{$f[0]}}, $f[1];
    }
    croak "No samples found in pedigree $ped_file." if not keys %{$self->{samples}};
    close $IN;
}

sub _checkSampleExists{
    my ($self, $sample) = @_;
    croak "Sample $sample does not exist in pedigree." if not exists $self->{samples}->{$sample};
}

sub _checkFamilyExists{
    my ($self, $family) = @_;
    croak "Family $family does not exist in pedigree." if not exists $self->{families}->{$family};
}

sub getAllFamilies{
    my ($self) = @_;
    return keys %{$self->{families}};
}

sub getAllSamples{
    my ($self) = @_;
    return keys %{$self->{samples}};
}

sub getSamplesFromFamily{
    my ($self, $family) = @_;
    $self->_checkFamilyExists($family);
    my @samples = ();
    foreach my $member (@{$self->{families}->{$family}}){
        push @samples, $member;
    }
    return @samples;
}

sub getFather{
    my ($self, $sample) = @_;
    $self->_checkSampleExists($sample);
    return $self->{samples}->{$sample}->{father};
}

sub getMother{
    my ($self, $sample) = @_;
    $self->_checkSampleExists($sample);
    return $self->{samples}->{$sample}->{mother};
}

sub getSex{
    my ($self, $sample) = @_;
    $self->_checkSampleExists($sample);
    return $self->{samples}->{$sample}->{sex};
}

sub getPhenotype{
    my ($self, $sample) = @_;
    $self->_checkSampleExists($sample);
    return $self->{samples}->{$sample}->{phenotype};
}

sub isAffected{
#return 1 if phenotype is known and is affected
#return 0 if phenotype is not known or is unaffected
    my ($self, $sample) = @_;
    $self->_checkSampleExists($sample);
    return 1 if $self->{samples}->{$sample}->{phenotype} == 2;
    return 0;
}

sub isUnaffected{
#return 1 if phenotype is known and is unaffected
#return 0 if phenotype is not known or is affected
    my ($self, $sample) = @_;
    $self->_checkSampleExists($sample);
    return 1 if $self->{samples}->{$sample}->{phenotype} == 1;
    return 0;
}

sub isMale{
    my ($self, $sample) = @_;
    $self->_checkSampleExists($sample);
    return 1 if $self->{samples}->{$sample}->{sex} == 1;
    return 0;
}

sub isFemale{
    my ($self, $sample) = @_;
    $self->_checkSampleExists($sample);
    return 1 if $self->{samples}->{$sample}->{sex} == 2;
    return 0;
}

sub getFamilyId{
    my ($self, $sample) = @_;
    $self->_checkSampleExists($sample);
    return $self->{samples}->{$sample}->{family};
}

sub getParents{
    my ($self, $sample) = @_;
    $self->_checkSampleExists($sample);
    my @parents = ();
    if ($self->{samples}->{$sample}->{father}){
        push @parents, $self->getFather($sample);
    }
    if ($self->{samples}->{$sample}->{mother}){
        push @parents, $self->getMother($sample);
    }
    return @parents if wantarray;
    return @parents if defined wantarray;
    carp "getParents function called in void context.";
}

sub getAllSiblings{
#return IDs of all siblings including half siblings
    my ($self, $sample) = @_;
    $self->_checkSampleExists($sample);
    my @sibs = ();
    foreach my $p ($self->getParents($sample)){
        push @sibs, $self->getChildren($p);
    }
    return @sibs;
}

sub getFullSiblings{
#return IDs of all FULL siblings
    my ($self, $sample) = @_;
    $self->_checkSampleExists($sample);
    my @sibs = ();
    my %par = ();
    foreach my $p ($self->getParents($sample)){
        push @{$par{$p}}, $self->getChildren($p);
    }
    if (keys(%par) != 2){
        carp "Can't find full siblings for $sample, no. parents = ".
            (keys(%par)) ." not 2.";
        return;
    }
    my %intersect = map {$_ => undef} @{$par{(keys%par)[0]}};
    @sibs = grep {exists ($intersect{$_}) } @{$par{(keys%par)[1]}};
    return @sibs;
}

sub getChildren{
    my ($self, $sample) = @_;
    $self->_checkSampleExists($sample);
    my @children = ();
    foreach my $member ($self->getSamplesFromFamily($self->getFamilyId($sample))){
        if ($self->getMother($member) eq $sample or $self->getFather($member) eq $sample){
            push @children, $member;
        }
    }
    return @children;
}

sub isObligateCarrierRecessive{
#assuming recessive inheritance
#checks if $sample has an affected child
    my ($self, $sample) = @_;
    $self->_checkSampleExists($sample);
    return 1 if $self->isAffected($sample);
    foreach my $c ($self->getChildren($sample)){
        return 1 if $self->isAffected($c);
    }
    return 0;
}

sub getAffectedsFromFamily{
#returns sample IDs of all affected members
    my ($self, $family) = @_;
    $self->_checkFamilyExists($family);
    my @affecteds = ();
    foreach my $member ($self->getSamplesFromFamily($family)){
       push @affecteds, $member if $self->isAffected($member);
    }
    return @affecteds;
}

sub getUnaffectedsFromFamily{
#returns sample IDs of all unaffected members
    my ($self, $family) = @_;
    $self->_checkFamilyExists($family);
    my @unaffecteds = ();
    foreach my $member ($self->getSamplesFromFamily($family)){
       push @unaffecteds, $member if $self->isUnaffected($member);
    }
    return @unaffecteds;

}

sub getAllAffecteds{
#returns sample IDs of all affected members
    my ($self) = @_;
    my @affecteds = ();
    foreach my $sample (keys %{$self->{samples}}){
       push @affecteds, $sample if $self->isAffected($sample);
    }
    return @affecteds;
}

sub getAllUnaffecteds{
#returns sample IDs of all unaffected members
    my ($self) = @_;
    my @unaffecteds = ();
    foreach my $sample (keys %{$self->{samples}}){
       push @unaffecteds, $sample if $self->isUnaffected($sample);
    }
    return @unaffecteds;
}



1; 

=head1 NAME

ParsePed.pm - reads and parses .ped pedigree files.

=head1 SYNOPSIS


$ped = ParsePed -> new(file => $ped);
#(initialise object with a ped file)

my @affected = $ped->getAffectedsFromFamily($family);
#(get IDs of all affected family members)

my @parents = $ped->getParents($sample);

my @siblings = $ped->getFullSiblings($sample);

=head1 DESCRIPTION

=head2 Overview

This module is an attempt at creating a simple object-oriented interface for retrieving information about samples in a PED file.

From the PLINK documentation:

The PED file is a white-space (space or tab) delimited file: the first six columns are mandatory:

     Family ID
     Individual ID
     Paternal ID
     Maternal ID
     Sex (1=male; 2=female; other=unknown)
     Phenotype


Affection status, by default, should be coded:

    -9 missing 
     0 missing
     1 unaffected
     2 affected


=head2 Constructor and initialization

To initialise call the 'new' method specifying a .ped file to read.

$ped = ParsePed -> new(file => "fam1.ped");

The information in the .ped file will be read and stored in $ped for access with the methods detailed below. 

=head2 Class and object methods


=head3 Accessors and mutators

The following features can be accessed using get_[feature] or set using set_[feature], substituting [feature] for the feature of interest. 

=over 8

=item B<file>

filename of .ped.

=back

=head3 Methods

=over 8

=item B<getAllSamples>

Returns all family IDs present in the .ped file.

my @families = $ped->getAllFamilies();

=item B<getAllSamples>

Returns all sample IDs present in the .ped file.

my @samples = $ped->getAllSamples();

=item B<getSamplesFromFamily>

Returns all sample IDs from given family. Takes a family ID as argument.

my @samples = $ped->getSamplesFromFamily($family);

=item B<getFather>

Returns sample ID of given sample's father. Requires sample ID as argument.

my $dad = $ped->getFather($sample);

=item B<getMother>

Returns sample ID of given sample's mother. Requires sample ID as argument.

my $mum = $ped->getMother($sample);

=item B<getSex>

Returns gender code for given sample. Requires sample ID as argument.

my $sex = $ped->getSex($sample);

=item B<getPhenotype>

Returns phenotype code for given sample. Requires sample ID as argument.

my $phenotype = $ped->getPhenotype($sample);

=item B<getFamilyId>

Returns family ID for given sample. Requires sample ID as argument.

my $family = $ped->getFamilyId($sample);

=item B<getParents>

Returns sample's parent's IDs. Requires sample ID as argument. Returns empty array if not parents are found.

my @parents = $ped->getParents($sample);

=item B<getAllSiblings>

Returns full and half siblings for sample. Requires sample ID as argument.

my @siblings = $ped->getAllSiblings($sample);

=item B<getFullSiblings>

Returns full for sample. Requires sample ID as argument.

my @siblings = $ped->getFullSiblings($sample);

=item B<getChildren>

Returns children of sample. Requires sample ID as argument.

my @kids = $ped->getChildren($sample);

=item B<getAffectedsFromFamily>

Returns the IDs of all affected members of a family. Requires family ID as argument.

my @affected = $ped->getAffectedsFromFamily($family);

=item B<getUnaffectedsFromFamily>

Returns the IDs of all unaffected members of a family. Requires family ID as argument.

my @unaffected = $ped->getUnaffectedsFromFamily($family);

=item B<getAllAffecteds>

Returns the IDs of all affected members of a pedigree.

my @affected = $ped->getAffecteds();

=item B<getAllUnaffecteds>

Returns the IDs of all unaffected members of a pedigree.

my @unaffected = $ped->getUnaffecteds();

=item B<isMale>

Returns 1 if gender of sample is known and is male. Returns 0 if gender is unknown or female.  Requires sample ID as argument.

my $is_male = $ped->isMale($sample);

=item B<isFemale>

Returns 1 if gender of sample is known and is female. Returns 0 if gender is unknown or male.  Requires sample ID as argument.

my $is_male = $ped->isFemale($sample);

=item B<isAffected>

Returns 1 if sample is affected, returns 0 if sample is unaffected or value is missing. Requires sample ID as argument.

my $is_affected = $ped->isAffected($sample);

=item B<isUnaffected>

Returns 1 if sample is unaffected, returns 0 if sample is affected or value is missing. Requires sample ID as argument.

my $is_unaffected = $ped->isUnaffected($sample);

=item B<isObligateCarrierRecessive>

Returns 1 if sample has an affected child (or if sample is affected). Requires sample ID as argument.

my $is_carrier = $ped->isObligateCarrierRecessive($sample);

=back


=head1 AUTHOR

David A. Parry

=head1 COPYRIGHT AND LICENSE

Copyright 2014  David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.


