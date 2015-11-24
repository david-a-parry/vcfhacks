package IdParser;

use strict;
use warnings;
use Carp;
use HTTP::Tiny;
use JSON;
our $AUTOLOAD;

##################################################

{
        my $_count = 0;
        my %_attrs = (
            _currentId      => ["", "read"],
            _isEnsemblId    => [0, "read"],
            _isTranscript   => [0, "read"],
            _isProtein      => [0, "read"],
            _identifierType => ["Gene Symbol or Other", "read"],
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

##################################################

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


##################################################

sub DESTROY{
    my ($self) = @_;
    $self -> _decr_count( );
}

##################################################

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
    return $self;
}

##################################################

sub parseId {
    my $self = shift;
    my $id = shift;
    if ($id =~ /ENS\w*G\d{11}(.\d+)*/){
        $self->{_isTranscript} = 0;
        $self->{_isEnsemblId} = 1;
        $self->{_isProtein} = 0;
        $self->{_identifierType} = "Ensembl Gene ID";
    }elsif ($id =~ /ENS\w*T\d{11}(.\d+)*/){
        $self->{_isTranscript} = 1;
        $self->{_isEnsemblId} = 1;
        $self->{_isProtein} = 0;
        $self->{_identifierType} = "Ensembl Transcript ID";
    }elsif ($id =~ /CCDS\d+(.\d+)*/){
        $self->{_isTranscript} = 1;
        $self->{_isEnsemblId} = 0;
        $self->{_isProtein} = 0;
        $self->{_identifierType} = "CCDS ID";
    }elsif ($id =~ /ENS\w*P\d{11}(.\d+)*/){
        $self->{_isEnsemblId} = 1;
        $self->{_isTranscript} = 0;
        $self->{_isProtein} = 1;
        $self->{_identifierType} = "Ensembl Protein ID";
    }elsif ($id =~ /[XN][MR]_\d+(.\d+)*/){
        $self->{_isTranscript} = 1;
        $self->{_isEnsemblId} = 0;
        $self->{_isProtein} = 0;
        $self->{_identifierType} = "RefSeq mRNA ID";
    }elsif ($id =~ /[XN]P_\d+(.\d+)*/){
        $self->{_isTranscript} = 0;
        $self->{_isEnsemblId} = 0;
        $self->{_isProtein} = 1;
        $self->{_identifierType} = "RefSeq Protein ID";
    }elsif ($id =~ /^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$/){
        $self->{_isTranscript} = 1;
        $self->{_isEnsemblId} = 0;
        $self->{_isProtein} = 1;
        $self->{_identifierType} = "Uniprot ID";
    }elsif ($id =~ /^\d+$/){
        $self->{_isTranscript} = 0;
        $self->{_isEnsemblId} = 0;
        $self->{_isProtein} = 0;
        $self->{_identifierType} = "Entrez Gene ID";
    }else{
        $self->{_isTranscript} = 0;
        $self->{_isEnsemblId} = 0;
        $self->{_identifierType} = "Gene Symbol or other";
    }
    $self->{_currentId} = $id; 
}

##################################################

1;
