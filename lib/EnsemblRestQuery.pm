package EnsemblRestQuery;

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
        _server         => ["http://rest.ensembl.org", "read/write"],
        _default_server => ["http://rest.ensembl.org", "read"],
        _grch37_server  => ["http://grch37.rest.ensembl.org", "read"],
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
    $self->{_http} = HTTP::Tiny->new();
    $class -> _incr_count();
    return $self;
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

sub ensRestQuery{
    my $self = shift;
    my $url = shift;
    my $response = $self->{_http}->get($url, {
          headers => { 'Content-type' => 'application/json' }
    });
    if (not $response->{success}){
        carp "Ensembl REST query ('$url') failed!\n" ;
        return;
    }
    if(length $response->{content}) {
      return decode_json($response->{content});
    }
    carp "No content for Ensembl REST query ('$url')!\n";
    return;
}
 
##################################################

sub getViaXreg{
    my ($self, $id, $species, $object_type) = @_;
    my $endpoint = "/xrefs/symbol/$species/$id?object_type=$object_type";
    my $url = $self->{_server} . $endpoint;
    return $self->ensRestQuery($url); 
}

##################################################

sub getTranscriptViaXreg{
    my ($self, $id, $species) = @_;
    return $self->getViaXreg($id, $species, 'transcript');
}
 
##################################################

sub getGeneViaXreg{
    my ($self, $id, $species) = @_;
    return $self->getViaXreg($id, $species, 'gene');
}
    

##################################################

sub lookUpEnsId{
    my ($self, $id, $expand) = @_;
    my $endpoint = "/lookup/id/$id";
    if ($expand){
        $endpoint .= "?expand=1";
    }
    my $url = $self->{_server} . $endpoint;
    return $self->ensRestQuery($url); 
} 

##################################################

sub getParent{
    my ($self, $id, $expand) = @_;
    my $endpoint = "/lookup/id/$id";
    my $url = $self->{_server} . $endpoint;
    my $hash = $self->ensRestQuery($url);
    return if not $hash;
    if ($hash->{Parent}){
        return $self->lookUpEnsId($hash->{Parent}, $expand);
    }else{
        return;
    }
}

##################################################

sub transcriptFromEnsp{
    my ($self, $id) = @_;
    $self->getParent($id);
}

##################################################

sub geneFromTranscript{
    my ($self, $id) = @_;
    $self->getParent($id);
}

##################################################

sub getXrefs{
    my ($self, %args) = @_;
    if (not $args{id}){
        carp "'id' argument is required for getXrefs method\n";
        return;
    }
    my $endpoint = "/xrefs/id/$args{id}";
    my @extra_args = ();
    if ($args{all_levels}){
        push @extra_args, "all_levels=$args{all_levels}";
    }
    if ($args{db}){
        push @extra_args, "external_db=$args{db}";
    }
    $endpoint .= "?" . join(";", @extra_args); 
    my $url = $self->{_server} . $endpoint;
    return $self->ensRestQuery($url); 
}

##################################################
sub proteinPosToGenomicPos{
    my ($self, %args) = @_;
    if (not $args{id}){
        carp "'id' argument is required for method\n";
        return;
    }
    if (not $args{start}){
        carp "'start' argument is required for method\n";
        return;
    }
    $args{end} ||= $args{start}; 
    my $endpoint = "/map/translation/$args{id}/$args{start}..$args{end}";
    my $url = $self->{_server} . $endpoint;
    return $self->ensRestQuery($url); 
}

##################################################
sub useGRCh37Server{
    my $self = shift;
    $self->{_server} = $self->{_grch37_server};
}

##################################################
sub useDefaultServer{
    my $self = shift;
    $self->{_server} = $self->{_default_server};
}
