use strict;
use warnings;
use Module::CoreList;
use Module::Extract::Use;
use Data::Dumper;

my @core = Module::CoreList->find_modules();

while (my $file = shift){

    my $extor = Module::Extract::Use->new;
    if( $extor->error ) { 
        die $extor->error;
    }
    my @modules = $extor->get_modules( $file );
    foreach my $m (@modules){
        if (grep {$_ eq $m } @core){
            print "Core module: $m\n";
        }else{
            print "Non-core module: $m\n";
        }
    }
    #print "Modules are:\n\n". join("\n", @modules) . "\n";
    my $details = $extor->get_modules_with_details( $file );
    foreach my $detail ( @$details ) {
        my $v = $detail->version || '';
        if ( @{ $detail->imports }){
            printf "%s %s imports %s\n",
                    $detail->module, $v,
                    join ' ', @{ $detail->imports }
        }
    }
}
