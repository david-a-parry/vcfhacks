#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy::Recursive qw(rcopy);
use File::Path qw (remove_tree);
use FindBin qw($Bin);

my $dir = "$Bin/../";
opendir (my $DIR, $dir) or die "Cannot read current directory $dir: $!\n";
chdir("$dir/par") or die "$!";
my $pp = "pp";
my $par = "vcfhacks_lib_$^O.par";
if ($^O eq 'darwin'){
    $pp = "pp5.16";
    $pp .= " --lib /opt/local/lib/perl5/site_perl/5.16.3/darwin-thread-multi-2level/";
}
my $pp_cmd = "$pp -p -x dummy.pl -o $par";
print STDERR "Making PAR with command: $pp_cmd\n";
system($pp_cmd); 
if ($?){
    die "ERROR - $pp_cmd exited with status $?\n";
}else{
    print STDERR "Done.\n"; 
}
print STDERR "Testing PAR file...\n";
system("perl test.pl"); 
if ($?){
    die "ERROR - test.pl exited with status $?\n";
}else{
    print STDERR "Done.\n"; 
}

