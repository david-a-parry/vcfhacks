use strict;
use warnings;
use Test::More tests => 68;
use FindBin qw($RealBin);
my $n_tests = 0;

opendir (my $DIR, "./") or die "Cannot read current directory: $!\n";
my @files = readdir($DIR); 
close $DIR;
opendir (my $ADIR, "accessories") or die "Cannot read accessories dir: $!\n";
push @files, map {"accessories/$_"} readdir($ADIR);
close $ADIR;

@files = grep {$_ !~ /(readme|examples)/ } 
         grep { -f $_ }
@files;

foreach my $f (@files){
    ok( -x $f, "$f is executable");
    my $output = `./$f -h 2>&1`;
    my $ex = check_exit($?);
    is( 1, $ex, "$f executes");
}

sub check_exit{
    my $e = shift;
    if($e == -1) {
        print "Failed to execute: $!\n";
        return $e;
    }elsif($e & 127) {
        printf "Child died with signal %d, %s coredump\n",
        ($e & 127),  ($e & 128) ? 'with' : 'without';
        return 0;
    }elsif($e != 0){
        printf "Child exited with value %d\n", $e >> 8;
        return 1;
    }
    return 1;
}
