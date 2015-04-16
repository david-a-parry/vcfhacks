#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy::Recursive qw(rcopy);
use File::Path qw (remove_tree);
use FindBin qw($Bin);

my $dir = "$Bin/../";
opendir (my $DIR, $dir) or die "Cannot read current directory $dir: $!\n";
my @files = readdir($DIR); 
my @scripts = ();
close $DIR;
@files = grep {$_ !~ /for_binaries|make_binaries|scripts/ } @files;
chdir $dir or die "Cannot move to directory $dir: $!\n";
my $out_dir = "scripts_$^O";
mkdir($out_dir) or die "$!\n";
foreach my $f (@files){
    next if $f =~ /^\./;
    if ($f =~ /\.pl$/){
        push @scripts, $f;
        print STDERR "Making refactored copy of $f...\n";
        (my $exe = $f) =~ s/\.pl$//;
        my $out = "$out_dir/$f";
        open (my $IN, $f) or die "Can't read file $f: $!\n";
        open (my $OUT, ">$out") or die "Can't open $out for writing: $!\n";
        while (my $line = <$IN>){
            if ($line =~ /^use strict/){
                print $OUT $line;
                print $OUT "use PAR;\nuse lib 'par/vcfhacks_lib.par';\n";
            }elsif($line =~ /^use lib|^use FindBin/){
                #don't print FindBin lines
            }else{
                $line =~ s/pod2usage\s*\(/pod2usage(noperldoc => 1, /;
                print $OUT $line;
            }
        }
        close $IN;
        close $OUT;
    }else{
        print STDERR "Copying $f...\n";
        rcopy($f, "$out_dir/$f") or die "error copying $f: $!\n"; 
    }
}
#we only make binaries now because we need to be sure that 
#the libs folder has already been copied
chddir("par") or die "$!";
my $pp = "pp";
if ($^O eq 'darwin'){
    $pp = "pp5.16";
    $pp .= " --lib /opt/local/lib/perl5/site_perl/5.16.3/darwin-thread-multi-2level/";
}
my $pp_cmd = "$pp -p -x dummy.pl -o vcfhacks_lib.par";
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

chdir("..") or die "$!";
print STDERR "Cleaning up $out_dir...\n";
chdir($out_dir);
foreach my $f (@files){
    next if $f =~ /^\./;
    if ($f !~ /\.pl$/){
        if ( -d $f and $f ne 'par'){
            print STDERR "Recursively removing directory $f.\n";
            remove_tree($f, {verbose => 1} );
        }else{
            if (-e $f){
                if ($f eq 'examples_bin.md' or $f eq 'readme_binaries.md'){
                    print STDERR "Removing file $f.\n";
                    unlink $f or warn "ERROR removing $f: $!\n";
                }
            }
        }
    }
}
print STDERR "Checking script syntax...\n";
foreach my $s (@scripts){
    system("perl -c $s");
}

