#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy::Recursive qw(rcopy);
use File::Path qw (remove_tree);
use FindBin qw($Bin);

my @needs_tabix =  
qw (
    annotateSnps.pl
    filterOnEvsMaf.pl
    filterOnSample.pl
    filterVcfOnVcf.pl
    getVariantsByLocation.pl
    rankOnCaddScore.pl
    sampleCallsToInfo.pl
);
my @needs_sort_external = 
qw (
    sortVcf.pl
    rankOnCaddScore.pl
);

my $dir = "$Bin/../";
opendir (my $DIR, $dir) or die "Cannot read current directory $dir: $!\n";
my @files = readdir($DIR); 
close $DIR;
@files = grep {$_ !~ /for_binaries|make_binaries/ } @files;
chdir $dir or die "Cannot move to directory $dir: $!\n";
my $bin_dir = "for_binaries_$^O";
mkdir($bin_dir) or die "$!\n";
foreach my $f (@files){
    next if $f =~ /^\./;
    if ($f =~ /\.pl$/){
        print STDERR "Making refactored copy of $f...\n";
        (my $exe = $f) =~ s/\.pl$//;
        my $out = "$bin_dir/$f";
        open (my $IN, $f) or die "Can't read file $f: $!\n";
        open (my $OUT, ">$out") or die "Can't open $out for writing: $!\n";
        while (my $line = <$IN>){
            $line =~ s/$f/$exe/g; 
            $line =~ s/pod2usage\s*\(/pod2usage(noperldoc => 1, /;
            print $OUT $line;
        }
        close $IN;
        close $OUT;
    }else{
        print STDERR "Copying $f...\n";
        rcopy($f, "$bin_dir/$f") or die "error copying $f: $!\n"; 
    }
}
#we only make binaries now because we need to be sure that 
#the libs folder has already been copied
foreach my $f (@files){
    next if $f =~ /^\./;
    if ($f =~ /\.pl$/){
        (my $exe = $f) =~ s/\.pl$//;
        chdir($bin_dir);
        my $pp = "pp --lib lib/ --lib lib/Bioperl --lib lib/BioASN1EntrezGene/lib";
        if ($^O eq 'darwin'){
            $pp = "pp5.16";
        }
        my $pp_cmd = "$pp -c $f -o $exe";
        if (grep {$_ eq $f} @needs_tabix){
            $pp_cmd .= " -M Tabix";
        }
        if (grep {$_ eq $f} @needs_sort_external){
            $pp_cmd .= " -M Sort::External";
            if ($^O eq 'darwin'){
                $pp_cmd .= " --lib /opt/local/lib/perl5/site_perl/5.16.3/darwin-thread-multi-2level/";
            
            }
        }
        if ($f eq "geneAnnotator.pl"){
            $pp_cmd .= " -M Bio::SeqIO::entrezgene" . 
                       " -M HTTP::Tiny -M JSON";
        } 
        print STDERR "Making binary with command: $pp_cmd\n";
        system($pp_cmd); 
        if ($?){
            print STDERR "WARNING - $pp_cmd exited with status $?\n";
        }else{
            print STDERR "Done.\n"; 
        }
        unlink $f or warn "Error removing $f from $bin_dir: $!\n"; 
        chdir("..");
    }
}
print STDERR "Cleaning up $bin_dir...\n";
chdir($bin_dir);
foreach my $f (@files){
    next if $f =~ /^\./;
    if ($f !~ /\.pl$/){
        if ( -d $f and $f ne 'data' ){
            print STDERR "Recursively removing directory $f.\n";
            remove_tree($f, {verbose => 1} );
        }else{
            if (-e $f){
                next if ($f eq 'examples_bin.md' or $f eq 'readme_binaries.md');
                print STDERR "Removing file $f.\n";
                unlink $f or warn "ERROR removing $f: $!\n";
            }
        }
    }
}
if (-e 'examples_bin.md'){
    rename 'examples_bin.md', 'examples.md';
    rename 'readme_binaries.md', 'readme.md';
}
