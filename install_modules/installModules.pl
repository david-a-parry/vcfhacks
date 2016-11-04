#!/usr/bin/env perl
use strict;
use warnings;
use CPAN;

use CPAN::Shell;
my @modules = qw /
    Parallel::ForkManager
    Sys::CPU
    Term::ProgressBar
    LWP::Simple 
    HTTP::Tiny 
    JSON 
    Excel::Writer::XLSX 
    Bio::Perl
/ ;

print STDERR "Attempting to start and configure CPAN...\n";

#CPAN::HandleConfig->load;
CPAN::HandleConfig->load
(
    doit => 1, 
    autoconfig => 1,
    prerequisites_policy => "follow",
    build_requires_install_policy => "yes",
);
CPAN::HandleConfig->commit;
CPAN::Shell::setup_output;
CPAN::Index->reload;
#CPAN::Shell->o("conf prerequisites_policy follow");
#CPAN::Shell->o("conf build_requires_install_policy yes");

print STDERR "Attempting to install the following modules with CPAN:\n" . 
 join("\n", @modules) . "\n";

my @fails = ();
foreach my $m (@modules){
    print STDERR "Attempting to install $m with CPAN::Shell...\n";
    my $mod = CPAN::Shell->expand("Module", $m);
    $mod->install;
    if (not $mod->uptodate){
        warn "Error installing $m\n";
        push @fails, $m;
    }
}


doNonCpan();

sub doNonCpan{
    my $mod = CPAN::Shell->expand("Module", "Bio::DB::HTS"); 
    chomp (my $get = `which wget`) ;
    if (not $get){
        chomp ($get = `which curl`) ;
        $get .= ' -O ' if ($get);
    }
    die "ERROR: could not find either curl or wget executables in your PATH.\n" 
     if not $get;
    if (not $mod->uptodate){
        print STDERR "Attempting to install Bio::DB::HTS...\n";
        print STDERR "Retrieving tarball from cpan...\n";
        my $success  = 0;
        if (runSystemCommand("$get http://search.cpan.org/CPAN/authors/id/R/RI/RISHIDEV/Bio-DB-HTS-2.5.tar.gz")){
            if (runSystemCommand("tar xvf Bio-DB-HTS-2.5.tar.gz")){
                chdir "Bio-DB-HTS-2.5" or die "Could not cd to Bio-DB-HTS-2.5 dir: $!\n";
                $success = runSystemCommand("perl INSTALL.pl");
            }
        }
        push @fails, "Bio::DB::HTS" if not $success;
    }

    $mod = CPAN::Shell->expand("Module", "Bio::DB::Sam"); 
    if (not $mod->uptodate){
        my $success  = 0;
        if (runSystemCommand("$get http://search.cpan.org/CPAN/authors/id/L/LD/LDS/Bio-SamTools-1.43.tar.gz")){
            if (runSystemCommand("tar xvf Bio-SamTools-1.43.tar.gz")){
                chdir "Bio-SamTools-1.43" or die "Could not cd to Bio-SamTools-1.43 dir: $!\n";
                $success = runSystemCommand("perl INSTALL.pl");
            }
        }
        push @fails, "Bio::DB::Sam" if not $success;
    }
}




if (@fails){
    warn scalar @fails . " modules failed to install. You may need to install ".
     "these modules manually and investigate the source of the error.\n";
    warn "The modules that failed were:\n" . join("\n", @fails) . "\n";
}

sub runSystemCommand{
    my $cmd = shift;
    system($cmd);
    return checkExit($?);
}

sub checkExit{
    my $e = shift;
    if($e == -1) {
        print "Failed to execute: $!\n";
        return 0;
    }elsif($e & 127) {
        printf "Child died with signal %d, %s coredump\n",
        ($e & 127),  ($e & 128) ? 'with' : 'without';
        return 0;
    }elsif($e != 0){
        printf "Child exited with value %d\n", $e >> 8;
        return 0;
    }
    return 1;
}
