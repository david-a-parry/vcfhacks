#!/usr/bin/env perl
use strict;
use warnings;
use CPAN;
use CPAN::Shell;
use File::Path qw (remove_tree);

BEGIN {
    $ENV{'AUTOMATED_TESTING'}      = 1;
    $ENV{'NONINTERACTIVE_TESTING'} = 1;
    $ENV{'PERL_MM_USE_DEFAULT'}    = 1;
    $ENV{'PERL_MM_NONINTERACTIVE'} = 1;
    $ENV{'BAIL_ON_FAIL'}           = 1;
}    

my @modules = qw /
    Parallel::ForkManager
    Term::ProgressBar
    List::MoreUtils
    LWP::Simple 
    HTTP::Tiny 
    JSON 
    Excel::Writer::XLSX 
    Sort::External
    Statistics::R
/ ;

print STDERR "Attempting to start and configure CPAN...\n";

CPAN::HandleConfig->load
(
    doit => 1, 
    autoconfig => 1,
    prerequisites_policy => "follow",
    build_requires_install_policy => "yes",
);
#CPAN::HandleConfig->load;
CPAN::HandleConfig->commit;
CPAN::Shell::setup_output;
CPAN::Index->reload;
#CPAN::Shell->o("conf prerequisites_policy follow");
#CPAN::Shell->o("conf build_requires_install_policy yes");

print STDERR "Attempting to install the following modules with CPAN:\n" . 
 join("\n", @modules) . "\n";

my @fails = ();
foreach my $m (@modules){
    my $mod = CPAN::Shell->expand("Module", $m);
    if (defined $mod and $mod->uptodate){
        print STDERR "$m is up to date. Skipping.\n";
        next;
    }
    print STDERR "Attempting to install $m with CPAN...\n";
    #system("perl -MCPAN -e 'install $m'"); 
    CPAN::Shell->install($m);
    #$mod->install;
    $mod = CPAN::Shell->expand("Module", $m);
    if (not defined $mod or not $mod->uptodate){
        warn "Error installing $m\n";
        push @fails, $m;
    }
}

#require Bio::Perl for remaining modules - but let's not be fussy about version
my $bpver = eval
{
    require Bio::Perl ; 
    Bio::Perl->import();
    $Bio::Perl::VERSION 
};
if ($bpver){
    print STDERR "Bio::Perl version $bpver installed - skipping...\n"
}else{
    CPAN::Shell->install("Bio::Perl");
    if (not checkExit($?)){
        push @fails, 'Bio::Perl';
    }
}

doNonCpan();

if (@fails){
    warn scalar @fails . " modules failed to install. You may need to install ".
     "these modules manually and investigate the source of the error.\n";
    warn "The modules that failed were:\n" . join("\n", @fails) . "\n";
}else{
    print STDERR "Done installing required modules for vcfhacks.\n";
}

##################################################
sub doNonCpan{
    chomp (my $get = `which wget`) ;
    if (not $get){
        chomp ($get = `which curl`) ;
        $get .= ' -O ' if ($get);
    }
    die "ERROR: could not find either curl or wget executables in your PATH.\n" 
     if not $get;
    my $ok = eval
    {
        require Bio::DB::HTS ; 
        Bio::DB::HTS->import();
        $Bio::DB::HTS::VERSION >= 2.5
    };
    if ( $ok){
        print STDERR "Bio::DB::HTS version 2.5 or higher already installed - skipping.\n";
    }else{
        print STDERR "Attempting to install Bio::DB::HTS...\n";
        print STDERR "Retrieving tarball from cpan...\n";
        my $success  = 0;
        if (runSystemCommand("$get http://search.cpan.org/CPAN/authors/id/R/RI/RISHIDEV/Bio-DB-HTS-2.7.tar.gz")){
            if (runSystemCommand("tar xvf Bio-DB-HTS-2.7.tar.gz")){
                chdir "Bio-DB-HTS-2.7" or die "Could not cd to Bio-DB-HTS-2.7 dir: $!\n";
                $success = runSystemCommand("perl INSTALL.pl");
                chdir "..";
                remove_tree("Bio-DB-HTS-2.7") 
                 or warn "Error removing Bio-DB-HTS-2.7 directory: $!\n";
            }
            unlink("Bio-DB-HTS-2.7.tar.gz") or warn "Error removing Bio-DB-HTS-2.7.tar.gz: $!\n";
        }
        push @fails, "Bio::DB::HTS" if not $success;
    }

    $ok = eval 
    {
        require Bio::DB::Sam ; 
        Bio::DB::Sam->import();
        $Bio::DB::Sam::VERSION >= 1.39
    };
    if ($ok){
        print STDERR "Bio::DB::Sam version 1.39 or higher already installed - skipping.\n";
    }else{
        my $success  = 0;
        if (runSystemCommand("$get http://search.cpan.org/CPAN/authors/id/L/LD/LDS/Bio-SamTools-1.43.tar.gz")){
            if (runSystemCommand("tar xvf Bio-SamTools-1.43.tar.gz")){
                chdir "Bio-SamTools-1.43" or die "Could not cd to Bio-SamTools-1.43 dir: $!\n";
                $success = runSystemCommand("perl INSTALL.pl");
                chdir "..";
                remove_tree("Bio-SamTools-1.43") 
                 or warn "Error removing Bio-SamTools-1.43 directory: $!\n";
            }
            unlink("Bio-SamTools-1.43.tar.gz") or warn "Error removing Bio-SamTools-1.43.tar.gz: $!\n";
        }
        push @fails, "Bio::DB::Sam" if not $success;
    }
}

##################################################
sub runSystemCommand{
    my $cmd = shift;
    system($cmd);
    return checkExit($?);
}

##################################################
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
