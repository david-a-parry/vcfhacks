#!/usr/bin/env perl
use strict;
use warnings;
use CPAN;

use CPAN::Shell;
my @modules = qw /
    Parallel::ForkManager
    Sys::CPU
    Term::ProgressBar
    Bio::DB::HTS::Tabix
    LWP::Simple 
    HTTP::Tiny 
    JSON 
    Excel::Writer::XLSX 
/ ;

print STDERR "Attempting to start and configure CPAN...\n";

CPAN::HandleConfig->load;
CPAN::Shell::setup_output;
CPAN::Index->reload;

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

if (@fails){
    warn scalar @fails . " modules failed to install. You may need to install ".
     "these modules manually and investigate the source of the error.\n";
    warn "The modules that failed were:\n" . join("\n", @fails) . "\n";
}
