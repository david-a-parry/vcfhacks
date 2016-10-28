#!/usr/bin/env perl
use warnings;
use strict;
use Parallel::ForkManager;
use Getopt::Long qw(:config no_ignore_case);
use Sys::CPU;
use Pod::Usage;
use Term::ProgressBar;
use Data::Dumper;
use List::MoreUtils qw(first_index);
use POSIX qw/strftime/;
use FindBin qw($RealBin);
use lib "$RealBin/../lib";
use VcfReader;

my %opts = ();
GetOptions(
    \%opts,
    "output|o=s",
    "input|i=s",
    "progress",
    "help|h|?",
    "manual",
    "c|clinvar_vcf=s",
) or pod2usage( -exitval => 2, -message => "Syntax error" );
pod2usage( -verbose => 2 ) if $opts{manual};
pod2usage( -verbose => 1 ) if $opts{help};
pod2usage( -exitval => 2, -message => "Syntax error" )
  if not $opts{input}
  or not $opts{c} ;

my $time = strftime( "%H:%M:%S", localtime );
print STDERR "[$time] INFO - Initializing input VCF...\n";
my $progressbar;
my $next_update      = 0;
my $total_vcf        = 0;
my $n = 0;
my ($header, $first_var, $VCF)  = VcfReader::getHeaderAndFirstVariant($opts{input});
die "Header not ok for input ($opts{input}) "
    if not VcfReader::checkHeader( header => $header );
if ( defined $opts{progress} and not -p $opts{input} and  $opts{input} ne '-') {
    $total_vcf = VcfReader::countVariants( $opts{input} );
    $time = strftime( "%H:%M:%S", localtime );
    print STDERR "\n[$time] INFO - $opts{input} has $total_vcf variants.\n";
}elsif(defined $opts{progress}){
    $time = strftime( "%H:%M:%S", localtime );
    print STDERR "\n[$time] INFO - Input is from STDIN or pipe - will report progress per 10000 variants.\n";
}

$time = strftime( "%H:%M:%S", localtime );
print STDERR "[$time] INFO - Finished initializing input VCF\n";


if ( defined $opts{progress}){
    my $name = "Annotating";
    $progressbar = Term::ProgressBar->new(
        { 
          name => $name,
          count => ($total_vcf), 
          ETA => "linear" 
        } 
    );
    #$progressbar->minor(0);
}

my %search_args = VcfReader::getSearchArguments( $opts{c});
print_header();
processLine($first_var);
VAR: while (my $line = <$VCF> ) {
    processLine($line);
    checkProgress();
}

close $VCF;



################################################
#####################SUBS#######################
################################################
sub processLine{
    my $line = shift;
    return if $line =~ /^#/;
    $n++;
    chomp $line;
    $line =~ s/^chr//;
    my $vcf_line = [split( "\t", $line )];
    my %min_vars   = VcfReader::minimizeAlleles($vcf_line);
    foreach my $allele ( sort { $a <=> $b } keys %min_vars ) {
        if ($min_vars{$allele}->{ORIGINAL_ALT} eq '*'){
            #filter new asterisk allele notation
            next;
        }
        my @hits = VcfReader::searchForPosition
        (
                    %search_args,
                    chrom => $min_vars{$allele}->{CHROM},
                    pos   => $min_vars{$allele}->{POS}
        );
       foreach my $snp_line (@hits) {
            #check whether the snp line(s) match our variant
            my @snp_split = split( "\t", $snp_line );
            if ( my $match = checkVarMatches( $min_vars{$allele}, \@snp_split ) ) {
                my %info = map {$_ => VcfReader::getVariantInfoField( \@snp_split, $_) } 
                qw/
                    CLNALLE
                    CLNHGVS
                    CLNDSDBID
                    CLNSIG
                    CLNDBN
                /;
                #make sure our matching allele is the clinically associated one
                #clinvar lists var details in order of alleles specified by CLNALLE field
                my $al_index = first_index {$_ == $match} split(",", $info{CLNALLE});
                next if $al_index < 0;
                my $sig = (split ",", $info{CLNSIG})[$al_index];
                next if not ( grep {$_ == 4 or $_ == 5} split(/\|/, $sig ) ); #4 or 5 == prob/pathognic
                $min_vars{$allele}->{CLNSIG} = $sig;
                $min_vars{$allele}->{CLNDSDBID} = (split ",", $info{CLNDSDBID})[$al_index];
                $min_vars{$allele}->{CLNDBN}    = (split ",", $info{CLNDBN})[$al_index];
                $min_vars{$allele}->{CLNHGVS}   = (split ",", $info{CLNHGVS})[$al_index];
            }
        }
    }
    my %cln_info = ();
    my $path = 0; 
    foreach my $allele ( sort { $a <=> $b } keys %min_vars ) {
        foreach my $k (
        qw/
            CLNALLE
            CLNHGVS
            CLNDSDBID
            CLNSIG
            CLNDBN
        /
        ){
            if (exists $min_vars{$allele}->{$k}){
                push @{$cln_info{$k}}, $min_vars{$allele}->{$k};
                $path++;
            }else{  
                push @{$cln_info{$k}}, ".";
            }
        }
    }
    return if not $path; 
    foreach my $k ( sort keys %cln_info ) {
        $vcf_line = VcfReader::addVariantInfoField 
        (
            line => $vcf_line,
            id   => "$k",
            value => join(",", @{ $cln_info{$k} } ), 
        );
    }       
    print join("\t", @$vcf_line) . "\n";       
 
}

#################################################
sub checkVarMatches {
    my ( $min_allele, $snp_line ) = @_;
    my %snp_min = VcfReader::minimizeAlleles($snp_line);
    foreach my $snp_allele ( keys %snp_min ) {
        $min_allele->{CHROM} =~ s/^chr//;
        $snp_min{$snp_allele}->{CHROM} =~ s/^chr//;
        next if $min_allele->{CHROM} ne $snp_min{$snp_allele}->{CHROM};
        next if $min_allele->{POS} ne $snp_min{$snp_allele}->{POS};
        next if $min_allele->{REF} ne $snp_min{$snp_allele}->{REF};
        next if $min_allele->{ALT} ne $snp_min{$snp_allele}->{ALT};
        return $snp_allele;
    }
    return 0;
}


################################################
sub print_header{
    my $meta_head = join("\n", grep {/^##/} @$header);
    print "$meta_head\n";
    print <<EOT
##INFO=<ID=CLNHGVS,Number=.,Type=String,Description="Variant names from HGVS.    The order of these variants corresponds to the order of the info in the other clinical  INFO tags.">
##INFO=<ID=CLNALLE,Number=.,Type=Integer,Description="Variant alleles from REF or ALT columns.  0 is REF, 1 is the first ALT allele, etc.  This is used to match alleles with other corresponding clinical (CLN) INFO tags.  A value of -1 indicates that no allele was found to match a corresponding HGVS allele name.">
##INFO=<ID=CLNSRC,Number=.,Type=String,Description="Variant Clinical Chanels">
##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Variant Clinical Significance, 0 - Uncertain significance, 1 - not provided, 2 - Benign, 3 - Likely benign, 4 - Likely pathogenic, 5 - Pathogenic, 6 - drug response, 7 - histocompatibility, 255 - other">
##INFO=<ID=CLNDSDBID,Number=.,Type=String,Description="Variant disease database ID">
##INFO=<ID=CLNDBN,Number=.,Type=String,Description="Variant disease name">
EOT
;   
    print "$header->[-1]\n";
}

################################################
sub checkProgress{
    return if not $progressbar;
    my $do_count_check = shift;
    if ($total_vcf > 0){
        $next_update = $progressbar->update($n) if $n >= $next_update;
    }elsif($do_count_check){#input from STDIN/pipe
        if (not $n % 10000) {
            my $time = strftime( "%H:%M:%S", localtime );
            $progressbar->message( "[INFO - $time] $n variants read" );
        }
    }
}
