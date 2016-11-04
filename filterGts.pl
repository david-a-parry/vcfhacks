#!/usr/bin/env perl
use strict;
use warnings;
use POSIX qw/strftime/;
use Data::Dumper;
use FindBin qw($RealBin);
use List::Util qw(first sum max);
use Term::ProgressBar;
use Getopt::Long;
use Pod::Usage;
use lib "$RealBin/lib";
use lib "$RealBin/lib/dapPerlGenomicLib";
use VcfReader 0.3;
use VcfhacksUtils;
my @ab = (); 
my %opts = (a => \@ab);
GetOptions(
    \%opts,
    "i|input=s",
    "o|output=s",
    "a|ab=f{,}",
    "d|depth=i",
    "s|supporting_reads=i",
    "p|progress",
    "h|help|?",  
    "m|manual",  
) or pod2usage( -exitval => 2, -message => "Syntax error" );
pod2usage( -verbose => 2, -exitval => 0 ) if $opts{m};
pod2usage( -verbose => 1, -exitval => 0 ) if $opts{h};
pod2usage( -exitval => 2, -message => "--input is required" ) if not $opts{i};
pod2usage
( 
    -exitval => 2, 
    -message => "--ab and/or --supporting_reads and/or --depth are required" 
) if not @ab and not $opts{d} and not $opts{s};

informUser("Checking input VCF\n");
my ($header, $first_var, $VCF)  = #this method means we can read from STDIN/pipe
 VcfReader::getHeaderAndFirstVariant($opts{i});
die "VCF header not OK for $opts{i}\n" 
 if not VcfReader::checkHeader(header => $header);

my %samples_to_columns = VcfReader::getSamples
(
    header => $header, 
    get_columns => 1
);
my @samp_order = sort {$samples_to_columns{$a} <=> $samples_to_columns{$b}} keys %samples_to_columns;
die "No samples found in VCF!\n" if not %samples_to_columns;
my $OUT;
if ( $opts{o} ) {
    open( $OUT, ">$opts{o}" ) || die "Can't open $opts{o} for writing: $!";
}else {
    $OUT = *STDOUT;
}

my $progressbar;
my $next_update      = 0;
my $total_vcf        = 0;
my $time = strftime( "%H:%M:%S", localtime );

if ( defined $opts{p} and not -p $opts{i} and  $opts{i} ne '-') {
    informUser("Counting variants for progress reporting...\n");
    $total_vcf = VcfReader::countVariants( $opts{i} );
    informUser("$opts{i} has $total_vcf variants.\n");
    $progressbar = Term::ProgressBar->new(
        { name => "Analyzing", 
          count => ($total_vcf), 
          ETA => "linear" 
        } 
    );
}elsif(defined $opts{p}){
    $time = strftime( "%H:%M:%S", localtime );
    informUser("Input is from STDIN or pipe - will report progress per 10000 variants.\n");
}else{
    print STDERR "Processing variants...\n";
}

print_header();
print $OUT filterGts($first_var);
my $n = 1;
if ($progressbar){
    $next_update = $progressbar->update($n) if $n >= $next_update;
}
while (my $l = <$VCF>){
    next if $l =~ /^#/;
    $n++;
    chomp $l;
    print $OUT filterGts($l);
    if ($progressbar){
        $next_update = $progressbar->update($n) if $n >= $next_update;
    }
    
} 
close $VCF;
close $OUT;
################################################
sub filterGts{
    my $line = shift;
    my @split = split( "\t", $line );
    my @gts = ();
    my @calls = VcfReader::getAllSampleVariants(\@split);
    my $format = VcfReader::getVariantField(\@split, "FORMAT");
SAMP: for (my $i = 0; $i < @samp_order; $i++){
        my @ads = VcfReader::getSampleAlleleDepths 
        (
              line => \@split,
              sample => $samp_order[$i],
              sample_to_columns => \%samples_to_columns,
        );
        @ads = map { $_ =~ /^\d+$/ ? $_ : 0 } @ads;
        my $gt = VcfReader::getSampleGenotypeField
        (
            line => \@split, 
            field => "GT",
            sample => $samp_order[$i],
            sample_to_columns => \%samples_to_columns,
        );
        my @alts = split(/[\/\|]/, $gt);
        
        my $depth = sum(@ads);
        if (not $depth){
            #might be absence of AD field, but DP field might be present
            $depth = VcfReader::getSampleGenotypeField
            (
                line => \@split, 
                field  => 'DP',
                sample => $samp_order[$i],
                sample_to_columns => \%samples_to_columns,
            );
        }
        $depth ||= 0;
        if ($opts{d} and $opts{d} > $depth){
            push @gts, filterGenotype($calls[$i], $format);
            next SAMP;
        }
        if (@ab and $depth > 0){
            my $al_b = 0;
            if ($alts[0] eq '0' and $alts[1] =~ /[1-9]/){ #only analyze het samples where one allele is REF
                $al_b = $ads[$alts[1]]/$depth;
                if ($al_b < $ab[0]){
                    push @gts, filterGenotype($calls[$i], $format);
                    next SAMP;
                }elsif (@ab > 1 and $al_b > $ab[1]){
                    push @gts, filterGenotype($calls[$i], $format);
                    next SAMP;
                }
            }
        }
        if ($opts{s}){
            my @non_ref = grep { $_ ne '0' and $_ ne '.' } @alts;
            my $x = 0;
            foreach my $non (@non_ref){
                $x++ if ($ads[$non] < $opts{s});#not enough reads to support allele
            
            }
            if (@non_ref and $x >= @non_ref){#all ALT alleles have to few supporting reads
                push @gts, filterGenotype($calls[$i], $format);
                next SAMP;
            }
        }
        push @gts, $calls[$i]; 
    }
    return join("\t", @split[0..8], @gts) . "\n";
}

################################################
sub filterGenotype{
    my ($call, $format) = @_;
    my @c = split(":", $call); 
    my @f = split(":", $format); 
    my $g = 0;
    $g++ until $f[$g] eq 'GT' or $g > $#f;
    if ($g > $#f){
        warn "Variant missing 'GT' field!\n";
        return $call;
    }
    $c[$g] = './.';
    return join(":", @c); 
}
################################################
sub print_header{
    my $meta_head = join("\n", grep {/^##/} @$header);
    my $headstring = '';
    my $optstring = VcfhacksUtils::getOptsVcfHeader(%opts) . "\n".$header->[-1] ."\n" ; 
    print $OUT "$meta_head\n";
    print $OUT $optstring; 
}
#################################################
sub informUser{
    my $msg = shift;
    my $time = strftime( "%H:%M:%S", localtime );
    if ($progressbar){
        $progressbar->message( "[INFO - $time] $msg" );
    }else{
        print STDERR "[INFO - $time] $msg";
    }
}

###########################################################
=head1 NAME

filterGts.pl - set genotypes to no-calls if below specified allele depths/balance

=head1 SYNOPSIS

        filterGts.pl -i [vcf file] -a [allele balance cutoff]
        filterGts.pl -i [vcf file] -d [minimum total allele depth]
        filterGts.pl -i [vcf file] -s [minimum supporting reads]
        filterGts.pl -h (display help message)
        filterGts.pl -m (display manual page)

=cut

=head1 ARGUMENTS

=over 8

=item B<-i    --input>

Input VCF file. Required.

=item B<-o    --output>

Output file.

=item B<-d    --depth>

Minimum depth. Genotypes with depth below this value will be converted to 
no-calls.

=item B<-b    --allele_balance>

Minimum and optional maximum alt allele ratio per sample call. If one value is
provided this will act as the minimum allele balance cutoff. If a second value
is provided this will be used as a maximum allele balance cutoff. Only 
heterozygous variants where one allele is the reference nucleotide will be 
filtered using this option. Valid values are between 0.00 and 1.00.

=item B<-s     --supporting_reads>

Minimum number of supporting reads for alternative alleles. If all 
non-reference alleles for a genotype have fewer than this number of reads 
supporting them, the genotype will be changed to a no-call.

=item B<-p,    --progress>

Use this flag to show a progress bar while this program is running.

=item B<-h    --help>

Display help message.

=item B<-m    --manual>

Show manual page.

=back 

=cut

=head1 DESCRIPTION

Changes called genotypes in a VCF to no-calls ('./.') while retaining other 
format fields according to user-specified allele depth options.

=head1 AUTHOR

David A. Parry

=head1 COPYRIGHT AND LICENSE

Copyright 2016, David A. Parry

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version. This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along
with this program. If not, see <http://www.gnu.org/licenses/>.

=cut



