#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;
use Term::ProgressBar;
use Data::Dumper;
use POSIX qw/strftime/;
use List::Util qw(sum);
use FindBin;
use lib "$FindBin::Bin";
use VcfReader;

my $minGQ = 0;
my %opts = ();
GetOptions(
    \%opts,
    "output=s",     
    "input=s",
    "genotype_quality=i",
    "keep_calls",
    "help",
    "manual",  
    "progress",
) or pod2usage( -exitval => 2, -message => "Syntax error" );
pod2usage( -verbose => 2 ) if $opts{manual};
pod2usage( -verbose => 1 ) if $opts{help};
pod2usage( -exitval => 2, -message => "Syntax error" )  if not $opts{input};
if ($opts{genotype_quality}){
    $minGQ = $opts{genotype_quality} 
}
my $OUT;
if ( $opts{output} ) {
    open( $OUT, ">$opts{output}" )
      || die "Can't open $opts{output} for writing: $!";
}
else {
    $OUT = *STDOUT;
}

my $progressbar;
my $next_update = 0;
my $total_vcf   = 0;
my $time = strftime( "%H:%M:%S", localtime );
print STDERR "[$time] Initializing input VCF... ";
die "Header not ok for input ($opts{input}) "
    if not VcfReader::checkHeader( vcf => $opts{input} );
if ( defined $opts{progress} ) {
    $total_vcf = VcfReader::countVariants( $opts{input} );
    print STDERR "$opts{input} has $total_vcf variants. ";
}


$time = strftime( "%H:%M:%S", localtime );
print STDERR "\n[$time] Finished initializing input VCF\n";

my @meta_head = VcfReader::getMetaHeader( $opts{input} );
my %info_fields = VcfReader::getInfoFields(header => \@meta_head);
print $OUT join("\n", @meta_head) ."\n";

if (not exists $info_fields{AC}){
    print $OUT "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes, for each ALT allele, in the same order as listed\">\n";
}
if (not exists $info_fields{AF}){
    print $OUT "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency, for each ALT allele, in the same order as listed\">\n";
}

print $OUT 
    "##INFO=<ID=PGTS,Number=G,Type=String,Description="
    ."\"Possible Genotype Call Codes (for ease of reference).\">\n"
    ."##INFO=<ID=GTC,Number=G,Type=String,Description="
    ."\"Genotype counts in order of genotypes listed by PGTS field.\">\n";
print $OUT "##sampleCallsToInfo.pl=\"";

my @opt_string = ();
foreach my $k ( sort keys %opts ) {
    if ( not ref $opts{$k} ) {
        push @opt_string, "$k=$opts{$k}";
    }
    elsif ( ref $opts{$k} eq 'SCALAR' ) {
        if ( defined ${ $opts{$k} } ) {
            push @opt_string, "$k=${$opts{$k}}";
        }
        else {
            push @opt_string, "$k=undef";
        }
    }
    elsif ( ref $opts{$k} eq 'ARRAY' ) {
        if ( @{ $opts{$k} } ) {
            push @opt_string, "$k=" . join( ",", @{ $opts{$k} } );
        }
        else {
            push @opt_string, "$k=undef";
        }
    }
}
print $OUT join( " ", @opt_string ) . "\"\n";
if (defined $opts{keep_calls}){
    my $col_header = VcfReader::getColumnHeader($opts{input});
    print $OUT "$col_header\n";
}else{
    print $OUT "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
}
$time = strftime( "%H:%M:%S", localtime );
print STDERR "[$time] Conversion starting\n";

if ( defined $opts{progress} and $total_vcf ) {
    $progressbar = Term::ProgressBar->new(
        { name => "Converting", 
          count => ($total_vcf), 
          ETA => "linear" 
        } 
    );
}

my $n                = 0;
my $VCF              = VcfReader::_openFileHandle( $opts{input} );

VAR: while ( my $line = <$VCF> ) {
    next if $line =~ /^#/;
    $n++;
    chomp $line;
    my @split_line = split( "\t", $line );
    my $l = convertCallsToInfo( \@split_line);
    print $OUT "$l\n";
    if ($progressbar) {
        $next_update = $progressbar->update($n) if $n >= $next_update;
    }
}
close $VCF;
close $OUT;
if ( defined $opts{progress} ) {
    $progressbar->update($total_vcf) if $total_vcf >= $next_update;
}
$time = strftime( "%H:%M:%S", localtime );
print STDERR "\nTime finished: $time\n";


################################################
#####################SUBS#######################
################################################

sub convertCallsToInfo{
    my ($line) = @_;
    my $info = VcfReader::getVariantField($line, "INFO");
    my $ac = VcfReader::getVariantInfoField($line, "AC");
    my $af = VcfReader::getVariantInfoField($line, "AF");
    #my $an = VcfReader::getVariantInfoField($line, "AN");
    if (not $ac or not $af){
        my %counts = VcfReader::countAlleles
                    (
                        line => $line,
                        minGQ => $minGQ,
                    );
        my @counts = ();
        foreach my $k (sort {$a <=> $b} keys %counts){
            push @counts, $counts{$k};
        }
        if (@counts < 2){#no ALT alleles(?)
            $info .= ";AF=0";
            $info .= ";AC=0";
        }
        my @alt_counts = @counts[1..$#counts];
        if (not $ac){
            $info .= ";AC=" . join(",", @alt_counts) ;
        }
        if (not $af){
            my @freqs = ();
            foreach my $alt (@alt_counts){
                my $f = "0.00";
                eval{
                    $f = $alt/sum(@counts);
                    $f = sprintf("%.3f", $f);
                };
                push @freqs, $f;
            }
            $info .= ";AF=" . join(",", @freqs);
        } 
    }

    my @pgts = VcfReader::getAllPossibleGenotypeCodes($line);    
    $info .= ";PGTS=" . join(",", @pgts);
    my %gts = VcfReader::countGenotypes
                    (
                        line => $line,
                        minGQ => $minGQ,
                    );
    my @gtcs = ();
    foreach my $g (@pgts){
        if (exists $gts{$g}){
            push @gtcs, $gts{$g};
        }else{
            push @gtcs, 0;
        }
    }
    $info .= ";GTC=" . join(",", @gtcs);
    my $l = VcfReader::replaceVariantField($line, 'INFO', $info);
    if (defined $opts{keep_calls}){
        return join("\t", @$l);
    }else{
        return join("\t", @$l[0..7]);
    }
}



=head1 NAME

sampleCallsToInfo.pl - remove sample calls from a VCF and write genotype and allele metrics to INFO fields instead

=head1 SYNOPSIS

sampleCallsToInfo.pl -i <input vcf file> [options]

sampleCallsToInfo.pl -h (show help message)

sampleCallsToInfo.pl -m (show manual page)


=head1 ARGUMENTS

=over 8 

=item B<-i    --input>

Input VCF file.

=item B<-o    --output>

Output file. Optional - defaults to STDOUT.

=item B<-g    --genotype_quality>

Minimum genotype quality.  Sample alleles will only be counted if they have a genotype quality score equal to or higher than this value. Default is 0 (i.e. no filtering on genotype quality).

=item B<-k    --keep_calls>

Use this flag to keep FORMAT and sample genotype fields in the output. Default is to remove FORMAT and all sample genotype fields.

=item B<-p    --progress>

Show a progress bar.

=item B<-h    --help>

Show help message.

=item B<-m    --manual>

Show manual page.

=back

=head1 DESCRIPTION

Converts genotype calls from a VCF into INFO fields in order to reduce file size and/or speed up searching (e.g. with filterVcfOnVcf.pl). Possible genotypes for a given variant are written to the PGTS field, genotype counts (one per PGTS value) are written to the GTC field and if the standard AC or AF fields are not already written to the input these are also added. 

=head1 AUTHOR

David A. Parry
University of Leeds


=head1 COPYRIGHT AND LICENSE

Copyright 2014  David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut
