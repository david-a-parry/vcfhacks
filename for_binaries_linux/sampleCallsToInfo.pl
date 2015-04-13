#!/usr/bin/perl
use warnings;
use strict;
use Parallel::ForkManager;
use Getopt::Long;
use Pod::Usage;
use Term::ProgressBar;
use Data::Dumper;
use POSIX qw/strftime/;
use Sys::CPU;
use List::Util qw(sum);
use FindBin;
use lib "$FindBin::Bin/lib";
use VcfReader;

my $minGQ = 0;
my $cpus  = Sys::CPU::cpu_count();
my $forks = 0;
my $buffer_size;
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
    "forks=i" => \$forks,
    "cache=i" => \$buffer_size,
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
if ( $forks < 2 ) {
    $forks       = 0;    #no point having overhead of forks for one fork
}
else {
    if ( $forks > $cpus ) {
        print STDERR
"[Warning]: Number of forks ($forks) exceeds number of CPUs on this machine ($cpus)\n";
    }
    if ( not $buffer_size ) {
        $buffer_size = 10000 > $forks * 1000 ? 10000 : $forks * 1000;
        ;
    }
    print STDERR
"[INFO] Processing in batches of $buffer_size variants split among $forks forks.\n";
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
my %contigs = ();
if ($forks){
    $total_vcf *= 3;
    %contigs = VcfReader::getContigOrder($opts{input});
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
print $OUT "##sampleCallsToInfo=\"";

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

my $n   = 0;
my $VCF = VcfReader::_openFileHandle( $opts{input} );

my @lines_to_process = ();
VAR: while ( my $line = <$VCF> ) {
    next if $line =~ /^#/;
    $n++;
    if ($forks ){
        push @lines_to_process, $line;
        if ( @lines_to_process >= $buffer_size ) {
            process_buffer();
            @lines_to_process = ();
        }
    }else{
        chomp $line;
        my @split_line = split( "\t", $line );
        my $l = convertCallsToInfo( \@split_line);
        print $OUT join("\t", @$l) ."\n";
    }
    if ($progressbar) {
        $next_update = $progressbar->update($n) if $n >= $next_update;
    }
}
if ($forks){
    process_buffer();
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
sub process_buffer {
    return if not @lines_to_process;
    my @lines_to_print;
    my $lines_per_slice = @lines_to_process;
    if ( $forks > 0 ) {
        $lines_per_slice =
            int( @lines_to_process / $forks ) > 1
          ? int( @lines_to_process / $forks )
          : 1;
    }
    my @batch = ();

    #get a batch for each thread
    for ( my $i = 0 ; $i < @lines_to_process ; $i += $lines_per_slice ) {
        my $last =
          ( $i + $lines_per_slice - 1 ) < $#lines_to_process
          ? $i + $lines_per_slice - 1
          : $#lines_to_process;
        if ( $i + $lines_per_slice >= @lines_to_process ) {
            $last = $#lines_to_process;
        }
        my @temp = @lines_to_process[ $i .. $last ];
        push @batch, \@temp;
    }
    my $pm = Parallel::ForkManager->new($forks);
    $pm->run_on_finish(    # called BEFORE the first call to start()
        sub {
            my ( $pid, $exit_code, $ident, $exit_signal, $core_dump,
                $data_structure_reference )
              = @_;

            if ( defined($data_structure_reference) ) {
                my %res = %{$data_structure_reference};
                if ( ref $res{keep} eq 'ARRAY' ) {
                    push @lines_to_print, \@{ $res{keep} } if @{ $res{keep} };
                }else{
                    die "Unexpected return from child process $pid:\n". Dumper %res;
                }
                if ($progressbar) {
                    $n += $res{batch_size};
                    $next_update = $progressbar->update($n)
                      if $n >= $next_update;
                }
            }
            else {
                die "ERROR: no message received from child process $pid!\n";
            }
        }
    );
    foreach my $b (@batch) {
        $pm->start() and next;
        my %results = process_batch($b);
        $pm->finish( 0, \%results );
    }
    $pm->wait_all_children;
    #print them
    @lines_to_print =
      sort { VcfReader::by_first_last_line( $a, $b, \%contigs ) }
      @lines_to_print;
    if (@lines_to_print) {
        my $incr_per_batch = @lines_to_process / @lines_to_print;
        foreach my $batch (@lines_to_print) {
            my $incr_per_line = $incr_per_batch / @$batch;
            foreach my $l (@$batch) {
                if ( ref $l eq 'ARRAY' ) {
                    print $OUT join( "\t", @$l ) . "\n";
                }
                else {
                    print $OUT "$l\n";
                }
                if ($progressbar) {
                    $n += $incr_per_line;
                    $next_update = $progressbar->update($n)
                      if $n >= $next_update;
                }
            }
        }
    }
    else {
        if ($progressbar) {
            $n += @lines_to_process;
            $next_update = $progressbar->update($n) if $n >= $next_update;
        }
    }
}

################################################
sub process_batch {

    #filter a set of lines
    my ($batch) = @_;
    my %results = ( batch_size => scalar(@$batch) );
    foreach my $line ( @{$batch} ) {
        chomp $line;
        my @split = split("\t", $line);
        my $l = convertCallsToInfo(\@split);
        push @{ $results{keep} }, $l;
    }
    return %results;
}

################################################

sub convertCallsToInfo{
    my ($line) = @_;
    #debug
    if (ref $line ne 'ARRAY'){
        die "$line is not an array ref ";
    }
    #debug
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
            if (not $af){
                $info .= ";AF=0";
            }
            if (not $ac){
                $info .= ";AC=0";
            }
        }else{
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
                        $f = sprintf("%g", $f);
                    };
                    push @freqs, $f;
                }
                $info .= ";AF=" . join(",", @freqs);
            } 
        } 
        $line = VcfReader::replaceVariantField($line, 'INFO', $info);
    }
    my @pgts = VcfReader::getAllPossibleGenotypeCodes($line);    
    $line = VcfReader::addVariantInfoField
        (
            line => $line,
            id => "PGTS",
            value => join(",", @pgts),
        );
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
    $line = VcfReader::addVariantInfoField
        (
            line => $line,
            id => "GTC",
            value => join(",", @gtcs),
        );
    if (defined $opts{keep_calls}){
        return $line;
    }else{
        my @ar = @$line[0..7];
        return \@ar;
    }
}



=head1 NAME

sampleCallsToInfo - remove sample calls from a VCF and write genotype and allele metrics to INFO fields instead

=head1 SYNOPSIS

sampleCallsToInfo -i <input vcf file> [options]

sampleCallsToInfo -h (show help message)

sampleCallsToInfo -m (show manual page)


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

=item B<--forks>

Number of forks to create for parallelising your analysis. By default no forking is done. To speed up your analysis you may specify the number of parallel processes to use here. (N.B. forking only occurs if a value of 2 or more is given here as creating 1 fork only results in increased overhead with no performance benefit).

=item B<--cache>

Cache size. Variants are processed in batches to allow for efficient parallelisation. When forks are used the default is to process up to 10,000 variants at once or 1,000 x no. forks if more than 10 forks are used. If you find this program comsumes too much memory when forking you may want to set a lower number here. When using forks you may get improved performance by specifying a higher cache size, however the increase in memory usage is proportional to your cache size multiplied by the number of forks.

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
