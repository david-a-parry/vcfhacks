#!/usr/bin/perl

use warnings;
use strict; 
use POSIX qw/strftime/;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use Term::ProgressBar;
use Tabix;
use FindBin;
use lib "$FindBin::Bin/lib";
use ParseVCF;

#TO DO - allow multiple cadd files, do not annotate if exact variant not found
#and optionally output not found variants to separate file

my @cadd_files = ();
my %opts = (cadd_file => \@cadd_files);
GetOptions(
        \%opts,
        "input=s",
        "output=s",
        "cadd_file=s{,}",
        "not_found=s",
        "do_not_rank",
        "filter=f",
        "progress",
        "help",
        "manual",
        );

pod2usage(-verbose => 2) if $opts{manual};
pod2usage(-verbose => 1) if $opts{help};
pod2usage(-message => "--input argument is required", exitval => 2) if not $opts{input};
pod2usage(-message => "--cadd_file argument is required", exitval => 2) if not @{$opts{cadd_file}};
pod2usage(-message => "--filter argument cannot be less than 0", exitval => 2) if defined $opts{filter} && $opts{filter} < 0;

my $sortex;
if (not $opts{do_not_rank}){
    if (-s $opts{input} > (1024**2 * 50)){#use Sort::External if file size bigger than 50 MB
        eval "use Sort::External; 1" or print STDERR "Sort::External module was not found - ".
            "will attempt to sort in memory. For huge files it is recommended to install Sort::External.\n";
        if (not $@){
            my $scheme = sub {
                substr($Sort::External::b, 0, 6) <=> 
                substr($Sort::External::a, 0, 6) };#rev numeric sort CADD scores
            $sortex = Sort::External->new
            (
                mem_threshold => 1024**2 * 24, 
                sortsub => $scheme,
            );
        }
    }
}
my %cadd_iters = ();
foreach my $cadd_file (@{$opts{cadd_file}}){
    checkCaddFile($cadd_file);
    $cadd_iters{$cadd_file} = Tabix->new(-data =>  $cadd_file, -index => "$cadd_file.tbi");
}

my $vcf_obj = ParseVCF->new(file=> $opts{input});
my $total_vcf = $vcf_obj->countLines("variants") if defined $opts{progress};
print STDERR "$opts{input} has $total_vcf variants\n" if defined $opts{progress};
my $NOT_FOUND;
if ($opts{not_found}){
    open ($NOT_FOUND, ">$opts{not_found}") or die "Can't open $opts{not_found} for writing: $!\n";
}
my $OUT;
if ($opts{output}){
    open ($OUT, ">$opts{output}") or die "Can't open $opts{output} for writing: $!\n";
}else{
    $OUT = \*STDOUT;
}
print $OUT  $vcf_obj->getHeader(0);
print $OUT "##INFO=<ID=CaddPhredScore,Number=A,Type=Float,Description=\"Variant site CADD phred score. ".
    "Score provided for each allele and given as '.' if not found.\">\n";
print $OUT "##rankOnCaddScore.pl=\"";

my @opt_string = ();
foreach my $k (sort keys %opts){
    if (not ref $opts{$k}){
        push @opt_string, "$k=$opts{$k}";
    }elsif (ref $opts{$k} eq 'SCALAR'){
        if (defined ${$opts{$k}}){
            push @opt_string, "$k=${$opts{$k}}";
        }else{
            push @opt_string, "$k=undef";
        }
    }elsif (ref $opts{$k} eq 'ARRAY'){
        if (@{$opts{$k}}){
            push @opt_string, "$k=" .join(",", @{$opts{$k}});
        }else{
            push @opt_string, "$k=undef";
        }
    }
}
my @variants = ();
print $OUT join(" ", @opt_string) . "\"\n" .  $vcf_obj->getHeader(1);


my $time = strftime("%H:%M:%S", localtime);
print STDERR "CADD annotation commencing: $time\n";
my $filtered = 0;
my $not_found = 0;
my $n = 0;
my $progressbar;
my $next_update = 0;
if (defined $opts{progress} and $total_vcf){
    $progressbar = Term::ProgressBar->new({name => "Annotating", count => $total_vcf, ETA => "linear", });
}
while (my $line = $vcf_obj->readLine){
    $n++;
    my %min_vars = $vcf_obj->minimizeAlleles();
    #get score for each allele or '-' if it can't be found
    my @cadd_scores = findCaddScore(\%min_vars, \%cadd_iters);
    if (grep {/^\.$/} @cadd_scores){
        $not_found++;
        if ($opts{not_found}){
            my @alts = getAltsWithoutScore(\@cadd_scores, \%min_vars);
            foreach my $al (@alts){
                print $NOT_FOUND "$min_vars{$al}->{CHROM}\t$min_vars{$al}->{POS}".
                    "\t.\t$min_vars{$al}->{REF}\t$min_vars{$al}->{ALT}\n";
            }
        }
    }
    my $max = ( sort {$b <=> $a} 
                grep {! /^\.$/ } @cadd_scores)[0];
    #if there is no score set $max to -1 so it goes to bottom of the pile
    $max = -1 if not defined $max;
    if ($opts{filter}){
        if (not defined $max){
            $filtered++;
            next;
        }elsif($max < $opts{filter}){
            $filtered++;
            next;
        }
    }
    $line = $vcf_obj->addVariantInfoField
        (
            'CaddPhredScore', 
            join(",", @cadd_scores),
        );
    if ($opts{do_not_rank}){
        print $OUT "$line\n";
    }else{
        push @variants, sprintf("%-6s%s", $max, $line);
        if ($sortex){#no Sort::External , sort in memory
            if (@variants > 99999){
                $sortex->feed(@variants);
                @variants = ();
            }
        }
    }
    if (defined $opts{progress}){
           $next_update = $progressbar->update($n) if $n >= $next_update;
    }
}
if (defined $opts{progress} and $total_vcf){
    $progressbar->update($total_vcf) if $total_vcf >= $next_update;
}
print STDERR "\nDone annotating CADD scores.\n";
if (not $opts{do_not_rank}){
    $n = 0;
    $next_update = 0;
    if ($sortex){
        $sortex->feed(@variants) if @variants;
        $sortex->finish;
        print STDERR "Writing output...\n";
        if (defined $opts{progress} and $total_vcf){
            $progressbar = Term::ProgressBar->new({name => "Writing", count => $total_vcf, ETA => "linear", });
        }
        while ( defined( $_ = $sortex->fetch ) ) {
            print $OUT substr($_, 6) ."\n";
            $n++;
            if (defined $opts{progress}){
               $next_update = $progressbar->update($n) if $n >= $next_update;
            }
        }
    }else{
        print STDERR "Performing sort in memory...\n";
        @variants = sort{
            substr($b, 0, 6) <=> substr($a, 0, 6) 
        } @variants;#rev numeric sort CADD scores
        print STDERR "Writing output...\n";
        if (defined $opts{progress} and $total_vcf){
            $progressbar = Term::ProgressBar->new({name => "Writing", count => $total_vcf, ETA => "linear", });
        }
        foreach my $var (@variants){
            print $OUT substr($var, 6) ."\n";
            $n++;
            if (defined $opts{progress}){
               $next_update = $progressbar->update($n) if $n >= $next_update;
            }
        }
    }
    if (defined $opts{progress} and $total_vcf){
        $progressbar->update($total_vcf) if $total_vcf >= $next_update;
    }
}
$time = strftime("%H:%M:%S", localtime);
print STDERR "Time finished: $time\n";
if ($opts{filter}){
    print "$filtered variants filtered on CADD score.";
    if ($not_found){
        print " $not_found variants not found.";
    }
    print "\n";
}
##########################
sub findCaddScore{
    my ($vars, $tabix_iter) = @_;
    my @scores = ();#return score for each allele
ALLELE:    foreach my $al (sort {$a<=>$b} keys %{$vars}){
        foreach my $iter (keys %{$tabix_iter}){
            my $it = $tabix_iter->{$iter}->query
                (
                $vars->{$al}->{CHROM},
                $vars->{$al}->{POS} -1, 
                $vars->{$al}->{POS},
                );
            next if not defined $it->{_};#not found
            while (my $result = $tabix_iter->{$iter}->read($it)){
                chomp($result);
                my @res = split("\t", $result);
                next if $res[0] ne $vars->{$al}->{CHROM};
                my ($pos, $ref, $alt) = reduceRefAlt($res[1], $res[2], $res[3]);
                next if ($vars->{$al}->{POS} != $pos);
                next if ($vars->{$al}->{REF} ne $ref); #should we error here?
                next if ($vars->{$al}->{ALT} ne $alt); #diff alt allele
                push @scores, $res[5];
                next ALLELE;
            }
        }#not found in any of our tabix iterators
        push @scores, '.';
    }
    return @scores;
}

##########################
sub checkCaddFile{
    my ($file) = @_;
    my $index = "$file.tbi";
    if (not -e $index ){
        die "Can't find tabix index for $file - please ensure $file is bgzip ".
            "compressed and indexed with tabix.\n";
    }
    if ($file !~ /\.gz$/){
        die "CADD file $file does not have a '.gz' extension and is ".
            "presumably not bgzip compressed. Please ensure to use a ".
            "bgzip compressed CADD file.\n";
    }
    my $FH = new IO::Uncompress::Gunzip $file or 
      die "IO::Uncompress::Gunzip failed while opening $file ".
      "for reading: \n$GunzipError";
    while (my $line = <$FH>){
        if ($line =~ /^##/){
            next;
        }elsif ($line =~ /^#/){
            if ($line !~ /^#Chrom\tPos\tRef\tAlt\t\w+Score\tPHRED/i){
                die "Invalid header line:\n\n$line\n\n".
                "Expected to find a header as follows:\n\n".
                "#Chrom\tPos\tRef\tAlt\tRawScore\tPHRED\n";
            }else{
                return;
            }
        }else{
            last;
        }
    }
    die "No header found for CADD file $file!\n";
}

##########################
sub reduceRefAlt{
    #reduce a single ref/alt pair to their simplest representation
    my ($pos, $ref, $alt) = @_;
    if (length($ref) > 1 and length($alt) > 1){
        #can only reduce if both REF and ALT are longer than 1
        my @r = split('', $ref);
        my @al = split('', $alt);
        while ($r[-1] eq $al[-1] and @r > 1 and @al > 1){
            #remove identical suffixes
            pop @r;
            pop @al;
        }
        while ($r[0] eq $al[0] and @r > 1 and @al > 1){
            #remove identical prefixes
            #increment position accordingly
            shift @r;
            shift @al;
            $pos++;
        }
        $ref = join('', @r);
        $alt = join('', @al);
    }
    return ($pos, $ref, $alt);
}


##########################
sub getAltsWithoutScore{
    my ($cadd_scores, $vars) = @_;
    my @alts = ();
    for (my $i = 0; $i < @$cadd_scores; $i++ ){
        if ($cadd_scores->[$i] eq '.'){
            #allele 1 will be first (i.e. pos 0) in array etc. 
            #so increment by 1
            push @alts, $i + 1;
        }
    }
    return @alts;
}
=head1 NAME

rankOnCaddScore.pl - annotate, rank and/or filter variants using CADD PHRED scores.

=head1 SYNOPSIS

        rankOnCaddScore.pl -i [vcf file] -c [tabix indexed cadd file(s)] [options]
        rankOnCaddScore.pl -h (display help message)
        rankOnCaddScore.pl -m (display manual page)

=cut

=head1 ARGUMENTS

=over 8

=item B<-i    --input>

Input VCF file.

=item B<-c    --cadd_file>

No, not a murder-solving monk, but one or more bgzip compressed and tabix indexed file of CADD scores. CADD score files are available from http://cadd.gs.washington.edu/download and can also be generated by uploading your own data to http://cadd.gs.washington.edu/score.

=item B<-o    --output>

Output CADD score annotated/filtered file. Optional - default is STDOUT.

=item B<-n    --not_found>

Output file for variants that can't be found in CADD files.  This output can be uploaded to the CADD server (http://cadd.gs.washington.edu/score) for scoring. The output from CADD server can then be added to your CADD files used in analysing future data if desired.
        
=item B<-d    --do_not_rank>

Use this flag to annotate variants but not reorder your VCF.

=item B<-f    --filter>

CADD Phred score cut-off. Filter any variant that has a CADD score below this score. 

=item B<-p    --progress>

Use this flag to print % progress to STDERR.

=item B<-h    --help>

Display help message.

=item B<-m    --manual>

Show manual page.

=back 

=cut


=head1 EXAMPLES

  rankOnCaddScore.pl -i input.vcf -o cadd_ranked.vcf -c cadd_file1.tsv.gz cadd_file2.tsv.gz
  (rank input.vcf on CADD score and output to cadd_ranked.vcf)

  rankOnCaddScore.pl -i input.vcf -o cadd_ranked.vcf -c cadd_file1.tsv.gz cadd_file2.tsv.gz -f 15
  (rank input.vcf on CADD score and filter anything with a CADD PHRED score below 15)

  rankOnCaddScore.pl -i input.vcf -o cadd_annotated.vcf -c cadd_file1.tsv.gz cadd_file2.tsv.gz -d
  (annotate but don't rank)

  rankOnCaddScore.pl -i input.vcf -o cadd_annotated.vcf -c cadd_file1.tsv.gz cadd_file2.tsv.gz -d -f 15
  (annotate but don't rank and filter anything  with a CADD PHRED score below 15)

  rankOnCaddScore.pl -i input.vcf -o cadd_annotated.vcf -c cadd_file1.tsv.gz cadd_file2.tsv.gz -d -f 15 -n not_found.tsv
  (as above but outputting variants that don't have CADD scores to not_found.tsv)


=head1 DESCRIPTION

This program reads a VCF file and searches the provided CADD score files for matching variants. Upon finding a corresponding variant in a CADD file it annotates the variant INFO field with "CaddPhredScore=" and a score for each allele (or . when a score is not found for an allele). By default variants are ordered by CADD PHRED score (or the highest CADD PHRED score for a multiallelic site) but this can be avoided using the --do_not_sort option. You may also filter on CADD score, but be aware that this is likely to be somewhat arbitrary. Advice from the authors suggests a PHRED score of 15 may be a good value to filter from because it is the median value for all possible canonical splice site changes and non-synonymous variants.

CADD files for all possible SNVs in the human genome can be downloaded from http://cadd.gs.washington.edu/download as well as indels present in 1000 genomes samples and ESP samples. When variants are not found you can use the --not_found option to output these to a file which can then be uploaded to http://cadd.gs.washington.edu/score to score these variants.  These can then be added to your database for future analyses after tabix indexing the output (tabix -s 1 -b 2 -e 2 cadd_output.tsv.gz). I would suggest collecting all your various CADD score files in one folder so that you can pass '/data/cadd_score/*.gz' to the --cadd_file option to analyze all the files in that folder.

=cut

=head1 AUTHOR

David A. Parry

University of Leeds

=head1 COPYRIGHT AND LICENSE

Copyright 2014  David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut

