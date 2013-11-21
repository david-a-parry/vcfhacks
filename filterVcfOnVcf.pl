#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;
use Term::ProgressBar;
use Data::Dumper;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use FindBin;
use lib "$FindBin::Bin";
use ParseVCF;


=head1 NAME

filterVcfOnVcf.pl - filter a VCF file using variants from one or more other VCF files.

=head1 SYNOPSIS

filterVcfOnVcf.pl -i <input vcf file> -f <vcfs to use to filter input> [options]
filterVcfOnVcf.pl -i <input vcf file> -d <directory containing vcfs to use to filter input> [options]
filterVcfOnVcf.pl -h (show help message)
filterVcfOnVcf.pl -m (show manual page)

=head1 ARGUMENTS

=over 8 

=item B<-i    --input>

Input VCF file.

=item B<-o    --output>

Output file. Optional - defaults to STDOUT.

=item B<-f    --filter>

One or more VCF files to use as to filter input. Variants that match from the input VCF will be filtered from the output.

=item B<-d    --directories>

One or more directories containing VCF files to use to filter input. Only files with a ".vcf" or ".vcf.gz" extension will be used.

=item B<-x    --expression>

Perl style regular expression to use to match VCF files when using --directories argument.

=item B<-s    --samples>

Samples from input file to check. If specified only variant alleles from these samples will be used to compare against variants in filter VCFs. Default is to look at all alleles.

=item B<-r    --reject>

Samples from filter VCF files to use for filtering.  If specified only variant alleles from these samples will be used to compare with the input VCF. Default is to look at all alleles.

=item B<-t    --threshold>

Minimum number of samples containing at least one matching variant allele before a variant is considered a match (optional).

=item B<-q    --quality>

Minimum variant quality as defined by the QUAL field.  Variants in the input will only be printed if they have an equal or higher variant quality score. Variants in filter VCFs will only be used for filtering if they have an equal or higher variant quality score.

=item B<-g    --genotype_quality>

Minimum genotype quality.  Only applicable if --samples,  --reject or --threshold arguments are used. Sample alleles will only be counted if they have a genotype quality score equal to or higher than this value.

=item B<-p    --print_matching>

Use this flag to reverse the script so that only matching lines are printed.

=item B<-b    --bar>

Show a progress bar.

=item B<-h    --help>

Show help message.

=item B<-m    --manual>

Show manual page.

=back

=head1 DESCRIPTION

Filter variants from a VCF using one or more other VCFs.  You may specify samples or a number of samples with matching variants to filter with.  You can also use minimum genotype or variant qualities to filter with. 

=head1 AUTHOR

David A. Parry
University of Leeds


=head1 COPYRIGHT AND LICENSE

Copyright 2013  David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut



my $vcf;
my @filter_vcfs;
my @dirs;
my @samples;#samples in vcf input to check calls for - filter only if they contain same allele as in filter_vcf. Default is to simply check alleles in ALT field and ignore samples.
my @reject;#if specified will only check alleles for these samples in filter_vcfs
my $threshold = 0;#only filter if we see the allele this many times in filter_vcfs
my $print_matching = 0;#flag telling the script to invert so we print matching lines and filter non-matching lines
my $out;
my $min_qual = 0;
my $minGQ = 0;
my $help;
my $man;
my $progress;
my $regex; #match this regex if looking in dir
GetOptions('x|expression=s' => \$regex, 'input=s' => \$vcf, 'output=s' => \$out, 'filter=s{,}' => \@filter_vcfs, 'directories=s{,}' => \@dirs, 'samples=s{,}' => \@samples, 'reject=s{,}' => \@reject, 'threshold=i' => \$threshold, 'quality=f' => \$min_qual,'genotype_quality=f' => \$minGQ, 'bar' => \$progress, 'help' => \$help, 'manual' => \$man, 'print_matching' => \$print_matching) or pod2usage(-message => "Syntax error.", exitval => 2);

pod2usage(-verbose => 2) if $man;
pod2usage(-verbose => 1) if $help;
pod2usage(-message => "Syntax error", exitval => 2) if (not $vcf or (not @filter_vcfs and not @dirs));

push @filter_vcfs, get_vcfs_from_directories(\@dirs, $regex) if @dirs;
die "No VCFs found to use as filters.\n" if not @filter_vcfs;
print STDERR "Initializing input VCF\n";
my $vcf_obj = ParseVCF->new(file=> $vcf);
my @filter_objs = ();
my $current_f = 0;
foreach my $filter(@filter_vcfs){
	$current_f++;
	my $filter_obj;
	print STDERR "Initializing filter VCF $current_f of " . scalar@filter_vcfs . " ($filter)\n";
	if ($filter =~ /\.gz$/){
        $filter_obj = ParseVCF->new(file=> $filter, noLineCount => 1);
    }else{
	    $filter_obj = ParseVCF->new(file=> $filter);
    }
	push @filter_objs, $filter_obj;
}
print STDERR "Done intializing filter VCFs.\n";
my $OUT;
if ($out){
	open ($OUT, ">$out") or die "Can't open $out for writing: $!\n";
}else{
	$OUT = \*STDOUT;
}
my $prev_chrom = 0;
my $progressbar;
if ($progress){
	$progressbar = Term::ProgressBar->new({name => "Filtering", count => $vcf_obj->countLines("variants"), ETA => "linear", });
}
my $next_update = 0;
my $n = 0;
print $OUT  $vcf_obj->getHeader(0);
print $OUT "##filterVcfOnVcf.pl=\"filter_vcfs=" .join(",", @filter_vcfs) ." genotype_quality=$minGQ quality=$min_qual samples=@samples reject=@reject threshold=$threshold print_matching=$print_matching\n";
print $OUT  $vcf_obj->getHeader(1);
my $lines_filtered = 0;
my $lines_passed = 0;
LINE: while (my $line = $vcf_obj->readLine){
	$n++;
	if ($progress){
                $next_update = $progressbar->update($n) if $n >= $next_update;
        }
	my $qual = $vcf_obj->getVariantField('QUAL');
	my $chrom = $vcf_obj->getVariantField('CHROM');
	my $pos = $vcf_obj->getVariantField('POS');
	my $ref = $vcf_obj->getVariantField('REF');
	next LINE if $qual < $min_qual;
	my @alts = ();
	if (@samples){
		@alts = $vcf_obj->getSampleActualGenotypes(multiple => \@samples, return_alleles_only => 1);
	}else{
		@alts = $vcf_obj->readAlleles(alt_alleles=>1);
	}
	my $threshold_counts = 0;#count samples in all filter_vcfs i.e. don't reset after each file
FILTER:	foreach my $filter_obj(@filter_objs){
		if ($filter_obj->searchForPosition(chrom=> $chrom, pos => $pos)){
FILTER_LINE:	while (my $filter_line = $filter_obj->readPosition()){
				my $filter_chrom = $filter_obj->getVariantField('CHROM');
				my $filter_pos = $filter_obj->getVariantField('POS');
				my $filter_qual = $filter_obj->getVariantField('QUAL');
				next FILTER_LINE if $filter_qual < $min_qual;
				my $f_ref = $filter_obj->getVariantField('REF');
				next FILTER_LINE if ($chrom ne $filter_chrom or $pos ne $filter_pos); #this might happen if tabix finds a deletion
				next FILTER_LINE if $ref ne $f_ref;
				my @f_alts = ();
				if (@reject){
					@f_alts = $filter_obj->getSampleActualGenotypes(multiple => \@reject, return_alleles_only => 1);
				}else{
					@f_alts = $filter_obj->readAlleles(alt_alleles=>1);
				}
				my $alts_match = 0;
				foreach my $alt (@alts){
					$alts_match++ if grep{/^$alt$/i} @f_alts;
					last if $alts_match;
				}
				next FILTER_LINE if not $alts_match;
				if ($threshold){
					my @t_samples = ();
					if (@reject){
						@t_samples = @reject;
					}else{
						@t_samples = $filter_obj->getSampleNames();
					}
					foreach my $t_samp(@t_samples){
						my @t_alleles = $filter_obj->getSampleActualGenotypes(sample => $t_samp, return_alleles_only => 1);
						foreach my $t_allele (@t_alleles){
							if (grep {/^$t_allele$/} @alts){
							#we're counting no. of samples with at least one matching allele
							#not no. matching alleles
								$threshold_counts++;
								last;
							}
						}
					}
					next FILTER_LINE if $threshold_counts < $threshold;
				}
				#we have a match - no need to look at other filter_vcfs
				$lines_filtered++;
				if ($print_matching){
					print $OUT "$line\n";
					next LINE;
				}else{
					next LINE;
				}
			}
		}
	}#no matching line
	$lines_passed++;
	print $OUT "$line\n" if not $print_matching; 
}
if ($progressbar){
        $progressbar->update($vcf_obj->countLines("variants")) if $vcf_obj->countLines("variants") >= $next_update;
}
if ($print_matching){
	print STDERR "$lines_filtered matching variants printed, $lines_passed filtered ";
}else{
	print STDERR "$lines_filtered matching variants filtered, $lines_passed printed ";
}
print STDERR "(" .$vcf_obj->countLines("variants") . " total)\n";

##################
sub get_vcfs_from_directories{
	my ($dirs, $regex) = @_;
	my @vcfs = ();
	foreach my $d (@$dirs){
		$d =~ s/\/$//;
		opendir (my $DIR, $d) or die "Can't open directory $d: $!\n";
		#my @dir_vcfs = grep {/\.vcf(\.gz)*$/i} readdir $DIR;
		my @dir_vcfs = grep {/\.vcf$/i} readdir $DIR;#can't use gzipped vcfs with ParseVCF search functions
		if (not @dir_vcfs){
			print STDERR "WARNING - no VCF files in directory $d\n";
		}elsif ($regex){
			@dir_vcfs = grep {/$regex/} @dir_vcfs;
			if (not @dir_vcfs){
				print STDERR "WARNING - no VCF files matching regex /$regex/ in directory $d\n" if not @dir_vcfs;
			}else{
				foreach my $v (@dir_vcfs){
					push @vcfs, "$d/$v";
				}
			}
		}else{
			foreach my $v (@dir_vcfs){
				push @vcfs, "$d/$v";
			}
		}
	}
	return @vcfs;
}
