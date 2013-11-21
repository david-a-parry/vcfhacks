#!/usr/bin/perl
#David Parry August 2011
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
my $vcf;
my $out;
my @samples;
my $check_presence_only;
my $ignore_non_existing;#don't exit if a sample is not found
my %sm_index;# this hash to index position of each sample in vcf file
my %rej_index;
my @reject = ();
my $threshold;
my $only_sample_alleles;
my $remove_calls;
my $help;
my $manual;

GetOptions('existing' => \$ignore_non_existing, 'input=s' =>\$vcf, "output=s" => \$out, 'k|samples=s{,}' => \@samples, 'reject=s{,}' => \@reject, 'threshold=i' => \$threshold, 'presence' => \$check_presence_only, 'x|sample_alleles_only' => \$only_sample_alleles, 'calls' => \$remove_calls, 'help' => \$help, "manual" => \$manual) or pod2usage(-message=> "syntax error.\n");
pod2usage(-verbose => 2) if $manual;
pod2usage(-verbose => 1) if $help;
pod2usage(-message=> "syntax error - required argument not found.\n") if (not $vcf or not @samples);
my $IN;
if ($vcf =~ /\.gz/){
	open ($IN, "gzip -dc $vcf |") || die "Can't open gzipped $vcf\n";
}else{
	open ($IN, $vcf) || die "Can't open $vcf\n";
}
my $OUT;
if ($out){
	open ($OUT, ">$out") || die "Can't open $out for writing: $!\n";
}else{
	$OUT = \*STDOUT;
}
my $format = 0;
my $ref = 0;
my $alt = 0;
my @sample_order;#for use if removing other sample calls
LINE: while (<$IN>){
	chomp;
	if (/^##/){
		print $OUT "$_\n";
		next;
	}
	if (/^#[A-Z]/){
		my @header = split("\t");
		foreach my $sample (@samples){
			my $n = 0;	
			no warnings 'uninitialized';
			$n++ until $header[$n] eq $sample or $n > $#header;
			if ($n > $#header){
				if (not $ignore_non_existing){
					die "Can't find sample $sample in header line.\nHeader:\n$_\n";
				}else{
					print STDERR "Warning - can't find sample $sample in header line.\n";
				}
			}else{
				$sm_index{$sample} = $n;
			}
		}
		if (@reject){
			foreach my $reject (@reject){
				my $n = 0;
				++$n until $header[$n] eq $reject or $n > $#header;
				if ($n > $#header){
					if (not $ignore_non_existing){
						die "Can't find sample $reject in header line.\nHeader:\n$_\n";
					}else{
						print STDERR "Warning - can't find sample $reject in header line.\n";
					}
				}else{
					$rej_index{$reject} = $n;
				}
			}
		}
		$ref++ until $header[$ref] eq "REF" or $ref > $#header;
		die "Can't find REF field in header line.\nHeader:\n$_\n" if ($ref > $#header);
		$alt++ until $header[$alt] eq "ALT" or $alt > $#header;
		die "Can't find ALT field in header line.\nHeader:\n$_\n" if ($alt > $#header);
		$format++ until $header[$format] eq "FORMAT" or $format > $#header;
		die "Can't find FORMAT field in header line.\nHeader:\n$_\n" if ($format > $#header);
		if ($remove_calls){
			@sample_order = ();#clear this array in case if we've stuck VCFs together and have come across a new header
			print $OUT join("\t", @header[0..$format]) . "\t";
			foreach my $k (sort {$sm_index{$a} <=> $sm_index{$b}} keys %sm_index){#if removing calls let's retain the original order
				push (@sample_order, $k);
			}
			print $OUT join("\t", @sample_order) ."\n";
		}else{
			print $OUT "$_\n"; 
		}
		next;
	}
	my @line = split("\t");
	my @form = split(":", $line[$format]);
	my $skip = 0;
	my $gt = 0;
	$gt++ until $form[$gt] eq "GT" or $gt > $#form;
	die "Can't find genotype field in variant line:\n$_\n" if ($gt > $#form);
	my %alleles;
	my $sample_no_calls = 0;
	my $sample_refs = 0;
SAMPLE:	foreach my $sample (keys %sm_index){
		my @call = split(":", $line[$sm_index{$sample}]); 
		if ($call[$gt] =~ /([\d])[\/\|]([\d])/){
			if ($1 == 0 and $2 == 0){
				if ($check_presence_only){
					next SAMPLE;
				}else{
					next LINE;#by default only keep variants present in all samples
				}
			}else{
				push (@{$alleles{$sample}}, $1, $2);
				#samples only get added to %alleles hash if they are called and not reference
				#because we compare the samples in %alleles hash only this allows variants 
				#with no call to go through
			}
		}else{
			next LINE unless $check_presence_only;
		}
	}
	next LINE if (keys%alleles < 1);#i.e. if only reference (0) or no calls were found
	my %reject_alleles;	
	if (@reject){
		foreach my $reject (keys %rej_index){
			my @call = split(":", $line[$rej_index{$reject}]); 
			if ($call[$gt] =~ /([\d])[\/\|]([\d])/){
				$reject_alleles{$1}++; #store alleles from rejection samples as keys of %reject_alleles
				$reject_alleles{$2}++;
			}
		}
	}
	my %var_call;
	my %breached;
	my %count;
	

	foreach my $samp (keys%alleles){
	#if we're looking for alleles that match in ALL samples than we only need to check a single hash entry
ALLELE:		foreach my $allele (@{$alleles{$samp}}){
			next ALLELE if ($allele == 0);
			next ALLELE if (exists $reject_alleles{$allele});
			$count{$allele}++; #count number of unrejected alleles for comparison with threshold by storing as key of %count hash
			if ($threshold){
                my $all_call_string = join("\t", @line[($format+1)..$#line]);
				my $t = 0;
				$t++ while $all_call_string =~ /(($allele)[\/\|]\d|\d[\/\|]$allele)/g;#searching all calls 
				if ($t > $threshold){
					$breached{$allele}++;
					next ALLELE;
				}
			}
			my $allele_matches = 0;
			foreach my $sample (keys %alleles){
				$allele_matches++ if (grep (/$allele/, @{$alleles{$sample}})); #we're counting the no. of samples with matching allele, not no. of occcurences (zygosity) of allele
			}
			if ($check_presence_only){#don't worry if all samples have variant or not if checking presence only
				$var_call{$allele}++ ;
			}else{
				$var_call{$allele}++ if ($allele_matches == keys%alleles); #i.e. if all of our calleda '--keep' sample genotypes match this allele in either het or homo state
			}
		}
	}
	next LINE if (keys %var_call < 1);#if we don't have at least one valid variant allele
	if (keys %count and keys %breached){
		next LINE if (keys %count == keys %breached);
	}
	if ($only_sample_alleles){
		my @alts = split(",", $line[$alt]);
		unshift @alts, $line[$ref];#add reference call to beginning of @alts array
		my @alt_ref=();
		if (@alts > 2){#i.e. more than one alt call
			foreach my $sample (keys %sm_index){
				foreach my $call (@{$alleles{$sample}}){
					push (@alt_ref, $alts[$call]) unless $call == 0;
				}
			}
			my %seen = ();
			@alt_ref = grep { !$seen{$_}++ } @alt_ref;#remove dups
			foreach my $sample (keys %sm_index){
				my @call = split(":", $line[$sm_index{$sample}]); 
				if ($call[$gt] =~ /([\d])[\/\|]([\d])/){
					#need to adjust allele numbering if only keeping sample alleles
					my @new_call = (); 
					foreach my $call ($1, $2){
						if ($call == 0){
							push (@new_call, 0);
						}else{
							my $actual_allele = $alts[$call];
							my $i = 0;
							++$i until $alt_ref[$i] eq $actual_allele or $i >$#alt_ref;
							die "Error identifying alternate alleles for line:\n$_" if $i > $#alt_ref;
							push (@new_call, $i+1);
						}
					}
					splice (@call, $gt, 1, join("/", @new_call));#replace genotype call
					splice (@line, $sm_index{$sample}, 1, join(":", @call));#replace sample genotype field
				}
			}
			splice (@line, $alt, 1, join(",", @alt_ref));#replace ALT field
		}
	}
	if ($remove_calls){
		print $OUT join("\t", @line[0..$format]) ."\t";
		my @rest_of_line = ();
		foreach my $sample (@sample_order){
			push (@rest_of_line, $line[$sm_index{$sample}]);
		}
		print $OUT join("\t", @rest_of_line) ."\n";
		
	}else{
		print $OUT "$_\n";
	}
}

=head1 NAME

filterOnSample.pl - filter variants in vcf that belong to specific samples.

=cut

=head1 SYNOPSIS

    filterOnSample.pl --input [var.vcf] --keep [samples to keep variants if present in all] --reject [samples to reject variants from if present in any] 

=head1 ARGUMENTS

=over 8 

=item B<-i    --input>

vcf file input.

=item B<-o    --output>

output filename.

=item B<-k    --samples>

IDs of samples to keep variants from. Variants will be kept only if present in ALL of these samples in either heterozygous or homozygous states unless --presence flag is set.  Samples must be entered as contained in the vcf header.

=item B<-p    --presence>

Use this flag to print variants present in any sample specified by the --keep option rather than variants present in all.

=item B<-r    --reject>

IDs of samples to reject variants from. Variants will be rejected if present in ANY of these samples.

=item B<-t    --threshold>

Reject variants present in more than this number of samples.

=item B<-x    --sample_alleles_only>

Use this flag to remove alternative alleles from other samples from output vcf.

=item B<-c    --calls>

Use this flag to remove calls from other samples from output vcf.

=item B<-e    --existing>

Use this flag to cause the program to ignore non-existing samples rather than exiting.

=item B<-h    --help>

Display help message.

=item B<-m    --manual>

Show manual page

=back
=cut

=head1 EXAMPLES

 filterOnSample.pl -i [var.vcf] -k Sample1 -r Sample2 Sample3
 (look for variants in Sample1 but reject variants also present in Sample2 or Sample3)
 
 filterOnSample.pl -i [var.vcf] -k Sample1 -r Sample2 Sample3 -t 4
 (same but also reject variants present in 4 or more samples)
 
 filterOnSample.pl -i [var.vcf] -k Sample1 -r Sample2 Sample3 -s 
 (same but don't print variant alleles from other samples for each variant if more than one variant allele exists)

=head1 DESCRIPTION


=cut

=head1 AUTHOR

David A. Parry

University of Leeds

=head1 COPYRIGHT AND LICENSE

Copyright 2013  David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.


=cut
