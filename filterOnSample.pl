#!/usr/bin/perl
#David Parry August 2011
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use Term::ProgressBar;
use FindBin;
use lib "$FindBin::Bin";
use ParseVCF;
my $vcf;
my $out;
my @samples;
my $check_presence_only;
my $ignore_non_existing;#don't exit if a sample is not found
my %sm_index;# this hash to index position of each sample in vcf file
my %rej_index;
my @reject = ();#reject if allele is present in these samples
my @reject_except = (); #reject all except these samples
my $threshold;
my $quality = 20;
my $num_matching;
my $help;
my $manual;
my $progress;
my %opts = ('existing' => \$ignore_non_existing,
        'input' =>\$vcf,
        "output" => \$out,
        'samples' => \@samples,
        'reject' => \@reject,
        'reject_all_except' => \@reject_except,
        'threshold' => \$threshold,
        'presence' => \$check_presence_only,
        'quality' => \$quality,
        'num_matching' => \$num_matching,
        'help' => \$help,
        "manual" => \$manual,
        'progress' => \$progress);

GetOptions(\%opts,
        'existing' => \$ignore_non_existing,
        'input=s' =>\$vcf,
        "output=s" => \$out,
        'samples=s{,}' => \@samples,
        'r|reject=s{,}' => \@reject,
        'x|reject_all_except:s{,}' => \@reject_except,
        'threshold=i' => \$threshold,
        'p|presence' => \$check_presence_only,
        'quality=i' => \$quality,
        'num_matching' => \$num_matching,
        'help' => \$help,
        "manual" => \$manual,
        'b|progress' => \$progress) or pod2usage(-message=> "syntax error.\n");
pod2usage(-verbose => 2) if $manual;
pod2usage(-verbose => 1) if $help;
pod2usage(-message=> "syntax error: --input (-i) argument is required.\n") if not $vcf;
pod2usage(-message=> "syntax error: you must specify samples using at least one of the arguments --samples (-s), --reject (-r) or --reject_all_except (-x).\n") if not @samples and not @reject and not @reject_except;
print STDERR "Warning - --num_matching has no effect when --presence flag is set.\n" if $check_presence_only and $num_matching;

print STDERR "Initializing VCF input ($vcf)\n";
my $vcf_obj = ParseVCF->new(file=> $vcf);
my @not_found = ();
my @samples_checked = ();
my @reject_checked = ();
   
my ($samples_found, $samples_not_found) = check_samples(\@samples);
my ($reject_found, $reject_not_found) = check_samples(\@reject);
if (@$samples_not_found or @$reject_not_found){
    my @not_found = (@$samples_not_found, @$reject_not_found);
    if (not $ignore_non_existing){
        print STDERR "Warning - could not find the following samples in VCF:\n".join("\n", @not_found)."\n";
        if (@samples and not @$samples_found){
            print STDERR "Warning - no samples specified by --samples identified in VCF.\n";
        }
        if (@reject and not @$reject_found){
            print STDERR "Warning - no samples specified by --reject identified in VCF.\n";
        }
        @samples = @$samples_found;
        @reject = @$reject_found;
    }else{
        die "Could not find the following samples in VCF:\n".join("\n", @not_found)."\n";
    }
}
if (@reject_except){
    my @all = $vcf_obj->getSampleNames();
    push @reject_except, @samples; 
    my %subtract = map {$_ => undef} @reject_except;
    @all = grep {!exists $subtract{$_} } @all;
    push @reject, @all;
    my %seen = ();
    @reject = grep { ! $seen{$_}++} @reject;
}

if (not @reject and not @samples){
    print STDERR "Warning - no samples from --samples (-s), --reject (-r) or --reject_all_except (-x) argument found to filter. Your output will remain unchanged.\n";
}

my $OUT;
if ($out){
    open ($OUT, ">$out") || die "Can't open $out for writing: $!\n";
}else{
    $OUT = \*STDOUT;
}
my $progressbar;
if ($progress){
    if ($vcf eq "-"){
        print STDERR "Can't use --progress option when input is from STDIN\n";
        $progress = 0;
    }else{
        $progressbar = Term::ProgressBar->new({name => "Filtering", count => $vcf_obj->countLines("variants"), ETA => "linear", });
    }
}
my $next_update = 0;
print $OUT  $vcf_obj->getHeader(0) ."##filterOnSample.pl\"";
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
print $OUT join(" ", @opt_string) . "\"\n" .  $vcf_obj->getHeader(1);
my $line_count = 0;
LINE: while (my $line = $vcf_obj->readLine){
    $line_count++;
    if ($progress){
        $next_update = $progressbar->update($line_count) if $line_count >= $next_update;
    }
    my %alleles = ();
    #do samples first for efficiency (if they don't have a variant allele)
    if (@samples){
SAMPLE: foreach my $sample (@samples){
            my $call = $vcf_obj->getSampleCall(sample => $sample, minGQ => $quality);
            if ($call =~ /(\d+)[\/\|](\d+)/){
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
    }#otherwise we'll collect --reject alleles and see if there are any alts that aren't in %reject_alleles

    my %reject_alleles;    
    if (@reject){
        foreach my $reject (@reject){
            my $call = $vcf_obj->getSampleCall(sample => $reject, minGQ => $quality);
            if ($call =~ /(\d+)[\/\|](\d+)/){
                $reject_alleles{$1}++; #store alleles from rejection samples as keys of %reject_alleles
                $reject_alleles{$2}++;
            }
        }
    }
    if (not @samples){
        my @ref_alt = $vcf_obj->readAlleles();
        for (my $i = 1; $i < @ref_alt; $i++){
            push @{$alleles{alt}}, $i if not exists $reject_alleles{$i};
        }
        next LINE if (keys%alleles < 1);#i.e. if only reference (0) or no calls were found
    }
    my %var_call;
    my %breached;
    my %count;
    my %genotypes = $vcf_obj->countGenotypes();

    foreach my $samp (keys%alleles){
    #if we're looking for alleles that match in ALL samples than we only need to check a single hash entry
ALLELE: foreach my $allele (@{$alleles{$samp}}){
            next ALLELE if ($allele == 0);
            next ALLELE if (exists $reject_alleles{$allele});
            $count{$allele}++; #count number of unrejected alleles for comparison with threshold by storing as key of %count hash
            if ($threshold){
                my $t = 0;
                foreach my $k (keys %genotypes){
                    $t += $genotypes{$k} if ($k =~ /(^($allele)[\/\|]\d+|\d+[\/\|]$allele)$/);
                }
                if ($t > $threshold){
                    $breached{$allele}++;
                    next ALLELE;
                }
            }
            my $allele_matches = 0;
            foreach my $sample (keys %alleles){
                $allele_matches++ if (grep (/^$allele$/, @{$alleles{$sample}})); #we're counting the no. of samples with matching allele, not no. of occcurences (zygosity) of allele
            }
            if ($check_presence_only){#don't worry if all samples have variant or not if checking presence only
                $var_call{$allele}++ ;
            }elsif($num_matching){
                $var_call{$allele}++ if $allele_matches >= $num_matching;
            }else{
                $var_call{$allele}++ if ($allele_matches == keys%alleles); #i.e. if all of our called '--keep' sample genotypes match this allele in either het or homo state
            }
        }
    }
    next LINE if (keys %var_call < 1 );#if we don't have at least one valid variant allele
    if (keys %count and keys %breached){
        next LINE if (keys %count == keys %breached);
    }
    print $OUT "$line\n";
}

if ($progressbar){
        $progressbar->update($vcf_obj->countLines("variants")) if $vcf_obj->countLines("variants") >= $next_update;
}

#################################################
sub check_samples{
    my ($sample_ref) = @_;
    my @found;
    my @not_found;
    foreach my $s (@$sample_ref){
        if (not $vcf_obj->checkSampleInVcf($s)){
            push @not_found, $s;
        }else{
            push @found, $s;
        }
    }
    return (\@found, \@not_found);
}
 
#################################################

=head1 NAME

filterOnSample.pl - filter variants in vcf that belong to specific samples.

=cut

=head1 SYNOPSIS

    filterOnSample.pl --input [var.vcf] --samples [samples to keep variants if present in all] --reject [samples to reject variants from if present in any] 

=head1 ARGUMENTS

=over 8 

=item B<-i    --input>

vcf file input.

=item B<-o    --output>

output filename.

=item B<-s    --samples>

IDs of samples to keep variants from. Variants will be kept only if present in ALL of these samples in either heterozygous or homozygous states unless --presence or --num_matching flags are set.  Samples must be entered as contained in the vcf header.

=item B<-p    --presence>

Use this flag to print variants present in any sample specified by the --samples option rather than variants present in all.

=item B<-n    --num_matching>

Use this flag to print variants present in at least this many samples rather than only variants present in all.

=item B<-r    --reject>

IDs of samples to reject variants from. Variants will be rejected if present in ANY of these samples.

=item B<-x    --reject_all_except>

Reject variants present in all samples except these. If used without an argument all samples in VCF that are not specified by --samples argument will be used to reject variants. If one or more samples are given as argument to this option then all samples in VCF that are not specified by --samples argument or this argument will be used to reject variants.

=item B<-t    --threshold>

Reject variants present in more than this number of samples in the VCF. Counts all samples in the VCF irrespective of whether they are specified by the --samples or any other argument.

=item B<-q    --quality>

Minimum phred-like genotype quality to consider.  All calls below this quality will be considered no calls. Default is 20.

=item B<-e    --existing>

Use this flag to cause the program to ignore non-existing samples rather than exiting.

=item B<-b    --progress>

Show a progress bar.

=item B<-h    --help>

Display help message.

=item B<-m    --manual>

Show manual page

=back
=cut

=head1 EXAMPLES

 filterOnSample.pl -i [var.vcf] -s Sample1
 (only print variants if Sample1's genotype has a variant allele)

 filterOnSample.pl -i [var.vcf] -k Sample1 -r Sample2 Sample3
 (look for variants in Sample1 but reject variants also present in Sample2 or Sample3)
 
 filterOnSample.pl -i [var.vcf] -k Sample1 -r Sample2 Sample3 -t 4
 (same but also reject variants present in 4 or more samples)
 
 filterOnSample.pl -i [var.vcf] -k Sample1 -x 
 (look for variants in Sample1 and reject variants if present in any other sample in the VCF)
 
 filterOnSample.pl -i [var.vcf] -k Sample1 -x Sample2 
 (look for variants in Sample1 and reject variants if present in any other sample in the VCF except for Sample2)
 

=head1 DESCRIPTION

This program reads a VCF file and filters variants depending on which samples contain variant. Samples to keep variants from can be specified using --samples (-s) and samples to reject variants from can be specified with --reject (-r). 

=cut

=head1 AUTHOR

David A. Parry

University of Leeds

=head1 COPYRIGHT AND LICENSE

Copyright 2013  David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.


=cut
