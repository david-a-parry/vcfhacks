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

One or more VCF files to use as to filter input. Variants that match from the input VCF will be filtered from the output. By default, if any sample matches a variant it will be filtered but combinations of --reject, --not_samples, --allele_frequency_filter, --threshold and --filter_homozygotes can be used to modify this behaviour.

=item B<-d    --directories>

One or more directories containing VCF files to use to filter input. Only files with a ".vcf" or ".vcf.gz" extension will be used.

=item B<-e    --expression>

Perl style regular expression to use to match VCF files when using --directories argument.

=item B<-s    --samples>

Samples from input file to check. If specified only variant alleles from these samples will be used to compare against variants in filter VCFs. Any line in your input without a variant call in these samples will be filtered. Default is to look at all alleles for each variant.

=item B<-r    --reject>

Samples from filter VCF files to use for filtering.  If specified only variant alleles from these samples will be used to compare with the input VCF. Default is to look at all alleles.

=item B<-x    --not_samples>

Samples from filter VCF files to ignore.  If specified variant alleles from all samples except these will be used to compare with the input VCF. Default is to look at all alleles.

=item B<-y    --allele_frequency_filter>

Reject variants if the allele frequency in all --filter VCFs is equal to or greater than this value. The value must be a float between 0 and 1. 

=item B<-t    --threshold>

Filter lines if the number of samples containing at least one matching variant allele is equal to or greater than this value.

=item B<-z    --filter_homozygotes>

Filter lines if a sample is homozygous for a variant. If --reject or --not_samples options are used only the relevant samples will be checked. If this option is used with --allele_frequency_filter or --threshold options lines will be filtered if the variant is homozygous in a sample regardless of allele frequency; if there is no homozygous sample lines will be filtered n frequency as normal. If used without --allele_frequency_filter or --threshold options lines will only be filtered if the variant is homozygous in at least one sample.

=item B<-q    --quality>

Minimum variant quality as defined by the QUAL field.  Variants in the input will only be printed if they have an equal or higher variant quality score. Variants in filter VCFs will only be used for filtering if they have an equal or higher variant quality score.

=item B<-g    --genotype_quality>

Minimum genotype quality.  Only applicable if --samples,  --reject or --threshold arguments are used. Sample alleles will only be counted if they have a genotype quality score equal to or higher than this value. Default is 0 (i.e. no filtering on genotype quality).

=item B<-a    --aff_quality>

Minimum genotype qualities to consider for samples specified by --samples argument only. Any sample call with a genotype quality below this threshold will be considered a no call. Default is 0.

=item B<-u    --un_quality>

Minimum genotype qualities to consider for samples specified by --reject argument only. Any sample call with a genotype quality below this threshold will be considered a no call. Default is 0.

=item B<-p    --print_matching>

Use this flag to reverse the script so that only matching lines are printed.

=item B<-b    --progress>

Show a progress bar.

=item B<-h    --help>

Show help message.

=item B<-m    --manual>

Show manual page.

=back

=head1 DESCRIPTION

Filter variants from a VCF using one or more other VCFs.  You may specify samples or a number of samples with matching variants to filter with.  You can also use minimum genotype or variant qualities to filter with. 

Each variant from the given --input file will be assessed and if all alternative/variant alleles for a given variant are represented in the --filter VCFs provided, the variant will be filtered. You may use --reject or --not_samples arguments to only filter variant alleles if present in specific samples in the --filter VCFs. You may also use the --samples argument to only compare variants from specific samples in your --input VCF and filter any variants that do not have alleles represented by these samples. See above for details of all available options. 

=head1 AUTHOR

David A. Parry
University of Leeds


=head1 COPYRIGHT AND LICENSE

Copyright 2013,2014  David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut



my $vcf;
my @filter_vcfs;
my @dirs;
my @samples;#samples in vcf input to check calls for - filter only if they contain same allele as in filter_vcf. Default is to simply check alleles in ALT field and ignore samples.
my @reject;#if specified will only check alleles for these samples in filter_vcfs
my @ignore_samples; #these samples will be ignored in either VCF 
my $threshold = 0;#only filter if we see the allele this many times in filter_vcfs
my $filter_homozygotes; #flag to filter if any of the reject samples are homozygous
my $maf = 0;
my $print_matching = 0;#flag telling the script to invert so we print matching lines and filter non-matching lines
my $out;
my $min_qual = 0;
my $minGQ = 0;
my $aff_quality ;#will convert to $minGQ value if not specified
my $unaff_quality ;#will convert to $minGQ value if not specified
my $help;
my $man;
my $progress;
my $regex; #match this regex if looking in dir
my %opts = (
        'not_samples' => \@ignore_samples,
        'expression' => \$regex,
        'input' => \$vcf,
        'output' => \$out,
        'filter' => \@filter_vcfs,
        'directories' => \@dirs,
        'samples=s' => \@samples,
        'reject' => \@reject,
        'allele_frequency_filter' => \$maf,
        'threshold' => \$threshold,
        'filter_homozygotes' => \$filter_homozygotes,
        'quality' => \$min_qual,
        'genotype_quality' => \$minGQ,
        'aff_quality' => \$aff_quality,
        'un_quality' => \$unaff_quality,
        'progress' => \$progress,
        'help' => \$help,
        'manual' => \$man,
        'print_matching' => \$print_matching
); 

GetOptions(
        \%opts,
        'x|not_samples=s{,}' => \@ignore_samples,
        'expression=s' => \$regex,
        'input=s' => \$vcf,
        'output=s' => \$out,
        'f|filter=s{,}' => \@filter_vcfs,
        'directories=s{,}' => \@dirs,
        'samples=s{,}' => \@samples,
        'reject=s{,}' => \@reject,
        'y|allele_frequency_filter=f' => \$maf,
        'threshold=i' => \$threshold,
        'z|filter_homozygotes' => \$filter_homozygotes,
        'quality=f' => \$min_qual,
        'genotype_quality=f' => \$minGQ,
        'a|aff_quality=i' => \$aff_quality,
        'un_quality=i' => \$unaff_quality,
        'bar|progress' => \$progress,
        'help' => \$help,
        'manual' => \$man,
        'p|print_matching' => \$print_matching) 
or pod2usage(
        -message => "Syntax error.",
        -exitval => 2
);

pod2usage(-verbose => 2) if $man;
pod2usage(-verbose => 1) if $help;
pod2usage(-message => "Syntax error", exitval => 2) if (not $vcf or (not @filter_vcfs and not @dirs));
pod2usage(-message => "Syntax error - cannot use --reject and --not_samples argument together", exitval => 2) if (@reject and @ignore_samples);
pod2usage(-message => "Variant quality scores must be 0 or greater.\n", -exitval => 2) if ($min_qual < 0 );
pod2usage(-message => "Genotype quality scores must be 0 or greater.\n", -exitval => 2) if ($minGQ < 0 );
pod2usage(-message => "--allele_frequency_filter (-y) value must be between 0 and 1.\n", -exitval => 2) if ($maf < 0 or $maf > 1);
if (defined $aff_quality){
    pod2usage(-message => "Genotype quality scores must be 0 or greater.\n", -exitval => 2) 
        if $aff_quality < 0;
}else{
    $aff_quality = $minGQ;
}
if (defined $unaff_quality){
    pod2usage(-message => "Genotype quality scores must be 0 or greater.\n", -exitval => 2) 
        if $unaff_quality < 0;
}else{
    $unaff_quality = $minGQ;
}
my %ignores = ();
if (@ignore_samples){
    %ignores = map {$_ => 1 } @ignore_samples;
}
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
    if (not $vcf_obj->get_inputIsStdin){
        $progressbar = Term::ProgressBar->new({name => "Filtering", count => $vcf_obj->countLines("variants"), ETA => "linear", });
    }else{
        print STDERR "Can't print a progress bar when input is from STDIN\n";
        $progress = 0;
    }
}
my $next_update = 0;
my $n = 0;
print $OUT  $vcf_obj->getHeader(0);
print $OUT "##filterVcfOnVcf.pl=\"";
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
my $lines_filtered = 0;
my $lines_passed = 0;
LINE: while (my $line = $vcf_obj->readLine){
    $n++;
    if ($progress){
        $next_update = $progressbar->update($n) if $n >= $next_update;
    }
    my $qual = $vcf_obj->getVariantField('QUAL');
    my $chrom = $vcf_obj->getVariantField('CHROM');
    next LINE if $qual < $min_qual;

    #process each allele separately in case we have MNVs/deletions that need to be simplified
    
    my %min_vars = $vcf_obj->minimizeAlleles();
    my %sample_alleles = ();
    if (@samples){
        %sample_alleles = map{$_ => undef}  $vcf_obj->getSampleCall(multiple => \@samples, return_alleles_only => 1, minGQ => $aff_quality);
        delete $sample_alleles{0};
        delete $sample_alleles{'.'};
        next LINE if not keys %sample_alleles;#filter if we don't have any variants in our samples 
    }else{
    #if no samples specified use all alleles for %sample_alleles
        %sample_alleles = map {$_ => undef} keys %min_vars;
    }
    my %sample_matches = ();#check each allele matches in all filter_vcfs but don't reset after each file
    my %thresh_counts = ();#count samples in all filter_vcfs i.e. don't reset after each file
    my %af_counts  ;#count allele occurences and total alleles to calculate allele frequency
    my %f_genos = ();#store genotypes as keys if we're using $filter_homozygotes
    
FILTER: foreach my $filter_obj(@filter_objs){
        my @temp_reject = ();
        if (@ignore_samples){
            @temp_reject = $filter_obj->getSampleNames() ;
            @temp_reject = grep {! $ignores{$_} } @temp_reject;
        }else{
            push @temp_reject, @reject;
        }
ALLELE:   foreach my $allele (keys %sample_alleles){
            if ($filter_obj->searchForPosition(chrom => $min_vars{$allele}->{CHROM}, pos => $min_vars{$allele}->{POS})){
                my %f_alts = ();
                #get genotype call codes (0, 1, 2 etc.) for filter samples
FILTER_LINE:    while (my $filter_line = $filter_obj->readPosition()){
                    my $filter_qual = $filter_obj->getVariantField('QUAL');
                    next FILTER_LINE if $filter_qual < $min_qual;
                    my %filter_min = $filter_obj->minimizeAlleles();
                    if (@temp_reject){
                        #%f_alts = map {$_ => undef} $filter_obj->getSampleActualGenotypes(multiple => \@temp_reject, return_alleles_only => 1, minGQ => $unaff_quality);
                        %f_alts = map{$_ => undef}  $filter_obj->getSampleCall(multiple => \@temp_reject, return_alleles_only => 1, minGQ => $unaff_quality);
                        if ($filter_homozygotes){
                            my %genos = $filter_obj->getSampleCall(multiple => \@temp_reject, minGQ => $unaff_quality);
                            %genos = reverse %genos;
                            foreach my $k (keys %genos){
                                $f_genos{$k}++;
                            }
                        }
                    }elsif(not @reject and not @ignore_samples){
                        #%f_alts = map {$_ => undef} $filter_obj->readAlleles(alt_alleles=>1, minGQ => $unaff_quality);
                        %f_alts = map{$_ => undef}  $filter_obj->getSampleCall(all => 1, return_alleles_only => 1, minGQ => $unaff_quality);
                        if ($filter_homozygotes){
                            my %genos = $filter_obj->getSampleCall(all => 1, minGQ => $unaff_quality);
                            %genos = reverse %genos;
                            foreach my $k (keys %genos){
                                $f_genos{$k}++;
                            }
                        }
                    }
                    my $filter_match = '';#if one of the filter's ALTs matches store the ALT allele code here
ALT:                foreach my $alt (keys %f_alts){
                        next if $alt eq '.';
                        next if $alt == 0;
                        next if $min_vars{$allele}->{POS} ne $filter_min{$alt}->{POS};
                        next if $min_vars{$allele}->{REF} ne $filter_min{$alt}->{REF};
                        next if $min_vars{$allele}->{ALT} ne $filter_min{$alt}->{ALT};
                        $min_vars{$allele}->{CHROM} =~ s/^chr//;
                        $filter_min{$alt}->{CHROM} =~ s/^chr//;
                        next if $min_vars{$allele}->{CHROM} ne $filter_min{$alt}->{CHROM};
                        if ($filter_homozygotes and not $threshold and not $maf){
                        #is using $filter_homozygotes on its own we only consider something a
                        #'match' if it's homozygous
                            if (exists $f_genos{"$alt/$alt"} or exists $f_genos{"$alt|$alt"}){
                                $filter_match = $alt;
                                $sample_matches{$allele}++;
                                last ALT;
                            }
                        }else{
                            $filter_match = $alt;
                            $sample_matches{$allele}++;
                            last ALT;
                        }
                    }
                    if (not $filter_match){
                        next FILTER_LINE;
                    }
                        
                    if ($threshold){
                        my @t_samples = ();
                        if (@temp_reject){
                            @t_samples = @temp_reject;
                        }elsif(not @reject and not @ignore_samples){
                            @t_samples = $filter_obj->getSampleNames();
                        }
                        foreach my $t_samp(@t_samples){
                            my %t_alleles = map { $_ => undef } $filter_obj->getSampleCall(sample => $t_samp, return_alleles_only => 1, minGQ => $unaff_quality);
                            if (exists $t_alleles{$filter_match}){
                                $thresh_counts{$allele}++;
                            }
                        }
                        #foreach my $alt (@alts){
                         #   next FILTER_LINE if not exists $thresh_counts{$alt};
                          #  next FILTER_LINE if $thresh_counts{$alt} < $threshold;
                        #}
                    }
                    if ($maf){
                        my %f_allele_counts = $filter_obj->countAlleles(minGQ => $unaff_quality);
                        foreach my $f_al (keys %f_allele_counts){
                            if ($f_al eq $filter_match){
                                $af_counts{$allele}->{counts} += $f_allele_counts{$f_al};
                            }
                            $af_counts{$allele}->{total} += $f_allele_counts{$f_al};
                        }
                    }
                }#read pos
            }#search
        }#foreach allele
        #done each allele - see if we've got enough data to filter this line
        #and can skip other filter vcfs for win
        next FILTER if keys %sample_alleles != keys %sample_matches;#can't skip - not all alleles accounted for
        my $homozygous_alleles = 0;
        if ($filter_homozygotes){
            foreach my $allele (keys %sample_alleles){
                if (exists $f_genos{"$allele/$allele"} or exists $f_genos{"$allele|$allele"}){
                    $homozygous_alleles++;
                }
            }
            if ($homozygous_alleles == keys %sample_alleles){
                $lines_filtered++;
                if ($print_matching){
                    print $OUT "$line\n";
                    next LINE;
                }else{
                    next LINE;
                }
            }
        }

        foreach my $allele (keys %sample_alleles){
            if ($threshold){
                next FILTER if $thresh_counts{$allele} < $threshold;   
            }
        }#if we haven't gone to next FILTER and 
        #are not filtering on allele frequency we can filter this line
        #no need to look at other filter_vcfs
        if (not $maf){
            $lines_filtered++;
            if ($print_matching){
                print $OUT "$line\n";
                next LINE;
            }else{
                next LINE;
            }
        }
    }#foreach filter vcf
    #no line meeting criteria for allele in any filter vcf
    #now check allele frequency if given
    if ($maf){
        my $alleles_over_maf = 0;
COUNTS: foreach my $allele (keys %sample_alleles){
            last COUNTS if not exists $af_counts{$allele};
            if (($af_counts{$allele}->{counts}/$af_counts{$allele}->{total}) >= $maf){
                $alleles_over_maf++ 
            }else{
                last COUNTS;
            }
        }
        if ($alleles_over_maf == keys %sample_alleles){
            $lines_filtered++;
            if ($print_matching){
                print $OUT "$line\n";
                next LINE;
            }else{
                next LINE;
            }
        }
    }
    $lines_passed++;
    print $OUT "$line\n" if not $print_matching; 
}#readline
if ($progressbar){
        $progressbar->update($vcf_obj->countLines("variants")) if $vcf_obj->countLines("variants") >= $next_update;
}
if ($print_matching){
    print STDERR "$lines_filtered matching variants printed, $lines_passed filtered ";
}else{
    print STDERR "$lines_filtered matching variants filtered, $lines_passed printed ";
}
if (not $vcf_obj->get_inputIsStdin){
    print STDERR "(" .$vcf_obj->countLines("variants") . " total)";
}
print STDERR "\n";

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
