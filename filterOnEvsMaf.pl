#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;
use Term::ProgressBar;
use Data::Dumper;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use POSIX qw/strftime/;
use FindBin;
use lib "$FindBin::Bin";
use ParseVCF;
my @samples;
my @evs;
my $freq;
my $quiet;
my $strict;
my %opts = (samples => \@samples, esp_file => \@evs, freq => \$freq, quiet=>\$quiet, Strict => \$strict);

GetOptions(\%opts,
        "output=s",
        "input=s",
        "samples=s{,}",
        "esp_file=s{,}",
        "dir=s",
        "help",
        "manual",
        "freq=f",
        "quiet",
        "Progress",
        ) or pod2usage(-exitval => 2, -message => "Syntax error") ;
pod2usage (-verbose => 2) if $opts{manual};
pod2usage (-verbose => 1) if $opts{help};

pod2usage(-exitval => 2, -message => "Syntax error") if not $opts{input} or (not $opts{dir} and not @{$opts{esp_file}});
if (defined $freq){
        pod2usage(-exitval => 2, -message => "value for --freq argument must be greater than 0 and less than 50") if $freq > 50 or $freq <= 0;
}

my $OUT;
if ($opts{output}){
        open ($OUT, ">$opts{output}") || die "Can't open $opts{output} for writing: $!";
}else{
        $OUT = *STDOUT;
}
my $progressbar;
my $next_update = 0;
my $prev_percent = 0;
my $n = 0;
my $kept = 0; #variants not filtered
my $filtered = 0; #variants filtered
my $pathogenic_snps = 0;
my $found = 0; #variants that match a known SNP
my $lines = 0; 
my $time = strftime("%H:%M:%S", localtime);
print STDERR "[$time] Initializing input VCF\n";
my $vcf_obj = ParseVCF->new( file=> $opts{input});
my $total_vcf = $vcf_obj->countLines("variants") if defined $opts{Progress};
print STDERR "$opts{input} has $total_vcf variants\n" if defined $opts{Progress};
$time = strftime("%H:%M:%S", localtime);
my @snp_objs = ();
my @snp_headers;
if ($opts{dir}){
    opendir (my $DIR, $opts{dir}) or die "Can't read directory $opts{dir}: $!\n";
    my @files = grep {/\.vcf(\.gz)*$/}readdir($DIR);
    die "No vcf files found in $opts{dir}\n" if not @files;
    foreach my $file (@files){
        push @evs, "$opts{dir}/$file";
    }
}
for (my $i = 0; $i < @evs; $i++){
    print STDERR "[$time] Initializing $evs[$i] ESP VCF " .($i + 1) ." of " .scalar(@evs) ."\n";
    my $snp_obj;
    if ($evs[$i] =~ /\.gz$/){
        $snp_obj = ParseVCF->new( file => $evs[$i], noLineCount =>1);
    }else{
        $snp_obj = ParseVCF->new( file => $evs[$i], noLineCount =>0);
    }
    push @snp_objs, $snp_obj;
    if (not -e $snp_obj->get_index){
        $snp_obj->indexVcf();
    }
    #my $total_snp = $snp_obj->countLines("variants") ;
    #print STDERR "$dbsnp[$i] has $total_snp variants\n" ;
    push @snp_headers, $snp_obj->getHeader(0);
}

if ($freq){
    if (not grep {/##INFO=<ID=MAF/} @snp_headers ){
        die  "Can't find MAF fields in evs file headers.\n";
    }
}
print $OUT  $vcf_obj->getHeader(0) ."##filterOnEvsMaf.pl=\"";
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

$time = strftime("%H:%M:%S", localtime);
print STDERR "[$time] EVS filter starting\n";

if (defined $opts{Progress} and $total_vcf){
        if ($opts{build} || $freq){
            $progressbar = Term::ProgressBar->new({name => "Filtering", count => $total_vcf, ETA => "linear", });
        }else{
            $progressbar = Term::ProgressBar->new({name => "Annotating", count => $total_vcf, ETA => "linear", });
        }
        #$progressbar->minor(0);
}

VAR: while (my $line = $vcf_obj->readLine){
    $n++;
    my $snp_meets_filter_criteria = 0; 
    my $snp_should_not_be_filtered = 0; #flag in case pathogenic flag is set or somethinga and we don't want to filter this snp no matter what
    my $is_known_snp = 0;
    my $replaced_already = 0; #use this flag to allow appending of multiple IDs found in our reference file even if --replace option is in place.
    if (defined $opts{Progress}){
           $next_update = $progressbar->update($n) if $n >= $next_update;
    }
    my $qual = $vcf_obj->getVariantField('QUAL');
    my $chrom = $vcf_obj->getVariantField('CHROM');
    my $pos = $vcf_obj->getVariantField('POS');
    my $ref = $vcf_obj->getVariantField('REF');
    foreach my $snp_obj (@snp_objs){
        if ($snp_obj->searchForPosition(chrom=> $chrom, pos => $pos)){
            while (my $snp_line = $snp_obj->readPosition()){#need to edit the below so we traverse through potentially multiple SNPs and still append/filter as necessarry
                #check whether the snp line(s) match our variant
                if (checkSnpMatches($vcf_obj, $snp_obj, $opts{samples})){
                    $is_known_snp++;
                    #perform filtering if fiters are set
                    if ($freq){
                        my $snp_info = $snp_obj->getVariantField("INFO");
                        my @info = split(';', $snp_info);
                        foreach my $info (@info){
                            if ($info =~ /^MAF=(.*)/){
                                my @mafs = split(',', $1);
                                foreach my $maf (@mafs){
                                    if ($maf >= $freq){
                                        $filtered++;
                                        $found++ ;
                                        next VAR;
                                    }
                                }
                            }
                        }
                    }else{#default behaviour is to filter everything that matches
                        $filtered++;
                        $found++ ;
                        next VAR ;

                    }
                }
            }
        }
    }#not found
    print $OUT "$line\n";
    $kept++;
    $found++ if $is_known_snp;
    if (defined $opts{Progress}){
        $next_update = $progressbar->update($lines) if $lines >= $next_update;
    }
}
if (defined $opts{Progress}){
        $progressbar->update($total_vcf) if $total_vcf >= $next_update;
}
$time = strftime("%H:%M:%S", localtime);
print STDERR "Time finished: $time\n";
print STDERR "$found matching variants identified.\n";
print STDERR "$filtered variants filtered, $kept variants retained.\n";


################################################
#####################SUBS#######################
sub evaluate_snp{
    #returns two values - first is 1 if snp should be filtered and 0 if not, 
    #the second is 1 if shouldn't be filtered under any circumstance (at the moment only if
    #pathogenic flag is set and snp has a SCS or CLNSIG value of 4 or 5
    my ($build, $path, $freq, $vcf_obj, $snp_obj, $samp_ar) = @_;
    my $matches = checkSnpMatches($vcf_obj, $snp_obj, $samp_ar);
    if (not $matches){
        return (0, 0);
    }
    my $snp_info = $snp_obj->getVariantField("INFO");
    if ($path){
        if ($snp_info =~ /SCS=([\d])/){
#SCS=4 indicates probable-pathogenic, SCS=5 indicates pathogenic, print regardless of freq or build if --pathogenic option in use
            if ($1 == 4 or $1 ==  5){
                return (0, 1);
            }
        }else{
            #die "No SCS field in snp line $snp_line" if $strict; #the SCS field appears to be dropped from the GATK bundle now
            #print STDERR "Warning, no SCS field in snp line $snp_line\n" unless $quiet;
        }
        if ($snp_info =~ /CLNSIG=([\d])/){
            if ($1 == 4 or $1 ==  5){
                return (0, 1);
            }
        }
    }
    if ($build){
        if ($snp_info =~ /dbSNPBuildID=([\d]+)/){
            if ($build >= $1){
                return (1, 0);
            }
        }else{
            die "No dbSNPBuildID field in snp line " .$snp_obj->get_currentLine ."\n"  if $strict;
            print STDERR "Warning, no dbSNPBuildID field in snp line " . $snp_obj->get_currentLine  . "\n" unless $quiet;
        }
    }
    if ($freq){
        my @info = split(/\;/, $snp_info);
        if ($freq <= 0.05){
            #G5 = minor allele freq > 5 % in at least 1 pop
            #G5A = minor allele freq > 5 % in all pops
            if (grep{/^(G5|G5A)/}@info){ 
                return (1, 0);
            }
        }
        my (@af) = grep {/^AF=[\d]\.[\d]+/} @info;
        if (@af){
            foreach my $af (@af){
                $af =~ s/AF=//;
                if ($freq <= $af){
                    return (1, 0);
                }
            }
        }
        my (@gmaf) = grep {/^GMAF=[\d]\.[\d]+/} @info;
        if (@gmaf){
            foreach my $gmaf (@gmaf){
                $gmaf =~ s/GMAF=//;
                if ($freq <= $gmaf){
                    return (1, 0);
                }
            }
        }
    }
    return (0, 0);
}
#################################################
sub checkSnpMatches{
    my ($obj, $snp, $samples_ref) = @_;

    my $chrom = $obj->getVariantField("CHROM");
    my $pos = $obj->getVariantField("POS");
    my $snp_chrom = $snp->getVariantField("CHROM");
    my $snp_pos = $snp->getVariantField("POS");
    my $snp_ref = $snp->getVariantField("REF");
    #my $ref = $obj->getVariantField("REF");
    #my $ref = $obj->getVariantField("REF");
    #my @alts = split(",", $obj->getVariantField("ALT"));
    $chrom =~ s/^chr//;
    $snp_chrom =~ s/^chr//;
    return 0 if $chrom ne $snp_chrom;
    return 0 if $pos ne $snp_pos;
    my @snp_alleles = $snp->readAlleles();    
    my @alleles = ();
    if (@$samples_ref){
        @alleles = $obj->getSampleActualGenotypes(multiple => $samples_ref, return_alleles_only => 1);
    }else{
        @alleles = $obj->readAlleles();    
    }
    #if none of our samples has a variant here then 
    #we might as well skip if snp meets our filters
    #so let's say it matches as default behaviour
    return 1 if not @alleles;    
    foreach my $allele (@alleles){
        return 0 if not grep { /^$allele$/ } @snp_alleles;
    }
    return 1;
}        

=head1 NAME

filterOnEvsMaf.pl - filter NHLBI ESP variants from a VCF file 

=head1 SYNOPSIS

        filterOnEvsMaf.pl -i [vcf file] -d [directory containing ESP VCF files] [options]
        filterOnEvsMaf.pl -i [vcf file] -e [ESP VCF file(s)] [options]
        filterOnEvsMaf.pl -h (display help message)
        filterOnEvsMaf.pl -m (display manual page)

=cut

=head1 ARGUMENTS

=over 8

=item B<-i    --input>

Input VCF file 

=item B<-o    --output>

Output snp filtered file. Optional - default is STDOUT.

=item B<-d    --dir>

Directory containing ESP VCF files.

=item B<-e    --esp_file>

One or more ESP VCF files to use.

=item B<--samples>

One or more samples to check variants for.  Default is to check all variants specified by the ALT field.

=item B<-f    --freq>

Percent minor allele frequency to filter from. Only variants with minor alleles equal to or over this frequency will be removed. Default is to remove any matching variant.

=item B<-p    --progress>

Use this flag to print % progress to STDERR.

=item B<-q    --quiet>

Use this flag to supress warnings if any variants in the SNP reference file do not have IDs.

=item B<--strict>

Use this flag to complain and exit if any variants in the SNP reference file do not have IDs.

=item B<-h    --help>

Display help message.

=item B<-m    --manual>

Show manual page.

=back 

=cut


=head1 DESCRIPTION

This program will filter a VCF file to remove variants matching variants in VCF files from NHLBI ESP.  If --freq is set only variants in ESP with a MAF equal to or greater than --freq will be filtered.

=cut

=head1 AUTHOR

David A. Parry

University of Leeds

=head1 COPYRIGHT AND LICENSE

Copyright 2013  David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut

