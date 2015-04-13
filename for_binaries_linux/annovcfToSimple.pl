#!/usr/bin/perl
#David Parry August 2011
#
=head1 NAME

annovcfToSimple - takes a vcf file and outputs a simplified version in Excel's .xlsx format.

=head1 SYNOPSIS

        annovcfToSimple -i [ file] [options]
        
        annovcfToSimple -v -i [VEP annotated VCF file] [options]
        
        annovcfToSimple -g -v -i [VEP and ensemblGeneAnnotator.pl annotated VCF file] [options]
        
        annovcfToSimple -h (display help message)
        
        annovcfToSimple -m (display manual page)

=cut

=head1 ARGUMENTS

=over 8

=item B<-i    --input>

Input VCF file, optionally annotated with variant_effect_predictor.pl and ensemblGeneAnnotator.pl. 

=item B<-o    --output>

Output file name. Defaults to [input name].xlsx (or [input name].txt if --text_output option is used).

=item B<-u    --summarise>

Use this flag to summarise the number of alleles and genotypes found in samples rather than outputting genotype columns for each individual sample. If sample IDs are specified with --samples or --pedigree options only these samples will be counted in the summarised allele and genotype counts. This can be overriden by adding the word 'all' after this argument, in which case all samples in the VCF will be summarised and any samples specified by --samples or --pedigree options will still have their individual genotypes written to the output.

=item B<--contains_variant>

Use this flag alongside the -u/--summarise option to add a columns summarising which samples contain variant genotypes. 

=item B<-s    --samples>

Sample ids to include in output. By default, genotypes for all samples in VCF will be written to the output, but if one or more samples are specified here only genotypes for these samples will be included in the output. Note, that this is ignored if --do_not_simplify argument is used unless using the --summarise flag.

=item B<-p    --pedigree>

If one or more valid PED files are provided here, samples found in the PED files will be included in the output if they exist in the VCF file (as opposed to the default of including all samples in the VCF). Can be used instead of or in conjunction with --samples argument.

PED format - from the PLINK documentation:

       The PED file is a white-space (space or tab) delimited file: the first six columns are mandatory:

            Family ID
            Individual ID
            Paternal ID
            Maternal ID
            Sex (1=male; 2=female; other=unknown)
            Phenotype

=item B<-v    --vep>

Use this flag to output annotations from Ensembl's variant_effect_predictor.pl at the beginning of each line (assuming your input has already been annotated by the VEP).

=item B<-g    --gene_anno>

Use this flag if annotated with ensemblGeneAnnotator.pl to output gene annotation information. This script cannot distinguish between gene annotation information belonging to canonical transcripts or particular variant classes, these distinctions have to be made when running ensemblGeneAnnotator.pl prior to this script.

=item B<--canonical_only>

Use this flag to only print consequences from canonical transcripts when using --vep option. 

=item B<--functional>

Use this flag to only annotate standard 'functional' variant classes (transcript_ablation, splice_donor_variant, splice_acceptor_variant, splice_region_variant, stop_gained, frameshift_variant, stop_lost, initiator_codon_variant, inframe_insertion, inframe_deletion, missense_variant, transcript_amplification, TFBS_ablation, TFBS_amplification, regulatory_region_ablation, regulatory_region_amplification).

=item B<-c    --classes>

Use this to specify a set of VEP variant classes to print in output. Overrides --functional option.

=item B<-a    --additional_classes>

Use this to specify additional VEP variant classes to output as well as those used by --functional.

=item B<--fields>

Specify one or more VEP fields to output. Use of --vep without this flag outputs the following fields:

            symbol
            gene
            feature
            allele
            consequence
            cds_position
            protein_position
            amino_acids
            codons
            existing_variation
            exon
            intron
            splice_consensus
            sift
            polyphen
            condel
            gmaf
            aa_maf
            ea_maf
            afr_maf
            amr_maf
            asn_maf
            eur_maf

If you wish to use the above fields in addition to fields specified here add 'default' to your list of fields.
       
=item B<--all>

Use with --vep option to output all available VEP fields.

=item B<-n    --info_fields>

One or more INFO field IDs from to output as columns. These are case sensitive and must appear exactly as defined in the VCF header. By default, if --info_fields is not specified, CaddPhredScore INFO fields will be written to the output automatically if found.

=item B<-d    --do_not_simplify>

Use this flag to output all standard VCF fields as they appear in the original VCF, but when used in conjunction with --vep still provides information for VEP annotations in a user-friendly manner. Genotypes for all samples in the VCF will be printed when this option is used regardless of --samples or --pedigree settings. The only exception is if used with the --summarise flag, in which case summary allele/genotype counts will be output instead of indvidual genotypes, but the genotypes will use the numeric call codes as they appear in the VCF rather than the actual REF/ALT alleles.

=item B<-t    --text_output>

Use this flag to output tab-delimited text instead of Excel format. 

=item B<-h    --help>

Show help message.

=item B<-m    --manual>

Show manual page.


=back 

=cut

=head1 DESCRIPTION

Reads a VCF file and outputs a simplified version in Excel's .xlsx format. Useful when using a VCF annotated with Ensembl's variant_effect_predictor.pl to output a more easily readable version (use the -v option), and particularly useful for VCF files annotated with ensemblGeneAnnotator.pl (use the -g option).

=cut

=head1 AUTHOR

David A. Parry
University of Leeds


=head1 COPYRIGHT AND LICENSE

Copyright 2013,2014  David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut




use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Excel::Writer::XLSX;
use Excel::Writer::XLSX::Utility;
use Pod::Usage;
use File::Basename;
use FindBin;
use lib "$FindBin::Bin/lib";
use ParseVCF;
use ParsePedfile;
use TextToExcel;


############

my $vcf;
my @samples;
my @vep_classes;
my @vep_add;
my @peds;
my $output;
my $config = {};
my $summarise_counts;
my $summarise_samples_with_variants;
my $help;
my $man;
GetOptions($config,
    'c|classes=s{,}' => \@vep_classes,
    'a|additional_classes=s{,}' => \@vep_add,
    'functional',
    'gene_anno',
    'text_output',
    'do_not_simplify',
    'canonical_only',
    'vep',
    'fields=s{,}', => \@{$config->{fields}}, 
    'all',
    'n|info_fields=s{,}', => \@{$config->{info_fields}},
    'i|input=s' =>\$vcf,
    'output=s' => \$output,
    's|samples=s{,}' => \@samples,
    'u|summarise:s' => \$summarise_counts,
    'contains_variant' => \$summarise_samples_with_variants,
    'pedigree=s{,}' => \@peds,
    'manual' => \$man,
    'help' => \$help) or pod2usage(-exitval => 2, -message => "Syntax error.\n");
pod2usage( -verbose => 2 ) if $man;
pod2usage( -verbose => 1 ) if $help;
pod2usage( -exitval => 2, -message => "--input is required" ) if (not $vcf);

if (@vep_classes or @vep_add){
    $config->{functional} = 1;
}


my @vep_valid = qw (transcript_ablation
                splice_donor_variant
                splice_acceptor_variant
                stop_gained
                frameshift_variant
                stop_lost
                initiator_codon_variant
                inframe_insertion
                inframe_deletion
                missense_variant
                transcript_amplification
                splice_region_variant
                incomplete_terminal_codon_variant
                synonymous_variant
                stop_retained_variant
                coding_sequence_variant
                mature_miRNA_variant
                5_prime_UTR_variant
                3_prime_UTR_variant
                intron_variant
                NMD_transcript_variant
                non_coding_exon_variant
                nc_transcript_variant
                upstream_gene_variant
                downstream_gene_variant
                TFBS_ablation
                TFBS_amplification
                TF_binding_site_variant
                regulatory_region_variant
                regulatory_region_ablation
                regulatory_region_amplification
                feature_elongation
                feature_truncation
                intergenic_variant);

if (not @vep_classes){
        @vep_classes = qw (transcript_ablation
                splice_donor_variant
                splice_acceptor_variant
                splice_region_variant
                stop_gained
                frameshift_variant
                stop_lost
                initiator_codon_variant
                inframe_insertion
                inframe_deletion
                missense_variant
                transcript_amplification
                TFBS_ablation
                TFBS_amplification
                regulatory_region_ablation
                regulatory_region_amplification);
}
push (@vep_classes, @vep_add) if (@vep_add);
#check VEP classes
foreach my $class (@vep_classes){
        die "Error - variant class '$class' not recognised.\n" if not grep {/$class/i} @vep_valid;
}
#get all samples from any peds specified
my @ped_samples = ();
if (@peds){
    foreach my $pedigree (@peds){
        my $ped_obj  = ParsePedfile->new(file => $pedigree);
        push @ped_samples, $ped_obj->getAllSamples();
    }
    @ped_samples = sort @ped_samples;
    my %seen = ();
    @ped_samples = grep {! $seen{$_}++ } @ped_samples;
}


#excel related global variables
my $workbook ;
my $worksheet ;
my $header_formatting;    
my $std_formatting;     
my $url_format;


my @suffixes = (".vcf", ".txt");
my $out_ext = 'xlsx';
$out_ext = 'txt' if defined $config->{text_output};
if (not $output){
    my ($out, $dir, $extension) = fileparse($vcf, @suffixes);
    $output = "$dir/$out.$out_ext";
}else{
    $output .= ".$out_ext" if $output !~ /\.$out_ext$/;
}
my @fields = ();
my @info_fields = ();
my @gene_anno_fields = ();
my $vcf_obj = ParseVCF->new( file=> $vcf);
my %info_fields = $vcf_obj->getInfoFields();
if (defined $config->{gene_anno}){
    if (not exists $info_fields{GeneAnno}){
        die "GeneAnno INFO field not found in VCF header - please annotate ".
            "your VCF with ensemblGeneAnnotator.pl or rerun this program without ".
            "the --gene_anno option.\n";
    }
    if ($info_fields{GeneAnno}->{Description} =~ /Format: ([\w+\|]+)\"/){
        @gene_anno_fields = split(/\|/, $1);
    }else{
        die "Could not parse GeneAnno description field - please annotate ".
            "your VCF with ensemblGeneAnnotator.pl or report this error if it recurs.\n";
    }
}
if (defined $config->{vep}){
    my $vep_header = $vcf_obj->readVepHeader();
    if (not @{$config->{fields}} and not $config->{all}){
        unshift @{$config->{fields}}, getDefaultVepFields($vep_header);
    }elsif(grep { /default/i } @{$config->{fields}}){
        @{$config->{fields}} = grep {! /^default$/i } @{$config->{fields}};
        unshift @{$config->{fields}}, getDefaultVepFields($vep_header);
    }
    @fields = map { lc($_) } @{$config->{fields}} if defined $config->{fields};
    my %seen = ();
    @fields = grep {! $seen{$_}++ } @fields;
    if (@fields){
        if(defined $config->{canonical_only}){
            push @fields, 'canonical' if not grep{/canonical/i} @fields;
        }
        foreach my $csq (@fields){
            if (not exists $vep_header->{($csq)}){
                die "Couldn't find '$csq' VEP field in header - please ensure your VCF is annotated with " .
                "Ensembl's variant effect precictor specifying the appropriate annotations.\n";
            }
        }
    }else{
        if(defined $config->{canonical_only}){
            if (not exists $vep_header->{canonical}){
                print STDERR "Canonical field not found in VEP header - all transcripts will be assessed.\n";
                $config->{canonical_only} = undef;
            }
        }
        foreach my $csq (sort {$vep_header->{$a} <=> $vep_header->{$b}} keys %{$vep_header}){
            push @fields, $csq;
        }
    }
}
if (@{$config->{info_fields}}){
    foreach my $f (@{$config->{info_fields}}){
        if (exists $info_fields{$f}){
            push @info_fields, $f;
        }else{
            print STDERR "INFO field $f not found in header and will not be included in output.\n";
        }
    }
}elsif(exists $info_fields{CaddPhredScore}){
#include CaddPhredScore annotation by default if it exists
    push @info_fields, "CaddPhredScore";
}


if (@ped_samples){
    my @sample_found = ();
    foreach my $s (@ped_samples){
        if ($vcf_obj->checkSampleInVcf($s)){
            push @sample_found, $s;
        }
    }
    if (not @sample_found){
        die "No samples from ped file(s) found in VCF.\n";
    }
    print STDERR "Found " . scalar(@sample_found) . " samples of " . scalar(@ped_samples) . 
        " samples in ped files in VCF.\n";
    push @samples, @sample_found;
}
#check samples or get them if none specified
if (not @samples){
    eval{
    @samples = $vcf_obj->getSampleNames();
    };
    warn "$@" if $@;
}else{
    my $header_string = $vcf_obj->getHeader(1);
    foreach my $sample (@samples){
        die "Can't find sample $sample in header line.\nHeader:\n$header_string\n" 
            if not $vcf_obj->checkSampleInVcf($sample);
    }
}

my $OUT;
my $xl_obj;
#old way
#my ($header, $gene_anno_cols) = get_header();
#new way
my @header = get_header();
if (defined $config->{text_output}){
    if ($output){
        open ($OUT, ">$output") or die "Can't open $output: $!\n";
    }else{
        $OUT = \*STDOUT;
    }
    foreach my $h (@header){
        $h =~ s/ /_/g;
    }
    print $OUT join("\t", @header) ."\n";
}else{
    $xl_obj = TextToExcel->new( file=> $output);
    $header_formatting = $xl_obj->createFormat(bold => 1);
    $xl_obj->writeLine(line => \@header, format => $header_formatting);
}

LINE: while (my $vcf_line = $vcf_obj->readLine){
    my @preceding = ();
    my @line = ();
    my @following = ();
    #Get preceding VEP columns (there may be multiple rows per line)
    if (defined $config->{vep}){
        my ($prec, $canonical_found, $functional_found, $sample_found) = get_vep_fields() ; 
        if (@$prec > 0){
            @preceding = @$prec;
        }else{
            my $e_string;
            my @temp_prec;
            if (not $canonical_found){
                $e_string = "No canonical transcript variant";
            }elsif (not $functional_found){
                $e_string = "No valid functional variant";
            }elsif (not $sample_found){
                $e_string = "No valid variant in samples";
            }
            push @temp_prec, $e_string;
            for (my $i = 1; $i < @fields; $i++){
                push @temp_prec, "-";
            }
            push @preceding, \@temp_prec;
        }

    }
    if (defined $config->{gene_anno}){
        my $g_anno = $vcf_obj->getVariantInfoField('GeneAnno');
        my @g_annots = split(",", $g_anno);
        foreach my $g (@g_annots){
            my @temp_follow = ();
            #Parse each field from @gene_anno_fields
            #convert back to normal text without ensemblGeneAnnotator replacements
            $g =~ tr/^`_/;, /;
            my @annots = split(/\|/, $g);
            foreach my $ann (@annots){
                $ann =~ s/::/\|/g;#use | as multivalue separator (ensemblGeneAnnotator uses ::)
                push @temp_follow, $ann;
            }
            push @following, \@temp_follow;
        }
    }
    #Get values for line
    if (defined $config->{do_not_simplify} and not defined $summarise_counts){
        @line = get_vcf_fields();    
    }else{
        @line = get_simplified_fields();    
    }
    #old way
#    if (ref $gene_anno_cols eq 'ARRAY'){
#        foreach my $field (@$gene_anno_cols){
#            push @line, $vcf_obj->getCustomField($field);
#        }
#    }
    #now write to text or excel as appropriate
    if (defined $config->{text_output}){
        my @fol_text = ();
        if (@following){
        #for text format we'll keep the gene annotations the same
        #for each line, enclosing values for different genes with []
            for (my $i = 0; $i < @{$following[0]}; $i++){
                my @temp_fol = ($following[0]->[$i]);
                for (my $j = 1; $j < @following; $j++){
                    push @temp_fol, $following[$j]->[$i];
                }
                if (@following > 1){
                    foreach my $t (@temp_fol){
                        $t = "[$t]";
                    }
                }
                push @fol_text, join("", @temp_fol);
            }
        }
        #@preceding will often span multiple rows
        if (@preceding){
            foreach my $p (@preceding){
                print $OUT join("\t", @$p) ."\t";
                print $OUT join("\t", @line);
                if (@fol_text){
                    print $OUT "\t" . join("\t", @fol_text);
                }
                print $OUT "\n";
            }
        }else{
            print $OUT join("\t", @line);
            if (@fol_text){
                print $OUT "\t" . join("\t", @fol_text);
            }
            print $OUT "\n";
        }
    }else{
        $xl_obj->writeLine(
                line => \@line, 
                preceding => \@preceding, 
                succeeding => \@following, 
                );
    }
}
if ($xl_obj){
    $xl_obj->DESTROY();
}
if ($OUT){
    close $OUT;
}
$vcf_obj->DESTROY();

####################################################
sub getDefaultVepFields{
    my ($vep_header) = @_;
    my @def = ();
    push @def, 
        qw ( 
            symbol
            gene
            feature
            allele
            consequence
            cds_position
            protein_position
            amino_acids
            codons
            existing_variation
            exon
            intron
            gmaf 
            aa_maf
            ea_maf
            afr_maf
            amr_maf
            asn_maf
            eur_maf
            sift
            polyphen
        );
    push @def, "condel" if exists $vep_header->{condel};
    push @def, "splice_consensus" if exists $vep_header->{splice_consensus};
    push @def, "hgvsc" if exists $vep_header->{hgvsc};
    push @def, "hgvsp" if exists $vep_header->{hgvsp};
    return @def;
}
####################################################
sub get_vcf_fields{
    my @vcf_values = ();
    foreach my $field (split("\t", $vcf_obj->get_currentLine)){
        push @vcf_values, $field;
    }
    #Add user-specified INFO fields after normal VCF fields
    foreach my $f (@info_fields){
        my $value = $vcf_obj->getVariantInfoField($f);
        if (defined $value){
            push @vcf_values, $value;
        }else{
            push @vcf_values, '.';
        }
    }
    return @vcf_values;
}

####################################################
sub get_simplified_fields{
    my @vcf_values = ();
    #Get preceding INFO fields (these will be one per line)
    foreach my $f (qw (CHROM POS ID QUAL FILTER REF ALT)){ 
        push @vcf_values, $vcf_obj->getVariantField($f);
    }
    my @all_alleles = $vcf_obj->readAlleles();
    if (defined $summarise_counts){
        my %allele_counts = ();
        my %genotype_counts = ();
        if ($summarise_counts eq 'all'){
            %allele_counts = $vcf_obj->countAlleles();
            %genotype_counts = $vcf_obj->countGenotypes();
        }else{
            %allele_counts = $vcf_obj->countAlleles(samples => \@samples);
            %genotype_counts = $vcf_obj->countGenotypes(samples => \@samples);
        }
            
        my @al_count_string = ();
        my @gt_count_string = ();
        my @sample_var_string = ();
        foreach my $allele (sort keys %allele_counts){
            if (defined $config->{do_not_simplify}){
                push @al_count_string, "$allele=$allele_counts{$allele}";
            }else{
                push @al_count_string, "$all_alleles[$allele]=$allele_counts{$allele}";
            }
        }
        foreach my $genotype (sort keys %genotype_counts){
            my @gt = split(/[\/\|]/, $genotype);
            my @actual_gt = ();
            if (defined $config->{do_not_simplify}){
                push @gt_count_string, "$genotype=$genotype_counts{$genotype}";
            }else{
                foreach my $g (@gt){
                    if ($g =~ /^\d+$/){
                        push @actual_gt, $all_alleles[$g];
                    }else{
                        push @actual_gt, $g; 
                    }
                }
                my $actual_geno_string = join("/", @actual_gt);
                push @gt_count_string, "$actual_geno_string=$genotype_counts{$genotype}";
            }
        }
        if (defined $summarise_samples_with_variants){
            my %calls = ();
            my %var_to_samples = ();
            if ($summarise_counts eq 'all'){
                %calls = $vcf_obj->getSampleCall(all => 1);
            }else{
                %calls = $vcf_obj->getSampleCall(multiple => \@samples);
            }
            foreach my $sample (sort keys %calls){
                next if ($calls{$sample} =~/^0[\|\/]0$/);
                next if ($calls{$sample} =~/^\.[\|\/]\.$/);
                push @{$var_to_samples{$calls{$sample}}}, $sample;
            }
            foreach my $genotype (sort keys %var_to_samples){
                my $sum_string = join(",", @{$var_to_samples{$genotype}});
                if (! defined $config->{do_not_simplify}){
                    my @gt = split(/[\/\|]/, $genotype);
                    my @actual_gt = ();
                    foreach my $g (@gt){
                        if ($g =~ /^\d+$/){
                            push @actual_gt, $all_alleles[$g];
                        }else{
                            push @actual_gt, $g; 
                        }
                    }
                    $genotype = join("/", @actual_gt);
                }
                push @sample_var_string, "$genotype=$sum_string";
                
            }
        }
        push @vcf_values, join(";", @al_count_string);
        push @vcf_values, join(";", @gt_count_string);
        push @vcf_values, join(";", @sample_var_string) if defined $summarise_samples_with_variants;
    }
    if (not defined $summarise_counts or (defined $summarise_counts and $summarise_counts eq 'all')){
        my @sample_calls = ();
        my @sample_allele_depths = ();
        my @sample_genotype_quality = ();
        foreach my $sample (@samples){
            my $genotype = $vcf_obj->getSampleActualGenotypes(sample=>$sample);
            if (defined $genotype){ 
                push (@sample_calls, $genotype);
            }else{
                push @sample_calls, "-";
            }
            my $ad = $vcf_obj->getSampleGenotypeField(sample => $sample, field => "AD"); 
            my $allele_depths = '-';
            if (defined $ad && $ad =~ /(\d+,*)+/){
                my @ads = split(",", $ad);
                my @allele_depth = ();
                for (my $n = 0; $n < @ads and $n < @all_alleles; $n++){
                    push (@allele_depth, "$all_alleles[$n]=$ads[$n]");
                }
                $allele_depths = join("/", @allele_depth);
            }
            push (@sample_allele_depths, $allele_depths);
            my $genotype_quality = "-";
            my $gq = $vcf_obj->getSampleGenotypeField(sample => $sample, field => "GQ");
            if (defined $gq && $gq =~ /(\d+\.*\d*)+/){
        #     $genotype_quality = 100 - (100 * 10**(-$gq/10));
                $genotype_quality = $gq;
            }
            push (@sample_genotype_quality, $genotype_quality);
        }
        foreach my $call (@sample_calls){
            push @vcf_values, $call;
        }
        foreach my $sample_depth (@sample_allele_depths){
            push @vcf_values, $sample_depth;
        }
        foreach my $sample_gq (@sample_genotype_quality){
            push @vcf_values, $sample_gq;
        }
    }
    #Add user-specified INFO fields after normal VCF fields
    foreach my $f (@info_fields){
        my $value = $vcf_obj->getVariantInfoField($f);
        if (defined $value){
            push @vcf_values, $value;
        }else{
            push @vcf_values, '.';
        }
    }
    return @vcf_values;
}

####################################################
sub get_vep_fields{
    #returns a 2D array of VEP values
    my ($row, $col) = @_;
    my @csq = $vcf_obj->getVepFields(\@fields);
    my ($canonical_found, $functional_found, $sample_found) = (0, 0, 0);
    my %vep_allele = ();
    my @vep_values = ();
    if ($config->{samples}){
        my @sample_alleles = $vcf_obj->getSampleActualGenotypes(multiple => $config->{samples}, return_alleles_only => 1);
        foreach my $allele (@sample_alleles){
            my $v_allele = $vcf_obj->altsToVepAllele(alt => $allele);
            $vep_allele{$v_allele} = $allele;
        }
     }
ANNOT:  foreach my $annot (@csq){
        my @csq_values = ();
        if ($config->{canonical_only}){
              next if (not $annot->{'canonical'} );
        }
        $canonical_found++;
        if ($config->{functional}){
            my $match = 0;
CLASS:      foreach my $class (@vep_classes){
                my @anno_csq = split(/\&/, $annot->{consequence});
                if (grep {/NMD_transcript_variant/i} @anno_csq){
                    next ANNOT;
                }else{
                    foreach my $ac (@anno_csq){
                        if (lc$ac eq lc$class){
                            $match++;
                            last CLASS;
                        }
                    }
                }
            }
            next ANNOT if not $match;
        }
        $functional_found++;
        if ($config->{samples}){
            next if not exists $vep_allele{$annot->{allele}};
        }
        $sample_found++;
        foreach my $vep  (@fields){
            my $value = $annot->{$vep};
            if (defined $value){
                push @csq_values, $value;
            }else{
                push @csq_values, '-';
            }
        }
        push @vep_values, \@csq_values;
    }   
    return (\@vep_values, $canonical_found, $functional_found, $sample_found);
}
####################################################################
sub get_header{
    my @head = ();
    #first deal with VEP fields
    my $vep_header;
    if (defined $config->{vep}){
        foreach my $csq (@fields){
            push @head, $csq;
        }
    }
    #GET HEADER INFO
    my $header_string = $vcf_obj->getHeader(1);
    my %header_columns = $vcf_obj->getHeaderColumns();
    my @head_split= split("\t", $header_string);
    my @anno_columns = ();
    if (defined $config->{do_not_simplify}){
        foreach my $h (@head_split){
            push @head, $h;
        }
        #Put user-specified INFO fields after normal vcf fields
        if (@info_fields){
            foreach my $f (@info_fields){
                push @head, $f;
            }
        }
    }else{
        foreach my $h ("Chrom", "Pos", "SNP ID", "Variant Quality", "Filters", "Genomic Ref", "Alt Alleles" ){
            push @head, $h;
        }
        if (defined $summarise_counts){
            push @head, "Allele Counts";
            push @head, "Genotype Counts";
            if (defined $summarise_samples_with_variants){
                push @head, "Samples with variant";
            }
        }
        if (not defined $summarise_counts or (defined $summarise_counts and $summarise_counts eq 'all')){
            foreach my $sample (@samples){
                push @head, $sample;
            }

            foreach my $sample (@samples){
                push @head, "$sample Allele Depth";
            }
            foreach my $sample (@samples){
                push @head, "$sample Phred Genotype Confidence";
            }
        }
        #Put user-specified INFO fields after normal vcf fields
        if (@info_fields){
            foreach my $f (@info_fields){
                push @head, $f;
            }
        }

        if (defined $config->{gene_anno}){
        #old way
#            my $i = 0; 
#            $i++ until $head_split[$i] =~ /CHROM/;
#            @anno_columns = @head_split[0..$i-1];
#            foreach my $anno_col (@anno_columns){
#                push @head, $anno_col;
#            }
        #new way
            foreach my $g (@gene_anno_fields){
                push @head, $g;
            }
        }

    }
#    return \@head, \@anno_columns;
    return @head;
}



