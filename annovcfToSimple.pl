#!/usr/bin/perl
#David Parry August 2011
#
=head1 NAME

annovcfToSimple.pl - takes a vcf file and outputs a simplified version in Excel's .xlsx format.

=head1 SYNOPSIS

        annovcfToSimple.pl -i [ file] [options]
        
        annovcfToSimple.pl -v -i [VEP annotated VCF file] [options]
        
        annovcfToSimple.pl -g -v -i [VEP and ensemblGeneAnnotator.pl annotated VCF file] [options]
        
        annovcfToSimple.pl -h (display help message)
        
        annovcfToSimple.pl -m (display manual page)

=cut

=head1 ARGUMENTS

=over 8

=item B<-i    --input>

Input VCF file, optionally annotated with variant_effect_predictor.pl and ensemblGeneAnnotator.pl. 

=item B<-o    --output>

Output file name. Defaults to [input name].xlsx (or [input name].txt if --text_output option is used).

=item B<-s    --samples>

Sample ids to include in output. By default, genotypes for all samples in VCF will be written to the output, but if one or more samples are specified here only genotypes for these samples will be included in the output. Note, that this is ignored if --do_not_simplify argument is used.

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

=item B<--classes>

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
       
=item B<--all>

Use with --vep option to output all available VEP fields.

=item B<-d    --do_not_simplify>

Use this flag to output all standard VCF fields to Excel as they appear in the original VCF, but when used in conjunction with --vep still provides information for VEP annotations in a user-friendly manner. Genotypes for all samples in the VCF will be printed when this option is used regardless of --samples or --pedigree settings.

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
use lib "$FindBin::Bin";
use ParseVCF;
use ParsePedfile;


############

my $vcf;
my @samples;
my @classes;
my @add;
my @peds;
my $output;
my $config = {};
my $help;
my $man;
GetOptions($config,
    'classes=s{,}' => \@classes,
    'additional_classes=s{,}' => \@add,
    'functional',
    'gene_anno',
    'text_output',
    'do_not_simplify',
    'canonical_only',
    'vep',
    'fields=s{,}', => @{$config->{fields}}, 
    'all',
    'input=s' =>\$vcf,
    'output=s' => \$output,
    'samples=s{,}' => \@samples,
    'pedigree=s{,}' => \@peds,
    'manual' => \$man,
    'help' => \$help) or die "Syntax error.\n";
pod2usage( -verbose => 2 ) if $man;
pod2usage( -verbose => 1 ) if $help;
pod2usage( -exitval => 2, -message => "--input is required" ) if (not $vcf);

if (@classes or @add){
    $config->{functional} = 1;
}

if (not defined @{$config->{fields}} and not $config->{all}){
    push @{$config->{fields}}, 
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
            );
}



my @valid = qw (transcript_ablation
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

if (not @classes){
        @classes = qw (transcript_ablation
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
push (@classes, @add) if (@add);
#check VEP classes
foreach my $class (@classes){
        die "Error - variant class '$class' not recognised.\n" if not grep {/$class/i} @valid;
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
my @fields;
my $vcf_obj = ParseVCF->new( file=> $vcf);
if (defined $config->{vep}){
    my $vep_header = $vcf_obj->readVepHeader();
    @fields = @{$config->{fields}} if defined $config->{fields};
    if (@fields){
        if(defined $config->{canonical_only}){
            push @fields, 'canonical' if not grep{/canonical/i} @fields;
        }
        foreach my $csq (@fields){
            if (not exists $vep_header->{$csq}){
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
    print "Found " . scalar(@sample_found) . " samples of " . scalar(@ped_samples) . 
        " samples in ped files in VCF.\n";
    push @samples, @sample_found;
}
#check samples or get them if none specified
if (not @samples){
    @samples = $vcf_obj->getSampleNames();
}else{
    my $header_string = $vcf_obj->getHeader(1);
    foreach my $sample (@samples){
        die "Can't find sample $sample in header line.\nHeader:\n$header_string\n" 
            if not $vcf_obj->checkSampleInVcf($sample);
    }
}


my $OUT;
if (defined $config->{text_output}){
    open ($OUT, ">$output") or die "Can't open $output: $!\n";
    process_as_text();
}else{
    process_as_xlsx();
}




###########################
sub process_as_text{
    open (my $OUT, ">$output") or die "Can't open $output: $!\n";
    if (defined $config->{vep}){
        print $OUT "#";
        print $OUT join("\t", @fields) ."\t";
    }
    #GET HEADER INFO
    my $header_string = $vcf_obj->getHeader(1);
    my %header_columns = $vcf_obj->getHeaderColumns();
    my @head_split= split("\t", $header_string);
    my @anno_columns  = ();
    if (defined $config->{do_not_simplify}){
        print $OUT join("\t", @head_split);
        
    }else{
        #print $OUT join("\t", ("Chrom", "Pos", "SNP_ID", "Variant Percent Confidence", "Filters", "Genomic Ref", ));
        print $OUT join("\t", ("Chrom", "Pos", "SNP_ID", "Variant Phred Quality", "Filters", "Genomic Ref", ));
        print $OUT join("\t", @samples) ."\t" ;
        print $OUT join(" Allele Depth\t", @samples) ." Allele Depth\t" ;
        print $OUT join(" Phred Genotype Confidence\t", @samples) ." Phred Genotype Confidence\t" ;
 #       print $OUT join(" Scaled Percent Probably Liklihoods\t", @samples) ." Scaled Percent Probably Liklihoods" ;
        if (defined $config->{gene_anno}){
            my $i = 0; 
            $i++ until $head_split[$i] =~ /CHROM/;
            @anno_columns = @head_split[0..$i-1];
            print $OUT join("\t", @anno_columns);
        }
        print $OUT "\n";
    }


    LINE: while (my $vcf_line = $vcf_obj->readLine){
        if (defined $config->{do_not_simplify}){
            write_line_to_text()
        }else{
            write_simplified_vcf_to_text(\@anno_columns);
        }
    }

}    

####################################################
sub write_line_to_text{
    my @output = split("\t", $vcf_obj->get_currentLine);
    if (defined $config->{vep}){
        my @csq = ();
        @csq = $vcf_obj->getVepFields(\@fields);
        my %vep_allele ; 
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
                  next if (not $annot->{canonical} );
              }
            if ($config->{functional}){
                my $match = 0;
CLASS:          foreach my $class (@classes){
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
            if ($config->{samples}){
                next if not exists $vep_allele{$annot->{allele}};
            }
            foreach my $vep  (@fields){
                my $value = $annot->{$vep};
                if (defined $value && $value ne ''){
                    push @csq_values, $value;
                }else{
                    push @csq_values, ".";
                }
            }
            print $OUT join("\t", @csq_values) ."\t";
            print $OUT join("\t", @output) ."\n";
        }
    }else{
        print $OUT join("\t", @output) ."\n";
    }
}

####################################################
sub write_simplified_vcf_to_text{
    my ($gene_anno_cols) = @_;
    my @output = ();
        
    foreach my $varfield ($vcf_obj->getVariantField("CHROM"), $vcf_obj->getVariantField("POS"), $vcf_obj->getVariantField("ID"),){ 
        push @output, $varfield;
    }
    #push @output, 100 - (100 * (10**(-($vcf_obj->getVariantField("QUAL"))/10)));
    push @output, $vcf_obj->getVariantField("QUAL");
    push @output, $vcf_obj->getVariantField("FILTER");
    push @output, $vcf_obj->getVariantField("REF");
    my @form = split(":", $vcf_obj->getVariantField("FORMAT"));
    my @all_alleles = ($vcf_obj->getVariantField("REF"), split(",", $vcf_obj->getVariantField("ALT")));
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
            for (my $n = 0; $n < @ads; $n++){
                push (@allele_depth, "$all_alleles[$n]=$ads[$n]");
            }
            $allele_depths = join("/", @allele_depth);
        }
        push (@sample_allele_depths, $allele_depths);
        my $genotype_quality = "-";
        my $gq = $vcf_obj->getSampleGenotypeField(sample => $sample, field => "GQ");
        if (defined $gq && $gq =~ /(\d+\.*\d*)+/){
            #$genotype_quality = 100 - (100 * 10**(-$gq/10));
            $genotype_quality = $gq;
        }
        push (@sample_genotype_quality, $genotype_quality);
    }
    foreach my $call (@sample_calls){
        push @output, $call;
    }
    foreach my $sample_depth (@sample_allele_depths){
        push @output, $sample_depth;
    }
    foreach my $sample_gq (@sample_genotype_quality){
        push @output, $sample_gq;
    }
    if ($config->{vep}){
        my @csq = $vcf_obj->getVepFields(\@fields);
        my %vep_allele ; 
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
                next if (not $annot->{canonical} );
            }
            if ($config->{functional}){
                my $match = 0;
CLASS:          foreach my $class (@classes){
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

            if ($config->{samples}){
                next if not exists $vep_allele{$annot->{allele}};
            }
            foreach my $vep  (@fields){
                my $value = $annot->{$vep};
                if (defined $value && $value ne ''){
                    push @csq_values, $value;
                }else{
                    push @csq_values, ".";
                }
                }
            print $OUT join("\t", @csq_values) ."\t";
            print $OUT join("\t", @output);
        }
    }else{
        print $OUT join("\t", @output);
    }
    if (ref $gene_anno_cols eq 'ARRAY'){
        my @geneanno;
        foreach my $field (@$gene_anno_cols){
            push @geneanno, $vcf_obj->getCustomField($field);
        }
        print $OUT join("\t", @geneanno);
    }
    print $OUT "\n";
}

    
###########################

sub process_as_xlsx{
    $workbook  = Excel::Writer::XLSX->new($output);
    $worksheet = $workbook->add_worksheet();
    #my ($ent, $go_acc, $go_desc, $rif, $sum) = (0,0,0,0,0);
    my $smpl;
    my $smend;
    my $col = 0;
    my $row = 0;
    $header_formatting = $workbook->add_format(bold => 1);
    $std_formatting = $workbook->add_format();
    $url_format = $workbook->add_format(color => 'blue', underline => 1, );
#    $worksheet->set_column( 0, 0, 10);
#    $worksheet->set_column( 1, 1, 12);
#    $worksheet->set_column( 2, 2, 25);
#    $worksheet->set_column( 4, 4, 9);
#    $worksheet->set_column( 10, 10, 11);
#    $worksheet->set_column( 11, 11, 13);
#    $worksheet->set_column( 12, 15, 20);
    my $vep_header;
    if (defined $config->{vep}){
        foreach my $csq (@fields){
            $worksheet->write($row, $col++, $csq, $header_formatting);
        }
    }
    #GET HEADER INFO
    my $header_string = $vcf_obj->getHeader(1);
    my %header_columns = $vcf_obj->getHeaderColumns();
    my @head_split= split("\t", $header_string);
    my @anno_columns = ();
    if (defined $config->{do_not_simplify}){
        foreach (@head_split){
            $worksheet->write($row, $col++, $_, $header_formatting);
        }
    }else{
        foreach ("Chrom", "Pos", "SNP_ID", "Variant Percent Confidence", "Filters", "Genomic Ref", ){
            $worksheet->write($row, $col++, $_, $header_formatting);
        }

        foreach my $sample (@samples){
            $worksheet->write($row, $col++, $sample, $header_formatting);
        }

        foreach my $sample (@samples){
            $worksheet->write($row, $col++, "$sample Allele Depth", $header_formatting);
        }
        foreach my $sample (@samples){
            $worksheet->write($row, $col++, "$sample Phred Genotype Confidence", $header_formatting);
        }
        if (defined $config->{gene_anno}){
            my $i = 0; 
            $i++ until $head_split[$i] =~ /CHROM/;
            @anno_columns = @head_split[0..$i-1];
            foreach my $anno_col (@anno_columns){
                $worksheet->write($row, $col++, $anno_col, $header_formatting);
            }
        }
    }


    LINE: while (my $vcf_line = $vcf_obj->readLine){
        $row++;
        $col = 0;
        if (defined $config->{do_not_simplify}){
            $row += write_line_to_excel($row, $col, );
        }else{
            $row += write_simplified_vcf_to_excel($row, $col, \@anno_columns);
        }
    }
    $workbook -> close();

}    


####################################################
sub write_vep_fields{
    my ($row, $col) = @_;
    my @csq = $vcf_obj->getVepFields(\@fields);
    my $lines = 0;
    my ($canonical_found, $functional_found, $sample_found) = (0, 0, 0);
    my $temp_col;
    my %vep_allele = ();
    if ($config->{samples}){
        my @sample_alleles = $vcf_obj->getSampleActualGenotypes(multiple => $config->{samples}, return_alleles_only => 1);
        foreach my $allele (@sample_alleles){
            my $v_allele = $vcf_obj->altsToVepAllele(alt => $allele);
            $vep_allele{$v_allele} = $allele;
        }
     }
ANNOT:  foreach my $annot (@csq){
        if ($config->{canonical_only}){
              next if (not $annot->{'canonical'} );
        }
        $canonical_found++;
        if ($config->{functional}){
            my $match = 0;
CLASS:      foreach my $class (@classes){
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
        $temp_col = $col;
        foreach my $vep  (@fields){
            my $value = $annot->{$vep};
            $worksheet->write($row, $temp_col++, $value);
        }
        $lines++;
        $row++;
    }   
    return ($lines, $temp_col, $canonical_found, $functional_found, $sample_found);
}
####################################################
sub write_line_to_excel{
#return no. of merged rows in order to increment $row
    my ($row, $col, ) = @_;
    my $lines = 1;
    my ($canonical_found, $functional_found, $sample_found) = (0, 0, 0);
    if (defined $config->{vep}){
        ($lines, $col, $canonical_found, $functional_found, $sample_found) = write_vep_fields($row, $col);
    }
    if ($lines < 1){
    #this happens if --canonical_only argument is in effect and we haven't got a canonical transcript for this variant
    #or if we haven't found a functional variant 
    #or if we haven't got a variant in our sample
    #we won't have printed anything for the consequence fields
        my $e_string ;
        if (not $canonical_found){
            $e_string = "No canonical variant";
        }elsif (not $functional_found){
           $e_string = "No valid functional variant";
        }elsif (not $sample_found){
            $e_string = "No valid variant in samples";
        }
        $worksheet->write($row, $col++, $e_string);
        for (my $i = 1; $i < @fields; $i++){
            $worksheet->write($row, $col++, "-");
        }
        $lines = 1;
    }
    foreach my $field (split("\t", $vcf_obj->get_currentLine)){
        write_worksheet($lines, $field, $worksheet, $row, $col++, $std_formatting);
    }
    #for (my $i = 0; $i < $lines; $i++){
    #    my $temp_col = $col;
    #    foreach my $field (split("\t", $vcf_obj->get_currentLine)){
    #        $worksheet->write($row, $temp_col++, $field);
    #    }
    #    $row++;
    #}
    return $lines - 1 ;

}
####################################################
sub write_simplified_vcf_to_excel{
#return no. of merged rows in order to increment $row
    my ($row, $col, $gene_anno_cols) = @_;
    my $lines = 1;
    my ($canonical_found, $functional_found, $sample_found) = (0, 0, 0);
    if (defined $config->{vep}){
        ($lines, $col, $canonical_found, $functional_found, $sample_found) = write_vep_fields($row, $col);
    }
    if ($lines < 1){
    #this happens if --canonical_only argument is in effect and we haven't got a canonical transcript for this variant
    #or if we haven't found a functional variant 
    #or if we haven't got a variant in our sample
    #we won't have printed anything for the consequence fields
        my $e_string ;
        if (not $canonical_found){
            $e_string = "No canonical variant";
        }elsif (not $functional_found){
           $e_string = "No valid functional variant";
        }elsif (not $sample_found){
            $e_string = "No valid variant in samples";
        }
        $worksheet->write($row, $col++, $e_string);
        for (my $i = 1; $i < @fields; $i++){
            $worksheet->write($row, $col++, "-");
        }
        $lines = 1;
    }
    foreach my $varfield ($vcf_obj->getVariantField("CHROM"), $vcf_obj->getVariantField("POS"), $vcf_obj->getVariantField("ID"),){ 
        write_worksheet($lines, $varfield, $worksheet, $row, $col++, $std_formatting);
    }
    #write_worksheet($lines, 100 - (100 * (10**(-($vcf_obj->getVariantField("QUAL"))/10))), $worksheet, $row, $col++, $std_formatting, );
    write_worksheet($lines, $vcf_obj->getVariantField("QUAL"), $worksheet, $row, $col++, $std_formatting, );
    write_worksheet($lines, $vcf_obj->getVariantField("FILTER"), $worksheet, $row, $col++, $std_formatting);
    write_worksheet($lines, $vcf_obj->getVariantField("REF"), $worksheet, $row, $col++, $std_formatting);
    my @all_alleles = $vcf_obj->readAlleles();
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
            for (my $n = 0; $n < @ads; $n++){
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
        write_worksheet($lines, $call, $worksheet, $row, $col++, $std_formatting);
    }
    foreach my $sample_depth (@sample_allele_depths){
        write_worksheet($lines, $sample_depth, $worksheet, $row, $col++, $std_formatting);
    }
    foreach my $sample_gq (@sample_genotype_quality){
        write_worksheet($lines, $sample_gq, $worksheet, $row, $col++, $std_formatting, );
    }
    if (ref $gene_anno_cols eq 'ARRAY'){
        foreach my $field (@$gene_anno_cols){
            write_worksheet($lines, $vcf_obj->getCustomField($field), $worksheet, $row, $col++, $std_formatting);
        }
    }
    return $lines - 1 ;#return additional rows added due to multiple annotations
}



####################################################
sub get_possible_genotypes{
    my ($allele_ref) = @_;
    my @combinations = ();
    for (my $n = 0; $n < @$allele_ref; $n++){
        for (my $m = 0; $m <= $n; $m++){
            push (@combinations, "$$allele_ref[$m]/$$allele_ref[$n]");
        }
    }
    return @combinations;
}

####################################################
sub write_worksheet{
#if we have multiple geneIDs then we need to write as a merged cell
    my ($gene_ids, $string, $worksheet, $row, $col, $formatting, $type) = @_;
    $type = 'string' if not $type;
    if ($string =~ /^\d+(\.\d+)*]$/){
        $type = 'number';
    }
    if ($gene_ids > 1){
        my $top_cell = xl_rowcol_to_cell($row, $col);
            my $bottom_cell = xl_rowcol_to_cell($row + $gene_ids-1, $col);
        $worksheet->merge_range_type($type, "$top_cell:$bottom_cell", $string, $formatting);
    }else{
        $worksheet->write($row, $col, $string);
    }
}
