#!/usr/bin/perl
#David Parry August 2011
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Excel::Writer::XLSX;
use Excel::Writer::XLSX::Utility;
use File::Basename;
use FindBin;
use lib "$FindBin::Bin";
use ParseVCF;

############
sub usage{
    my ($help) = @_;
    if ($help){
        print "\nTakes an annotated vcf file optionally annotated with Ensembl's variant_effect_predictor.pl and outputs a simplified version in Excel's .xlsx format.\n\n" .
        "Usage:  $0 -i [annotated vcf file] [options]\n\n" .
        "Arguments:\n\n" .
        "--input [annotated vcf input]\n". 
        "--output [output file]\n". 
        "--samples [sample ids to output, optional. Does not currently work with --do_not_simplify argument.]\n". 
        "--vep [use this flag to output annotations from Ensembl's variant_effect_predictor.pl at the beginning of each line]\n". 
        "--functional [use this flag to only output standard functional variant classes]\n".
        "--classes [use this to specify VEP variant classes to output]\n".
        "--additional_classes [use this to specify additional VEP variant classes to output as well as those used by --functional]\n".
        "--fields [specify one or more VEP fields to output. Use of --vep without this flag outputs all available fields]\n". 
        "--canonical_only [use this flag to only print consequences from canonical transcripts when using --vep flag]\n". 
        "--do_not_simplify [use this flag to output all fields to Excel as they appear in the original VCF. Useful if you have a customised VCF but doesn't allow you to pick out individual samples]\n". 
        "--text_output [use this flag to output tab-delimited text instead of Excel format]\n". 
        "--quiet [use this flag to supress warning messages if AD, GQ or PL fields can't be found in FORMAT field for a variant]\n". 
        "--gene_anno [use this flag if annotated with ensembl_gene_annotator.pl]\n--help [use this flag to print this help message]\n\n";
        exit;
    }else{
        die "--input argument is required\n";
    }
}
############

my $vcf;
my @samples;
my @classes;
my @add;
my $output;
my $config = {};
my $help;
my $quiet;
GetOptions($config,
    'classes=s{,}' => \@classes,
    'additional_classes=s{,}' => \@add,
    'functional',
    'gene_anno',
    'text_output',
    'do_not_simplify',
    'quiet', => \$quiet,
    'canonical_only',
    'vep',
    'fields=s{,}', => @{$config->{fields}}, 
    'input=s' =>\$vcf,
    'output=s' => \$output,
    'samples=s{,}' => \@samples,
    'help' => \$help) or die "Syntax error.\n";
usage($help) if $help;
usage() if (not $vcf);

if (@classes or @add){
    $config->{functional} = 1;
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
foreach my $class (@classes){
        die "Error - variant class '$class' not recognised.\n" if not grep {/$class/i} @valid;
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
    if (@fields){
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
        print $OUT join("\t", ("Chrom", "Pos", "SNP_ID", "Variant Percent Confidence", "Filters", "Genomic Ref", ));
        if (not @samples){
            @samples = $vcf_obj->getSampleNames();
        }else{
            foreach my $sample (@samples){
                die "Can't find sample $sample in header line.\nHeader:\n$header_string\n" if not grep{/^$sample$/} $vcf_obj->getSampleNames();
            }
        }
        print $OUT join("\t", @samples) ."\t" ;
        print $OUT join(" Allele Depth\t", @samples) ." Allele Depth\t" ;
        print $OUT join(" Genotype Confidence\t", @samples) ." Genotype Confidence\t" ;
        print $OUT join(" Scaled Percent Probably Liklihoods\t", @samples) ." Scaled Percent Probably Liklihoods" ;
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
    if (@fields){
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
    push @output, 100 - (100 * (10**(-($vcf_obj->getVariantField("QUAL"))/10)));
    push @output, $vcf_obj->getVariantField("FILTER");
    push @output, $vcf_obj->getVariantField("REF");
    my @form = split(":", $vcf_obj->getVariantField("FORMAT"));
    my @all_alleles = ($vcf_obj->getVariantField("REF"), split(",", $vcf_obj->getVariantField("ALT")));
    my $gt = 0;
    my $ad = 0;
    my $gq = 0;
    my $pl = 0;
    {
    no warnings 'uninitialized';
    $gt++ until $form[$gt] eq "GT" or $gt > $#form;
    $ad++ until $form[$ad] eq "AD" or $ad > $#form;
    $gq++ until $form[$gq] eq "GQ" or $gq > $#form;
    $pl++ until $form[$pl] eq "PL" or $pl > $#form;
    }
    die "Can't find genotype field in variant line:\n" . $vcf_obj->get_currentLine() ."\n" if ($gt > $#form);
    print STDERR "Can't find allele depth field in variant line:\n" . $vcf_obj->get_currentLine() ."\n" if ($ad > $#form and not $quiet);
    print STDERR "Can't find genotype quality field in variant line:\n" . $vcf_obj->get_currentLine() ."\n" if ($gq > $#form and not $quiet);
    print STDERR "WARNING - Can't find probable likelihoods quality field in variant line:\n" . $vcf_obj->get_currentLine() ."\n" if ($pl > $#form and not $quiet);
    my @sample_calls = ();
    my @sample_allele_depths = ();
    my @sample_genotype_quality = ();
    my @sample_probable = ();
    foreach my $sample (@samples){
        my $genotype;
        my @call = split(":", $vcf_obj->getSampleVariant($sample));
        if (@call > $gt and $call[$gt] =~ /(\d+)[\/\|](\d+)/){
            $genotype = "$all_alleles[$1]/$all_alleles[$2]";
        }else{
            $genotype = "-";
        }
        push (@sample_calls, $genotype);
        my $allele_depths;
        if (@call > $ad){
            if ($call[$ad] =~ /(\d+,*)+/){
                my @ads = split(",", $call[$ad]);
                my @allele_depth = ();
                for (my $n = 0; $n < @ads; $n++){
                    push (@allele_depth, "$all_alleles[$n]=$ads[$n]");
                }
                $allele_depths = join("/", @allele_depth);
            }else{
                $allele_depths = "-";
            }
        }else{
            $allele_depths = "-";
        }
        push (@sample_allele_depths, $allele_depths);
        my $genotype_quality;
        if (@call > $gq){
            if ($call[$gq] =~ /(\d+\.*\d*)+/){
                $genotype_quality = 100 - (100 * 10**(-$1/10));
            }else{
                $genotype_quality = "-";
            }
        }else{
            $genotype_quality = "-";
        }
        push (@sample_genotype_quality, $genotype_quality);
        my $probable_likelihoods;
        if (@call > $pl){
            if ($call[$pl] =~ /(\d+,\d+,\d+)/){
                my @prob = split(",", $call[$pl]);
                my @possible_genotypes = get_possible_genotypes(\@all_alleles);
                if (scalar@possible_genotypes != scalar@prob){
                    print STDERR "ERROR - number of possible genotypes doesn't match genotype probabilities for line:\n" . $vcf_obj->get_currentLine() ."\n";
                    $probable_likelihoods = "-";
                }else{    
                    my @probables = ();
                    for (my $n = 0; $n < @prob; $n++){
                        my $percent = 100 * (10**(-$prob[$n]/10));
                        push (@probables, "$possible_genotypes[$n]=$percent");
                    }
                    $probable_likelihoods =  join(", ", @probables);
                }
            }else{
                $probable_likelihoods = "-";
            }
        }else{
            $probable_likelihoods = "-";
        }
        push (@sample_probable, $probable_likelihoods);
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
    foreach my $probable_percent (@sample_probable){
        push @output, $probable_percent;
    }
    if (@fields){
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
=cut$worksheet->set_column( 0, 0, 10);
    $worksheet->set_column( 1, 1, 12);
    $worksheet->set_column( 2, 2, 25);
    $worksheet->set_column( 4, 4, 9);
    $worksheet->set_column( 10, 10, 11);
    $worksheet->set_column( 11, 11, 13);
    $worksheet->set_column( 12, 15, 20);
=cut
    my $vep_header;
    if (@fields){
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

        if (not @samples){
            @samples = $vcf_obj->getSampleNames();
        }else{
            foreach my $sample (@samples){
                die "Can't find sample $sample in header line.\nHeader:\n$header_string\n" if not grep{/^$sample$/} $vcf_obj->getSampleNames();
            }
        }
        foreach my $sample (@samples){
            $worksheet->write($row, $col++, $sample, $header_formatting);
        }

        foreach my $sample (@samples){
            $worksheet->write($row, $col++, "$sample Allele Depth", $header_formatting);
        }
        foreach my $sample (@samples){
            $worksheet->write($row, $col++, "$sample Genotype Confidence", $header_formatting);
        }
        foreach my $sample (@samples){
            $worksheet->write($row, $col++, "$sample Scaled Percent Probably Liklihoods", $header_formatting);
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
        if ($config->{samples}){
            next if not exists $vep_allele{$annot->{allele}};
        }
        $temp_col = $col;
        foreach my $vep  (@fields){
            my $value = $annot->{$vep};
            $worksheet->write($row, $temp_col++, $value);
        }
        $lines++;
        $row++;
    }   
    return ($lines, $temp_col);
}
####################################################
sub write_line_to_excel{
#return no. of merged rows in order to increment $row
    my ($row, $col, ) = @_;
    my $lines = 1;
    if (@fields){
        ($lines, $col) = write_vep_fields($row, $col);
    }
    if ($lines < 1){
    #this happens if --canonical_only argument is in effect and we haven't got a canonical transcript for this variant
    #we won't have printed anything for the consequence fields
        $worksheet->write($row, $col++, "No canonical variant found for this site");
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
    if (@fields){
        ($lines, $col) = write_vep_fields($row, $col);
    }
    if ($lines < 1){
    #this happens if --canonical_only argument is in effect and we haven't got a canonical transcript for this variant
    #we won't have printed anything for the consequence fields
        $worksheet->write($row, $col++, "No canonical variant found for this site");
        for (my $i = 1; $i < @fields; $i++){
            $worksheet->write($row, $col++, "-");
        }
        $lines = 1;
    }
    foreach my $varfield ($vcf_obj->getVariantField("CHROM"), $vcf_obj->getVariantField("POS"), $vcf_obj->getVariantField("ID"),){ 
        #$worksheet->write($row, $col++, $varfield);
        #push (@shared_values, $varfield);
        write_worksheet($lines, $varfield, $worksheet, $row, $col++, $std_formatting);
    }
    #push (@shared_values, 100 - (100 * (10**(-$split_line[$qual]/10))), $split_line[$header_hash->{FILTER}], $split_line[$header_hash->{REF}]);
    write_worksheet($lines, 100 - (100 * (10**(-($vcf_obj->getVariantField("QUAL"))/10))), $worksheet, $row, $col++, $std_formatting, );
    write_worksheet($lines, $vcf_obj->getVariantField("FILTER"), $worksheet, $row, $col++, $std_formatting);
    write_worksheet($lines, $vcf_obj->getVariantField("REF"), $worksheet, $row, $col++, $std_formatting);
    #$worksheet->write($row, $col++, 100 - (100 * (10**(-$line[$qual]/10))));
    #$worksheet->write($row, $col++, $line[$filter]);
    #$worksheet->write($row, $col++, $line[$ref]);
    my @form = split(":", $vcf_obj->getVariantField("FORMAT"));
    my @all_alleles = ($vcf_obj->getVariantField("REF"), split(",", $vcf_obj->getVariantField("ALT")));
    my $gt = 0;
    my $ad = 0;
    my $gq = 0;
    my $pl = 0;
    {
    no warnings 'uninitialized';
    $gt++ until $form[$gt] eq "GT" or $gt > $#form;
    $ad++ until $form[$ad] eq "AD" or $ad > $#form;
    $gq++ until $form[$gq] eq "GQ" or $gq > $#form;
    $pl++ until $form[$pl] eq "PL" or $pl > $#form;
    }
    die "Can't find genotype field in variant line:\n" . $vcf_obj->get_currentLine() ."\n" if ($gt > $#form);
    print STDERR "Can't find allele depth field in variant line:\n" . $vcf_obj->get_currentLine() ."\n" if ($ad > $#form and not $quiet);
    print STDERR "Can't find genotype quality field in variant line:\n" . $vcf_obj->get_currentLine() ."\n" if ($gq > $#form and not $quiet);
    print STDERR "WARNING - Can't find probable likelihoods quality field in variant line:\n" . $vcf_obj->get_currentLine() ."\n" if ($pl > $#form and not $quiet);
    my @sample_calls = ();
    my @sample_allele_depths = ();
    my @sample_genotype_quality = ();
    my @sample_probable = ();
    foreach my $sample (@samples){
        my $genotype;
        my @call = split(":", $vcf_obj->getSampleVariant($sample));
        if (@call > $gt and $call[$gt] =~ /(\d+)[\/\|](\d+)/){
            $genotype = "$all_alleles[$1]/$all_alleles[$2]";
        }else{
            $genotype = "-";
        }
        push (@sample_calls, $genotype);
        my $allele_depths;
        if (@call > $ad){
            if ($call[$ad] =~ /(\d+,*)+/){
                my @ads = split(",", $call[$ad]);
                my @allele_depth = ();
                for (my $n = 0; $n < @ads; $n++){
                    push (@allele_depth, "$all_alleles[$n]=$ads[$n]");
                }
                $allele_depths = join("/", @allele_depth);
            }else{
                $allele_depths = "-";
            }
        }else{
            $allele_depths = "-";
        }
        push (@sample_allele_depths, $allele_depths);
        my $genotype_quality;
        if (@call > $gq){
            if ($call[$gq] =~ /(\d+\.*\d*)+/){
                $genotype_quality = 100 - (100 * 10**(-$1/10));
            }else{
                $genotype_quality = "-";
            }
        }else{
            $genotype_quality = "-";
        }
        push (@sample_genotype_quality, $genotype_quality);
        my $probable_likelihoods;
        if (@call > $pl){
            if ($call[$pl] =~ /(\d+,\d+,\d+)/){
                my @prob = split(",", $call[$pl]);
                my @possible_genotypes = get_possible_genotypes(\@all_alleles);
                if (scalar@possible_genotypes != scalar@prob){
                    print STDERR "ERROR - number of possible genotypes doesn't match genotype probabilities for line:\n" . $vcf_obj->get_currentLine() ."\n";
                    $probable_likelihoods = "-";
                }else{    
                    my @probables = ();
                    for (my $n = 0; $n < @prob; $n++){
                        my $percent = 100 * (10**(-$prob[$n]/10));
                        push (@probables, "$possible_genotypes[$n]=$percent");
                    }
                    $probable_likelihoods =  join(", ", @probables);
                }
            }else{
                $probable_likelihoods = "-";
            }
        }else{
            $probable_likelihoods = "-";
        }
        push (@sample_probable, $probable_likelihoods);
    }
    foreach my $call (@sample_calls){
        #$worksheet->write($row, $col++, $call);
        write_worksheet($lines, $call, $worksheet, $row, $col++, $std_formatting);
        #push (@shared_values, $call);
    }
    foreach my $sample_depth (@sample_allele_depths){
        #$worksheet->write($row, $col++,$sample_depth);
        write_worksheet($lines, $sample_depth, $worksheet, $row, $col++, $std_formatting);
        #push (@shared_values, $sample_depth);
    }
    foreach my $sample_gq (@sample_genotype_quality){
        #$worksheet->write($row, $col++,$sample_gq);
        write_worksheet($lines, $sample_gq, $worksheet, $row, $col++, $std_formatting, );
        #push (@shared_values, $sample_gq);
    }
    foreach my $probable_percent (@sample_probable){
        #$worksheet->write($row, $col++,$probable_percent);
        write_worksheet($lines, $probable_percent, $worksheet, $row, $col++, $std_formatting);
        #push (@shared_values, $probable_percent);
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
