#!/usr/bin/env perl
#David Parry August 2011
#
=head1 NAME

annovcfToSimple.pl - takes a vcf file and outputs a simplified version in Excel's .xlsx format.

=head1 SYNOPSIS

        annovcfToSimple.pl -i [ file] [options]
        
        annovcfToSimple.pl -s -i [SnpEff annotated VCF file] [options]
        
        annovcfToSimple.pl -v -i [VEP annotated VCF file] [options]
        
        annovcfToSimple.pl -g -v -i [VEP and geneAnnotator.pl annotated VCF file] [options]
        
        annovcfToSimple.pl -s -i [SnpEff annotated VCF file] [options]

        annovcfToSimple.pl -h (display help message)
        
        annovcfToSimple.pl -m (display manual page)

=cut

=head1 ARGUMENTS

=over 8

=item B<-i    --input>

Input VCF file, optionally annotated with variant_effect_predictor.pl/SnpEff and geneAnnotator.pl. 

=item B<-o    --output>

Output file name. Defaults to [input name].xlsx (or [input name].txt if --text_output option is used).

=item B<-u    --summarise>

Use this flag to summarise the number of alleles and genotypes found in samples rather than outputting genotype columns for each individual sample. If sample IDs are specified with --samples or --pedigree options only these samples will be counted in the summarised allele and genotype counts. This can be overriden by adding the word 'all' after this argument, in which case all samples in the VCF will be summarised and any samples specified by --samples or --pedigree options will still have their individual genotypes written to the output. Alternatively you may specify the word 'and' after this argument if you are not using the --samples argument in order to print both a summary of all samples in addition to columns for the individual genotypes, allele depths etc. for every sample in the VCF.

=item B<--contains_variant>

Use this flag alongside the -u/--summarise option to add a column summarising which samples contain variant genotypes, plus two additional columns (ADs and GQs) to give allele depths and genotype quality scores for samples with variants.

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

=item B<--snpeff>

Use this flag to output annotations from SnpEff at the beginning of each line (assuming your input has already been annotated by SnpEff). Can not be used in conjunction with --vep option.

=item B<-g    --gene_anno>

Use this flag if annotated with geneAnnotator.pl to output gene annotation information. This script cannot distinguish between gene annotation information belonging to canonical transcripts or particular variant classes, these distinctions have to be made when running geneAnnotator.pl prior to this script.

=item B<--canonical_only>

Use this flag to only print consequences from canonical transcripts when using --vep option. 

=item B<-f    --functional>

Use this flag to only annotate standard 'functional' variant classes. VEP/SnpEff results not matching these classes will not be reported. Default 'functional' variant classes are:

B<VEP:>

        TFBS_ablation
        TFBS_amplification
        frameshift_variant
        inframe_deletion
        inframe_insertion
        initiator_codon_variant
        missense_variant
        protein_altering_variant
        regulatory_region_ablation
        regulatory_region_amplification
        splice_acceptor_variant
        splice_donor_variant
        splice_region_variant
        stop_gained
        stop_lost
        transcript_ablation
        transcript_amplification

B<SnpEff:>

        chromosome
        coding_sequence_variant
        inframe_insertion
        disruptive_inframe_insertion
        inframe_deletion
        disruptive_inframe_deletion
        exon_loss_variant
        frameshift_variant
        missense_variant
        rare_amino_acid_variant
        splice_acceptor_variant
        splice_donor_variant
        splice_region_variant
        stop_lost
        5_prime_UTR_premature_start_codon_gain_variant
        start_lost
        stop_gained
        exon_loss

Available classes that can be chosen instead of (or in addition to - see below) these classes can be found in the data/vep_classes.tsv and data/snpeff_classes.tsv files respectively using the -c/--classes. 



=item B<-c    --classes>

Use this to specify a set of VEP variant classes to print in output. Overrides --functional option.

=item B<-a    --additional_classes>

Use this to specify additional VEP variant classes to output as well as those used by --functional.

=item B<-e    --fields>

Specify one or more VEP/SnpEff fields to output. Use of --vep or --snpeff options without this flag outputs the following fields (assuming they are present in your VCF):

B<VEP:>
       
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
            eas_maf
            eur_maf
            sas_maf
B<SnpEff:>
            allele
            annotation
            annotation_impact
            gene_name
            gene_id
            feature_type
            feature_id
            transcript_biotype
            hgvs.c
            hgvs.p
            rank
            cds.pos / cds.length
            aa.pos / aa.length

If you wish to use the above fields in addition to fields specified here add 'default' to your list of fields.
       
=item B<--all>

Use with --vep or --snpeff option to output all available VEP/SnpEff fields.

=item B<-n    --info_fields>

One or more INFO field IDs from to output as columns. These are case sensitive and must appear exactly as defined in the VCF header. By default, if --info_fields is not specified, CaddPhredScore INFO fields will be written to the output automatically if found.

=item B<-b    --biotype_filters>

When using the -f/--functional or -c/--classes options, features/transcripts with the following biotypes are ignored by default:

    IG_C_pseudogene
    IG_J_pseudogene
    IG_V_pseudogene
    nonsense_mediated_decay
    non_stop_decay
    processed_pseudogene
    pseudogene
    retained_intron
    transcribed_processed_pseudogene
    transcribed_unprocessed_pseudogene
    TR_J_pseudogene
    TR_V_pseudogene
    unitary_pseudogene
    unprocessed_pseudogene

You may override this filter by specifying biotypes to filter with this option or prevent biotype filtering using the --no_biotype_filtering option (see below). 

The 'data/biotypes.tsv' file contains a list of valid biotypes and the default behaviour of this program (i.e. 'keep' or 'filter'). 

=item B<--no_biotype_filtering>

Use this flag to consider consequences affecting ALL biotypes when using -f/--functional or -c/--classes options.


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

Reads a VCF file and outputs a simplified version in Excel's .xlsx format. Useful when using a VCF annotated with Ensembl's variant_effect_predictor.pl or SnpEff to output a more easily readable version (use the --vep or --snpeff option), and particularly useful for VCF files annotated with geneAnnotator.pl (use the -g option).

=cut

=head1 AUTHOR

David A. Parry


=head1 COPYRIGHT AND LICENSE

Copyright 2013, 2014, 2015  David A. Parry

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
use FindBin qw($RealBin);
use lib "$RealBin/lib";
use VcfReader;
use VcfhacksUtils;
use ParsePedfile;
use TextToExcel;


############

my $vcf;
my @samples;
my @csq_classes;
my @csq_add;
my @biotypes = ();
my @peds;
my $output;
my $config = {};
my $summarise_counts;
my $summarise_samples_with_variants;
my $help;
my $man;
GetOptions($config,
    'all',
    'a|additional_classes=s{,}' => \@csq_add,
    'b|biotypes=s{,}' => \@biotypes,
    'canonical_only',
    'contains_variant' => \$summarise_samples_with_variants,
    'c|classes=s{,}' => \@csq_classes,
    'do_not_simplify',
    'fields|e=s{,}', => \@{$config->{fields}}, 
    'functional|f',
    'gene_anno',
    'i|input=s' =>\$vcf,
    'manual' => \$man,
    'no_biotype_filtering',
    'n|info_fields=s{,}', => \@{$config->{info_fields}},
    'o|output=s' => \$output,
    'pedigree=s{,}' => \@peds,
    'snpeff',
    's|samples=s{,}' => \@samples,
    'text_output',
    'u|summarise:s' => \$summarise_counts,
    'vep',
    'help' => \$help, 
) or pod2usage(-exitval => 2, -message => "Syntax error.\n");
pod2usage( -verbose => 2 ) if $man;
pod2usage( -verbose => 1 ) if $help;
pod2usage( -exitval => 2, -message => "--input is required" ) if (not $vcf);

if (@csq_classes or @csq_add){
    $config->{functional} = 1;
}

#get and check header
my @header = VcfReader::getHeader($vcf);
die "ERROR: Invalid VCF header in $vcf\n" 
  if not VcfReader::checkHeader(header => \@header);

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
my @info_fields = ();
my @gene_anno_fields = ();

my %info_fields = VcfReader::getInfoFields(header => \@header);
if (defined $config->{gene_anno}){
    if (not exists $info_fields{GeneAnno}){
        die "GeneAnno INFO field not found in VCF header - please annotate ".
            "your VCF with geneAnnotator.pl or rerun this program without ".
            "the --gene_anno option.\n";
    }
    if ($info_fields{GeneAnno}->{Description} =~ /Format: ([\w+\|]+)\"/){
        @gene_anno_fields = split(/\|/, $1);
    }else{
        die "Could not parse GeneAnno description field - please annotate ".
            "your VCF with geneAnnotator.pl or report this error if it recurs.\n";
    }
}


my %csq_header = getAndCheckCsqHeader();
my @csq_fields = getCsqFields();#consequence fields to retrieve from VCF
@csq_classes = getAndCheckClasses();
my %class_filters = map { $_ => undef } @csq_classes;
my %biotype_filters = map { $_ => undef } getAndCheckBiotypes();


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


my %sample_to_col = VcfReader::getSamples
(
    header => \@header,
    get_columns => 1,
);
if (@ped_samples){
    my @sample_found = ();
    foreach my $s (@ped_samples){
        if (exists $sample_to_col{$s}){
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
    if (not defined $summarise_counts or $summarise_counts ne 'all'){
        @samples = VcfReader::getSamples
        (
            header => \@header,
        );
    }
}else{
    foreach my $s (@samples){
        die "Can't find sample $s in header line.\nHeader:\n$header[-1]\n" 
            if not exists $sample_to_col{$s};
    }
}

my $OUT;
my $xl_obj;
my @col_header = get_header();
if (defined $config->{text_output}){
    if ($output){
        open ($OUT, ">$output") or die "Can't open $output: $!\n";
    }else{
        $OUT = \*STDOUT;
    }
    foreach my $h (@col_header){
        $h =~ s/ /_/g;
    }
    print $OUT join("\t", @col_header) ."\n";
}else{
    $xl_obj = TextToExcel->new( file=> $output);
    $header_formatting = $xl_obj->createFormat(bold => 1);
    $url_format = $xl_obj->createFormat(
        color     => 'blue',
        underline => 1,
    );
    $xl_obj->writeLine(line => \@col_header, format => $header_formatting);
}

my $VCF = VcfReader::openVcf($vcf); 
LINE: while (my $var = <$VCF>){
    next if $var =~ /^#/;
    chomp $var;
    my @split = split("\t", $var);
    my @preceding = ();
    my @line = ();
    my @following = ();
    #Get preceding VEP columns (there may be multiple rows per line)
    if (defined $config->{vep} or defined $config->{snpeff}){
        my ($prec, $canonical_found, $functional_found, $sample_found) 
            = get_csq_fields(\@split); 
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
            for (my $i = 1; $i < @csq_fields; $i++){
                push @temp_prec, "-";
            }
            push @preceding, \@temp_prec;
        }

    }
    if (defined $config->{gene_anno}){
        my $g_anno = VcfReader::getVariantInfoField
        (
            \@split,
            'GeneAnno'
        );
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
        @line = get_vcf_fields(\@split);    
    }else{
        @line = get_simplified_fields(\@split);    
    }
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
        my $row = $xl_obj->writeLine(
                    line => \@line, 
                    preceding => \@preceding, 
                    succeeding => \@following, 
                    );
    }
}
close $VCF;
if ($xl_obj){
    $xl_obj->DESTROY();
}
if ($OUT){
    close $OUT;
}

#################################################
sub getAndCheckClasses{
    return if not $config->{vep} and not $config->{snpeff}; 
    my %all_classes = ();
    if ($config->{vep} ){
        %all_classes =  VcfhacksUtils::readVepClassesFile();
    }else{
        %all_classes =  VcfhacksUtils::readSnpEffClassesFile();
    }
    if (not @csq_classes){
        @csq_classes = grep { $all_classes{$_} eq 'default' } keys %all_classes;
        push @csq_classes, 'splice_region_variant';
    }
    push @csq_classes, @csq_add if (@csq_add);
    @csq_classes = map { lc($_) } @csq_classes; 
    @csq_classes = VcfhacksUtils::removeDups(@csq_classes);
    foreach my $class (@csq_classes) {
        die "Error - variant class '$class' not recognised.\n"
          if not exists $all_classes{lc($class)} ;
    }
    return @csq_classes;
}

####################################################
sub getAndCheckCsqHeader{
    return if not $config->{vep} and not $config->{snpeff}; 
    if ($config->{vep} and $config->{snpeff}){ 
        die "Please supply only one of --vep or --snpeff options.\n";
    }
    my %csq_head = ();
    if ($config->{vep}){
        eval { 
            %csq_head = VcfReader::readVepHeader
            (
                header => \@header
            ); 
        };
        if ($@){
            die "ERROR: Could not find VEP header in input. Please ".
                "annotate your input with Ensembl's VEP and try again.\n";
        } ;
    }elsif ($config->{snpeff}){
        eval { 
            %csq_head = VcfReader::readSnpEffHeader
            (
                header => \@header
            ); 
        } ;
        if ($@){
            die "ERROR: Could not find SnpEff header in input. Please ".
                "annotate your input with SnpEff and try again.\n";
        }
    }
    return %csq_head;
}

####################################################
sub getCsqFields{
    return if not $config->{vep} and not $config->{snpeff}; 
    my @default_fields = ();
    if ($config->{vep}){
       @default_fields =  
        qw(
            allele
            gene
            feature
            feature_type
            consequence
            symbol
            biotype
        );
        foreach my $f (qw / 
            amino_acids
            codons
            existing_variation
            exon
            intron
            condel
            splice_consensus
            hgvsc
            hgvsp
            cds_position
            protein_position
            sift
            polyphen
            condel
            gmaf
            aa_maf
            ea_maf
            afr_maf
            amr_maf
            eas_maf
            eur_maf
            sas_maf
            lof
            lof_filter
            lof_flags
            lof_info
        /){
            push @default_fields, $f 
                if exists $csq_header{$f};
        }
    }else{
        @default_fields = 
        qw(
            allele
            annotation
            annotation_impact
            gene_name
            gene_id
            feature_type
            feature_id
            transcript_biotype
            hgvs.c
            hgvs.p
            rank
        );
        push @default_fields, 
            "cds.pos / cds.length", 
            "aa.pos / aa.length",
            "errors / warnings / info";
    }

    if (not @{$config->{fields}} and not $config->{all}){
        unshift @{$config->{fields}}, @default_fields;
    }elsif($config->{all}){
         push @{$config->{fields}},  
            sort {$csq_header{$a} <=> $csq_header{$b}} keys %csq_header;
    }elsif(grep { /default/i } @{$config->{fields}}){
        @{$config->{fields}} = grep {! /^default$/i } @{$config->{fields}};
        unshift @{$config->{fields}}, @default_fields;
    }
    my @required = ( "allele" );
    if ($config->{vep}){
        push @required, "consequence";
        push @required, "biotype";
    }else{
        push @required, "annotation";
        push @required, "transcript_biotype";
    }
    foreach my $req (@required){
        if (not grep { $_ eq $req } @{$config->{fields}} ){
            push @{$config->{fields}}, $req;
        }
    }
    @csq_fields = map { lc($_) } @{$config->{fields}} if defined $config->{fields};
    @csq_fields = VcfhacksUtils::removeDups(@csq_fields);
    if($config->{vep} and $config->{canonical_only}){
        push @csq_fields, 'canonical' if not grep{$_ eq 'canonical'} @csq_fields;
    }
    my @temp_fields = (); 
    foreach my $csq (@csq_fields){
        if (not exists $csq_header{$csq}){
            if ( grep {$_ eq $csq} @required ){
                die "Couldn't find '$csq' consequence field in header - ".
                    "please ensure your VCF is annotated with " .
                    "the appropriate annotations.\n";
            }else{
                warn "Couldn't find '$csq' consequence field in header - will not output this annotaiton.\n";
            }
        }else{
            push @temp_fields, $csq;
        }
    }
    return @temp_fields;
}

####################################################
sub get_vcf_fields{
    my $line = shift;
    my @vcf_values = ();
    foreach my $field (@$line){
        push @vcf_values, $field;
    }
    #Add user-specified INFO fields after normal VCF fields
    foreach my $f (@info_fields){
        my $value = VcfReader::getVariantInfoField
        (   
            $line,
            $f
        );
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
    my $line = shift;
    my @vcf_values = ();
    #Get preceding INFO fields (these will be one per line)
    push @vcf_values, VcfReader::getMultipleVariantFields
    (
        $line, 
        qw (
            CHROM 
            POS 
            ID 
            QUAL 
            FILTER 
            REF 
            ALT
        )
    ); 
    my @all_alleles = VcfReader::readAlleles(line => $line);
    if (defined $summarise_counts){
        my %allele_counts = ();
        my %genotype_counts = ();
        if ($summarise_counts eq 'all' or $summarise_counts eq 'and'){
            %allele_counts = VcfReader::countAlleles(line => $line);
            %genotype_counts = VcfReader::countGenotypes(line => $line);
        }else{
            %allele_counts = VcfReader::countAlleles
            (
                line => $line,
                samples => \@samples,
                sample_to_columns => \%sample_to_col,
            );
            %genotype_counts = VcfReader::countGenotypes
            (
                line => $line,
                samples => \@samples,
                sample_to_columns => \%sample_to_col,
            );
        }
            
        my @al_count_string = ();
        my @gt_count_string = ();
        my @sample_var_string = ();
        my @sample_ad_string = ();
        my @sample_gq_string = ();
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
                %calls = VcfReader::getSampleCall
                (
                    line => $line, 
                    all => 1,
                    sample_to_columns => \%sample_to_col,
                );
            }else{
                %calls = VcfReader::getSampleCall
                (
                    line => $line, 
                    multiple => \@samples,
                    sample_to_columns => \%sample_to_col,
                );
            }
            foreach my $sample (sort keys %calls){
                next if ($calls{$sample} =~/^0[\|\/]0$/);
                next if ($calls{$sample} =~/^\.[\|\/]\.$/);
                push @{$var_to_samples{$calls{$sample}}}, $sample;
                my @ads = VcfReader::getSampleAlleleDepths 
                (
                    sample =>$sample,
                    line => $line,
                    sample_to_columns => \%sample_to_col,
                );
                my @allele_depth = ();
                for (my $n = 0; $n < @ads and $n < @all_alleles; $n++){
                    push (@allele_depth, "$all_alleles[$n]=$ads[$n]");
                }
                my $allele_depths = join("/", @allele_depth);
                push @sample_ad_string, "$sample:$allele_depths";
                my $gq = VcfReader::getSampleGenotypeField
                (
                    sample  => $sample, 
                    line => $line,
                    sample_to_columns => \%sample_to_col,
                    field => "GQ",
                );
                push @sample_gq_string, "$sample=$gq";
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
        if (defined $summarise_samples_with_variants){
            push @vcf_values, join(";", @sample_var_string);
            push @vcf_values, join(";", @sample_ad_string);
            push @vcf_values, join(";", @sample_gq_string);
        }
    }
    if (not defined $summarise_counts 
        or ($summarise_counts eq 'all' and @samples)
        or $summarise_counts eq 'and'
    ){
        my @sample_calls = ();
        my @sample_allele_depths = ();
        my @sample_genotype_quality = ();
        foreach my $sample (@samples){
            my $genotype = VcfReader::getSampleActualGenotypes 
            (
                sample =>$sample,
                line => $line,
                sample_to_columns => \%sample_to_col,
            );
            if (defined $genotype){ 
                push (@sample_calls, $genotype);
            }else{
                push @sample_calls, "-";
            }
            my @ads = VcfReader::getSampleAlleleDepths 
            (
                sample =>$sample,
                line => $line,
                sample_to_columns => \%sample_to_col,
            );
            my $allele_depths = '-';
            if (@ads){
                my @allele_depth = ();
                for (my $n = 0; $n < @ads and $n < @all_alleles; $n++){
                    push (@allele_depth, "$all_alleles[$n]=$ads[$n]");
                }
                $allele_depths = join("/", @allele_depth);
            }
            push (@sample_allele_depths, $allele_depths);
            my $genotype_quality = "-";
            my $gq = VcfReader::getSampleGenotypeField
            (
                sample  => $sample, 
                line => $line,
                sample_to_columns => \%sample_to_col,
                field => "GQ",
            );
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
        my $value = VcfReader::getVariantInfoField
        (
            $line,
            $f,
        ); 
        if (defined $value){
            push @vcf_values, $value;
        }else{
            push @vcf_values, '.';
        }
    }
    return @vcf_values;
}

#################################################
sub getAndCheckBiotypes{
    return if $config->{no_biotype_filtering}; 
    my %all_biotypes = VcfhacksUtils::readBiotypesFile();
    if (not @biotypes){
        @biotypes = grep { $all_biotypes{$_} eq 'filter' } keys %all_biotypes;
    }
    @biotypes = map { lc($_) } @biotypes; 
    @biotypes = VcfhacksUtils::removeDups(@biotypes);
    foreach my $biotyp (@biotypes) {
        die "Error - biotype '$biotyp' not recognised.\n"
          if not exists $all_biotypes{lc($biotyp)} ;
    }
    return @biotypes;
}
 
####################################################
sub get_csq_fields{
    my $line = shift;
    if ($config->{vep}){
        return get_vep_fields($line);
    }else{
        return get_snpeff_fields($line);
    }
}

####################################################
sub get_snpeff_fields{
    #returns a 2D array of VEP values
    my $line = shift;
    my ($row, $col) = @_;
    my @csq = VcfReader::getSnpEffFields
    ( 
        line          => $line,
        field         => \@csq_fields,
        snpeff_header => \%csq_header,
    );
    my ($canonical_found, $functional_found, $sample_found) = (0, 0, 0);
    my @snpeff_values = ();
    my %sample_alleles = ();
    if ($config->{samples}){
        %sample_alleles = map {$_ => undef} 
        VcfReader::getSampleActualGenotypes
        (
            line                => $line,
            multiple            => $config->{samples}, 
            return_alleles_only => 1
        );
    }
ANNOT:  foreach my $annot (@csq){
        my @csq_values = ();
        $canonical_found++;#no canonical check for snpeff
        if ($config->{functional}){
            next ANNOT if not check_functional($annot);
        }
        $functional_found++;
        if ($config->{samples}){
            next if not exists $sample_alleles{$annot->{allele}};
        }
        $sample_found++;
        foreach my $csq  (@csq_fields){
            my $value = $annot->{$csq};
            if (defined $value){
                push @csq_values, $value;
            }else{
                push @csq_values, '-';
            }
        }
        push @snpeff_values, \@csq_values;
    }
    return (\@snpeff_values, $canonical_found, $functional_found, $sample_found);
}
####################################################
sub get_vep_fields{
    #returns a 2D array of VEP values
    my $line = shift;
    my ($row, $col) = @_;
    my @csq = VcfReader::getVepFields
    ( 
        line        => $line,
        field       => \@csq_fields,
        vep_header  => \%csq_header,
    );
    my ($canonical_found, $functional_found, $sample_found) = (0, 0, 0);
    my %vep_allele = ();
    my @vep_values = ();
    my @all_alleles = VcfReader::readAlleles(line => $line);
    if ($config->{samples}){
        my @sample_alleles = VcfReader::getSampleActualGenotypes
        (
            line                => $line,
            multiple            => $config->{samples}, 
            return_alleles_only => 1
        );
        my @v_alleles = VcfReader::altsToVepAllele
        (
            line => $line,
        );
        @vep_allele{@v_alleles} = @sample_alleles;
     }
ANNOT:  foreach my $annot (@csq){
        my @csq_values = ();
        if ($config->{canonical_only}){
              next if (not $annot->{'canonical'} );
        }
        $canonical_found++;
        if ($config->{functional}){
            next ANNOT if not check_functional($annot);
        }
        $functional_found++;
        if ($config->{samples}){
            next if not exists $vep_allele{$annot->{allele}};
        }
        $sample_found++;
        foreach my $vep  (@csq_fields){
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
sub check_functional{
    my $annot = shift;
    if ($config->{vep}){
        return check_vep_consequence($annot);
    }else{
        return check_snpeff_consequence($annot);
    }
}

########################################
sub check_snpeff_consequence {
    my $annot = shift;
    #skip variants with undef features (intergenic variants)
    return 0 if not defined $annot->{feature_id};
    #skip unwanted biotypes
    return 0 if exists $biotype_filters{lc $annot->{transcript_biotype} };
    my @anno_csq = split( /\&/, $annot->{annotation} );
ANNO: foreach my $ac (@anno_csq){
        $ac = lc($ac);#we've already converted %class_filters to all lowercase
        if ( exists $class_filters{$ac} ){
            return 1;
        }
    }
    return 0;#no annotation matching %class_filters
}

########################################
sub check_vep_consequence {
    my $annot = shift;
    #intergenic variants have no feature associated with them - skip
    return 0 if $annot->{consequence} eq "intergenic_variant";
    #skip unwanted biotypes
    return 0 if (exists $biotype_filters{$annot->{biotype}}) ;
    #skip non-canonical transcripts if --canonical_only selected
    if ($config->{canonical_only}) {
        return 0 if ( not $annot->{canonical} );
    }
    
    my @anno_csq = split( /\&/, $annot->{consequence} );

ANNO: foreach my $ac (@anno_csq){
        $ac = lc($ac);#we've already converted %class_filters to all lowercase
        if ( exists $class_filters{$ac} ){
            return 1;
        }
    }
    return 0;#no annotation matching %class_filters
}


####################################################################
sub get_header{
    my @head = ();
    #first deal with VEP/SnpEff fields
    foreach my $csq (@csq_fields){
        push @head, $csq;
    }

    #GET HEADER INFO
    my @head_split= split("\t", $header[-1]);
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
                push @head, "Samples with variant", "ADs", "GQs";
            }
        }
        if (not defined $summarise_counts 
            or $summarise_counts eq 'all' 
            or $summarise_counts eq 'and'
         ){
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
            foreach my $g (@gene_anno_fields){
                push @head, $g;
            }
        }

    }
    return @head;
}



