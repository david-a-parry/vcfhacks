#vcfhacks examples#

###Getting started###

This is the examples document for the binary versions of vcfhacks. If you are using the perl scripts (i.e. files with a '.pl' extension) see the original examples documentation [here](https://github.com/gantzgraf/vcfhacks/blob/master/examples.md). 

These programs are all command line utilities. To run these programs you simply need to download and extract the tar.bz2 file for your plaform (currently only 64 bit Linux Mac OS X systems are supported) and change into the newly created vcfhacks_binaries directory. Each program can be run from this directory using the command *./[program_name]*  (e.g. *./annotateSnps*) or you may prefer to move the programs somewhere in your $PATH or add the new vcfhacks_binaries to your $PATH to be able to run these programs from any directory. 

These programs are written primarily with a view to assisting mutation discovery in Mendelian disease. They have been tested on data generated using GATK variant callers (UnifiedGenotyper and HaplotypeCaller) but should work with data generate in VCF format from other variant callers. When using tools such as findBiallelicVep or filterOnSample it is assumed that you are using a multisample VCF as input.  Help information is available for most programs by running the program with either '--help' or '--manual' options. 

The following programs support parallel execution using forks:

    annotateSnps
    filterOnEvsMaf
    filterVcfOnVcf
    filterOnSample
    sampleCallsToInfo 

Use the *--forks*  option with these programs to specify the number of forks to use (if desired). It is probably a bad idea to specify more forks than there are CPU cores on your machine and these programs will warn you if this happens.  

Examples are given below of what are envisaged to be typical uses of these programs. They are only intended as guidelines for use of some of these programs in order to get people started with these tools. 

###FILTERING VARIANTS:###

Typically we may first want to remove common variation from a VCF. The example below is typical for processing data relating to a rare recessive condiation. It removes variants in the **dbSnp VCF file** (dbSnp138.b37.vcf.gz) present at **build 129** (generally considered the last build before 'contamination' with data from large sequencing projects) or with a **minor allele frequency (MAF) of 1 % or higher** in later build versions.  However, any variant with a **'Pathogenic'** or **'Probably Pathogenic'** annotation in the **ClinVar** VCF file (clinvar_20130506.vcf) will be retained regardless of frequency or dbSNP build. We optionally use the --progress flag to show a progressbar as variants are processed.

    ./annotateSnps -d dbSnp138.b37.vcf.gz clinvar_20130506.vcf -b 129 -f 1 --pathogenic -i input.vcf -o input_snpfiltered.vcf --progress

However, these days, we probably *do not need to use the fudged tactic of removing SNPs from dbSNP129 or earlier* thanks to the wealth of frequency data no available both in dbSNP and other databases, so generally I would omit the --build option such as in the example below:

    ./annotateSnps -d dbSnp138.b37.vcf.gz clinvar_20130506.vcf -f 1 --pathogenic -i input.vcf -o input_snpfiltered.vcf --progress

Similarly, we may want to use data from the **NHLBI Exome Sequencing Project (ESP)** to filter common variants. Below we filter any variant with a MAF of 1 % or higher in the ESP VCF files present in the directory passed to the -d argument:

    ./filterOnEvsMaf -i input_snpfiltered.vcf -o input_snpfiltered_evsfiltered.vcf -f 1 -d ~/ESP_vcfs --progress

Having filtered common variation from dbSNP and ESP databases we may want to **filter variants present in other VCFs** (e.g. control samples sequenced locally).

    ./filterVcfOnVcf -i input_snpfiltered_evsfiltered.vcf -f controls.vcf -o input_snpfiltered_evsfiltered_vcffiltered.vcf

We may only want to filter using these VCFs if variants are present above a certain **allele frequency**, in which case we can pass the **--allele_frequency_filter/-y** option. The following command filters any variant present at an allele frequency of 1 % or higher in the controls.vcf file: 

    ./filterVcfOnVcf -i input_snpfiltered_evsfiltered.vcf -f controls.vcf -o input_snpfiltered_evsfiltered_vcffiltered.vcf -y 0.01 

We can also filter using VCFs such as those available from the **ExAC** project using allele frequency (AF) annotations in the VCF INFO field if we pass the **--info_filter/-w** option: 

    ./filterVcfOnVcf -i input_snpfiltered_evsfiltered.vcf -f ExAC.r0.1.sites.vep.vcf.gz -o input_snpfiltered_evsfiltered_vcffiltered.vcf -y 0.01  -w 

[Also see notes on *sampleCallsToInfo* and filterVcfOnVcf under the *MISC USEFUL TOOLS* section.]

If we have a multisample VCF we may want to only keep variants present in our sample of interest or family members:

    ./filterOnSample -i input_snpfiltered_evsfiltered_vcffiltered.vcf -s Sample1 Sample2 -x Sample3 Sample4 -o input_snpfiltered_evsfiltered_vcffiltered_samplefiltered.vcf

In the example above Sample1 and Sample2 might be affected individuals while Sample3 and Sample4 could be unaffected family members who are potential carriers. The **filterOnSample** program will only output variants present in Sample1 and Sample2 but will ignore anything that is also present in any other sample in the VCF with the exception of Sample3 and Sample4. Alternatively you can explicitly specify which samples to use for filtering, such as the example below where variants are only kept if present in Sample1 and Sample2 but rejected if Control1 or Control2 contain the same variant:

    ./filterOnSample -i input_snpfiltered_evsfiltered_vcffiltered.vcf -s Sample1 Sample2 -r Control1 Control2 -o input_snpfiltered_evsfiltered_vcffiltered_samplefiltered.vcf

The above might be useful for identifying mutations with dominant inheritance patterns by identifying mutations only present in affected samples. See the full documentation for filterOnSample 

At this point, you may want to identify variants that have potentially pathogenic functional consequences on gene products. In order to do so you will need to have annotated your VCF with functional consequences using either Ensembl's *Variant Effect Predictor* (http://www.ensembl.org/info/docs/tools/vep/index.html) or *SnpEff* (http://snpeff.sourceforge.net/).  You can filter this VCF using getFunctionalVariantsVep/getFunctionalVariantsSnpEff or **findBiallelicVep/findBiallelicSnpEff** depending on your needs.  Below is an example using findBiallelicVep to identify potential compound heterozygous or homozygous 'functional' variants present in Sample1 and Sample2 (they are related so we are  using the -e flag to specify that we expect both to contain the same causative mutation) but not in Sample3 and Sample4 (unaffected family members):

    ./findBiallelicVep -i vep.input_snpfiltered_evsfiltered_vcffiltered.vcf -s Sample1 Sample2 -r Sample3 Sample4 -e -o output.vcf 

The above might be better acheived by using a *PED file* of the family instead:

    ./findBiallelicVep -i vep.input_snpfiltered_evsfiltered_vcffiltered.vcf -f family.ped -o output.vcf 

You might instead be looking at several unrelated samples with the same condition and want to identify shared genes with potential compound heterozygous or homozygous 'functional' variants (not necessarily the same causative variants in each sample). In the example below we look at four unrelated samples with the same condition and output variants for genes that have potential biallelic variation in 3 or more of the 4 samples:

    ./findBiallelicVep -i vep.input_snpfiltered_evsfiltered_vcffiltered.vcf -s Sample1 Sample2 Sample3 Sample4 -n 3 -o output.vcf -l shared_genes.txt

See the various options in the help documentation for findBiallelicVep for more details on how to filter variants on VEP fields such as 1000 genomes allele frequency or PolyPhen predictions.

The **getFunctionalVariantsVep** program identifies variants matching 'functional' criteria and can also be used to provide a list of gene IDs to ignore (e.g. lof tolerant genes).  You may want to use this program to filter genes on functional consequence when not looking for biallelic variation or to prefilter with a list of gene IDs before running findBiallelicVep. You may also/instead want to use rankOnCaddScore to rank variants ignoring functional annotation using Combined Annotation Dependent Depletion (CADD) scores.

    ./rankOnCaddScore -i input.vcf -o cadd_ranked.vcf -c ~/cadd_score_files/*.gz -n not_found_cadd_scores.tsv

You may then want to submit your 'not_found_cadd_scores.tsv' to http://cadd.gs.washington.edu/score and rerun once you have downloaded and tabix indexed your missing results.

###ANNOTATING YOUR OUTPUT:###

Having filtered/ranked you variant calls you may want to annotate gene information using ensemblGeneAnnotator . You may either create your own database with ensemblGeneAnnotator or download a pre-build database from: 

    http://sourceforge.net/projects/vcfhacks/files/ensAnnotatorDatabase/

Remember to unzip the database if you have downloaded the pre-built version. Annotation requires ensembl gene IDs annotated using Ensembl's variant effect predictor. 

    ./ensemblGeneAnnotator -i input.vcf -d ~/ensAnnotator -o annotated.vcf

You may only want to annotate gene information for 'functional' variants and omit information for overlapping genes with synonymous/intronic/UTR variants, in which case use a command like the the following:

    ./ensemblGeneAnnotator -i input.vcf -d ~/ensAnnotator -o annotated.vcf --functional

At this point you may wish to produce a spreadsheet of your results using annovcfToSimple. The following should create an Excel spreadsheet (.xlsx) from a VEP and ensemblGeneAnnotator annotated file:

    ./annovcfToSimple -g -v -i input.vcf -o input.xlsx

If you have many samples in your VCF you may want to specify only a few using the --samples option or alternatively use the --summarise option to summarise the alleles and genotypes in your VCF rather than giving information for each sample.


###GETTING VARIANTS FROM COORDINATES:###

You may quickly obtain variants that lie within genomic regions (specified in a BED file or on the command line), overlap single genomic coordinates or match those of another VCF using getVariantsByLocation. If an index does not already exist one will be created, which may take some time for large files, but once created retrieval of regions should be quick. Note that your input must be coordinate sorted for this to work, otherwise the program will error.

To get variants from vars.vcf that lie within regions in bedfile regions.bed: 
 
    ./getVariantsByLocation -i vars.vcf -b regions.bed 
 
To get variants from vars.vcf that lie within the region 1:2000000-50000000 and output to file filtered.vcf
 
    ./getVariantsByLocation -i vars.vcf -r 1:2000000-50000000 -o filtered.vcf

To get variants from vars.vcf that lie within regions in bedfile regions.bed or the region 1:2000000-50000000 and output to file filtered.vcf
 
    ./getVariantsByLocation -i vars.vcf -b regions.bed -r 1:2000000-50000000 -o filtered.vcf

To get variants from vars.vcf that overlap any variant in other.vcf 
 
    ./getVariantsByLocation -i vars.vcf -v other.vcf
    
To get variants from vars.vcf that match a variant in other.vcf
 
    ./getVariantsByLocation -i vars.vcf -v other.vcf --matching

To get variants from vars.vcf that overlap the coordinates 1:2000000 or 1:30000000: 
    
    ./getVariantsByLocation -i vars.vcf -r 1:2000000 1:30000000 -o filtered.vcf

If you have a VCF that is no longer coordinate sorted (e.g. if you've ranked your variants on CADD scores), you may use filterVcfOnLocation to output variants that lie within specific genomic regions. This may be significantly slower compared to using getVariantsByLocation on an indexed VCF but is useful if you do not want to coordinate sort your file.

To get variants from vars.vcf that lie within regions in bedfile regions.bed: 
    
    ./filterVcfOnLocation -i vars.vcf -b regions.bed

To get variants from vars.vcf that lie within the region 1:2000000-50000000 and output to file filtered.vcf

    ./filterVcfOnLocation -i vars.vcf -r 1:2000000-50000000 -o filtered.vcf


###MISC USEFUL TOOLS###

If you need to sort a VCF in coordinate order use sortVcf:

    ./sortVcf -i input.vcf -o output.vcf

To get the names of all samples and the total number of samples in a VCF use getSampleNames

    ./getSampleNames -i vars.vcf

To get the total number of variants in a VCF:

    ./countVariants vars.vcf

To get a summary of the number of observed alleles and genotypes in a VCF use countVcfCalls:
    
    ./countVcfCalls -i vars.vcf

If you want to see the sample names with ALT genotypes as well as the number of observed genotypes use the -r flag:
    
    ./countVcfCalls -i vars.vcf -r 

You can pipe to countVcfCalls from STDIN by specifying '-' as the input file. So, if you have a particular variant you want compare counts for at position 1:2000000 you may use the following command:

    ./getVariantsByLocation -i vars.vcf -r 1:2000000 | perl countVcfCalls -i -

If you have a VCF with many hundreds to thousands of samples that you frequently use for filtering on allele frequency using filterVcfOnVcf, you can summarise the allele counts with sampleCallsToInfo:

    ./sampleCallsToInfo -i vars.vcf -o summarised.vcf 

For data from a VCF with many hundreds to thousands of samples, using a summarised VCF with filterVcfOnVcf is likely to be a great deal quicker than using the original VCF. To filter on allele frequency or genotypes use filterVcfOnVcf  with the --info_filter flag

    ./filterVcfOnVcf -i input.vcf  -f summarised.vcf  -o filtered.vcf -y 0.01 --info_filter



