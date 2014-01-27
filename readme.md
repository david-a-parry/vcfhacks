__vcfhacks__

This project comprises a set of perl scripts and modules that may be useful for VCF manipulation when trying to discover disease-causing variants in sequencing data. Although more biologists are taking to the commandline many do not have the time or inclination to really get to grips with the more obtuse programs used in the field of bioinformatics. The aim of this project is to provide conceptually simple programs that perform common useful tasks. 

__UPDATE__

VERSION 0.7.4:

27/1/14

-added getFunctionalVariantsSnpEff.pl program and DbnsfpVcfFilter.pm module.
-added SnpEff annotation parsing functions to ParseVcf.pm.
-added features to getFunctionalVariantsVep.pl and findBiallelicVep.pl enabling filtering of splice_region_variants using annotations from my SpliceConsensus VEP plugin. In practice this means that you can select variants on a slightly stricter definition of the splice consensus region (3 bp before the exon to the first 3 bp of the exon or the last bp of the exon to 6 bp after the exon).
-following report of the Bio::SeqIO::entrezgene bug the bioperl team have fixed the issue, but at the time of writing the fix will not be in the current version of bioperl.  Instead you will need to replace the Bio::SeqIO::entrezgene.pm module in your bioperl installation with the one here: https://github.com/bioperl/bioperl-live/blob/master/Bio/SeqIO/entrezgene.pm if you wish to update your local database for ensemblGeneAnnotator.pl.
-various documentation tidy-ups
-various code tidy-ups, minor optimisations


VERSION 0.7.3:

16/12/13 

-added sortVcf.pl program
-fixed issues with splitMultiallelicVariants.pl always printing to STDOUT and missing out final header line. 
-Updated location of HMD_human5.rpt reference file for ensemblGeneAnnotator.pl.  Bio::SeqIO::entrezgene no longer parses all entrez gene asn.1 records correctly for creating the database - bug reported but an acceptable workaround needed so it is recommended to only use existing databases at this time.

__INSTALLATION__

Unzip the downloaded file and ensure you keep the .pl scripts in the same directory as the .pm module files.  The exceptions to this are the two Variant Effect Predictor (VEP) plugin modules (SpliceConsensus.pm and SpliceConsensusFilter.pm) which should be installed in your VEP cache 'Plugins' folder if you want to use them. The SpliceConsensus.pm can be used for variant filtering purposes in getFunctionalVariantsVep.pl and findBiallelicVep.pl.

You will also need to install the Tabix.pm perl module by Heng Li. Download the latest version of tabix from github (https://github.com/samtools/tabix) or clone the git repository (git clone https://github.com/samtools/tabix.git). If you downloaded the zipped version unzip to create the tabix directory. Next, cd into the new tabix directory and run 'make', cd into the 'perl' subdirectory, run 'perl Makefile.PL' and 'make test'. If tests succeed run '[sudo] make install' to complete Tabix.pm installation. If you get the error 'Subroutine Tabix::tabix_open redefined...' this is harmless and can be removed by replacing: 

<pre><code>require XSLoader;
XSLoader::load('Tabix', $VERSION);
</code></pre>

with 

<pre><code>{no warnings 'redefine';
require XSLoader;
XSLoader::load('Tabix', $VERSION);
}
</code></pre>

Other perl modules required by these scripts are installable via CPAN - perl will complain that they are not available in "@INC" when you attempt to run these programs if they are not on your system. Please see http://www.cpan.org/modules/INSTALL.html for instructions on how to install these modules.

__CODE__

These scripts were written __David A. Parry__, a geneticist/molecular biologist at the University of Leeds.  These were written and modified over a period of self-tutelage in which my programming styles/skills have been changeable so anyone browsing code may see varying states of readability and documentation. This is an issue I hope to tidy up over time. 

__CREDIT__

These scripts were written __David A. Parry__, a geneticist/molecular biologist at the University of Leeds.

__If you find these programs useful and use them for a paper__ ***please cite the URL <https://sourceforge.net/projects/vcfhacks/> *** (although, especially if you are local, an authorship is always much appreciated).

__COPYRIGHT AND LICENSE__

Copyright 2013  David A. Parry

These programs are free software: you can redistribute them and/or modify them under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. These programs are distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.



__Thanks!__
