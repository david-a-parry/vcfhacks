# vcfhacks

This project comprises a set of perl programs and modules for VCF manipulation when trying to discover disease-causing variants in sequencing data. This project provides conceptually simple programs that can perform variant filtering and annotation to aid identification of disease-causing variants.

Usage examples are included in the examples.md markdown documents which can also be viewed at https://github.com/gantzgraf/vcfhacks/blob/master/examples.md or https://github.com/gantzgraf/vcfhacks/blob/master/examples_bin.md

## INSTALLATION

You may either download these programs as perl scripts or as precompiled binaries for 64-bit Linux or Mac OS X. The binary executables are provided for ease so that users do not have to install the various perl modules required by these scripts. Users that are comfortable installing perl modules may prefer to use the scripts, which take up less space, do not suffer the same lag in startup and provide a scrollable manual page in most cases. The scripts also allow you to receive more regular updates using git (see below).

### BINARIES

These programs are all command line utilities. To run these programs you simply need to download and extract the tar.bz2 file for your plaform from https://github.com/gantzgraf/vcfhacks/releases (currently only 64 bit Linux Mac OS X systems are supported) and change into the newly created vcfhacks_binaries directory. Each program can be run from this directory using the command *./[program_name]*  (e.g. *./annotateSnps*) or you may prefer to move the programs somewhere in your $PATH or add the new vcfhacks_binaries to your $PATH to be able to run these programs from any directory. You may need to make the programs executable before they will run (e.g. by running 'chmod +x *' from within your vcfhacks_binaries directory). 

To check whether the binaries are functioning and compatible with your system type 'prove' from within the downloaded and extracted directory to run a small series of tests.

### SCRIPTS

You may either download and unpack a script bundle from a specific release (https://github.com/gantzgraf/vcfhacks/releases), or to stay up to date, you may clone this repository: 

    git clone --recursive https://github.com/gantzgraf/vcfhacks.git

...and regularly run the following commands to receive updates. 
    
    git pull --recurse-submodules 
    git submodule update --recursive --remote 

The .pl scripts must remain in the same directory as the 'lib' and 'data' folders but you may create symlinks to the scripts (e.g. in your ~/bin folder) if you wish.  The vep_plugins folder also contains two Variant Effect Predictor (VEP) plugin modules (SpliceConsensus.pm and SpliceConsensusFilter.pm) which should be installed in your VEP cache 'Plugins' folder if you want to use them. The SpliceConsensus.pm annotations can be used for variant filtering purposes in getFunctionalVariants.pl and findBiallelic.pl.

To run the scripts either add the enclosing directory to your PATH and make sure the scripts are executable (e.g. 'chmod +x vcfhacks/*.pl'), run the scripts directly from the folder (e.g. ./annotateSnps.pl) or run each script using 'perl [location of script]' followed by the required arguments. Help information is available for most scripts by running the script with either '--help' or '--manual' options. 

If you want to use scripts to search bgzip compressed VCFs or use rankOnCaddScore.pl, as of version 0.2.2 you will also need to build and install the Bio::DB::HTS::Tabix. The easiest way to acheive this is to download the tarball from:

    http://search.cpan.org/dist/Bio-DB-HTS/lib/Bio/DB/HTS/Tabix.pm
    
Unpack the tarball cd into the newly created directory and run the 'INSTALL.pl' script. For example:

    wget http://search.cpan.org/CPAN/authors/id/R/RI/RISHIDEV/Bio-DB-HTS-2.7.tar.gz
    tar xvf Bio-DB-HTS-2.7.tar.gz
    cd Bio-DB-HTS-2.7
    perl INSTALL.pl

Alternatively, see below for an attempt to automatically install the required perl modules.

Other perl modules required by these scripts are installable via CPAN - perl will complain that they are not available in "@INC" when you attempt to run these programs if they are not on your system. Please see http://www.cpan.org/modules/INSTALL.html for instructions on how to install these modules. Below is a list of these non-core modules that you are likely to need install:

    List::MoreUtils
    Parallel::ForkManager
    Term::ProgressBar
    Bio::DB::HTS::Tabix
    List::MoreUtils
    LWP::Simple (geneAnnotator.pl only)
    LWP::Protocol::https (for geneAnnotator.pl download of some database files only)
    HTTP::Tiny (geneAnnotator.pl for remote retrieval of gene IDs only)
    JSON (geneAnnotator.pl for remote retrieval of gene IDs only)
    Excel::Writer::XLSX (annovcfToSimple.pl only)
    Bio::DB::Sam (hgmdMartToVcf.pl only)
    Statistics::R (geneBurdenTest.pl only)

An experimental script is also provided to attempt to automate the installation of the required modules. This requires developer tools to be installed (e.g. xcode on Mac OS X or base-devel on debian/ubuntu systems) and for your CPAN installation to be functioning. To try automated install of these modules try the following (prefix with sudo if you have sudo priveleges and require system-wide installation): 
    
    perl install_modules/installModules.pl

#### Testing 

You can test your installation (scripts or binaries) by simply running the following command from within the vcfhacks directory:
   
    prove
    
The 'prove' command will run all tests in the 't' folder. The tests are not exhaustive, but should give you a good idea about whether there is an issue with your installation.

In order for all tests to complete successfully, you will need the geneAnnotator database, which can be downloaded from: https://github.com/gantzgraf/vcfhacks/releases/download/v.0.2.1/geneAnnotatorDb.zip. Unzip this file and place the 'geneAnnotatorDb' folder in vcfhacks 'data' directory'.

__CREDIT__

These programs were written __David A. Parry__, currently at the University of Edinburgh. 

__If you find these programs useful and use them for a paper__ ***please cite the URL <https://github.com/gantzgraf/vcfhacks> *** 

__COPYRIGHT AND LICENSE__

Copyright 2013,2014,2015  David A. Parry

These programs are free software: you can redistribute them and/or modify them under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. These programs are distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.



