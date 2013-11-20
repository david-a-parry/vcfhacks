__vcfhacks__

This project comprises a set of perl scripts and modules that may be useful for VCF manipulation when trying to discover disease-causing variants in sequencing data. Although more biologists are taking to the commandline many do not have the time or inclination to really get to grips with the more obtuse programs used in the field of bioinformatics. The aim of this project is to provide conceptually simple programs that perform common useful tasks. 

__INSTALLATION__

Unzip the downloaded file and ensure you keep the .pl scripts in the same directory as the .pm module files.  The exceptions to this are the two Variant Effect Predictor (VEP) modules (SpliceConsensus.pm and SpliceConsensusFilter.pm) which should be installed in your VEP cache 'Plugins' folder (assuming you have the VEP installed). The SpliceConsensus.pm annotations are not yet integrated into the other scripts but shall be shortly.

You will also need to install the Tabix.pm perl module by Heng Li. Download the latest version of tabix from https://sourceforge.net/projects/samtools/files/tabix/, unpack the tabix-x.x.x.tar.bz2 file and run 'make' in the newly created 'tabix-x.x.x' directory. Next cd into the 'perl' subdirectory, run 'perl Makefile.PL' followed by 'make test'. If tests succeed run 'sudo make install' to complete Tabix.pm installation.

__CODE__

These scripts were written __David A. Parry__, a geneticist/molecular biologist at the University of Leeds.  These were written and modified over a period of self-tutelage in which my programming styles/skills have been changeable so anyone browsing code may see varying states of readability and documentation. This is an issue I hope to tidy up over time. 

__CREDIT__

These scripts were written __David A. Parry__, a geneticist/molecular biologist at the University of Leeds.

__If you find these programs useful and use them for a paper__ ***please cite the URL <https://sourceforge.net/projects/vcfhacks/> *** (although, especially if you are local, an authorship is always much appreciated).

__COPYRIGHT AND LICENSE__

Copyright 2013  David A. Parry

These programs are free software: you can redistribute them and/or modify them under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. These programs are distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.



__Thanks!__
