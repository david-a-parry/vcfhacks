use strict;
use warnings;
use Test::More tests => 9;
use FindBin qw($RealBin);
use lib "$RealBin/../lib";
use lib "$RealBin/../lib/Bioperl";
use lib "$RealBin/../lib/BioASN1EntrezGene/lib";

#test loading of 9 local modules
require_ok("VcfReader");
require_ok("ClinVarReader");
require_ok("ParsePedfile");
require_ok("VcfhacksUtils");
require_ok("TextToExcel");
require_ok("Tabix");
require_ok("SortGenomicCoordinates");
require_ok("IdParser");
require_ok("Bio::SeqIO::entrezgene");

