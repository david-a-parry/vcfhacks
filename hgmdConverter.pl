#!/usr/bin/env perl 
use strict;
use warnings;

my $in = shift; 
my $FH;
if ($in =~ /\.gz$/){
    open ($FH, "gzip -dc $in |") or die "Can't open $in via gzip: $!\n";
}else{
    open ($FH, $in ) or die "Can't open $in for reading: $!\n";
}
my %hgmd_columns = ();
my $header = <$FH>;
chomp (my @head = split("\t", $header)); 
my $n = 0;
%hgmd_columns = map { lc($_) => $n++ } map { (my $trim = $_) =~ s/^#+//; $trim } @head;
foreach my $req
    (
    "hgmd id",
    "disease",
    "variant class",
    "gene symbol",
    "chromosome",
    "coordinate start",
    "coordinate end",
    "strand",
    "hgvs",
){
    if (not exists $hgmd_columns{$req}){
        die "Required column ($req) not found in HGMD file $in\nFound the following columns:\n" 
            .join("\n", @head) . "\n";
    }
}
while (my $line = <$FH>){
    chomp($line);
    my @out = ();
    my @details = ();
    my @fields = split("\t", $line); 
    my $strand = $fields[ $hgmd_columns{strand} ] ;
    my ($ref, $alt) = getHgmdAlleles(\@fields); 
    next if (not $ref and not $alt);
    foreach my $f ( "chromosome", "coordinate start" ){
        push @out, $fields[ $hgmd_columns{$f} ];
    }
    if (length($alt) > length($ref)){#insertion
        #VEP represents insertions as having an end -1 of start
        push @out, $fields[ $hgmd_columns{"coordinate start"} ] -1 ;
    }else{
        push @out, $fields[ $hgmd_columns{"coordinate end"} ] ;
    }
    $ref = $ref ? $ref : "-";
    $alt = $alt ? $alt : "-";
    push @out, "$ref/$alt";
    push @out, $strand;
    foreach my $d ( "hgmd id", "disease", "variant class", "gene symbol", "hgvs"){
        (my $no_space = $fields[ $hgmd_columns{$d} ]) =~ s/\s/_/g;
        push @details, $no_space;
    }
    
    push @out, join("|", @details);
    print join("\t", @out) . "\n";
}

###########################################################
sub getHgmdAlleles{
    my $hgmd = shift;#array ref to split line from HGMD file
    my $hgvs = $hgmd->[ $hgmd_columns{hgvs} ] ;
    my $ref = "";
    my $alt = "";
    if ($hgvs =~ /c\.\d+([\+-]\d+)*([ATGCN])>([ATGCN])$/i){ #SNV
        $ref = $2;
        $alt = $3;
    }elsif ($hgvs =~ /c\.\d+([\+-]\d+)*(_\d+([\+-]\d+)*)*del([ATGCN]+)$/i){
        $ref = $4;
    }elsif ($hgvs =~ /c\.\d+([\+-]\d+)*(_\d+([\+-]\d+)*)*(ins|dup)([ATGCN]+)$/i){
        $alt = $5;
    }elsif ($hgvs =~ /c\.\d+([\+-]\d+)*(_\d+([\+-]\d+)*)*del([ATGCN]+)ins([ATGCN]+)$/i){
        $ref = $4;
        $alt = $5;
    }
    return ($ref, $alt);
}
