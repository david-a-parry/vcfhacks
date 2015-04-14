
=head1 LICENSE
                                                                                                                     
 Copyright (c) David A. Parry
 University of Leeds
                                                                                                                     
=head1 CONTACT                                                                                                       

 David Parry <d.a.parry@leeds.ac.uk>
    
=cut

=head1 NAME

 SpliceConsensusFilter

=head1 SYNOPSIS

 mv SpliceConsensusFilter.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin SpliceConsensusFilter

=head1 DESCRIPTION

 A simple example VEP filter plugin that prints splice_region_variants only if they match certain criteria (ignores other variant classes). The plugin finds variants with a 'splice_region_variant' annotation and asseses whether the variants lie within a stricter definition of the splice consensus sequence (that is, 3 bp before the exon to the first 3 bp of the exon or the last bp of the exon to 6 bp after the exon) and outputs only these variants.

=cut

package SpliceConsensusFilter;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepFilterPlugin);

sub feature_types {
    return ['Transcript'];
}

sub include_line {
    my ($self, $tva) = @_;

    #we define splice consensus as...

    my $t = $tva->transcript;
    my $tv = $tva->transcript_variation;
    my $vf = $tva->variation_feature;
    #foreach my $c ($vf->consequence_type){
    foreach my $c (@{$tv->consequence_type}){
        if ($c eq 'splice_region_variant'){
            my @exons =  sort {$a->start() <=> $b->start()}  @{$t -> get_all_Exons()};
            my @introns =  sort {$a->start() <=> $b->start()}  @{$t -> get_all_Introns()};
            #check if it really is splice consensus
            if (my $intron_no = $tv->intron_number){#is intronic
                #find nearest exon and see if in splice consensus
                my ($i, $total) =  split(/\//, $intron_no); 
                if ($t->strand < 0){ #on - strand
                    $i = 1 + $total - $i;
                }
                my $down_dist = $vf->start - $exons[$i-1]->end();
                my $up_dist = $exons[$i]->start() -  $vf->end;
                if ($t->strand > 0){ #on + strand
                    return 1 if $up_dist <= 3;#in last 3 nt of intron
                    return 1 if $down_dist <= 6;#in first 6 nt of intron
                }else{
                    return 1 if $up_dist <= 6;
                    return 1 if $down_dist <= 3;
                }
            }
            if (my $exon_no = $tv->exon_number){
                #find nearest intron and see if in splice consensus
                my ($i, $total) =  split(/\//, $exon_no); 
                if ($t->strand < 0){ #on - strand
                    $i = 1 + $total - $i;
                }
                my $down_dist = $vf->start - $exons[$i-1]->start();
                my $up_dist = $exons[$i-1]->end() -  $vf->end;
                if ($t->strand > 0){ #on + strand
                    return 1 if $up_dist <= 0; #in last nt of exon
                    return 1 if $down_dist <= 2;#in first three nt of exon
 
                }else{
                    return 1 if $up_dist <= 2;
                    return 1 if $down_dist <= 0;
                }
            }
        }
    }
    return 0;#filter line
}

1;

