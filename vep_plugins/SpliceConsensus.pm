=head1 LICENSE
                                                                                                                     
 Copyright (c) David A. Parry
 University of Leeds
                                                                                                                     
=head1 CONTACT                                                                                                       

 David Parry <d.a.parry@leeds.ac.uk>
    
=cut

=head1 NAME

 SpliceConsensus

=head1 SYNOPSIS

 mv SpliceConsensus.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin SpliceConsensus

=head1 DESCRIPTION

 A VEP plugin that adds a splice_consensus annotation field. Finds variants with a 'splice_region_variant' annotation and asseses whether the variants lie within a stricter definition of the splice consensus sequence (that is, 3 bp before the exon to the first 3 bp of the exon or the last bp of the exon to 6 bp after the exon). Could do with better indel handling.
=cut

package SpliceConsensus;

    use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);
    
    sub feature_types {
        return ['Transcript'];
    }

    sub get_header_info {
        return {
            SPLICE_CONSENSUS => "Splice consensus as defined by 3 bp before, 6 bp after or first 3 bp of exon or last bp of exon"
        };
    }

sub run {
    my ($self, $tva) = @_;

    my $t = $tva->transcript;
    my $tv = $tva->transcript_variation;
    my $vf = $tva->variation_feature;
    foreach my $c (@{$tv->consequence_type}){
        if ($c eq 'splice_region_variant'){
            my @exons =  sort {$a->start() <=> $b->start()}  @{$t -> get_all_Exons()};
            my @introns =  sort {$a->start() <=> $b->start()}  @{$t -> get_all_Introns()};
            #check if it really is splice consensus
            if (my $intron_no = $tv->intron_number){#is intronic
                #find nearest exons and see if in splice consensus
                my ($i, $total) =  split(/\//, $intron_no); 
                if ($t->strand < 0){ #on - strand
                    $i = 1 + $total - $i;
                }
                my $down_dist = $vf->start - $exons[$i-1]->end();
                my $up_dist = $exons[$i]->start() -  $vf->end;
                if ($t->strand > 0){ #on + strand
                    if ($up_dist <= 3){
                    #in last 3 nt of intron
                        return {SPLICE_CONSENSUS => 
                            "SPLICE_CONSENSUS-Intronic-$up_dist"};
                    }
                    if ($down_dist <= 6){
                    #in first 6 nt of intron
                        return {SPLICE_CONSENSUS => 
                            "SPLICE_CONSENSUS-Intronic+$down_dist"};
                    }
                }else{
                    if ($up_dist <= 6){
                    #in first 6 nt of intron
                        return {SPLICE_CONSENSUS => 
                            "SPLICE_CONSENSUS-Intronic+$up_dist"};
                    }
                    if($down_dist <= 3){
                    #in last 3 nt of intron
                        return {SPLICE_CONSENSUS => 
                            "SPLICE_CONSENSUS-Intronic-$down_dist"};
                    }
                }
            }
            if (my $exon_no = $tv->exon_number){
                #find nearest exon and see if in splice consensus
                my ($i, $total) =  split(/\//, $exon_no); 
                if ($t->strand < 0){ #on - strand
                    $i = 1 + $total - $i;
                }
                my $down_dist = $vf->start - $exons[$i-1]->start();
                my $up_dist = $exons[$i-1]->end() -  $vf->end;
                if ($t->strand > 0){ #on + strand
                    if ($up_dist <= 0){ #in last nt of exon
                        $up_dist++;
                        return {SPLICE_CONSENSUS => 
                            "SPLICE_CONSENSUS-Exonic-$up_dist"};
                    }
                    if ($down_dist <= 2){#in first three nt of exon
                        $down_dist++;
                        return {SPLICE_CONSENSUS => 
                            "SPLICE_CONSENSUS-Exonic-$down_dist"};
                    }
                }else{
                    if ($up_dist <= 2){
                        $up_dist++;
                        return {SPLICE_CONSENSUS => 
                            "SPLICE_CONSENSUS-Exonic-$up_dist"};
                    }
                    if ($down_dist <= 0){
                        $down_dist++;
                        return {SPLICE_CONSENSUS => 
                            "SPLICE_CONSENSUS-Exonic-$down_dist"};
                    }
                }
            }
        }
    }
    return {};
}

    1;
