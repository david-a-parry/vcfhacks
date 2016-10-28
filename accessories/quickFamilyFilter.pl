#!/usr/bin/env perl
#David Parry January 2014
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Term::ProgressBar;
use File::Temp qw/ tempfile / ;
use File::Basename;
use Data::Dumper;
use FindBin qw ( $RealBin ) ;
use lib "$FindBin::Bin/../lib";
use ParsePedfile;

my $script_prefix = "perl $RealBin/..";
my @input = ();
my @peds = ();
my @args = ();
my $eddie = 0;

my $cadd_dir = "$ENV{HOME}/ref/GRCh37/cadd/v1.3/";
my $gatk = "$ENV{HOME}/programs/GATK/v3.6/GenomeAnalysisTK.jar";
my $fasta = "$ENV{HOME}/ref/GRCh37/hs37d5.fa";

if ($ENV{HOSTNAME} and $ENV{HOSTNAME} =~ /^login\d+.ecdf.ed.ac.uk$/){
    $cadd_dir = "/exports/igmm/eddie/aitman-lab/ref/GRCh37/cadd/v1.3/";
    $gatk = "$ENV{HOME}/GATK/v3.6/GenomeAnalysisTK.jar";
    $fasta = "/exports/igmm/eddie/aitman-lab/ref/GRCh37/hs37d5.fa";
    $eddie = 1;
}

my %opts = 
(
    i        => \@input, 
    p        => \@peds, 
    a        => \@args,
    fasta    => $fasta,
    cadd_dir => $cadd_dir,
    gatk_jar => $gatk,
);

sub usage{
    my $msg = shift;
    print STDERR "ERROR: $msg\n" if $msg;
    print <<EOT

    Usage: $0 -i input.vcf -p pedfile.ped [options] 
           $0 -i input_1pc_filtered.vcf input_0pt1pc_filtered.vcf -p pedfile.ped [options]

    Options:
    
        -i,--input
            Input VCF file(s). You may specify two VCF files, where the first is to be provided for recessive analysis while the second is provided for dominant analysis (e.g. you might have filtered on a lower allele frequency for dominant/de novo analysis).

        -p,--ped
            One or more PED files to process

        -d,--do_not_exclude
            Use this flag to disable the default excluding of variants in all samples not in a given family

        -u,--unaffecteds_file
            A file containing one sample ID per line (or separated by whitespace) to use for excluding variants in all analyses

        -x,--exclude_except_file
            A file containing one sample ID per line (or separated by whitespace) to not be used for excluding variants in all analyses, but instead all samples not in a given family and not in this file will be used for filtering.

        -y,--sample_frequency
            Value between 0 and 1 to use for frequency filtering in filterOnSample.pl commands. Default is to filter if present in any non-familial unaffected sample. 

        -c,--cadd_score
            Cadd score input before filtering

        --cadd_dir
            Directory containing CADD score files - default = $cadd_dir

        -o,--denovo_only
            Only look for de novo variants

        -g,--gatk_denovos
            Use GATK to select hiConfDeNovos

        --gatk_jar
            Location of GATK jar file - default = $gatk

        --fasta
            Location of fasta for GATK commands - default = $fasta

        -s,--seg_only
            Only look for dominant or recessive variants - do not look for de novo variants

        -l,--lentient
            Also perform 'lenient' findBiallelic.pl commands (i.e. allowing for de novo occurence of alleles or hemizygosity masking a change).

        -z,--denovo_mode
            Also perform '--denovo_biallelic_mode' findBiallelic.pl commands (looking only for biallelic variants where an allele has occured de novo).

        -b,--progress
            Use a progress bar for each command.

        -f,--forks
            Use this many forks for the scripts that support it.

        -n,--no_exec
            Print but do not execute commands.

        -q,--qsub
            Submit jobs through qsub scripts. If a name of a directory is given here, the scripts will be output into this directory. 

        -h,-?,--help
            Show this message and exit

EOT
;
    exit 1 if $msg;
    exit;
}

GetOptions(
        \%opts,
        'i|input=s',
        'p|ped=s{,}',
        'd|do_not_exclude',
        'z|denovo_mode',
        'u|unaffecteds_file=s',
        'x|exclude_except_file=s',
        'a|args=s{,}',
        'y|sample_frequency=f',
        'l|lenient',
        'b|progress',
        'f|forks=i',#forks for filterOnSample.pl commands
        'n|no_exec',
        'c|cadd_score:f',
        's|seg_only',
        'o|denovo_only',
        'q|qsub:s',
        'g|gatk_denovos',
        'gatk_jar=s',
        'fasta=s',
        'cadd_dir=s',
        'h|?|help',
) or die "Syntax error in option spec!\n";
usage() if $opts{h};
usage("-i,--input argument is required\n") if not @input;
usage("Please supply only one or two input files, not " . scalar(@input) . "\n")
  if @input > 2;
usage("-p,--ped argument is required\n") if not $opts{p};
my $extra_args = join(" ", @args);

my $script_dir = '.';
if ($opts{q}){
    $script_dir = $opts{q};
    if (not -d $script_dir){
        mkdir($script_dir) or die "Could not create subscript directory ".
                                  "'$script_dir': $!\n"
    }
}

my %commands = (); 
for (my $i = 0; $i < @input; $i++){
    if (defined $opts{c}){
        ( my $out_stub = $input[$i] ) =~ s/\.vcf(\.(b)*gz)*$//; 
        my $scored = "$out_stub.caddscored.vcf";
        my $cadd_cmd = "$script_prefix/rankOnCaddScore.pl -r $cadd_dir -i $input[$i] -o $scored -n not_found.$scored -d ";
        $cadd_cmd .= " --progress " if $opts{b};
        push @{$commands{cadd}}, $cadd_cmd;
        $input[$i] = $scored;
    }
}
    
foreach my $ped (@{$opts{p}}){
    my $p = ParsePedfile->new( file => $ped );
    foreach my $f ( sort  $p->getAllFamilies ) {
        my $pf = fileparse($ped);
        push @{$commands{"$pf-$f"}}, getSegCommands($f, $p);
    }
}

foreach my $k (keys %commands){
    @{$commands{$k}} = map { /^java/ ? $_ : "perl $RealBin/$_" } @{$commands{$k}}; 
}

if ($opts{n}){
    printCommands();
}elsif (defined $opts{q}){
    qsubCommands();
}else{
    executeCommands();
}

#####################################################################
sub execute{
    my $cmd = shift;
    print STDERR "[Executing] $cmd\n";
    system($cmd); 
    die "Error $? / $@\n" if $?;
}

#####################################################################
sub qsubCommands{
    my @hold_ids = ();
    if (exists $commands{cadd}){
        foreach my $cmd (@{$commands{cadd}}){
            push @hold_ids, qsub($cmd);
        }
    }
    foreach my $k (keys %commands){
        next if $k eq 'cadd';
        my $cmds = join("\n\n", (@{$commands{$k}}));
        qsub($cmds, $k, \@hold_ids);
    }
}

#####################################################################
sub qsub{
    my $cmd = shift;
    my $id  = shift;
    my $hold = shift;
    my $hold_string = '';
    if ($hold){
        if (@$hold){
            $hold_string = " -hold_jid " . join(",", @$hold);
        }
    }
    my $script = makeScript($cmd, $id); 
    my $sub = "qsub$hold_string $script";
    print STDERR "Executing: $sub\n";
    my $output = `$sub`;
    die "Error $? / $@\n" if $?;
    if ($output =~ /Your job (\d+) .* has been submitted/){
        return $1;
    }else{
        die "Error parsing qsub output for $cmd!\nOutput was: $output";
    }
}

#####################################################################
sub makeScript{
    my $cmds = shift;
    my $id = shift;
    if ($cmds =~ /-p (tmp_ped_\S+)/){
        $cmds .= "\n\nrm $1\n";
    }
    my $s = "$script_dir/$id.sh";
    my $runtime = 2;
    if ($id eq 'cadd'){
        $runtime = 12;
    }
    open (my $SCRIPT, ">", $s) or die "Could not write to $s: $!\n";
    print $SCRIPT <<EOT
 
#\$ -cwd
#\$ -V
#\$ -e $s.stderr
#\$ -o $s.stdout
#\$ -l h_rt=$runtime:00:00
#\$ -l h_vmem=8G
. /etc/profile.d/modules.sh

$cmds

EOT
;
    close $SCRIPT;
    return $s;
}
#####################################################################
sub executeCommands{
    if (exists $commands{cadd}){
        foreach my $cmd (@{$commands{cadd}}){
            execute($cmd);
        }
    }
    foreach my $k (keys %commands){
        next if $k eq 'cadd';
        foreach my $cmd (@{$commands{$k}}){
            execute($cmd);
        }
    }
}

#####################################################################
sub printCommands{
    if (exists $commands{cadd}){
        print join("\n", @{$commands{cadd}}) . "\n";
    }
    foreach my $k (keys %commands){
        next if $k eq 'cadd';
        print join("\n", @{$commands{$k}}) . "\n";
    }
}

#####################################################################

sub getFamilyFilterCommands{
    my $f = shift;
    my $ped = shift;
    my $in = shift;
    my @aff = $ped->getAffectedsFromFamily($f);
    ( my $out_stub = $in ) =~ s/\.vcf(\.(b)*gz)*$//; 
    my $output = "$out_stub". "_$f" . "_samplefilter.vcf";
    if (not @aff){
        warn "ERROR - No affecteds found for family $f!\n";
        return;
    }
    my @un = $ped->getUnaffectedsFromFamily($f);
    my $filter_cmd = 
        sprintf("$script_prefix/filterOnSample.pl -i %s -o %s -s %s ", 
                $in,
                $output,
                join(" ", @aff),
        );
    if (not $opts{d}){
        $filter_cmd .= " -x " . join(" ", @un); 
    }else{
        if ($opts{u}){ 
            $filter_cmd .= " -r \$(cat \"$opts{u}\") ";
        }
        if ($opts{x}){ 
            $filter_cmd .= " -x \$(cat \"$opts{x}\") ";
        }
    }
    $filter_cmd .= " $extra_args" if $extra_args; 
    $filter_cmd .= " --fork $opts{f} " if $opts{f};
    $filter_cmd .= " --frequency $opts{y} --min_alleles_for_freq 100 " if $opts{y};
    $filter_cmd .= " --progress " if $opts{b};
    return $filter_cmd, $output;
}
    
#####################################################################
sub getSegCommands{
    my $f = shift;
    my $ped = shift;
    my $u = $eddie ? 0 : 1;
    my ($TEMPPED, $temp_ped) = tempfile("tmp_ped_XXXXX", UNLINK => $u);
    foreach my $fam_samp ($ped->getSamplesFromFamily($f)){
        print $TEMPPED join
                    ("\t",
                            $f,
                            $fam_samp,
                            $ped->getFather($fam_samp),
                            $ped->getMother($fam_samp),
                            $ped->getSex($fam_samp),
                            $ped->getPhenotype($fam_samp),
                    ) . "\n";
    }
    close $TEMPPED;

    my @seg_commands = ();
    my @filtered = ();
    my @stubs = ();
    foreach my $in (@input){
        my ($f_cmd, $f_out) = getFamilyFilterCommands($f, $ped, $in);
        return if not $f_cmd;#no affecteds
        push @seg_commands, $f_cmd;
        push @filtered, $f_out;
        ( my $f_out_stub = $f_out ) =~ s/\.vcf(\.(b)*gz)*$//; 
        push @stubs, $f_out_stub; 
    }
    if (@stubs < 2){
        push @stubs, $stubs[0];
        push @filtered, $filtered[0];
    }
    my @aff = $ped->getAffectedsFromFamily($f);
    my @un = $ped->getUnaffectedsFromFamily($f);
    my $check_dominant = 0;
    foreach my $af (@aff){
        foreach my $p ($ped->getParents($af)){
            $check_dominant++ if $ped->isAffected($p);
        }
    }
    if ($check_dominant and not $opts{o}){
        my $seg_out = $stubs[1]. "_$f" . "_dominantfilter.vcf",
        print STDERR "Checking dominant inheritance in family $f...\n";
        my $seg_cmd = sprintf("$script_prefix/filterOnSample.pl -i %s -o %s -s %s ", 
                $filtered[1],
                $seg_out,
                join(" ", @aff),
        );
        if ($opts{d}){
            if (@un){
                $seg_cmd .= "-r " . join(" ", @un); 
            }
            if ($opts{u}){ 
                $seg_cmd .= " -r \$(cat \"$opts{u}\") ";
            }
            if ($opts{x}){ 
                $seg_cmd .= " -x \$(cat \"$opts{x}\") ";
            }
        }else{
            $seg_cmd .= " -x ";
        }
            
        $seg_cmd .= " $extra_args" if $extra_args; 
        $seg_cmd .= " --fork $opts{f} " if $opts{f};
        $seg_cmd .= " --progress " if $opts{b};
        push @seg_commands, $seg_cmd;
        my $rnk_cmd = "$script_prefix/rankOnCaddScore.pl -r ~/ref/GRCh37/cadd/v1.3/ -i $seg_out -o $seg_out.caddranked -n not_found.$seg_out.caddranked ";
        $rnk_cmd .= " --progress " if $opts{b};
        push @seg_commands, $rnk_cmd;
        my $ens_cmd = "$script_prefix/geneAnnotator.pl --functional -i $seg_out.caddranked -o $seg_out.caddranked.geneanno"; 
        push @seg_commands, $ens_cmd;
        my $fnc_cmd = "$script_prefix/getFunctionalVariants.pl --input $seg_out.caddranked.geneanno -o functional.$seg_out.caddranked.geneanno  --add splice_region_variant -d all -k --af 0.01 --pass ";
        push @seg_commands, $fnc_cmd;
        my $anno_cmd = "$script_prefix/annovcfToSimple.pl -v -g -u all -p $temp_ped -i functional.$seg_out.caddranked.geneanno "; 
        push @seg_commands, $anno_cmd;
    }elsif(not $opts{o}){
        print STDERR "Checking recessive inheritance in family $f...\n";
        my @out_files = ();
        my @bial_commands = ();
        my $seg_out ="biallelic_$f.$stubs[0].vcf";
        #create temp ped file of just this family
        push @bial_commands, sprintf("$script_prefix/findBiallelic.pl --af 0.01 -g 0.05 -c 10 --add splice_region_variant -i %s -o %s -l %s.genelist -f %s ", 
                                $filtered[0], 
                                $seg_out,
                                $seg_out,
                                $temp_ped,
                            );
        push @out_files, $seg_out;
        if ($opts{l}){
            my $lenient_out = "biallelic_lenient_$f.$seg_out.vcf";
            push @bial_commands, sprintf("$script_prefix/findBiallelic.pl --af 0.01 -c 10 --add splice_region_variant -i %s -o %s -l %s.genelist -f %s -t", 
                                    $filtered[0], 
                                    $lenient_out,
                                    $lenient_out,
                                    $temp_ped,
                                );
            push @out_files, $lenient_out; 
        }
        if ($opts{z}){
            my $denovo_b_out = "biallelic_denovo_$f.$seg_out.vcf";
            push @bial_commands, sprintf("$script_prefix/findBiallelic.pl --af 0.01 -g 0.05 -c 10 --add splice_region_variant -i %s -o %s -l %s.genelist -f %s --denovo_biallelic_mode", 
                                    $filtered[0], 
                                    $denovo_b_out,
                                    $denovo_b_out,
                                    $temp_ped,
                                );
            push @out_files, $denovo_b_out; 
        }
        foreach my $cmd (@bial_commands){
            if ($opts{d}){
                if ($opts{u}){ 
                    $cmd .= " -r \$(cat \"$opts{u}\") ";
                }
                if ($opts{x}){ 
                    $cmd .= " -x \$(cat \"$opts{x}\") ";
                }
            }else{
                $cmd .= " -x " ;
            }
            $cmd .= " $extra_args" if $extra_args; 
            $cmd .= " --progress " if $opts{b};
            push @seg_commands, $cmd;
        }
        foreach my $out (@out_files){
            my $ens_cmd = "$script_prefix/geneAnnotator.pl --functional -i $out -o $out.geneanno"; 
            push @seg_commands, $ens_cmd;
            my $anno_cmd = "$script_prefix/annovcfToSimple.pl -v -g -u all -p $temp_ped -i $out.geneanno"; 
            if ($out =~ /biallelic_(denovo|lenient)_/){
                $anno_cmd .= " -n hiConfDeNovo loConfDeNovo CaddPhredScore" ;
            }
            push @seg_commands, $anno_cmd;
        }
    }
    if (not $check_dominant and not $opts{s}){
        print STDERR "Checking de novo occurence in family $f...\n";
        my $denovo_cmd = sprintf("$script_prefix/filterOnSample.pl -i %s -o %s -s %s -r %s --confirm -z 0.2 ", 
                $filtered[1],
                "denovo_$f.$filtered[1]",
                join(" ", @aff),
                join(" ", @un),
        );
        $denovo_cmd .= " $extra_args" if $extra_args; 
        $denovo_cmd .= " --progress " if $opts{b};
        push @seg_commands, $denovo_cmd;
        my $rnk_cmd = "$script_prefix/rankOnCaddScore.pl -r ~/ref/GRCh37/cadd/v1.3/  -i denovo_$f.$filtered[1] -o denovo_$f.$filtered[1].caddranked -n not_found.denovo_$f.$filtered[1].caddranked ";
        $rnk_cmd .= " --progress " if $opts{b};
        push @seg_commands, $rnk_cmd;
        my $ens_cmd = "$script_prefix/geneAnnotator.pl --functional -i denovo_$f.$filtered[1].caddranked -o denovo_$f.$filtered[1].caddranked.geneanno"; 
        push @seg_commands, $ens_cmd;
        my $fnc_cmd = "$script_prefix/getFunctionalVariants.pl --input denovo_$f.$filtered[1].caddranked.geneanno -o functional.denovo_$f.$filtered[1].caddranked.geneanno  --add splice_region_variant -c 10 --af 0.01 --pass ";
        push @seg_commands, $fnc_cmd;
        my $anno_cmd = "$script_prefix/annovcfToSimple.pl -n hiConfDeNovo loConfDeNovo CaddPhredScore -v -g -u all -p $temp_ped -i functional.denovo_$f.$filtered[1].caddranked.geneanno"; 
        push @seg_commands, $anno_cmd;
    }
    if (not $check_dominant and not $opts{s} and $opts{g}){
        foreach my $d_type ( qw /hiConfDeNovo loConfDeNovo / ){ 
            my $gatk_cmd = "java -Xmx4g -jar $gatk -R $fasta -T SelectVariants -select \"vc.hasAttribute('$d_type') "; 
            foreach my $af (@aff){
                $gatk_cmd .= "&& vc.getAttribute('$d_type').contains('$af') ";
            }
            my $gatk_out =  join("_", @aff) . "_$d_type.vcf";
            $gatk_cmd .= "\" -V $filtered[1] -o $gatk_out";
            push @seg_commands, $gatk_cmd;
            my $rnk_cmd = "$script_prefix/rankOnCaddScore.pl -r ~/ref/GRCh37/cadd/v1.3/  -i $gatk_out -o $gatk_out.caddranked -n not_found.$gatk_out.caddranked ";
            $rnk_cmd .= " --progress " if $opts{b};
            push @seg_commands, $rnk_cmd;
            my $ens_cmd = "$script_prefix/geneAnnotator.pl --functional -i $gatk_out.caddranked -o $gatk_out.caddranked.geneanno"; 
            push @seg_commands, $ens_cmd;
            my $fnc_cmd = "$script_prefix/getFunctionalVariants.pl --input $gatk_out.caddranked.geneanno -o functional.$gatk_out.caddranked.geneanno  --add splice_region_variant -c 10 --af 0.01 --pass ";
            push @seg_commands, $fnc_cmd;
            my $anno_cmd = "$script_prefix/annovcfToSimple.pl -n hiConfDeNovo loConfDeNovo CaddPhredScore -v -g -u all -p $temp_ped -i functional.$gatk_out.caddranked.geneanno"; 
            push @seg_commands, $anno_cmd;
        }
    }   
    return @seg_commands;
}


