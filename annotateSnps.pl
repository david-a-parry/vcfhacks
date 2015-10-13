#!/usr/bin/env perl
use warnings;
use strict;
use Parallel::ForkManager;
use Getopt::Long qw(:config no_ignore_case);
use Sys::CPU;
use Pod::Usage;
use Term::ProgressBar;
use Data::Dumper;
use POSIX qw/strftime/;
use FindBin;
use lib "$FindBin::Bin/lib";
use VcfReader;
use ClinVarReader;

my @samples;
my @dbsnp ;
my $cvar;
my $freq;
my $quiet;
my $strict;
my $cpus = Sys::CPU::cpu_count();
my $forks = 0;
my $buffer_size;
my %opts = (
    CACHE      => \$buffer_size,
    forks      => \$forks,
    samples    => \@samples,
    dbsnp_file => \@dbsnp,
    freq       => \$freq,
    quiet      => \$quiet,
    strict     => \$strict
);

GetOptions(
    \%opts,    "output|o=s",     "input|i=s",         "known_snps|k=s",
    "replace|r", "samples|s=s{,}", "dbsnp_file|d=s{,}", "help|h|?",
    "manual",  "build|b=i", "clinvar_file|c=s",
    "f|freq=f"  => \$freq,
    "t|forks=i" => \$forks,
    "CACHE|C=i", "pathogenic", "quiet|q", "progress",
    "VERBOSE|V", "no_common_tag|n",
) or pod2usage( -exitval => 2, -message => "Syntax error" );
pod2usage( -verbose => 2 ) if $opts{manual};
pod2usage( -verbose => 1 ) if $opts{help};
pod2usage( -exitval => 2, -message => "Syntax error" )
  if not $opts{input}
  or not @{ $opts{dbsnp_file} };

if ( $forks < 2 ) {
    $forks = 0;    #no point having overhead of forks for one fork
}
else {
    if ($forks > $cpus){
        print STDERR "Warning - Number of forks ($forks) exceeds number of CPUs on this machine ($cpus)\n";
    }
    if ( not $buffer_size ) {
        $buffer_size = 10000 > $forks * 1000 ? 10000 : $forks * 1000;
    }
    print STDERR
"INFO - Processing in batches of $buffer_size variants split among $forks forks.\n";
}
if ( defined $freq ) {
    pod2usage(
        -exitval => 2,
        -message =>
          "value for --freq argument must be greater than 0 and less than 50"
    ) if $freq > 50 or $freq <= 0;
}
pod2usage(
    -exitval => 2,
    -message => "value for --build argument must be equal to or greater than 0"
) if defined $opts{build} && $opts{build} < 0;

my $OUT;
if ( $opts{output} ) {
    open( $OUT, ">$opts{output}" )
      || die "Can't open $opts{output} for writing: $!";
}
else {
    $OUT = *STDOUT;
}
my $KNOWN;
if ( $opts{known_snps} ) {
    open( $KNOWN, ">$opts{known_snps}" )
      || die "Can't open $opts{known_snps} for writing: $!";
}

if (not defined $opts{no_common_tag}){
    $opts{no_common_tag} = 0;
}
my $progressbar;
my $next_update      = 0;
my $kept             = 0;    #variants not filtered
my $filtered         = 0;    #variants filtered
my $pathogenic_snps  = 0;
my $found            = 0;    #variants that match a known SNP
my $printed_to_known = 0;
my $total_vcf        = 0;
my $time = strftime( "%H:%M:%S", localtime );
print STDERR "[$time] INFO - Initializing input VCF... ";

die "Header not ok for input ($opts{input}) "
    if not VcfReader::checkHeader( vcf => $opts{input} );
if ( defined $opts{progress} ) {
    $total_vcf = VcfReader::countVariants( $opts{input} );
    $time = strftime( "%H:%M:%S", localtime );
    print STDERR "\n[$time] INFO - $opts{input} has $total_vcf variants. ";
}
my %contigs = ();
if ($forks){
    %contigs       = VcfReader::getContigOrder( $opts{input} );
}
my %sample_to_col = ();
if (@samples) {
    %sample_to_col = VcfReader::getSamples(
        vcf         => $opts{input},
        get_columns => 1,
    );
}

$time = strftime( "%H:%M:%S", localtime );
print STDERR "\n[$time] INFO - Finished initializing input VCF\n";

if ($opts{clinvar_file}){
    $time = strftime( "%H:%M:%S", localtime );
    print STDERR "[$time] INFO - Checking ClinVar file, $opts{clinvar_file}...\n";
    my ( $cv, $col ) = ClinVarReader::checkClinVarFile($opts{clinvar_file});
    $opts{clinvar_file} = $cv;
    $time = strftime( "%H:%M:%S", localtime );
    print STDERR "[$time] INFO - ClinVar file is OK.\n";
}

if (@dbsnp > 1){
    $time = strftime( "%H:%M:%S", localtime );
    print STDERR "[$time] INFO - Multiple dbSNP files in use - "
      . "annotations will be added in following order of precedence: "
      . join (" > ", @dbsnp) . "\n";
}

my %dbsnp_to_info;#key is file, value is ref to hash of info fields
my %dbsnp_to_index = ();
my $dbpm           = Parallel::ForkManager->new($forks);
for ( my $i = 0 ; $i < @dbsnp ; $i++ ) {
    $dbpm->run_on_finish(    # called BEFORE the first call to start()
        sub {
            my ( $pid, $exit_code, $ident, $exit_signal, $core_dump,
                $data_structure_reference )
              = @_;

            if ( defined($data_structure_reference) ){
                if ( ref $data_structure_reference eq 'ARRAY' ) {
                    $dbsnp_to_info{$data_structure_reference->[2]} = $data_structure_reference->[0];
                    $dbsnp_to_index{ $data_structure_reference->[2] } =
                      $data_structure_reference->[1];
                }
                else {
                    die
                      "Unexpected return from dbSNP reference initialization:\n"
                      . Dumper $data_structure_reference;
                }
            }
            else{ # problems occuring during storage or retrieval 
                die "No message received from child process $pid!\n";
            }
        }
    );
    $time = strftime( "%H:%M:%S", localtime );
    print STDERR "[$time] INFO - Initializing $dbsnp[$i] dbSNP reference VCF "
      . ( $i + 1 ) . " of "
      . scalar(@dbsnp) . "\n";
    $dbpm->start() and next;
    my @info_and_index = initializeDbsnpVcfs( $dbsnp[$i] );
    push @info_and_index, $dbsnp[$i];
    $time = strftime( "%H:%M:%S", localtime );
    print STDERR
      "[$time] INFO - Finished initializing $dbsnp[$i] dbSNP reference VCF.\n";
    $dbpm->finish( 0, \@info_and_index );
}
print STDERR "INFO - Waiting for children...\n" if $opts{VERBOSE};
$dbpm->wait_all_children;
my @add_head = ();
if ( $opts{pathogenic} ) {
    my $inf = checkAndAddHeaders("CLNSIG");
    push @add_head, $inf if $inf;
    $inf = checkAndAddHeaders("SCS");
    push @add_head, $inf if $inf;
    if ( not @add_head and not $opts{clinvar_file}) {
        print STDERR
"WARNING - can't find CLNSIG or SCS fields in dbSNP file headers, your SNP reference files probably don't have pathogenic annotations.\n";
    }
}

if ( $opts{build} ) {
    my $inf = checkAndAddHeaders("dbSNPBuildID");
    if ( not $inf ) {
        print STDERR
"WARNING - can't find dbSNPBuildID fields in dbSNP file headers, your SNP reference files probably don't have readable dbSNP build annotations.\n";
    }
    else {
        push @add_head, $inf;
    }
}
if ( $opts{freq} ) {
    my @freq_head;
    foreach my $f ( qw / AF CAF G5A G5 COMMON / ) { 
        my $inf = checkAndAddHeaders("$f");
        push @freq_head, $inf if $inf;
    }
    if ( not @freq_head ) {
        print STDERR
"WARNING - can't find allele frequency fields (GMAF, CAF, AF, G5A, G5 or COMMON) in dbSNP file headers, your SNP reference files probably don't have readable frequency data.\n";
    }
    else {
        push @add_head, @freq_head;
    }
}

my $meta_head = VcfReader::getMetaHeader( $opts{input} );
print $OUT "$meta_head\n";
print $KNOWN "$meta_head\n" if $KNOWN;
my $headstring = '';
if (@add_head) {
    $headstring .= join("\n", @add_head) . "\n";
}
if ($opts{clinvar_file}){
    $headstring .= <<EOT
##INFO=<ID=ClinVarPathogenic,Number=A,Type=Integer,Description="For each allele, a value of 1 is given if the variant has ever been asserted "Pathogenic" or "Likely pathogenic" by any submitter for any phenotype, and 0 if present in ClinVar but does not meet this criteria">
##INFO=<ID=ClinVarConflicted,Number=A,Type=Integer,Description="For each allele, a value of 1 is given if the variant has ever been asserted "Pathogenic" or "Likely pathogenic" by any submitter for any phenotype, and has also been asserted "Benign" or "Likely benign" by any submitter for any phenotype, and 0 if present in ClinVar but does not meet this criteria">
##INFO=<ID=ClinVarTraits,Number=A,Type=String,Description="Any traits associated with this allele in ClinVar">
##INFO=<ID=ClinVarClinicalSignificance,Number=A,Type=String,Description="Clinical significance terms for this allele given in ClinVar">
EOT
    ;
}
if ($headstring){
    print $OUT $headstring;
    print $KNOWN $headstring if $KNOWN;
}

print $OUT "##annotateSnps.pl=\"";
print $KNOWN "##annotateSnps.pl=\"" if $KNOWN;

my @opt_string = ();
foreach my $k ( sort keys %opts ) {
    if ( not ref $opts{$k} ) {
        push @opt_string, "$k=$opts{$k}";
    }
    elsif ( ref $opts{$k} eq 'SCALAR' ) {
        if ( defined ${ $opts{$k} } ) {
            push @opt_string, "$k=${$opts{$k}}";
        }
        else {
            push @opt_string, "$k=undef";
        }
    }
    elsif ( ref $opts{$k} eq 'ARRAY' ) {
        if ( @{ $opts{$k} } ) {
            push @opt_string, "$k=" . join( ",", @{ $opts{$k} } );
        }
        else {
            push @opt_string, "$k=undef";
        }
    }
}
my $col_header = VcfReader::getColumnHeader( $opts{input} );
print $OUT join( " ", @opt_string ) . "\"\n" . $col_header . "\n";
print $KNOWN join( " ", @opt_string ) . "\"\n" . $col_header . "\n"
  if $KNOWN;

$time = strftime( "%H:%M:%S", localtime );
print STDERR "[$time] INFO - SNP annotation starting\n";

$freq /= 100 if ($freq);
if ( $opts{build} || $freq ) {
    $time = strftime( "%H:%M:%S", localtime );
    print STDERR "[$time] INFO - Filtering variants on following criteria:\n";
    print STDERR ">in dbSNP$opts{build} or previous\n"
      if defined $opts{build} && $opts{build};
    print STDERR
      ">with a minor allele frequency equal to or greater than $freq ("
      . $freq * 100 . " %)\n"
      if defined $freq && $freq;
    print STDERR
      ">without \"pathogenic\" or \"probably pathogenic\" annotations\n"
      if defined $opts{pathogenic} && $opts{pathogenic};
}
elsif ( $opts{pathogenic} ) {
    if ($KNOWN) {
        print STDERR
"INFO - Annotating output and printing variants with \"pathogenic\" or \"probably pathogenic\" annotations to $opts{known_snps}\n";
    }
    else {
        print STDERR
"WARNING: --pathogenic flag acheives nothing if neither --build or --freq filtering is in effect and --known_out output is not specified\n";
    }
}
my $prog_total;
if ( defined $opts{progress} and $total_vcf ) {
    my $x_prog = 3;
    $x_prog = 4 if $KNOWN;
    $prog_total = $total_vcf * $x_prog;
    if ( $opts{build} || $freq ) {
        $progressbar = Term::ProgressBar->new(
            { name => "Filtering", 
              count => ($prog_total), 
              ETA => "linear" 
            } 
        );
    }
    else {
        $progressbar = Term::ProgressBar->new(
            #{ name => "Annotating", count => ($prog_total), ETA => "linear", }
            { name => "Annotating", 
              count => ($prog_total), 
              ETA => "linear" 
            } 
        );
    }

    #$progressbar->minor(0);
}

my $n                = 0;
my @lines_to_process = ();
my $VCF              = VcfReader::_openFileHandle( $opts{input} );
my %no_fork_args = ();
my %no_fork_cvar_args    = ();

if ($forks < 2){
    foreach my $d (@dbsnp) {
        my %s = VcfReader::getSearchArguments( $d, $dbsnp_to_index{$d} );
        $no_fork_args{$d} = \%s;
    }
    if ($opts{clinvar_file}){
        %no_fork_cvar_args = ClinVarReader::getClinVarSearchArgs($opts{clinvar_file});
    }
}
VAR: while ( my $line = <$VCF> ) {
    next if $line =~ /^#/;
    $n++;
    if ($progressbar) {
        $next_update = $progressbar->update($n) if $n >= $next_update;
    }
    if ($forks > 1){
        push @lines_to_process, $line;
        if ( @lines_to_process >= $buffer_size ) {
            process_buffer();
            @lines_to_process = ();
        }
    }else{
        chomp $line;
        #our VcfReader methods should be more efficient on pre-split lines
        my @split_line = split( "\t", $line );
        my %res = filterSnps( \@split_line, \%no_fork_args, \%no_fork_cvar_args);
        if ($res{keep}){
            print $OUT join("\t", @{$res{keep}}) ."\n"; 
            $kept++;
        }
        $filtered++ if $res{filter};
        $found++ if $res{found};
        $pathogenic_snps++ if  $res{pathogenic} ;
        $n += 2;
        if ($KNOWN){
            $n++;
             if ($res{known}){
                print $KNOWN join("\t", @{$res{known}}) ."\n"; 
                $printed_to_known++;
            }
        }
        if ($progressbar) {
            $next_update = $progressbar->update($n) if $n >= $next_update;
        }
    }
}
close $VCF;
process_buffer() if $forks > 1;
close $OUT;
close $KNOWN if $KNOWN;
if ( defined $opts{progress} ) {
    $progressbar->update($prog_total) if $prog_total >= $next_update;
}
$time = strftime( "%H:%M:%S", localtime );
print STDERR "\nTime finished: $time\n";
print STDERR "$found known SNPs identified.\n";
print STDERR "$filtered variants filtered, $kept variants retained.\n";
print STDERR
  "$pathogenic_snps pathogenic or probable pathogenic variants identified.\n"
  if $pathogenic_snps;
print STDERR 
  "$printed_to_known variants matching user-specified criteria printed to $opts{known_snps}.\n" 
    if $KNOWN;

################################################
#####################SUBS#######################
################################################

sub process_buffer {

    #this array is an array of refs to batches
    # so that we can quickly sort our batches rather than
    # performing a big sort on all lines
    return if not @lines_to_process;
    my @lines_to_print;
    my @known;
    my $i               = 0;
    my $t               = 0;
    my $lines_per_slice = @lines_to_process;
    if ( $forks > 0 ) {
        $lines_per_slice = int( @lines_to_process / $forks ) > 1 ? int( @lines_to_process / $forks ) : 1;
    }
    my @batch = ();

    #get a batch for each thread
    for ( my $i = 0 ; $i < @lines_to_process ; $i += $lines_per_slice ) {
        my $last =
          ( $i + $lines_per_slice - 1 ) < $#lines_to_process
          ? $i + $lines_per_slice - 1
          : $#lines_to_process;
        if ($i + $lines_per_slice >= @lines_to_process){
            $last = $#lines_to_process;
        }
        my @temp = @lines_to_process[ $i .. $last ];
        push @batch, \@temp;
    }
    my $pm = Parallel::ForkManager->new($forks);
    $pm->run_on_finish(    # called BEFORE the first call to start()
        sub {
            my ( $pid, $exit_code, $ident, $exit_signal, $core_dump,
                $data_structure_reference )
              = @_;

            if ( defined($data_structure_reference) ){
                my %res = %{$data_structure_reference}
                  ;        
                if ( ref $res{keep} eq 'ARRAY' ) {
                    push @lines_to_print, \@{ $res{keep} } if @{ $res{keep} };
                }
                if ( $res{filter} ) {
                    $filtered += $res{filter} ;
                }
                if ( $res{pathogenic} ) {
                    $pathogenic_snps += $res{pathogenic} ;
                }
                if ( ref $res{known} eq 'ARRAY' ) {
                    push @known, \@{ $res{known} } if @{ $res{known} };
                }
                if ( $res{found}  ) {
                    $found += $res{found} ;
                }
                if ($progressbar) {
                    $n += $res{batch_size};
                    $next_update = $progressbar->update($n)
                      if $n >= $next_update;
                }
            }else{
                die "ERROR: no message received from child process $pid!\n";
            }
        }
    );
    foreach my $b (@batch) {
        $pm->start() and next;
        my %results = process_batch($b);
        $pm->finish( 0, \%results );
    }
    $pm->wait_all_children;

    #print them
    if (@lines_to_print){
        @lines_to_print = sort {
            VcfReader::by_first_last_line($a, $b, \%contigs) 
            } @lines_to_print;
        my $incr_per_batch = @lines_to_process / @lines_to_print;
        foreach my $batch (@lines_to_print) {
            my $incr_per_line = $incr_per_batch / @$batch;
            foreach my $l (@$batch) {
                if ( ref $l eq 'ARRAY' ) {
                    print $OUT join( "\t", @$l ) . "\n";
                }
                else {
                    print $OUT "$l\n";
                }
                $kept++;
                if ($progressbar) {
                    $n += $incr_per_line;
                    $next_update = $progressbar->update($n) if $n >= $next_update;
                }
            }
        }
    }else{
        if ($progressbar) {
            $n += @lines_to_process;
            $next_update = $progressbar->update($n) if $n >= $next_update;
        }
    }
    if ($KNOWN) {
        if (@known){
            @known = sort {
                VcfReader::by_first_last_line($a, $b, \%contigs) 
                } @known;
            my $incr_per_batch = @lines_to_process / @known;
            foreach my $k (@known) {
                my $incr_per_line = $incr_per_batch / @$k;
                foreach my $l (@$k) {
                    if ( ref $l eq 'ARRAY' ) {
                        print $KNOWN join( "\t", @$l ) . "\n";
                    }
                    else {
                        print $KNOWN "$l\n";
                    }
                    $printed_to_known++;
                    if ($progressbar) {
                        $n += $incr_per_line;
                        $next_update = $progressbar->update($n)
                        if $n >= $next_update;
                    }
                }
            }
        }else{
            if ($progressbar) {
                $n += @lines_to_process;
                $next_update = $progressbar->update($n) if $n >= $next_update;
            }
        }
    }
    else {
        foreach my $batch (@known) {
            $printed_to_known += @$batch;
        }
    }
}

################################################
sub process_batch {

    #filter a set of lines
    my ($batch) = @_;
    my %results = ( batch_size => scalar(@$batch) );
    my %sargs;
    my %cvargs;
    foreach my $d (@dbsnp) {

        #WE'VE FORKED AND COULD HAVE A RACE CONDITION HERE -
        # WHICH IS WHY WE DO A FIRST PASS WITH OUR initializeDbsnpVcfs
        # METHOD TOWARDS THE START OF THE PROGRAM
        my %s = VcfReader::getSearchArguments( $d, $dbsnp_to_index{$d} );
        $sargs{$d} = \%s;
    }
    if ($opts{clinvar_file}){
        %cvargs = ClinVarReader::getClinVarSearchArgs($opts{clinvar_file});
    }
    foreach my $line ( @{$batch} ) {
        chomp $line;
        #our VcfReader methods should be more efficient on pre-split lines
        my @split_line = split( "\t", $line );
        my %res = filterSnps( \@split_line, \%sargs, \%cvargs );
        push @{ $results{keep} },   $res{keep}   if $res{keep};
        $results{filter}++  if $res{filter};
        push @{ $results{known} },  $res{known}  if $res{known};
        $results{pathogenic}++  if $res{pathogenic};
        $results{found}++  if $res{found};
    }
    foreach my $d (keys %sargs){
        if (exists $sargs{$d}->{file_handle}){
            close $sargs{$d}->{file_handle} if $sargs{$d}->{file_handle};
        }
    }
    return %results;
}

################################################
sub filterSnps {
    my ( $vcf_line, $search_args, $cvar_args ) = @_;
    my ( $keep, $filter, $known, $pathogenic );
    my $line_should_not_be_filtered = 0
      ; #flag in case pathogenic flag is set or something and we don't want to filter this snp no matter what
    my $is_known_snp     = 0;
    my $replaced_already = 0
      ; #use this flag to allow appending of multiple IDs found in our reference file even if --replace option is in place.

#process each allele separately in case we have MNVs/deletions that need to be simplified
#
    my %min_vars       = VcfReader::minimizeAlleles($vcf_line);
    my %sample_alleles = ();
    if (@samples) {
        %sample_alleles = map { $_ => undef } VcfReader::getSampleCall(
            multiple            => \@samples,
            sample_to_columns   => \%sample_to_col,
            return_alleles_only => 1,
            line                => $vcf_line
        );
    }
    my @info;
    foreach my $allele ( sort { $a <=> $b } keys %min_vars ) {
        if ($min_vars{$allele}->{ORIGINAL_ALT} eq '*'){
            #filter new asterisk allele notation
            $min_vars{$allele}->{filter_snp} = 1;
            next;
        }
        $min_vars{$allele}->{filter_snp} = 0;
        if (@samples) {

            #doesn't exist in any sample...
            next if not exists $sample_alleles{$allele};
        }
        foreach my $k ( @dbsnp ) {#do each file in order so that annotations from file 1 take precedence over others
        #foreach my $k ( keys %{$search_args} ) {
            if (
                my @snp_hits = VcfReader::searchForPosition(
                    %{ $search_args->{$k} },

                    #vcf => $k,
                    chrom => $min_vars{$allele}->{CHROM},
                    pos   => $min_vars{$allele}->{POS}
                )
              )
            {
                foreach my $snp_line (@snp_hits) {

                    #check whether the snp line(s) match our variant
                    my @snp_split = split( "\t", $snp_line );
                    if ( my $match = checkVarMatches( $min_vars{$allele}, \@snp_split ) ) {
                        $is_known_snp++;

                        #replace or append to ID field
                        if (
                            VcfReader::getVariantField( $vcf_line, 'ID' ) eq '.'
                            || ( $opts{replace} && not $replaced_already ) )
                        {
                            $replaced_already++;
                            $vcf_line = VcfReader::replaceVariantField(
                                $vcf_line,
                                'ID',
                                VcfReader::getVariantField( \@snp_split, 'ID' )
                            );
                        }
                        else {
                            my $id =
                              VcfReader::getVariantField( $vcf_line, 'ID' );
                            my @split_id = split( ";", $id );

                            #check it's not already correctly annotated
                            my $new_id = $id;
                            foreach my $sid (
                                split(
                                    ";",
                                    VcfReader::getVariantField(
                                        \@snp_split, 'ID'
                                    )
                                )
                              )
                            {
                                if ( not grep { /^$sid$/i } @split_id ) {
                                    $new_id .= ";" . $sid;
                                }
                            }
                            if ( $new_id ne $id ) {
                                $vcf_line =
                                  VcfReader::replaceVariantField( $vcf_line,
                                    'ID', $new_id );
                            }
                        }

                        #get snp info and perform 
                        #filtering if fiters are set
                        $line_should_not_be_filtered += 
                          evaluate_snp( $min_vars{$allele},  \@snp_split, $match );
                    }
                }
            }
        }
        if (%$cvar_args){#perform ClinVar filtering with ClinVar.tsv file
            $line_should_not_be_filtered += evaluate_clinvar( $min_vars{$allele}, $cvar_args );
        }
    }
    my $filter_count = 0;
    my %snp_info = ();
    my %clinvar_info = ();
    #TO DO
    #Get all snp_info fields present in %min_vars and create undef hash entries
    foreach my $allele ( keys %min_vars ) {
        map { $snp_info{$_} = undef }  keys %{ $min_vars{$allele}->{snp_info} } ;
    }
    #then work through each allele, adding value to array or '.' in if field is unavailable
    foreach my $allele ( sort { $a <=> $b } keys %min_vars ) {
        foreach my $k (keys %snp_info){
            if (exists $min_vars{$allele}->{snp_info}->{$k}){
                push @{$snp_info{$k}}, $min_vars{$allele}->{snp_info}->{$k};
            }else{  
                push @{$snp_info{$k}}, ".";
            }
        }
        #add clinvar info for allele
        foreach my $k ( keys %{ $min_vars{$allele}->{cvar_info} } ) {
            push @{$clinvar_info{$k}}, $min_vars{$allele}->{cvar_info}->{$k};
        }
    }
    foreach my $k ( sort keys %snp_info ) {
        $vcf_line = VcfReader::addVariantInfoField 
            (
                line => $vcf_line,
                id   => "AS_$k",
                value => join(",", @{ $snp_info{$k} } ), 
            );
    }
    foreach my $k ( sort keys %clinvar_info ) {
        $vcf_line = VcfReader::addVariantInfoField 
            (
                line => $vcf_line,
                id   => $k,
                value => join(",", @{ $clinvar_info{$k} } ), 
            );
    }
    foreach my $allele ( keys %min_vars ) {
        $filter_count++ if ( $min_vars{$allele}->{filter_snp} );
    }
    if ( $opts{build} || $opts{pathogenic} || $freq )
    {    #if using filtering only put filtered in $KNOWN
        if ($line_should_not_be_filtered)
        { #at the moment this means it's pathogenic flagged and $opts{pathogenic} is set
             #so if --build or --freq are set keep, if not send to both @$keep and @$known
            $keep = $vcf_line;
            $pathogenic = $vcf_line;
            if ( not $opts{build} && not $freq ) {
                $known = $vcf_line;
            }
        }
        elsif ( $filter_count == keys(%min_vars) )
        {    #filter if all alleles meet criteria
            $filter = $vcf_line;
            $known  = $vcf_line;
        }
        else {    # don't filter
            $keep = $vcf_line;
        }
    }
    else
    { #otherwise put all identified dbSNP variants in $@$known and all lines in @$keep
        $keep = $vcf_line;
        if ($is_known_snp) {
            $known = $vcf_line;
        }
    }

    #print STDERR "DEBUG: RETURNING FROM THREAD $thread\n" if $opts{VERBOSE};
    return ( keep => $keep, known => $known, found => $is_known_snp, filter => $filter, pathogenic => $pathogenic );
}

################################################
sub initializeDbsnpVcfs {
    my ($snpfile) = @_;
    #we getSearchArguments here simply to prevent a race condition later
    my @head = VcfReader::getHeader($snpfile);
    die "Header not ok for $snpfile " 
        if not VcfReader::checkHeader( header => \@head );
    my %info = VcfReader::getInfoFields( header => \@head);
    #my %sargs = VcfReader::getSearchArguments($snpfile);
    my %index = VcfReader::readIndex($snpfile);
    return ( \%info, \%index );
}

#################################################
sub evaluate_clinvar {
    my $min_allele = shift;
    my $cvar_args = shift;
    my @matches = ClinVarReader::searchForMatchingVariant
    (
        %$cvar_args,
        chrom => $min_allele->{CHROM},
        pos   => $min_allele->{POS},
        ref   => $min_allele->{REF},
        alt   => $min_allele->{ALT},
    );
    my $isPathogenic = 0;
    my $isConflicted = 0;
    my @traits = ();
    my @sig = ();
    my @path = ();
    my @conf = ();
    foreach my $m (@matches){
        push @traits, ClinVarReader::getColumnValue($m, 'all_traits', $cvar_args->{col_hash});
        push @sig, ClinVarReader::getColumnValue($m, 'clinical_significance', $cvar_args->{col_hash});
        my $c =  ClinVarReader::getColumnValue($m, 'conflicted', $cvar_args->{col_hash});
        my $p = ClinVarReader::getColumnValue($m, 'pathogenic', $cvar_args->{col_hash});
        push @path, $p;
        push @conf, $c;
        $isPathogenic += $p;
        $isConflicted += $c;
    }
    $min_allele->{cvar_info}->{ClinVarPathogenic} = $isPathogenic ? 1 : 0;
    $min_allele->{cvar_info}->{ClinVarConflicted} = $isConflicted ? 1 : 0;
    if (@traits){
        @traits = map { VcfReader::convertTextForInfo($_) } @traits;
        $min_allele->{cvar_info}->{ClinVarTraits} = join("/", @traits);
    }
    if (@sig){
        @sig = map { VcfReader::convertTextForInfo($_) } @sig;
        $min_allele->{cvar_info}->{ClinVarClinicalSignificance} = join("/", @sig);
    }
    return $isPathogenic;
}
#################################################
sub evaluate_snp {

#returns 1 if shouldn't be filtered under any circumstance (at the moment only if
#pathogenic flag is set and snp is clinically associated) and 0 otherwise
    my ( $min_allele, $snp_line, $snp_alt ) = @_;
    my %info_values = ();
    my @al =  VcfReader::readAlleles(line => $snp_line);
    foreach my $f (qw (SCS CLNSIG dbSNPBuildID G5 G5A GMAF CAF AF COMMON)) {
        my $value = VcfReader::getVariantInfoField( $snp_line, $f );
        if ( defined $value ) {
            if (not exists $min_allele->{snp_info}->{$f} 
              or $min_allele->{snp_info}->{$f} eq '.' ){ 
              #we only keep the first field we come accross in our dbSNP files
                if ($f eq 'CAF'){#we just take the alt value for CAF
                    my @caf = split(",", $value); 
                    $min_allele->{snp_info}->{$f} = $caf[$snp_alt];
                    $info_values{$f} = $caf[$snp_alt];
                }else{
                    #Only add COMMON/G5/G5A if there is only one ALT allele for this line
                    #Otherwise it is not clear whether match is for our alt allele
                    if ($f =~ /^(G5|COMMON|GMAF)/){
                        next if (@al > 2);
                    }
                    $min_allele->{snp_info}->{$f} = $value;
                    $info_values{$f} = $value;
                }
            }
        }
    }
    if ($opts{pathogenic}) {
        if ( exists $info_values{SCS} ) {
            my @scs = split( /[\,\|]/, $info_values{SCS} );
            foreach my $s (@scs) {

#SCS=4 indicates probable-pathogenic, SCS=5 indicates pathogenic, print regardless of freq or build if --pathogenic option in use
                if ( $s eq '4' or $s eq '5' ) {
                    return 1;
                }
            }
        }
        else {
#die "No SCS field in snp line $snp_line" if $strict; #the SCS field appears to be dropped from the GATK bundle now
#print STDERR "Warning, no SCS field in snp line $snp_line\n" unless $quiet;
        }
        if ( exists $info_values{CLNSIG} ) {
            my @scs = split( /[\,\|]/, $info_values{CLNSIG} );
            foreach my $s (@scs) {

#SCS=4 indicates probable-pathogenic, SCS=5 indicates pathogenic, print regardless of freq or build if --pathogenic option in use
                if ( $s eq '4' or $s eq '5' ) {
                    return 1;
                }
            }
        }
    }
    if ($opts{build}) {
        if ( exists $info_values{dbSNPBuildID} ) {
            if ( $opts{build} >= $info_values{dbSNPBuildID} ) {
                $min_allele->{filter_snp}++;
                return 0 ;
            }
        }
        else {
            die "No dbSNPBuildID field in snp line "
              . join( "\t", @$snp_line ) . "\n"
              if $strict;
            print STDERR "Warning, no dbSNPBuildID field in snp line "
              . join( "\t", @$snp_line ) . "\n"
              unless $quiet;
        }
    }
    if ($freq) {
        if ( $freq <= 0.05 ) {

            #G5 = minor allele freq > 5 % in at least 1 pop
            #G5A = minor allele freq > 5 % in all pops
            if ( $info_values{G5} ) {
                $min_allele->{filter_snp}++;
                return 0;
            }
            if ( $info_values{G5A} ) {
                $min_allele->{filter_snp}++;
                return 0;
            }
        }
        if ( $freq <= 0.01 && not $opts{no_common_tag}) {
              #only use COMMON tag if user has specifically asked for it 
              #as it seems to be quite inaccurate (small n?)
            if ( $info_values{COMMON} ) {
                $min_allele->{filter_snp}++;
                return 0;
            }
        }
        if ( exists $info_values{GMAF} ) {
            if ( $freq <= $info_values{GMAF} ) {
                $min_allele->{filter_snp}++;
                return 0;
            }
        }
        if ( exists $info_values{CAF} ) {
            if ($info_values{CAF} ne '.'){
                if ( $freq <= $info_values{CAF} ) {
                    $min_allele->{filter_snp}++;
                    return 0;
                }
            }
        }
        if ( exists $info_values{AF} ) {
            my $c = $info_values{AF};
            $c =~ s/^\[//;
            $c =~ s/\]$//;
            my @af = split( ',', $c );
            my %snp_min = VcfReader::minimizeAlleles($snp_line);
            foreach my $al ( keys %snp_min ) {
                next if $min_allele->{CHROM} ne $snp_min{$al}->{CHROM};
                next if $min_allele->{POS} ne $snp_min{$al}->{POS};
                next if $min_allele->{REF} ne $snp_min{$al}->{REF};
                next if $min_allele->{ALT} ne $snp_min{$al}->{ALT};
                die "AF values don't match no. alternate alleles for "
                  . " SNP line:\n"
                  . join( "\t", @$snp_line ) . "\n"
                  if ( @af < $al );
                next if $af[ $al - 1 ] eq '.';
                if ( $freq <= $af[ $al - 1 ] ) {
                    $min_allele->{filter_snp}++;
                    return 0;
                }
            }
        }
    }
    return 0;
}
#################################################
sub checkVarMatches {
    my ( $min_allele, $snp_line ) = @_;
    my %snp_min = VcfReader::minimizeAlleles($snp_line);
    foreach my $snp_allele ( keys %snp_min ) {
        $min_allele->{CHROM} =~ s/^chr//;
        $snp_min{$snp_allele}->{CHROM} =~ s/^chr//;
        next if $min_allele->{CHROM} ne $snp_min{$snp_allele}->{CHROM};
        next if $min_allele->{POS} ne $snp_min{$snp_allele}->{POS};
        next if $min_allele->{REF} ne $snp_min{$snp_allele}->{REF};
        next if $min_allele->{ALT} ne $snp_min{$snp_allele}->{ALT};
        return $snp_allele;
    }
    return 0;
}

#################################################

sub checkAndAddHeaders{
    my $field = shift; 
    my $inf_head;
    foreach my $d (@dbsnp){
        $time = strftime( "%H:%M:%S", localtime );
        if ( $dbsnp_to_info{$d}->{$field} ){
            print STDERR "[$time] INFO - $field field found in $d...\n";
            (my $desc = $dbsnp_to_info{$d}->{$field}->{Description}) =~ s/\"//;
       
            $inf_head = "##INFO=<ID=AS_$field,Number=A,Type=$dbsnp_to_info{$d}->{$field}->{Type},Description=\"This annotation has been altered by annotateSnps.pl to ".
                        "report consequences per ALT allele. Original description was as follows: $desc>";
        }
    }
    return $inf_head;
}

#################################################


=head1 NAME

annotateSnps.pl - annotate and optionally filter SNPs from a VCF file 

=head1 SYNOPSIS

        annotateSnps.pl -i [vcf file] -d [reference SNP vcf file] [options]
        annotateSnps.pl -h (display help message)
        annotateSnps.pl -m (display manual page)

=cut

=head1 ARGUMENTS

=over 8

=item B<-i    --input>

Input VCF file 

=item B<-o    --output>

Output snp annotated/filtered file. Optional - default is STDOUT.

=item B<-k    --known_out>

Optional output file to print identified SNPs to.  If --build or --freq arguments are in use only SNPs meeting filter criteria will be printed to this file.

=item B<-d    --dbsnp_file>

SNP reference VCF file(s). IDs from these files will be used to annotate/filter the input VCF file.  If an ID already exists IDs from matching variants will be appended to the ID field. Your dbSNP file MUST USE THE SAME REFERENCE as your input VCF.

SNP vcf files for use can be downloaded from the NCBI FTP site (ftp://ftp.ncbi.nih.gov/snp/) or from the Broad Institutes FTP site (e.g. ftp://ftp.broadinstitute.org/bundle/1.5/b37/). Clinically annotated VCFs are available from the NCBI FTP site. If more than one dbSNP file is used annotations from the first file will take precedence over subsequent files if the same INFO fields exist for a given variant. 

Your are advised to use bgzip compressed and tabix indexed files as your SNP reference VCFs.

=item B<-c    --clinvar_file>

Optional ClinVar TSV file from to identify pathogenic variants as obtained from https://github.com/macarthur-lab/clinvar. If used the following INFO fields will be added:

=over 8
    
=item ClinVarPathogenic 

For each allele, a value of 1 is given if the variant has ever been asserted "Pathogenic" or "Likely pathogenic" by any submitter for any phenotype, and 0 if present in ClinVar but does not meet this criteria

=item ClinVarConflicted 

For each allele, a value of 1 is given if the variant has ever been asserted "Pathogenic" or "Likely pathogenic" by any submitter for any phenotype, and has also been asserted "Benign" or "Likely benign" by any submitter for any phenotype, and 0 if present in ClinVar but does not meet this criteria

=item ClinVarTraits 

Any traits associated with an allele in ClinVar

=item ClinVarClinicalSignificance

Clinical significance terms for alleles given in ClinVar

=back

=item B<-r    --replace>

Use this option to replace the ID field with IDs from your SNP reference VCF file rather than appending to existing IDs.

=item B<-b    --build>

Build number to filter from (e.g. 129).  SNPs from this build or before will be filtered from output regardless of any value given for --freq. Must be an integer.

=item B<-f    --freq>

Percent SNP minor allele frequency to filter from. SNPs with minor alleles equal to or over this frequency will be removed from output regardless of value given for --build.

=item B<--pathogenic>

When using --build or --freq filtering use this flag to print SNPs with "pathogenic" or "probably pathogenic" annotations (as indicated by SCS or CLNSIG flags in the dbSNP/ClinVar VCF) to your filtered output regardless of build or frequency settings.  If used on its own the program will only print SNPs with "pathogenic" or "probably pathogenic" annotations to your --known_out file (if specified).

=item B<-n    --no_common_tag>

If filtering on a minor allele frequency of 1 % or lower, use of this flag will disable the 'COMMON' tag for variant filtering. The 'COMMON' tag in dbSNP denotes variants that have "at least one 1000Genomes population with a minor allele of frequency >= 1% and for which 2 or more founders contribute to that minor allele frequency", however due to the relatively small sample sizes it might be desirable in some cases not to use this annotation. 

=item B<-s    --samples>

One or more samples to check variants for.  Default is to check all variants specified by the ALT field. If specified, SNPs will not be filtered unless variant is present in at least one of the samples specified by this argument.

=item B<-t    --forks>

Number of forks to create for parallelising your analysis. By default no forking is done. To speed up your analysis you may specify the number of parallel processes to use here. (N.B. forking only occurs if a value of 2 or more is given here as creating 1 fork only results in increased overhead with no performance benefit).

=item B<-C    --CACHE>

Cache size. Variants are processed in batches to allow for efficient parallelisation. When forks are used the default is to process up to 10,000 variants at once or 1,000 x no. forks if more than 10 forks are used. If you find this program comsumes too much memory when forking you may want to set a lower number here. When using forks you may get improved performance by specifying a higher cache size, however the increase in memory usage is proportional to your cache size multiplied by the number of forks.

=item B<--progress>

Use this flag to show a progress bar while this program is running.

=item B<-q    --quiet>

Use this flag to supress warnings if any variants in the SNP reference file do not have IDs.

=item B<-h    --help>

Display help message.

=item B<-m    --manual>

Show manual page.

=back 

=cut

=head1 DESCRIPTION

This program will annotate a VCF file with SNP IDs, frequency and dbSNP build information from a given dbSNP VCF file and can optionally filter SNPs from a vcf file on user-specified criteria. In its simplest form this program writes ID fields from files specified using the --dbsnp argument to matching variants in input files and adds any dbSNP build, frequency or clinical significance information to the INFO field. However, it can also be used to filter variants based on their presence in different builds of dbSNP or on allele frequency. Use the --build, --freq and --pathogenic arguments to set up your filtering parameters as detailed above (assuming the relevant annotations are present in the dbSNP files used). 

For example:

    annotateSnps.pl -d dbSnp138.b37.vcf.gz clinvar_20130506.vcf -b 129 -f 1 --pathogenic -i input.vcf -o input_filtered.vcf

The above command will remove variants if they were present in dbSNP build 129 or earlier dbSNP builds or if they are in later builds but have an allele frequency equal to or greater than 1 %. However, any variant with a 'pathogenic' or 'probably pathogenic' annotation will not be filtered regardless of frequency of dbSNP build.

dbSNP and ClinVar VCF files are available from NCBI's ftp site the Broad Institutes FTP site. Make sure you are using a file with the correct genome version for your data. A more consistently annotated collection of ClinVar variants is available from https://github.com/macarthur-lab/clinvar for use with the -c option instead of or as well as a ClinVar VCF.

=cut

=head1 AUTHOR

David A. Parry

=head1 COPYRIGHT AND LICENSE

Copyright 2012, 2013, 2014, 2015  David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut

