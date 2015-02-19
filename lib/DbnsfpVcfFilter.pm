package DbnsfpVcfFilter;
use strict;
use warnings;
use Carp;
use Scalar::Util qw(looks_like_number);
require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(checkAndParseDbnsfpExpressions getFrequencyFilters getAnyDamagingFilters getAllDamagingFilters filterDbnsfpExpressions);
my %dbnsfp_short = ( #selection of shortcuts for the dbNSFP fields
    "1000Gp1_AC" => "dbNSFP_1000Gp1_AC",
    "1000Gp1_AF" => "dbNSFP_1000Gp1_AF",
    "1000Gp1_AFR_AC" => "dbNSFP_1000Gp1_AFR_AC",
    "1000Gp1_AFR_AF" => "dbNSFP_1000Gp1_AFR_AF",
    "1000Gp1_AMR_AC" => "dbNSFP_1000Gp1_AMR_AC",
    "1000Gp1_AMR_AF" => "dbNSFP_1000Gp1_AMR_AF",
    "1000Gp1_ASN_AC" => "dbNSFP_1000Gp1_ASN_AC",
    "1000Gp1_ASN_AF" => "dbNSFP_1000Gp1_ASN_AF",
    "1000Gp1_EUR_AC" => "dbNSFP_1000Gp1_EUR_AC",
    "1000Gp1_EUR_AF" => "dbNSFP_1000Gp1_EUR_AF",
    "29way_logOdds" => "dbNSFP_29way_logOdds",
    "29way_pi" => "dbNSFP_29way_pi",
    "AA" => "dbNSFP_ESP6500_AA_AF",
    "AA_AF" => "dbNSFP_ESP6500_AA_AF",
    "AC" => "dbNSFP_1000Gp1_AC",
    "AF" => "dbNSFP_1000Gp1_AF",
    "AFR" => "dbNSFP_1000Gp1_AFR_AC",
    "AFR" => "dbNSFP_1000Gp1_AFR_AF",
    "AFR_AC" => "dbNSFP_1000Gp1_AFR_AC",
    "AFR_AF" => "dbNSFP_1000Gp1_AFR_AF",
    "AMR" => "dbNSFP_1000Gp1_AMR_AC",
    "AMR" => "dbNSFP_1000Gp1_AMR_AF",
    "AMR_AC" => "dbNSFP_1000Gp1_AMR_AC",
    "AMR_AF" => "dbNSFP_1000Gp1_AMR_AF",
    "ASN" => "dbNSFP_1000Gp1_ASN_AC",
    "ASN" => "dbNSFP_1000Gp1_ASN_AF",
    "ASN_AC" => "dbNSFP_1000Gp1_ASN_AC",
    "ASN_AF" => "dbNSFP_1000Gp1_ASN_AF",
    "EA" => "dbNSFP_ESP6500_EA_AF",
    "EA_AF" => "dbNSFP_ESP6500_EA_AF",
    "ESP6500_AA_AF" => "dbNSFP_ESP6500_AA_AF",
    "ESP6500_EA_AF" => "dbNSFP_ESP6500_EA_AF",
    "EUR" => "dbNSFP_1000Gp1_EUR_AC",
    "EUR" => "dbNSFP_1000Gp1_EUR_AF",
    "EUR_AC" => "dbNSFP_1000Gp1_EUR_AC",
    "EUR_AF" => "dbNSFP_1000Gp1_EUR_AF",
    "Ensembl_geneid" => "dbNSFP_Ensembl_geneid",
    "Ensembl_transcriptid" => "dbNSFP_Ensembl_transcriptid",
    "FATHMM" => "dbNSFP_FATHMM_pred",
    "FATHMM_pred" => "dbNSFP_FATHMM_pred",
    "FATHMM_score" => "dbNSFP_FATHMM_score",
    "FATHMM_score_converted" => "dbNSFP_FATHMM_score_converted",
    "GERP" => "dbNSFP_GERP++_RS",
    "GERP++" => "dbNSFP_GERP++_RS",
    "GERP++_NR" => "dbNSFP_GERP++_NR",
    "GERP++_RS" => "dbNSFP_GERP++_RS",
    "Interpro_domain" => "dbNSFP_Interpro_domain",
    "LRT" => "dbNSFP_LRT_pred",
    "LRT_Omega" => "dbNSFP_LRT_Omega",
    "LRT_pred" => "dbNSFP_LRT_pred",
    "LRT_score" => "dbNSFP_LRT_score",
    "LRT_score_converted" => "dbNSFP_LRT_score_converted",
    "MutationAssessor" => "dbNSFP_MutationAssessor_pred",
    "MutationAssessor_pred" => "dbNSFP_MutationAssessor_pred",
    "MutationAssessor_score" => "dbNSFP_MutationAssessor_score",
    "MutationAssessor_score_converted" => "dbNSFP_MutationAssessor_score_converted",
    "MutationTaster" => "dbNSFP_MutationTaster_pred",
    "MutationTaster_pred" => "dbNSFP_MutationTaster_pred",
    "MutationTaster_score" => "dbNSFP_MutationTaster_score",
    "MutationTaster_score_converted" => "dbNSFP_MutationTaster_score_converted",
    "Polyphen2" => "dbNSFP_Polyphen2_HVAR_pred",
    "Polyphen2_HDIV_pred" => "dbNSFP_Polyphen2_HDIV_pred",
    "Polyphen2_HDIV_score" => "dbNSFP_Polyphen2_HDIV_score",
    "Polyphen2_HVAR" => "dbNSFP_Polyphen2_HVAR_pred",
    "Polyphen2_HVAR_pred" => "dbNSFP_Polyphen2_HVAR_pred",
    "Polyphen2_HVAR_score" => "dbNSFP_Polyphen2_HVAR_score",
    "Polyphen2_pred" => "dbNSFP_Polyphen2_HVAR_pred",
    "Polyphen2_score" => "dbNSFP_Polyphen2_HVAR_score",
    "SIFT" => "dbNSFP_SIFT_pred",
    "SIFT_pred" => "dbNSFP_SIFT_pred",
    "SIFT_score" => "dbNSFP_SIFT_score",
    "SIFT_score_converted" => "dbNSFP_SIFT_score_converted",
    "UniSNP_ids" => "dbNSFP_UniSNP_ids",
    "Uniprot_aapos" => "dbNSFP_Uniprot_aapos",
    "Uniprot_acc" => "dbNSFP_Uniprot_acc",
    "Uniprot_id" => "dbNSFP_Uniprot_id",
    "aapos" => "dbNSFP_aapos",
    "cds_strand" => "dbNSFP_cds_strand",
    "codonpos" => "dbNSFP_codonpos",
    "gene" => "dbNSFP_Ensembl_geneid",
    "geneid" => "dbNSFP_Ensembl_geneid",
    "genename" => "dbNSFP_genename",
    "name" => "dbNSFP_genename",
    "phyloP" => "dbNSFP_phyloP",
    "refcodon" => "dbNSFP_refcodon",
    "strand" => "dbNSFP_cds_strand",
    "transcript" => "dbNSFP_Ensembl_transcriptid",
    "transcriptid" => "dbNSFP_Ensembl_transcriptid",
);

my %operators = (
    '&' => 'and',
    '&&' => 'and',
    'and' => 'and',
    '|' => 'or',
    '||' => 'or',
    'or' => 'or',
    '^' => 'xor',
    'xor' => 'xor',
);

my %num_comparators = (
    '!' => '!=',
    '!=' => '!=',
    '<' => '<',
    '<=' => '<=',
    '=' => '==',
    '==' => '==',
    '>' => '>',
    '>=' => '>=',
);
my %char_comparators = (
    '!' => 'ne',
    '!=' => 'ne',
    '<' => 'lt',
    '<=' => 'le',
    '=' => 'eq',
    '==' => 'eq',
    '>' => 'gt',
    '>=' => 'ge',
    'eq' => 'eq',
    'ge' => 'ge',
    'gt' => 'gt',
    'le' => 'le',
    'lt' => 'lt',
    'ne' => 'ne',
);

my %benign_prediction_progs = (dbNSFP_MutationAssessor_pred => ['L', 'N'],
                        dbNSFP_LRT_pred => ['N'],
                        dbNSFP_MutationTaster_pred => ['N', 'P'],
                        dbNSFP_Polyphen2_HVAR_pred => ['B'],
                        dbNSFP_SIFT_pred => ['T'],
                    );
my %path_prediction_progs = (dbNSFP_MutationAssessor_pred => ['H', 'M'],
                        dbNSFP_LRT_pred => ['D'],
                        dbNSFP_MutationTaster_pred => ['A', 'D'],
                        dbNSFP_Polyphen2_HVAR_pred => ['P', 'D'],
                        dbNSFP_SIFT_pred => ['D'],
                    );

my @frequency_fields = qw(
                dbNSFP_1000Gp1_AF
                dbNSFP_1000Gp1_AFR_AF
                dbNSFP_1000Gp1_AMR_AF
                dbNSFP_1000Gp1_ASN_AF
                dbNSFP_1000Gp1_EUR_AF
                dbNSFP_ESP6500_AA_AF
                dbNSFP_ESP6500_EA_AF
                );


##################################################
#returns array of filters like checkAndParseDbnsfpExpressions
#containing all population frequency fields
sub getFrequencyFilters{
    #$comparator tells whether we want to filter variants that are greater/less/equal etc. 
    #to $freq.  $operator defines whether we want any (or) or all (and) or just one (xor)
    #frequencies to match expression in order to filter
    my ($freq, $comparator, $operator) = @_;
    my @filters;
    #check arguments
    croak("getFrequencyFilters method requires frequency to filter on")
         if not defined $freq;
    croak("getFrequencyFilters method frequency must be between 0.00 and 1.00") 
        if $freq > 1.00 or $freq < 0.00;
    $comparator = '>=' if not $comparator;
    croak("comparator \"$comparator\" passed to getFrequencyFilters method is "
        . "not recognised") if not exists $num_comparators{$comparator};
    $operator = 'or' if not $operator;
    croak("Operator \"$operator\" passed to getFrequencyFilters subroutine not recognised "
        ."- please use either 'and', 'or' or 'xor'.") if not exists $operators{$operator};
    #finished checking arguments
    #create filter expressions
    my $exp_hash = {field => [], operator => [], expression => []};
    for (my $i = 0; $i < @frequency_fields; $i++){
        push (@{$exp_hash->{field}}, $frequency_fields[$i]);
        push (@{$exp_hash->{expression}}, "$num_comparators{$comparator}$freq");
        push (@{$exp_hash->{operator}}, $operators{$operator}) if $i < (@frequency_fields -1);
    }
    return ($exp_hash) if wantarray;
    return $exp_hash if defined wantarray;
    carp ("getAllDamagingFilters called in void context.");
}

##################################################
sub getAnyDamagingFilters{
    #return anon hash like those in checkAndParseDbnsfpExpressions
    #containing all standard damaging predictions by default
    #returns anon hash if called in scalar context, array containing a single anon hash
    #if called in array context
    my ($operator, $progs) = @_;
    #by default all programs must give benign prediction in order to filter
    # so $operator is 'and' by default
    my $exp_hash = {field => [], operator => [], expression => []};
    $operator = 'and' if not $operator;
    croak("Operator \"$operator\" passed to getAnyDamagingFilters subroutine not recognised "
        ."- please use either 'and', 'or' or 'xor'.") if not exists $operators{$operator};
    croak("Second argument to getAnyDamagingFilters subroutine must be an ARRAY "
        ."reference to a list of program names or UNDEF.") if defined $progs and ref $progs ne 'ARRAY';

    my @filters = ();
    my @pr;
    if (not defined $progs){
        #return all prediction programs
        push @$progs, keys %benign_prediction_progs;
    }
    for (my $i = 0; $i < @$progs;  $i++){
        if (not exists $benign_prediction_progs{$progs->[$i]}){
            if (not exists $dbnsfp_short{$progs->[$i]}){
                croak("Unrecognised prediction program $progs->[$i] passed to getAnyDamagingFilters method.");
            }elsif(not exists $benign_prediction_progs{$dbnsfp_short{$progs->[$i]}}){
                croak("Unrecognised prediction program $progs->[$i] passed to getAnyDamagingFilters method.");
            }else{
                $progs->[$i] = $dbnsfp_short{$progs->[$i]};
            }
        }
        for (my $j = 0; $j < @{$benign_prediction_progs{$progs->[$i]}}; $j++){
            if ($j == 0){
                #enclose each prediction program's values in brackets in case of multiple values
                push (@{$exp_hash->{field}}, "($progs->[$i]");
            }else{
                push (@{$exp_hash->{field}}, $progs->[$i]);
            }
            if ($j == (@{$benign_prediction_progs{$progs->[$i]}} - 1)){
                #enclose each prediction program's values in brackets in case of multiple values
                push (@{$exp_hash->{expression}}, "eq \"$benign_prediction_progs{$progs->[$i]}->[$j]\")");
            }else{
                push (@{$exp_hash->{expression}}, "eq \"$benign_prediction_progs{$progs->[$i]}->[$j]\"");
            }
            if ($j < (@{$benign_prediction_progs{$progs->[$i]}} - 1)){
                #multiple benign values given for a prediction program 
                #(e.g. 'L' and 'N' for dbNSFP_MutationAssessor_pred)
                #should be separated with 'or' regardless of value of $operator
                push (@{$exp_hash->{operator}}, 'or') 
            }
        }
        if ($i < @$progs -1){
            push (@{$exp_hash->{operator}}, $operator) 
        }
    }
    return ($exp_hash) if wantarray;
    return $exp_hash if defined wantarray;
    carp ("getAnyDamagingFilters called in void context.");
}



##################################################
sub getAllDamagingFilters{
    #return anon hash like those in checkAndParseDbnsfpExpressions
    #containing all standard damaging predictions by default
    #returns anon hash if called in scalar context, array containing a single anon hash
    #if called in array context
    my ($operator, $progs) = @_;
    #by default all programs must give damaging prediction or else we filter
    # so $operator is 'and' by default
    my $exp_hash = {field => [], operator => [], expression => []};
    $operator = 'and' if not $operator;
    croak("Operator \"$operator\" passed to getAllDamagingFilters subroutine not recognised "
        ."- please use either 'and', 'or' or 'xor'.") if not exists $operators{$operator};
    croak("Second argument to getAllDamagingFilters subroutine must be an ARRAY "
        ."reference to a list of program names or UNDEF.") if defined $progs and ref $progs ne 'ARRAY';

    my @filters = ();
    my @pr;
    if (not defined $progs){
        #return all prediction programs
        push @$progs, keys %path_prediction_progs;
    }
    for (my $i = 0; $i < @$progs;  $i++){
        if (not exists $path_prediction_progs{$progs->[$i]}){
            if (not exists $dbnsfp_short{$progs->[$i]}){
                croak("Unrecognised prediction program $progs->[$i] passed to getAllDamagingFilters method.");
            }elsif(not exists $path_prediction_progs{$dbnsfp_short{$progs->[$i]}}){
                croak("Unrecognised prediction program $progs->[$i] passed to getAllDamagingFilters method.");
            }else{
                $progs->[$i] = $dbnsfp_short{$progs->[$i]};
            }
        }
        for (my $j = 0; $j < @{$path_prediction_progs{$progs->[$i]}}; $j++){
            if ($j == 0){
                #enclose each prediction program's values in brackets in case of multiple values
                push (@{$exp_hash->{field}}, "($progs->[$i]");
            }else{
                push (@{$exp_hash->{field}}, $progs->[$i]);
            }
            if ($j == (@{$path_prediction_progs{$progs->[$i]}} - 1)){
                #enclose each prediction program's values in brackets in case of multiple values
                push (@{$exp_hash->{expression}}, "ne \"$path_prediction_progs{$progs->[$i]}->[$j]\")");
            }else{
                push (@{$exp_hash->{expression}}, "ne \"$path_prediction_progs{$progs->[$i]}->[$j]\"");
            }
            if ($j < (@{$path_prediction_progs{$progs->[$i]}} - 1)){
                #multiple benign values given for a prediction program 
                #(e.g. 'M' and 'H' for dbNSFP_MutationAssessor_pred)
                #should be separated with 'and' regardless of value of $operator
                push (@{$exp_hash->{operator}}, 'and') 
            }
        }
        if ($i < @$progs -1){
            push (@{$exp_hash->{operator}}, $operator) 
        }
    }
    return ($exp_hash) if wantarray;
    return $exp_hash if defined wantarray;
    carp ("getAllDamagingFilters called in void context.");
}

##################################################
sub filterDbnsfpExpressions{
    #return 1 if line should be filtered
    my ($expressions, $vcf_line, $multi_value_operator) = @_;
    #$expressions is reference to an array of filter expressions
    #as returned by checkAndParseDbnsfpExpressions subroutine
    #$multi_value_operator variable tells how to handle fields with multiple values 
    # - default is to only filter if all match condition
    $multi_value_operator = 'and' if not $multi_value_operator;
    croak("Operator \"$multi_value_operator\" passed to filterDbnsfpExpressions subroutine not recognised "
        ."- please use either 'and', 'or' or 'xor'.") if not exists $operators{$multi_value_operator};
    croak("First value passed to filterDbnsfpExpressions must be a "
        ."reference to an array") if ref $expressions ne 'ARRAY';
    my @line = split("\t", $vcf_line);
    my @info = split(";", $line[7]);
    my %values = ();
    foreach my $inf (@info){
        my ($field, $value) = split("=", $inf);
        $values{$field} = $value if defined $value;
    }
    my $n = 0;
EXPR: foreach my $e (@$expressions){
        $n++;
        #sanity checks first
        croak("Array reference passed to filterDbnsfpExpressions must contain "
            ."HASH references only.") if ref $e ne 'HASH';
        foreach (qw(field expression operator)){
            if (not exists ($e->{$_})){
                croak("Missing required hash entry: $_ in expression $n passed to "
                    ."filterDbnsfpExpressions subroutine.");
            }
            if (ref $e->{$_} ne 'ARRAY'){
                my $ref = ref $e->{$_};
                if ($ref){
                    croak("Hash entry $_ in expression $n passed to filterDbnsfpExpressions "
                        ."subroutine must be an ARRAY reference but $ref reference found instead.");
                }else{
                    croak("Hash entry $_ in expression $n passed to filterDbnsfpExpressions "
                        ."subroutine must be an ARRAY reference but instead has value $e->{$_}.");
                }
            }
        }
        #checks done - create and test @eval expression
        my @eval_full = ();
        for (my $i = 0; $i < @{$e->{field}}; $i++){
            my @eval = ();
            my $f = $e->{field}->[$i];
            $f =~ s/^\(//;#we may have used ( to specify precedence
            my $x = $e->{expression}->[$i];
            $x =~ s/\)$//;#we may have used ) to specify precedence
            if (not exists $values{$f}){#not in INFO field, can't evaluate expression
                #make this eval 0 if empty so that and/or/xor operator functions
                #for other parts of this expression are still honoured
                push @eval_full, 0; 
                push @eval_full, $e->{operator}->[$i] if $i < @{$e->{operator}};
                next;
            }
            #multiple entries may be separated by '|' or ','
            my @vals = split(/[\|\,]/, $values{$f});
            foreach (my $j = 0; $j < @vals -1; $j++){
                next if ($vals[$j] eq '.');
                
                #only add multi value operator if we have a previous value
                push @eval, $operators{$multi_value_operator} if $j > 0 and @eval > 0;
                push @eval, "\"$vals[$j]\"";
                push @eval, "$x";
            }
            unless ($vals[-1] eq '.' ){
                #only add multi value operator if we have a previous value
                push @eval, $operators{$multi_value_operator} if @eval > 0;
                push @eval, "\"$vals[-1]\"";
                push @eval, "$x";
            }
            #if @eval is empty then all our values are '.' (i.e. not defined)
            if (@eval == 0){
                #make this eval 0 if empty so that and/or/xor operator functions
                #for other parts of this expression are still honoured
                push @eval_full, 0; 
                push @eval_full, $e->{operator}->[$i] if $i < @{$e->{operator}};
                next;
            }
            if ($e->{field}->[$i] =~/^\(/){
                unshift @eval, "(";#we may have used ( to specify precedence
            }
            if ($e->{expression}->[$i] =~/\)$/){
                push @eval, ")";#we may have used ( to specify precedence
            }
            push @eval, $e->{operator}->[$i] if $i < @{$e->{operator}};
            push @eval_full, @eval;
        }
        #DEBUG
      #  my $res = eval join(" ", @eval_full);
      #  my $error = $@;
      #  print "DEBUG: " . join(" ", @eval_full) . " : $res : ";
      #  print "\n";
      #  print $error if $error;
        #end DEBUG
        my $ev =eval join(" ", @eval_full);
        carp "$@\n" if $@;
        return 1 if $ev;
    }
    return 0;
}

##################################################
sub checkAndParseDbnsfpExpressions{
    my ($expressions, $head) = @_;
    #$expressions should be an array of strings in format 
    # "1000Gp1_EUR_AF<=0.01" 
    # or multiple of these expressions joined by and/or/xor 
    # e.g. "dbNSFP_Polyphen2_HVAR_pred=B  or SIFT_pred=T"
    #$head is a vcf header containing at a minimum the relevant dbNSFP INFO lines
    #returns
    my @filters = (); #array of anon hashes
                      # each anon hash (
                      # field => [dbNSFP_fieldX, dbNSFP_fieldX, dbNSFP_fieldY]
                      # operator => [ or, and] - must equal size of fields -1
                      # expression => [>1, <5, >1] - must equal size of fields
    my %types = (); #key is field, value is types
    foreach my $info ( split("\n", $head) ){
        if ($info =~ /##INFO=<ID=(dbNSFP_\S+),Number=\S+,Type=(\w+),Description=.*>$/){
            $types{$1} = $2;
        }
    }
	
    foreach my $exp (@$expressions){
        my $open_brackets = () = $exp =~ /\(/g;
        my $close_brackets = () = $exp =~ /\)/g;
        if ($open_brackets > $close_brackets){
            croak "Unclosed brackets in expression $exp\n";
        }elsif ($open_brackets < $close_brackets){
            croak "Trailing brackets in expression $exp\n";
        }
        $exp =~ s/^\s+//;
        my $exp_hash = {field => [], operator => [], expression => []};
        #split on white space. 
        my @split = split(/\s/, $exp); 
        #every second array element must be an operator (and/or/xor)
        for(my $i = 0; $i<@split; $i++){
        #so if index % 2 == 0 check expression
            if ($i %2 ==0){
                if ($split[$i] =~ /([\(\w\+]+)([=<>!]{1,2})(\S+)/){
                    my ($field, $comparator, $value) = ($1, $2, $3);
                    (my $f = $field) =~ s/^\(//;#we may have used ( to specify precedence
                    (my $v = $value) =~ s/\)$//;#we may have used ( to specify precedence
                    #Check $field exists in header 
                    my $t;
                    if (exists $types{$f}){
                       $t = $types{$f};
                    }elsif(exists $dbnsfp_short{$f}){
                        if (exists $types{$dbnsfp_short{$f}}){
                            $f = $dbnsfp_short{$f};
                            $t = $types{$f};
                        }else{
                            croak("Can't find $dbnsfp_short{$f} ($f) dbNSFP INFO field "
                                ."in header - please ensure you have annotated your VCF with "
                                ."appropriate dbNSFP fields using SnpSift.jar\n");
                        }
                    }else{
                        croak("Can't find $f dbNSFP INFO field in header"
                                ." - please ensure you annotate your VCF with "
                                ."appropriate dbNSFP fields using SnpSift.jar\n");
                    }
                    #Check comparator is valid and value matches $field type
                    if ($t eq 'Integer' or $t eq 'Float'){
                        if (exists $num_comparators{$comparator}){
                            $comparator = $num_comparators{$comparator};
                        }else{
                            croak("Unrecognised comparator ($comparator) for "
                                ."dbNSFP INFO field type $t\n");
                        }
                        if (not looks_like_number($v)){
                            croak("dbNSFP INFO field for $field is of type $t but value "
                                ."($v) in expression ($split[$i]) does not "
                                ."look like a number.\n");
                        }
                    }elsif($t eq 'Character' or $t eq 'String'){
                        if (exists $char_comparators{$comparator}){
                            if (looks_like_number($v)){
                                print STDERR "Warning: dbNSFP INFO field for $field is of "
                                    ."type $t but value ($v) in expression ($split[$i])"
                                    ." look like a number.\n";
                                if (exists $num_comparators{$comparator}){
                                    print STDERR "Interpretting as a number.\n";
                                    print STDERR "Changing comparator ($comparator) to "
                                        ."$num_comparators{$comparator} in expression ($split[$i]).\n"
                                        if $num_comparators{$comparator} ne $comparator;
                                    $comparator = $num_comparators{$comparator};
                                }else{
                                    print STDERR "Interpretting as a character/string.\n";
                                    $comparator = $char_comparators{$comparator};
                                }
                            }else{
                                $comparator = $char_comparators{$comparator};
                            }
                        }else{
                            croak("Unrecognised comparator ($comparator) for "
                                ."dbNSFP INFO field type $t\n");
                        }
                    }else{
                        if ($t eq 'Flag'){
                            croak("dbNSFP Flag fields not currently supported.\n"
                                ."Offending expression: $split[$i]\n");
                        }
                        croak("Unrecognised INFO type ($t) in dbNSFP field.\n"
                            ."Offending expression: $split[$i]\n");
                    }
                    #Create a valid perl expression from user input
                    # $exp_hash = {field => [], operator => [], expression => []};
                    if ($field =~ /^\(/){
                        push (@{$exp_hash->{field}}, "($f");
                    }else{
                        push (@{$exp_hash->{field}}, $f);
                    }
                    if ($value =~ s/\)$//){
                        push (@{$exp_hash->{expression}}, "$comparator \"$v\")");
                    }else{
                        push (@{$exp_hash->{expression}}, "$comparator \"$v\"");
                    }
                }else{
                    croak("Invalid dbNSFP expression: $exp\n");
                }
            #otherwise check operator
            }else{
                if (exists $operators{lc$split[$i]}){
                    push (@{$exp_hash->{operator}}, $operators{lc$split[$i]});
                }else{
                    my $pos = $i +1 ;
                    croak("Expected operator (or/and/xor) at position $pos of dbNSFP expression "
                        ."$exp but found $split[$i].\n");
                }
            }
        }
        #sanity check for $exp_hash
        # field => [dbNSFP_fieldX, dbNSFP_fieldX, dbNSFP_fieldY]
        # expression => [>1, <5, >1] - must equal size of fields
        # operator => [ or, and] - must equal size of fields -1
        if (@{$exp_hash->{field}} != @{$exp_hash->{expression}}){
            croak("Number of fields does not match number of expressions "
                . "for expression $exp\n"); #shouldn't be possible
        }
        if (@{$exp_hash->{operator}} != (@{$exp_hash->{field}} -1) ){
            if (@{$exp_hash->{operator}} > (@{$exp_hash->{field}} -1) ){
                croak("Trailing operators in dbNSFP expression $exp\n");
            }else{
                my $f = @{$exp_hash->{field}};
                my $o = @{$exp_hash->{operator}};
                croak("Invalid number of fields vs. operators in dbNSFP expression $exp  - found $f fields "
                    ." so expected " . ($f -1 ) . " operators but found $o.\n"); #this shouldn't be possible
            }
        }
        push @filters, $exp_hash;
        
    }
    return @filters if defined wantarray;
    carp ("checkAndParseDbnsfpExpressions called in void context.");
}

####################################################################
sub getDbnsfpFieldsFromHeader{
    my ($head) = @_;
    my @fields = ();
    #$head is a vcf header containing at a minimum the relevant dbNSFP INFO lines
    foreach my $info ( split("\n", $head) ){
        if ($info =~ /##INFO=<ID=(dbNSFP_\S+),Number=\S+,Type=(\w+),Description=.*>$/){
            push @fields, $1;
        }
    }
    return @fields if defined wantarray;
    carp ("checkAndParseDbnsfpExpressions called in void context.");
}

####################################################################
sub getDbnsfpValue{
    my ($vcf_line, $f) = @_;
    #$f is the name of the dbnsfp field we want the value for
    #$vcf_line is a complete VCF line containing dbNSFP annotations
    #returns the value for given field if found, otherwise returns nothing
    my @line = split("\t", $vcf_line);
    my @info = split(";", $line[7]);
    my %values = ();
    foreach my $inf (@info){
        my ($field, $value) = split("=", $inf);
        $values{$field} = $value if defined $value;
    }
    return $values{$f} if (exists $values{$f});
    return;
}

####################################################################
sub getDbnsfpAnnotations{
    my ($vcf_line) = @_;
    #$vcf_line is a complete VCF line containing dbNSFP annotations
    #returns hash of field names to values
    my @line = split("\t", $vcf_line);
    my @info = split(";", $line[7]);
    my %values = ();
    foreach my $inf (@info){
        my ($field, $value) = split("=", $inf);
        $values{$field} = $value if defined $value;
    }
    return %values;
}

1;


=head1 NAME

DbnsfpVcfFilter.pm - read and filter dbNSFP fields added by SnpSift.jar

=head1 SYNOPSIS

use DbnsfpVcfFilter; 

my @expressions = ("1000Gp1_EUR_AFE<lt>=0.01", "Polyphen2=B or SIFT=T");

my @filters = checkAndParseDbnsfpExpressions(\@expressions, $vcf_header_string);

filterDbnsfpExpressions(\@filters, $vcf_line);

=head1 DESCRIPTION

Check/create filter expressions to test against VCF lines containing dbNSFP annotations added using SnpSift.jar (http://snpeff.sourceforge.net/SnpSift.html). Also contains a few simple functions for returning dbNSFP annotations.


=head1 FUNCTIONS

=head2 checkAndParseDbnsfpExpressions

The purpose of this function is to convert strings expressing a filtering condition according to dbNSFP annotations into formats accepted by this module's filterDbnsfpExpressions function and to exit if the strings are not valid.  It returns an array of anonymous hashes that can be passed to filterDbnsfpExpressions to perform the desired filtering.

The first argument passed to this function must be a reference to an array of strings describing filter expressions in the format described below.  The second argument must be the header string from your VCF.

Each expression consists of a dbNSFP annotation (e.g. "dbNSFP_GERP++_RS"), a comparator (e.g. "<") and a value (e.g. "2"). Multiple expressions can be joined with 'and', 'or' or 'xor' operators (given as strings) if desired in which case these operators behave as they would in standard perl code. Valid comparators are:

=over 8

= (or ==): equal to

! (or !=): not equal to

<: less than

<= less than or equal to

>: more than

>=: more than or equal to

=back 

Refer to dbNSFP documentation (https://sites.google.com/site/jpopgen/dbNSFP) for details of the available dbNSFP fields. 

B<Examples of filter expressions:>

=over

"dbNSFP_GERP++_RSE<lt>2"

"dbNSFP_SIFT_pred=T or dbNSFP_Polyphen2_HVAR_pred=B"

"dbNSFP_ESP6500_EA_AFE<gt>=0.01 and dbNSFP_ESP6500_AA_AFE<lt>=0.02"

=back

In the first example any line with a GERP Rejected Substitution score less than 2 will return 1 when processed by filterDbnsfpExpressions. In the second example any variant that SIFT predicts as 'tolerated' (i.e. the dbNSFP_SIFT_pred annotation is 'T') or predicted by Polyphen to be benign (i.e. the dbNSFP_Polyphen2_HVAR_pred annotation is 'B') will be filtered. In the third example a variant is filtered if the NHLBI ESP6500 European American allele frequency is greater than or equal to 0.01 (1 %) and the NHLBI ESP6500 African American allele frequency if less than or equal to 0.02 (2 %). 

This module accepts the following shortened names for fields:

=over 8

    1000Gp1_AC: dbNSFP_1000Gp1_AC
    1000Gp1_AF: dbNSFP_1000Gp1_AF
    1000Gp1_AFR_AC: dbNSFP_1000Gp1_AFR_AC
    1000Gp1_AFR_AF: dbNSFP_1000Gp1_AFR_AF
    1000Gp1_AMR_AC: dbNSFP_1000Gp1_AMR_AC
    1000Gp1_AMR_AF: dbNSFP_1000Gp1_AMR_AF
    1000Gp1_ASN_AC: dbNSFP_1000Gp1_ASN_AC
    1000Gp1_ASN_AF: dbNSFP_1000Gp1_ASN_AF
    1000Gp1_EUR_AC: dbNSFP_1000Gp1_EUR_AC
    1000Gp1_EUR_AF: dbNSFP_1000Gp1_EUR_AF
    29way_logOdds: dbNSFP_29way_logOdds
    29way_pi: dbNSFP_29way_pi
    AA: dbNSFP_ESP6500_AA_AF
    AA_AF: dbNSFP_ESP6500_AA_AF
    AC: dbNSFP_1000Gp1_AC
    AF: dbNSFP_1000Gp1_AF
    AFR: dbNSFP_1000Gp1_AFR_AC
    AFR: dbNSFP_1000Gp1_AFR_AF
    AFR_AC: dbNSFP_1000Gp1_AFR_AC
    AFR_AF: dbNSFP_1000Gp1_AFR_AF
    AMR: dbNSFP_1000Gp1_AMR_AC
    AMR: dbNSFP_1000Gp1_AMR_AF
    AMR_AC: dbNSFP_1000Gp1_AMR_AC
    AMR_AF: dbNSFP_1000Gp1_AMR_AF
    ASN: dbNSFP_1000Gp1_ASN_AC
    ASN: dbNSFP_1000Gp1_ASN_AF
    ASN_AC: dbNSFP_1000Gp1_ASN_AC
    ASN_AF: dbNSFP_1000Gp1_ASN_AF
    EA: dbNSFP_ESP6500_EA_AF
    EA_AF: dbNSFP_ESP6500_EA_AF
    ESP6500_AA_AF: dbNSFP_ESP6500_AA_AF
    ESP6500_EA_AF: dbNSFP_ESP6500_EA_AF
    EUR: dbNSFP_1000Gp1_EUR_AC
    EUR: dbNSFP_1000Gp1_EUR_AF
    EUR_AC: dbNSFP_1000Gp1_EUR_AC
    EUR_AF: dbNSFP_1000Gp1_EUR_AF
    Ensembl_geneid: dbNSFP_Ensembl_geneid
    Ensembl_transcriptid: dbNSFP_Ensembl_transcriptid
    FATHMM: dbNSFP_FATHMM_pred
    FATHMM_pred: dbNSFP_FATHMM_pred
    FATHMM_score: dbNSFP_FATHMM_score
    FATHMM_score_converted: dbNSFP_FATHMM_score_converted
    GERP: dbNSFP_GERP++_RS
    GERP++: dbNSFP_GERP++_RS
    GERP++_NR: dbNSFP_GERP++_NR
    GERP++_RS: dbNSFP_GERP++_RS
    Interpro_domain: dbNSFP_Interpro_domain
    LRT: dbNSFP_LRT_pred
    LRT_Omega: dbNSFP_LRT_Omega
    LRT_pred: dbNSFP_LRT_pred
    LRT_score: dbNSFP_LRT_score
    LRT_score_converted: dbNSFP_LRT_score_converted
    MutationAssessor: dbNSFP_MutationAssessor_pred
    MutationAssessor_pred: dbNSFP_MutationAssessor_pred
    MutationAssessor_score: dbNSFP_MutationAssessor_score
    MutationAssessor_score_converted: dbNSFP_MutationAssessor_score_converted
    MutationTaster: dbNSFP_MutationTaster_pred
    MutationTaster_pred: dbNSFP_MutationTaster_pred
    MutationTaster_score: dbNSFP_MutationTaster_score
    MutationTaster_score_converted: dbNSFP_MutationTaster_score_converted
    Polyphen2: dbNSFP_Polyphen2_HVAR_pred
    Polyphen2_HDIV_pred: dbNSFP_Polyphen2_HDIV_pred
    Polyphen2_HDIV_score: dbNSFP_Polyphen2_HDIV_score
    Polyphen2_HVAR: dbNSFP_Polyphen2_HVAR_pred
    Polyphen2_HVAR_pred: dbNSFP_Polyphen2_HVAR_pred
    Polyphen2_HVAR_score: dbNSFP_Polyphen2_HVAR_score
    Polyphen2_pred: dbNSFP_Polyphen2_HVAR_pred
    Polyphen2_score: dbNSFP_Polyphen2_HVAR_score
    SIFT: dbNSFP_SIFT_pred
    SIFT_pred: dbNSFP_SIFT_pred
    SIFT_score: dbNSFP_SIFT_score
    SIFT_score_converted: dbNSFP_SIFT_score_converted
    UniSNP_ids: dbNSFP_UniSNP_ids
    Uniprot_aapos: dbNSFP_Uniprot_aapos
    Uniprot_acc: dbNSFP_Uniprot_acc
    Uniprot_id: dbNSFP_Uniprot_id
    aapos: dbNSFP_aapos
    cds_strand: dbNSFP_cds_strand
    codonpos: dbNSFP_codonpos
    gene: dbNSFP_Ensembl_geneid
    geneid: dbNSFP_Ensembl_geneid
    genename: dbNSFP_genename
    name: dbNSFP_genename
    phyloP: dbNSFP_phyloP
    refcodon: dbNSFP_refcodon
    strand: dbNSFP_cds_strand
    transcript: dbNSFP_Ensembl_transcriptid
    transcriptid: dbNSFP_Ensembl_transcriptid

=back

=head2 filterDbnsfpExpressions

Takes a reference to an array of filter conditions (in the form of anonymous hashes described below) as its first argument and a vcf line as its second argument. Optionally takes a third argument which must be an 'and', 'or' or 'xor' operator to tell it how to handle annotations with multiple values (e.g. Polyphen2 will often have multiple predictions for each protein isoform encoded). By default this value is 'and' such that all values for an annotation must match the condition. Returns 1 if the line meets a filter condition and 0 if not.

The anonymous hashes contained in the array of filter conditions should follow the following structure where each key points to an anonymous array:

$exp_hash = {field => [], expression => [], operator => []};

Expressions are constructed by finding the value of $exp_hash->{field}->[$i] in the vcf line provided and testing it against $exp_hash->{expression}->[$i]. There must be an equal number of items in @{$exp_hash->{field}} as in @{$exp_hash->{expression}}.  If there is more than one value for these arrays each part of this expression must be joined so that the condition tested by $exp_hash->{field}->[$i] and $exp_hash->{expression}->[$i] are followed by $exp_hash->{operator}->[$i] while $i < @{$exp_hash->{field}} -1. 

However, the best way to construct these data structures is using the checkAndParseDbnsfpExpressions function of this module. In this manner "dbNSFP_GERP++_RSE<lt>2" passed as an expression to checkAndParseDbnsfpExpressions becomes: 

=over 4

B<{field =E<gt> ["dbNSFP_GERP++_RSE"], expression =E<gt> ["<2"], operator =E<gt> []}>

=back 

"dbNSFP_SIFT_pred=T or dbNSFP_Polyphen2_HVAR_pred=B" becomes:

=over 4

B<{field =E<gt> ["dbNSFP_SIFT_pred", dbNSFP_Polyphen2_HVAR_pred,], expression =E<gt> ["eqT", "eqB"], operator =E<gt> ["or"]}>

=back 

"dbNSFP_ESP6500_EA_AFE<gt>=0.01 and dbNSFP_ESP6500_AA_AFE<lt>=0.02" becomes:

=over 4

B<{field =E<gt> ["dbNSFP_ESP6500_EA_AF", "dbNSFP_ESP6500_AA_AF",], expression =E<gt> ["E<gt>=0.01", "<=0.02"], operator =E<gt> ["and"]}>

=back


=head2 getFrequencyFilters

Constructs a set of filter expressions used to compare against allele frequencies for the following allele frequency annotations:

=over 8

    dbNSFP_1000Gp1_AF
    dbNSFP_1000Gp1_AFR_AF
    dbNSFP_1000Gp1_AMR_AF
    dbNSFP_1000Gp1_ASN_AF
    dbNSFP_1000Gp1_EUR_AF
    dbNSFP_ESP6500_AA_AF
    dbNSFP_ESP6500_EA_AF

=back

The first value passed to this method must be the frequency value to filter against (e.g. 0.01 for 1 % allele frequency value), the second (optional - defaults to ">=") argument must be a numerical comparator ("<", "<=", ">", ">=", "==" or "!="), and the final argument is an optional operator which can be either 'and', 'or' or 'xor' to determine whether all, any, or just one of the allele frequency values must match the given condition.

getFrequencyFilters(0.05, "<=", "or") is equivelent to constructing an expression from the following string with checkAndParseDbnsfpExpressions:

=over 4

"dbNSFP_1000Gp1_AFE<gt>0.05 or dbNSFP_1000Gp1_AFR_AFE<gt>0.05 or dbNSFP_1000Gp1_AMR_AFE<gt>0.05 or dbNSFP_1000Gp1_ASN_AFE<gt>0.05 or dbNSFP_1000Gp1_EUR_AFE<gt>0.05 or dbNSFP_ESP6500_AA_AFE<gt>0.05 or dbNSFP_ESP6500_EA_AFE<gt>0.05"

=back

Returns a single element array containing the appropriate anonymous hash for use with filterDbnsfpExpressions if called in an array context, otherwise returns the anonymous hash.

=head2 getAnyDamagingFilters

Constructs a set of filter expressions used to compare against prediction values given by missense prediction programs in order to filter variants with benign predictions. The first argument should be an operator ('and', 'or', or 'xor') given as a string to tell the function how to compare multiple programs. This defaults to 'and' such that ALL given programs must give a benign prediction in order to filter a line. The second argument is an array giving names of the prediction programs to be used. These can be:

=over 8

"dbNSFP_LRT_pred"

"dbNSFP_MutationAssessor_pred"

"dbNSFP_MutationTaster_pred"

"dbNSFP_Polyphen2_HVAR_pred"

"dbNSFP_SIFT_pred"

=back

If this argument is not given all of these programs will be used. 

As such, getAnyDamagingFilters() is equivelent to constructing an expression from the following string with checkAndParseDbnsfpExpressions:

=over 4

"dbNSFP_LRT_pred=N and (dbNSFP_MutationAssessor_pred=L or dbNSFP_MutationAssessor_pred=N) and (dbNSFP_MutationTaster_pred=N or dbNSFP_MutationTaster_pred=P) and dbNSFP_Polyphen2_HVAR_pred=B and dbNSFP_SIFT_pred=T"

=back 

Returns a single element array containing the appropriate anonymous hash for use with filterDbnsfpExpressions if called in an array context, otherwise returns the anonymous hash.

=head2 getAllDamagingFilters

Constructs a set of filter expressions used to compare against prediction values given by missense prediction programs in order to filter variants without damaging predictions. The first argument should be an operator ('and', 'or', or 'xor') given as a string to tell the function how to compare multiple programs. This defaults to 'and' such that ALL given programs must NOT give a damaging prediction in order to filter a line. The second argument is an array giving names of the prediction programs to be used as for getAnyDamagingFilters. If this argument is not given all of these programs will be used.  As such the following:

getAllDamagingFilters() is equivelent to constructing an expression from the following string with checkAndParseDbnsfpExpressions:

=over 4

"dbNSFP_LRT_pred!=D and (dbNSFP_MutationAssessor_pred!=H and dbNSFP_MutationAssessor_pred=!M) and (dbNSFP_MutationTaster_pred!=A and dbNSFP_MutationTaster_pred!=D) and ((dbNSFP_Polyphen2_HVAR_pred!=P and dbNSFP_Polyphen2_HVAR_pred!=D) and dbNSFP_SIFT_pred!=D"

=back 

Returns a single element array containing the appropriate anonymous hash for use with filterDbnsfpExpressions if called in an array context, otherwise returns the anonymous hash.


=head2 getDbnsfpFieldsFromHeader

Returns an array of the names of all dbNSFP fields found in the dbNSFP info field in VCF header. Takes a VCF header string as an argument - e.g.

=over 4

my @dbnsfp_fields = DbnsfpVcfFilter::getDbnsfpFieldsFromHeader($header);

=back 

=head2 getDbnsfpValue

Takes a VCF line string and the name of a dbNSFP field and returns the value for that field if found.

=over 4

my $ma_pred = DbnsfpVcfFilter::getDbnsfpValue($line, "dbNSFP_MutationAssessor_pred");

=back 

=head2 getDbnsfpAnnotations

Given an annotated VCF line, returns hash of field names to values.

=over 4

my %dbnsfp = DbnsfpVcfFilter::getDbnsfpAnnotations($line);

=back


=head1 EXPORTS

checkAndParseDbnsfpExpressions
filterDbnsfpExpressions
getFrequencyFilters
getAllDamagingFilters
getAnyDamagingFilters

=head1 AUTHOR

David A. Parry
University of Leeds


=head1 COPYRIGHT AND LICENSE

Copyright 2014 David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut


