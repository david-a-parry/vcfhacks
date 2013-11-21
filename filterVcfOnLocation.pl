#!/usr/bin/perl
#David Parry University of Leeds April 2011

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.If not, see <http://www.gnu.org/licenses/>.


=head1 NAME

filterVcfOnLocation.pl - print variants from a vcf file that lie within a list of genomic regions


=head1 SYNOPSIS
 
filterVcfOnLocation.pl -i [vars.vcf] -b [regions.bed] -o [output.vcf]

=cut

=head1 ARGUMENTS
   

=over 8

=item B<-i    --input>

VCF format file containing list of variants to filter (required). 

=item B<-b    --bed>

Bed file containing regions to filter on (required unless using --regions argument). Can also be a file containing a list of intervals, one interval per line.

=item B<-r    --regions>

Specify one or more regions to filter on seperated with spaces in format "chr1:1-2000". Can be used in conjunction with or instead of --bed argument.

=item B<-f    --flanks>

Add this number of bases to the 5' and 3' flanks of regions.

=item B<-e    --exclude>

Reverses the purpose of this script - i.e. prints variants that do NOT lie within bed regions.

=item B<-o    --output>

File to print output (optional). Will print to STDOUT by default.

=item B<-p    --progress>

If redirecting STDOUT or using B<--output> argument this flag can be used to print percent progress to STDERR. 

=item B<-c    --offset>

Column offset for chr field of vcf (e.g. if annotated vcf).

=item B<-h    --help>

Show this script's help information.

=item B<-m    --manual>

Show this script's manual page.

=back

=cut


=head1 EXAMPLES

filterVcfOnLocation.pl -i [vars.vcf] -b [regions.bed] 

filterVcfOnLocation.pl -i [vars.vcf] -r chr1:2000000-50000000 -o [filtered.vcf]

   
=cut


=head1 DESCRIPTION

Takes a bed file or a list of regions and filters variants from a vcf file that lie in these regions.


=cut

=head1 AUTHOR

David A. Parry

University of Leeds


=head1 COPYRIGHT AND LICENSE

Copyright 2011, 2012, 2013  David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut

use warnings;
use strict;
use Getopt::Long qw(:config no_ignore_case);
use Data::Dumper;
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin";
use POSIX qw/strftime/;
use SortGenomicCoordinates;
use ParseVCF;
use Term::ProgressBar;

my $help;
my $manual;
my $vcf_file;
my $outfile;
my @bedfile; 
my @reg;
my $offset;
my $exclude;
my $progress;
my $total_lines = 0;
my $f = 0;
GetOptions('exclude' => \$exclude, 'regions=s{,}' => \@reg, 'c|offset=i' => \$offset, 'flanks=i' => \$f, 'o|output=s' => \$outfile, 'progress' => \$progress, 'input=s' => \$vcf_file, 'bed=s{,}' => \@bedfile, 'help' => \$help, 'manual' => \$manual) or pod2usage(-message => "Syntax error.", -exitval => 2);
pod2usage(-verbose => 2) if ($manual);
pod2usage(-verbose => 1) if ($help);
pod2usage(-message => "Syntax error.", -exitval => 2) if (not $vcf_file or (not @bedfile and not @reg));
if ($progress and $vcf_file eq "-"){
    print STDERR "Can't use progress option when input is from STDIN.\n";
    $progress = "";
}
my $time = strftime("%H:%M:%S", localtime);
print STDERR "Time started: $time\n";
$offset = 0 if not $offset;
my $vcf_obj = ParseVCF->new( file=> $vcf_file);
my $progressbar ;
if ($progress){
    $total_lines = $vcf_obj->countLines("variants");
    print STDERR "vcf file has $total_lines variants\n";
        $progressbar = Term::ProgressBar->new({name => "Filtering", count => $total_lines, ETA => "linear", });    
}
my @regions = ();

if (@bedfile){
    foreach my $bedfile (@bedfile){
        open (BED, $bedfile) || die "can't open $bedfile";
        my @temp = ();
        while (my $temp =  <BED>){
            chomp $temp;
            next if not $temp;
            $temp =~ s/[:-]/\t/g;
            $temp =~ s/\,//g;
            push (@temp, $temp);
        }
    #    my @temp = (<BED>) =~ s/[:-]/\t/g;
        close BED;
        my $invalid = grep {!/^(chr)*[\dXYMU][\w\d]*\t[\d]+\t[\d]+/} @temp;
        warn "$invalid invalid region(s) found in $bedfile" if $invalid;
        push (@regions, @temp);
    }
}
if (@reg){
    foreach my $reg (@reg){
        $reg =~ s/\,//g;
        die "invalid region specied: $reg" if $reg !~ /^(chr)*[\dXYMU][\w\d]*:[\d]+-[\d]+$/; 
        $reg =~ s/[\:\-]/\t/g;
        push (@regions, $reg);
    }
}

my $reg_obj = SortGenomicCoordinates -> new(array => \@regions, type => 'bed', col => 1);
$reg_obj->prep();
my $OUT;
if ($outfile){
    open ($OUT, ">$outfile") || die "Can't open $outfile for writing\n";
}else{
    $OUT = *STDOUT;
}
my $lines = 0;
my $next_update = 0;
my $prev_percent = 0;
my $printed = 0;
my $filtered = 0;
print $OUT $vcf_obj->getHeader ;
while (my $varLine = $vcf_obj->readLine){
    my @line = split("\t", $varLine);
    if ($reg_obj->locate(chrom => $line[$offset+0], position => $line[$offset+1]) > -1){
        if (not $exclude){
            print $OUT "$varLine\n";
            $printed++;
        }else{
            $filtered++;
        }
    }else{
        if ($exclude){
            print $OUT "$varLine\n";
            $printed++;
        }else{
            $filtered++;
        }
    }
    $lines++;
    if ($progress){
        $next_update = $progressbar->update($lines) if $lines >= $next_update;
    }
}
if ($progress){
    $progressbar->update($total_lines) if $total_lines >= $next_update;
}
$time = strftime("%H:%M:%S", localtime);
print STDERR "Time finished: $time\n";
print STDERR "$filtered variants filtered, $printed variants retained.\n";

#############################################
        
############################################
sub percentcomplete {
    my ($processedLines, $totalLines, $previousPercentComplete) = @_;
    my $percentComplete = int(($processedLines/$totalLines)*100);
    if ($percentComplete != $$previousPercentComplete){
        if ($percentComplete % 10 == 0 or $processedLines == $totalLines){
            $$previousPercentComplete = $percentComplete;
            return sprintf ("%7d out of %7d processed... %3d %% complete\n", $processedLines, $totalLines, $percentComplete);
        }
        else{
            $$previousPercentComplete = $percentComplete;
            return 0;
        }
    }else{
        return 0;
    }
}
