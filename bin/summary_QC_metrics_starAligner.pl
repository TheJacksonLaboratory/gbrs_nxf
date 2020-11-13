#!/usr/bin/perl
use strict;
use warnings;

open(FILEIN3, $ARGV[0]) || die "cannot open the file"; ####*Log.final.out
open(FILEIN4, $ARGV[1]) || die "cannot open the file"; ####picard stats

my ($value1, $value2, $value3, $flag);

$flag =  0;

print "##################\n";
print "STAR Aligner Report\n";
print "##################\n";

while(my $readFile3 = <FILEIN3>)
{
   if(($readFile3 =~ /^\s*Number of input reads.*$/) && ($flag == 0))
    {
        chomp $readFile3;
        $readFile3 =~ s/^\s+//;
        $readFile3 =~ s/|//;
        print "$readFile3\n";
        $flag = 1;
        next;
    }
   elsif(($flag == 1))
    {
        chomp $readFile3;
        $readFile3 =~ s/^\s+//;
        $readFile3 =~ s/\|//;
        print "$readFile3\n";
    }

}

print "\n";

$flag = 0;

my @splitHeader = ();
my @splitValue  = ();

print "##################\n";
print "picard collectRNAseqMetrics Report\n";
print "##################\n";

while(my $readFile4 = <FILEIN4>)
{

   if(($readFile4 =~ /^\s*##\s+METRICS\s+CLASS\s+picard.analysis.RnaSeqMetrics.*$/) && ($flag == 0))
    {
           $flag = 1;
           next;
    }
   elsif(($readFile4 =~ /^\s*PF_BASES.*$/) && ($flag == 1))
    {
           $flag = 2;
           chomp $readFile4;
           @splitHeader = split("\t", $readFile4);
           next;
    }
   elsif($flag == 2)
    {
           $flag = 3;
           chomp $readFile4;
           @splitValue = split("\t", $readFile4);
           next;
    }

}

for (my $i =0; $i<=$#splitHeader-3; $i++)
{
   print "$splitHeader[$i]\t$splitValue[$i]\n";

}
