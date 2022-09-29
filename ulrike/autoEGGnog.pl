#! /usr/bin/perl -w
use strict;
#author ulrike.loeber@mdc-berlin.de


#create log
my $date=localtime();
open (LOG, '>>', 'CelticFireAnnotationEggnog.log') or die $!;
print LOG "$date\n";
close LOG;
open(STDERR, '>>', 'CelticFireAnnotationEggnog.log') or die "Can't open log";

open (INFILE,$ARGV[0]) or die $!;
my @sample=<INFILE>;
close INFILE;
#Cutibacterium_acnes_20925_1_38/annotation/Cutibacterium_acnes_20925_1_38.faa
my $i=1;
	 open (DELFILE,">QDELeggnog.sh");

foreach my $dir(@sample){
	open (OUTFILE,">eggnog.$i.sh") or die $!;

        chomp $dir;
        my $file=$dir;
        $file="/annotation/".$dir.".faa";


        my $in=$file;
        my $out="/".$dir.".eggnog";


  print OUTFILE   "#!/bin/bash

#\$ -N eggnog.$i
#\$ -cwd
#\$ -l h_vmem=246G
#\$ -l h_rt=40:00:00
python2.7 emapper.py -i $dir$in --output $dir$out -m diamond -d bact --override"; 

        
        close OUTFILE;
        print DELFILE "qdel eggnog.$i\n";
        $i++;
}
close DELFILE;

for (my $j=1;$j<=$i;$j++){
        system "qsub eggnog.$j.sh";
        print "qsub eggnog.$j.sh\n";
}

