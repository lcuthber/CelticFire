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
	 open (DELFILE,">QDELhmmer.sh");

foreach my $dir(@sample){
	open (OUTFILE,">hmmer.$i.sh") or die $!;

        chomp $dir;
        my $file=$dir;
        $file="/annotation/".$dir.".faa";


        my $in=$file;
        my $out="/".$dir.".pfam35.txt";


  print OUTFILE   "#!/bin/bash

#\$ -N hmmer.$i
#\$ -cwd
#\$ -l h_vmem=100G
#\$ -l h_rt=10:00:00
source /home/user/.bashrc
hmmsearch --tblout $dir$out -E 1e-5 --cpu 1 Pfam-A.hmm $dir$in";
        
        close OUTFILE;
        print DELFILE "qdel hmmer.$i\n";
        $i++;
}
close DELFILE;

for (my $j=1;$j<=$i;$j++){
        system "qsub hmmer.$j.sh";
        print "qsub hmmer.$j.sh\n";
}
