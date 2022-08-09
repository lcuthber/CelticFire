#! /home/uloeber/.guix-profile/bin/perl -w
use strict;

#author: Ulrike Löber (ulrike.loeber@mdc-berlin.de)
#create log
my $date=localtime();
open (LOG, '>>', 'bactopiaeach.log') or die $!;
print LOG "$date\n";
close LOG;
open(STDERR, '>>', 'bactopiaeach.log') or die "Can't open log";

open (INFILE,$ARGV[0]) or die $!;
my @sample=<INFILE>;
close INFILE;
open (DELFILE,">QDELbactopia.sh");
my $i=1;

foreach my $dir(@sample){
        chomp $dir;
        my $file=$dir;
        $file=~s/\//_contigs.fa/g;

        my $prefix=$dir;
        $prefix=~s/\///g;
        my $fwd="RAW/".$prefix.".1.fastq.gz";
        my $rev="RAW/".$prefix.".2.fastq.gz";
        my $name="job".$prefix;
        open (OUTFILE,">bactopiaeach.$i.sh") or die $!;

        print OUTFILE   "#!/bin/bash


#\$ -l os=centos7
#\$ -N $name
#\$ -cwd
#\$ -l h_vmem=100G
#\$ -l h_rt=10:00:00
source ~/anaconda3/etc/profile.d/conda.sh
conda activate bactopia
bactopia --R1 $fwd --R2 $rev --sample $prefix --datasets datasets/ --outdir BACTOPIA

";

        close OUTFILE;
        print DELFILE "qdel $name\n";
        $i++;
}
close DELFILE;

for (my $j=1;$j<=($i-1);$j++){
        system "qsub bactopiaeach.$j.sh";
        print "qsub bactopiaeach.$j.sh\n";
}



