#! /usr/bin/perl -w
use strict;

#author: Ulrike Löber (ulrike.loeber@mdc-berlin.de)
#create log
my $date=localtime();

open (LOG, '>>', 'mergeISO.log') or die $!;
#print LOG "$date\n";



#head sample_species.txt 
#ID      species
#20925_1_10      Streptococcus parasanguinis
#20925_1_22      Streptococcus sanguinis
#20925_1_23      Neisseria subflava
#20925_1_24      Streptococcus mitis
#20925_1_25      Actinomyces bouchesdurhonensis
#20925_1_26      Streptococcus pseudopneumoniae
#20925_1_62      Staphylococcus epidermidis
#20925_1_63      Staphylococcus capitis
#20925_1_65      Gemella sanguinis

my %specieshash;

open (SPECIN,"sample_species.txt") or die $!;
while (<SPECIN>){
	chomp $_;
	my ($key, $value) = split /\t/, $_;
	next LINE if not $key;
	$value=~s/ /_/g;
	$specieshash{$key} = $value;

}
close SPECIN;
open (DELFILE,">delBACTOPIA.sh") or die $!;
#$ ls RAW/
#16234_1_10.1.fastq.gz  20925_1_29.1.fastq.gz  20925_1_59.1.fastq.gz  20925_1_89.1.fastq.gz   27098_8_123.1.fastq.gz  27098_8_153.1.fastq.gz  27098_8_183.1.fastq.gz  27098_8_42.2.fastq.gz  27098_8_72.2.fastq.gz
#16234_1_10.2.fastq.gz  20925_1_29.2.fastq.gz  20925_1_59.2.fastq.gz  20925_1_89.2.fastq.gz   27098_8_123.2.fastq.gz  27098_8_153.2.fastq.gz  27098_8_183.2.fastq.gz  27098_8_4.2.fastq.gz   27098_8_7.2.fastq.gz
#16234_1_11.1.fastq.gz  20925_1_30.1.fastq.gz  20925_1_60.1.fastq.gz  20925_1_90.1.fastq.gz   27098_8_124.1.fastq.gz  27098_8_154.1.fastq.gz  27098_8_184.1.fastq.gz  27098_8_43.1.fastq.gz  27098_8_73.1.fastq.gz

#head clusters.txt 

#20925_1_10      20925_1_17      20925_1_35      20925_1_39      20925_1_57      20925_1_71      27098_8_10      27098_8_114     27098_8_20      27098_8_26      27098_8_31      27098_8_3       27098_8_5       27098_8_8  27098_8_9
#20925_1_11      20925_1_18      27098_8_11
#20925_1_13
#20925_1_15
#20925_1_16      20925_1_82
#20925_1_19
#20925_1_21      20925_1_47      20925_1_54      27098_8_21      27098_8_24      27098_8_28
#20925_1_22
#20925_1_23      20925_1_40      20925_1_67      20925_1_73      20925_1_79
my $i=0;
open (CLUSTER, "clusters.txt") or die $!;
while (<CLUSTER>){
	chomp $_;
	my @members = split /\t/, $_;
	my $species=$specieshash{$members[0]};
	my $repres=$members[0];
	my @R1string;
	my @R2string;
	my $outR1="reassemble220920/".$species."_".$repres.".R1.fastq.gz";
	my $outR2="reassemble220920/".$species."_".$repres.".R2.fastq.gz";
	my $prefix=$species."_".$repres;
	foreach my $isolate(@members){
		if($species !~ /$specieshash{$isolate}/){print "WARNING $isolate $species !~ $specieshash{$isolate}\n";}
		my $fileR1="RAW/".$isolate.".1.fastq.gz";		
		my $fileR2="RAW/".$isolate.".2.fastq.gz";		
		push (@R1string,$fileR1);
		push (@R2string,$fileR2);
	}
#	print "cat @R1string >$outR1\n";
#	system "cat @R1string >$outR1";
#	print "cat @R2string >$outR2\n";
#	system "cat @R2string >$outR2";

        my $fwd=$outR1;
        my $rev=$outR2;
        my $name=$prefix.".log";
	$species=~s/_/ /g;

        open (OUTFILE,">bactopiaeach.$i.sh") or die $!;

        print OUTFILE   "#!/bin/bash


#\$ -l os=centos7
#\$ -N $name
#\$ -cwd
#\$ -l h_vmem=100G
#\$ -l h_rt=10:00:00
source ~/.bashrc
conda activate bactopia
bactopia --R1 $fwd --R2 $rev --sample $prefix --datasets datasets/  --outdir BACTOPIA141120mergedIso --species '$species'

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



close CLUSTER;
close LOG;

