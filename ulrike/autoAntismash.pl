use strict;

#author: Ulrike Löber (ulrike.loeber@mdc-berlin.de)
#create log
my $date=localtime();
open (LOG, '>>', 'antismasheach.log') or die $!;
print LOG "$date\n";
close LOG;
open(STDERR, '>>', 'antismasheach.log') or die "Can't open log";

open (INFILE,$ARGV[0]) or die $!;
my @sample=<INFILE>;
close INFILE;
open (DELFILE,">QDELantismash.sh");
my $i=1;

foreach my $dir(@sample){
        chomp $dir;
	my $in=$dir."/annotation/".$dir.".gbk";
	my $out=$dir."/antismash";
	my $name="antismash".$dir;
        open (OUTFILE,">antismasheach.$i.sh") or die $!;

        print OUTFILE   "#!/bin/bash


#\$ -l os=centos7
#\$ -N $name
#\$ -cwd
#\$ -l h_vmem=20G
#\$ -pe smp 4
#\$ -l h_rt=01:00:00
source ~/.bashrc
conda activate antismash
antismash --clusterblast --subclusterblast --asf --full-hmmer --smcogs --tta --cpus 4 $in --outputfolder $out
";

        close OUTFILE;
        print DELFILE "qdel $name\n";
        $i++;
}
close DELFILE;

for (my $j=2;$j<=($i-1);$j++){
        system "qsub antismasheach.$j.sh";
        print "qsub antismasheach.$j.sh\n";
}



