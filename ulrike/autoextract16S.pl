#! /usr/bin/perl -w
use strict;
#author ulrike.loeber@mdc-berlin.de

open (INFILE,$ARGV[0]) or die $!;
my @sample=<INFILE>;
close INFILE;

foreach my $dir(@sample){

        chomp $dir;
	my $in=$dir;
	$in=$dir."/assembly/".$in.".fna";
	my $out=$in;
	$out=~s/\.fna/_16S\.metaxa2/g;
	print "metaxa2 -i $in -o $out\n";
	system "metaxa2 -i $in -o $out";
}
