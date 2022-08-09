#! /usr/bin/perl -w
use strict;
#author: Ulrike Löber (ulrike.loeber@mdc-berlin.de)
#head binVFDB.txt 
#TAG	keyword	low level
#adsA	Immune evasion	Immune evasion
#bexA	Toxin	Toxin
#bexB	Toxin	Toxin
#cap8B	Toxin; Membrane-damaging; Pore-forming toxin	Toxin
#cap8C	Toxin; Membrane-damaging; Pore-forming toxin	Toxin
#cap8D	Toxin; Membrane-damaging; Pore-forming toxin	Toxin
#cap8E	Toxin; Membrane-damaging; Pore-forming toxin	Toxin
my %virumap;
my %categories;
open (VFDBmap, "binVFDB.txt") or die $!;
my $line=0;
while(<VFDBmap>){
	chomp $_;
	if($line==0){next;}
	else{
		my @tmp=split(/\t/,$_);
		$virumap{$tmp[0]}=$tmp[2];		
		$categories{$tmp[2]}="";
	}
}
close VFDBmap;

#head ariba-vfdb_core-summary.renamed.txt 
#sample_name	adsA	bexA	bexB	cap8B	cap8C	cap8D	cap8E	cap8F	cap8G	cap8L	cap8M	cap8N	cap8O	cap8P	cps4D	cps4J	cps4K	cps4L	ctrD	ebp	esaA	esaB	esaG1	essA	essB	esxA	geh	hld	hlgA	hlgB	hlgC	hly_hla	hscA	hscB	hysA	hysA_1	icaA	icaB	icaC	icaD	icaR	isdA	isdB	isdC	isdD	isdE	isdF	isdG	katA_1	lbpA	liplipA	lipB	lytB	lytC	msrA_B_pilB_	nanB	pavA	pce	ply	psaA	rffG	sbi	scn	sdrE	sipA	spa	srtB_1	sspB	sspC
#Actinomyces_bouchesdurhonensis_20925_1_25	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	00	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	00	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
#Actinomyces_naeslundii_27098_8_156	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	00	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	00	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
#Actinomyces_oris_27098_8_157	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	00	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	00	0	0	0	0	0	0	0	0	0

open (aribaVFDB, "ariba-vfdb_core-summary.renamed.txt") or die $!;
$line=0;
my @genes;
my %vfcount;
while(<aribaVFDB>){
	chomp $_;
	if ($line==0){
		@genes=split(/\t/,$_);
		$line++;
	}
	else{
		my @temp=split(/\t/,$_);
		for(my $e = 1;$e>@genes;$e++){
			if($temp[$e]>=1){
#			$vfcount{$temp[0]}{$genes[$e]}=$temp[$e];
			$vfcount{$temp[0]}{$virumap{$genes[$e]}}{$genes[$e]}=$temp[$e];
		
			}

		}
	}
}
close aribaVFDB;

foreach my $key(sort keys %vfcount){
	foreach my $cat(sort keys %categories){
		my $count=0;
		foreach my $gene(sort keys %{$vfcount{$key}{$cat}}){
			$count=$count+$vfcount{$key}{$cat}{$gene};
		}
		print "$key\t$cat\t$count\n";
	}
}
#print OUTFILE "
#DATASET_TEXT
##In text datasets, each ID is associated to text label, which can be displayed directly on the node branch, or outside the tree
##lines starting with a hash are comments and ignored during parsing
##=================================================================#
##                    MANDATORY SETTINGS                           #
##=================================================================#
##select the separator which is used to delimit the data below (TAB,SPACE or COMMA).This separator must be used throughout this file.
#SEPARATOR TAB

##label is used in the legend table (can be changed later)
#DATASET_LABEL,example text dataset

##dataset color (can be changed later)
#COLOR,#ff0000

##=================================================================#
##                    OPTIONAL SETTINGS                            #
##=================================================================#

##=================================================================#
##     all other optional settings can be set or changed later     #
##           in the web interface (under 'Datasets' tab)           #
##=================================================================#

##left margin, used to increase/decrease the spacing to the next dataset. Can be negative, causing datasets to overlap. Used only for text labels which are displayed on the outside
#MARGIN,0

##applies to external text labels only; if set, text labels associated to internal nodes will be displayed even if these nodes are not collapsed. It could cause overlapping in the dataset display.
#SHOW_INTERNAL,0

##Rotate all labels by the specified angle
#LABEL_ROTATION,0

##By default, internal labels will be placed above the branches. If LABELS_BELOW is set to 1, labels will be below the branches
#LABELS_BELOW,1

##Shift internal labels vertically by this amount of pixels (positive or negative)
#VERTICAL_SHIFT,0

##If set to 1, tree rotation will not influence the individual label rotation
#STRAIGHT_LABELS,0

##applies to external text labels only; If set to 1, labels will be displayed in arcs aligned to the tree (in circular mode) or vertically (in normal mode). All rotation parameters (global or individual) will be ignored.
#ALIGN_TO_TREE,0

##font size factor; For external text labels, default font size will be slightly less than the available space between leaves, but you can set a multiplication factor here to increase/decrease it (values from 0 to 1 will decrease it, values above 1 will increase it)
#SIZE_FACTOR,1

##add extra horizontal shift to the external labels. Useful in unrooted display mode to shift text labels further away from the node labels.
#EXTERNAL_LABEL_SHIFT,0

##Internal tree nodes can be specified using IDs directly, or using the 'last common ancestor' method described in iTOL help pages
##=================================================================#
##       Actual data follows after the "DATA" keyword              #
##=================================================================#
##the following fields are possible for each node:
##ID,label,position,color,style,size_factor,rotation

##position defines the position of the text label on the tree:
##  -1 = external label
##  a number between 0 and 1 = internal label positioned at the specified value along the node branch (for example, position 0 is exactly at the start of node branch, position 0.5 is in the middle, and position 1 is at the end)
##style can be 'normal',''bold','italic' or 'bold-italic'
##size factor will be multiplied with the standard font size

#DATA
##Examples

##node 9606 will have an external label 'Homo sapiens' in bold red and twice the size of standard labels
##9606,Homo sapiens,-1,#ff0000,bold,2,0

##node 4530 will have an internal label 'Oryza sativa' in bold italic blue, starting directly over the node
##4530,Oryza sativa,0,#0000ff,bold-italic,1
#";
