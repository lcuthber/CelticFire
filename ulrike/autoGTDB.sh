#!/bin/bash


#$ -l os=centos7
#$ -N GTDB
#$ -cwd
#$ -l h_vmem=500G
#$ -l h_rt=96:00:00
source ~/.bashrc
conda activate bactopia
bactopia tools gtdb --bactopia ./ --gtdb GTDBrelease202 -resume


