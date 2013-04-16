#!/bin/bash

#$ -cwd
#$ -pe mpi 4
#$ -S /bin/bash
#$ -v PATH

CONFIG_PATH=~/hla/hlaforest/scripts/config.sh
source $CONFIG_PATH

OUT_DIR=$1
READ_1=$2
READ_2=$3

OUT_PREFIX=$OUT_DIR/sim
FILTERED_1=$OUT_PREFIX\_1.fa
FILTERED_2=$OUT_PREFIX\_2.fa
SIM_1=$OUT_DIR/sim_1.fa
SIM_2=$OUT_DIR/sim_2.fa
TMP_DIR=$OUT_DIR/tmp

mkdir -p $TMP_DIR

#sh $SCRIPT_PATH/build-forest-from-fastq.sh $READ_1 $READ_2 $OUT_DIR

$SCRIPT_PATH/prune2.pl -threshold 0.05 -f `echo $TMP_DIR/*.forest | tr " " ","` -t 2 > $TMP_DIR/t2.txt
$SCRIPT_PATH/prune2.pl -threshold 0.05 -f `echo $TMP_DIR/*.forestpruned | tr " " ","` -t 3 > $TMP_DIR/t3.txt
$SCRIPT_PATH/prune2.pl -threshold 0.05 -f `echo $TMP_DIR/*.forestprunedpruned | tr " " ","` -t 4 > $TMP_DIR/t4.txt
$SCRIPT_PATH/prune2.pl -threshold 0.05 -f `echo $TMP_DIR/*.forestprunedprunedpruned | tr " " ","` -t 5 > $TMP_DIR/t5.txt

$SCRIPT_PATH/call-haplotypes.pl -threshold .3 -2 $TMP_DIR/t2.txt -3 $TMP_DIR/t3.txt -4 $TMP_DIR/t4.txt -5 $TMP_DIR/t5.txt > $OUT_DIR/haplotypes.txt

