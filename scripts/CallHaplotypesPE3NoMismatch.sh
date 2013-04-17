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

sh $SCRIPT_PATH/build-forest-from-fastq-no-mismatch.sh $READ_1 $READ_2 $OUT_DIR

$SCRIPT_PATH/prune3.pl -threshold 0.05 -f `echo $TMP_DIR/*.forest | tr " " ","` > $OUT_DIR/haplotypes.txt

