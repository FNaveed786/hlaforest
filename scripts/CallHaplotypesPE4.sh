#!/bin/bash

#$ -cwd
#$ -pe mpi 1
#$ -S /bin/bash
#$ -v PATH

CONFIG_PATH=~/hla/hlaforest/scripts/config.sh
source $CONFIG_PATH

OUT_DIR=$1
READ_1=$2
READ_2=$3

OUT_PREFIX=$OUT_DIR/sim
SIM_1=$OUT_DIR/sim_1.fa
SIM_2=$OUT_DIR/sim_2.fa
TMP_DIR=$OUT_DIR/tmp
FILTERED_PREFIX=$TMP_DIR/tmp_hla_hits
ALIGNED_PAIRS=$TMP_DIR/hla_aligned_pairs.sam
FOREST_PAIRS=$TMP_DIR/hla_aligned_pairs

mkdir -p $TMP_DIR

#sh $SCRIPT_PATH/build-forest-from-fastq.sh $READ_1 $READ_2 $OUT_DIR
# filter reads
bowtie	-p $NUM_THREADS --al $FILTERED_PREFIX $BOWTIE_INDEX -1 $READ_1 -2 $READ_2> /dev/null
# align reads
bowtie -p $NUM_THREADS --sam --all $BOWTIE_INDEX -1 $FILTERED_PREFIX\_1 -2 $FILTERED_PREFIX\_2 > $ALIGNED_PAIRS

$SCRIPT_PATH/build-forest.pl -reads $ALIGNED_PAIRS -o $FOREST_PAIRS -b $TREES_IN_FOREST $MAX_TREES_FLAG

$SCRIPT_PATH/prune4.pl -threshold 0.05 -f `echo $TMP_DIR/*.forest | tr " " ","` > $OUT_DIR/haplotypes.txt

