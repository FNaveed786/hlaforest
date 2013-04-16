#!/bin/bash
#$ -V
#$ -cwd 
#$ -S /bin/bash
#$ -pe mpi $NUM_THREADS

READ_ONE=$1
READ_TWO=$2
OUT_DIR=$3

TMP_DIR=$OUT_DIR/tmp
FILTERED_PREFIX=$TMP_DIR/tmp_hla_hits
ALIGNED_PAIRS=$TMP_DIR/hla_aligned_pairs.sam
#TREE_PAIRS=$TMP_DIR/hla_aligned_pairs.tree
FOREST_PAIRS=$TMP_DIR/hla_aligned_pairs

mkdir -p $TMP_DIR

bowtie -p $NUM_THREADS --al $FILTERED_PREFIX $BOWTIE_INDEX -1 $READ_ONE -2 $READ_TWO > /dev/null
# Exact match results in less required memory
bowtie -n 0 -e 1 -p $NUM_THREADS --sam --all $BOWTIE_INDEX -1 $FILTERED_PREFIX\_1 -2 $FILTERED_PREFIX\_2 > $ALIGNED_PAIRS

# Build the forest from the sam alignment
$SCRIPT_PATH/build-forest.pl -reads $ALIGNED_PAIRS -o $FOREST_PAIRS -b $TREES_IN_FOREST $MAX_TREES_FLAG
