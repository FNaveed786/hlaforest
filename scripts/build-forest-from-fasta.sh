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

MAX_INSERT=1000

mkdir -p $TMP_DIR

bowtie -X $MAX_INSERT -f -p $NUM_THREADS --al $FILTERED_PREFIX $BOWTIE_INDEX -1 $READ_ONE -2 $READ_TWO > /dev/null
# This was used to generate uniquely mapping reads
#bowtie -X $MAX_INSERT -f -n 0 -e 1 -p 8 -S -a $BOWTIE_INDEX -1 $FILTERED_PREFIX\_1 -2 $FILTERED_PREFIX\_2 > $ALIGNED_PAIRS
# Generate all reads
bowtie -X $MAX_INSERT -f -n 0 -e 1 -p $NUM_THREADS -S -a $BOWTIE_INDEX -1 $FILTERED_PREFIX\_1 -2 $FILTERED_PREFIX\_2 > $ALIGNED_PAIRS

$SCRIPT_PATH/build-forest.pl -reads $ALIGNED_PAIRS -o $FOREST_PAIRS -verbose -n $TREES_IN_FOREST
