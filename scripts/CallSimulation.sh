#!/bin/bash

#$ -cwd
#$ -pe mpi $NUM_THREADS
#$ -S /bin/bash
#$ -v PATH

CONFIG_PATH=~/hla/hlaforest/scripts/config.sh
source $CONFIG_PATH

OUT_DIR=$1
READ_LENGTH=$2
NUM_READS=$3
INSERT_SIZE=$4
GENE_LIST=$5
ERROR_RATE=$6

OUT_PREFIX=$OUT_DIR/sim
FILTERED_1=$OUT_PREFIX\_1.fa
FILTERED_2=$OUT_PREFIX\_2.fa
SIM_1=$OUT_DIR/sim_1.fa
SIM_2=$OUT_DIR/sim_2.fa
TMP_DIR=$OUT_DIR/tmp


$SCRIPT_PATH/SimReads.pl -pe -l $READ_LENGTH -n $NUM_READS -f /home/hyjkim/hla/reference/hla_nuc_nospace_filtered.fasta -i $INSERT_SIZE -g $GENE_LIST -o $OUT_PREFIX -e $ERROR_RATE

sh $SCRIPT_PATH/build-forest-from-fasta.sh $SIM_1 $SIM_2 $OUT_DIR

$SCRIPT_PATH/prune.pl -f `echo $TMP_DIR/*.forest | tr " " ","` -t 2 > $TMP_DIR/t2.txt
$SCRIPT_PATH/prune.pl -f `echo $TMP_DIR/*.forestpruned | tr " " ","` -t 3 > $TMP_DIR/t3.txt
$SCRIPT_PATH/prune.pl -f `echo $TMP_DIR/*.forestprunedpruned | tr " " ","` -t 4 > $TMP_DIR/t4.txt
$SCRIPT_PATH/prune.pl  -f `echo $TMP_DIR/*.forestprunedprunedpruned | tr " " ","` -t 5 > $TMP_DIR/t5.txt

$SCRIPT_PATH/call-haplotypes.pl -2 $TMP_DIR/t2.txt -3 $TMP_DIR/t3.txt -4 $TMP_DIR/t4.txt -5 $TMP_DIR/t5.txt > $OUT_DIR/haplotypes.txt

$SCRIPT_PATH/score-sim-call.pl -t $OUT_DIR/sim_chosen_haplotypes.txt -p $OUT_DIR/haplotypes.txt > $OUT_DIR/sim-score.txt
