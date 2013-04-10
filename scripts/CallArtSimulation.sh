#!/bin/bash

#$ -cwd
#$ -pe mpi 4
#$ -S /bin/bash
#$ -v PATH

CONFIG_PATH=~/hla/hlaforest/scripts/config.sh
source $CONFIG_PATH

OUT_DIR=$1
IMGT_REF_FA=$2
MIN_REF_LENGTH=$3
#GENE_LIST=$4

OUT_PREFIX=$OUT_DIR/sim
SIM_REF_FA=$OUT_DIR/ref.fa
SIM_REF_LOG=$OUT_DIR/sim_chosen_haplotypes.txt
SIM_1=$OUT_PREFIX\1.fq
SIM_2=$OUT_PREFIX\2.fq
TMP_DIR=$OUT_DIR/tmp

# Select haplotypes

$SCRIPT_PATH/SelectHaplotypes.pl -f $IMGT_REF_FA -minRefLength $MIN_REF_LENGTH -o $SIM_REF_FA -l $SIM_REF_LOG

# Simulate reads with ART
~/bin/art_illumina -i $SIM_REF_FA -o $OUT_PREFIX  --len 50 --fcov 100 -sam -p --mflen 250 --sdev 25 --insRate 0 --delRate 0 --insRate2 0 --delRate2 0 

# Call haplotypes with existing script
$SCRIPT_PATH/CallHaplotypesPE.sh $OUT_DIR $SIM_1 $SIM_2

# Score the called haplotypes with the true haplotypes
perl $SCRIPT_PATH/score-sim-call.pl -t $OUT_DIR/sim_chosen_haplotypes.txt -p $OUT_DIR/haplotypes.txt > $OUT_DIR/sim-score.txt
