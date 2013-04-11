#!/bin/bash

#$ -cwd
#$ -pe mpi 4
#$ -S /bin/bash
#$ -v PATH

CONFIG_PATH=~/hla/hlaforest/scripts/config.sh
source $CONFIG_PATH

OUT_DIR=$1
IMGT_REF_FA=$2
READ_LENGTH=$3
FRAG_LENGTH=$4
FRAG_LENGTH_SD=$5
COVERAGE=$6
Q_SHIFT=$7
GENE_LIST=$8

if [ "$GENE_LIST" ]; then
GENE_FLAG="-g $GENE_LIST"
fi

OUT_PREFIX=$OUT_DIR/sim
SIM_REF_FA=$OUT_DIR/ref.fa
SIM_REF_LOG=$OUT_DIR/sim_chosen_haplotypes.txt
SIM_1=$OUT_PREFIX\1.fq
SIM_2=$OUT_PREFIX\2.fq
ALN_1=$OUT_PREFIX\1.aln
ALN_2=$OUT_PREFIX\2.aln
TMP_DIR=$OUT_DIR/tmp

SUBRATE_1=$OUT_PREFIX\1\_sub_rate.tab
SUBRATE_2=$OUT_PREFIX\2\_sub_rate.tab

# Select haplotypes

$SCRIPT_PATH/SelectHaplotypes.pl -f $IMGT_REF_FA -minRefLength $MIN_REF_LENGTH -o $SIM_REF_FA -l $SIM_REF_LOG $GENE_FLAG

# Simulate reads with ART
/bin/art_illumina -i $SIM_REF_FA -o $OUT_PREFIX  --len $READ_LENGTH --fcov $COVERAGE -sam -p --mflen $FRAG_LENGTH --sdev $FRAG_LENGTH_SD --insRate 0 --delRate 0 --insRate2 0 --delRate2 0 -qs $Q_SHIFT -qs2 $Q_SHIFT

# Get substitution rate of aln files
perl $SCRIPT_PATH/aln2subrate.pl < $ALN_1 > $SUBRATE_1
perl $SCRIPT_PATH/aln2subrate.pl < $ALN_2 > $SUBRATE_2

# Call haplotypes with existing script
$SCRIPT_PATH/CallHaplotypesPE.sh $OUT_DIR $SIM_1 $SIM_2

# Score the called haplotypes with the true haplotypes
perl $SCRIPT_PATH/score-sim-call.pl -t $OUT_DIR/sim_chosen_haplotypes.txt -p $OUT_DIR/haplotypes.txt > $OUT_DIR/sim-score.txt

# Cleanup tmp dir
rm -rf $OUT_DIR/tmp
