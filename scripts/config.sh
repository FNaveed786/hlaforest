#!/bin/bash

# User modifiable variables
HLAFOREST_HOME=/home/hyjkim/hla/hlaforest
SCRIPT_PATH=$HLAFOREST_HOME/scripts # installation directory of your scripts
BOWTIE_INDEXES=$HLAFOREST_HOME/reference # path to bowtie references directory
NUM_THREADS=8 # the number of threads bowtie should use
TREES_IN_FOREST=50000 # the number of trees to output in each forest file


PERL5LIB=$PERL5LIB:$SCRIPT_PATH
BOWTIE_INDEX=hla_nuc_nospace_filtered


export PERL5LIB
export BOWTIE_INDEXES
export BOWTIE_INDEX
export NUM_THREADS
export SCRIPT_PATH
export TREES_IN_FOREST
