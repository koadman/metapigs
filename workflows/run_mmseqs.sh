#!/bin/bash

##############################################################################################
#
# Usage: As a qsub submission with variables METAPIGS_REPO and RUN_TABLE
#
# > qsub -l select=1:mem=2g:ncpus=1 -v METAPIGS_REPO=$HOME/metapigs,READ_TABLE=reads.tsv
#
##############################################################################################

#
# Switch to working directory
#
cd $PBS_O_WORKDIR

#
# Check for env variable
#
if [ -z $METAPIGS_REPO ]
then
    echo "METAPIGS_REPO has not been set."
    exit 1
fi

#
#  Run nextflow on an execution host.
#
#  It is possible to test run using "--debug" option.
#
#  In testing, the Phix genome and adapter sequences were being sourced from BBMAP. However,
#  as these tools are satisfied using conda, the path to these files is cumbersome to determine.
#  It would be better to place these sequences at some location in the repo. I have presently
#  put the BBMAP Phix genome in workflows.
#
#  NOTE: double hyphens are user parameters, single hyphen args are to nextflow itself.
#
#  Some settings are defined in nextflow.config
#   - default process memory and concurrency
#   - queue size
#   - conda environment to install
#   - conda cache dir
#
nextflow -C $METAPIGS_REPO/workflows/mmseqs.config run $METAPIGS_REPO/workflows/mmseqs.nf \
	-resume \
	-with-report \
        -profile cluster \
        --targetDB $TARGET_DB \
        --out_dir $OUT \
        --read_table $READ_TABLE \
	--raw_dir $RAWDIR
