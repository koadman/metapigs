#!/bin/bash

##############################################################################################
#
# Usage: As a qsub submission with variables METAPIGS_REPO and RUN_TABLE
#
# > qsub -l select=1:mem=2g:ncpus=1 -v METAPIGS_REPO=$HOME/metapigs,RUN_TABLE=mytable.csv
#
# The run table is expected to be tab delimited -- despite the suffix implying otherwise
#
# Headers are expected to be:
#
#        *sample_name     - unique sample name
#        r1_filename      - path to R1
#        r2_filename      - path to R2
#        isolation_source - source (pig) id
#
# Although a pain, absolute paths to R1/R2 will save file-not-found headaches if the user
# fines relative paths confusing.
#
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
# Set scratch location for nextflow processes, though this shouldn't be necessary here
#
export METAPIGS_TMP=/scratch/work/
export TMP=$METAPIGS_TMP
export TMPDIR=$TMP

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
nextflow -C $METAPIGS_REPO/workflows/initial.config run $METAPIGS_REPO/workflows/initial.nf \
	-with-report \
        -profile cluster \
        --run_table $RUN_TABLE \
        --adapters $METAPIGS_REPO/source_data/custom_adapters.fa \
        --phix $METAPIGS_REPO/source_data/phix174_ill.ref.fa
