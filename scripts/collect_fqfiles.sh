#!/bin/bash

SEQ_FOLDER=MON5838

# collect a flat list of all FastQ files in sequencing run folder
# the extract the plate and well as new columns 
find $SEQ_FOLDER -type f -name '*.fastq.gz' | \
    awk '{match($1,/plate_([0-9]+)_([a-zA-Z]+[0-9]+)_/,a); print "P" a[1],a[2],$0}' | sort -gk1,2 > fq_well_and_plate.txt
