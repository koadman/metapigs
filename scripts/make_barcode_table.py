#!/usr/bin/env python3
import numpy as np
import sys
from Bio.Seq import Seq

i5 = np.genfromtxt(sys.argv[1],skip_header=1,delimiter='\t',dtype=None,encoding=None)
i7 = np.genfromtxt(sys.argv[2],skip_header=1,delimiter='\t',dtype=None,encoding=None)

bc_table = open('phylosift_barcodes.tsv','w')
for i in range(i5.shape[0]):
    bc_table.write(i5[i,4]+'\ti5_'+i5[i,0]+"\t0\t0\t0\t0\n")

for i in range(i7.shape[0]):
    bc_table.write(i7[i,4]+'\ti7_'+i7[i,0]+"\t0\t0\t0\t0\n")
bc_table.close()

sample_table = open('phylosift_samples.tsv','w')
sample_table.write("#plate_well\ti7_name\ti7_barcode\ti5_name\ti5_barcode\n")
for plate in range(10):
    # iterate through the plates, rotating the i7 barcode plate column
    for i in range(12):
        for j in range(8):
            i5_bc = j * 12 + i
            i7_bc = j * 12 + (i + 12 - plate) % 12
            i7_seq = Seq(i7[i7_bc,4])
            i7_rc = i7_seq.reverse_complement()
            sample_table.write("plate_"+str(plate+1)+"_"+i5[i5_bc,0]+"\ti7_bc_"+i7[i7_bc,0]+"\t"+str(i7_rc)+"\ti5_bc_"+i5[i5_bc,0]+"\t"+i5[i5_bc,4]+"\n")
sample_table.close()

