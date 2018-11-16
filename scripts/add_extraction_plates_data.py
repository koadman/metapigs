#!/usr/bin/env python3
import argparse
import pandas
import sys

parser = argparse.ArgumentParser(description='Adding extraction plates info to NCBI matrix')
parser.add_argument('NCBI_table', help='Excel spreadsheet from NCBI')
parser.add_argument('extraction_table', help='Excel spreadsheet of extraction plates')
args = parser.parse_args()

#this below is working:
#it's going to grab automatically the first file given in the command in the terminal to run this piece of code
NCBI_table = pandas.read_excel(args.NCBI_table,sheet_name=0,skiprows=14,skip_footer=0,index_col=None)
extraction_table = pandas.read_excel(args.extraction_table,sheet_name=0,index_col=None)

#dictionary can also be written this way: dict() = {}
sample_map=dict()
for row in range (extraction_table.shape[0]):
    if not extraction_table.iat[row,6] in sample_map:
        sample_map[extraction_table.iat[row,6]]=[(extraction_table.iat[row,0], extraction_table.iat[row,1])]
    else:
        sys.stderr.write("WARNING: found extra DNA sample for poop "+extraction_table.iat[row,6]+"\n")
        sample_map[extraction_table.iat[row,6]].append((extraction_table.iat[row,0], extraction_table.iat[row,1]))

print('\t'.join(NCBI_table.columns)+'\tDNA_plate\tDNA_well')
used = {}
for row in range(NCBI_table.shape[0]):
        used[NCBI_table.iat[row,0]]=1
        if not NCBI_table.iat[row,0] in sample_map:
            sys.stderr.write("WARNING: DNA sample not found for poop "+NCBI_table.iat[row,0]+"\n")
        else:
            for i in range(len(sample_map[NCBI_table.iat[row,0]])):
                _t = sample_map[NCBI_table.iat[row,0]][i]
                ncbi_row = NCBI_table.iloc[row,:]
                if i>0:
                    ncbi_row[0]=ncbi_row[0]+'.'+str(i+1)
                print('\t'.join(map(str,ncbi_row))+'\t'+"\t".join(_t))
#        print("\t".join(sample_map[NCBI_table.iat[row,0]]))

for s in sample_map.items():
    if not s[0] in used:
        sys.stderr.write("WARNING: DNA sample for poop "+s[0]+" is missing metadata\n")
