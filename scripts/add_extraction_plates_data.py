#!/usr/bin/env python
import argparse
import pandas

parser = argparse.ArgumentParser(description='Adding extraction plates info to NCBI matrix')
parser.add_argument('NCBI_table', help='Excel spreadsheet from NCBI')
parser.add_argument('extraction_table', help='Excel spreadsheet of extraction plates')
args = parser.parse_args()

#this below is working:
#it's going to grab automatically the first file given in the command in the terminal to run this piece of code
NCBI_table = pandas.read_excel(args.NCBI_table,sheet_name=None,skiprows=14,skip_footer=0,index_col=None)
print(NCBI_table.iloc[:,0])

#converters = {col: str for col in column_list}
#extraction_table = pandas.read_excel(args.extraction_table,sheet_name=0,skiprows=1,skip_footer=0,index_col=None,converters=converters)
#print(extraction_table)

column_list = []
df_column = pandas.read_excel(args.extraction_table,sheet_name=0,index_col=None).columns
for i in df_column:
    column_list.append(i)
converter = {col: str for col in column_list}
extraction_table = pandas.read_excel(args.extraction_table,sheet_name=0,index_col=None,converters=converter)
print(extraction_table)

#NCBI_matrix = NCBI_table.as_matrix()
#extraction_matrix = extraction_table.as_matrix()

#dictionary can also be written this way: dict() = {}
#sample_map=dict()
#for row in range (extraction_table.shape[0]):
#    sample_map[extraction_table.iat[row,6]]=(extraction_table.iat[row,0], extraction_table.iat[row,1])
#print (sample_map)

#for row in range(NCBI_table.shape[0]):
#        print("\t".join(sample_map[NCBI_table.iat[row,0]]))
