#!/usr/bin/env python3
import os
import sys
import gzip

m8file = open(sys.argv[1],'r')
eval_threshold=float(sys.argv[2])
fq_file=gzip.open(sys.argv[3])
sample_id=sys.argv[4]

fq_lines=0
for line in fq_file:
	fq_lines += 1

prev_q = ""
gene_hits={}
gene_kbp={}
for line in m8file:
	d = line.rstrip().split('\t')
	if d[0]==prev_q:
		continue	# we've already seen this query and it can only be used once
	if float(d[10])>eval_threshold:
		continue	# hit not strong enough
	if not d[1] in gene_hits:
		gene_hits[d[1]]=0	# hit a new DB sequence. initialize to 0
		gene_kbp[d[1]]=1.1 # HACK - this needs to get read in from the DB
	gene_hits[d[1]]+= 1.0/gene_kbp[d[1]]

millions = fq_lines/4000000

for gene in gene_hits:
	print(sample_id+'\t'+gene+'\t'+str(gene_hits[gene]/millions))
