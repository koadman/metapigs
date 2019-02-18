#!/usr/bin/env python3
import os
import sys
import gzip
import argparse

parser = argparse.ArgumentParser(description='Process a BLAST m8 file to rpkm')
parser.add_argument('query_file', type=str, nargs='+',
                    help='query sequence file')
parser.add_argument('m8_file', type=str, nargs='+',
                    help='DB hits in m8 format')
parser.add_argument('evalue_threshold', type=float, nargs='+',
                    help='significance threshold for inclusion of hit')
parser.add_argument('fastq_file', type=str, nargs='+',
                    help='Sequence reads FastQ file')
parser.add_argument('sample_id', type=str, nargs='+',
                    help='Name to tag results with in tabular output')
args = parser.parse_args()

query_file = open(args.query_file, 'r')
m8file = open(args.m8_file, 'r')
eval_threshold = args.evalue_threshold
fq_file = gzip.open(args.fastq_file)
sample_id = args.sample_id

fq_lines = 0
for line in fq_file:
    fq_lines += 1

prev_q = ""
gene_hits = {}
gene_kbp = {}
for line in m8file:
    d = line.rstrip().split('\t')
    if d[0] == prev_q:
        continue  # we've already seen this query and it can only be used once
    if float(d[10]) > eval_threshold:
        continue  # hit not strong enough
    if not d[1] in gene_hits:
        gene_hits[d[1]] = 0  # hit a new DB sequence. initialize to 0
        gene_kbp[d[1]] = 1.1  # HACK - this needs to get read in from the DB
    gene_hits[d[1]] += 1.0 / gene_kbp[d[1]]

millions = fq_lines / 4000000

for gene in gene_hits:
    print(sample_id + '\t' + gene + '\t' + str(gene_hits[gene] / millions))
