#!/usr/bin/env python3
import os
import sys
import gzip
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Process a BLAST m8 file to rpkm')
parser.add_argument('--db-file', type=str,
                    help='DB sequence file')
parser.add_argument('--m8-file', type=str,
                    help='DB hits in m8 format')
parser.add_argument('--evalue-threshold', type=float,
                    help='significance threshold for inclusion of hit')
parser.add_argument('--fastq-file', type=str,
                    help='Sequence reads FastQ file')
parser.add_argument('--sample-id', type=str,
                    help='Name to tag results with in tabular output')
args = parser.parse_args()

db_file = open(args.db_file, 'r')
m8file = open(args.m8_file, 'r')
eval_threshold = args.evalue_threshold
fq_file = gzip.open(args.fastq_file)
sample_id = args.sample_id


gene_kbp = {}
for record in SeqIO.parse(args.db_file, "fasta"):
    gene_kbp[record.id]=len(record)/1000

fq_lines = 0
for line in fq_file:
    fq_lines += 1

prev_q = ""
gene_hits = {}
for line in m8file:
    d = line.rstrip().split('\t')
    if d[0] == prev_q:
        continue  # we've already seen this query and it can only be used once
    if float(d[10]) > eval_threshold:
        continue  # hit not strong enough
    if not d[1] in gene_hits:
        gene_hits[d[1]] = 0  # hit a new DB sequence. initialize to 0
    gene_hits[d[1]] += 1.0 / gene_kbp[d[1]]

millions = fq_lines / 4000000

for gene in gene_hits:
    print(sample_id + '\t' + gene + '\t' + str(gene_hits[gene] / millions))
