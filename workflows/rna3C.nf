#!/usr/bin/env nextflow

params.out_dir = 'output'
params.cpus = 4
params.mem = '8G'

/**
 * Build up the read sets to munge from a CSV run_table
 **/
read_sets = Channel.fromPath(params.run_table)
    .splitCsv(header: true, sep: '\t', strip: true)
    .map{[file("${params.raw_dir}/${it['r1_filename']}"), file("${params.raw_dir}/${it['r2_filename']}")]}

process sortmerna {
	publishDir params.out_dir, mode: 'copy' 
    cpus = params.cpus
    memory = params.mem
    container = 'quay.io/biocontainers/sortmerna:2.1b--he860b03_4'

    input:
    set r1, r2 from read_sets
    output:
    file("*.rrna.fq") into rna_reads

    script:
"""
zcat ${r1} ${r2} > read_pairs.fq
BNAME=`basename ${r1} _R1_001.fastq.gz`
DB=${params.db_path}
RNADB=\$DB/rfam-5.8s-database-id98.fasta,\$DB/rfam-5.8s-database-id98.idx:\$DB/rfam-5s-database-id98.fasta,\$DB/rfam-5s-database-id98.idx:\$DB/silva-arc-16s-id95.fasta,\$DB/silva-arc-16s-id95.idx:\$DB/silva-arc-23s-id98.fasta,\$DB/silva-arc-23s-id98.idx:\$DB/silva-bac-16s-id90.fasta,\$DB/silva-bac-16s-id90.idx:\$DB/silva-bac-23s-id98.fasta,\$DB/silva-bac-23s-id98.idx:\$DB/silva-euk-18s-id95.fasta,\$DB/silva-euk-18s-id95.idx:\$DB/silva-euk-28s-id98.fasta,\$DB/silva-euk-28s-id98.idx
sortmerna --reads read_pairs.fq --ref \$RNADB --fastx --aligned \$BNAME.rrna -v -a ${params.cpus}
rm -f read_pairs.fq.gz
"""
}

