#!/usr/bin/env nextflow

/**
	nextflow run mmseqs.nf 
		--targetDB = /shared/homes/s1/pig_microbiome/mmseqs_nextflow/aminoglycosideDB
		--out_dir = /shared/homes/s1/pig_microbiome/mmseqs_nextflow/
		--run_table = /shared/homes/s1/pig_microbiome/mmseqs_nextflow/reads.tsv 
		--raw_dir = /shared/homes/s1/pig_microbiome/MON5838
**/

/**
 *
 * Build up the read sets to munge from a CSV run_table
 *
 **/
read_sets = Channel.fromPath(params.run_table)
    .splitCsv(header: true, sep: '\t', strip: true)
    .map{[file("${params.raw_dir}/${it['r1_filename']}"), file("${params.raw_dir}/${it['r2_filename']}")]}


target_db = Channel.fromPath(params.targetDB + "*").collect().map{it.toSorted{a,b -> a.name <=> b.name}}

read_sets = read_sets.combine(target_db).map{it -> [it[0], it[1], it[2..-1]]} //.subscribe{println it}


process mmseqs_all {
	cpus 8
	scratch '/scratch/work/'
	stageInMode 'copy'
	publishDir params.out_dir, mode: 'copy' 

	input:
	set file(r1), file(r2), file('*') from read_sets

	output:
	set file("*.querydb*"), file("*.resultdb*") into query_out
	

	script:
	run = r1.baseName

	"""
	mkdir tmp
	cat $r1 $r2 > both.fq.gz
	mmseqs createdb both.fq.gz ${run}.querydb
	mmseqs search ${run}.querydb agly.db ${run}.resultdb tmp --threads 8
	mmseqs convertalis ${run}.querydb agly.db ${run}.resultdb ${run}.resultdb.m8 --threads 8
	"""
}
