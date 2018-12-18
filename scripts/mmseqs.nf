#!/usr/bin/env nextflow

/**
activate conda environment with mmseqs and nextflow
nextflow run mmseqs.nf 
--targetDB = aminoglycosideDB
--out_dir = /shared/homes/12705859/mmseqs_nextflow
--run_table = /shared/homes/12705859/mmseqs_nextflow/reads.tsv 
--raw_dir = /shared/homes/s1/pig_microbiome/MON5838
**/


/**
 * usage: mmseqs.nf --run_table=table.tsv
 **/
params.ncpu = 8
params.out_dir = 'out'
params.raw_dir = '.'
BNAME='basename{r1}'

/**
 *
 * Build up the read sets to munge from a CSV run_table
 *
 **/
Channel.fromPath(params.run_table)
    .splitCsv(header: true, sep: '\t', strip: true)
    .map{[file("${params.raw_dir}/${it['r1_filename']}"), file("${params.raw_dir}/${it['r2_filename']}")]}
    .into{test_exists; read_sets}

process TestExistence {
        input:
        set r1, r2 from test_exists

        exec:
        assert r1.exists(), "The file $r1 did not exist"
        assert r2.exists(), "The file $r2 did not exist"
}

/**
 * compute mmseqs
 **/
 process mmseqs_create_query_search_convertalis {
  publishDir params.out_dir, mode: 'copy'

  input:
  set r1, r2 from read_sets

  output:
  file("*.queryDB") 
  file("*resultDB.m8") 

  """
  cat ${r1} ${r2} > both.fq.gz
  mmseqs createdb both.fq.gz ${r1}.queryDB
  mmseqs search $BNAME.queryDB ${params.targetDB} $BNAME.resultDB --threads 8 
  mmseqs convertalis $BNAME.queryDB ${params.targetDB} $BNAME.resultDB $BNAME.resultDB.m8
  """
}