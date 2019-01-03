#!/usr/bin/env nextflow

/**
activate conda environment with mmseqs and nextflow
nextflow run mmseqs.nf 
--targetDB = /shared/homes/12705859/mmseqs_nextflow/aminoglycosideDB
--out_dir = /shared/homes/12705859/mmseqs_nextflow
--run_table = /shared/homes/12705859/mmseqs_nextflow/reads.tsv 
--raw_dir = /shared/homes/s1/pig_microbiome/MON5838
**/


/**
 * usage: mmseqs.nf --run_table=table.tsv
 **/
params.ncpu = 4
params.out_dir = 'out'
params.raw_dir = '.'

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
 process mmseqs {
  publishDir params.out_dir, mode: 'copy' 

  cpus 4

  input:
  set r1, r2 from read_sets

  output:
  file('*.queryDB') 
  file('*.dbtype') 
  file('*.index')
  file('*.lookup')
  file('*.resultDB')  
  file('*.resultDB.m8')
   

  """
  cat ${r1} ${r2} > both.fq.gz
  mmseqs createdb both.fq.gz  ${r1.baseName}.queryDB
  mmseqs search ${r1.baseName}.queryDB ${params.targetDB} ${r1.baseName}.resultDB tmp.${r1.baseName} --threads 4 
  mmseqs convertalis ${r1.baseName}.queryDB ${params.targetDB} ${r1.baseName}.resultDB ${r1.baseName}.resultDB.m8
  """
}
