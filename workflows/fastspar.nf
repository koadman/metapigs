#!/usr/bin/env nextflow

/**
 * process a biom table with fastspar
 **/

otu_channel = Channel.fromPath(params.otu_table)
otu_channel2 = Channel.fromPath(params.otu_table)
params.bootstraps = 1000

process fastspar_setup_bootstrap {
    cpus = 2
    memory = '4G'
    container = 'quay.io/biocontainers/fastspar:0.0.10--hd41b482_0'

    input:
	file(otutable) from otu_channel

    output:
    file("otu_data*.tsv") into otu_replicates mode flatten

    script:
"""
fastspar_bootstrap --otu_table ${otutable} --number ${params.bootstraps} --prefix otu_data
"""
}

process fastspar_bootstrap {
	publishDir params.out_dir, mode: 'copy'
    cpus = 2
    memory = '4G'
    container = 'quay.io/biocontainers/fastspar:0.0.10--hd41b482_0'

    input:
    file(otutable) from otu_replicates
    output:
    file("cor_*")
    file("cov_*")

    script:
"""
fastspar --threads 2 --otu_table ${otutable} --correlation cor_${otutable} --covariance cov_${otutable} -i 5
"""
}


