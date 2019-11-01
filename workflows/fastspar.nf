#!/usr/bin/env nextflow

/**
 * process a biom table with fastspar
 **/

params.bootstraps = 1000
params.out_dir = "output"
params.mem = '4G'
params.cpus = 2

process fastspar_setup_bootstrap {
    cpus = params.cpus
    memory = params.mem
    container = 'quay.io/biocontainers/fastspar:0.0.10--hd41b482_0'

    input:
	file(otutable) from Channel.fromPath(params.otu_table)

    output:
    file("otu_data*.tsv") into otu_replicates mode flatten

    script:
"""
fastspar_bootstrap --otu_table ${otutable} --number ${params.bootstraps} --prefix otu_data
"""
}

process fastspar_bootstrap {
    cpus = params.cpus
    memory = params.mem
    container = 'quay.io/biocontainers/fastspar:0.0.10--hd41b482_0'

    input:
    file(otutable) from otu_replicates
    output:
    file("cor_*") into correlations
    file("cov_*")

    script:
"""
export OMP_NUM_THREADS=${params.cpus}
fastspar -y --threads ${params.cpus} --otu_table ${otutable} --correlation cor_${otutable} --covariance cov_${otutable} -i 5
"""
}

all_cors = correlations.collect()

process fastspar {
	publishDir params.out_dir, mode: 'copy'
    cpus = params.cpus
    memory = params.mem
    container = 'quay.io/biocontainers/fastspar:0.0.10--hd41b482_0'

    input:
	file(otutable) from Channel.fromPath(params.otu_table)
    output:
    file('median_covariance.tsv')
    file('median_correlation.tsv') into median_cor

    script:
"""
export OMP_NUM_THREADS=${params.cpus}
fastspar -y --threads ${params.cpus} --otu_table ${otutable}  --correlation median_correlation.tsv --covariance median_covariance.tsv 
"""
}

process fastspar_pvalue {
	publishDir params.out_dir, mode: 'copy'
    cpus = params.cpus
    memory = params.mem
    container = 'quay.io/biocontainers/fastspar:0.0.10--hd41b482_0'

    input:
	file(otutable) from Channel.fromPath(params.otu_table)
    file('*') from all_cors
    file('median_correlation.tsv') from median_cor
    output:
    file('pvalues.tsv')
    file('median_correlation.tsv')

    script:
"""
fastspar_pvalues --otu_table ${otutable} --correlation median_correlation.tsv  --prefix cor_ --permutations ${params.bootstraps} -o pvalues.tsv 
"""
}
