#!/usr/bin/env nextflow
params.ncpu = 2
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
 * run phylosift on a read set
 **/
process phylosift {
    cpus params.ncpu
    time '20h'
    memory '48 GB'
    publishDir params.out_dir, mode: 'copy'

    input:
    set r1, r2 from read_sets

    output:
    file("PS_temp/*/*.jplace") into jplace

    """
    ~/software/phylosift_v1.0.1/bin/phylosift all --disable_updates --chunks 1 --chunk_size 100000 --paired ${r1} ${r2}
    """
}
