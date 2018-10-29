#!/usr/bin/env nextflow
/**
 * usage: phylosift.nf --run_table=table.tsv
 **/
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
    phylosift all --disable_updates --chunks 1 --chunk_size 100000 --paired ${r1} ${r2}
    """
}

jplace_all = jplace.collect()
(jplace_alpha, jplace_squash, jplace_beta) = jplace_all.into(3)

/**
 * compute per-sample alpha diversity with various metrics
 **/
process alpha_diversity {
  publishDir params.out_dir, mode: 'copy'

  input:
  file("*.jplace") from jplace_alpha
  file("metadata.tsv") from Channel.fromPath(params.run_table)

  output:
  file("*.alphadiv") into alphadiv

  """
  guppy fpd *.jplace > all.alphadiv
  """
}

/**
 * cluster the samples
 **/
process squash_clust {
  publishDir params.out_dir, mode: 'copy'

  input:
  file("*.jplace") from jplace_squash
  file("metadata.tsv") from Channel.fromPath(params.run_table)

  output:
  file("*.clust") into squash

  """
  guppy squash *.jplace > all.clust
  """
}

/**
 * edge PCA to explore variation in community composition among samples
 **/
process beta_diversity {
  publishDir params.out_dir, mode: 'copy'

  input:
  file("*.jplace") from jplace_beta
  file("metadata.tsv") from Channel.fromPath(params.run_table)

  output:
  file("pca*") into betadiv

  """
  guppy epca --prefix pca *.jplace
  """
}
