/**
 * Runtime Variables:
 * 
 * r1/r2 = read sets in split files
 * adapter = fasta of adapter sequence
 * phix = fasta of PhiX genome
 * 
 * Expects:
 * 
 * bbduk.sh on the path
 * 
 **/

params.ncpu = 1
params.out_dir = 'out'
params.raw_dir = '.'

/**
 * 
 * Build up the read sets to munge from a CSV run_table
 * 
 **/
Channel.fromPath(params.run_table)
    .splitCsv(header: true, sep: ',', strip: true)
    .map{[it['run_id'],  file("${params.raw_dir}/${it['r1_filename']}"), file("${params.raw_dir}/${it['r2_filename']}"), it['source_id']]}
    .into{test_exists; read_sets}

process TestExistence {
        input:
        set run_id, r1, r2, source_id from test_exists
        
        exec:
        assert r1.exists(), "The file $r1 did not exist"
        assert r2.exists(), "The file $r2 did not exist"
}

                
/**
 * Clean up a readset using bbduk from BBTools
 * 
 * A three stage process:
 *
 *  1. removal of adapter sequences
 *  2. quality trimming
 *  3. removal of any PhiX contamination
 * 
 * Cleaned reads, those reads removed and statistics for each readset are published into per source directories
 **/

process CleanUp {
    cpus params.ncpu
    publishDir params.out_dir, mode: 'copy', saveAs: {fn -> "${source_id}/reads/${fn}" }
    
    input:
    set run_id, r1, r2, source_id from read_sets
    each file(adapters) from Channel.fromPath(params.adapters)
    each file(phix) from Channel.fromPath(params.phix)
    
    output:
    set source_id, run_id, file("${run_id}_cleaned_paired.fq.gz"), file('*matched.fq.gz'), file('*_stats.txt') into cleaned_reads
    
    """
    bbduk.sh t=${task.cpus} k=23 hdist=1 tpe tbo mink=11 ktrim=r ref=$adapters \
        in=$r1 in2=$r2 out=stdout.fq outm=${run_id}_adapter_matched.fq.gz stats=${run_id}_adapter_stats.txt | \
    bbduk.sh t=${task.cpus} ftm=0 qtrim=r trimq=10 \
        in=stdin.fq out=stdout.fq stats=${run_id}.quality_stats.txt |
    bbduk.sh t=${task.cpus} k=31 hdist=1 ref=$phix \
        in=stdin.fq out=${run_id}_cleaned_paired.fq.gz outm=${run_id}_phix_matched.fq.gz stats=${run_id}_phix_stats.txt
    """
}

cleaned_reads = cleaned_reads.map{it -> [it[0], it[1], it[2]]}.groupTuple()

/**
 * 
 * Assemble a pooled data set, where pooling was by source id
 * 
 * The resulting output will be published to out/source_id/asm
 *
 **/

process Assembly {
    cpus params.ncpu
    publishDir params.out_dir, mode: 'copy', saveAs: {fn -> "${source_id}/asm/${fn}" }
    
    input:
    set source_id, run_id, reads from cleaned_reads

    output:
    set file("megahit_out/${source_id}.contigs.fa"), file("megahit_out/${source_id}.log"), file("megahit_out/opts.txt") into assembly
    
    """
    megahit -t ${task.cpus} -m 0.33 -o megahit_out --out-prefix $source_id --12 ${reads.join(",")}
    """
}
