/**
Runtime Variables:

r1/r2 = read sets in split files
adapter = fasta of adapter sequence
phix = fasta of PhiX genome

Expects:

bbduk.sh on the path

**/

params.ncpu = 1
params.out_dir = 'out'
params.raw_dir = '.'

read_sets = Channel.fromPath(params.run_table)
                .splitCsv(header: true, sep: ',', strip: true)
                .map{[it['run_id'],  file("${params.raw_dir}/${it['r1_filename']}"), file("${params.raw_dir}/${it['r2_filename']}"), it['source_id']]}


/**
 Clean Up a readset using bbduk from BBTools
 
 A three stage process, removing first adapter sequences, 
 then quality trimmed and finally PhiX contamination.
 
 Removed reads and statistics for each readset are also published into per source directories
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

cleaned_reads = cleaned_reads.map{it -> [it[0], it[1], it[2]]}.groupTuple() //.subscribe{println it}

process Assembly {
    cpus params.ncpu
    publishDir params.out_dir, mode: 'copy', saveAs: {fn -> "${source_id}/asm/${fn}" }
    
    input:
    set source_id, run_id, reads from cleaned_reads

    output:
    set file("megahit_out/${source_id}.contigs.fa"), file("megahit_out/${source_id}.log"), file("megahit_out/opts.txt") into assembly
    
    """
    megahit -t ${task.cpus} -o megahit_out --out-prefix $source_id --12 ${reads.join(",")}
    """
}

/**
 A fine-grained version of clean-up, which probably costs too much in disk storage to be of value.
 **/
 
/*
process RemoveAdapters {
    cpus params.ncpu

    input:
    set run_id, r1, r2, source_id from read_sets
    each file(adapters) from Channel.fromPath(params.adapters)
    
    output:
    set source_id, run_id, file('stage1.fq.gz') into stage1

    """
    bbduk.sh t=${task.cpus} k=23 hdist=1 tpe tbo mink=11 ktrim=r ref=$adapters \
        in=$r1 in2=$r2 out=stage1.fq.gz outm=adapter_matched.fq.gz stats=stage1_adapter.stats
    """
}

process QualityTrim {
    cpus params.ncpu

    input:
    set source_id, run_id, reads from stage1
    
    output:
    set source_id, run_id, file('stage2.fq.gz') into stage2

    """
    bbduk.sh t=${task.cpus} ftm=0 qtrim=r trimq=10 \
        in=$reads out=stage2.fq.gz stats=stage2_quality.stats 
    """
}

process RemovePhiX {
    cpus params.ncpu

    input:
    set source_id, run_id, reads from stage2
    each file(phix) from Channel.fromPath(params.phix)
    
    output:
    set source_id, run_id, file('cleaned_*R1.fq.gz'), file('cleaned_*R2.fq.gz') into cleaned_reads

    """
    bbduk.sh t=${task.cpus} k=31 hdist=1 ref=$phix \
        in=$reads out=cleaned_${run_id}_R1.fq.gz out2=cleaned_${run_id}_R2.fq.gz outm=phix_matched.fq.gz stats=stage3_phix.stats
    """
}*/

