# versions: 
#bbmap                     38.22
#bwa                       0.7.17
#megahit                   1.1.3
#metabat2                  2.12.1
#samtools                  1.9

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

params.debug = false
params.seed = 1234

/**
 * Build up the read sets to munge from a CSV run_table
 **/
test_exists = Channel.fromPath(params.run_table)
    .splitCsv(header: true, sep: '\t', strip: true)
    .map{[it['alt_sample_name'],  file(it['r1_filename']),
          file(it['r2_filename']), it['isolation_source']]}

/**
 * Test that each paired readset mentioned in the run table exists. The process will stop
 * the workflow of a either part of any pair is missing.
 *
 * The output of this process is the same as the input. Doing this, rather than cloning
 * the channel above, makes testing existence a prerequisite of the clean-up process. Therefore,
 * the existence of all readsets must be established before attempting any cleanups.
 */
process TestExistence {

    // directives
    cpus 1
    conda false

    input:
    set run_id, r1, r2, source_id from test_exists

    output:
    set run_id, r1, r2, source_id into valid_sets

    exec:
    if (params.debug) {
        "everything is fine Dave and all my circuits are functioning perfectly"
    }
    else {
        assert r1.exists(), "The file $r1 did not exist"
        assert r2.exists(), "The file $r2 did not exist"
        }
}


/**
 * Clean up a readset using bbduk from BBTools and publish the results
 *
 * A three stage process:
 *
 *  1. removal of adapter sequences
 *  2. quality trimming
 *  3. removal of any PhiX contamination
 *
 * Cleaned reads, those reads removed and statistics for each readset are published into per source directories. The
 * process uses pipes to avoid creating many intermediate files, which would become a significant resource requirement
 * and IO load.
 *
 * Reads which match either adapter or PhiX are written to files and kept. Statistics on removal of reads and base
 * trimming are kept in the ".stats.txt" files.
 *
 * The final cleaned readset is interleaved, as are the matched "bad" read files.
 **/

process CleanUp {

    // directives
    scratch '/scratch/work/'
    stageInMode 'copy'
    publishDir params.out_dir, mode: 'symlink', saveAs: {fn -> "${source_id}/reads/${fn}" }

    memory { 30.GB * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
    set run_id, r1, r2, source_id from valid_sets
    each file(adapters) from Channel.fromPath(params.adapters)
    each file(phix) from Channel.fromPath(params.phix)

    output:
    set source_id, run_id, file("${run_id}_cleaned_paired.fq.gz"), file('*_matched.fq.gz'), file('*_stats.txt') into cleaned_reads

    script:
    if (params.debug) {
        """
        echo "$source_id $run_id" > ${run_id}_cleaned_paired.fq.gz
        echo "$source_id $run_id" > ${run_id}_cleaned_paired.fq.gz
        echo "$source_id $run_id" > ${run_id}_adapters_matched.fq.gz
        echo "$source_id $run_id" > ${run_id}_adapters_stats.txt
        echo "$source_id $run_id" > ${run_id}_phix_matched.fq.gz
        echo "$source_id $run_id" > ${run_id}_phix_stats.txt
        echo "$source_id $run_id" > ${run_id}_quality_stats.txt
        """
    }
    else {
        """
        bbduk.sh t=1 k=23 hdist=1 tpe tbo mink=11 ktrim=r ref=$adapters \
            int=t in=$r1 in2=$r2 out=stdout.fq outm=${run_id}_adapter_matched.fq.gz stats=${run_id}_adapter_stats.txt | \
        bbduk.sh t=1 ftm=0 qtrim=r trimq=20 \
            int=t in=stdin.fq out=stdout.fq stats=${run_id}_quality_stats.txt |
        bbduk.sh t=2 k=31 hdist=1 ref=$phix \
            int=t in=stdin.fq out=${run_id}_cleaned_paired.fq.gz outm=${run_id}_phix_matched.fq.gz stats=${run_id}_phix_stats.txt
        """
    }
}

// as we need to consume it twice, create a temporary copy of the cleaned_reads
(cleaned_reads, cleaned_reads_copy) = cleaned_reads.into(2)

// pooled assemblies require reads grouped by source "pig"
asm_tasks = cleaned_reads.map{it -> [it[0], it[1], it[2]]}.groupTuple()

/**
 *
 * Assemble a pooled data set, where pooling was by source id
 *
 * The resulting output will be published to out/source_id/asm
 *
 **/

process PooledAssembly {

    // directives
    cpus 16
    memory { 80.GB * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
    maxRetries 3

    scratch '/scratch/work/'
    stageInMode 'copy'
    publishDir params.out_dir, mode: 'symlink', saveAs: {fn -> "${source_id}/asm/${fn}" }

    input:
    set source_id, run_id, reads from asm_tasks

    output:
    set source_id, file("megahit_out/${source_id}.contigs.fa"), file("megahit_out/${source_id}.log"), file("megahit_out/opts.txt") into assemblies

    script:
    if (params.debug) {
        """
        mkdir megahit_out
        echo "$source_id $reads" > megahit_out/${source_id}.contigs.fa
        echo "$source_id $reads" > megahit_out/${source_id}.log
        echo "$source_id $reads" > megahit_out/opts.txt
        """
    }
    else {
        """
        megahit -t ${task.cpus} --memory ${task.memory.size} --min-contig-len 5 --prune-level 3 --max-tip-len 280 -o megahit_out --out-prefix $source_id --12 ${reads.join(",")}
        """
    }
}


index_tasks = assemblies.map{it[0..1]}

process CreateContigsIndexes {

    cpus 1
    scratch '/scratch/work/'
    stageInMode 'copy'
    publishDir params.out_dir, mode: 'symlink', saveAs: {fn -> "${source_id}/asm/${fn}" }

    input:
    set source_id, file('contigs.fa') from index_tasks

    output:
    set source_id, file("contigs.fa"), file("contigs.fa.*") into indexed_assemblies

    script:
    if (params.debug) {
        """
        echo "$source_id contigs" > contigs.fa.bwt
        echo "$source_id contigs" > contigs.fa.pac
        echo "$source_id contigs" > contigs.fa.ann
        echo "$source_id contigs" > contigs.fa.amb
        echo "$source_id contigs" > contigs.fa.sa
        """
    }
    else {
        """
        bwa index contigs.fa
        """
    }
}

// create a channel which combines each readset with the relevant pooled assembly
// unwrap the lists and take only relevant variables
(map_tasks, map_tasks_copy) = indexed_assemblies.cross(cleaned_reads_copy).map{it[0][0..1] + it[1][1..2]}.into(2)

process MapReadsToPooled {

    memory { 8.GB * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
    maxRetries 3

    scratch '/scratch/work/'
    stageInMode 'copy'
    publishDir params.out_dir, mode: 'symlink', saveAs: {fn -> "${source_id}/mapped/${fn}" }

    input:
    set source_id, contigs, run_id, reads from map_tasks

    output:
    set source_id, contigs, file("*.bam"), file("*.bam.bai") into bam_files

    script:
    if (params.debug) {
        """
        echo "$source_id $run_id" > ${source_id}_${run_id}.bam
        echo "$source_id $run_id" > ${source_id}_${run_id}.bam.bai
        """
    }
    else {
        """
        bwa mem -p -t${task.cpus} $contigs $reads | samtools view -@${task.cpus} -uS | samtools sort -@${task.cpus} -o ${source_id}_${run_id}.bam -
        samtools index -@${task.cpus} ${source_id}_${run_id}.bam
        """
    }
}

timeseries_bams = bam_files.groupTuple().map{it ->
    source_id = it[0]
    contigs = it[1][0]
    bams = it[2].toSorted{a, b -> a.name <=> b.name}
    [source_id, contigs, bams]}

process MetaBat2 {

    cpus 16
    memory { 80.GB * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
    maxRetries 3

    scratch '/scratch/work/'
    stageInMode 'copy'
    publishDir params.out_dir, mode: 'symlink', saveAs: {fn -> "${source_id}/metabat/${fn}" }

    input:
    set source_id, contigs, bams from timeseries_bams

    output:
    set source_id, file("depth.txt"), file("bins*")

    script:
    if (params.debug) {
        """
        echo "jgi_summarize_bam_contig_depths --outputDepth depth.txt ${bams.join(' ')}" > depth.txt
        echo "metabat2 --seed $params.seed -t ${task.cpus} -i $contigs -a depth.txt -o metabat.bins --minContig 2500" > bins
        """
    }
    else {
        """
        jgi_summarize_bam_contig_depths --outputDepth depth.txt ${bams.join(' ')}
        metabat2 --seed $params.seed -t ${task.cpus} -i $contigs -a depth.txt -o bins --minContig 2500 --saveCls
        """
    }
}
