#!/usr/bin/env nextflow

/**
 * process Hi-C timeseries data
 **/

// Build up the read sets from the supplied table, then combine each read-pair set
// with the target DB files. This allows staging of the database for local access
// on execution hosts

hiC_reads = Channel.fromPath(params.hic_table)
		.splitCsv(header: true, sep: '\t', strip: true)
		.map{[file("${params.raw_dir}/${it['hic_r1_filename']}"), 
            file("${params.raw_dir}/${it['hic_r2_filename']}"),
            file("${it['asm_filename']}")]}
        // the map organises the row as a list of 3 elements: r1 file, r2 file and then a list of db files
		.map{it -> [it[0], it[1], it[2]]}


process clean_hic {
	publishDir params.out_dir, mode: 'copy' 
    cpus = 3
    memory = '8G'
    container = 'quay.io/biocontainers/bbmap:38.68--h516909a_0'

    input:
	set file(r1), file(r2), file(assembly) from hiC_reads
    output:
    set file("cleaned.*"), file(assembly) into cleanedreads

    script:
"""
    bbduk.sh threads=2 -Xmx7g ref=/usr/local/opt/bbmap-38.68-0/resources/adapters.fa minlength=80 qtrim=r trimq=10 in=${r1} in2=${r2} out=cleaned.${r1} -minavgquality=15 2> ${r1}.bbduk.metrics

"""
}

process map_hic {
    memory = '12G'
    cpus = 8
    container = 'quay.io/biocontainers/bwa:0.7.17--hed695b0_6'

    input:
    set file(cleaned), file(assembly) from cleanedreads
    output:
    set file("${cleaned}.hic.sam"), file(assembly)  into mappedreads

    script:
"""
    bwa index ${assembly}
    bwa mem -t 8 -5 -S -P -p ${assembly} ${cleaned} > ${cleaned}.hic.sam
"""
}

process bamsort_hic {
    memory = '20G'
    cpus = 8
    container = 'quay.io/biocontainers/samtools:1.9--h10a08f8_12'

    input:
    set file(sam), file(assembly)  from mappedreads
    output:
    set file("*.bam"), file("*.bai"), file(assembly) into bam

    script:
"""
    BNAME=`basename ${sam} .hic.sam`
    samtools view -S -b *.hic.sam | samtools sort -@ 6 -m 2G - > \$BNAME.hic.bam
    samtools index \$BNAME.hic.bam
"""
}

process markdup_hic {
	publishDir params.out_dir, mode: 'copy' 
	stageInMode 'copy'
    memory = '12G'
    cpus = 2
    container = 'quay.io/biocontainers/gatk4:4.1.3.0--0'

    input:
    set file(bam), file(bai), file(assembly) from bam
    output:
    set file("*.hic.dm.bam"), file(bai), file(assembly) into dupmarked
    file("*.dupmark_metrics.txt") into dm_metrics

    script:
"""
    BNAME=`basename ${bam} .hic.bam`
    gatk MarkDuplicates -I \$BNAME.hic.bam -O \$BNAME.hic.dm.bam -M \$BNAME.dupmark_metrics.txt
"""

}

