@Grab('org.codehaus.groovy:groovy-ant') // not part of Nextflow setup
import groovy.util.FileNameFinder
import java.io.File
import java.nio.file.Paths

params.refDB = 'gtdb_kraken' # '/shared/homes/s1/databases/gtdb_kraken'
params.threads = 2
params.debug = false
params.out_dir = 'kraken2_out' # '/shared/homes/s1/pig_microbiome/kraken2/kraken2_out'
params.mem = '4G'

fnf = new FileNameFinder()

/**
 * Read a list of folders each pertaining to a metagenome and containing the binnned fasta files
 *
 * Each non-empty line of the file is converted to a path object
 */
genome_dirs = Channel.fromPath(params.genomes)
        .splitText{it.trim()}       // remove leading and trailing whitespace/newline
        .filter{s -> !s.isEmpty()}  // remove empty strings
        .map{s -> file(s)}          // convert to file objects
        .filter{ p ->               // remove any entries that are not directories or do not contain fasta files
            numFasta = fnf.getFileNames(p.toString(), '*.fa').size()
            p.isDirectory() && numFasta > 0
        }


process run_kraken2 {

    scratch '/scratch/work/'
    stageInMode 'copy'
    publishDir params.out_dir, mode: 'copy', saveAs: {fn -> "${dir.parent.name}/${fn}" }
    memory = params.mem
    container = 'quay.io/biocontainers/kraken2:2.0.8_beta--pl526h6bb024c_0'

    input:
    val dir from genome_dirs

    output:
    file('unclassified.tsv'), file('classified.tsv'), file('kraken2_out.tsv') into kraken2_out_files mode flatten

    script:
    if (params.debug) {
        """
        echo $dir > unclassified.tsv
        echo $dir > classified.tsv
        """
    }
    else {
        """
        kraken2 --db ${params.refDB} $dir --threads ${params.threads} --unclassified-out $dir unclassifed.tsv --classified-out $dir classified.tsv --output $dir kraken2_out.tsv --confidence 0.99 --memory mapping
        """
    }

}
