@Grab('org.codehaus.groovy:groovy-ant') // not part of Nextflow setup
import groovy.util.FileNameFinder
import java.io.File
import java.nio.file.Paths

# --genomes /shared/homes/s1/pig_microbiome/kraken2/paths_to_bins.txt
# --refDB /shared/homes/s1/databases/gtdb_kraken
# --out_dir /shared/homes/s1/pig_microbiome/kraken2/kraken2_out

params.debug = false
params.refDB = 'gtdb_kraken' 
params.threads = 2
params.out_dir = 'kraken2_out' 
params.mem = '150G'

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


process remove_contig_headers {
	
	input:
	val dir from genome_dirs
	
	output:
	bin to bins_channel
	
	script:
		"""
		sed '/^>/d' dir > bin
		"""
}


process run_kraken2 {

    memory = params.mem
    container = 'quay.io/biocontainers/kraken2:2.0.8_beta--pl526h6bb024c_0'

    input:
    val bin from bins_channel

    output:
    file('kraken2_out.tsv') into kraken2_out_files mode flatten

    script:
        """
        kraken2 --db ${params.refDB} $dir --threads ${params.threads} --unclassified-out $dir unclassifed.tsv --classified-out $dir classified.tsv --output $dir kraken2_out.tsv --confidence 0.99 --memory mapping
        """

}

all_krakens = kraken2_out_files.collect()

process concatenate_kraken2_outputs {

	input:
	file('*') from all_krakens
	
	output:
	file('all_kraken2_out.tsv') into all_out 
	
	script:
		"""
		cat file('*') > all_kraken2_out.tsv
		"""
		
