@Grab('org.codehaus.groovy:groovy-ant') // not part of Nextflow setup
import groovy.util.FileNameFinder
import java.io.File
import java.nio.file.Paths

params.debug = false
params.out_dir = 'checkm_out'

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


process CheckM {

    scratch '/scratch/work/'
    stageInMode 'copy'
    publishDir params.out_dir, mode: 'copy', saveAs: {fn -> "${dir.parent.name}/${fn}" }

    input:
    val dir from genome_dirs

    output:
    set file('checkm_qa.tsv'), file('tree_qa.tsv') into checkm_files

    script:
    if (params.debug) {
        """
        echo $dir > tree_qa.tsv
        echo $dir > checkm_qa.tsv
        """
    }
    else {
        // the same as lineage_wf but with extended output
        """
        checkm tree -t 2 --extension fa $dir checkm_work
        checkm tree_qa --out_format 2 -f tree_qa.tsv --tab_table checkm_work
        checkm lineage_set checkm_work checkm_work/lineage.ms
        checkm analyze -t 2 --extension fa checkm_work/lineage.ms $dir checkm_work
        checkm qa -t 2 --tab_table --out_format 2 -f checkm_qa.tsv checkm_work/lineage.ms checkm_work
        """
    }

}
