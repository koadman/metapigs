/**
 * This should be run from run_mmseqs.sh
 **/

// The target DB will be incorporated into the channel
db_files = Channel.fromPath(params.targetDB + "*").collect().map{it.toSorted{a,b -> a.name <=> b.name}}
db = file(params.targetDB)

// Build up the read sets from the supplied table, then combine each read-pair set
// with the target DB files. This allows staging of the database for local access
// on execution hosts
read_sets = Channel.fromPath(params.read_table)
		.splitCsv(header: true, sep: '\t', strip: true)
		.map{[file("${params.raw_dir}/${it['r1_filename']}"), file("${params.raw_dir}/${it['r2_filename']}")]}
        // the map organises the row as a list of 3 elements: r1 file, r2 file and then a list of db files
		.combine(db_files).map{it -> [it[0], it[1], it[2..-1]]}


// Carry out the steps to search using mmseqs2, all in one process working in local scratch space
// and returning only the desired files.
process mmseqs_all {
	scratch '/scratch/work/'
	stageInMode 'copy'
	publishDir params.out_dir, mode: 'copy' 

	input:
    // the asterisk allows all the files pertaining to the DB to be staged without explicitly mentioning them all.
	set file(r1), file(r2), file('*') from read_sets

	output:
    file("*.resultdb*") into query_out
	
	script:

	// a local variable to simplify the base code.
	run = r1.baseName

	"""
	mkdir tmp
	cat $r1 $r2 > both.fq.gz
	mmseqs createdb both.fq.gz ${run}.querydb
	mmseqs search ${run}.querydb $db ${run}.resultdb tmp --threads 8
	mmseqs convertalis ${run}.querydb $db ${run}.resultdb ${run}.resultdb.m8 --threads 8
	"""
}
