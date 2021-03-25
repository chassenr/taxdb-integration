rule fix_ncbi_taxpath:
	input:
		gtdb = config["rdir"] + "/tax_combined/gtdb_derep_taxonomy.txt",
		checkv = config["rdir"] + "/tax_combined/checkv_derep_taxonomy.txt",
		ncbi = expand(config["rdir"] + "/tax_combined/{library_highres}_derep_taxonomy.txt", library_highres = LIBRARY_HIGHRES)
	output:
		tax_combined = config["rdir"] + "/tax_combined/derep_taxonomy_combined.txt"
	params:
		outdir = config["rdir"] + "/tax_combined",
		script = config["wdir"] + "/scripts/fix_ncbi_taxpath.R"
	conda:
		config["wdir"] + "/envs/r.yaml"
	log:
		config["rdir"] + "/logs/fix_ncbi_taxpath.log"
	shell:
		"""
		cat {input} > "{params.outdir}/tmp"
		{params.script} -i "{params.outdir}/tmp" -o "{output.tax_combined}" &>> {log}
		rm "{params.outdir}/tmp"
		"""

rule format_taxonomy:
	input:
		tax_combined = config["rdir"] + "/tax_combined/derep_taxonomy_combined.txt"
	output:
		tax_good = config["rdir"] + "/tax_combined/derep_taxonomy_good.txt",
		nodes = config["rdir"] + "/kraken2_db/taxonomy/nodes.dmp",
		names = config["rdir"] + "/kraken2_db/taxonomy/names.dmp",
		file_list = config["rdir"] + "/kraken2_genomes/file_names_derep_genomes.txt"
	params:
		krakendir = config["rdir"] + "/kraken2_genomes/genome_files",
		genomedir = config["rdir"] + "/derep_combined",
		tax_script = config["tax_script"]
	conda:
		config["wdir"] + "/envs/biopython.yaml"
	log:
		config["rdir"] + "/logs/format_taxonomy.log"
	shell:
		"""
		cut -f1,2 {input.tax_combined} > {output.tax_good}
		{params.tax_script} --gtdb {output.tax_good} --assemblies {params.genomedir} --nodes {output.nodes} --names {output.names} --kraken_dir {params.krakendir} &>> {log}
		# replace 'domain' with 'superkingdom (required for kaiju) to ensure that the same nodes.dmp and names.dmp files can be used for both databases
		sed -i -e 's/domain/superkingdom/g' {output.nodes}
		find {params.krakendir} -type f -name '*.fa' > {output.file_list}
		"""
# depending on the number and size of the genomes, it may be required to delete the zipped fna in the derep directory

# optional: run conterminator
if config["kingdoms_highres"]:
	rule cat_library:
		input:
			fasta_ncbi = expand(config["rdir"] + "/kraken2_db/tmp/{library_highres}_library.fna", library_highres = LIBRARY_HIGHRES),
			fasta_gtdb = config["rdir"] + "/kraken2_db/tmp/gtdb_library.fna",
			fasta_checkv = config["rdir"] + "/kraken2_db/library/checkv/library.fna",
			map_ncbi = expand(config["rdir"] + "/kraken2_db/tmp/{library_highres}_prelim_map.txt", library_highres = LIBRARY_HIGHRES),
			map_gtdb = config["rdir"] + "/kraken2_db/tmp/gtdb_prelim_map.txt",
			map_checkv = config["rdir"] + "/kraken2_db/library/checkv/prelim_map.txt"
		output:
			tmp_fna = config["rdir"] + "/kraken2_db/tmp/library.fna",
			tmp_map = config["rdir"] + "/kraken2_db/tmp/prelim_map.txt"
		shell:
			"""
			cat {input.fasta_ncbi} {input.fasta_gtdb} {input.fasta_checkv} > {output.tmp_fna}
			cat {input.map_ncbi} {input.map_gtdb} {input.map_checkv} > {output.tmp_map}
			"""

	localrules: detect_contam_highres
	
	rule detect_contam_highres:
		input:
			tmp_fna = config["rdir"] + "/kraken2_db/tmp/library.fna",
			tmp_map = config["rdir"] + "/kraken2_db/tmp/prelim_map.txt",
			nodes = config["rdir"] + "/kraken2_db/taxonomy/nodes.dmp",
			names = config["rdir"] + "/kraken2_db/taxonomy/names.dmp"
		output:
			delnodes = config["rdir"] + "/kraken2_db/taxonomy/delnodes.dmp",
			merged = config["rdir"] + "/kraken2_db/taxonomy/merged.dmp",
			kstring = config["rdir"] + "/decontamination/conterminator_string.txt",
			xstring = config["rdir"] + "/decontamination/conterminator_blacklist.txt",
			cmap = config["rdir"] + "/decontamination/cmap.txt",
			contam = config["rdir"] + "/decontamination/highres_db_conterm_prediction"
		params:
			script = config["wdir"] + "/scripts/get_kingdoms_conterminator.R",
			tmpdir = config["rdir"] + "/decontamination/tmp",
			taxdir = config["rdir"] + "/kraken2_db/taxonomy/",
			prefix = config["rdir"] + "/decontamination/highres_db",
			cmem = config["cmem"],
			kingdoms = config["kingdoms_highres"]
		conda:
			config["wdir"] + "/envs/r.yaml"
		threads: config["masking_threads"]
		log:
			config["rdir"] + "/logs/highres_conterminator.log"
		shell:
			"""
			# prepare fasta header mapping file for conterminator
			cut -f2,3 {input.tmp_map} > {output.cmap}
			# create dummy delnodes and merged files in the taxonomic directory for compatibility with conterminator
			touch {output.delnodes}
			touch {output.merged}
			# parse taxid string for conterminator kingdoms parameter
			{params.script} -t {params.taxdir} -k "{params.kingdoms}" -s "{params.taxdir}/accessionTaxa.sql" -o {output.kstring} -x {output.xstring} &>> {log}
			# run conterminator
			KSTR=$(cat {output.kstring})
			XSTR=$(cat {output.xstring})
			conterminator dna {input.tmp_fna} {output.cmap} {params.prefix} {params.tmpdir} --mask-lower-case 1 --ncbi-tax-dump {params.taxdir} --threads {threads} --split-memory-limit {params.cmem} --blacklist $XSTR --kingdoms $KSTR &>> {log}
			"""

	rule parse_contam_highres:
		input:
			contam = config["rdir"] + "/decontamination/highres_db_conterm_prediction"
		output:
			contam_filt = config["rdir"] + "/decontamination/highres_db_conterm_prediction_filt",
			id_contam = config["rdir"] + "/decontamination/contam_id.accnos"
		shell:
			"""
			awk -v FS="\\t" -v OFS="\\t" '$5 >= 0 && $6 >= 0' {input.contam} > {output.contam_filt}
			cut -f2 {output.contam_filt} | sort | uniq > {output.id_contam}
			"""

rule build_krakendb:
	input:
		gtdb_map = config["rdir"] + "/kraken2_db/library/gtdb/prelim_map.txt",
		gtdb_fasta = config["rdir"] + "/kraken2_db/library/gtdb/library.fna",
		ncbi_map = expand(config["rdir"] + "/kraken2_db/library/{library_highres}/prelim_map.txt", library_highres = LIBRARY_HIGHRES),
		ncbi_fasta = expand(config["rdir"] + "/kraken2_db/library/{library_highres}/library.fna", library_highres = LIBRARY_HIGHRES),
		checkv_map = config["rdir"] + "/kraken2_db/library/checkv/prelim_map.txt",
		checkv_fasta = config["rdir"] + "/kraken2_db/library/checkv/library.fna"
	output:
		hash = config["rdir"] + "/kraken2_db/hash.k2d",
		opts = config["rdir"] + "/kraken2_db/opts.k2d",
		map  = config["rdir"] + "/kraken2_db/seqid2taxid.map",
		taxo = config["rdir"] + "/kraken2_db/taxo.k2d"
	params:
		dbdir = config["rdir"] + "/kraken2_db",
		kmer_len = config["kmer_len"],
		min_len = config["minimizer_len"],
		min_spaces = config["minimizer_spaces"],
		max_dbsize = config["max_dbsize"]
	threads: config["krakenbuild_threads"]
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	log:
		config["rdir"] + "/logs/build_kraken.log"
	shell:
		"""
		kraken2-build --build --threads {threads} --db {params.dbdir} --kmer-len {params.kmer_len} --minimizer-len {params.min_len} --minimizer-spaces {params.min_spaces} --max-db-size {params.max_dbsize} &>> {log}
		"""

