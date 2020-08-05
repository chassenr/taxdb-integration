rule format_taxonomy:
	input:
		config["rdir"] + "/tax_combined/ncbi_derep_taxonomy.txt",
		config["rdir"] + "/tax_combined/gtdb_derep_taxonomy.txt"
	output:
		tax_combined = config["rdir"] + "/tax_combined/derep_taxonomy_combined.txt",
		nodes = config["rdir"] + "/kraken2_db/taxonomy/nodes.dmp",
		names = config["rdir"] + "/kraken2_db/taxonomy/names.dmp",
		file_list = config["rdir"] + "/kraken2_genomes/file_names_derep_genomes.txt"
	params:
		krakendir = config["rdir"] + "/kraken2_genomes/genome_files",
		genomedir = config["rdir"] + "/derep_combined",
		tax_script = config["tax_script"]
	conda:
		config["wdir"] + "/envs/derep.yaml"
	log:
                config["rdir"] + "/logs/format_taxonomy.log"
	shell:
		"""
		cat {input} > {output.tax_combined}
		{params.tax_script} --gtdb {output.tax_combined} --assemblies {params.genomedir} --nodes {output.nodes} --names {output.names} --kraken_dir {params.krakendir} &>> {log}
		# replace 'domain' with 'superkingdom (required for kaiju) to ensure that the same nodes.dmp and names.dmp files can be used for both databases
		sed -i -e 's/domain/superkingdom/g' {output.nodes}
		find {params.krakendir} -type f -name '*.fa' > {output.file_list}
		"""
# depending on the number and size of the genomes, it may be required to delete the zipped fna in the derep directory

rule add_krakendb_ncbi:
	input:
		file_list = config["rdir"] + "/kraken2_genomes/file_names_derep_genomes.txt",
		ncbi = config["rdir"] + "/tax_combined/{library_name}_derep_taxonomy.txt",
		nodes = config["rdir"] + "/kraken2_db/taxonomy/nodes.dmp",
		names = config["rdir"] + "/kraken2_db/taxonomy/names.dmp"
	output:
		add_complete = config["rdir"] + "/kraken2_genomes/added/{library_name}.done"
	params:
		cat_complete = config["rdir"] + "/kraken2_genomes/added/{library_name}_cat.fa",
		dbdir = config["rdir"] + "/kraken2_db"
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	log:
                config["rdir"] + "/logs/add_krakendb_{library_name}.log"
	shell:
		"""
		cut -f1 {input.ncbi} | grep -F -f - {input.file_list} | xargs cat > {params.cat_complete}
		kraken2-build --db {params.dbdir} --add-to-library {params.cat_complete} &>> {log}
		touch {output.add_complete}
		rm {params.cat_complete}
		"""

rule add_krakendb_gtdb:
	input:
		file_list = config["rdir"] + "/kraken2_genomes/file_names_derep_genomes.txt",
		gtdb = config["rdir"] + "/tax_combined/gtdb_derep_taxonomy.txt",
		nodes = config["rdir"] + "/kraken2_db/taxonomy/nodes.dmp",
		names = config["rdir"] + "/kraken2_db/taxonomy/names.dmp"
	output:
		add_complete = config["rdir"] + "/kraken2_genomes/added/gtdb.done"
	params:
		cat_complete = config["rdir"] + "/kraken2_genomes/added/gtdb_cat.fa",
		dbdir = config["rdir"] + "/kraken2_db"
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	log:
                config["rdir"] + "/logs/add_krakendb_gtdb.log"
	shell:
		"""
		cut -f1 {input.gtdb} | grep -F -f - {input.file_list} | xargs cat > {params.cat_complete}
		kraken2-build --db {params.dbdir} --add-to-library {params.cat_complete} &>> {log}
		touch {output.add_complete}
		rm {params.cat_complete}
		"""
# depending on the number and size of the genomes, it may be required to delete the output of the tax_from_gtdb.py script (i.e. the kraken2_genomes directory)
# by adding these genomes to the kraken database, they will be copied to the krakendir anyway

# to avoid ftp issue, recreate kraken2 code for adding UniVec files
# https://github.com/DerrickWood/kraken2/blob/561cc73fababe1dfd996e553e36ea1aff5642ef8/scripts/download_genomic_library.sh#L102-L117
if config["univec"]:
	rule add_univec:
		output:
			fasta = config["rdir"] + "/kraken2_db/library/" + config["univec"] + "/library.fna",
			map = config["rdir"] + "/kraken2_db/library/" + config["univec"] + "/prelim_map.txt"
		params:
			dbdir = config["rdir"] + "/kraken2_db",
			uv_name = config["univec"],
			scan_script = config["scan_script"],
			mask_script = config["mask_script"]
		conda:
			config["wdir"] + "/envs/kraken2.yaml"
		log:
			config["rdir"] + "/logs/add_krakendb_univec.log"
		shell:
			"""
			wget -O "{params.dbdir}/library/{params.uv_name}/tmp.fna" "ftp://ftp.ncbi.nlm.nih.gov/pub/UniVec/{params.uv_name}"
			sed -e 's/^>/>kraken:taxid|123456|/' "{params.dbdir}/library/{params.uv_name}/tmp.fna" > {output.fasta}
			{params.scan_script} {output.fasta} > {output.map}
			{params.mask_script} {output.fasta}
			"""

rule build_krakendb:
	input:
		config["rdir"] + "/kraken2_genomes/added/gtdb.done",
		expand(config["rdir"] + "/kraken2_genomes/added/{library_name}.done",  library_name = LIBRARY_NAME)
	output:
		hash = config["rdir"] + "/kraken2_db/hash.k2d",
		opts = config["rdir"] + "/kraken2_db/opts.k2d",
		map  = config["rdir"] + "/kraken2_db/seqid2taxid.map",
		taxo = config["rdir"] + "/kraken2_db/taxo.k2d"
	params:
		dbdir = config["rdir"] + "/kraken2_db",
		kmer_len = config["kmer_len"],
		min_len = config["minimizer_len"],
		min_spaces = config["minimizer_spaces"]
	threads: config["krakenbuild_threads"]
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	log:
                config["rdir"] + "/logs/build_kraken.log"
	shell:
		"""
		kraken2-build --build --threads {threads} --db {params.dbdir} --kmer-len {params.kmer_len} --minimizer-len {params.min_len} --minimizer-spaces {params.min_spaces} &>> {log}
		"""

