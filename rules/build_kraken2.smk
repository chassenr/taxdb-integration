checkpoint format_taxonomy:
	input:
		expand(config["rdir"] + "/tax_combined/{library_name}_derep_taxonomy.txt", library_name = LIBRARY_NAME),
		config["rdir"] + "/tax_combined/gtdb_derep_taxonomy.txt"
	output:
		tax_combined = config["rdir"] + "/tax_combined/derep_taxonomy_combined.txt",
		nodes = config["rdir"] + "/kraken2_db/taxonomy/nodes.dmp",
		names = config["rdir"] + "/kraken2_db/taxonomy/names.dmp",
		krakendir = directory(config["rdir"] + "/kraken2_genomes")
	params:
		krakendir = config["rdir"] + "/kraken2_genomes",
		genomedir = config["rdir"] + "/derep_combined",
		tax_script = config["tax_script"]
	conda:
		config["wdir"] + "/envs/derep.yaml"
	shell:
		"""
		cat {input} > {output.tax_combined}
		{params.tax_script} --gtdb {output.tax_combined} --assemblies {params.genomedir} --nodes {output.nodes} --names {output.names} --kraken_dir {params.krakendir}
		# replace 'domain' with 'superkingdom (required for kaiju) to ensure that the same nodes.dmp and names.dmp files can be used for both databases
		sed -i -e 's/domain/superkingdom/g' {output.nodes}
		"""

# adapted from: https://github.com/leylabmpi/Struo/blob/f8fdf3d6f04678502fb8d6b094cb4135b7c361e3/bin/kraken2/Snakefile
rule add_krakendb:
	input:
		fasta = config["rdir"] + "/kraken2_genomes/{genome}.gz",
		nodes = config["rdir"] + "/kraken2_db/taxonomy/nodes.dmp",
		names = config["rdir"] + "/kraken2_db/taxonomy/names.dmp"
	output:
		config["rdir"] + "/kraken2_genomes/added/{genome}.done"
	params:
		dbdir = config["rdir"] + "/kraken2_db"
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	shell:
		"""
		kraken2-build --db {params.dbdir} --add-to-library {input.fasta}
		touch {output}
		"""

def aggregate_input(wildcards):
	checkpoint_output = checkpoints.format_taxonomy.get(**wildcards).output.krakendir
	genome_names = expand(config["rdir"] + "/kraken2_genomes/added/{genome}.done",
	genome = glob_wildcards(os.path.join(checkpoint_output,"/added/{genome}.done")).genome) 
	return genome_names

rule build_krakendb:
	input:
		aggregate_input
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
	shell:
		"""
		kraken2-build --build --threads {threads} --db {params.dbdir} --kmer-len {params.kmer_len} --minimizer-len {params.min_len} --minimizer-spaces {params.min_spaces}
		"""

