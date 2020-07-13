rule format_taxonomy:
	input:
		config["rdir"] + "/tax_combined/ncbi_derep_taxonomy.txt",
		config["rdir"] + "/tax_combined/gtdb_derep_taxonomy.txt"
	output:
		tax_combined = config["rdir"] + "/tax_combined/derep_taxonomy_combined.txt",
		nodes = config["rdir"] + "/kraken2_db/taxonomy/nodes.dmp",
		names = config["rdir"] + "/kraken2_db/taxonomy/names.dmp"
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
# depending on the number and size of the genomes, it may be required to delete the zipped fna in the derep directory

# adapted from: https://github.com/leylabmpi/Struo/blob/f8fdf3d6f04678502fb8d6b094cb4135b7c361e3/bin/kraken2/Snakefile
GENOMES, = glob_wildcards(config["rdir"] + "/kraken2_genomes/{genome}.fa")
rule add_krakendb:
	input:
		fasta = config["rdir"] + "/kraken2_genomes/{genome}.fa",
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

# depending on the number and size of the genomes, it may be required to delete the output of the tax_from_gtdb.py script (i.e. the kraken2_genomes directory)
# by adding these genomes to the kraken database, they will be copied to the krakendir anyway
rule build_krakendb:
	input:
		expand(config["rdir"] + "/kraken2_genomes/added/{genome}.done", genome = GENOMES)
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

