rule fix_ncbi_taxpath:
	input:
		ncbi = expand(config["rdir"] + "/tax_combined/{library_name}_derep_taxonomy.txt", library_name = LIBRARY_NAME),
		gtdb = config["rdir"] + "/tax_combined/gtdb_derep_taxonomy.txt",
		checkv = config["rdir"] + "/tax_combined/checkv_derep_taxonomy.txt"
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
                cat {input.ncbi} {input.checkv} > "{params.outdir}/tmp"
                {params.script} -i "{params.outdir}/tmp" -g {input.gtdb} -o {output.tax_combined} &>> {log}
                rm "{params.outdir}/tmp"
                """


rule format_taxonomy:
	input:
		tax_combined = config["rdir"] + "/tax_combined/derep_taxonomy_combined.txt"
	output:
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
		{params.tax_script} --gtdb {output.tax_combined} --assemblies {params.genomedir} --nodes {output.nodes} --names {output.names} --kraken_dir {params.krakendir} &>> {log}
		# replace 'domain' with 'superkingdom (required for kaiju) to ensure that the same nodes.dmp and names.dmp files can be used for both databases
		sed -i -e 's/domain/superkingdom/g' {output.nodes}
		find {params.krakendir} -type f -name '*.fa' > {output.file_list}
		"""
# depending on the number and size of the genomes, it may be required to delete the zipped fna in the derep directory

rule masking_ncbi:
	input:
		file_list = config["rdir"] + "/kraken2_genomes/file_names_derep_genomes.txt",
		ncbi = config["rdir"] + "/tax_combined/{library_name}_derep_taxonomy.txt",
		nodes = config["rdir"] + "/kraken2_db/taxonomy/nodes.dmp",
		names = config["rdir"] + "/kraken2_db/taxonomy/names.dmp"
	output:
		fasta = config["rdir"] + "/kraken2_db/library/{library_name}/library.fna"
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	threads: config["masking_threads"]
	shell:
		"""
		cut -f1 {input.ncbi} | grep -F -f - {input.file_list} | parallel -j{threads} 'dustmasker -in {{}} -outfmt fasta' | sed -e '/^>/!s/[a-z]/x/g' >> {output.fasta}
		"""

rule masking_gtdb:
	input:
		file_list = config["rdir"] + "/kraken2_genomes/file_names_derep_genomes.txt",
		gtdb = config["rdir"] + "/tax_combined/gtdb_derep_taxonomy.txt",
		nodes = config["rdir"] + "/kraken2_db/taxonomy/nodes.dmp",
		names = config["rdir"] + "/kraken2_db/taxonomy/names.dmp"
	output:
		fasta = config["rdir"] + "/kraken2_db/library/gtdb/library.fna"
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	threads: config["masking_threads"]
	shell:
		"""
		cut -f1 {input.gtdb} | grep -F -f - {input.file_list} | parallel -j{threads} 'dustmasker -in {{}} -outfmt fasta' | sed -e '/^>/!s/[a-z]/x/g' >> {output.fasta}
		"""

rule prelim_map_ncbi:
	input:
		fasta = config["rdir"] + "/kraken2_db/library/{library_name}/library.fna"
	output:
		map = config["rdir"] + "/kraken2_db/library/{library_name}/prelim_map.txt"
	params:
		libdir = config["rdir"] + "/kraken2_db/library/{library_name}"
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	shell:
		"""
		LC_ALL=C grep '^>' {input.fasta} | sed 's/^>//' > {params.libdir}/tmp.accnos
		NSEQ=$(wc -l {params.libdir}/tmp.accnos | cut -d' ' -f1)
		printf 'TAXID\\n%.0s' $(seq 1 $NSEQ) | paste - {params.libdir}/tmp.accnos | paste - <(cut -d'|' -f3 {params.libdir}/tmp.accnos) > {output.map}
		rm {params.libdir}/tmp.accnos
		"""

rule prelim_map_gtdb:
	input:  
		fasta = config["rdir"] + "/kraken2_db/library/gtdb/library.fna"
	output:
		map = config["rdir"] + "/kraken2_db/library/gtdb/prelim_map.txt"
	params: 
		libdir = config["rdir"] + "/kraken2_db/library/gtdb"
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	shell:
		"""
		LC_ALL=C grep '^>' {input.fasta} | sed 's/^>//' > {params.libdir}/tmp.accnos
		NSEQ=$(wc -l {params.libdir}/tmp.accnos | cut -d' ' -f1)
		printf 'TAXID\\n%.0s' $(seq 1 $NSEQ) | paste - {params.libdir}/tmp.accnos | paste - <(cut -d'|' -f3 {params.libdir}/tmp.accnos) > {output.map}
		rm {params.libdir}/tmp.accnos
		"""

# to avoid ftp issue, recreate kraken2 code for adding UniVec files
# https://github.com/DerrickWood/kraken2/blob/561cc73fababe1dfd996e553e36ea1aff5642ef8/scripts/download_genomic_library.sh#L102-L117
if config["univec"]:
	rule add_univec:
		output:
			fasta = config["rdir"] + "/kraken2_db/library/" + config["univec"] + "/library.fna",
			map = config["rdir"] + "/kraken2_db/library/" + config["univec"] + "/prelim_map.txt"
		params:
			ncbi_server = config["ncbi_server"],
			uv_name = config["univec"],
			libdir = config["rdir"] + "/kraken2_db/library/" + config["univec"]
		conda:
			config["wdir"] + "/envs/kraken2.yaml"
		log:
			config["rdir"] + "/logs/add_krakendb_univec.log"
		shell:
			"""
			wget -O "{params.libdir}/tmp.fna" "{params.ncbi_server}/pub/UniVec/{params.uv_name}"
			# choosing random artificial taxid (this taxid must not exist elsewhere in the database)
			sed -i 's/^>/>kraken:taxid|1234567|/' "{params.libdir}/tmp.fna"
			dustmasker -in "{params.libdir}/tmp.fna" -outfmt fasta | sed -e '/^>/!s/[a-z]/x/g' > {output.fasta}
			rm "{params.libdir}/tmp.fna"
			grep '^>' {output.fasta} | sed 's/^>//' > {params.libdir}/tmp.accnos
			NSEQ=$(wc -l {params.libdir}/tmp.accnos | cut -d' ' -f1)
			printf 'TAXID\\n%.0s' $(seq 1 $NSEQ) | paste - {params.libdir}/tmp.accnos | paste - <(cut -d'|' -f3 {params.libdir}/tmp.accnos) > {output.map}
			rm {params.libdir}/tmp.accnos
			"""

rule build_krakendb:
	input:
		gtdb_map = config["rdir"] + "/kraken2_db/library/gtdb/prelim_map.txt",
		gtdb_fasta = config["rdir"] + "/kraken2_db/library/gtdb/library.fna",
		ncbi_map = expand(config["rdir"] + "/kraken2_db/library/{library_name}/prelim_map.txt", library_name = LIBRARY_NAME),
		ncbi_fasta = expand(config["rdir"] + "/kraken2_db/library/{library_name}/library.fna", library_name = LIBRARY_NAME),
		univec_map = config["rdir"] + "/kraken2_db/library/" + config["univec"] + "/prelim_map.txt" if config["univec"] else [],
		univec_fasta = config["rdir"] + "/kraken2_db/library/" + config["univec"] + "/library.fna" if config["univec"] else []
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

