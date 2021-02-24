rule format_taxonomy_coarse:
	input:
		tax_all_coarse = config["cdir"] + "/tax_coarse_all.txt"
	output:
		nodes = config["cdir"] + "/kraken2_db/taxonomy/nodes.dmp",
		names = config["cdir"] + "/kraken2_db/taxonomy/names.dmp",
		file_list = config["cdir"] + "/kraken2_genomes/file_names_coarse_genomes.txt"
	params:
		krakendir = config["cdir"] + "/kraken2_genomes/genome_files",
		genomedir = config["cdir"] + "/genomes_select",
		tax_script = config["tax_script"]
	conda:
		config["wdir"] + "/envs/biopython.yaml"
	log:
		config["rdir"] + "/logs/format_taxonomy_coarse.log"
	shell:
		"""
		{params.tax_script} --gtdb {input.tax_all_coarse} --assemblies {params.genomedir} --nodes {output.nodes} --names {output.names} --kraken_dir {params.krakendir} &>> {log}
		# replace 'domain' with 'superkingdom (required for kaiju) to ensure that the same nodes.dmp and names.dmp files can be used for both databases
		sed -i -e 's/domain/superkingdom/g' {output.nodes}
		find {params.krakendir} -type f -name '*.fa' > {output.file_list}
		"""

rule masking_coarse:
	input:
		tax_all_coarse = config["cdir"] + "/tax_coarse_all.txt",
		file_list = config["cdir"] + "/kraken2_genomes/file_names_coarse_genomes.txt",
		nodes = config["cdir"] + "/kraken2_db/taxonomy/nodes.dmp",
		names = config["cdir"] + "/kraken2_db/taxonomy/names.dmp"
	output:
		fasta = config["cdir"] + "/kraken2_db/library/coarse/library.fna"
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	threads: config["masking_threads"]
	shell:
		"""
		cut -f1 {input.ncbi} | grep -F -f - {input.file_list} | parallel -j{threads} 'dustmasker -in {{}} -outfmt fasta' | sed -e '/^>/!s/[a-z]/x/g' >> {output.fasta}
		"""

rule prelim_map_coarse:
	input:
		fasta = config["cdir"] + "/kraken2_db/library/coarse/library.fna"
	output:
		map = config["cdir"] + "/kraken2_db/library/coarse/prelim_map.txt"
	params:
		libdir = config["cdir"] + "/kraken2_db/library/coarse"
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	shell:
		"""
		LC_ALL=C grep '^>' {input.fasta} | sed 's/^>//' > {params.libdir}/tmp.accnos
		NSEQ=$(wc -l {params.libdir}/tmp.accnos | cut -d' ' -f1)
		printf 'TAXID\\n%.0s' $(seq 1 $NSEQ) | paste - {params.libdir}/tmp.accnos | paste - <(cut -d'|' -f3 {params.libdir}/tmp.accnos) > {output.map}
		rm {params.libdir}/tmp.accnos
		"""

localrules: detect_contamination

# make sure that conterminator is available in PATH. Manual installation required as conda version is not up to date
rule detect_contamination:
	input:
		fasta = config["cdir"] + "/kraken2_db/library/coarse/library.fna",
		map = config["cdir"] + "/kraken2_db/library/coarse/prelim_map.txt",
		nodes = config["cdir"] + "/kraken2_db/taxonomy/nodes.dmp",
                names = config["cdir"] + "/kraken2_db/taxonomy/names.dmp",
		tax_all_coarse = config["cdir"] + "/tax_coarse_all.txt"
	output:
		delnodes = config["cdir"] + "/kraken2_db/taxonomy/delnodes.dmp",
		merged = config["cdir"] + "/kraken2_db/taxonomy/merged.dmp",
		cstring = config["cdir"] + "/decontamination/conterminator_string.txt",
		cmap = config["cdir"] + "/decontamination/cmap.txt",
		contam = config["cdir"] + "/decontamination/coarse_db_conterm_prediction"
	params:
		script = config["wdir"] + "/scripts/get_kingdoms_conterminator.R",
		tmpdir = config["cdir"] + "/decontamination/tmp",
		taxdir = config["cdir"] + "/kraken2_db/taxonomy/",
		prefix = config["cdir"] + "/decontamination/coarse_db",
		cmem = config["cmem"]
	conda:
		config["wdir"] + "/envs/r.yaml"
	threads: config["masking_threads"]
	log:
		config["rdir"] + "/logs/coarse_conterminator.log"
	shell:
		"""
		# prepare fasta header mapping file for conterminator
		cut -f2,3 {input.map} > {output.cmap}
		# create dummy delnodes and merged files in the taxonomic directory for compatibility with conterminator
		touch {output.delnodes}
		touch {output.merged}
		# parse taxid string for conterminator kingdoms parameter
		{params.script} -i {input.tax_all_coarse} -t {params.taxdir} -s "{params.taxdir}/accessionTaxa.sql" -o {output.cstring}
		# run conterminator
		CSTR=$(cat {output.cstring})
		conterminator dna {input.fasta} {output.cmap} {params.prefix} {params.tmpdir} --mask-lower-case 1 --ncbi-tax-dump {params.taxdir} --threads {threads} --split-memory-limit {params.cmem} --blacklist '5' --kingdoms $CSTR
		"""

rule remove_contamination:
	input:
		contam = config["cdir"] + "/decontamination/coarse_db_conterm_prediction",
		fasta = config["cdir"] + "/kraken2_db/library/coarse/library.fna",
		map = config["cdir"] + "/kraken2_db/library/coarse/prelim_map.txt"
	output:
		cleaned_fasta = config["cdir"] + "/kraken2_db/library/coarse/library.fna",
		cleaned_map = config["cdir"] + "/kraken2_db/library/coarse/prelim_map.txt"
	params:
		script = config["wdir"] + "/scripts/remove_contamination.R"
	conda:
		config["wdir"] + "/envs/r.yaml"
	log:
		config["rdir"] + "/logs/coarse_cleaning.log"
	shell:
		"""
		"""

if config["univec"]:
	rule add_univec:
		output:
			fasta = config["cdir"] + "/kraken2_db/library/" + config["univec"] + "/library.fna",
			map = config["cdir"] + "/kraken2_db/library/" + config["univec"] + "/prelim_map.txt"
		params:
			ncbi_server = config["ncbi_server"],
			uv_name = config["univec"],
			libdir = config["cdir"] + "/kraken2_db/library/" + config["univec"]
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

rule build_kraken_coarse:
	input:
	output:
	params:
	threads: config["krakenbuild_threads"]
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	log:
		config["rdir"] + "/logs/build_kraken.log"
	shell:
		"""
		kraken2-build --build --threads {threads} --db {params.dbdir} --kmer-len {params.kmer_len} --minimizer-len {params.min_len} --minimizer-spaces {params.min_spaces} &>> {log}
		"""

