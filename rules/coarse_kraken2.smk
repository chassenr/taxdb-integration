rule format_taxonomy_coarse:
	input:
		tax_all_coarse = config["cdir"] + "/tax_coarse_all.txt"
	output:
		tax_good = config["cdir"] + "/tax_coarse_good.txt",
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
		cut -f1,2 {input.tax_all_coarse} > {output.tax_good}
		{params.tax_script} --gtdb {output.tax_good} --assemblies {params.genomedir} --nodes {output.nodes} --names {output.names} --kraken_dir {params.krakendir} &>> {log}
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
		fasta = config["cdir"] + "/kraken2_db/tmp/library.fna"
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	threads: config["masking_threads"]
	shell:
		"""
		cut -f1 {input.tax_all_coarse} | grep -F -f - {input.file_list} | parallel -j{threads} 'dustmasker -in {{}} -outfmt fasta' | sed -e '/^>/!s/[a-z]/x/g' >> {output.fasta}
		"""

rule prelim_map_coarse:
	input:
		fasta = config["cdir"] + "/kraken2_db/tmp/library.fna"
	output:
		map = config["cdir"] + "/kraken2_db/tmp/prelim_map.txt"
	params:
		libdir = config["cdir"] + "/kraken2_db/tmp"
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
		fasta = config["cdir"] + "/kraken2_db/tmp/library.fna",
		map = config["cdir"] + "/kraken2_db/tmp/prelim_map.txt",
		nodes = config["cdir"] + "/kraken2_db/taxonomy/nodes.dmp",
                names = config["cdir"] + "/kraken2_db/taxonomy/names.dmp"
	output:
		delnodes = config["cdir"] + "/kraken2_db/taxonomy/delnodes.dmp",
		merged = config["cdir"] + "/kraken2_db/taxonomy/merged.dmp",
		kstring = config["cdir"] + "/decontamination/conterminator_string.txt",
		xstring = config["cdir"] + "/decontamination/conterminator_blacklist.txt",
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
		{params.script} -t {params.taxdir} -k "prokaryotes,fungi,protists,plants,metazoa" -s "{params.taxdir}/accessionTaxa.sql" -o {output.kstring} -x {output.xstring}
		# run conterminator
		KSTR=$(cat {output.kstring})
		XSTR=$(cat {output.xstring})
		conterminator dna {input.fasta} {output.cmap} {params.prefix} {params.tmpdir} --mask-lower-case 1 --ncbi-tax-dump {params.taxdir} --threads {threads} --split-memory-limit {params.cmem} --blacklist $XSTR --kingdoms $KSTR &>> {log}
		"""

rule filter_contamination:
	input:
		contam = config["cdir"] + "/decontamination/coarse_db_conterm_prediction",
		fasta = config["cdir"] + "/kraken2_db/tmp/library.fna",
		map = config["cdir"] + "/kraken2_db/tmp/prelim_map.txt"
	output:
		contam_filt = config["cdir"] + "/decontamination/coarse_db_conterm_prediction_filt",
		id_contam = config["cdir"] + "/decontamination/contam_id.accnos",
		fasta_noncontam = config["cdir"] + "/kraken2_db/library/coarse/library.fna",
		fasta_contam = config["cdir"] + "/decontamination/library_contam.fna",
		map_noncontam = config["cdir"] + "/kraken2_db/library/coarse/prelim_map.txt"
	conda:
		config["wdir"] + "/envs/bbmap.yaml"
	log:
		config["rdir"] + "/logs/coarse_contam_filter.log"
	shell:
		"""
		awk -v FS="\\t" -v OFS="\\t" '$5 >= 0 && $6 >= 0' {input.contam} > {output.contam_filt}
		cut -f2 {output.contam_filt} | sort | uniq > {output.id_contam}
		filterbyname.sh in={input.fasta} out={output.fasta_contam} names={output.id_contam} include=t ow=t &>> {log}
		filterbyname.sh in={input.fasta} out={output.fasta_noncontam} names={output.id_contam} include=f ow=t &>> {log}
		grep -v -F -f {output.id_contam} {input.map} > {output.map_noncontam}
		"""

rule remove_contamination:
	input:
		contam = config["cdir"] + "/decontamination/coarse_db_conterm_prediction_filt",
		fasta_contam = config["cdir"] + "/decontamination/library_contam.fna",
		map_noncontam = config["cdir"] + "/kraken2_db/library/coarse/prelim_map.txt",
		fasta_noncontam = config["cdir"] + "/kraken2_db/library/coarse/library.fna",
	output:
		cleaned_fasta = config["cdir"] + "/decontamination/cleaned.fna",
		cleaned_map = config["cdir"] + "/decontamination/cleaned_map.txt"
	params:
		script = config["wdir"] + "/scripts/remove_contam_contigs.R",
		contam_dir = config["cdir"] + "/decontamination"
	conda:
		config["wdir"] + "/envs/r.yaml"
	log:
		config["rdir"] + "/logs/coarse_contam_cleaning.log"
	shell:
		"""
		{params.script} -i {input.fasta_contam} -c {input.contam} -o {output.cleaned_fasta} &>> {log}
		LC_ALL=C grep '^>' {output.cleaned_fasta} | sed 's/^>//' > "{params.contam_dir}/tmp.accnos"
		NSEQ=$(wc -l "{params.contam_dir}/tmp.accnos" | cut -d' ' -f1)
		printf 'TAXID\\n%.0s' $(seq 1 $NSEQ) | paste - "{params.contam_dir}/tmp.accnos" | paste - <(cut -d'|' -f3 "{params.contam_dir}/tmp.accnos") > {output.cleaned_map}
		rm "{params.contam_dir}/tmp.accnos"
		cat {output.cleaned_map} >> {input.map_noncontam}
		cat {output.cleaned_fasta} >> {input.fasta_noncontam}
		# to be implemented later: remove tmp folder in krakendb to save disk space
		"""

# to avoid ftp issue, recreate kraken2 code for adding UniVec files
# https://github.com/DerrickWood/kraken2/blob/561cc73fababe1dfd996e553e36ea1aff5642ef8/scripts/download_genomic_library.sh#L102-L117
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
		printf 'TAXID\\n%.0s' $(seq 1 $NSEQ) | paste - {params.libdir}/tmp.accnos | paste - <(cut -d'|' -f2 {params.libdir}/tmp.accnos) > {output.map}
		rm {params.libdir}/tmp.accnos
		"""

rule build_kraken_coarse:
	input:
		univec_fasta = config["cdir"] + "/kraken2_db/library/" + config["univec"] + "/library.fna",
		univec_map = config["cdir"] + "/kraken2_db/library/" + config["univec"] + "/prelim_map.txt",
		coarse_fasta = config["cdir"] + "/kraken2_db/library/coarse/library.fna",
		coarse_map = config["cdir"] + "/kraken2_db/library/coarse/prelim_map.txt"
	output:
		hash = config["cdir"] + "/kraken2_db/hash.k2d",
		opts = config["cdir"] + "/kraken2_db/opts.k2d",
		map  = config["cdir"] + "/kraken2_db/seqid2taxid.map",
		taxo = config["cdir"] + "/kraken2_db/taxo.k2d"
	params:
		dbdir = config["cdir"] + "/kraken2_db",
		kmer_len = config["kmer_len"],
		min_len = config["minimizer_len"],
		min_spaces = config["minimizer_spaces"],
		max_dbsize = config["max_dbsize"]
	threads: config["krakenbuild_threads"]
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	log:
		config["rdir"] + "/logs/build_kraken_coarse.log"
	shell:
		"""
		kraken2-build --build --threads {threads} --db {params.dbdir} --kmer-len {params.kmer_len} --minimizer-len {params.min_len} --minimizer-spaces {params.min_spaces} --max-db-size {params.max_dbsize} &>> {log}
		"""

