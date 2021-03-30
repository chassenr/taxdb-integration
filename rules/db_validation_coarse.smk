rule ps_build_db_coarse_nolim:
	input:
		univec_fasta = config["cdir"] + "/kraken2_db/library/" + config["univec"] + "/library.fna",
		univec_map = config["cdir"] + "/kraken2_db/library/" + config["univec"] + "/prelim_map.txt",
		coarse_fasta = config["cdir"] + "/kraken2_db/library/coarse/library.fna",
		coarse_map = config["cdir"] + "/kraken2_db/library/coarse/prelim_map.txt"
	output:
		hash = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_nolim/hash.k2d",
		opts = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_nolim/opts.k2d",
		map  = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_nolim/seqid2taxid.map",
		taxo = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_nolim/taxo.k2d"
	params:
		libdir = config["cdir"] + "/kraken2_db/library/",
		taxdir = config["cdir"] + "/kraken2_db/taxonomy/",
		dbdir = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_nolim",
		kmer_len = "{kl}",
		min_len = "{ml}",
		min_spaces = "{ms}"
	threads: config["krakenbuild_threads"]
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	log:
		config["sdir"] + "/logs/db_build_coarse_k{kl}_m{ml}_s{ms}_nolim.log"
	shell:
		"""
		mkdir -p {params.dbdir}
		cp -r {params.libdir} {params.dbdir}/
		cp -r {params.taxdir} {params.dbdir}/
		/usr/bin/time -v kraken2-build --build --threads {threads} --db {params.dbdir} --kmer-len {params.kmer_len} --minimizer-len {params.min_len} --minimizer-spaces {params.min_spaces} &>> {log}
		"""

rule ps_build_db_coarse_mem:
	input:
		univec_fasta = config["cdir"] + "/kraken2_db/library/" + config["univec"] + "/library.fna",
		univec_map = config["cdir"] + "/kraken2_db/library/" + config["univec"] + "/prelim_map.txt",
		coarse_fasta = config["cdir"] + "/kraken2_db/library/coarse/library.fna",
		coarse_map = config["cdir"] + "/kraken2_db/library/coarse/prelim_map.txt"
	output:
		hash = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_mem{dbsize}/hash.k2d",
		opts = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_mem{dbsize}/opts.k2d",
		map  = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_mem{dbsize}/seqid2taxid.map",
		taxo = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_mem{dbsize}/taxo.k2d"
	params:
		libdir = config["cdir"] + "/kraken2_db/library/",
		taxdir = config["cdir"] + "/kraken2_db/taxonomy/",
		dbdir = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_mem{dbsize}",
		kmer_len = "{kl}",
		min_len = "{ml}",
		min_spaces = "{ms}",
		max_dbsize = "{dbsize}" + "000000000"
	threads: config["krakenbuild_threads"]
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	log:
		config["sdir"] + "/logs/db_build_coarse_k{kl}_m{ml}_s{ms}_mem{dbsize}.log"
	shell:
		"""
		mkdir -p {params.dbdir}
		cp -r {params.libdir} {params.dbdir}/
		cp -r {params.taxdir} {params.dbdir}/
		/usr/bin/time -v kraken2-build --build --threads {threads} --db {params.dbdir} --kmer-len {params.kmer_len} --minimizer-len {params.min_len} --minimizer-spaces {params.min_spaces} --max-db-size {params.max_dbsize} &>> {log}
		"""


