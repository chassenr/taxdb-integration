rule ps_build_db_highres_nolim_euk:
	input:
		ncbi_map = expand(config["rdir"] + "/kraken2_db/library/{library_highres}/prelim_map.txt", library_name = LIBRARY_HIGHRES),
		ncbi_fasta = expand(config["rdir"] + "/kraken2_db/library/{library_highres}/library.fna", library_name = LIBRARY_HIGHRES)
	output:
		hash = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_nolim_euk/hash.k2d",
		opts = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_nolim_euk/opts.k2d",
		map  = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_nolim_euk/seqid2taxid.map",
		taxo = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_nolim_euk/taxo.k2d"
	params:
		libdir = config["rdir"] + "/kraken2_db_euk/library/",
		taxdir = config["rdir"] + "/kraken2_db_euk/taxonomy/",
		dbdir = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_nolim_euk",
		kmer_len = "{kl}",
		min_len = "{ml}",
		min_spaces = "{ms}"
	threads: config["krakenbuild_threads"]
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	log:
		config["sdir"] + "/logs/db_build_highres_k{kl}_m{ml}_s{ms}_nolim_euk.log"
	shell:
		"""
		mkdir -p {params.dbdir}
		cp -r {params.libdir} {params.dbdir}/
		cp -r {params.taxdir} {params.dbdir}/
		/usr/bin/time -v kraken2-build --build --threads {threads} --db {params.dbdir} --kmer-len {params.kmer_len} --minimizer-len {params.min_len} --minimizer-spaces {params.min_spaces} &>> {log}
		"""

rule ps_build_db_highres_mem_euk:
	input:
		ncbi_map = expand(config["rdir"] + "/kraken2_db/library/{library_highres}/prelim_map.txt", library_name = LIBRARY_HIGHRES),
		ncbi_fasta = expand(config["rdir"] + "/kraken2_db/library/{library_highres}/library.fna", library_name = LIBRARY_HIGHRES)
	output:
		hash = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_euk/hash.k2d",
		opts = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_euk/opts.k2d",
		map  = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_euk/seqid2taxid.map",
		taxo = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_euk/taxo.k2d"
	params:
		libdir = config["rdir"] + "/kraken2_db_euk/library/",
		taxdir = config["rdir"] + "/kraken2_db_euk/taxonomy/",
		dbdir = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_euk",
		kmer_len = "{kl}",
		min_len = "{ml}",
		min_spaces = "{ms}",
		max_dbsize = "{dbsize}" + "000000000"
	threads: config["krakenbuild_threads"]
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	log:
		config["sdir"] + "/logs/db_build_highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_euk.log"
	shell:
		"""
		mkdir -p {params.dbdir}
		cp -r {params.libdir} {params.dbdir}/
		cp -r {params.taxdir} {params.dbdir}/
		/usr/bin/time -v kraken2-build --build --threads {threads} --db {params.dbdir} --kmer-len {params.kmer_len} --minimizer-len {params.min_len} --minimizer-spaces {params.min_spaces} --max-db-size {params.max_dbsize} &>> {log}
		"""


