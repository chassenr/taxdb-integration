rule ps_build_db_highres_nolim_pro:
	input:
		gtdb_map = config["rdir"] + "/kraken2_db_pro/library/gtdb/prelim_map.txt",
		gtdb_fasta = config["rdir"] + "/kraken2_db_pro/library/gtdb/library.fna",
		checkv_map = config["rdir"] + "/kraken2_db_pro/library/checkv/prelim_map.txt",
		checkv_fasta = config["rdir"] + "/kraken2_db_pro/library/checkv/library.fna"
	output:
		hash = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_nolim_pro/hash.k2d",
		opts = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_nolim_pro/opts.k2d",
		map  = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_nolim_pro/seqid2taxid.map",
		taxo = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_nolim_pro/taxo.k2d"
	params:
		libdir = config["rdir"] + "/kraken2_db_pro/library/",
		taxdir = config["rdir"] + "/kraken2_db_pro/taxonomy/",
		dbdir = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_nolim_pro",
		kmer_len = "{kl}",
		min_len = "{ml}",
		min_spaces = "{ms}"
	threads: config["krakenbuild_threads"]
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	log:
		config["sdir"] + "/logs/db_build_highres_k{kl}_m{ml}_s{ms}_nolim_pro.log"
	shell:
		"""
		mkdir -p {params.dbdir}
		cp -r {params.libdir} {params.dbdir}/
		cp -r {params.taxdir} {params.dbdir}/
		/usr/bin/time -v kraken2-build --build --threads {threads} --db {params.dbdir} --kmer-len {params.kmer_len} --minimizer-len {params.min_len} --minimizer-spaces {params.min_spaces} &>> {log}
		"""

rule ps_build_db_highres_mem_pro:
	input:
		gtdb_map = config["rdir"] + "/kraken2_db_pro/library/gtdb/prelim_map.txt",
		gtdb_fasta = config["rdir"] + "/kraken2_db_pro/library/gtdb/library.fna",
		checkv_map = config["rdir"] + "/kraken2_db_pro/library/checkv/prelim_map.txt",
		checkv_fasta = config["rdir"] + "/kraken2_db_pro/library/checkv/library.fna"
	output:
		hash = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_pro/hash.k2d",
		opts = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_pro/opts.k2d",
		map  = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_pro/seqid2taxid.map",
		taxo = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_pro/taxo.k2d"
	params:
		libdir = config["rdir"] + "/kraken2_db_pro/library/",
		taxdir = config["rdir"] + "/kraken2_db_pro/taxonomy/",
		dbdir = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_pro",
		kmer_len = "{kl}",
		min_len = "{ml}",
		min_spaces = "{ms}",
		max_dbsize = "{dbsize}" + "000000000"
	threads: config["krakenbuild_threads"]
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	log:
		config["sdir"] + "/logs/db_build_highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_pro.log"
	shell:
		"""
		mkdir -p {params.dbdir}
		cp -r {params.libdir} {params.dbdir}/
		cp -r {params.taxdir} {params.dbdir}/
		/usr/bin/time -v kraken2-build --build --threads {threads} --db {params.dbdir} --kmer-len {params.kmer_len} --minimizer-len {params.min_len} --minimizer-spaces {params.min_spaces} --max-db-size {params.max_dbsize} &>> {log}
		"""


