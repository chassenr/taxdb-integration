rule ps_classify_coarse_nolim:
	input:
		hash = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_nolim/hash.k2d",
		fq = config["sdir"] + "/sim_out/sim_out_70bp_SE/1/R1.fq"
	output:
		kraken = config["sdir"] + "/sim_out/sim_out_70bp_SE/coarse_k{kl}_m{ml}_s{ms}_nolim_c{conf}.kraken",
		kreport = config["sdir"] + "/sim_out/sim_out_70bp_SE/coarse_k{kl}_m{ml}_s{ms}_nolim_c{conf}.kreport"
	params:
		dbdir = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_nolim",
		conf = "{conf}"
	threads: config["kraken_classify_threads"]
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	log:
		config["sdir"] + "/logs/coarse_classify_k{kl}_m{ml}_s{ms}_nolim_c{conf}.log"
	shell:
		"""
		/usr/bin/time -v kraken2 --db {params.dbdir} --threads {threads} --confidence {params.conf} --report {output.kreport} --output {output.kraken} --report-minimizer-data {input.fq} &>> {log}
		"""

rule ps_classify_coarse_mem220:
	input:
		hash = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_mem220/hash.k2d",
		fq = config["sdir"] + "/sim_out/sim_out_70bp_SE/1/R1.fq"
	output:
		kraken = config["sdir"] + "/sim_out/sim_out_70bp_SE/coarse_k{kl}_m{ml}_s{ms}_mem220_c{conf}.kraken",
		kreport = config["sdir"] + "/sim_out/sim_out_70bp_SE/coarse_k{kl}_m{ml}_s{ms}_mem220_c{conf}.kreport"
	params:
		dbdir = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_mem220",
		conf = "{conf}"
	threads: config["kraken_classify_threads"]
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	log:
		config["sdir"] + "/logs/coarse_classify_k{kl}_m{ml}_s{ms}_mem220_c{conf}.log"
	shell:
		"""
		/usr/bin/time -v kraken2 --db {params.dbdir} --threads {threads} --confidence {params.conf} --report {output.kreport} --output {output.kraken} --report-minimizer-data {input.fq} &>> {log}
		"""

rule ps_classify_highres_nolim:
	input:
		hash = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_nolim/hash.k2d",
		fq = config["sdir"] + "/sim_out/sim_out_70bp_SE/1/R1.fq"
	output:
		kraken = config["sdir"] + "/sim_out/sim_out_70bp_SE/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf}.kraken",
		kreport = config["sdir"] + "/sim_out/sim_out_70bp_SE/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf}.kreport"
	params:
		dbdir = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_nolim",
		conf = "{conf}"
	threads: config["kraken_classify_threads"]
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	log:
		config["sdir"] + "/logs/highres_classify_k{kl}_m{ml}_s{ms}_nolim_c{conf}.log"
	shell:
		"""
		/usr/bin/time -v kraken2 --db {params.dbdir} --threads {threads} --confidence {params.conf} --report {output.kreport} --output {output.kraken} --report-minimizer-data {input.fq} &>> {log}
		"""

rule ps_classify_highres_mem220:
	input:
		hash = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_mem220/hash.k2d",
		fq = config["sdir"] + "/sim_out/sim_out_70bp_SE/1/R1.fq"
	output:
		kraken = config["sdir"] + "/sim_out/sim_out_70bp_SE/highres_k{kl}_m{ml}_s{ms}_mem220_c{conf}.kraken",
		kreport = config["sdir"] + "/sim_out/sim_out_70bp_SE/highres_k{kl}_m{ml}_s{ms}_mem220_c{conf}.kreport"
	params:
		dbdir = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_mem220",
		conf = "{conf}"
	threads: config["kraken_classify_threads"]
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	log:
		config["sdir"] + "/logs/highres_classify_k{kl}_m{ml}_s{ms}_mem220_c{conf}.log"
	shell:
		"""
		/usr/bin/time -v kraken2 --db {params.dbdir} --threads {threads} --confidence {params.conf} --report {output.kreport} --output {output.kraken} --report-minimizer-data {input.fq} &>> {log}
		"""

rule conifer_coarse_nolim:
	input:
		kraken = config["sdir"] + "/sim_out/sim_out_70bp_SE/coarse_k{kl}_m{ml}_s{ms}_nolim_c{conf}.kraken",
		taxo = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_nolim/taxo.k2d
	output:
		conifer = config["sdir"] + "/sim_out/sim_out_70bp_SE/coarse_k{kl}_m{ml}_s{ms}_nolim_c{conf}.conifer_out
	shell:
		"""
		conifer --all --both_scores -i {input.kraken} -d {input.taxo} > {output.conifer}
		"""

rule conifer_coarse_mem220:
	input:
		kraken = config["sdir"] + "/sim_out/sim_out_70bp_SE/coarse_k{kl}_m{ml}_s{ms}_mem220_c{conf}.kraken",
		taxo = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_nolim/taxo.k2d
	output:
		conifer = config["sdir"] + "/sim_out/sim_out_70bp_SE/coarse_k{kl}_m{ml}_s{ms}_mem220_c{conf}.conifer_out
	shell:
		"""
		conifer --all --both_scores -i {input.kraken} -d {input.taxo} > {output.conifer}
		"""

rule conifer_highres_nolim:
	input:
		kraken = config["sdir"] + "/sim_out/sim_out_70bp_SE/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf}.kraken",
		taxo = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_nolim/taxo.k2d
	output:
		conifer = config["sdir"] + "/sim_out/sim_out_70bp_SE/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf}.conifer_out
	shell:
		"""
		conifer --all --both_scores -i {input.kraken} -d {input.taxo} > {output.conifer}
		"""

rule conifer_highres_mem220:
	input:
		kraken = config["sdir"] + "/sim_out/sim_out_70bp_SE/highres_k{kl}_m{ml}_s{ms}_mem220_c{conf}.kraken",
		taxo = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_mem220/taxo.k2d
	output:
		conifer = config["sdir"] + "/sim_out/sim_out_70bp_SE/highres_k{kl}_m{ml}_s{ms}_mem220_c{conf}.conifer_out
	shell:
		"""
		conifer --all --both_scores -i {input.kraken} -d {input.taxo} > {output.conifer}
		"""

rule filter_reads_nolim:
	input:
		conifer = config["sdir"] + "/sim_out/sim_out_70bp_SE/coarse_k{kl}_m{ml}_s{ms}_nolim_c{conf}.conifer_out,
		names = ,
		nodes = ,
	output:
		
	params:
		script = ,
	conda:
		config["wdir"] + "/envs/r.yaml"
	shell:
		"""
		
		"""
