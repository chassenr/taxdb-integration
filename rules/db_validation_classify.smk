rule ps_classify_coarse_nolim:
	input:
		hash = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_nolim/hash.k2d",
		fq_R1 = config["sdir"] + "/sim_out/sim_PE_{rl}_out/1/R1.fq",
		fq_R2 = config["sdir"] + "/sim_out/sim_PE_{rl}_out/1/R2.fq"
	output:
		kraken = config["sdir"] + "/sim_out/sim_PE_{rl}_out/coarse_k{kl}_m{ml}_s{ms}_nolim_c{conf}.kraken",
		kreport = config["sdir"] + "/sim_out/sim_PE_{rl}_out/coarse_k{kl}_m{ml}_s{ms}_nolim_c{conf}.kreport"
	params:
		dbdir = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_nolim",
		conf = "{conf}"
	threads: config["kraken_classify_threads"]
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	log:
		config["sdir"] + "/logs/coarse_classify_sim_PE_rl{rl}_k{kl}_m{ml}_s{ms}_nolim_c{conf}.log"
	shell:
		"""
		/usr/bin/time -v kraken2 --db {params.dbdir} --threads {threads} --confidence {params.conf} --report {output.kreport} --output {output.kraken} --report-minimizer-data --paired {input.fq_R1} {input.fq_R2} &>> {log}
		"""

rule ps_classify_coarse_mem220:
	input:
		hash = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_mem220/hash.k2d",
		fq_R1 = config["sdir"] + "/sim_out/sim_PE_{rl}_out/1/R1.fq",
		fq_R2 = config["sdir"] + "/sim_out/sim_PE_{rl}_out/1/R2.fq"
	output:
		kraken = config["sdir"] + "/sim_out/sim_PE_{rl}_out/coarse_k{kl}_m{ml}_s{ms}_mem220_c{conf}.kraken",
		kreport = config["sdir"] + "/sim_out/sim_PE_{rl}_out/coarse_k{kl}_m{ml}_s{ms}_mem220_c{conf}.kreport"
	params:
		dbdir = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_mem220",
		conf = "{conf}"
	threads: config["kraken_classify_threads"]
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	log:
		config["sdir"] + "/logs/coarse_classify_sim_PE_rl{rl}_k{kl}_m{ml}_s{ms}_mem220_c{conf}.log"
	shell:
		"""
		/usr/bin/time -v kraken2 --db {params.dbdir} --threads {threads} --confidence {params.conf} --report {output.kreport} --output {output.kraken} --report-minimizer-data --paired {input.fq_R1} {input.fq_R2} &>> {log}
		"""

rule ps_classify_highres_nolim:
	input:
		hash = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_nolim/hash.k2d",
		fq_R1 = config["sdir"] + "/sim_out/sim_PE_{rl}_out/1/R1.fq",
		fq_R2 = config["sdir"] + "/sim_out/sim_PE_{rl}_out/1/R2.fq"
	output:
		kraken = config["sdir"] + "/sim_out/sim_PE_{rl}_out/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf}.kraken",
		kreport = config["sdir"] + "/sim_out/sim_PE_{rl}_out/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf}.kreport"
	params:
		dbdir = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_nolim",
		conf = "{conf}"
	threads: config["kraken_classify_threads"]
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	log:
		config["sdir"] + "/logs/highres_classify_sim_PE_rl{rl}_k{kl}_m{ml}_s{ms}_nolim_c{conf}.log"
	shell:
		"""
		/usr/bin/time -v kraken2 --db {params.dbdir} --threads {threads} --confidence {params.conf} --report {output.kreport} --output {output.kraken} --report-minimizer-data --paired {input.fq_R1} {input.fq_R2} &>> {log}
		"""

rule ps_classify_highres_mem220:
	input:
		hash = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_mem220/hash.k2d",
		fq_R1 = config["sdir"] + "/sim_out/sim_PE_{rl}_out/1/R1.fq",
		fq_R2 = config["sdir"] + "/sim_out/sim_PE_{rl}_out/1/R2.fq"
	output:
		kraken = config["sdir"] + "/sim_out/sim_PE_{rl}_out/highres_k{kl}_m{ml}_s{ms}_mem220_c{conf}.kraken",
		kreport = config["sdir"] + "/sim_out/sim_PE_{rl}_out/highres_k{kl}_m{ml}_s{ms}_mem220_c{conf}.kreport"
	params:
		dbdir = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_mem220",
		conf = "{conf}"
	threads: config["kraken_classify_threads"]
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	log:
		config["sdir"] + "/logs/highres_classify_sim_PE_rl{rl}_k{kl}_m{ml}_s{ms}_mem220_c{conf}.log"
	shell:
		"""
		/usr/bin/time -v kraken2 --db {params.dbdir} --threads {threads} --confidence {params.conf} --report {output.kreport} --output {output.kraken} --report-minimizer-data --paired {input.fq_R1} {input.fq_R2} &>> {log}
		"""

rule conifer_coarse_nolim:
	input:
		kraken = config["sdir"] + "/sim_out/sim_PE_{rl}_out/coarse_k{kl}_m{ml}_s{ms}_nolim_c{conf}.kraken",
		taxo = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_nolim/taxo.k2d"
	output:
		conifer = config["sdir"] + "/sim_out/sim_PE_{rl}_out/coarse_k{kl}_m{ml}_s{ms}_nolim_c{conf}.conifer_out"
	shell:
		"""
		conifer --all --both_scores -i {input.kraken} -d {input.taxo} > {output.conifer}
		"""

rule conifer_coarse_mem220:
	input:
		kraken = config["sdir"] + "/sim_out/sim_PE_{rl}_out/coarse_k{kl}_m{ml}_s{ms}_mem220_c{conf}.kraken",
		taxo = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_nolim/taxo.k2d"
	output:
		conifer = config["sdir"] + "/sim_out/sim_PE_{rl}_out/coarse_k{kl}_m{ml}_s{ms}_mem220_c{conf}.conifer_out"
	shell:
		"""
		conifer --all --both_scores -i {input.kraken} -d {input.taxo} > {output.conifer}
		"""

rule conifer_highres_nolim:
	input:
		kraken = config["sdir"] + "/sim_out/sim_PE_{rl}_out/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf}.kraken",
		taxo = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_nolim/taxo.k2d"
	output:
		conifer = config["sdir"] + "/sim_out/sim_PE_{rl}_out/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf}.conifer_out"
	shell:
		"""
		conifer --all --both_scores -i {input.kraken} -d {input.taxo} > {output.conifer}
		"""

rule conifer_highres_mem220:
	input:
		kraken = config["sdir"] + "/sim_out/sim_PE_{rl}_out/highres_k{kl}_m{ml}_s{ms}_mem220_c{conf}.kraken",
		taxo = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_mem220/taxo.k2d"
	output:
		conifer = config["sdir"] + "/sim_out/sim_PE_{rl}_out/highres_k{kl}_m{ml}_s{ms}_mem220_c{conf}.conifer_out"
	shell:
		"""
		conifer --all --both_scores -i {input.kraken} -d {input.taxo} > {output.conifer}
		"""

rule confusion_matrix_coarse_nolim:
	input:
		conifer = config["sdir"] + "/sim_out/sim_PE_{rl}_out/coarse_k{kl}_m{ml}_s{ms}_nolim_c{conf}.conifer_out",
		names = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_nolim/taxonomy/names.dmp",
		nodes = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_nolim/taxonomy/nodes.dmp",
	output:
		stats_taxongroups = config["sdir"] + "/sim_out/sim_PE_{rl}_out/coarse_k{kl}_m{ml}_s{ms}_nolim_c{conf}_rtl{rtl}_stats_taxongroups.txt",
		stats_domain = config["sdir"] + "/sim_out/sim_PE_{rl}_out/coarse_k{kl}_m{ml}_s{ms}_nolim_c{conf}_rtl{rtl}_stats_domain.txt",
		cm_taxongroups = config["sdir"] + "/sim_out/sim_PE_{rl}_out/coarse_k{kl}_m{ml}_s{ms}_nolim_c{conf}_rtl{rtl}_cm_taxongroups.txt",
		cm_domain = config["sdir"] + "/sim_out/sim_PE_{rl}_out/coarse_k{kl}_m{ml}_s{ms}_nolim_c{conf}_rtl{rtl}_cm_domain.txt"
	params:
		taxdir = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_nolim/taxonomy",
		script = config["wdir"] + "/scripts/sim_results_confusion_coarse.R",
		rtl = "{rtl}",
		prefix = config["sdir"] + "/sim_out/sim_PE_{rl}_out/coarse_k{kl}_m{ml}_s{ms}_nolim_c{conf}_rtl{rtl}"
	conda:
		config["wdir"] + "/envs/r.yaml"
	shell:
		"""
		{params.script} -i {input.conifer} -t {params.taxdir} -s "{params.taxdir}/accessionTaxa.sql" -r {params.rtl} -m "PE" -o {params.prefix}
		"""

rule confusion_matrix_coarse_mem220:
	input:
		conifer = config["sdir"] + "/sim_out/sim_PE_{rl}_out/coarse_k{kl}_m{ml}_s{ms}_mem220_c{conf}.conifer_out",
		names = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_mem220/taxonomy/names.dmp",
		nodes = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_mem220/taxonomy/nodes.dmp",
	output:
		stats_taxongroups = config["sdir"] + "/sim_out/sim_PE_{rl}_out/coarse_k{kl}_m{ml}_s{ms}_mem220_c{conf}_rtl{rtl}_stats_taxongroups.txt",
		stats_domain = config["sdir"] + "/sim_out/sim_PE_{rl}_out/coarse_k{kl}_m{ml}_s{ms}_mem220_c{conf}_rtl{rtl}_stats_domain.txt",
		cm_taxongroups = config["sdir"] + "/sim_out/sim_PE_{rl}_out/coarse_k{kl}_m{ml}_s{ms}_mem220_c{conf}_rtl{rtl}_cm_taxongroups.txt",
		cm_domain = config["sdir"] + "/sim_out/sim_PE_{rl}_out/coarse_k{kl}_m{ml}_s{ms}_mem220_c{conf}_rtl{rtl}_cm_domain.txt"
	params:
		taxdir = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_mem220/taxonomy",
		script = config["wdir"] + "/scripts/sim_results_confusion_coarse.R",
		rtl = "{rtl}",
		prefix = config["sdir"] + "/sim_out/sim_PE_{rl}_out/coarse_k{kl}_m{ml}_s{ms}_mem220_c{conf}_rtl{rtl}"
	conda:
		config["wdir"] + "/envs/r.yaml"
	shell:
		"""
		{params.script} -i {input.conifer} -t {params.taxdir} -s "{params.taxdir}/accessionTaxa.sql" -r {params.rtl} -m "PE" -o {params.prefix}
		"""

rule ps_classify_SE_coarse_nolim:
	input:
		hash = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_nolim/hash.k2d",
		fq = config["sdir"] + "/sim_out/sim_SE_{serl}_out/1/R1.fq"
	output:
		kraken = config["sdir"] + "/sim_out/sim_SE_{serl}_out/coarse_k{kl}_m{ml}_s{ms}_nolim_c{conf}.kraken",
		kreport = config["sdir"] + "/sim_out/sim_SE_{serl}_out/coarse_k{kl}_m{ml}_s{ms}_nolim_c{conf}.kreport"
	params:
		dbdir = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_nolim",
		conf = "{conf}"
	threads: config["kraken_classify_threads"]
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	log:
		config["sdir"] + "/logs/coarse_classify_sim_SE_rl{serl}_k{kl}_m{ml}_s{ms}_nolim_c{conf}.log"
	shell:
		"""
		/usr/bin/time -v kraken2 --db {params.dbdir} --threads {threads} --confidence {params.conf} --report {output.kreport} --output {output.kraken} --report-minimizer-data {input.fq} &>> {log}
		"""

rule ps_classify_SE_coarse_mem220:
	input:
		hash = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_mem220/hash.k2d",
		fq = config["sdir"] + "/sim_out/sim_SE_{serl}_out/1/R1.fq"
	output:
		kraken = config["sdir"] + "/sim_out/sim_SE_{serl}_out/coarse_k{kl}_m{ml}_s{ms}_mem220_c{conf}.kraken",
		kreport = config["sdir"] + "/sim_out/sim_SE_{serl}_out/coarse_k{kl}_m{ml}_s{ms}_mem220_c{conf}.kreport"
	params:
		dbdir = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_mem220",
		conf = "{conf}"
	threads: config["kraken_classify_threads"]
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	log:
		config["sdir"] + "/logs/coarse_classify_sim_SE_rl{serl}_k{kl}_m{ml}_s{ms}_mem220_c{conf}.log"
	shell:
		"""
		/usr/bin/time -v kraken2 --db {params.dbdir} --threads {threads} --confidence {params.conf} --report {output.kreport} --output {output.kraken} --report-minimizer-data {input.fq} &>> {log}
		"""

rule conifer_SE_coarse_nolim:
	input:
		kraken = config["sdir"] + "/sim_out/sim_SE_{serl}_out/coarse_k{kl}_m{ml}_s{ms}_nolim_c{conf}.kraken",
		taxo = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_nolim/taxo.k2d"
	output:
		conifer = config["sdir"] + "/sim_out/sim_SE_{serl}_out/coarse_k{kl}_m{ml}_s{ms}_nolim_c{conf}.conifer_out"
	shell:
		"""
		conifer --all --both_scores -i {input.kraken} -d {input.taxo} > {output.conifer}
		"""

rule conifer_SE_coarse_mem220:
	input:
		kraken = config["sdir"] + "/sim_out/sim_SE_{serl}_out/coarse_k{kl}_m{ml}_s{ms}_mem220_c{conf}.kraken",
		taxo = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_nolim/taxo.k2d"
	output:
		conifer = config["sdir"] + "/sim_out/sim_SE_{serl}_out/coarse_k{kl}_m{ml}_s{ms}_mem220_c{conf}.conifer_out"
	shell:
		"""
		conifer --all --both_scores -i {input.kraken} -d {input.taxo} > {output.conifer}
		"""

rule confusion_matrix_SE_coarse_nolim:
	input:
		conifer = config["sdir"] + "/sim_out/sim_SE_{serl}_out/coarse_k{kl}_m{ml}_s{ms}_nolim_c{conf}.conifer_out",
		names = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_nolim/taxonomy/names.dmp",
		nodes = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_nolim/taxonomy/nodes.dmp",
	output:
		stats_taxongroups = config["sdir"] + "/sim_out/sim_SE_{serl}_out/coarse_k{kl}_m{ml}_s{ms}_nolim_c{conf}_rtl{rtl}_stats_taxongroups.txt",
		stats_domain = config["sdir"] + "/sim_out/sim_SE_{serl}_out/coarse_k{kl}_m{ml}_s{ms}_nolim_c{conf}_rtl{rtl}_stats_domain.txt",
		cm_taxongroups = config["sdir"] + "/sim_out/sim_SE_{serl}_out/coarse_k{kl}_m{ml}_s{ms}_nolim_c{conf}_rtl{rtl}_cm_taxongroups.txt",
		cm_domain = config["sdir"] + "/sim_out/sim_SE_{serl}_out/coarse_k{kl}_m{ml}_s{ms}_nolim_c{conf}_rtl{rtl}_cm_domain.txt"
	params:
		taxdir = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_nolim/taxonomy",
		script = config["wdir"] + "/scripts/sim_results_confusion_coarse.R",
		rtl = "{rtl}",
		prefix = config["sdir"] + "/sim_out/sim_SE_{serl}_out/coarse_k{kl}_m{ml}_s{ms}_nolim_c{conf}_rtl{rtl}"
	conda:
		config["wdir"] + "/envs/r.yaml"
	shell:
		"""
		{params.script} -i {input.conifer} -t {params.taxdir} -s "{params.taxdir}/accessionTaxa.sql" -r {params.rtl} -m "SE" -o {params.prefix}
		"""

rule confusion_matrix_SE_coarse_mem220:
	input:
		conifer = config["sdir"] + "/sim_out/sim_SE_{serl}_out/coarse_k{kl}_m{ml}_s{ms}_mem220_c{conf}.conifer_out",
		names = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_mem220/taxonomy/names.dmp",
		nodes = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_mem220/taxonomy/nodes.dmp",
	output:
		stats_taxongroups = config["sdir"] + "/sim_out/sim_SE_{serl}_out/coarse_k{kl}_m{ml}_s{ms}_mem220_c{conf}_rtl{rtl}_stats_taxongroups.txt",
		stats_domain = config["sdir"] + "/sim_out/sim_SE_{serl}_out/coarse_k{kl}_m{ml}_s{ms}_mem220_c{conf}_rtl{rtl}_stats_domain.txt",
		cm_taxongroups = config["sdir"] + "/sim_out/sim_SE_{serl}_out/coarse_k{kl}_m{ml}_s{ms}_mem220_c{conf}_rtl{rtl}_cm_taxongroups.txt",
		cm_domain = config["sdir"] + "/sim_out/sim_SE_{serl}_out/coarse_k{kl}_m{ml}_s{ms}_mem220_c{conf}_rtl{rtl}_cm_domain.txt"
	params:
		taxdir = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_mem220/taxonomy",
		script = config["wdir"] + "/scripts/sim_results_confusion_coarse.R",
		rtl = "{rtl}",
		prefix = config["sdir"] + "/sim_out/sim_SE_{serl}_out/coarse_k{kl}_m{ml}_s{ms}_mem220_c{conf}_rtl{rtl}"
	conda:
		config["wdir"] + "/envs/r.yaml"
	shell:
		"""
		{params.script} -i {input.conifer} -t {params.taxdir} -s "{params.taxdir}/accessionTaxa.sql" -r {params.rtl} -m "SE" -o {params.prefix}
		"""

