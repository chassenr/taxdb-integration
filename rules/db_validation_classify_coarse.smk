rule ps_classify_PE_coarse_nolim:
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

rule ps_classify_PE_coarse_mem:
	input:
		hash = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_mem{dbsize}/hash.k2d",
		fq_R1 = config["sdir"] + "/sim_out/sim_PE_{rl}_out/1/R1.fq",
		fq_R2 = config["sdir"] + "/sim_out/sim_PE_{rl}_out/1/R2.fq"
	output:
		kraken = config["sdir"] + "/sim_out/sim_PE_{rl}_out/coarse_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf}.kraken",
		kreport = config["sdir"] + "/sim_out/sim_PE_{rl}_out/coarse_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf}.kreport"
	params:
		dbdir = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_mem{dbsize}",
		conf = "{conf}"
	threads: config["kraken_classify_threads"]
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	log:
		config["sdir"] + "/logs/coarse_classify_sim_PE_rl{rl}_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf}.log"
	shell:
		"""
		/usr/bin/time -v kraken2 --db {params.dbdir} --threads {threads} --confidence {params.conf} --report {output.kreport} --output {output.kraken} --report-minimizer-data --paired {input.fq_R1} {input.fq_R2} &>> {log}
		"""

rule conifer_PE_coarse_nolim:
	input:
		kraken = config["sdir"] + "/sim_out/sim_PE_{rl}_out/coarse_k{kl}_m{ml}_s{ms}_nolim_c{conf}.kraken",
		taxo = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_nolim/taxo.k2d"
	output:
		conifer = config["sdir"] + "/sim_out/sim_PE_{rl}_out/coarse_k{kl}_m{ml}_s{ms}_nolim_c{conf}.conifer_out"
	shell:
		"""
		conifer --all --both_scores -i {input.kraken} -d {input.taxo} > {output.conifer}
		"""

rule conifer_PE_coarse_mem:
	input:
		kraken = config["sdir"] + "/sim_out/sim_PE_{rl}_out/coarse_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf}.kraken",
		taxo = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_mem{dbsize}/taxo.k2d"
	output:
		conifer = config["sdir"] + "/sim_out/sim_PE_{rl}_out/coarse_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf}.conifer_out"
	shell:
		"""
		conifer --all --both_scores -i {input.kraken} -d {input.taxo} > {output.conifer}
		"""

rule confusion_matrix_PE_coarse_nolim:
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

rule confusion_matrix_PE_coarse_mem:
	input:
		conifer = config["sdir"] + "/sim_out/sim_PE_{rl}_out/coarse_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf}.conifer_out",
		names = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_mem{dbsize}/taxonomy/names.dmp",
		nodes = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_mem{dbsize}/taxonomy/nodes.dmp",
	output:
		stats_taxongroups = config["sdir"] + "/sim_out/sim_PE_{rl}_out/coarse_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf}_rtl{rtl}_stats_taxongroups.txt",
		stats_domain = config["sdir"] + "/sim_out/sim_PE_{rl}_out/coarse_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf}_rtl{rtl}_stats_domain.txt",
		cm_taxongroups = config["sdir"] + "/sim_out/sim_PE_{rl}_out/coarse_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf}_rtl{rtl}_cm_taxongroups.txt",
		cm_domain = config["sdir"] + "/sim_out/sim_PE_{rl}_out/coarse_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf}_rtl{rtl}_cm_domain.txt"
	params:
		taxdir = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_mem{dbsize}/taxonomy",
		script = config["wdir"] + "/scripts/sim_results_confusion_coarse.R",
		rtl = "{rtl}",
		prefix = config["sdir"] + "/sim_out/sim_PE_{rl}_out/coarse_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf}_rtl{rtl}"
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

rule ps_classify_SE_coarse_mem:
	input:
		hash = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_mem{dbsize}/hash.k2d",
		fq = config["sdir"] + "/sim_out/sim_SE_{serl}_out/1/R1.fq"
	output:
		kraken = config["sdir"] + "/sim_out/sim_SE_{serl}_out/coarse_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf}.kraken",
		kreport = config["sdir"] + "/sim_out/sim_SE_{serl}_out/coarse_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf}.kreport"
	params:
		dbdir = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_mem{dbsize}",
		conf = "{conf}"
	threads: config["kraken_classify_threads"]
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	log:
		config["sdir"] + "/logs/coarse_classify_sim_SE_rl{serl}_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf}.log"
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

rule conifer_SE_coarse_mem:
	input:
		kraken = config["sdir"] + "/sim_out/sim_SE_{serl}_out/coarse_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf}.kraken",
		taxo = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_mem{dbsize}/taxo.k2d"
	output:
		conifer = config["sdir"] + "/sim_out/sim_SE_{serl}_out/coarse_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf}.conifer_out"
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

rule confusion_matrix_SE_coarse_mem:
	input:
		conifer = config["sdir"] + "/sim_out/sim_SE_{serl}_out/coarse_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf}.conifer_out",
		names = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_mem{dbsize}/taxonomy/names.dmp",
		nodes = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_mem{dbsize}/taxonomy/nodes.dmp",
	output:
		stats_taxongroups = config["sdir"] + "/sim_out/sim_SE_{serl}_out/coarse_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf}_rtl{rtl}_stats_taxongroups.txt",
		stats_domain = config["sdir"] + "/sim_out/sim_SE_{serl}_out/coarse_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf}_rtl{rtl}_stats_domain.txt",
		cm_taxongroups = config["sdir"] + "/sim_out/sim_SE_{serl}_out/coarse_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf}_rtl{rtl}_cm_taxongroups.txt",
		cm_domain = config["sdir"] + "/sim_out/sim_SE_{serl}_out/coarse_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf}_rtl{rtl}_cm_domain.txt"
	params:
		taxdir = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_mem{dbsize}/taxonomy",
		script = config["wdir"] + "/scripts/sim_results_confusion_coarse.R",
		rtl = "{rtl}",
		prefix = config["sdir"] + "/sim_out/sim_SE_{serl}_out/coarse_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf}_rtl{rtl}"
	conda:
		config["wdir"] + "/envs/r.yaml"
	shell:
		"""
		{params.script} -i {input.conifer} -t {params.taxdir} -s "{params.taxdir}/accessionTaxa.sql" -r {params.rtl} -m "SE" -o {params.prefix}
		"""


