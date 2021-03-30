rule ps_classify_PE_highres_nolim_pro:
	input:
		hash = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_nolim_pro/hash.k2d",
		fq_R1 = config["sdir"] + "/sim_out/sim_PE_{rl}_out/1/R1.fq",
		fq_R2 = config["sdir"] + "/sim_out/sim_PE_{rl}_out/1/R2.fq"
	output:
		kraken = config["sdir"] + "/sim_out/sim_PE_{rl}_out/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf}_pro.kraken",
		kreport = config["sdir"] + "/sim_out/sim_PE_{rl}_out/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf}_pro.kreport"
	params:
		dbdir = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_nolim_pro",
		conf = "{conf}"
	threads: config["kraken_classify_threads"]
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	log:
		config["sdir"] + "/logs/highres_classify_sim_PE_rl{rl}_k{kl}_m{ml}_s{ms}_nolim_c{conf}_pro.log"
	shell:
		"""
		/usr/bin/time -v kraken2 --db {params.dbdir} --threads {threads} --confidence {params.conf} --report {output.kreport} --output {output.kraken} --report-minimizer-data --paired {input.fq_R1} {input.fq_R2} &>> {log}
		"""

rule ps_classify_PE_highres_mem_pro:
	input:
		hash = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_pro/hash.k2d",
		fq_R1 = config["sdir"] + "/sim_out/sim_PE_{rl}_out/1/R1.fq",
		fq_R2 = config["sdir"] + "/sim_out/sim_PE_{rl}_out/1/R2.fq"
	output:
		kraken = config["sdir"] + "/sim_out/sim_PE_{rl}_out/highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf}_pro.kraken",
		kreport = config["sdir"] + "/sim_out/sim_PE_{rl}_out/highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf}_pro.kreport"
	params:
		dbdir = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_pro",
		conf = "{conf}"
	threads: config["kraken_classify_threads"]
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	log:
		config["sdir"] + "/logs/highres_classify_sim_PE_rl{rl}_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf}_pro.log"
	shell:
		"""
		/usr/bin/time -v kraken2 --db {params.dbdir} --threads {threads} --confidence {params.conf} --report {output.kreport} --output {output.kraken} --report-minimizer-data --paired {input.fq_R1} {input.fq_R2} &>> {log}
		"""

rule conifer_PE_highres_nolim_pro:
	input:
		kraken = config["sdir"] + "/sim_out/sim_PE_{rl}_out/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf}_pro.kraken",
		taxo = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_nolim_pro/taxo.k2d"
	output:
		conifer = config["sdir"] + "/sim_out/sim_PE_{rl}_out/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf}_pro.conifer_out"
	shell:
		"""
		conifer --all --both_scores -i {input.kraken} -d {input.taxo} > {output.conifer}
		"""

rule conifer_PE_highres_mem_pro:
	input:
		kraken = config["sdir"] + "/sim_out/sim_PE_{rl}_out/highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf}_pro.kraken",
		taxo = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_pro/taxo.k2d"
	output:
		conifer = config["sdir"] + "/sim_out/sim_PE_{rl}_out/highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf}_pro.conifer_out"
	shell:
		"""
		conifer --all --both_scores -i {input.kraken} -d {input.taxo} > {output.conifer}
		"""

rule confusion_matrix_PE_highres_nolim_pro:
	input:
		conifer = config["sdir"] + "/sim_out/sim_PE_{rl}_out/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf}_pro.conifer_out",
		coarse = config["sdir"] + "/sim_out/sim_SE_{serl}_out/" + config["best_coarse_PE"] + ".conifer_out",
		names = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_nolim_pro/taxonomy/names.dmp",
		nodes = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_nolim_pro/taxonomy/nodes.dmp",
	output:
		stats_taxlevels = config["sdir"] + "/sim_out/sim_PE_{rl}_out/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf}_rtl{rtl}_pro_stats_taxlevels.txt",
		stats_domain = config["sdir"] + "/sim_out/sim_PE_{rl}_out/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf}_rtl{rtl}_pro_stats_domain.txt",
		cm_taxlevels = config["sdir"] + "/sim_out/sim_PE_{rl}_out/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf}_rtl{rtl}_pro_cm_taxlevels.txt",
		cm_domain = config["sdir"] + "/sim_out/sim_PE_{rl}_out/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf}_rtl{rtl}_pro_cm_domain.txt"
	params:
		taxdir = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_nolim_pro/taxonomy",
		script = config["wdir"] + "/scripts/sim_results_confusion_highres_pro.R",
		rtl = "{rtl}",
		prefix = config["sdir"] + "/sim_out/sim_PE_{rl}_out/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf}_rtl{rtl}_pro"
	conda:
		config["wdir"] + "/envs/r.yaml"
	shell:
		"""
		{params.script} -i {input.conifer} -c {input.coarse} -t {params.taxdir} -s "{params.taxdir}/accessionTaxa.sql" -r {params.rtl} -m "PE" -o {params.prefix}
		"""

rule confusion_matrix_PE_highres_mem:
	input:
		conifer = config["sdir"] + "/sim_out/sim_PE_{rl}_out/highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf}_pro.conifer_out",
		coarse = config["sdir"] + "/sim_out/sim_SE_{serl}_out/" + config["best_coarse_PE"] + ".conifer_out",
		names = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_pro/taxonomy/names.dmp",
		nodes = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_pro/taxonomy/nodes.dmp",
	output:
		stats_taxlevels = config["sdir"] + "/sim_out/sim_PE_{rl}_out/highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf}_rtl{rtl}_pro_stats_taxlevels.txt",
		stats_domain = config["sdir"] + "/sim_out/sim_PE_{rl}_out/highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf}_rtl{rtl}_pro_stats_domain.txt",
		cm_taxlevels = config["sdir"] + "/sim_out/sim_PE_{rl}_out/highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf}_rtl{rtl}_pro_cm_taxlevels.txt",
		cm_domain = config["sdir"] + "/sim_out/sim_PE_{rl}_out/highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf}_rtl{rtl}_pro_cm_domain.txt"
	params:
		taxdir = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_pro/taxonomy",
		script = config["wdir"] + "/scripts/sim_results_confusion_highres_pro.R",
		rtl = "{rtl}",
		prefix = config["sdir"] + "/sim_out/sim_PE_{rl}_out/highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf}_rtl{rtl}_pro"
	conda:
		config["wdir"] + "/envs/r.yaml"
	shell:
		"""
		{params.script} -i {input.conifer} -c {input.coarse} -t {params.taxdir} -s "{params.taxdir}/accessionTaxa.sql" -r {params.rtl} -m "PE" -o {params.prefix}
		"""

rule ps_classify_SE_highres_nolim_pro:
	input:
		hash = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_nolim_pro/hash.k2d",
		fq = config["sdir"] + "/sim_out/sim_SE_{serl}_out/1/R1.fq"
	output:
		kraken = config["sdir"] + "/sim_out/sim_SE_{serl}_out/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf}_pro.kraken",
		kreport = config["sdir"] + "/sim_out/sim_SE_{serl}_out/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf}_pro.kreport"
	params:
		dbdir = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_nolim_pro",
		conf = "{conf}"
	threads: config["kraken_classify_threads"]
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	log:
		config["sdir"] + "/logs/highres_classify_sim_SE_rl{serl}_k{kl}_m{ml}_s{ms}_nolim_c{conf}_pro.log"
	shell:
		"""
		/usr/bin/time -v kraken2 --db {params.dbdir} --threads {threads} --confidence {params.conf} --report {output.kreport} --output {output.kraken} --report-minimizer-data {input.fq} &>> {log}
		"""

rule ps_classify_SE_highres_mem_pro:
	input:
		hash = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_pro/hash.k2d",
		fq = config["sdir"] + "/sim_out/sim_SE_{serl}_out/1/R1.fq"
	output:
		kraken = config["sdir"] + "/sim_out/sim_SE_{serl}_out/highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf}_pro.kraken",
		kreport = config["sdir"] + "/sim_out/sim_SE_{serl}_out/highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf}_pro.kreport"
	params:
		dbdir = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_pro",
		conf = "{conf}"
	threads: config["kraken_classify_threads"]
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	log:
		config["sdir"] + "/logs/highres_classify_sim_SE_rl{serl}_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf}_pro.log"
	shell:
		"""
		/usr/bin/time -v kraken2 --db {params.dbdir} --threads {threads} --confidence {params.conf} --report {output.kreport} --output {output.kraken} --report-minimizer-data {input.fq} &>> {log}
		"""

rule conifer_SE_highres_nolim_pro:
	input:
		kraken = config["sdir"] + "/sim_out/sim_SE_{serl}_out/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf}_pro.kraken",
		taxo = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_nolim_pro/taxo.k2d"
	output:
		conifer = config["sdir"] + "/sim_out/sim_SE_{serl}_out/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf}_pro.conifer_out"
	shell:
		"""
		conifer --all --both_scores -i {input.kraken} -d {input.taxo} > {output.conifer}
		"""

rule conifer_SE_highres_mem_pro:
	input:
		kraken = config["sdir"] + "/sim_out/sim_SE_{serl}_out/highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf}_pro.kraken",
		taxo = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_pro/taxo.k2d"
	output:
		conifer = config["sdir"] + "/sim_out/sim_SE_{serl}_out/highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf}_pro.conifer_out"
	shell:
		"""
		conifer --all --both_scores -i {input.kraken} -d {input.taxo} > {output.conifer}
		"""

rule confusion_matrix_SE_highres_nolim_pro:
	input:
		conifer = config["sdir"] + "/sim_out/sim_SE_{serl}_out/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf}_pro.conifer_out",
		coarse = config["sdir"] + "/sim_out/sim_SE_{serl}_out/" + config["best_coarse_SE"] + ".conifer_out",
		names = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_nolim_pro/taxonomy/names.dmp",
		nodes = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_nolim_pro/taxonomy/nodes.dmp",
	output:
		stats_taxlevels = config["sdir"] + "/sim_out/sim_SE_{serl}_out/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf}_rtl{rtl}_pro_stats_taxlevels.txt",
		stats_domain = config["sdir"] + "/sim_out/sim_SE_{serl}_out/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf}_rtl{rtl}_pro_stats_domain.txt",
		cm_taxlevels = config["sdir"] + "/sim_out/sim_SE_{serl}_out/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf}_rtl{rtl}_pro_cm_taxlevels.txt",
		cm_domain = config["sdir"] + "/sim_out/sim_SE_{serl}_out/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf}_rtl{rtl}_pro_cm_domain.txt"
	params:
		taxdir = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_nolim_pro/taxonomy",
		script = config["wdir"] + "/scripts/sim_results_confusion_highres_pro.R",
		rtl = "{rtl}",
		prefix = config["sdir"] + "/sim_out/sim_SE_{serl}_out/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf}_rtl{rtl}_pro"
	conda:
		config["wdir"] + "/envs/r.yaml"
	shell:
		"""
		{params.script} -i {input.conifer} -c {input.coarse} -t {params.taxdir} -s "{params.taxdir}/accessionTaxa.sql" -r {params.rtl} -m "SE" -o {params.prefix}
		"""

rule confusion_matrix_SE_highres_mem_pro:
	input:
		conifer = config["sdir"] + "/sim_out/sim_SE_{serl}_out/highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf}_pro.conifer_out",
		coarse = config["sdir"] + "/sim_out/sim_SE_{serl}_out/" + config["best_coarse_SE"] + ".conifer_out",
		names = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_pro/taxonomy/names.dmp",
		nodes = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_pro/taxonomy/nodes.dmp",
	output:
		stats_taxlevels = config["sdir"] + "/sim_out/sim_SE_{serl}_out/highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf}_rtl{rtl}_pro_stats_taxlevels.txt",
		stats_domain = config["sdir"] + "/sim_out/sim_SE_{serl}_out/highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf}_rtl{rtl}_pro_stats_domain.txt",
		cm_taxlevels = config["sdir"] + "/sim_out/sim_SE_{serl}_out/highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf}_rtl{rtl}_pro_cm_taxlevels.txt",
		cm_domain = config["sdir"] + "/sim_out/sim_SE_{serl}_out/highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf}_rtl{rtl}_pro_cm_domain.txt"
	params:
		taxdir = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_pro/taxonomy",
		script = config["wdir"] + "/scripts/sim_results_confusion_highres_pro.R",
		rtl = "{rtl}",
		prefix = config["sdir"] + "/sim_out/sim_SE_{serl}_out/highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf}_rtl{rtl}_pro"
	conda:
		config["wdir"] + "/envs/r.yaml"
	shell:
		"""
		{params.script} -i {input.conifer} -c {input.coarse} -t {params.taxdir} -s "{params.taxdir}/accessionTaxa.sql" -r {params.rtl} -m "SE" -o {params.prefix}
		"""


