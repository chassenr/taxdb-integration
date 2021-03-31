rule ps_classify_PE_highres_nolim_pro:
	input:
		hash = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_nolim_pro/hash.k2d",
		fq_R1 = config["sdir"] + "/sim_out/sim_PE_{rl}_out/1/R1.fq",
		fq_R2 = config["sdir"] + "/sim_out/sim_PE_{rl}_out/1/R2.fq"
	output:
		kraken = config["sdir"] + "/sim_out/sim_PE_{rl}_out/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf_highres}_pro.kraken",
		kreport = config["sdir"] + "/sim_out/sim_PE_{rl}_out/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf_highres}_pro.kreport"
	params:
		dbdir = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_nolim_pro",
		conf = "{conf_highres}"
	threads: config["kraken_classify_threads"]
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	log:
		config["sdir"] + "/logs/highres_classify_sim_PE_rl{rl}_k{kl}_m{ml}_s{ms}_nolim_c{conf_highres}_pro.log"
	shell:
		"""
		/usr/bin/time -v kraken2 --db {params.dbdir} --threads {threads} --confidence {params.conf} --report {output.kreport} --output {output.kraken} --report-minimizer-data --paired {input.fq_R1} {input.fq_R2} &>> {log}
		"""

#rule ps_classify_PE_highres_mem_pro:
#	input:
#		hash = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_pro/hash.k2d",
#		fq_R1 = config["sdir"] + "/sim_out/sim_PE_{rl}_out/1/R1.fq",
#		fq_R2 = config["sdir"] + "/sim_out/sim_PE_{rl}_out/1/R2.fq"
#	output:
#		kraken = config["sdir"] + "/sim_out/sim_PE_{rl}_out/highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf_highres}_pro.kraken",
#		kreport = config["sdir"] + "/sim_out/sim_PE_{rl}_out/highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf_highres}_pro.kreport"
#	params:
#		dbdir = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_pro",
#		conf = "{conf_highres}"
#	threads: config["kraken_classify_threads"]
#	conda:
#		config["wdir"] + "/envs/kraken2.yaml"
#	log:
#		config["sdir"] + "/logs/highres_classify_sim_PE_rl{rl}_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf_highres}_pro.log"
#	shell:
#		"""
#		/usr/bin/time -v kraken2 --db {params.dbdir} --threads {threads} --confidence {params.conf} --report {output.kreport} --output {output.kraken} --report-minimizer-data --paired {input.fq_R1} {input.fq_R2} &>> {log}
#		"""

rule conifer_PE_highres_nolim_pro:
	input:
		kraken = config["sdir"] + "/sim_out/sim_PE_{rl}_out/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf_highres}_pro.kraken",
		taxo = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_nolim_pro/taxo.k2d"
	output:
		conifer = config["sdir"] + "/sim_out/sim_PE_{rl}_out/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf_highres}_pro.conifer_out"
	shell:
		"""
		conifer --all --both_scores -i {input.kraken} -d {input.taxo} > {output.conifer}
		"""

#rule conifer_PE_highres_mem_pro:
#	input:
#		kraken = config["sdir"] + "/sim_out/sim_PE_{rl}_out/highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf_highres}_pro.kraken",
#		taxo = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_pro/taxo.k2d"
#	output:
#		conifer = config["sdir"] + "/sim_out/sim_PE_{rl}_out/highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf_highres}_pro.conifer_out"
#	shell:
#		"""
#		conifer --all --both_scores -i {input.kraken} -d {input.taxo} > {output.conifer}
#		"""

rule confusion_matrix_PE_highres_nolim_pro:
	input:
		conifer = config["sdir"] + "/sim_out/sim_PE_{rl}_out/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf_highres}_pro.conifer_out",
		coarse = config["sdir"] + "/sim_out/sim_PE_{rl}_out/coarse_k{kl}_m{ml}_s{ms}_nolim_c{conf_coarse}.conifer_out",
		names_highres = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_nolim_pro/taxonomy/names.dmp",
		nodes_highres = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_nolim_pro/taxonomy/nodes.dmp",
		names_coarse = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_nolim/taxonomy/names.dmp",
		nodes_coarse = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_nolim/taxonomy/nodes.dmp"
	output:
		cm_domain_filt = config["sdir"] + "/sim_out/sim_PE_{rl}_out/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf_highres}_rtl{rtl_highres}_pro_coarse_c{conf_coarse}_rtl{rtl_coarse}_cm_domain_filt.txt",
		cm_domain_nofilt = config["sdir"] + "/sim_out/sim_PE_{rl}_out/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf_highres}_rtl{rtl_highres}_pro_coarse_c{conf_coarse}_rtl{rtl_coarse}_cm_domain_nofilt.txt",
		stats_domain_filt = config["sdir"] + "/sim_out/sim_PE_{rl}_out/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf_highres}_rtl{rtl_highres}_pro_coarse_c{conf_coarse}_rtl{rtl_coarse}_stats_domain_filt.txt",
		stats_domain_nofilt = config["sdir"] + "/sim_out/sim_PE_{rl}_out/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf_highres}_rtl{rtl_highres}_pro_coarse_c{conf_coarse}_rtl{rtl_coarse}_stats_domain_nofilt.txt",
		summary_taxlevels_filt = config["sdir"] + "/sim_out/sim_PE_{rl}_out/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf_highres}_rtl{rtl_highres}_pro_coarse_c{conf_coarse}_rtl{rtl_coarse}_summary_taxlevels_filt.txt",
		summary_taxlevels_nofilt = config["sdir"] + "/sim_out/sim_PE_{rl}_out/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf_highres}_rtl{rtl_highres}_pro_coarse_c{conf_coarse}_rtl{rtl_coarse}_summary_taxlevels_nofilt.txt"
	params:
		taxdir_highres = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_nolim_pro/taxonomy",
		taxdir_coarse = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_nolim/taxonomy",
		script = config["wdir"] + "/scripts/sim_results_confusion_highres_pro.R",
		rtl_highres = "{rtl_highres}",
		rtl_coarse = "{rtl_coarse}",
		prefix = config["sdir"] + "/sim_out/sim_PE_{rl}_out/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf_highres}_rtl{rtl_highres}_pro_coarse_c{conf_coarse}_rtl{rtl_coarse}"
	conda:
		config["wdir"] + "/envs/r.yaml"
	shell:
		"""
		{params.script} -i {input.conifer} -t {params.taxdir_highres} -s "{params.taxdir_highres}/accessionTaxa.sql" -r {params.rtl_highres} -m "PE" -c {input.coarse} -T {params.taxdir_coarse} -S "{params.taxdir_coarse}/accessionTaxa.sql" -R {params.rtl_coarse} -o {params.prefix}
		"""

# as the highres DBs are all below the the mem limit, use nolim here in combination with mem limited coarse
rule confusion_matrix_PE_highres_mem_pro:
	input:
		conifer = config["sdir"] + "/sim_out/sim_PE_{rl}_out/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf_highres}_pro.conifer_out",
		coarse = config["sdir"] + "/sim_out/sim_PE_{rl}_out/coarse_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf_coarse}.conifer_out",
		names_highres = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_nolim_pro/taxonomy/names.dmp",
		nodes_highres = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_nolim_pro/taxonomy/nodes.dmp",
		names_coarse = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_mem{dbsize}/taxonomy/names.dmp",
		nodes_coarse = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_mem{dbsize}/taxonomy/nodes.dmp"
	output:
		cm_domain_filt = config["sdir"] + "/sim_out/sim_PE_{rl}_out/highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf_highres}_rtl{rtl_highres}_pro_coarse_c{conf_coarse}_rtl{rtl_coarse}_cm_domain_filt.txt",
		cm_domain_nofilt = config["sdir"] + "/sim_out/sim_PE_{rl}_out/highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf_highres}_rtl{rtl_highres}_pro_coarse_c{conf_coarse}_rtl{rtl_coarse}_cm_domain_nofilt.txt",
		stats_domain_filt = config["sdir"] + "/sim_out/sim_PE_{rl}_out/highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf_highres}_rtl{rtl_highres}_pro_coarse_c{conf_coarse}_rtl{rtl_coarse}_stats_domain_filt.txt",
		stats_domain_nofilt = config["sdir"] + "/sim_out/sim_PE_{rl}_out/highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf_highres}_rtl{rtl_highres}_pro_coarse_c{conf_coarse}_rtl{rtl_coarse}_stats_domain_nofilt.txt",
		summary_taxlevels_filt = config["sdir"] + "/sim_out/sim_PE_{rl}_out/highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf_highres}_rtl{rtl_highres}_pro_coarse_c{conf_coarse}_rtl{rtl_coarse}_summary_taxlevels_filt.txt",
		summary_taxlevels_nofilt = config["sdir"] + "/sim_out/sim_PE_{rl}_out/highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf_highres}_rtl{rtl_highres}_pro_coarse_c{conf_coarse}_rtl{rtl_coarse}_summary_taxlevels_nofilt.txt"
	params:
		taxdir_highres = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_nolim_pro/taxonomy",
		taxdir_coarse = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_mem{dbsize}/taxonomy",
		script = config["wdir"] + "/scripts/sim_results_confusion_highres_pro.R",
		rtl_highres = "{rtl_highres}",
		rtl_coarse = "{rtl_coarse}",
		prefix = config["sdir"] + "/sim_out/sim_PE_{rl}_out/highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf_highres}_rtl{rtl_highres}_pro_coarse_c{conf_coarse}_rtl{rtl_coarse}"
	conda:
		config["wdir"] + "/envs/r.yaml"
	shell:
		"""
		{params.script} -i {input.conifer} -t {params.taxdir_highres} -s "{params.taxdir_highres}/accessionTaxa.sql" -r {params.rtl_highres} -m "PE" -c {input.coarse} -T {params.taxdir_coarse} -S "{params.taxdir_coarse}/accessionTaxa.sql" -R {params.rtl_coarse} -o {params.prefix}
		"""

rule ps_classify_SE_highres_nolim_pro:
	input:
		hash = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_nolim_pro/hash.k2d",
		fq = config["sdir"] + "/sim_out/sim_SE_{serl}_out/1/R1.fq"
	output:
		kraken = config["sdir"] + "/sim_out/sim_SE_{serl}_out/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf_highres}_pro.kraken",
		kreport = config["sdir"] + "/sim_out/sim_SE_{serl}_out/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf_highres}_pro.kreport"
	params:
		dbdir = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_nolim_pro",
		conf = "{conf_highres}"
	threads: config["kraken_classify_threads"]
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	log:
		config["sdir"] + "/logs/highres_classify_sim_SE_rl{serl}_k{kl}_m{ml}_s{ms}_nolim_c{conf_highres}_pro.log"
	shell:
		"""
		/usr/bin/time -v kraken2 --db {params.dbdir} --threads {threads} --confidence {params.conf} --report {output.kreport} --output {output.kraken} --report-minimizer-data {input.fq} &>> {log}
		"""

#rule ps_classify_SE_highres_mem_pro:
#	input:
#		hash = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_pro/hash.k2d",
#		fq = config["sdir"] + "/sim_out/sim_SE_{serl}_out/1/R1.fq"
#	output:
#		kraken = config["sdir"] + "/sim_out/sim_SE_{serl}_out/highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf_highres}_pro.kraken",
#		kreport = config["sdir"] + "/sim_out/sim_SE_{serl}_out/highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf_highres}_pro.kreport"
#	params:
#		dbdir = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_pro",
#		conf = "{conf_highres}"
#	threads: config["kraken_classify_threads"]
#	conda:
#		config["wdir"] + "/envs/kraken2.yaml"
#	log:
#		config["sdir"] + "/logs/highres_classify_sim_SE_rl{serl}_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf_highres}_pro.log"
#	shell:
#		"""
#		/usr/bin/time -v kraken2 --db {params.dbdir} --threads {threads} --confidence {params.conf} --report {output.kreport} --output {output.kraken} --report-minimizer-data {input.fq} &>> {log}
#		"""

rule conifer_SE_highres_nolim_pro:
	input:
		kraken = config["sdir"] + "/sim_out/sim_SE_{serl}_out/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf_highres}_pro.kraken",
		taxo = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_nolim_pro/taxo.k2d"
	output:
		conifer = config["sdir"] + "/sim_out/sim_SE_{serl}_out/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf_highres}_pro.conifer_out"
	shell:
		"""
		conifer --all --both_scores -i {input.kraken} -d {input.taxo} > {output.conifer}
		"""

#rule conifer_SE_highres_mem_pro:
#	input:
#		kraken = config["sdir"] + "/sim_out/sim_SE_{serl}_out/highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf_highres}_pro.kraken",
#		taxo = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_pro/taxo.k2d"
#	output:
#		conifer = config["sdir"] + "/sim_out/sim_SE_{serl}_out/highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf_highres}_pro.conifer_out"
#	shell:
#		"""
#		conifer --all --both_scores -i {input.kraken} -d {input.taxo} > {output.conifer}
#		"""

rule confusion_matrix_SE_highres_nolim_pro:
	input:
		conifer = config["sdir"] + "/sim_out/sim_SE_{serl}_out/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf_highres}_pro.conifer_out",
		coarse = config["sdir"] + "/sim_out/sim_SE_{serl}_out/coarse_k{kl}_m{ml}_s{ms}_nolim_c{conf_coarse}.conifer_out",
		names_highres = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_nolim_pro/taxonomy/names.dmp",
		nodes_highres = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_nolim_pro/taxonomy/nodes.dmp",
		names_coarse = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_nolim/taxonomy/names.dmp",
		nodes_coarse = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_nolim/taxonomy/nodes.dmp"
	output:
		cm_domain_filt = config["sdir"] + "/sim_out/sim_SE_{serl}_out/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf_highres}_rtl{rtl_highres}_pro_coarse_c{conf_coarse}_rtl{rtl_coarse}_cm_domain_filt.txt",
		cm_domain_nofilt = config["sdir"] + "/sim_out/sim_SE_{serl}_out/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf_highres}_rtl{rtl_highres}_pro_coarse_c{conf_coarse}_rtl{rtl_coarse}_cm_domain_nofilt.txt",
		stats_domain_filt = config["sdir"] + "/sim_out/sim_SE_{serl}_out/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf_highres}_rtl{rtl_highres}_pro_coarse_c{conf_coarse}_rtl{rtl_coarse}_stats_domain_filt.txt",
		stats_domain_nofilt = config["sdir"] + "/sim_out/sim_SE_{serl}_out/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf_highres}_rtl{rtl_highres}_pro_coarse_c{conf_coarse}_rtl{rtl_coarse}_stats_domain_nofilt.txt",
		summary_taxlevels_filt = config["sdir"] + "/sim_out/sim_SE_{serl}_out/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf_highres}_rtl{rtl_highres}_pro_coarse_c{conf_coarse}_rtl{rtl_coarse}_summary_taxlevels_filt.txt",
		summary_taxlevels_nofilt = config["sdir"] + "/sim_out/sim_SE_{serl}_out/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf_highres}_rtl{rtl_highres}_pro_coarse_c{conf_coarse}_rtl{rtl_coarse}_summary_taxlevels_nofilt.txt"
	params:
		taxdir_highres = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_nolim_pro/taxonomy",
		taxdir_coarse = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_nolim/taxonomy",
		script = config["wdir"] + "/scripts/sim_results_confusion_highres_pro.R",
		rtl_highres = "{rtl_highres}",
		rtl_coarse = "{rtl_coarse}",
		prefix = config["sdir"] + "/sim_out/sim_SE_{serl}_out/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf_highres}_rtl{rtl_highres}_pro_coarse_c{conf_coarse}_rtl{rtl_coarse}"
	conda:
		config["wdir"] + "/envs/r.yaml"
	shell:
		"""
		{params.script} -i {input.conifer} -t {params.taxdir_highres} -s "{params.taxdir_highres}/accessionTaxa.sql" -r {params.rtl_highres} -m "SE" -c {input.coarse} -T {params.taxdir_coarse} -S "{params.taxdir_coarse}/accessionTaxa.sql" -R {params.rtl_coarse} -o {params.prefix}
		"""

# as the highres DBs are all below the the mem limit, use nolim here in combination with mem limited coarse
rule confusion_matrix_SE_highres_mem_pro:
	input:
		conifer = config["sdir"] + "/sim_out/sim_SE_{serl}_out/highres_k{kl}_m{ml}_s{ms}_nolim_c{conf_highres}_pro.conifer_out",
		coarse = config["sdir"] + "/sim_out/sim_SE_{serl}_out/coarse_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf_coarse}.conifer_out",
		names_highres = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_nolim_pro/taxonomy/names.dmp",
		nodes_highres = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_nolim_pro/taxonomy/nodes.dmp",
		names_coarse = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_mem{dbsize}/taxonomy/names.dmp",
		nodes_coarse = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_mem{dbsize}/taxonomy/nodes.dmp"
	output:
		cm_domain_filt = config["sdir"] + "/sim_out/sim_SE_{serl}_out/highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf_highres}_rtl{rtl_highres}_pro_coarse_c{conf_coarse}_rtl{rtl_coarse}_cm_domain_filt.txt",
		cm_domain_nofilt = config["sdir"] + "/sim_out/sim_SE_{serl}_out/highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf_highres}_rtl{rtl_highres}_pro_coarse_c{conf_coarse}_rtl{rtl_coarse}_cm_domain_nofilt.txt",
		stats_domain_filt = config["sdir"] + "/sim_out/sim_SE_{serl}_out/highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf_highres}_rtl{rtl_highres}_pro_coarse_c{conf_coarse}_rtl{rtl_coarse}_stats_domain_filt.txt",
		stats_domain_nofilt = config["sdir"] + "/sim_out/sim_SE_{serl}_out/highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf_highres}_rtl{rtl_highres}_pro_coarse_c{conf_coarse}_rtl{rtl_coarse}_stats_domain_nofilt.txt",
		summary_taxlevels_filt = config["sdir"] + "/sim_out/sim_SE_{serl}_out/highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf_highres}_rtl{rtl_highres}_pro_coarse_c{conf_coarse}_rtl{rtl_coarse}_summary_taxlevels_filt.txt",
		summary_taxlevels_nofilt = config["sdir"] + "/sim_out/sim_SE_{serl}_out/highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf_highres}_rtl{rtl_highres}_pro_coarse_c{conf_coarse}_rtl{rtl_coarse}_summary_taxlevels_nofilt.txt"
	params:
		taxdir_highres = config["sdir"] + "/db_highres_k{kl}_m{ml}_s{ms}_nolim_pro/taxonomy",
		taxdir_coarse = config["sdir"] + "/db_coarse_k{kl}_m{ml}_s{ms}_mem{dbsize}/taxonomy",
		script = config["wdir"] + "/scripts/sim_results_confusion_highres_pro.R",
		rtl_highres = "{rtl_highres}",
		rtl_coarse = "{rtl_coarse}",
		prefix = config["sdir"] + "/sim_out/sim_SE_{serl}_out/highres_k{kl}_m{ml}_s{ms}_mem{dbsize}_c{conf_highres}_rtl{rtl_highres}_pro_coarse_c{conf_coarse}_rtl{rtl_coarse}"
	conda:
		config["wdir"] + "/envs/r.yaml"
	shell:
		"""
		{params.script} -i {input.conifer} -t {params.taxdir_highres} -s "{params.taxdir_highres}/accessionTaxa.sql" -r {params.rtl_highres} -m "SE" -c {input.coarse} -T {params.taxdir_coarse} -S "{params.taxdir_coarse}/accessionTaxa.sql" -R {params.rtl_coarse} -o {params.prefix}
		"""


