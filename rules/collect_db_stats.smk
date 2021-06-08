rule stats_genomes_derep:
	input:
		gen2taxid = config["rdir"] + "/tax_combined/full_genome2taxid.txt",
		sql = config["rdir"] + "/DB_taxonomy/accessionTaxa.sql"
	output:
		genome_stats = config["wdir"] + "/stats/db_stats_genomes_derep.txt",
		taxid = config["rdir"] + "/tmp/tmp_derep_genomes2taxid.txt"
	params:
		script = config["wdir"] + "/scripts/tax_db_stats.R",
		gendir = config["rdir"] + "/derep_combined"
	conda:
		config["wdir"] + "/envs/r.yaml"
	shell:
		"""
		find {params.gendir} -name '*.gz' | grep -o -F -f <(cut -f1 {input.gen2taxid}) | grep -w -F -f - {input.gen2taxid} | cut -f1,2 > {output.taxid}
		{params.script} -t {output.taxid} -s {input.sql} -o {output.genome_stats}
		"""

rule stats_kraken_db:
	input:
		kraken2_select = config["rdir"] + "/" + config["db_name"] + "/kraken2_select_accessions.txt",
		sql = config["rdir"] + "/DB_taxonomy/accessionTaxa.sql"
	output:
		kraken_stats = config["wdir"] + "/stats/db_stats_" + config["db_name"] + ".txt",
		taxid = config["rdir"] + "/tmp/tmp_"  + config["db_name"] + "_genomes2taxid.txt"
	params:
		script = config["wdir"] + "/scripts/tax_db_stats.R"
	conda:
		config["wdir"] + "/envs/r.yaml"
	shell:
		"""
		cut -f1,3 {input.kraken2_select} > {output.taxid}
		{params.script} -t {output.taxid} -s {input.sql} -o {output.kraken_stats}
		"""

rule stats_proteins_all:
	input:
		acc2taxid = config["rdir"] + "/tax_combined/prot_accession2taxid.txt",
		sql = config["rdir"] + "/DB_taxonomy/accessionTaxa.sql"
	output:
		protein_stats = config["wdir"] + "/stats/db_stats_proteins_all.txt"
	params:
		script = config["wdir"] + "/scripts/tax_db_stats.R"
	conda:
		config["wdir"] + "/envs/r.yaml"
	shell:
		"""
		{params.script} -t {input.acc2taxid} -s {input.sql} -o {output.protein_stats}
		"""

rule stats_proteins_clustered:
	input:
		kaiju_select = config["rdir"] + "/kaiju_db/kaiju_select_accession2taxid.txt",
		sql = config["rdir"] + "/DB_taxonomy/accessionTaxa.sql"
	output:
		clustered_stats = config["wdir"] + "/stats/db_stats_proteins_clustered.txt"
	params:
		script = config["wdir"] + "/scripts/tax_db_stats.R"
	conda:
		config["wdir"] + "/envs/r.yaml"
	shell:
		"""
		{params.script} -t {input.kaiju_select} -s {input.sql} -o {output.clustered_stats}
		"""
