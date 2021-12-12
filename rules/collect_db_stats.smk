rule stats_genomes_derep_taxid:
	input:
		gen2taxid = config["rdir"] + "/tax_combined/full_genome2taxid.txt"
	output:
		tmp_taxid = config["rdir"] + "/tmp/tmp_derep_genomes_paths.txt"
	params:
		gendir = config["rdir"] + "/derep_combined",
		tmpdir = config["rdir"] + "/tmp",
		taxdir = config["rdir"] + "/DB_taxonomy"
	conda:
		config["wdir"] + "/envs/taxonkit.yaml"
	shell:
		"""
		find {params.gendir} -name '*.gz' | grep -o -F -f <(cut -f1 {input.gen2taxid}) | grep -w -F -f - {input.gen2taxid} | cut -f1,3 > {params.tmpdir}/tmp_derep_genomes_taxid.txt
		echo -e 'genome\\ttaxid\\tpath' > {output.tmp_taxid}
		cut -f2 {params.tmpdir}/tmp_derep_genomes_taxid.txt | taxonkit lineage --data-dir {params.taxdir} | cut -f2 | paste {params.tmpdir}/tmp_derep_genomes_taxid.txt - >> {output.tmp_taxid}
		rm {params.tmpdir}/tmp_derep_genomes_taxid.txt
		"""

rule stats_genomes_derep_summary:
	input:
		tmp_taxid = config["rdir"] + "/tmp/tmp_derep_genomes_paths.txt"
	output:
		genome_stats = config["wdir"] + "/stats/db_stats_genomes_derep.txt"
	params:
		script = config["wdir"] + "/scripts/tax_db_stats.R"
	conda:
		config["wdir"] + "/envs/r.yaml"
	shell:
		"""
		{params.script} -t {input.tmp_taxid} -o {output.genome_stats}
		"""

rule stats_kraken_db_taxid:
	input:
		kraken2_select = config["rdir"] + "/" + config["db_name"] + "/kraken2_select_accessions.txt"
	output:
		tmp_taxid = config["rdir"] + "/tmp/tmp_"  + config["db_name"] + "_paths.txt"
	params:
		taxdir = config["rdir"] + "/DB_taxonomy"
	conda:
		config["wdir"] + "/envs/taxonkit.yaml"
	shell:
		"""
		echo -e 'genome\\ttaxid\\tpath' > {output.tmp_taxid}
		cut -f3 {input.kraken2_select} | taxonkit lineage --data-dir {params.taxdir} | paste <(cut -f1 {input.kraken2_select}) - >> {output.tmp_taxid}
		"""

rule stats_kraken_db_summary:
	input:
		tmp_taxid = config["rdir"] + "/tmp/tmp_"  + config["db_name"] + "_paths.txt"
	output:
		kraken_stats = config["wdir"] + "/stats/db_stats_" + config["db_name"] + ".txt"
	params:
		script = config["wdir"] + "/scripts/tax_db_stats.R"
	conda:
		config["wdir"] + "/envs/r.yaml"
	shell:
		"""
		{params.script} -t {input.tmp_taxid} -o {output.kraken_stats}
		"""

rule stats_proteins_all_taxid:
	input:
		acc2taxid = config["rdir"] + "/tax_combined/prot_accession2taxid.txt"
	output:
		tmp_taxid = config["rdir"] + "/tmp/tmp_proteins_all_paths.txt"
	params:
		taxdir = config["rdir"] + "/DB_taxonomy"
	conda:
		config["wdir"] + "/envs/taxonkit.yaml"
	shell:
		"""
		echo -e 'genome\\tsequence\\ttaxid\\tpath' > {output.tmp_taxid}
		cut -f3 {input.acc2taxid} | taxonkit lineage --data-dir {params.taxdir} | cut -f2 | paste {input.acc2taxid} - >> {output.tmp_taxid}
		"""

rule stats_proteins_all_summary:
	input:
		tmp_taxid = config["rdir"] + "/tmp/tmp_proteins_all_paths.txt"
	output:
		protein_stats = config["wdir"] + "/stats/db_stats_proteins_all.txt"
	params:
		script = config["wdir"] + "/scripts/tax_db_stats.R"
	threads: config["fread_threads"]
	conda:
		config["wdir"] + "/envs/r.yaml"
	shell:
		"""
		{params.script} -t {input.tmp_taxid} -o {output.protein_stats} -c {threads}
		"""

rule stats_proteins_clustered_taxid:
	input:
		protein_accmap = config["rdir"] + "/kaiju_db/library/proteins_kaiju_accession_map.txt"
	output:
		tmp_taxid = config["rdir"] + "/tmp/tmp_proteins_clustered_paths.txt"
	params:
		taxdir = config["rdir"] + "/DB_taxonomy"
	conda:
		config["wdir"] + "/envs/taxonkit.yaml"
	shell:
		"""
		echo -e 'taxid\\tpath' > {output.tmp_taxid}
		cut -f1 {input.protein_accmap} | taxonkit lineage --data-dir {params.taxdir} >> {output.tmp_taxid}
		"""

rule stats_proteins_clustered_summary:
	input:
		tmp_taxid = config["rdir"] + "/tmp/tmp_proteins_clustered_paths.txt"
	output:
		clustered_stats = config["wdir"] + "/stats/db_stats_proteins_clustered.txt"
	params:
		script = config["wdir"] + "/scripts/tax_db_stats.R"
	threads: config["fread_threads"]
	conda:
		config["wdir"] + "/envs/r.yaml"
	shell:
		"""
		{params.script} -t {input.tmp_taxid} -o {output.clustered_stats} -c {threads}
		"""

