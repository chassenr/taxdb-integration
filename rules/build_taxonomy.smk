rule check_taxpath:
	input:
		gtdb = config["rdir"] + "/tax_combined/gtdb_derep_taxonomy.txt",
		# this is to make sure that I can add faa reps for taxa if no protein otherwise available (maybe not necessary anymore once prodigal integrated)
		gtdb_reps = config["rdir"] + "/gtdb/metadata/gtdb_reps_tax.txt",
		checkv = config["rdir"] + "/tax_combined/checkv_derep_taxonomy.txt",
		# add all reps also for checkv
		checkv_reps = config["rdir"] + "/checkv/checkv_reps_taxonomy.txt",
		ncbi = expand(config["rdir"] + "/tax_combined/{library_name}_derep_taxonomy.txt", library_name = LIBRARY_NAME),
		added_nuc_ncbi = config["rdir"] + "/tax_combined/euk_custom_post_derep_taxonomy.txt" if config["custom_ncbi_post_derep"] else [],
		added_nuc_gtdb = config["rdir"] + "/tax_combined/pro_custom_post_derep_taxonomy.txt" if config["custom_gtdb_post_derep"] else [],
		added_nuc_checkv = config["rdir"] + "/tax_combined/vir_custom_post_derep_taxonomy.txt" if config["custom_checkv_post_derep"] else [],
		added_prot_ncbi = config["custom_euk_prot"] if config["custom_euk_prot"] else [],
		added_prot_gtdb = config["custom_pro_prot"] if config["custom_pro_prot"] else [],
		added_prot_checkv = config["custom_vir_prot"] if config["custom_vir_prot"] else []
	output:
		tax_combined = config["rdir"] + "/tax_combined/full_taxonomy_combined.txt"
	params:
		outdir = config["rdir"] + "/tax_combined",
		script = config["wdir"] + "/scripts/fix_ncbi_taxpath.R"
	conda:
		config["wdir"] + "/envs/r.yaml"
	log:
		config["rdir"] + "/logs/fix_ncbi_taxpath.log"
	shell:
		"""
		cat {input.gtdb} {input.gtdb_reps} {input.checkv} {input.checkv_reps} {input.ncbi} {input.added_nuc_ncbi} {input.added_nuc_gtdb} {input.added_nuc_checkv} > "{params.outdir}/tmp"
		cat {input.added_prot_ncbi} {input.added_prot_gtdb} {input.added_prot_checkv} | cut -f1,2 >> "{params.outdir}/tmp"
		{params.script} -i "{params.outdir}/tmp" -o "{output.tax_combined}" &>> {log}
		rm "{params.outdir}/tmp"
		"""

rule build_taxonomy:
	input:
		tax_combined = config["rdir"] + "/tax_combined/full_taxonomy_combined.txt"
	output:
		tax_good = config["rdir"] + "/tax_combined/full_taxonomy_good.txt",
		nodes = config["rdir"] + "/DB_taxonomy/nodes.dmp",
		names = config["rdir"] + "/DB_taxonomy/names.dmp"
	params:
		script = config["gtdb_to_taxdump"],
		taxdir = config["rdir"] + "/DB_taxonomy"
	conda:
		config["wdir"] + "/envs/biopython.yaml"
	log:
		config["rdir"] + "/logs/build_common_taxonomy.log"
	shell:
		"""
		cut -f1,2 {input.tax_combined} > {output.tax_good}
		{params.script} -o {params.taxdir} {output.tax_good} &>> {log}
		sed '/\#/d' -i {params.taxdir}/merged.dmp
		sed '/\#/d' -i {params.taxdir}/delnodes.dmp 
		"""

# this needs to be collected separately to use the fixed taxpath and exclude accessions without proteins
rule full_protein_taxonomy:
	input:
		config["rdir"] + "/tax_combined/euk_custom_protein_taxonomy.txt" if config["custom_euk_prot"] else [],
		expand(config["rdir"] + "/tax_combined/{library_name}_protein_taxonomy.txt", library_name = LIBRARY_NAME),
		config["rdir"] + "/tax_combined/pro_custom_protein_taxonomy.txt" if config["custom_pro_prot"] else [],
		config["rdir"] + "/tax_combined/gtdb_protein_taxonomy.txt",
		config["rdir"] + "/tax_combined/vir_custom_protein_taxonomy.txt" if config["custom_vir_prot"] else [],
		config["rdir"] + "/tax_combined/checkv_protein_taxonomy.txt"
	output:
		config["rdir"] + "/tax_combined/protein_taxonomy_good.txt"
	shell:
		"""
		cat {input} > {output}
		"""

rule generate_genome2taxid:
	input:
		tax_good = config["rdir"] + "/tax_combined/full_taxonomy_good.txt"
	output:
		gen2taxid = config["rdir"] + "/tax_combined/full_genome2taxid.txt"
	params:
		dbdir = config["rdir"] + "/DB_taxonomy"
	conda:
		config["wdir"] + "/envs/taxonkit.yaml"
	shell:
		"""
		cut -f1 {input.tax_good} | taxonkit name2taxid --data-dir {params.dbdir} > "{params.dbdir}/tmp_taxid.txt"
		cut -f2 "{params.dbdir}/tmp_taxid.txt" | taxonkit lineage -t --data-dir {params.dbdir} | cut -f3 | cut -d';' -f7 | paste "{params.dbdir}/tmp_taxid.txt" - > {output.gen2taxid}
		rm "{params.dbdir}/tmp_taxid.txt"
		"""

rule generate_nucl_accesion2taxid:
	input:
		gen2taxid = config["rdir"] + "/tax_combined/full_genome2taxid.txt"
	output:
		acc2taxid = config["rdir"] + "/tax_combined/nucl_accession2taxid.txt"
	params:
		script = config["wdir"] + "/scripts/generate_accession2taxid.sh",
		gendir = config["rdir"] + "/derep_combined"
	conda:
		config["wdir"] + "/envs/parallel.yaml"
	threads: config["taxmap_threads"]
	shell:
		"""
		cut -f1,3 {input.gen2taxid} | parallel -j{threads} '{params.script} {{}} {params.gendir}' >> {output.acc2taxid}
		"""

rule generate_prot_accesion2taxid:
	input:
		gen2taxid = config["rdir"] + "/tax_combined/full_genome2taxid.txt"
	output:
		acc2taxid = config["rdir"] + "/tax_combined/prot_accession2taxid.txt"
	params:
		script = config["wdir"] + "/scripts/generate_accession2taxid.sh",
		gendir = config["rdir"] + "/proteins_all"
	conda:
		config["wdir"] + "/envs/parallel.yaml"
	shell:
		"""
		cut -f1,3 {input.gen2taxid} | parallel -j{threads} '{params.script} {{}} {params.gendir}' >> {output.acc2taxid}
		"""

