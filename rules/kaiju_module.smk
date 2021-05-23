rule extract_clustered_protein_accessions:
	input:
		acc2taxid = config["rdir"] + "/tax_combined/prot_accession2taxid.txt",
		protein_clustered = config["rdir"] + "/proteins_clustered/proteins_clustered.faa"
	output:
		acc2taxid_clustered = config["rdir"] + "/kaiju_db/kaiju_select_accession2taxid.txt"
	params:
		script = config["wdir"] + "/scripts/subset_prot_accession2taxid.R",
		outdir = config["rdir"] + "/kaiju_db"
	conda:
		config["wdir"] + "/envs/r.yaml"
	threads: config["parallel_threads"]
	shell:
		"""
		mkdir -p {params.outdir}
		grep '^>' {input.protein_clustered} | sed -e 's/^>//' -e 's/\s.*$//' > {params.outdir}/tmp
		{input.acc2taxid} > {output.acc2taxid_clustered}
		rm {params.outdir}/tmp
		"""

rule format_faa_kaiju:
	input:
		acc2taxid_clustered = config["rdir"] + "/kaiju_db/kaiju_select_accession2taxid.txt",
		protein_clustered = config["rdir"] + "/proteins_clustered/proteins_clustered.faa"
	output:
		protein_formatted = config["rdir"] + "/kaiju_db/library/proteins.faa"
	params:
		script = config["wdir"] + "/scripts/format_kaiju.sh"
	conda:
		config["wdir"] + "/envs/parallel.yaml"
	threads: config["parallel_threads"]
	shell:
		"""
		cut -f2 {input.acc2taxid_clustered} | sort | uniq | parallel -j{threads} '{params.script} {{}} {input.acc2taxid_clustered} {input.protein_clustered}' >> {output.protein_formatted}
		"""

# optional: run conterminator on protein sequences: TO BE IMPLEMENTED!!!

rule build_kaiju_db:
	input:
		faa = config["rdir"] + "/kaiju_db/library/proteins.faa"
	output:
		kaiju = config["rdir"] + "/kaiju_db/proteins.fmi"
	params:
		dbdir = config["rdir"] + "/kaiju_db"
	threads: config["kaiju_threads"]
	conda:
		config["wdir"] + "/envs/kaiju.yaml"
	log:
		config["rdir"] + "/logs/build_kaiju.log"
	shell:
		"""
		kaiju-mkbwt -n {threads} -a ACDEFGHIKLMNPQRSTVWY -o "{params.dbdir}/proteins" {input.faa}
		kaiju-mkfmi "{params.dbdir}/proteins"
		"""

