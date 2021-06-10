rule extract_clustered_protein_accessions:
	input:
		acc2taxid = config["rdir"] + "/tax_combined/prot_accession2taxid.txt",
		protein_clustered = config["rdir"] + "/proteins_clustered/proteins_clustered.faa"
	output:
		acc2taxid_clustered = config["rdir"] + "/kaiju_db/kaiju_select_accession2taxid.txt"
	params:
		script = config["wdir"] + "/scripts/filter_acc2taxid.R",
		outdir = config["rdir"] + "/kaiju_db"
	conda:
		config["wdir"] + "/envs/r.yaml"
	threads: config["parallel_threads"]
	shell:
		"""
		mkdir -p {params.outdir}
		grep '^>' {input.protein_clustered} | sed -e 's/^>//' -e 's/ $//' > {params.outdir}/tmp
		{params.script} -i {params.outdir}/tmp -a {input.acc2taxid} -o {output.acc2taxid_clustered} -c {threads}
		# rm {params.outdir}/tmp
		"""

rule format_faa_kaiju:
	input:
		acc2taxid_clustered = config["rdir"] + "/kaiju_db/kaiju_select_accession2taxid.txt",
		protein_clustered = config["rdir"] + "/proteins_clustered/proteins_clustered.faa"
	output:
		protein_formatted = config["rdir"] + "/kaiju_db/library/proteins.faa",
		protein_accmap = config["rdir"] + "/kaiju_db/library/proteins_kaiju_accession_map.txt"
	params:
		script = config["wdir"] + "/scripts/format_kaiju.R",
		libdir = config["rdir"] + "/kaiju_db/library"
	conda:
		config["wdir"] + "/envs/r.yaml"
	threads: config["parallel_threads"]
	shell:
		"""
		{params.script} -i {input.protein_clustered} -t {input.acc2taxid_clustered} -c {threads} -o "{params.libdir}/tmp.faa" -m {output.protein_accmap}
		seqkit seq -w 0 -j {threads} "{params.libdir}/tmp.faa" | tr 'BZ' 'DE' | sed '/^>/! s/[^ARNDCQEGHILKMFPSTWYV]//g' > {output.protein_formatted}
		rm "{params.libdir}/tmp.faa"
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
		kaiju-mkbwt -n {threads} -a ACDEFGHIKLMNPQRSTVWY -o "{params.dbdir}/proteins" {input.faa} &>> {log}
		kaiju-mkfmi "{params.dbdir}/proteins" &>> {log}
		"""

