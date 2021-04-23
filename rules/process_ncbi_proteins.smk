rule download_proteins_ncbi:
	input:
		url = config["rdir"] + "/{library_name}/assembly_url_genomic.txt",
		tax = config["rdir"] + "/tax_combined/{library_name}_derep_taxonomy.txt",
		tax_good = config["rdir"] + "/tax_combined/full_taxonomy_good.txt"
	output:
		tax_prot = config["rdir"] + "/{library_name}/assembly_protein_taxonomy.txt"
	params:
		outdir = config["rdir"] + "/{library_name}/proteins"
	threads: config["download_threads"]
	conda:
		config["wdir"] + "/envs/download.yaml"
	log:
		config["rdir"] + "/logs/download_proteins_ncbi_{library_name}.log"
	shell:
		"""
		awk '$3 == "TRUE"' {input.url} | cut -f2 | grep -F -f <(cut -f1 {input.tax}) > "{params.outdir}/links.list"
		aria2c -i "{params.outdir}/links.list" -c -l "{params.outdir}/links.log" --dir {params.outdir} --max-tries=20 --retry-wait=5 --max-connection-per-server=1 --max-concurrent-downloads={threads} &>> {log}
		# We need to verify all files are there
		cat "{params.outdir}/links.list" | xargs -n 1 basename | sort > "{params.outdir}/tmp1"
		find {params.outdir} -type f -name '*.gz' | xargs -n 1 basename | sort > "{params.outdir}/tmp2"
		if diff "{params.outdir}/tmp1" "{params.outdir}/tmp2"
		then
		  cut -d'_' -f1,2 "{params.outdir}/tmp1" | grep -F -f - {input.tax_good} > {output.tax_prot}
		fi
		rm "{params.outdir}/links.list" "{params.outdir}/links.log" "{params.outdir}/tmp1" "{params.outdir}/tmp2"
		"""

rule custom_ncbi_proteins_pre_derep:
	input:
		tax = config["rdir"] + "/tax_combined/{library_name}_derep_taxonomy.txt",
		tax_good = config["rdir"] + "/tax_combined/full_taxonomy_good.txt",
		tax_prot = config["rdir"] + "/{library_name}/assembly_protein_taxonomy.txt"
	output:
		tax_prot_added = config["rdir"] + "/tax_combined/{library_name}_protein_taxonomy.txt"
	params:
		add = lambda wildcards: config["custom_ncbi_pre_derep"][wildcards.library_name],
		outdir = config["rdir"] + "/{library_name}/proteins"
	shell:
		"""
		if [[ "{params.add}" != "" ]]
		then
		  awk '$4 != "NA"' {params.add} | grep -F -f <(cut -f1 {input.tax}) > "{params.outdir}/tmp"
		  cut -f4 "{params.outdir}/tmp" | while read line
		  do
		    ln -sf "$line" {params.outdir}
		  done
		  cat {input.tax_prot} > {output.tax_prot_added}
		  awk '$4 != "NA"' "{params.outdir}/tmp" | cut -f1 | grep -F -f - {input.tax_good} >> {output.tax_prot_added}
		else
		  cp {input.tax_prot} {output.tax_prot_added}
		fi
		"""

rule collect_ncbi_proteins:
	input:
		tax_prot_added = config["rdir"] + "/tax_combined/{library_name}_protein_taxonomy.txt"
	output:
		linked = config["rdir"] + "/{library_name}/proteins/linked"
	params:
		pdir = config["rdir"] + "/{library_name}/proteins",
		outdir = config["rdir"] + "/proteins_all/"
	shell:
		"""
		mkdir -p {params.outdir}
		find {params.pdir} -type f -name '*.gz' | while read line
		do
		  ln -sf "$line" {params.outdir}
		done
		touch {output.linked}
		"""

if config["custom_euk_prot"]:
	rule add_custom_euk_prot:
		input:
			tax_good = config["rdir"] + "/tax_combined/full_taxonomy_good.txt"
		output:
			prot_added = config["rdir"] + "/tax_combined/euk_custom_protein_taxonomy.txt"
		params:
			add = config["custom_euk_prot"],
			outdir = config["rdir"] + "/proteins_all/"
		shell:
			"""
			mkdir -p {params.outdir}
			cut -f3 {params.add} | while read line
			do
			  ln -sf "$line" {params.outdir}
			done
			cut -f1 {params.add} | grep -F -f - {input.tax_good} > {output.prot_added}
			"""

