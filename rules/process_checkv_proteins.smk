rule split_checkv_proteins:
	input:
		faa = config["rdir"] + "/checkv/checkv_full.faa"
	output:
		done = config["rdir"] + "/checkv/proteins/done"
	params:
		outdir = config["rdir"] + "/checkv/proteins/"
	conda:
		config["wdir"] + "/envs/parallel.yaml"
	threads: config["parallel_threads"]
	shell:
		"""
		mkdir -p {params.outdir}
		cd {params.outdir}
		# thanks to: https://bioinformatics.stackexchange.com/questions/3021/how-to-split-multifasta-based-on-partial-fasta-header
		awk -F'_' '{{if(/^>/){{sub("^>", ""); name=$1"_"$2; print > name".faa"}}else{{print > name".faa"}}}}' {input.faa}
		find . -type f -name '*.faa' | parallel -j {threads} gzip {{}}
		touch {output.done}
		"""

rule collect_checkv_proteins:
	input:
		done = config["rdir"] + "/checkv/proteins/done",
		tax = config["rdir"] + "/tax_combined/checkv_derep_taxonomy.txt"
	output:
		tax_prot = config["rdir"] + "/checkv/checkv_protein_taxonomy.txt"
	params:
		pdir = config["rdir"] + "/checkv/proteins",
		outdir = config["rdir"] + "/proteins_all/"
	shell:
		"""
		mkdir -p {params.outdir}
		find {params.pdir} -type f -name '*.gz' | grep -F -f <(cut -f1 {input.tax}) > "{params.pdir}/tmp"
		cat "{params.pdir}/tmp" | while read line
		do
		  ln -sf "$line" {params.outdir}
		done
		cat "{params.pdir}/tmp" | xargs -n 1 basename | sed 's/\.faa\.gz//' | grep -F -f - {input.tax} > {output.tax_prot}
		rm "{params.pdir}/tmp"
		"""

rule custom_checkv_proteins_pre_derep:
	input:
		tax = config["rdir"] + "/tax_combined/checkv_derep_taxonomy.txt",
		tax_prot = config["rdir"] + "/checkv/checkv_protein_taxonomy.txt"
	output:
		tax_prot_added = config["rdir"] + "/tax_combined/checkv_protein_taxonomy.txt"
	params:
		add = config["custom_checkv_pre_derep"],
		outdir = config["rdir"] + "/proteins_all/"
	shell:
		"""
		if [[ "{params.add}" != "" ]]
		then
		  awk -v FS="\\t" -v OFS="\\t" '$4 != "NA"' {params.add} | grep -F -f <(cut -f1 {input.tax}) > "{params.outdir}/tmp" 
		  cut -f4 "{params.outdir}/tmp" | while read line
		  do
		    ln -sf "$line" {params.outdir}
		  done
		  cat {input.tax_prot} > {output.tax_prot_added}
		  cut -f1,2 "{params.outdir}/tmp" >> {output.tax_prot_added}
		else
		  cp {input.tax_prot} {output.tax_prot_added}
		fi
		"""

if config["custom_vir_prot"]:
	rule add_custom_vir_prot:
		output:
			prot_added = config["rdir"] + "/tax_combined/vir_custom_protein_taxonomy.txt"
		params:
			add = config["custom_vir_prot"],
			outdir = config["rdir"] + "/proteins_all/"
		shell:
			"""
			mkdir -p {params.outdir}
			cut -f3 {params.add} | while read line
			do
			  ln -sf "$line" {params.outdir}
			done
			cut -f1,2 {params.add} > {output.prot_added}
			"""

