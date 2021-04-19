rule download_proteins_ncbi:
	input:
		url = config["rdir"] + "/{library_name}/assembly_url_genomic.txt",
		tax = config["rdir"] + "/tax_combined/{library_name}_derep_taxonomy.txt",
		tax_good = config["rdir"] + "/tax_combined/full_taxonomy_good.txt"
	output:
		download_complete = config["rdir"] + "/{library_name}/proteins/done",
		tax_prot = config["rdir"] + "/tax_combined/{library_name}_protein_taxonomy.txt"
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
		  touch "{params.outdir}/done"
		  cut -d'_' -f1,2 "{params.outdir}/tmp1" | grep -F -f - {input.tax_combined} > {output.tax_prot}
		fi
		rm "{params.outdir}/links.list" "{params.outdir}/links.log" "{params.outdir}/tmp1" "{params.outdir}/tmp2"
		"""

rule collect_ncbi_proteins:
	input:
		tax_prot = config["rdir"] + "/tax_combined/{library_name}_protein_taxonomy.txt"
	output:
		linked = config["rdir"] + "/{library_name}/proteins/linked"
	params:
		pdir = config["rdir"] + "/{library_name}/proteins",
		outdir = config["rdir"] + "/protein_combined/"
	shell:
		"""
		mkdir -p {params.outdir}
		find {params.pdir} -type f -name '*.gz' | while read line
		do
		  ln -sf "$line" {params.outdir}
		done
		touch {output.linked}

if config["custom_euk_prot"]:
        rule add_custom_euk_prot:
		input:
			tax_good = config["rdir"] + "/tax_combined/full_taxonomy_good.txt"
                output:
                        prot_added = config["rdir"] + "/tax_combined/euk_custom_protein_taxonomy.txt"
                params:
                        add = config["custom_euk_prot"],
			indir = config["custom_euk_prot_dir"],
                        outdir = config["rdir"] + "/protein_combined/"
                shell:
                        """
			ls -1 {params.indir} | grep -F -f <(cut -f1 {params.add}) > {params.indir}/tmp
			if [[ "$(wc -l < {params.indir}/tmp)" -eq "$(wc -l < {params.add})" ]]
			then
			  cat {params.indir}/tmp | while read line
			  do
			    ln -sf "$line" {params.outdir}
			  done
			cut -f1 {params.add} | grep -F -f - {input.tax_good} > {output.prot_added}
			fi
			rm {params.dir}/tmp
			"""

rule cluster_proteins:
	input:
		tax_prot = config["rdir"] + "/tax_combined/{library_name}_protein_taxonomy.txt",
		linked = config["rdir"] + "/{library_name}/proteins/linked",
		
	output:
		

