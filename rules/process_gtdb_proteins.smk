rule download_proteins_gtdb:
	input:
		gtdb_links = config["rdir"] + "/gtdb/metadata/gtdb_download_info.txt",
		tax_good = config["rdir"] + "/tax_combined/full_taxonomy_good.txt",
		tax = config["rdir"] + "/tax_combined/gtdb_derep_taxonomy.txt"
	output:
		tax_prot = config["rdir"] + "/gtdb/metadata/gtdb_protein_taxonomy.txt"
	params:
		outdir = config["rdir"] + "/gtdb/proteins"
	conda:
		config["wdir"] + "/envs/download.yaml"
	threads: config["download_threads"]
	log:
		config["rdir"] + "/logs/download_proteins_gtdb.log"
	shell:
		"""
		awk '$11 == "TRUE"' {input.gtdb_links} |  grep -F -f <(cut -f1 {input.tax} > "{params.outdir}/tmp"
		cut -f8,10 "{params.outdir}/tmp" | awk -v FS="\\t" -v OFS="\\t" '{{print $1"\\n out="$2}}' > "{params.outdir}/links"
		aria2c -i "{params.outdir}/links" -c -l "{params.outdir}/links.log" --dir {params.outdir} --max-tries=20 --retry-wait=5 --max-connection-per-server=1 --max-concurrent-downloads={threads} &>> {log}
		# We need to verify all files are there
		cut -f10 "{params.outdir}/tmp" | sort > "{params.outdir}/tmp1"
		find {params.outdir} -type f -name '*.gz' | xargs -n 1 basename | sort > "{params.outdir}/tmp2"
		if diff "{params.outdir}/tmp1" "{params.outdir}/tmp2"
		then
		  cut -f6,7 "{params.outdir}/tmp" > {output.tax_prot}
		fi
		rm "{params.outdir}/links" "{params.outdir}/links.log" "{params.outdir}/tmp1" "{params.outdir}/tmp2" "{params.outdir}/tmp"
		"""

rule custom_gtdb_proteins_pre_derep:
	input:
		tax = config["rdir"] + "/tax_combined/gtdb_derep_taxonomy.txt",
		tax_good = config["rdir"] + "/tax_combined/full_taxonomy_good.txt",
		tax_prot = config["rdir"] + "/gtdb/metadata/gtdb_protein_taxonomy.txt"
	output:
		tax_prot_added = config["rdir"] + "/gtdb/metadata/gtdb_protein_taxonomy_added.txt"
	params:
		add = config["custom_gtdb_pre_derep"],
		outdir = config["rdir"] + "/gtdb/proteins/"
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
		  cut -f1 "{params.pdir}/tmp" | grep -F -f - {input.tax_good} >> {output.tax_prot_added}
		else
		  cp {input.tax_prot} {output.tax_prot_added}
		fi
		"""

rule download_gtdb_faa_reps:
	output:
		gtdb_faa_done = config["rdir"] + "/gtdb/proteins/done"
	params:
		outdir = config["rdir"] + "/gtdb",
		gtdb_link = config["gtdb_link"]
	threads: config["download_onefile"]
	conda:
		config["wdir"] + "/envs/download.yaml"
	log:
		config["rdir"] + "/logs/download_faa_reps_gtdb.log"
	shell:
		"""
		aria2c -c -l "{params.outdir}/faa_reps_links.log" -d {params.outdir} --max-tries=20 --retry-wait=5 -x {threads} -j {threads} -s {threads} "{params.gtdb_link}/genomic_files_reps/gtdb_proteins_aa_reps.tar.gz" &>> {log}
		tar -C {params.outdir} -xzf "{params.outdir}/gtdb_proteins_aa_reps.tar.gz"
		mv {params.outdir}/gtdb_protein_aa_reps_* "{params.outdir}/reps_proteins/"
		rm "{params.outdir}/gtdb_protein_aa_reps.tar.gz"
		"""

rule gtdb_add_protein_reps:
	input:
		gtdb_faa_done = config["rdir"] + "/gtdb/proteins/done",
		tax_prot_added = config["rdir"] + "/gtdb/metadata/gtdb_protein_taxonomy_added.txt",
		gtdb_reps = config["rdir"] + "/gtdb/metadata/gtdb_reps_tax.txt",
		tax_good = config["rdir"] + "/tax_combined/full_taxonomy_good.txt"
	output:
		tax_prot = config["rdir"] + "/tax_combined/gtdb_protein_taxonomy.txt"
	params:
		metadir = config["rdir"] + "/gtdb/metadata",
		repdir = config["rdir"] + "/gtdb/reps_proteins",
		outdir = config["rdir"] + "/gtdb/proteins/"
	shell:
		"""
		cut -f2 {input.tax_prot_added} | sort -t$'\\t' | uniq | grep -v -F -f - {input.gtdb_reps} > "{params.metadir}/tmp_reps_prot.txt"
		cut -f1 "{params.metadir}/tmp_reps_prot.txt" | sed 's/$/_protein\.faa\.gz/' | while read line
		do
		  ln -sf "{params.repdir}/$line" {params.outdir}
		done
		cat {input.tax_prot_added} > {output.tax}
		cut -f1 "{params.metadir}/tmp_reps_prot.txt" | grep -F -f - {input.tax_good} >> {output.tax}
                rm "{params.metadir}/tmp_reps_prot.txt"
		"""

rule collect_gtdb_proteins:
	input:
		tax_prot = config["rdir"] + "/tax_combined/gtdb_protein_taxonomy.txt"
	output:
		linked = config["rdir"] + "/gtdb/proteins/linked"
	params:
		pdir = config["rdir"] + "/gtdb/proteins",
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

if config["custom_pro_prot"]:
	rule add_custom_pro_prot:
		input:
			tax_good = config["rdir"] + "/tax_combined/full_taxonomy_good.txt"
		output:
			prot_added = config["rdir"] + "/tax_combined/pro_custom_protein_taxonomy.txt"
		params:
			add = config["custom_pro_prot"],
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

