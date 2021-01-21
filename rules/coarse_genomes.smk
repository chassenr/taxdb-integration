rule select_ncbi_coarse:
	input:
		tax = config["rdir"] + "/{library_name}/assembly_taxonomy.txt",
		meta = config["rdir"] + "/{library_name}/metadata/genome_metadata.txt"
	output:
		tax_coarse = config["cdir"] + "/{library_name}/assembly_taxonomy_coarse.txt"
	params:
		script = config["wdir"] + "/scripts/coarse_ncbi_selection.R",
		rank = lambda wildcards: config["rank_coarse"][wildcards.library_name],
		sublineage = lambda wildcards: config["coarse_sublineage"][wildcards.library_name],
		rank_sublineage = lambda wildcards: config["rank_sublineage"][wildcards.library_name],
		nmax = config["nmax_coarse"]
	conda:
		config["wdir"] + "/envs/r.yaml"
	log:
		config["rdir"] + "/logs/select_coarse_{library_name}.log"
	shell:
		"""
		{params.script} -t {input.tax} -m {input.meta} -r {params.rank} -s "{params.sublineage}" -R {params.rank_sublineage} -n {params.nmax} -o {output.tax_coarse} &>>{log}
		"""

rule download_ncbi_coarse:
	input:
		tax_coarse = config["cdir"] + "/{library_name}/assembly_taxonomy_coarse.txt",
		url = config["rdir"] + "/{library_name}/assembly_url_genomic.txt"
	output:
		coarse_ncbi_download = config["cdir"] + "/{library_name}/genomes/done"
	params:
		outdir = config["cdir"] + "/{library_name}/genomes"
	threads: config["download_threads"]
	conda:
		config["wdir"] + "/envs/download.yaml"
	log:
		config["rdir"] + "/logs/download_coarse_ncbi_{library_name}.log"
	shell:
		"""
		cut -f1 {input.tax_coarse} | sed 's/$/_/' | grep -F -f - {input.url} > "{params.outdir}/links"
		aria2c -i "{params.outdir}/links" -c -l "{params.outdir}/links.log" --dir {params.outdir} --max-tries=20 --retry-wait=5 --max-connection-per-server=1 --max-concurrent-downloads={threads} &>> {log}
		# We need to verify all files are there
		cat "{params.outdir}/links" | xargs -n 1 basename | sort > "{params.outdir}/tmp1"
		find {params.outdir} -type f -name '*.gz' | xargs -n 1 basename | sort > "{params.outdir}/tmp2"
		if diff "{params.outdir}/tmp1" "{params.outdir}/tmp2" 
		then
		  touch "{params.outdir}/done"
		fi
		rm "{params.outdir}/links.log" "{params.outdir}/tmp1" "{params.outdir}/tmp2"
		"""

rule download_gtdb_reps:
	input:
		ar_tax = config["rdir"] + "/gtdb/metadata/ar_tax.tsv",
		bac_tax = config["rdir"] + "/gtdb/metadata/bac_tax.tsv",
		ar_meta = config["rdir"] + "/gtdb/metadata/ar_meta.tsv",
		bac_meta = config["rdir"] + "/gtdb/metadata/bac_meta.tsv"
	output:
		gtdb_reps = config["cdir"] + "/gtdb/gtdb_reps_tax.txt"
	params:
		outdir = config["cdir"] + "/gtdb",
		gtdb_link = config["gtdb_link"]
	threads: config["download_threads"]
	conda:
		config["wdir"] + "/envs/download.yaml"
	log:
		config["rdir"] + "/logs/download_coarse_gtdb.log"
	shell:
		"""
		aria2c -c -l "{params.outdir}/links.log" -d {params.outdir} --max-tries=20 --retry-wait=5 -x {threads} -j {threads} -s {threads} "{params.gtdb_link}/genomic_files_reps/gtdb_genomes_reps.tar.gz" &>> {log}
		# wget -P {params.outdir} "{params.gtdb_link}/genomic_files_reps/gtdb_genomes_reps.tar.gz"
		tar -xzvf "{params.outdir}/gtdb_genomes_reps.tar.gz"
		mv "{params.outdir}/gtdb_genomes_reps_*" "{params.outdir}/genomes/"
		awk -v FS="\\t" -v OFS="\\t" '$16 == "t"' {input.ar_meta} | cut -f1 | grep -F -f - {input.ar_tax} | sed 's/^[RG][SB]_//' > {output.gtdb_reps}
		awk -v FS="\\t" -v OFS="\\t" '$16 == "t"' {input.bac_meta} | cut -f1 | grep -F -f - {input.bac_tax} | sed 's/^[RG][SB]_//' >> {output.gtdb_reps}
		"""

rule download_checkv_reps:
	input:
		reps_metadata = config["rdir"] + "/checkv/checkv_reps_metadata.txt"
	output:
		checkv_tax = config["cdir"] + "/checkv/checkv_reps_tax.txt",
		fna = config["cdir"] + "/checkv/checkv_reps.fna",
		done = config["cdir"] + "/checkv/genomes/done"
	params:
		checkv_link = config["checkv_link"],
		outdir = config["cdir"] + "/checkv/genomes/"
	threads: config["download_threads"]
	conda:
		config["wdir"] + "/envs/kraken2.yaml" # this env has parallel
	log:
		config["rdir"] + "/logs/download_coarse_checkv.log"
	shell:
		"""
		wget -O {output.fna} "{params.checkv_link}/checkv_reps.fna" &>> {log}
		cut -f1,8 {input.reps_metadata} | sed '1d' > {output.checkv_tax}
		mkdir -p {params.outdir}
		cd {params.outdir}
		cat {output.fna} | awk '{{ if (substr($0, 1, 1)==">") {{filename=(substr($0,2) ".fa")}} print $0 > filename }}'
		find . -type f -name '*.fa' | parallel -j {threads} gzip {{}}
		touch done
		"""

rule collect_coarse:
	input:
		checkv_tax = config["cdir"] + "/checkv/checkv_reps_tax.txt",
		gtdb_reps = config["cdir"] + "/gtdb/gtdb_reps_tax.txt",
		tax_coarse = expand(config["cdir"] + "/{library_name}/assembly_taxonomy_coarse.txt", library_name = LIBRARY_NAME)
	output:
		tax_all = config["cdir"] + "/tax_coarse_all.txt",
		sn_done = config["cdir"] + "/genomes_select/done"
	params:
		script = config["wdir"] + "/scripts/fix_ncbi_taxpath.R",
		cdir = config["cdir"],
		outdir = config["cdir"] + "/genomes_select"
	log:
		config["rdir"] + "/logs/collect_coarse_genomes.log"
	conda:
		config["wdir"] + "/envs/r.yaml"
	shell:
		"""
		# fix NCBI path
		cat {input.tax_coarse} {input.checkv_tax} > "{params.outdir}/tmp"
		{params.script} -i "{params.outdir}/tmp" -g {input.gtdb_reps} -o {output.tax_all} &>> {log}
		rm "{params.outdir}/tmp"
		# create softlinks for genomes
		find {params.cdir}/*/genomes/ -type f -name '*.gz' | while read line
		do
		  ln -s "$line" {params.outdir}
		done
		touch {output.sn_done}
		"""

