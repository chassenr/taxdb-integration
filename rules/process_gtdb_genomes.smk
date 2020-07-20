rule get_gtdb_metadata:
	output:
		ar_tax = config["rdir"] + "/gtdb/metadata/ar_tax.tsv",
		bac_tax = config["rdir"] + "/gtdb/metadata/bac_tax.tsv",
		ar_meta = config["rdir"] + "/gtdb/metadata/ar_meta.tsv",
		bac_meta = config["rdir"] + "/gtdb/metadata/bac_meta.tsv",
		sp_cluster = config["rdir"] + "/gtdb/metadata/sp_cluster.tsv",
		genomes_refseq = config["rdir"] + "/gtdb/metadata/genomes_refseq.tsv",
		genomes_genbank = config["rdir"] + "/gtdb/metadata/genomes_genbank.tsv"
	params:
		outdir = config["rdir"] + "/gtdb/metadata",
		gtdb_link = config["gtdb_link"]
	log:
                config["rdir"] + "/logs/get_gtdb_metadata.log"
	shell:
		"""
		wget -O {output.ar_tax} "{params.gtdb_link}/ar122_taxonomy.tsv" &>> {log}
		wget -O {output.bac_tax} "{params.gtdb_link}/bac120_taxonomy.tsv" &>> {log}
		wget -O {output.ar_meta} "{params.gtdb_link}/ar122_metadata.tsv" &>> {log}
		wget -O {output.bac_meta} "{params.gtdb_link}/bac120_metadata.tsv" &>> {log}
		wget -O {output.sp_cluster} "{params.gtdb_link}/sp_clusters.tsv" &>> {log}
		wget -O "{params.outdir}/assembly_summary_refseq.txt" http://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt &>> {log}
		wget -O "{params.outdir}/assembly_summary_genbank.txt" http://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt &>> {log}
		# parse NCBI data
		awk -v FS="\\t" -v OFS="\\t" '$0 !~ /^#/{{print $1,$18,$20}}' "{params.outdir}/assembly_summary_refseq.txt" > {output.genomes_refseq}
		awk -v FS="\\t" -v OFS="\\t" '$0 !~ /^#/{{print $1,$18,$20}}' "{params.outdir}/assembly_summary_genbank.txt" > {output.genomes_genbank}
		"""

rule parse_gtdb_metadata:
	input:
		ar_tax = config["rdir"] + "/gtdb/metadata/ar_tax.tsv",
		bac_tax = config["rdir"] + "/gtdb/metadata/bac_tax.tsv",
		ar_meta = config["rdir"] + "/gtdb/metadata/ar_meta.tsv",
		bac_meta = config["rdir"] + "/gtdb/metadata/bac_meta.tsv",
		sp_cluster = config["rdir"] + "/gtdb/metadata/sp_cluster.tsv",
		genomes_refseq = config["rdir"] + "/gtdb/metadata/genomes_refseq.tsv",
		genomes_genbank = config["rdir"] + "/gtdb/metadata/genomes_genbank.tsv"
	output:
		gtdb_links = config["rdir"] + "/gtdb/metadata/gtdb_download_info.txt",
		uba = config["rdir"] + "/gtdb/metadata/gtdb_uba.txt"
	params:
		outdir = config["rdir"] + "/gtdb/metadata",
		script = config["wdir"] + "/scripts/prepare_files_gtdb.R"
	conda:
		config["wdir"] + "/envs/r.yaml"
	log:
                config["rdir"] + "/logs/parse_gtdb_metadata.log"
	shell:
		"""
		{params.script} -a {input.ar_tax} -b {input.bac_tax} -r {input.genomes_refseq} -g {input.genomes_genbank} -o {params.outdir} &>> {log}
		"""	
	
rule download_gtdb_ncbi:
	input:
		gtdb_links = config["rdir"] + "/gtdb/metadata/gtdb_download_info.txt"
	output:
		config["rdir"] + "/gtdb/genomes/done"
	params:
		outdir = config["rdir"] + "/gtdb/genomes"
	threads: config["download_threads"]
	conda:
		config["wdir"] + "/envs/download.yaml"
	log:
                config["rdir"] + "/logs/download_gtdb_ncbi.log"
	shell:
		"""
		cut -f2,4 {input.gtdb_links} | awk -v FS="\\t" -v OFS="\\t" '{{print $1"\\n out="$2}}' > "{params.outdir}/links"
		aria2c -i "{params.outdir}/links" -c -l "{params.outdir}/links.log" --dir {params.outdir} --max-tries=20 --retry-wait=5 --max-connection-per-server=1 --max-concurrent-downloads={threads} &>> {log}
		# We need to verify all files are there
		cut -f4 {input.gtdb_links} | sort > "{params.outdir}/tmp1" 
		find {params.outdir} -type f -name '*.gz' | xargs -n 1 basename | sort > "{params.outdir}/tmp2"
		if diff "{params.outdir}/tmp1" "{params.outdir}/tmp2" 
		then
		  touch "{params.outdir}/done"
		fi
		rm "{params.outdir}/links" "{params.outdir}/links.log" "{params.outdir}/tmp1" "{params.outdir}/tmp2"
		# In case we also want to get the md5 files
		# cut -f5,1 {input.gtdb_links} | awk '{{print $2"\\n out="$1".md5"}}' > "{params.outdir}/links_md5"
		"""

rule download_gtdb_uba:
	input:
		uba = config["rdir"] + "/gtdb/metadata/gtdb_uba.txt"
	output:
		config["rdir"] + "/gtdb/genomes_uba/done"
	params:
		outdir = config["rdir"] + "/gtdb",
		gtdb_link = config["gtdb_link"]
	log:
                config["rdir"] + "/logs/download_gtdb_uba.log"
	shell:
		"""
		wget -O "{params.outdir}/genomes_uba/gtdb_uba_mags_arc.tar.gz" "{params.gtdb_link}/gtdb_uba_mags_arc.tar.gz" &>> {log}
		wget -O "{params.outdir}/genomes_uba/gtdb_uba_mags.tar.gz" "{params.gtdb_link}/gtdb_uba_mags.tar.gz" &>> {log}
		tar -xzvf "{params.outdir}/genomes_uba/gtdb_uba_mags_arc.tar.gz" -C "{params.outdir}/genomes_uba" &>> {log}
		tar -xzvf "{params.outdir}/genomes_uba/gtdb_uba_mags.tar.gz" -C "{params.outdir}/genomes_uba" &>> {log}
		rm "{params.outdir}/genomes_uba/gtdb_uba_mags_arc.tar.gz" "{params.outdir}/genomes_uba/gtdb_uba_mags.tar.gz"
		# only keep those UBA listed in the metadata (some UBA genomes are available already on NCBI, and this will avoid duplication)
		cut -f1 {input.uba} | sed 's/$/\\./' | grep -v -F -f - <(find "{params.outdir}/genomes_uba" -name 'UBA*') | xargs -n 1 rm
		gzip {params.outdir}/genomes_uba/UBA*
		mv {params.outdir}/genomes_uba/*.gz {params.outdir}/genomes/ 
		touch "{params.outdir}/genomes_uba/done"
		"""

rule derep_gtdb:
	input:
		download_complete_ncbi = config["rdir"] + "/gtdb/genomes/done",
		download_complete_uba = config["rdir"] + "/gtdb/genomes_uba/done",
		tax_ncbi = config["rdir"] + "/gtdb/metadata/gtdb_download_info.txt",
		tax_uba = config["rdir"] + "/gtdb/metadata/gtdb_uba.txt"
	output:
		taxonomy = config["rdir"] + "/gtdb/metadata/gtdb_taxonomy.txt",
		derep_taxonomy = config["rdir"] + "/gtdb/metadata/gtdb_derep_taxonomy.txt"
	params:
		indir = config["rdir"] + "/gtdb/genomes",
		outdir = config["rdir"] + "/gtdb/derep_genomes",
		derep_script = config["derep_script"],
		derep_threshold = config["derep_threshold_gtdb"]
	threads: config["derep_threads"]
	conda:
		config["wdir"] + "/envs/derep.yaml"
	log:
                config["rdir"] + "/logs/derep_gtdb.log"
	shell:
		"""
		# parse taxonomy file
		cp {input.tax_uba} {output.taxonomy}
		cut -f6,7 {input.tax_ncbi} >> {output.taxonomy}
		{params.derep_script} --threads {threads} --threshold {params.derep_threshold} {params.indir} {params.outdir} {output.taxonomy} &>> {log}
		# only select dereplicated genomes from taxonomy table for further processing
		find {params.outdir} -type f -name '*.gz' | xargs -n1 basename | sed 's/\\.f.*//' | cut -d'_' -f1,2 | grep -w -F -f - {output.taxonomy} > {output.derep_taxonomy}
		# delete non-dereplicated genomes
		find {params.indir} -type f -name '*.gz' | xargs -n 1 -P {threads} rm
		"""

rule collect_gtdb_genomes:
 	input:
 		derep_taxonomy = config["rdir"] + "/gtdb/metadata/gtdb_derep_taxonomy.txt"
 	output:
 		tax = config["rdir"] + "/tax_combined/gtdb_derep_taxonomy.txt"
	params:
		indir = config["rdir"] + "/gtdb/derep_genomes",
		outdir = config["rdir"] + "/derep_combined/"
	shell:
		"""
		mkdir -p {params.outdir}
		find {params.indir} -type f -name '*.gz' | xargs -n 1 mv -t {params.outdir}
		cp {input.derep_taxonomy} {output.tax}
		"""

