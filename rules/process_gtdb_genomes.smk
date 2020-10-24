rule get_gtdb_metadata:
	output:
		ar_tax = config["rdir"] + "/gtdb/metadata/ar_tax.tsv",
		bac_tax = config["rdir"] + "/gtdb/metadata/bac_tax.tsv",
		genomes_refseq = config["rdir"] + "/gtdb/metadata/genomes_refseq.tsv",
		genomes_genbank = config["rdir"] + "/gtdb/metadata/genomes_genbank.tsv",
		ar_meta = config["rdir"] + "/gtdb/metadata/ar_meta.tsv",
		bac_meta = config["rdir"] + "/gtdb/metadata/bac_meta.tsv"
	params:
		outdir = config["rdir"] + "/gtdb/metadata",
		gtdb_link = config["gtdb_link"]
	log:
                config["rdir"] + "/logs/get_gtdb_metadata.log"
	shell:
		"""
		wget -O {output.ar_tax} "{params.gtdb_link}/ar122_taxonomy.tsv" &>> {log}
		wget -O {output.bac_tax} "{params.gtdb_link}/bac120_taxonomy.tsv" &>> {log}
		wget -O "{params.outdir}/tmp_ar.tar.gz" "{params.gtdb_link}/ar122_metadata.tar.gz" &>> {log}
		wget -O "{params.outdir}/tmp_bac.tar.gz" "{params.gtdb_link}/bac120_metadata.tar.gz" &>> {log}
		tar -xzf "{params.outdir}/tmp_ar.tar.gz"
		tar -xzf "{params.outdir}/tmp_bac.tar.gz"
		rm "{params.outdir}/tmp_ar.tar.gz" "{params.outdir}/tmp_bac.tar.gz"
		mv "{params.outdir}/ar122_metadata*" {output.ar_meta}
		mv "{params.outdir}/bac120_metadata*" {output.bac_meta}
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
		genomes_refseq = config["rdir"] + "/gtdb/metadata/genomes_refseq.tsv",
		genomes_genbank = config["rdir"] + "/gtdb/metadata/genomes_genbank.tsv",
		ar_meta = config["rdir"] + "/gtdb/metadata/ar_meta.tsv",
                bac_meta = config["rdir"] + "/gtdb/metadata/bac_meta.tsv"
	output:
		gtdb_links = config["rdir"] + "/gtdb/metadata/gtdb_download_info.txt",
		gtdb_meta = config["rdir"] + "/gtdb/metadata/gtdb_metadata.tsv"
	params:
		script = config["wdir"] + "/scripts/prepare_files_gtdb.R"
	conda:
		config["wdir"] + "/envs/r.yaml"
	log:
                config["rdir"] + "/logs/parse_gtdb_metadata.log"
	shell:
		"""
		{params.script} -a {input.ar_tax} -b {input.bac_tax} -A {input.ar_meta} -B {input.bac_meta} -r {input.genomes_refseq} -g {input.genomes_genbank} -o {output.gtdb_links} -m {output.gtdb_meta} &>> {log}
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

localrules: derep_gtdb

rule derep_gtdb:
	input:
		download_complete_ncbi = config["rdir"] + "/gtdb/genomes/done",
		tax_ncbi = config["rdir"] + "/gtdb/metadata/gtdb_download_info.txt"
	output:
		taxonomy = config["rdir"] + "/gtdb/metadata/gtdb_taxonomy.txt",
		# derep_db = config["rdir"] + "/gtdb/derep_genomes/gtdb_derep_db",
		derep_meta = config["rdir"] + "/gtdb/metadata/gtdb_derep_taxonomy_meta.txt"
	params:
		indir = config["rdir"] + "/gtdb/genomes",
		outdir = config["rdir"] + "/gtdb/derep_genomes",
		z_threshold = config["z_threshold_gtdb"],
		derep_slurm = config["wdir"] + "/config/cluster_derep.yaml",
		derep_chunks = config["gtdb_derep_chunks"]
	threads: config["derep_threads"]
	conda:
		config["wdir"] + "/envs/derep.yaml"
	log:
                # config["rdir"] + "/logs/derep_gtdb.log"
		config["rdir"] + "/logs/derep_new_gtdb.log"
	shell:
		"""
		# parse taxonomy file
		cut -f6,7 {input.tax_ncbi} > {output.taxonomy}
		cd {params.outdir}
		derepG --threads {threads} --in-dir {params.indir} --taxa {output.taxonomy} --tmp ./ --slurm-config {params.derep_slurm} --db gtdb_derep_db --threshold {params.z_threshold} --chunk-size {params.derep_chunks} --debug --slurm-arr-size 10000 &>> {log}
		mv *derep-genomes_results.tsv {output.derep_meta}
		# do not delete redundant genomes until DB workflow is finished, work with soft links for remaining steps
		"""

rule collect_gtdb_genomes:
 	input:
 		derep_meta = config["rdir"] + "/gtdb/metadata/gtdb_derep_taxonomy_meta.txt"
 	output:
 		tax = config["rdir"] + "/tax_combined/gtdb_derep_taxonomy.txt"
	params:
		outdir = config["rdir"] + "/derep_combined/"
	shell:
		"""
		mkdir -p {params.outdir}
		cut -f4 {input.derep_meta} | sed '1d' | while read line
		do
		  ln -s "$line" {params.outdir}
		done
		# find {params.indir} -type f -name '*.gz' | xargs -n 1 mv -t {params.outdir}
		awk -v FS="\\t" -v OFS="\\t" '{{print $2,$1}}' {input.derep_meta} | sed '1d' > {output.tax}
		"""

