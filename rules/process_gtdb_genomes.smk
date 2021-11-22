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
		mkdir -p {params.outdir}
		wget -O {output.ar_tax} "{params.gtdb_link}/ar122_taxonomy.tsv" &>> {log}
		wget -O {output.bac_tax} "{params.gtdb_link}/bac120_taxonomy.tsv" &>> {log}
		wget -O "{params.outdir}/tmp_ar.tar.gz" "{params.gtdb_link}/ar122_metadata.tar.gz" &>> {log}
		wget -O "{params.outdir}/tmp_bac.tar.gz" "{params.gtdb_link}/bac120_metadata.tar.gz" &>> {log}
		tar -xzf "{params.outdir}/tmp_ar.tar.gz" -C {params.outdir}
		tar -xzf "{params.outdir}/tmp_bac.tar.gz" -C {params.outdir}
		rm "{params.outdir}/tmp_ar.tar.gz" "{params.outdir}/tmp_bac.tar.gz"
		mv {params.outdir}/ar122_metadata* {output.ar_meta}
		mv {params.outdir}/bac120_metadata* {output.bac_meta}
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
	threads: config["download_threads"]
	log:
                config["rdir"] + "/logs/parse_gtdb_metadata.log"
	shell:
		"""
		{params.script} -a {input.ar_tax} -b {input.bac_tax} -A {input.ar_meta} -B {input.bac_meta} -r {input.genomes_refseq} -g {input.genomes_genbank} -o {output.gtdb_links} -m {output.gtdb_meta} -c {threads} &>> {log}
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

# also download gtdb reps separately to augment any missing species in the NCBI download
rule download_gtdb_reps:
	input:
		ar_tax = config["rdir"] + "/gtdb/metadata/ar_tax.tsv",
		bac_tax = config["rdir"] + "/gtdb/metadata/bac_tax.tsv",
		ar_meta = config["rdir"] + "/gtdb/metadata/ar_meta.tsv",
		bac_meta = config["rdir"] + "/gtdb/metadata/bac_meta.tsv"
	output:
		gtdb_reps = config["rdir"] + "/gtdb/metadata/gtdb_reps_tax.txt"
	params:
		outdir = config["rdir"] + "/gtdb",
		gtdb_link = config["gtdb_link"]
	threads: config["download_onefile"]
	conda:
		config["wdir"] + "/envs/download.yaml"
	log:
		config["rdir"] + "/logs/download_coarse_gtdb.log"
	shell:
		"""
		aria2c -c -l "{params.outdir}/reps_links.log" -d {params.outdir} --max-tries=20 --retry-wait=5 -x {threads} -j {threads} -s {threads} "{params.gtdb_link}/genomic_files_reps/gtdb_genomes_reps.tar.gz" &>> {log}
		# wget -P {params.outdir} "{params.gtdb_link}/genomic_files_reps/gtdb_genomes_reps.tar.gz"
		tar -C {params.outdir} -xzf "{params.outdir}/gtdb_genomes_reps.tar.gz"
		mv {params.outdir}/gtdb_genomes_reps_* "{params.outdir}/reps_genomes/"
		rm "{params.outdir}/gtdb_genomes_reps.tar.gz"
		awk -v FS="\\t" -v OFS="\\t" '$16 == "t"' {input.ar_meta} | cut -f1 | grep -F -f - {input.ar_tax} | sed 's/^[RG][SB]_//' | sed 's/d__Archaea;/d__Archaea;l__Archaea;k__Archaea;/' > {output.gtdb_reps}
		awk -v FS="\\t" -v OFS="\\t" '$16 == "t"' {input.bac_meta} | cut -f1 | grep -F -f - {input.bac_tax} | sed 's/^[RG][SB]_//' | sed 's/d__Bacteria;/d__Bacteria;l__Bacteria;k__Bacteria;/' >> {output.gtdb_reps}
		"""

# to add cusom assemblies, manually include genome files in the respective directories and provide the taxonomy file name in the config file
rule add_custom_gtdb_pre_derep:
	input:
		download_info = config["rdir"] + "/gtdb/metadata/gtdb_download_info.txt"
	output:
		tax_gtdb = config["rdir"] + "/gtdb/metadata/gtdb_taxonomy.txt",
		tax_added = config["rdir"] + "/gtdb/metadata/gtdb_taxonomy_added.txt"
	params:
		add = config["custom_gtdb_pre_derep"],
		gendir = config["rdir"] + "/gtdb/genomes"
	shell:
		"""
		# parse taxonomy file
		cut -f6,7 {input.download_info} > {output.tax_gtdb}
		# add additional genomes if provided
		if [[ "{params.add}" != "" ]]
		then
		  cut -f3 {params.add} | while read line
		  do
		    ln -sf "$line" {params.gendir}
		  done
		  cat {output.tax_gtdb} > {output.tax_added}
		  cut -f1,2 {params.add} >> {output.tax_added}
		else
		  cp {output.tax_gtdb} {output.tax_added}
		fi
		"""

localrules: derep_gtdb

rule derep_gtdb:
	input:
		download_complete_ncbi = config["rdir"] + "/gtdb/genomes/done",
		tax_added = config["rdir"] + "/gtdb/metadata/gtdb_taxonomy_added.txt"
	output:
		derep_meta = config["rdir"] + "/gtdb/metadata/gtdb-derep-genomes_results.tsv"
	params:
		indir = config["rdir"] + "/gtdb/genomes",
		outdir = config["rdir"] + "/gtdb/derep_genomes",
		z_threshold = config["z_threshold_gtdb"],
		m_threshold = config["m_threshold_gtdb"],
		derep_db = config["rdir"] + "/gtdb/derep_genomes/gtdb_derep_db",
		derep_slurm = config["wdir"] + "/config/cluster_derep.yaml",
		derep_chunks = config["gtdb_derep_chunks"]
	threads: config["derep_threads"]
	conda:
		config["wdir"] + "/envs/derep.yaml"
	log:
                config["rdir"] + "/logs/derep_gtdb.log"
	shell:
		"""
		mkdir -p {params.outdir}
		cd {params.outdir}
		derepG --threads {threads} --in-dir {params.indir} --taxa {input.tax_added} --tmp ./ --db {params.derep_db} --prefix ../metadata/gtdb --threshold {params.z_threshold} --mash-threshold {params.m_threshold} --debug --slurm-config {params.derep_slurm} --chunk-size {params.derep_chunks} --slurm-arr-size 10000 &>> {log}
		# do not delete redundant genomes until DB workflow is finished, work with soft links for remaining steps
		"""

rule collect_gtdb_genomes:
 	input:
 		derep_meta = config["rdir"] + "/gtdb/metadata/gtdb-derep-genomes_results.tsv",
		gtdb_reps = config["rdir"] + "/gtdb/metadata/gtdb_reps_tax.txt"
 	output:
 		tax = config["rdir"] + "/tax_combined/gtdb_derep_taxonomy.txt"
	params:
		metadir = config["rdir"] + "/gtdb/metadata",
		repdir = config["rdir"] + "/gtdb/reps_genomes",
		outdir = config["rdir"] + "/derep_combined/"
	shell:
		"""
		mkdir -p {params.outdir}
		cut -f4 {input.derep_meta} | sed '1d' | while read line
		do
		  ln -sf "$line" {params.outdir}
		done
		awk -v FS="\\t" -v OFS="\\t" '{{print $2,$1}}' {input.derep_meta} | sed '1d' > {output.tax}
		# for the species not in the dereplicated set (because there genomes could not be obtained from genbank/refseq) take at least the GTDB reps
		# very annoyingly, GTDB has again changed their directory structure
		cut -f1 {input.derep_meta} | sed '1d' | sort -t$'\\t' | uniq | grep -v -F -f - {input.gtdb_reps} > "{params.metadir}/tmp_reps.txt"
		find {params.repdir} -type f -name '*.gz' > "{params.metadir}/tmp_files.txt"
		cut -f1 "{params.metadir}/tmp_reps.txt" | grep -F -f - "{params.metadir}/tmp_files.txt" | while read line
		do
		  ln -sf $line {params.outdir}
		done
		cat "{params.metadir}/tmp_reps.txt" >> {output.tax}
		rm "{params.metadir}/tmp_reps.txt" "{params.metadir}/tmp_files.txt"
		"""

if config["custom_gtdb_post_derep"]:
	rule add_custom_gtdb_post_derep:
		output:
			tax_added = config["rdir"] + "/tax_combined/pro_custom_post_derep_taxonomy.txt"
		params:
			add = config["custom_gtdb_post_derep"],
			outdir = config["rdir"] + "/derep_combined/"
		shell:
			"""
			mkdir -p {params.outdir}
			cut -f4 {params.add} | sed '1d' | while read line
			do
			  ln -sf "$line" {params.outdir}
			done
			awk -v FS="\\t" -v OFS="\\t" '{{print $2,$1}}' {params.add} | sed '1d' > {output.tax}
			"""

