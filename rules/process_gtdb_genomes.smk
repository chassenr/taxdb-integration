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

# to add cusom assemblies, manually include genome files in the respective directories and provide the taxonomy file name in the config file for custom_ncbi
rule add_custom_gtdb:
	input:
		download_info = config["rdir"] + "/gtdb/metadata/gtdb_download_info.txt"
	output:
		tax_gtdb = config["rdir"] + "/gtdb/metadata/gtdb_taxonomy.txt",
		tax_added = config["rdir"] + "/gtdb/metadata/gtdb_taxonomy_added.txt"
	params:
		add = config["custom_gtdb"],
		dir = config["rdir"] + "/gtdb/genomes"
	shell:
		"""
		# parse taxonomy file
		cut -f6,7 {input.download_info} > {output.tax_gtdb}
		# add additional genomes if provided
		if [[ "{params.add}" != "" ]]
		then
		  ls -1 {params.dir} | grep -F -f <(cut -f1 {params.add}) > {params.dir}/tmp
		  if [[ "$(wc -l < {params.dir}/tmp)" -eq "$(wc -l < {params.add})" ]]
		  then
		    cat {output.tax_gtdb} {params.add} > {output.tax_added}
		  fi
		  rm {params.dir}/tmp
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
		derep_meta = config["rdir"] + "/gtdb/metadata/gtdb_derep_taxonomy_meta.txt"
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
                # config["rdir"] + "/logs/derep_gtdb.log"
		config["rdir"] + "/logs/derep_new_gtdb.log"
	shell:
		"""
		mkdir -p {params.outdir}
		cd {params.outdir}
		derepG --threads {threads} --in-dir {params.indir} --taxa {input.tax_added} --tmp ./ --slurm-config {params.derep_slurm} --db {params.derep_db} --threshold {params.z_threshold} --mash-threshold {params.m_threshold} --chunk-size {params.derep_chunks} --debug --slurm-arr-size 10000 &>> {log}
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
		awk -v FS="\\t" -v OFS="\\t" '{{print $2,$1}}' {input.derep_meta} | sed '1d' > {output.tax}
		"""

rule masking_gtdb:
	input:
		file_list = config["rdir"] + "/kraken2_genomes/file_names_derep_genomes.txt",
		gtdb = config["rdir"] + "/tax_combined/gtdb_derep_taxonomy.txt",
		nodes = config["rdir"] + "/kraken2_db/taxonomy/nodes.dmp",
		names = config["rdir"] + "/kraken2_db/taxonomy/names.dmp"
	output:
		fasta = config["rdir"] + "/kraken2_db/library/gtdb/library.fna"
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	threads: config["masking_threads"]
	shell:
		"""
		cut -f1 {input.gtdb} | grep -F -f - {input.file_list} | parallel -j{threads} 'dustmasker -in {{}} -outfmt fasta' | sed -e '/^>/!s/[a-z]/x/g' >> {output.fasta}
		"""

rule prelim_map_gtdb:
	input:  
		fasta = config["rdir"] + "/kraken2_db/library/gtdb/library.fna"
	output:
		map = config["rdir"] + "/kraken2_db/library/gtdb/prelim_map.txt"
	params: 
		libdir = config["rdir"] + "/kraken2_db/library/gtdb"
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	shell:
		"""
		LC_ALL=C grep '^>' {input.fasta} | sed 's/^>//' > {params.libdir}/tmp.accnos
		NSEQ=$(wc -l {params.libdir}/tmp.accnos | cut -d' ' -f1)
		printf 'TAXID\\n%.0s' $(seq 1 $NSEQ) | paste - {params.libdir}/tmp.accnos | paste - <(cut -d'|' -f3 {params.libdir}/tmp.accnos) > {output.map}
		rm {params.libdir}/tmp.accnos
		"""

if config["kingdoms"]:
	rule separate_contam_ncbi:
		input:
			id_contam = config["cdir"] + "/decontamination/contam_id.accnos",
			fasta = config["rdir"] + "/kraken2_db/library/gtdb/library.fna",
			map = config["rdir"] + "/kraken2_db/library/gtdb/prelim_map.txt"
		output:
			fasta_contam = config["cdir"] + "/decontamination/gtdb_library_contam.fna"
		conda:
			config["wdir"] + "/envs/bbmap.yaml"
		log:
			config["rdir"] + "/logs/gtdb_contam_filter.log"
		shell:
			"""
			filterbyname.sh in={input.fasta} out={output.fasta_contam} names={input.id_contam} include=t
			"""

	rule remove_contam_ncbi:
		input:
			contam = config["rdir"] + "/decontamination/highres_db_conterm_prediction",
			fasta_contam = config["cdir"] + "/decontamination/gtdb_library_contam.fna"
		output:
			cleaned_fasta = config["cdir"] + "/decontamination/gtdb_cleaned.fna",
			cleaned_map = config["cdir"] + "/decontamination/gtdb_cleaned_map.txt"
		params:
			script = config["wdir"] + "/scripts/remove_contamination.R",
			contamdir = config["cdir"] + "/decontamination"
		conda:
			config["wdir"] + "/envs/r.yaml"
		log:
			config["rdir"] + "/logs/gtdb_contam_remove.log"
		shell:
			"""
			{params.script} -i {input.fasta_contam} -c {input.contam} -o {output.cleaned_fasta}
			C_ALL=C grep '^>' {output.cleaned_fasta} | sed 's/^>//' > "{params.contam_dir}/tmp.accnos"
			NSEQ=$(wc -l "{params.contam_dir}/tmp.accnos" | cut -d' ' -f1)
			printf 'TAXID\n%.0s' $(seq 1 $NSEQ) | paste - "{params.contam_dir}/tmp.accnos" | paste - <(cut -d'|' -f3 "{params.contam_dir}/tmp.accnos") > {output.cleaned_map}
			rm "{params.contam_dir}/tmp.accnos"
			"""

	rule concat_cleaned_ncbi:
		input:
			id_contam = config["cdir"] + "/decontamination/contam_id.accnos",
			fasta = config["rdir"] + "/kraken2_db/library/gtdb/library.fna",
			map = config["rdir"] + "/kraken2_db/library/gtdb/prelim_map.txt",
			cleaned_fasta = config["cdir"] + "/decontamination/gtdb_cleaned.fna",
			cleaned_map = config["cdir"] + "/decontamination/gtdb_cleaned_map.txt"
		output:
			done = config["cdir"] + "/decontamination/gtdb_cleaning.done
		params:
			tmpdir = config["rdir"] + "/kraken2_db/tmp/gtdb"
		conda:
			config["wdir"] + "/envs/bbmap.yaml"
		shell:
			"""
			mv {input.fasta} "{params.tmpdir}/library.fna"
			filterbyname.sh in="{params.tmpdir}/library.fna" out={input.fasta} names={input.id_contam} include=f
			mv {input.map} "{params.tmpdir}/prelim_map.txt"
			grep -v -F -f {input.id_contam} "{params.tmpdir}/prelim_map.txt" > {input.map}
			cat {input.cleaned_fasta} >> {input.fasta}
			cat {input.cleaned_map} >> {input.map}
			rm "{params.tmpdir}/prelim_map.txt" "{params.tmpdir}/library.fna"
			touch {output}
			"""

