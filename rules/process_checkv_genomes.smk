rule download_genomes_checkv:
	output: 
		fna = config["rdir"] + "/checkv/checkv_full.fna",
		clusters = config["rdir"] + "/checkv/checkv_clusters.tsv",
		meta_genbank = config["rdir"] + "/checkv/checkv_genbank.tsv",
		meta_circular = config["rdir"] + "/checkv/checkv_circular.tsv",
		fna_reps = config["cdir"] + "/checkv/checkv_reps.fna",
		faa = config["rdir"] + "/checkv/checkv_full.faa",
		faa_reps = config["cdir"] + "/checkv/checkv_reps.faa"
	params:
		checkv_link = config["checkv_link"],
		outdir = config["rdir"] + "/checkv/"
	log:
		config["rdir"] + "/logs/download_genomes_checkv.log"
	shell:
		"""
		wget -P {params.outdir} {params.checkv_link} &>> {log}
		mkdir -p "{params.outdir}/checkv_tmp"
		tar -C "{params.outdir}/checkv_tmp" -xzf {params.outdir}/*.tar.gz
		mv {params.outdir}/checkv_tmp/checkv-db-*/genome_db/checkv_full.fna {output.fna}
		mv {params.outdir}/checkv_tmp/checkv-db-*/genome_db/checkv_clusters.tsv {output.clusters}
		mv {params.outdir}/checkv_tmp/checkv-db-*/genome_db/checkv_genbank.tsv {output.meta_genbank}
		mv {params.outdir}/checkv_tmp/checkv-db-*/genome_db/checkv_circular.tsv {output.meta_circular}
		mv {params.outdir}/checkv_tmp/checkv-db-*/genome_db/checkv_reps.fna {output.fna_reps}
		mv {params.outdir}/checkv_tmp/checkv-db-*/genome_db/checkv_full.faa {output.faa}
		mv {params.outdir}/checkv_tmp/checkv-db-*/genome_db/checkv_reps.faa {output.faa_reps}
		rm -rf "{params.outdir}/checkv_tmp"
		rm {params.outdir}/*.tar.gz
		"""

rule parse_taxa_checkv:
	input:
		meta_genbank = config["rdir"] + "/checkv/checkv_genbank.tsv",
                meta_circular = config["rdir"] + "/checkv/checkv_circular.tsv",
		clusters = config["rdir"] + "/checkv/checkv_clusters.tsv",
		nodes = config["rdir"] + "/ncbi_taxdump/nodes.dmp",
                names = config["rdir"] + "/ncbi_taxdump/names.dmp"
	output:
		checkv_taxonomy = config["rdir"] + "/checkv/checkv_taxonomy.txt",
		reps_metadata = config["rdir"] + "/checkv/checkv_reps_metadata.txt"
	params:
		script = config["wdir"] + "/scripts/parse_checkv_taxonomy.R",
		outdir = config["rdir"]
	conda:
		config["wdir"] + "/envs/r.yaml"
	log:
                config["rdir"] + "/logs/parse_taxa_checkv.log"
	shell:
		"""
		{params.script} -g {input.meta_genbank} -i {input.meta_circular} -c {input.clusters} -t "{params.outdir}/ncbi_taxdump" -s "{params.outdir}/ncbi_taxdump/accessionTaxa.sql" -o {output.checkv_taxonomy} -m {output.reps_metadata}
		"""

rule split_fasta:
	input:
		fna = config["rdir"] + "/checkv/checkv_full.fna"
	output:
		done = config["rdir"] + "/checkv/genomes/done"
	params:
		outdir = config["rdir"] + "/checkv/genomes/"
	shell:
		"""
		# thanks: https://gist.github.com/astatham/621901
		cd {params.outdir}
		cat {input} | awk '{{ if (substr($0, 1, 1)==">") {{filename=(substr($0,2) ".fa")}} print $0 > filename }}'
		touch done
		"""

rule add_custom_checkv:
	input:
		tax_checkv = config["rdir"] + "/checkv/checkv_taxonomy.txt"
	output:
		tax_added = config["rdir"] + "/checkv/checkv_taxonomy_added.txt"
	params:
		add = config["custom_checkv"],
		dir = config["rdir"] + "/checkv/genomes"
	shell:
		"""
		if [[ "{params.add}" != "" ]]
		then
		  ls -1 {params.dir} | grep -F -f <(cut -f1 {params.add}) > {params.dir}/tmp
		  if [[ "$(wc -l < {params.dir}/tmp)" -eq "$(wc -l < {params.add})" ]]
		  then
		    cat {input.tax_checkv} {params.add} > {output.tax_added}
		  fi
		  rm {params.dir}/tmp
		else
		  cp {input.tax_checkv} {output.tax_added}
		fi
		"""

localrules: derep_checkv

rule derep_checkv:
	input:
		split_done = config["rdir"] + "/checkv/genomes/done",
		checkv_taxonomy = config["rdir"] + "/checkv/checkv_taxonomy_added.txt"
	output:
		derep_meta = config["rdir"] + "/checkv/checkv_derep_taxonomy_meta.txt"
	params:
		indir = config["rdir"] + "/checkv/genomes",
		outdir = config["rdir"] + "/checkv/derep_genomes",
		z_threshold = config["z_threshold_checkv"],
		m_threshold = config["m_threshold_checkv"],
		derep_db = config["rdir"] + "/checkv/derep_genomes/checkv_derep_db",
		derep_slurm = config["wdir"] + "/config/cluster_derep.yaml",
		derep_chunks = config["checkv_derep_chunks"]
	threads: config["derep_threads"]
	conda:
		config["wdir"] + "/envs/derep.yaml"
	log:
		config["rdir"] + "/logs/derep_checkv.log"
	shell:
		"""
		mkdir -p {params.outdir}
		cd {params.outdir}
		derepG --threads {threads} --in-dir {params.indir} --taxa {input.checkv_taxonomy} --tmp ./ --db {params.derep_db} --threshold {params.z_threshold} --mash-threshold {params.m_threshold} --debug --slurm-config {params.derep_slurm} --chunk-size {params.derep_chunks} --slurm-arr-size 10000 &>> {log}
		mv *derep-genomes_results.tsv {output.derep_meta}
		# do not delete redundant genomes until DB workflow is finished, work with soft links for remaining steps
		"""

rule collect_checkv_genomes:
	input:
		derep_meta = config["rdir"] + "/checkv/checkv_derep_taxonomy_meta.txt",
		checkv_taxonomy = config["rdir"] + "/checkv/checkv_taxonomy_added.txt"
	output:
		tax = config["rdir"] + "/tax_combined/checkv_derep_taxonomy.txt"
	params:
		gendir = config["rdir"] + "/checkv/genomes",
		outdir = config["rdir"] + "/derep_combined/"
	shell:
		"""
		mkdir -p {params.outdir}
		cut -f4 {input.derep_meta} | sed '1d' | while read line
		do
		  ln -sf "$line" {params.outdir}
		done
		awk -v FS="\\t" -v OFS="\\t" '{{print $2,$1}}' {input.derep_meta} | sed '1d' > {output.tax}
		# assuming that the only reason that derepG jobs may have failed is because of too divergent assemblies,
		# include all assemblies for these taxonomic paths
		if [[ $(cut -f2 {input.checkv_taxonomy} | sort -t$'\\t' | uniq | wc -l) != $(cut -f1 {input.derep_meta} | sed '1d' | sort -t$'\\t' | uniq | wc -l) ]]
		then
		  cut -f2 {input.checkv_taxonomy} | sort -t$'\\t' | uniq | grep -v -F -f <(cut -f1 {input.derep_meta} | sed '1d' | sort -t$'\\t' | uniq) | grep -F -f - {input.checkv_taxonomy} > "{params.gendir}/../tmp"
		  cut -f1 "{params.gendir}/../tmp" | sed 's/$/\.fa/' | while read line
		  do
		    ln -sf "{params.gendir}/$line" {params.outdir}
		  done
		  cat "{params.gendir}/../tmp" >> {output.tax}
		  rm "{params.gendir}/../tmp"
		fi
		"""

rule masking_checkv:
	input:
		file_list = config["rdir"] + "/kraken2_genomes/file_names_derep_genomes.txt",
		checkv = config["rdir"] + "/tax_combined/checkv_derep_taxonomy.txt",
		nodes = config["rdir"] + "/kraken2_db_pro/taxonomy/nodes.dmp",
		names = config["rdir"] + "/kraken2_db_pro/taxonomy/names.dmp"
	output:
		fasta = config["rdir"] + "/kraken2_db_pro/library/checkv/library.fna"
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	threads: config["masking_threads"]
	shell:
		"""
		cut -f1 {input.checkv} | grep -F -f - {input.file_list} | parallel -j{threads} 'dustmasker -in {{}} -outfmt fasta' | sed -e '/^>/!s/[a-z]/x/g' >> {output.fasta}
		"""

rule prelim_map_checkv:
	input:  
		fasta = config["rdir"] + "/kraken2_db_pro/library/checkv/library.fna"
	output:
		map = config["rdir"] + "/kraken2_db_pro/library/checkv/prelim_map.txt"
	params: 
		libdir = config["rdir"] + "/kraken2_db_pro/library/checkv"
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	shell:
		"""
		LC_ALL=C grep '^>' {input.fasta} | sed 's/^>//' > {params.libdir}/tmp.accnos
		NSEQ=$(wc -l {params.libdir}/tmp.accnos | cut -d' ' -f1)
		printf 'TAXID\\n%.0s' $(seq 1 $NSEQ) | paste - {params.libdir}/tmp.accnos | paste - <(cut -d'|' -f3 {params.libdir}/tmp.accnos) > {output.map}
		rm {params.libdir}/tmp.accnos
		"""

