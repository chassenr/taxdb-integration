rule download_genomes_checkv:
	output: 
		fna = config["rdir"] + "/checkv/checkv_full.fna",
		clusters = config["rdir"] + "/checkv/checkv_clusters.tsv",
		meta_genbank = config["rdir"] + "/checkv/checkv_genbank.tsv",
		meta_circular = config["rdir"] + "/checkv/checkv_circular.tsv"
	params:
		library_dir = config["rdir"],
		checkv_link = config["checkv_link"]
	log:
		config["rdir"] + "/logs/download_genomes_checkv.log"
	shell:
		"""
		wget -O {output.fna} "{params.checkv_link}/checkv_full.fna" &>> {log}
		wget -O {output.clusters} "{params.checkv_link}/checkv_clusters.tsv" &>> {log}
		wget -O {output.meta_genbank} "{params.checkv_link}/checkv_genbank.tsv"
		wget -O {output.meta_circular} "{params.checkv_link}/checkv_circular.tsv"
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

localrules: derep_checkv

rule derep_checkv:
	input:
		split_done = config["rdir"] + "/checkv/genomes/done",
		checkv_taxonomy = config["rdir"] + "/checkv/checkv_taxonomy.txt"
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
		derepG --threads {threads} --in-dir {params.indir} --taxa {input.checkv_taxonomy} --tmp ./ --slurm-config {params.derep_slurm} --db {params.derep_db} --threshold {params.z_threshold} --mash-threshold {params.m_threshold} --chunk-size {params.derep_chunks} --debug --slurm-arr-size 10000 &>> {log}
		mv *derep-genomes_results.tsv {output.derep_meta}
		# do not delete redundant genomes until DB workflow is finished, work with soft links for remaining steps
		"""

rule collect_checkv_genomes:
	input:
		derep_meta = config["rdir"] + "/checkv/checkv_derep_taxonomy_meta.txt"
	output:
		tax = config["rdir"] + "/tax_combined/checkv_derep_taxonomy.txt"
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

