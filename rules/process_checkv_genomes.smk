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
		outdir = config["wdir"]
	conda:
		config["wdir"] + "/envs/r.yaml"
	log:
                config["rdir"] + "/logs/parse_taxa_checkv.log"
	shell:
		"""
		{params.script} -g {input.meta_genbank} -i {input.meta_circular} -c {input.clusters} -t "{params.outdir}/ncbi_taxdump" -s "{params.outdir}/ncbi_taxdump/accessionTaxa.sql" -o {output.checkv_taxonomy} -m {output.reps_metadata}
		"""

localrules: derep_checkv

rule derep_checkv:
	input:
		checkv_taxonomy = config["rdir"] + "/checkv/checkv_taxonomy.txt",
		fna = config["rdir"] + "/checkv/checkv_full.fna"
	output:
		derep_db = config["rdir"] + "/checkv/derep_genomes/gtdb_derep_db",
		derep_meta = config["rdir"] + "/checkv/checkv_derep_taxonomy_meta.txt"
	params:
		dir = config["rdir"] + "/checkv",
		outdir = config["rdir"] + "/checkv/derep_genomes",
		z_threshold = config["z_threshold_checkv"],
		derep_slurm = config["wdir"] + "/config/cluster_derep.yaml",
		derep_chunks = config["checkv_derep_chunks"]
	threads: config["derep_threads"]
	conda:
		config["wdir"] + "/envs/derep.yaml"
	log:
		config["rdir"] + "/logs/derep_checkv.log"
	shell:
		"""
		cd {params.outdir}
		derepG --threads {threads} --in-dir {params.indir} --taxa "{params.dir}/assembly_taxonomy_select.txt" --tmp ./ --slurm-config {params.derep_slurm} --db {output.derep_db} --threshold {params.z_threshold} --chunk-size {params.derep_chunks} --copy --out-dir {params.outdir} &>> {log}
		mv *derep-genomes_results.tsv {output.derep_meta}
		# do not delete redundant genomes until DB workflow is finished, work with soft links for remaining steps
		"""

rule collect_checkv_genomes:
	input:
		derep_meta = config["rdir"] + "/checkv/checkv_derep_taxonomy_meta.txt"
	output:
		tax = config["rdir"] + "/tax_combined/checkv_derep_taxonomy.txt"
	shell:
		"""
		awk -v FS="\\t" -v OFS="\\t" '{{print $2,$1}}' {input.derep_meta} | sed '1d' > {output.tax}
		"""

