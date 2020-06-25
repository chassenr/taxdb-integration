rule get_genomes_ncbi:
	output: 
		config["rdir"] + "/" + config["library_name"] + "/assembly_url_genomic.txt",
		config["rdir"] + "/" + config["library_name"] + "/assembly_summary_combined.txt"
	params:
		library_dir = config["rdir"],
		library_name = config["library_name"],
		assembly_level = config["assembly_level"],
		script = config["wdir"] + "/scripts/get_genomes_ncbi.sh"
	conda:
                config["wdir"] + "/envs/download.yaml"
	shell:
		"""
		{params.script} "{params.library_dir}" "{params.library_name}" "{params.assembly_level}"
		"""

rule download_genomes_ncbi:
	input:
		url = config["rdir"] + "/" + config["library_name"] + "/assembly_url_genomic.txt"
	output:
		download_complete = config["rdir"] + "/" + config["library_name"] + "/genomes/done"
	params:
		outdir = config["rdir"] + "/" + config["library_name"] + "/genomes"
	threads: config["download_threads"]
	conda:
		config["wdir"] + "/envs/download.yaml"
	shell:
		"""
		aria2c -i {input.url} -c -l links.log --dir {params.outdir} --max-tries=20 --retry-wait=5 --max-connection-per-server=1 --max-concurrent-downloads={threads}
		# We need to verify all files are there
		cat {input.url} | xargs -n 1 basename | sort > tmp1 
		ls -1 {params.outdir}/*.gz | xargs -n 1 basename > tmp2
		if diff tmp1 tmp2 
		then
		  touch "{params.outdir}/done"
		fi
		rm links.log tmp1 tmp2
		"""

rule get_taxdump:
	output:
		nodes = config["rdir"] + "/ncbi_taxdump/nodes.dmp",
		names = config["rdir"] + "/ncbi_taxdump/names.dmp"
	params:
		outdir = config["rdir"]
	conda:
		config["wdir"] + "/envs/download.yaml"
	shell:
		"""
		aria2c --max-tries=20 --retry-wait=5 --dir {params.outdir} http://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
		mkdir -p "{params.outdir}/ncbi_taxdump/"
		tar -xzvf "{params.outdir}/new_taxdump.tar.gz" --directory "{params.outdir}/ncbi_taxdump/"
		rm "{params.outdir}/new_taxdump.tar.gz"
		"""

rule get_taxpath:
	input:
		nodes = config["rdir"] + "/ncbi_taxdump/nodes.dmp",
		names = config["rdir"] + "/ncbi_taxdump/names.dmp",
		genomes = config["rdir"] + "/" + config["library_name"] + "/assembly_summary_combined.txt"
	output:
		taxonomy = config["rdir"] + "/" + config["library_name"] + "/assembly_taxonomy.txt"
	params:
		outdir = config["rdir"],
		script = config["wdir"] + "/scripts/get_taxpath.R"
	conda:
		config["wdir"] + "/envs/r.yaml"
	shell:
		"""
		{params.script} -i {input.genomes} -t "{params.outdir}/ncbi_taxdump" -s "{params.outdir}/ncbi_taxdump/accessionTaxa.sql" -o {output.taxonomy}
		"""

rule derep_genomes:
	input:
		download_complete = config["rdir"] + "/" + config["library_name"] + "/genomes/done",
		taxonomy = config["rdir"] + "/" + config["library_name"] + "/assembly_taxonomy.txt"
	output:
		derep_taxonomy = config["rdir"] + "/" + config["library_name"] + "/derep_assembly_taxonomy.txt"
	params:
		indir = config["rdir"] + "/" + config["library_name"] + "/genomes",
		outdir = config["rdir"] + "/" + config["library_name"] + "/derep_genomes",
		derep_script = config["derep_script"],
		derep_threshold = config["derep_threshold"]
	threads: config["derep_threads"]
	conda:
		config["wdir"] + "/envs/derep.yaml"
	shell:
		"""
		{params.derep_script} --threads {threads} --threshold {params.derep_threshold} {params.indir} {params.outdir} {input.taxonomy}
		# only select dereplicated genomes from taxonomy table for further processing
		ls -1 {params.outdir}/*.gz | xargs -n1 basename | sed 's/\\([0-9]\\)_.*/\\1/' | grep -F -f - {input.taxonomy} > {output.derep_taxonomy}
		rm {params.indir}/*.gz
		"""

# collect all genomes (from GTDB and NCBI) and associated taxonomy file in common directory
