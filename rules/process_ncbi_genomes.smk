rule get_genomes_ncbi:
	output: 
		links = config["rdir"] + "/{library_name}/assembly_url_genomic.txt",
		summary = config["rdir"] + "/{library_name}/assembly_summary_combined.txt"
	params:
		library_dir = config["rdir"],
		library_name = "{library_name}",
		assembly_level = lambda wildcards: config["assembly_level"][wildcards.library_name],
		script = config["wdir"] + "/scripts/prepare_files_ncbi.R",
		ncbi_server = config["ncbi_server"]
	conda:
                config["wdir"] + "/envs/r.yaml"
	threads: config["download_threads"]
	log:
		config["rdir"] + "/logs/get_genomes_ncbi_{library_name}.log"
	shell:
		"""
		wget -O "{params.library_dir}/{params.library_name}/assembly_summary_refseq.txt" "{params.ncbi_server}/genomes/refseq/{params.library_name}/assembly_summary.txt" &>> {log}
		wget -O "{params.library_dir}/{params.library_name}/assembly_summary_genbank.txt" "{params.ncbi_server}/genomes/genbank/{params.library_name}/assembly_summary.txt" &>> {log}
		{params.script} -r "{params.library_dir}/{params.library_name}/assembly_summary_refseq.txt" -g "{params.library_dir}/{params.library_name}/assembly_summary_genbank.txt" -a "{params.assembly_level}" -c {threads} -o "{output.links}" -s "{output.summary}" &>> {log}
		"""

rule download_genomes_ncbi:
	input:
		url = config["rdir"] + "/{library_name}/assembly_url_genomic.txt"
	output:
		download_complete = config["rdir"] + "/{library_name}/genomes/done"
	params:
		outdir = config["rdir"] + "/{library_name}/genomes"
	threads: config["download_threads"]
	conda:
		config["wdir"] + "/envs/download.yaml"
	log:
                config["rdir"] + "/logs/download_genomes_ncbi_{library_name}.log"
	shell:
		"""
		cut -f1 {input.url} > "{params.outdir}/links.list"
		aria2c -i "{params.outdir}/links.list" -c -l "{params.outdir}/links.log" --dir {params.outdir} --max-tries=20 --retry-wait=5 --max-connection-per-server=1 --max-concurrent-downloads={threads} &>> {log}
		# We need to verify all files are there
		cat "{params.outdir}/links.list" | xargs -n 1 basename | sort > "{params.outdir}/tmp1"
		find {params.outdir} -type f -name '*.gz' | xargs -n 1 basename | sort > "{params.outdir}/tmp2"
		if diff "{params.outdir}/tmp1" "{params.outdir}/tmp2" 
		then
		  touch "{params.outdir}/done"
		fi
		rm "{params.outdir}/links.list" "{params.outdir}/links.log" "{params.outdir}/tmp1" "{params.outdir}/tmp2"
		"""

rule calculate_contig_stats:
	input:
		download_genomes = config["rdir"] + "/{library_name}/genomes/done",
	output:
		contig_stats = config["rdir"] + "/{library_name}/metadata/genome_metadata.txt"
	params:
		genome_dir = config["rdir"] + "/{library_name}/genomes"
	conda:
		config["wdir"] + "/envs/bbmap.yaml"
	log:
		config["rdir"] + "/logs/calculate_contig_stats_{library_name}.log"
	shell:
		"""
		statswrapper.sh "{params.genome_dir}/*.gz" format=5 > "{params.genome_dir}/../metadata/bb_out.tmp" 2> {log}
		paste "{params.genome_dir}/../metadata/bb_out.tmp" <(cut -f11 "{params.genome_dir}/../metadata/bb_out.tmp" | sed 's/.*\///' | cut -d'_' -f1,2 | sed 's/filename/accession/') | sed 's/contig_bp/genome_size/'> {output.contig_stats}
		rm "{params.genome_dir}/../metadata/bb_out.tmp"
		"""

rule get_taxdump:
	output:
		nodes = config["rdir"] + "/ncbi_taxdump/nodes.dmp",
		names = config["rdir"] + "/ncbi_taxdump/names.dmp"
	params:
		outdir = config["rdir"]
	conda:
		config["wdir"] + "/envs/download.yaml"
	log:
                config["rdir"] + "/logs/get_taxdump.log"
	shell:
		"""
		aria2c --max-tries=20 --retry-wait=5 --dir {params.outdir} http://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz &>> {log}
		mkdir -p "{params.outdir}/ncbi_taxdump/"
		tar -xzvf "{params.outdir}/new_taxdump.tar.gz" --directory "{params.outdir}/ncbi_taxdump/" &>> {log}
		rm "{params.outdir}/new_taxdump.tar.gz"
		"""

rule get_taxpath:
	input:
		genomes = config["rdir"] + "/{library_name}/assembly_summary_combined.txt",
		nodes = config["rdir"] + "/ncbi_taxdump/nodes.dmp",
		names = config["rdir"] + "/ncbi_taxdump/names.dmp"
	output:
		all_ranks = config["rdir"] + "/{library_name}/assembly_taxonomy_all_ranks.txt"
	params:
		taxdir = config["rdir"] + "/ncbi_taxdump"
	conda:
		config["wdir"] + "/envs/taxonkit.yaml"
	shell:
		"""
		cut -f7 {input.genomes} | taxonkit lineage -R --data-dir {params.taxdir} | paste <(cut -f1 {input.genomes}) - > {output.all_ranks}
		"""

rule format_taxpath:
	input:
		all_ranks = config["rdir"] + "/{library_name}/assembly_taxonomy_all_ranks.txt"
	output:
		taxonomy = config["rdir"] + "/{library_name}/assembly_taxonomy.txt"
	params:
		script = config["wdir"] + "/scripts/format_ncbi_taxpath.R",
		library = "{library_name}"
	conda:
		config["wdir"] + "/envs/r.yaml"
	log:
                config["rdir"] + "/logs/get_taxpath_{library_name}.log"
	shell:
		"""
		{params.script} -i {input.all_ranks} -t {params.library} -o {output.taxonomy} &>> {log}
		"""

rule filter_assemblies:
	input:
		taxonomy = config["rdir"] + "/{library_name}/assembly_taxonomy.txt",
		contig_stats = config["rdir"] + "/{library_name}/metadata/genome_metadata.txt",
		assembly_summary = config["rdir"] + "/{library_name}/assembly_summary_combined.txt"
	output:
		tax_filt = config["rdir"] + "/{library_name}/assembly_taxonomy_filtered.txt"
	params:
		script = config["wdir"] + "/scripts/filter_genome_assemblies.R",
		fields = lambda wildcards: config["filter_ncbi_fields"][wildcards.library_name],
		abs_cut = lambda wildcards: config["filter_ncbi_abs"][wildcards.library_name],
		rel_cut = lambda wildcards: config["filter_ncbi_rel"][wildcards.library_name],
		cut_dif = lambda wildcards: config["filter_ncbi_dif"][wildcards.library_name],
		cut_sign = lambda wildcards: config["filter_ncbi_sign"][wildcards.library_name]
	conda:
		config["wdir"] + "/envs/r.yaml"
	log:
		config["rdir"] + "/logs/filter_ncbi_assemblies_{library_name}.log"
	shell:
		"""
		{params.script} -t "{input.taxonomy}" -m "{input.contig_stats}" -n "{input.assembly_summary}" -p "{params.fields}" --absolute_cutoff="{params.abs_cut}" -r "{params.rel_cut}" -d "{params.cut_dif}" -s "{params.cut_sign}" -o {output.tax_filt} &>> {log}
		"""

rule add_custom_ncbi_pre_derep:
	input:
		tax_ncbi = config["rdir"] + "/{library_name}/assembly_taxonomy_filtered.txt"
	output:
		tax_added = config["rdir"] + "/{library_name}/assembly_taxonomy_added.txt"
	params:
		add = lambda wildcards: config["custom_ncbi_pre_derep"][wildcards.library_name],
		gendir = config["rdir"] + "/{library_name}/genomes"
	shell:
		"""
		if [[ "{params.add}" != "" ]]
		then
		  cut -f3 {params.add} | while read line
		  do
		    ln -sf "$line" {params.gendir}
		  done
		  cat {input.tax_ncbi} > {output.tax_added}
		  cut -f1,2 {params.add} >> {output.tax_added}
		else
		  cp {input.tax_ncbi} {output.tax_added}
		fi
		"""
		
localrules: derep_ncbi

rule derep_ncbi:
	input:
		download_complete = config["rdir"] + "/{library_name}/genomes/done",
		taxonomy = config["rdir"] + "/{library_name}/assembly_taxonomy_added.txt"
	output:
		tax_select = config["rdir"] + "/{library_name}/assembly_taxonomy_select.txt",
		derep_meta = config["rdir"] + "/{library_name}/assembly-derep-genomes_results.tsv"
	params:
		indir = config["rdir"] + "/{library_name}/genomes",
		outdir = config["rdir"] + "/{library_name}/derep_genomes",
		z_threshold = lambda wildcards: config["z_threshold_ncbi"][wildcards.library_name],
		m_threshold = lambda wildcards: config["m_threshold_ncbi"][wildcards.library_name],
		ani_fraglen = lambda wildcards: config["ani_fraglen_ncbi"][wildcards.library_name],
		derep_db = config["rdir"] + "/{library_name}/derep_genomes/derep_db",
		derep_lineage = lambda wildcards: config["derep_lineage_exclude"][wildcards.library_name],
		derep_slurm = config["wdir"] + "/config/cluster_derep.yaml",
		derep_chunks = lambda wildcards: config["ncbi_derep_chunks"][wildcards.library_name]
	threads: config["derep_threads"]
	conda:
		config["wdir"] + "/envs/derep.yaml"
	log:
		config["rdir"] + "/logs/derep_ncbi_{library_name}.log"
	shell:
		"""
		mkdir -p {params.outdir}
		# remove exlcuded lineages from taxonomy table (alternatively supply file to derep command with taxa to keep)
		if [[ "{params.derep_lineage}" != "" ]]
		then
		  grep -v "{params.derep_lineage}" {input.taxonomy} > {output.tax_select}
		else
		  cp {input.taxonomy} {output.tax_select}
		fi
		cd {params.outdir}
		derepG --threads {threads} --in-dir {params.indir} --prefix ../assembly --taxa {output.tax_select} --tmp ./ --db {params.derep_db} --threshold {params.z_threshold} --mash-threshold {params.m_threshold} --ani-fraglen-fraction {params.ani_fraglen} --debug --slurm-config {params.derep_slurm} --chunk-size {params.derep_chunks} --slurm-arr-size 10000 &>> {log}
		# mv *-derep-genomes_results.tsv {output.derep_meta}
		# do not delete redundant genomes until DB workflow is finished, work with soft links for remaining steps
		"""

rule collect_ncbi_genomes:
	input:
		derep_meta = config["rdir"] + "/{library_name}/assembly-derep-genomes_results.tsv",
		taxonomy = config["rdir"] + "/{library_name}/assembly_taxonomy_select.txt",
		url = config["rdir"] + "/{library_name}/assembly_url_genomic.txt"
	output:
		tax = config["rdir"] + "/tax_combined/{library_name}_derep_taxonomy.txt"
	params:
		gendir = config["rdir"] + "/{library_name}/genomes",
		indir = config["rdir"] + "/{library_name}/derep_genomes",
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
		if [[ $(cut -f2 {input.taxonomy} | sort -t$'\\t' | uniq | wc -l) != $(cut -f1 {input.derep_meta} | sed '1d' | sort -t$'\\t' | uniq | wc -l) ]]
		then
		  cut -f2 {input.taxonomy} | sort -t$'\\t' | uniq | grep -v -F -f <(cut -f1 {input.derep_meta} | sed '1d' | sort -t$'\\t' | uniq) | grep -F -f - {input.taxonomy} > "{params.indir}/../tmp"
		  cut -f1 "{params.indir}/../tmp" | grep -F -f - <(cut -f1 {input.url}) | xargs -n1 basename | while read line
		  do
		    ln -sf "{params.gendir}/$line" {params.outdir}
		  done
		  cat "{params.indir}/../tmp" >> {output.tax}
		  rm "{params.indir}/../tmp"
		fi
		"""

if config["custom_ncbi_post_derep"]:
	rule add_custom_ncbi_post_derep:
		output:
			tax = config["rdir"] + "/tax_combined/euk_custom_post_derep_taxonomy.txt"
		params:
			add = config["custom_ncbi_post_derep"],
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
