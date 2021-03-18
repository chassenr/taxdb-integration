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
	log:
		config["rdir"] + "/logs/get_genomes_ncbi_{library_name}.log"
	shell:
		"""
		wget -O "{params.library_dir}/{params.library_name}/assembly_summary_refseq.txt" "{params.ncbi_server}/genomes/refseq/{params.library_name}/assembly_summary.txt" &>> {log}
		wget -O "{params.library_dir}/{params.library_name}/assembly_summary_genbank.txt" "{params.ncbi_server}/genomes/genbank/{params.library_name}/assembly_summary.txt" &>> {log}
		{params.script} -r "{params.library_dir}/{params.library_name}/assembly_summary_refseq.txt" -g "{params.library_dir}/{params.library_name}/assembly_summary_genbank.txt" -a "{params.assembly_level}" -o "{output.links}" -s "{output.summary}" &>> {log}
		"""

rule download_genomes_ncbi:
	input:
		url = config["rdir"] + "/{library_highres}/assembly_url_genomic.txt"
	output:
		download_complete = config["rdir"] + "/{library_highres}/genomes/done"
	params:
		outdir = config["rdir"] + "/{library_highres}/genomes"
	threads: config["download_threads"]
	conda:
		config["wdir"] + "/envs/download.yaml"
	log:
                config["rdir"] + "/logs/download_genomes_ncbi_{library_highres}.log"
	shell:
		"""
		aria2c -i {input.url} -c -l "{params.outdir}/links.log" --dir {params.outdir} --max-tries=20 --retry-wait=5 --max-connection-per-server=1 --max-concurrent-downloads={threads} &>> {log}
		# We need to verify all files are there
		cat {input.url} | xargs -n 1 basename | sort > "{params.outdir}/tmp1"
		find {params.outdir} -type f -name '*.gz' | xargs -n 1 basename | sort > "{params.outdir}/tmp2"
		if diff "{params.outdir}/tmp1" "{params.outdir}/tmp2" 
		then
		  touch "{params.outdir}/done"
		fi
		rm "{params.outdir}/links.log" "{params.outdir}/tmp1" "{params.outdir}/tmp2"
		"""

rule calculate_contig_stats:
	input:
		download_genomes = config["rdir"] + "/{library_highres}/genomes/done",
	output:
		contig_stats = config["rdir"] + "/{library_highres}/metadata/genome_metadata.txt"
	params:
		genome_dir = config["rdir"] + "/{library_highres}/genomes"
	conda:
		config["wdir"] + "/envs/bbmap.yaml"
	log:
		config["rdir"] + "/logs/calculate_contig_stats_{library_highres}.log"
	shell:
		"""
		statswrapper.sh "{params.genome_dir}/*.gz" format=5 > "{params.genome_dir}/../metadata/bb_out.tmp" 2> {log}
		paste "{params.genome_dir}/../metadata/bb_out.tmp" <(cut -f11 "{params.genome_dir}/../metadata/bb_out.tmp" | sed 's/.*\///' | cut -d'_' -f1,2 | sed 's/filename/accession/') | sed 's/contig_bp/genome_size/'> {output.contig_stats}
		rm "{params.genome_dir}/../metadata/bb_out.tmp"
		"""

rule download_assembly_stats:
	input:
		url = config["rdir"] + "/{library_coarse}/assembly_url_genomic.txt"
	output:
		stats_url = config["rdir"] + "/{library_coarse}/assembly_url_contig_stats.txt",
		download_complete = config["rdir"] + "/{library_coarse}/contig_stats/done"
	params:
		outdir = config["rdir"] + "/{library_coarse}/contig_stats"
	threads: config["download_threads"]
	conda:
		config["wdir"] + "/envs/download.yaml"
	log:
		config["rdir"] + "/logs/download_assembly_stats_{library_coarse}.log"
	shell:
		"""
		sed 's/_genomic\.fna\.gz/_assembly_stats\.txt/' {input.url} > {output.stats_url}
		aria2c -i {output.stats_url} -c -l "{params.outdir}/links.log" --dir {params.outdir} --max-tries=20 --retry-wait=5 --max-connection-per-server=1 --max-concurrent-downloads={threads} &>> {log}
		# We need to verify all files are there
		cat {output.stats_url} | xargs -n 1 basename | sort > "{params.outdir}/tmp1"
		find {params.outdir} -type f -name '*.txt' | xargs -n 1 basename | sort > "{params.outdir}/tmp2"
		if diff "{params.outdir}/tmp1" "{params.outdir}/tmp2"
		then
		  touch "{params.outdir}/done"
		fi
		rm "{params.outdir}/links.log" "{params.outdir}/tmp1" "{params.outdir}/tmp2"
		"""

rule parse_assembly_metadata:
	input:
		download_contig_stats = config["rdir"] + "/{library_coarse}/contig_stats/done",
		summary = config["rdir"] + "/{library_coarse}/assembly_summary_combined.txt"
	output:
		metadata = config["rdir"] + "/{library_coarse}/metadata/genome_metadata.txt"
	params:
		outdir = config["rdir"] + "/{library_coarse}/metadata",
		stats_dir = config["rdir"] + "/{library_coarse}/contig_stats",
		script = config["wdir"] + "/scripts/parse_assembly_metadata.R"
	log:
		config["rdir"] + "/logs/parse_assembly_metadata_{library_coarse}.log"
	shell:
		"""
		grep "total-length" {params.stats_dir}/*.txt | grep ":all" | sed 's/^.*\///' | sed -E 's/(GC[AF]_[0-9]+\.[0-9]+)_.*_assembly_stats\.txt:/\\1\\t/' > "{params.outdir}/contig_stats.tmp"
		echo -e "accession\\tgenome_size" > {output.metadata}
		cut -f1,7 "{params.outdir}/contig_stats.tmp" >> {output.metadata}
		rm "{params.outdir}/contig_stats.tmp"
		# optional: delete contig_stats directory again (not implemented at the moment).
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
		nodes = config["rdir"] + "/ncbi_taxdump/nodes.dmp",
		names = config["rdir"] + "/ncbi_taxdump/names.dmp",
		genomes = config["rdir"] + "/{library_name}/assembly_summary_combined.txt"
	output:
		taxonomy = config["rdir"] + "/{library_name}/assembly_taxonomy.txt"
	params:
		outdir = config["rdir"],
		script = config["wdir"] + "/scripts/get_taxpath.R"
	conda:
		config["wdir"] + "/envs/r.yaml"
	log:
                config["rdir"] + "/logs/get_taxpath_{library_name}.log"
	shell:
		"""
		{params.script} -i {input.genomes} -t "{params.outdir}/ncbi_taxdump" -s "{params.outdir}/ncbi_taxdump/accessionTaxa.sql" -o {output.taxonomy} &>> {log}
		"""

rule filter_assemblies:
	input:
		taxonomy = config["rdir"] + "/{library_highres}/assembly_taxonomy.txt",
		contig_stats = config["rdir"] + "/{library_highres}/metadata/genome_metadata.txt"
	output:
		tax_filt = config["rdir"] + "/{library_highres}/assembly_taxonomy_filtered.txt"
	params:
		script = config["wdir"] + "/scripts/filter_genome_assemblies.R",
		fields = lambda wildcards: config["filter_ncbi_fields"][wildcards.library_highres],
		abs_cut = lambda wildcards: config["filter_ncbi_abs"][wildcards.library_highres],
		rel_cut = lambda wildcards: config["filter_ncbi_rel"][wildcards.library_highres],
		cut_dif = lambda wildcards: config["filter_ncbi_dif"][wildcards.library_highres],
		cut_sign = lambda wildcards: config["filter_ncbi_sign"][wildcards.library_highres]
	conda:
		config["wdir"] + "/envs/r.yaml"
	log:
		config["rdir"] + "/logs/filter_ncbi_assemblies_{library_highres}.log"
	shell:
		"""
		{params.script} -t "{input.taxonomy}" -m "{input.contig_stats}" -p "{params.fields}" --absolute_cutoff="{params.abs_cut}" -r "{params.rel_cut}" -d "{params.cut_dif}" -s "{params.cut_sign}" -o {output.tax_filt} &>> {log}
		"""

# to add cusom assemblies, manually include genome files in the respective directories and provide the taxonomy file name in the config file for custom_ncbi
rule add_custom_assemblies:
	input:
		tax_ncbi = config["rdir"] + "/{library_highres}/assembly_taxonomy_filtered.txt"
	output:
		tax_added = config["rdir"] + "/{library_highres}/assembly_taxonomy_added.txt"
	params:
		add = lambda wildcards: config["custom_ncbi"][wildcards.library_highres],
		dir = config["rdir"] + "/{library_highres}/genomes"
	shell:
		"""
		if [[ "{params.add}" != "" ]]
		then
		  ls -1 {params.dir} | grep -F -f <(cut -f1 {params.add}) > {params.dir}/tmp
		  if [[ "$(wc -l < {params.dir}/tmp)" -eq "$(wc -l < {params.add})" ]]
		  then
		    cat {input.tax_ncbi} {params.add} > {output.tax_added}
		  fi
		  rm {params.dir}/tmp
		else
		  cp {input.tax_ncbi} {output.tax_added}
		fi
		"""
		
localrules: derep_ncbi

rule derep_ncbi:
	input:
		download_complete = config["rdir"] + "/{library_highres}/genomes/done",
		taxonomy = config["rdir"] + "/{library_highres}/assembly_taxonomy_added.txt"
	output:
		derep_meta = config["rdir"] + "/{library_highres}/derep_taxonomy_meta.txt"
	params:
		dir = config["rdir"] + "/{library_highres}",
		indir = config["rdir"] + "/{library_highres}/genomes",
		outdir = config["rdir"] + "/{library_highres}/derep_genomes",
		z_threshold = lambda wildcards: config["z_threshold_ncbi"][wildcards.library_highres],
		m_threshold = lambda wildcards: config["m_threshold_ncbi"][wildcards.library_highres],
		ani_fraglen = lambda wildcards: config["ani_fraglen_ncbi"][wildcards.library_highres],
		derep_db = config["rdir"] + "/{library_highres}/derep_genomes/derep_db",
		derep_lineage = lambda wildcards: config["derep_lineage_exclude"][wildcards.library_highres],
		derep_slurm = config["wdir"] + "/config/cluster_derep.yaml",
		derep_chunks = lambda wildcards: config["ncbi_derep_chunks"][wildcards.library_highres]
	threads: config["derep_threads"]
	conda:
		config["wdir"] + "/envs/derep.yaml"
	log:
		config["rdir"] + "/logs/derep_ncbi_{library_highres}.log"
	shell:
		"""
		mkdir -p {params.outdir}
		# remove exlcuded lineages from taxonomy table (alternatively supply file to derep command with taxa to keep)
		if [[ "{params.derep_lineage}" != "" ]]
		then
		  grep -v "{params.derep_lineage}" {input.taxonomy} > "{params.dir}/assembly_taxonomy_select.txt"
		else
		  cp {input.taxonomy} "{params.dir}/assembly_taxonomy_select.txt"
		fi
		cd {params.outdir}
		derepG --threads {threads} --in-dir {params.indir} --taxa "{params.dir}/assembly_taxonomy_select.txt" --tmp ./ --db {params.derep_db} --threshold {params.z_threshold} --mash-threshold {params.m_threshold} --ani-fraglen-fraction {params.ani_fraglen} --debug --slurm-config {params.derep_slurm} --chunk-size {params.derep_chunks} --slurm-arr-size 10000 &>> {log}
		mv *derep-genomes_results.tsv {output.derep_meta}
		# do not delete redundant genomes until DB workflow is finished, work with soft links for remaining steps
		"""

rule collect_ncbi_genomes:
	input:
		derep_meta = config["rdir"] + "/{library_highres}/derep_taxonomy_meta.txt",
		taxonomy = config["rdir"] + "/{library_highres}/assembly_taxonomy_select.txt",
		url = config["rdir"] + "/{library_highres}/assembly_url_genomic.txt"
	output:
		tax = config["rdir"] + "/tax_combined/{library_highres}_derep_taxonomy.txt"
	params:
		gendir = config["rdir"] + "/{library_highres}/genomes",
		indir = config["rdir"] + "/{library_highres}/derep_genomes",
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
		  cut -f1 "{params.indir}/../tmp" | grep -F -f - {input.url} | xargs -n1 basename | while read line
		  do
		    ln -sf "{params.gendir}/$line" {params.outdir}
		  done
		  cat "{params.indir}/../tmp" >> {output.tax}
		  rm "{params.indir}/../tmp"
		fi
		"""

rule masking_ncbi:
	input:
		file_list = config["rdir"] + "/kraken2_genomes/file_names_derep_genomes.txt",
		ncbi = config["rdir"] + "/tax_combined/{library_highres}_derep_taxonomy.txt",
		nodes = config["rdir"] + "/kraken2_db/taxonomy/nodes.dmp",
		names = config["rdir"] + "/kraken2_db/taxonomy/names.dmp"
	output:
		fasta = config["rdir"] + "/kraken2_db/tmp/{library_highres}_library.fna"
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	threads: config["masking_threads"]
	shell:
		"""
		cut -f1 {input.ncbi} | grep -F -f - {input.file_list} | parallel -j{threads} 'dustmasker -in {{}} -outfmt fasta' | sed -e '/^>/!s/[a-z]/x/g' >> {output.fasta}
		"""

rule prelim_map_ncbi:
	input:
		fasta = config["rdir"] + "/kraken2_db/tmp/{library_highres}_library.fna"
	output:
		map = config["rdir"] + "/kraken2_db/tmp/{library_highres}_prelim_map.txt"
	params:
		libdir = config["rdir"] + "/kraken2_db/tmp",
		libname = "{library_highres}"
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	shell:
		"""
		LC_ALL=C grep '^>' {input.fasta} | sed 's/^>//' > {params.libdir}/tmp_{params.libname}.accnos
		NSEQ=$(wc -l {params.libdir}/tmp_{params.libname}.accnos | cut -d' ' -f1)
		printf 'TAXID\\n%.0s' $(seq 1 $NSEQ) | paste - {params.libdir}/tmp_{params.libname}.accnos | paste - <(cut -d'|' -f3 {params.libdir}/tmp_{params.libname}.accnos) > {output.map}
		rm {params.libdir}/tmp_{params.libname}.accnos
		"""

if config["kingdoms_highres"]:
	rule filter_contam_ncbi:
		input:
			id_contam = config["rdir"] + "/decontamination/contam_id.accnos",
			fasta = config["rdir"] + "/kraken2_db/tmp/{library_highres}_library.fna",
			map = config["rdir"] + "/kraken2_db/tmp/{library_highres}_prelim_map.txt"
		output:
			fasta_contam = config["rdir"] + "/decontamination/{library_highres}_library_contam.fna",
			fasta_noncontam = config["rdir"] + "/kraken2_db/library/{library_highres}/library.fna",
			map_noncontam = config["rdir"] + "/kraken2_db/library/{library_highres}/prelim_map.txt"
		conda:
			config["wdir"] + "/envs/bbmap.yaml"
		log:
			config["rdir"] + "/logs/{library_highres}_contam_filter.log"
		shell:
			"""
			filterbyname.sh in={input.fasta} out={output.fasta_contam} names={input.id_contam} include=t
			filterbyname.sh in={input.fasta} out={output.fasta_noncontam} names={input.id_contam} include=f
			grep -v -F -f {input.id_contam} {input.map} > {output.map_noncontam}
			"""

	rule remove_contam_ncbi:
		input:
			contam = config["rdir"] + "/decontamination/highres_db_conterm_prediction_filt",
			fasta_contam = config["rdir"] + "/decontamination/{library_highres}_library_contam.fna",
			fasta_noncontam = config["rdir"] + "/kraken2_db/library/{library_highres}/library.fna",
			map_noncontam = config["rdir"] + "/kraken2_db/library/{library_highres}/prelim_map.txt",
			fasta_tmp = config["rdir"] + "/kraken2_db/tmp/{library_highres}_library.fna",
			map_tmp = config["rdir"] + "/kraken2_db/tmp/{library_highres}_prelim_map.txt"
		output:
			cleaned_fasta = config["rdir"] + "/decontamination/{library_highres}_cleaned.fna",
			cleaned_map = config["rdir"] + "/decontamination/{library_highres}_cleaned_map.txt"
		params:
			script = config["wdir"] + "/scripts/remove_contamination.R",
			contam_dir = config["rdir"] + "/decontamination",
			libname = "{library_highres}"
		conda:
			config["wdir"] + "/envs/r.yaml"
		log:
			config["rdir"] + "/logs/{library_highres}_contam_remove.log"
		shell:
			"""
			{params.script} -i {input.fasta_contam} -c {input.contam} -o {output.cleaned_fasta} &>> {log}
			LC_ALL=C grep '^>' {output.cleaned_fasta} | sed 's/^>//' > "{params.contam_dir}/tmp_{params.libname}.accnos"
			NSEQ=$(wc -l "{params.contam_dir}/tmp_{params.libname}.accnos" | cut -d' ' -f1)
			printf 'TAXID\\n%.0s' $(seq 1 $NSEQ) | paste - "{params.contam_dir}/tmp_{params.libname}.accnos" | paste - <(cut -d'|' -f3 "{params.contam_dir}/tmp_{params.libname}.accnos") > {output.cleaned_map}
			rm "{params.contam_dir}/tmp_{params.libname}.accnos"
			cat {output.cleaned_map} >> {input.map_noncontam}
			cat {output.cleaned_fasta} >> {input.fasta_noncontam}
			rm {input.fasta_tmp} {input.map_tmp}
			"""

else:
	rule move_library_ncbi:
		input:
			fasta = config["rdir"] + "/kraken2_db/tmp/{library_highres}_library.fna",
			map = config["rdir"] + "/kraken2_db/tmp/{library_highres}_prelim_map.txt"
		output:
			fasta = config["rdir"] + "/kraken2_db/library/{library_highres}/library.fna",
			map = config["rdir"] + "/kraken2_db/library/{library_highres}/prelim_map.txt"
		shell:
			"""
			mv {input.fasta} {output.fasta}
			mv {input.map} {output.map}
			"""

