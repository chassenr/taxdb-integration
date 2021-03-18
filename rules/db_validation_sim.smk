rule sim_select_genomes:
	input:
		tax_highres = config["rdir"] + "/tax_combined/derep_taxonomy_good.txt",
		tax_coarse = config["cdir"] + "/tax_coarse_good.txt",
		tax_gtdb = config["rdir"] + "/gtdb/metadata/gtdb_taxonomy.txt",
		tax_ncbi = expand(config["rdir"] + "/{library_name}/assembly_taxonomy.txt", library_name = LIBRARY_NAME),
		tax_checkv = config["rdir"] + "/checkv/checkv_taxonomy.txt"
	output:
		db_accnos = config["sdir"] + "/sim_out/db.accnos",
		tax_sim = config["sdir"] + "/sim_out/genome_candidates.txt",
		tax_select = config["sdir"] + "/sim_out/genomes_select.txt"
	params:
		simdir = config["sdir"] + "/sim_out"
	shell:
		"""
		cat {input.tax_highres} {input_coarse} | cut -f1 | sort | uniq > {output.db_accnos}
		cat {input.tax_gtdb} {input.tax_checkv} {input.tax_ncbi} | grep -v -F -f {output.db_accnos} | sort -t$'\\t' -k2,2 | sort -t$'\\t' -u -k2,2 --merge > {output.tax_sim}
		cut -f2 {output.tax_sim} | cut -d';' -f-2 | sort -t$'\\t' | uniq > "{params.simdir}/phylum_list.txt"
		while read line
		do
		  grep "${line}" {output.tax_sim} | shuf -n 3 >> "{params.simdir}/tmp_select.txt"
		done < "{params.simdir}/phylum_list.txt"
		sort -t$'\\t' -k2,2 "{params.simdir}/tmp_select.txt" > {output.tax_select}
		rm "{params.simdir}/tmp_select.txt"
		"""

rule sim_genome_download:
	input:
		tax_select = config["sdir"] + "/sim_out/genomes_select.txt",
		gtdb_url = config["rdir"] + "/gtdb/metadata/gtdb_download_info.txt",
		ncbi_url = expand(config["rdir"] + "/{library_name}/assembly_url_genomic.txt", library_name = LIBRARY_NAME),
		checkv_genomes = config["rdir"] + "/checkv/genomes/done",
		tax_checkv = config["rdir"] + "/checkv/checkv_taxonomy.txt"
	output:
		sim_url = config["sdir"] + "/sim_out/genome_links.txt",
		sim_genomes = config["sdir"] + "/sim_out/genomes/done"
	params:
		outdir = config["sdir"] + "/sim_out/genomes"
	threads: config["download_threads"]
	conda:
		config["wdir"] + "/envs/download.yaml"
	log:
		config["sdir"] + "logs/sim_genomes_download.log"
	shell:
		"""
		cut -f2 {input.gtdb_url} | cat - {input.ncbi_url} | grep -F -f <(cut -f1 {input.tax_select}) - > {output.sim_url}
		# append any gtdb genomes with mismatches in assembly versions
		sed 's/.*\///' {output.sim_url} | cut -d'_' -f2 | grep -v -F -f - {input.tax_select} | cut -f1 | cut -d'_' -f2 | sed 's/\.[0-9]\+//' | grep -F -f - <(cut -f2 {input.gtdb_url}) >> {output.sim_url}
		aria2c -i {output.sim_url} -c -l "{params.outdir}/links.log" --dir {params.outdir} --max-tries=20 --retry-wait=5 --max-connection-per-server=1 --max-concurrent-downloads={threads} &>> {log}
                # We need to verify all files are there
                cat {output.sim_url} | xargs -n 1 basename | sort > "{params.outdir}/tmp1"
                find {params.outdir} -type f -name '*.gz' | xargs -n 1 basename | sort > "{params.outdir}/tmp2"
                if diff "{params.outdir}/tmp1" "{params.outdir}/tmp2"
                then
                  touch {output.sim_genomes}
                fi
                rm "{params.outdir}/links.log" "{params.outdir}/tmp1" "{params.outdir}/tmp2"
		gunzip {params.outdir}/*.gz
		"""

rule sim_community:
	input:
		tax_select = config["sdir"] + "/sim_out/genomes_select.txt",
		sim_url = config["sdir"] + "/sim_out/genome_links.txt",
		sim_genomes = config["sdir"] + "/sim_out/genomes/done"
	output:
		sim_table = config["sdir"] + "/sim_out/sim_genomes_table.txt",
		sim_abund = config["sdir"] + "/sim_out/sim_PE_{RL}_out_abund.txt"
	params:
		gendir = config["sdir"] + "/sim_out/genomes/",
		prefix = config["sdir"] + "/sim_out/sim_PE_out"
	conda:
		config["wdir"] + "/envs/mgsim.yaml"
	log:
		
	shell:
		"""
		# format genome table for MGSIM
		paste <(cat {input.sim_url} | xargs -n 1 basename | sed -e 's/\.gz//' -e 's/^/{params.gendir}/') <(cat {input.sim_url} | xargs -n 1 basename | cut -d'_' -f2 | sed 's/\.[0-9]\+//') | sort -t$'\\t' -k2,2 > "{params.gendir}../tmp1" 
		cut -f1 {input.tax_select} | cut -d'_' -f2 | sed 's/\.[0-9]\+//' | paste {input.tax_select} - | sort -t$'\\t' -k3,3 > "{params.gendir}../tmp2"
		if diff <(cut -f2 "{params.gendir}../tmp1") <(cut -f3 "{params.gendir}../tmp2")
		then		
		  echo -e "Accnos\tTaxon\tFasta" > {output.sim_table}
		  cut -f2,1 tmp2 | paste - <(cut -f1 tmp1) >> {output.sim_table}
		fi
		rm "{params.gendir}../tmp1" "{params.gendir}../tmp2"
		# prepare community percentages
		MGSIM communities {output.sim_table} {params.prefix}
		"""

rule sim_add_random:
	input:
		sim_abund = config["sdir"] + "/sim_out/sim_PE_out_abund.txt"
	output:
		sim_abund_wr = config["sdir"] + "/sim_out/sim_PE_out_abund_wr.txt"
	params:
		base_gen = config["random_link"],
		simdir = config["sdir"] + "/sim_out",
		script = 
	conda:
		config["wdir"] + "/envs/r.yaml"
	log:
		
	shell:
		"""
		wget wget -P {params.outdir} {params.checkv_link} &>> {log}
		{params.script}
		"""

rule sim_reads_modern:
	input:
		sim_abund_wr = config["sdir"] + "/sim_out/sim_PE_out_abund_wr.txt",
		sim_table = config["sdir"] + "/sim_out/sim_genomes_table.txt"
	output:
		sim_fq_R1 = config["sdir"] + "/sim_out/sim_PE_{RL}_out/1/R1.fq".
		sim_fq_R2 = config["sdir"] + "/sim_out/sim_PE_{RL}_out/1/R2.fq"
	params:
		outdir = config["sdir"] + "/sim_out/sim_PE_{RL}_out"
		seq_depth = config["sim_seqdepth"],
		read_len = "{RL}"
	conda:
		config["wdir"] + "/envs/mgsim.yaml"
	threads: config["sim_threads"]
	log:
		
	shell:
		"""
		MGSIM reads --sr-seq-depth {params.seq_depth} --art-len={params.read_len} --art-mflen=0 -n {threads} {input.sim_table} {input.sim_abund_wr} {params.outdir}
		"""


