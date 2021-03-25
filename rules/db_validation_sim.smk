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
		cat {input.tax_highres} {input.tax_coarse} | cut -f1 | sort | uniq > {output.db_accnos}
		cat {input.tax_gtdb} {input.tax_ncbi} | grep -v -F -f {output.db_accnos} | sort -t$'\\t' -k2,2 | sort -t$'\\t' -u -k2,2 --merge > {output.tax_sim}
		sed 's/\\t/__accnos\\t/' {input.tax_checkv} | grep -v -F -f <(sed 's/$/__accnos/' {output.db_accnos}) | sed 's/__accnos\\t/\\t/' | sort -t$'\\t' -k2,2 | sort -t$'\\t' -u -k2,2 --merge >> {output.tax_sim}
		cut -f2 {output.tax_sim} | cut -d';' -f1,2 | sort -t$'\\t' | uniq | grep -v "p__Eukaryota" > "{params.simdir}/phylum_list.txt"
		while read line
		do
		  grep "$line" {output.tax_sim} | shuf -n 3 >> "{params.simdir}/tmp_select.txt"
		done < "{params.simdir}/phylum_list.txt"
		sort -t$'\\t' -k2,2 "{params.simdir}/tmp_select.txt" | uniq > {output.tax_select}
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
		gendir = config["rdir"] + "/checkv/genomes",
		outdir = config["sdir"] + "/sim_out/genomes"
	threads: config["download_threads"]
	conda:
		config["wdir"] + "/envs/download.yaml"
	log:
		config["sdir"] + "/logs/sim_genomes_download.log"
	shell:
		"""
		# copy checkv genomes
		mkdir -p {params.outdir}
		cut -f1 {input.tax_checkv} | sed 's/$/\.fa/' | grep -F -f <(cut -f1 {input.tax_select} | sed 's/$/\.fa/') > "{params.outdir}/tmp_checkv.txt" 
		while read line
		do
		  cp {params.gendir}/$line {params.outdir}/$line
		done < "{params.outdir}/tmp_checkv.txt"
		# re-download gtdb and ncbi (less dependent this way)
		cut -f2 {input.gtdb_url} | cat - {input.ncbi_url} | grep -F -f <(cut -f1 {input.tax_select}) > {output.sim_url}
		# append any gtdb genomes with mismatches in assembly versions
		cat {output.sim_url} | xargs -n 1 basename | cut -d'_' -f2 | grep -v -F -f - {input.tax_select} | cut -f1 | cut -d'_' -f2 | sed -e 's/\.[0-9]\+//' -e 's/^/_/' | grep -F -f - <(cut -f2 {input.gtdb_url}) >> {output.sim_url}
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
		cat "{params.outdir}/tmp_checkv.txt" >> {output.sim_url}
		rm "{params.outdir}/tmp_checkv.txt"
		"""

rule sim_community:
	input:
		tax_select = config["sdir"] + "/sim_out/genomes_select.txt",
		sim_url = config["sdir"] + "/sim_out/genome_links.txt",
		sim_genomes = config["sdir"] + "/sim_out/genomes/done"
	output:
		sim_table = config["sdir"] + "/sim_out/sim_genomes_table.txt",
		sim_abund = config["sdir"] + "/sim_out/sim_out_abund.txt"
	params:
		gendir = config["sdir"] + "/sim_out/genomes/",
		prefix = config["sdir"] + "/sim_out/sim_out"
	conda:
		config["wdir"] + "/envs/mgsim.yaml"
	log:
		config["sdir"] + "/logs/sim_community.log"
	shell:
		"""
		# format genome table for MGSIM
		DIRSTRING=$(echo "{params.gendir}" | sed 's/\//\\\\\\//g')
		paste <(cat {input.sim_url} | xargs -n 1 basename | sed -e "s/\.gz//" -e "s/^/$DIRSTRING/") <(cat {input.sim_url} | xargs -n 1 basename | cut -d'_' -f2 | sed -e 's/\.[0-9]\+//' -e 's/\.fa//') | sort -t$'\\t' -k2,2 > "{params.gendir}../tmp1" 
		cut -f1 {input.tax_select} | cut -d'_' -f2 | sed -e 's/\.[0-9]\+//' -e 's/\.fa//' | paste {input.tax_select} - | sort -t$'\\t' -k3,3 > "{params.gendir}../tmp2"
		if diff <(cut -f2 "{params.gendir}../tmp1") <(cut -f3 "{params.gendir}../tmp2")
		then		
		  echo -e "Accnos\\tTaxon\\tFasta" > {output.sim_table}
		  cut -f2,1 "{params.gendir}../tmp2" | paste - <(cut -f1 "{params.gendir}../tmp1") >> {output.sim_table}
		fi
		rm "{params.gendir}../tmp1" "{params.gendir}../tmp2"
		# prepare community percentages
		MGSIM communities {output.sim_table} {params.prefix} &>> {log}
		"""

rule sim_add_random:
	input:
		sim_abund = config["sdir"] + "/sim_out/sim_out_abund.txt",
		sim_table = config["sdir"] + "/sim_out/sim_genomes_table.txt"
	output:
		sim_abund_wr = config["sdir"] + "/sim_out/sim_out_abund_wr.txt",
		sim_table_wr = config["sdir"] + "/sim_out/sim_genomes_table_wr.txt",
		random_fasta = config["sdir"] + "/sim_out/genomes/random.fna"
	params:
		random_template = config["random_link"],
		simdir = config["sdir"] + "/sim_out/",
		script = config["wdir"] + "/scripts/sim_add_random.R"
	conda:
		config["wdir"] + "/envs/r.yaml"
	log:
		config["sdir"] + "/logs/sim_random.log"
	shell:
		"""
		wget -P {params.simdir} {params.random_template} &>> {log}
		DIRSTRING=$(echo "{params.simdir}" | sed 's/\//\\\\\\//g')
		RANDOMGEN=$(basename {params.random_template} | sed "s/^/$DIRSTRING/")
		{params.script} -a {input.sim_abund} -g {input.sim_table} -f $RANDOMGEN -A {output.sim_abund_wr} -G {output.sim_table_wr} -r {output.random_fasta} &>> {log}
		"""

rule sim_reads_modern:
	input:
		sim_abund_wr = config["sdir"] + "/sim_out/sim_out_abund_wr.txt",
		sim_table_wr = config["sdir"] + "/sim_out/sim_genomes_table_wr.txt"
	output:
		sim_fq_R1 = config["sdir"] + "/sim_out/sim_PE_{rl}_out/1/R1.fq",
		sim_fq_R2 = config["sdir"] + "/sim_out/sim_PE_{rl}_out/1/R2.fq"
	params:
		outdir = config["sdir"] + "/sim_out/sim_PE_{rl}_out",
		seq_depth = config["sim_seqdepth"],
		read_length = "{rl}",
		insert_size = config["insert_size"]
	conda:
		config["wdir"] + "/envs/mgsim.yaml"
	threads: config["sim_threads"]
	log:
		config["sdir"] + "/logs/sim_reads_PE_rl{rl}_modern.log"
	shell:
		"""
		MGSIM reads --sr-seq-depth {params.seq_depth} --art-paired --art-len={params.read_length} --art-mflen={params.insert_size} -n {threads} {input.sim_table_wr} {input.sim_abund_wr} {params.outdir} &>> {log}
		"""

rule sim_reads_short:
	input:
		sim_abund_wr = config["sdir"] + "/sim_out/sim_out_abund_wr.txt",
		sim_table_wr = config["sdir"] + "/sim_out/sim_genomes_table_wr.txt"
	output:
		sim_fq_SE = config["sdir"] + "/sim_out/sim_SE_{serl}_out/1/R1.fq"
	params:
		outdir = config["sdir"] + "/sim_out/sim_SE_{serl}_out",
		seq_depth = config["sim_seqdepth"],
		read_length = "{serl}"
	conda:
		config["wdir"] + "/envs/mgsim.yaml"
	threads: config["sim_threads"]
	log:
		config["sdir"] + "/logs/sim_reads_SE_rl{serl}_short.log"
	shell:
		"""
		MGSIM reads --sr-seq-depth {params.seq_depth} --art-len={params.read_length} --art-mflen=0 -n {threads} {input.sim_table_wr} {input.sim_abund_wr} {params.outdir} &>> {log}
		"""

