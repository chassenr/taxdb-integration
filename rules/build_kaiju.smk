rule download_predicted_proteins:
	input:
		tax_combined = config["rdir"] + "/tax_combined/derep_taxonomy_combined.txt",
		gtdb_links = config["rdir"] + "/gtdb/metadata/gtdb_download_info.txt",
		ncbi_links = expand(config["rdir"] + "/{library_name}/assembly_url_genomic.txt", library_name = LIBRARY_NAME)
	output:
		download_complete = config["rdir"] + "/proteins_combined/done",
		protein_links = config["rdir"] + "/proteins_metadata/faa_download_links.txt"
	params:
		outdir = config["rdir"] + "/proteins_combined"
	threads: config["download_threads"]
	conda:
		config["wdir"] + "/envs/download.yaml"
	log:
		config["rdir"] + "/logs/download_proteins.log"
	shell:
		"""
		cut -f2 {input.gtdb_links} | sed 's/_genomic\.fna\.gz/_protein\.faa\.gz/' > "{params.outdir}/tmp" 
		cat {input.ncbi_links} | sed 's/_genomic\.fna\.gz/_protein\.faa\.gz/' > "{params.outdir}/tmp" 
		cut -f1 {input.tax_combined} | grep -F -f - "{params.outdir}/tmp" > {output.protein_links}                                                           
		rm "{params.outdir}/tmp"
		# implement check that download link file is consistent with number of dereplicated genomes?
		aria2c -i {output.protein_links} -c -l "{params.outdir}/links.log" --dir {params.outdir} --max-tries=20 --retry-wait=5 --max-connection-per-server=1 --max-concurrent-downloads={threads} &>> {log}
		cat {output.protein_links} | xargs -n 1 basename | sort > "{params.outdir}/tmp1"
		find {params.outdir} -type f -name '*.gz' | xargs -n 1 basename | sort > "{params.outdir}/tmp2"
		if diff "{params.outdir}/tmp1" "{params.outdir}/tmp2"
		then
		  touch "{params.outdir}/done"
		fi
		rm "{params.outdir}/links.log" "{params.outdir}/tmp1" "{params.outdir}/tmp2"
		"""

rule download_checkv_faa:
	input:
	output:
	params:
	log:
	shell:
		"""
		"""

rule format_faa_headers:
	input:
		map_gtdb = config["rdir"] + "/kraken2_db/library/gtdb/prelim_map.txt",
		map_ncbi = expand(config["rdir"] + "/kraken2_db/library/{library_name}/prelim_map.txt", library_name = LIBRARY_NAME),
		download_complete = config["rdir"] + "/proteins_combined/done",
		protein_links = config["rdir"] + "/proteins_metadata/faa_download_links.txt"
	output:
		taxid_map = config["rdir"] + "/proteins_metadata/taxid_map.txt",
		faa = config["rdir"] + "/kaiju_db/library/proteins.faa"
	params:
		outdir = config["rdir"] + "/proteins_combined",
		script = config["format_faa_script"]
	threads: config["format_threads"]
	conda:
		config["wdir"] + "/envs/kraken2.yaml" # only because it has GNU parallel
	log:
		config["rdir"] + "/logs/format_proteins.log"
	shell:
		"""
		# extract taxid per accession from kraken seqmap
		cat {input.map_gtdb} {input.map_ncbi} | cut -f2,3 | sed -e 's/_/\t/g' -e 's/\t/_/' | cut -f1,4 | uniq | sort -t$'\t' -k1,1 > {output.taxid_map}
		find {params.outdir} -type f -name '*.gz' | parallel -j{threads} '{params.script} {{}} {output.taxid_map}' >> {output.faa} 
		"""

# rule check_aa_seqs

rule build_kaiju_db:
	input:
		faa = config["rdir"] + "/kaiju_db/library/proteins.faa"
	output:
		kaiju = config["rdir"] + "/kaiju_db/proteins.fmi"
	params:
		dbdir = config["rdir"] + "/kaiju_db"
	threads: config["kaiju_threads"]
	conda:
		config["wdir"] + "/envs/kaiju.yaml"
	log:
		config["rdir"] + "/logs/build_kaiju.log"
	shell:
		"""
		kaiju-mkbwt -n {threads} -a ACDEFGHIKLMNPQRSTVWY -o "{params.dbdir}/proteins" {input.faa}
		kaiju-mkfmi "{params.dbdir}/proteins"
		"""

