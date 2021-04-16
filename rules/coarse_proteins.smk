rule get_faa_ncbi:
	input:
		tax_coarse = config["cdir"] + "/{library_name}/assembly_taxonomy_coarse.txt",
		url = config["rdir"] + "/{library_name}/assembly_url_faa.txt"
	output:
		coarse_ncbi_faa_download = config["cdir"] + "/{library_name}/proteins/done"
	params:
		outdir = config["cdir"] + "/{library_name}/proteins"
	conda:
		config["wdir"] + "/envs/download.yaml"
	threads: config["download_threads"]
	log:
		config["rdir"] + "/logs/download_faa_coarse_ncbi_{library_name}.log"
	shell:
		"""
		cut -f1 {input.tax_coarse} | sed 's/$/_/' | grep -F -f - {input.url} > "{params.outdir}/links"
		aria2c -i "{params.outdir}/links" -c -l "{params.outdir}/links.log" --dir {params.outdir} --max-tries=20 --retry-wait=5 --max-connection-per-server=1 --max-concurrent-downloads={threads} &>> {log}
		# We need to verify all files are there
		cat "{params.outdir}/links" | xargs -n 1 basename | sort > "{params.outdir}/tmp1"
		find {params.outdir} -type f -name '*.gz' | xargs -n 1 basename | sort > "{params.outdir}/tmp2"
		if diff "{params.outdir}/tmp1" "{params.outdir}/tmp2"
		then
		  touch "{params.outdir}/done"
		fi
		rm "{params.outdir}/links.log" "{params.outdir}/tmp1" "{params.outdir}/tmp2"
		"""

rule get_faa_gtdb_reps:
	input:
		gtdb_reps = config["cdir"] + "/gtdb/gtdb_reps_tax.txt"
	output:
		download_done = config["cdir"] + "/gtdb/proteins/done"
	params:
		outdir = config["cdir"] + "/gtdb",
		gtdb_link = config["gtdb_link"]
	threads: config["download_onefile"]
	conda:
		config["wdir"] + "/envs/download.yaml"
	log:
		config["rdir"] + "/logs/download_coarse_gtdb_faa.log"
	shell:
		"""
		aria2c -c -l "{params.outdir}/links.log" -d {params.outdir} --max-tries=20 --retry-wait=5 -x {threads} -j {threads} -s {threads} "{params.gtdb_link}/genomic_files_reps/gtdb_proteins_aa_reps.tar.gz" &>> {log}
                tar -C {params.outdir} -xzf "{params.outdir}/gtdb_proteins_aa_reps.tar.gz"
                mv {params.outdir}/gtdb_proteins_aa_reps_* "{params.outdir}/proteins/"
                rm "{params.outdir}/gtdb_proteins_aa_reps.tar.gz"
		touch {output.download_done}
                """

