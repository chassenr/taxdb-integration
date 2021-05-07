rule group_by_taxpath:
	input:
		protein_tax = config["rdir"] + "/tax_combined/protein_taxonomy_good.txt"
	output:
		path_list = config["rdir"] + "/proteins_clustered/taxon_list.txt"
	params:
		script = config["wdir"] + "/scripts/protein_group_by_taxpath.R"
	conda:
		config["wdir"] + "/envs/r.yaml"
	log:
		config["rdir"] + "/logs/protein_group_by_taxpath.log"
	shell:
		"""
		{params.script} -i {input.protein_tax} -o {output.path_list} &>> {log}
		"""

rule cluster_proteins:
	input:
		protein_tax = config["rdir"] + "/tax_combined/protein_taxonomy_good.txt",
		path_list = config["rdir"] + "/proteins_clustered/taxon_list.txt"
	output:
		protein_clustered = config["rdir"] + "/proteins_clustered/proteins_clustered.faa"
	params:
		script = config["wdir"] + "/scripts/cluster_proteins_mmseqs.sh",
		cov = config["cov"],
		cov_mode = config["cov_mode"],
		min_seq_id = config["min_seq_id"],
		pdir = config["rdir"] + "/proteins_all",
		outdir = config["rdir"] + "/proteins_clustered"
	conda:
		config["wdir"] + "/envs/mmseqs2.yaml"
	threads: config["parallel_threads"]
	shell:
		"""
		mkdir -p {params.outdir}/tmp
		cat {input.path_list} | parallel -j{threads} '{params.script} {{}} {input.protein_tax} {params.pdir} {params.outdir}/tmp {params.cov} {params.cov_mode} {params.min_seq_id}' >> {output.protein_clustered}
		"""

rule clean_up_post_clustering:
	input:
		protein_clustered = config["rdir"] + "/proteins_clustered/proteins_clustered.faa"
	output:
		done = config["rdir"] + "/proteins_clustered/done"
	params:
		tmpdir = config["rdir"] + "/proteins_clustered/tmp"
	threads: config["parallel_threads"]
	shell:
		"""
		find {params.tmpdir} -type d -name 'tmp_*' | xargs -n 1 -P {threads} rm -rf
		find {params.tmpdir} -type f -name '*_all_seqs.fasta' | xargs -n 1 -P {threads} rm
		find {params.tmpdir} -type f -name '*_rep_seq.fasta' | xargs -n 1 -P {threads} rm
		touch {output.done} 
		"""

