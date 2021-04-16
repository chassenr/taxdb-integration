rule get_faa_links:
	input:
		links_genomic = config["rdir"] + "/{library_name}/assembly_url_genomic.txt"
	output:
		links_faa = config["rdir"] + "/{library_name}/assembly_url_faa.txt"
	shell:
		"""
		sed 's/_genomic\.fna\.gz/_protein\.faa\.gz/' {input.links_genomic} > {output.links_faa}
		"""

