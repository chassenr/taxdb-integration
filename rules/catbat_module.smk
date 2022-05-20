rule copy_taxonomy_cat:
	input:
		names = config["rdir"] + "/DB_taxonomy/names.dmp"
	output:
		names = config["rdir"] + "/CAT_db/CAT_taxonomy/names.dmp"
	params:
		taxdir = config["rdir"] + "/DB_taxonomy",
		cat_dir = config["rdir"] + "/CAT_db",
		cat_tax = config["rdir"] + "/CAT_db/CAT_taxonomy"
	shell:
		"""
		mkdir -p {params.cat_dir}
		cp {params.taxdir}/*.dmp {params.cat_tax}/
		"""

rule format_accession_taxid_cat:
	input:
		protein_accmap = config["rdir"] + "/kaiju_db/library/proteins_kaiju_accession_map.txt"
	output:
		acc2taxid = config["rdir"] + "/CAT_db/CAT_taxonomy/prot.accession2taxid.FULL.gz"
	params:
		cat_tax = config["rdir"] + "/CAT_db/CAT_taxonomy"
	shell:
		"""
		echo -e 'accession.version\\ttaxid' > {params.cat_tax}/prot.accession2taxid.FULL
		paste <(cut -f3 {input.protein_accmap} | sed 's/$/.1/') <(cut -f1 {input.protein_accmap}) >> {params.cat_tax}/prot.accession2taxid.FULL
		gzip {params.cat_tax}/prot.accession2taxid.FULL
		"""

rule format_protein_cat:
	input:
		protein_formatted = config["rdir"] + "/kaiju_db/library/proteins.faa"
	output:
		faa_gz = config["rdir"] + "/CAT_db/CAT_database/nr.gz"
	shell:
		"""
		sed '/^>/s/$/.1 x/' {input.protein_formatted} | gzip > {output.faa_gz}
		"""

# don't work with clustered protein file, but rather with kaiju library
# diamond does not work with unusual line endings (mmseqs2 bug)

#rule format_protein_cat:
#	input:
#		protein_clustered = config["rdir"] + "/proteins_clustered/proteins_clustered.faa"
#	output:
#		faa_gz = config["rdir"] + "/CAT_db/CAT_database/nr.gz"
#	shell:
#		"""
#		sed -e '/^>/s/^>[0-9]* />/' -e '/^>/s/ *$/.1/' {input.protein_clustered} | tr 'BZ' 'DE' | sed '/^>/! s/[^ARNDCQEGHILKMFPSTWYV]//g' | gzip > {output.faa_gz}
#		"""

#rule format_accession_taxid_cat:
#	input:
#		acc2taxid = config["rdir"] + "/tax_combined/prot_accession2taxid.txt"
#	output:
#		acc2taxid = config["rdir"] + "/CAT_db/CAT_taxonomy/prot.accession2taxid.FULL.gz"
#	params:
#		cat_tax = config["rdir"] + "/CAT_db/CAT_taxonomy"
#	shell:
#		"""
#		echo -e 'accession.version\\ttaxid' > {params.cat_tax}/prot.accession2taxid.FULL
#		sed -e 's/\\t/_/' -e 's/\\t/.1\\t/' {input.acc2taxid} >> {params.cat_tax}/prot.accession2taxid.FULL
#		gzip {params.cat_tax}/prot.accession2taxid.FULL
#		"""

rule build_CAT_db:
	input:
		acc2taxid = config["rdir"] + "/CAT_db/CAT_taxonomy/prot.accession2taxid.FULL.gz",
		faa_gz = config["rdir"] + "/CAT_db/CAT_database/nr.gz",
		names = config["rdir"] + "/CAT_db/CAT_taxonomy/names.dmp"
	output:
		done = config["rdir"] + "/CAT_db/done"
	params:
		cat_tax = config["rdir"] + "/CAT_db/CAT_taxonomy",
		cat_db = config["rdir"] + "/CAT_db/CAT_database"
	conda: 
		config["wdir"] + "/envs/cat.yaml"
	threads: config["parallel_threads"]
	log:
		config["rdir"] + "/logs/build_CAT_db.log"
	shell:
		"""
		CAT prepare --existing -t {params.cat_tax} -d {params.cat_db} -n {threads} &>> {log}
		FILE=$(ls -1 {params.cat_db}/*.dmnd)
		if [[ -f $FILE ]]
		then
		  touch {output.done}
		fi
		"""

