# the fasta reps from GTDB include more than just the actual representatives
# take those faa preferentially
rule download_gtdb_faa_reps:
	output:
		gtdb_faa_done = config["rdir"] + "/gtdb/reps_proteins/done"
	params:
		outdir = config["rdir"] + "/gtdb/reps_proteins",
		gtdb_link = config["gtdb_link"]
	threads: config["download_onefile"]
	conda:
		config["wdir"] + "/envs/download.yaml"
	log:
		config["rdir"] + "/logs/download_faa_reps_gtdb.log"
	shell:
		"""
		mkdir -p {params.outdir}
		aria2c -c -l "{params.outdir}/faa_reps_links.log" -d {params.outdir}/../ --max-tries=20 --retry-wait=5 -x {threads} -j {threads} -s {threads} "{params.gtdb_link}/genomic_files_reps/gtdb_proteins_aa_reps.tar.gz" &>> {log}
		tar -xzf "{params.outdir}/../gtdb_proteins_aa_reps.tar.gz" -C {params.outdir} --strip-components 2
		cd {params.outdir}
		find ./ -type f -name "*.faa" | xargs -n 1 basename | while read line
		do
		  newfile=$(echo $line | sed 's/^[RG][SB]_//')
		  mv $line $newfile
		done
		rm "{params.outdir}/../gtdb_proteins_aa_reps.tar.gz"
		touch {output.gtdb_faa_done}
		"""

rule gzip_gtdb_faa_reps:
	input:
		gtdb_faa_done = config["rdir"] + "/gtdb/reps_proteins/done"
	output:
		gtdb_faa_zipped = config["rdir"] + "/gtdb/reps_proteins/zipped"
	params:
		outdir = config["rdir"] + "/gtdb/reps_proteins"
	conda:
		config["wdir"] + "/envs/parallel.yaml"
	threads: config["parallel_threads"]
	shell:
		"""
		find {params.outdir} -type f -name '*.faa' | parallel -j{threads} gzip {{}}
		touch {output.gtdb_faa_zipped}
		"""

# link those faa which are available via gtdb reps
rule gtdb_add_protein_reps:
	input:
		gtdb_faa_zipped = config["rdir"] + "/gtdb/reps_proteins/zipped",
		tax = config["rdir"] + "/tax_combined/gtdb_derep_taxonomy.txt"
	output:
		faa_reps = config["rdir"] + "/gtdb/metadata/gtdb_faa_from_reps.accnos"
	params:
		outdir = config["rdir"] + "/gtdb/proteins",
		repdir = config["rdir"] + "/gtdb/reps_proteins"
	shell:
		"""
		mkdir -p {params.outdir}
		find {params.repdir} -type f -name '*.gz' | xargs -n 1 basename | cut -d'_' -f1,2 | grep -F -f <(cut -f1 {input.tax}) > {output.faa_reps}
		sed 's/$/_protein\.faa\.gz/' {output.faa_reps} | while read line
		do
		  ln -sf "{params.repdir}/$line" {params.outdir}
		done
		"""

# download the missing from ncbi
rule download_proteins_gtdb:
	input:
		faa_reps = config["rdir"] + "/gtdb/metadata/gtdb_faa_from_reps.accnos",
		gtdb_links = config["rdir"] + "/gtdb/metadata/gtdb_download_info.txt",
		tax = config["rdir"] + "/tax_combined/gtdb_derep_taxonomy.txt"
	output:
		faa_ncbi = config["rdir"] + "/gtdb/metadata/gtdb_faa_from_ncbi.accnos"
	params:
		outdir = config["rdir"] + "/gtdb/proteins"
	conda:
		config["wdir"] + "/envs/download.yaml"
	threads: config["download_threads"]
	log:
		config["rdir"] + "/logs/download_proteins_gtdb.log"
	shell:
		"""
		cut -f1 {input.tax} | grep -v -F -f {input.faa_reps} | grep -F -f - {input.gtdb_links} | awk -v FS="\\t" -v OFS="\\t" '$11 == "TRUE"' | cut -f6 > {output.faa_ncbi}
		grep -F -f {output.faa_ncbi} {input.gtdb_links} | cut -f8,10 | awk -v FS="\\t" -v OFS="\\t" '{{print $1"\\n out="$2}}' > "{params.outdir}/links"
		aria2c -i "{params.outdir}/links" -c -l "{params.outdir}/links.log" --dir {params.outdir} --max-tries=20 --retry-wait=5 --max-connection-per-server=1 --max-concurrent-downloads={threads} &>> {log}
		rm "{params.outdir}/links" "{params.outdir}/links.log"
		"""

# if no faa gtdb reps or download, predict with prodigal
rule predict_faa_gtdb:
	input:
		faa_reps = config["rdir"] + "/gtdb/metadata/gtdb_faa_from_reps.accnos",
		faa_ncbi = config["rdir"] + "/gtdb/metadata/gtdb_faa_from_ncbi.accnos",
		tax = config["rdir"] + "/tax_combined/gtdb_derep_taxonomy.txt",
		gtdb_links = config["rdir"] + "/gtdb/metadata/gtdb_download_info.txt"
	output:
		tax_prot = config["rdir"] + "/gtdb/metadata/gtdb_protein_taxonomy.txt"
	params:
		outdir = config["rdir"] + "/gtdb/proteins",
		cdsdir = config["rdir"] + "/gtdb/cds_proteins",
		gendir = config["rdir"] + "/gtdb/genomes",
		tmpdir = config["rdir"] + "/gtdb/tmp"
	conda:
		config["wdir"] + "/envs/prodigal.yaml"
	threads: config["parallel_threads"]
	log:
		config["rdir"] + "/logs/predict_proteins_gtdb.log"
	shell:
		"""
		mkdir -p {params.tmpdir}
		mkdir -p {params.cdsdir}
		cat {input.faa_ncbi} {input.faa_reps} | grep -v -F -f - {input.tax} | cut -f1 | grep -F -f - {input.gtdb_links} | cut -f4 | sed 's/_genomic\.fna\.gz//' | parallel -j{threads} 'zcat {params.gendir}/{{}}_genomic.fna.gz > {params.tmpdir}/{{}}_genomic.fna'
		cd {params.tmpdir}
		ls -1 | sed 's/_genomic\.fna//' | parallel -j{threads} 'prodigal -i {{}}_genomic.fna -o {params.cdsdir}/{{}}.gff -a {params.outdir}/{{}}_protein.faa -f gff'
		cd
		find {params.outdir} -type f -name '*.faa' | parallel -j{threads} 'gzip {{}}'
		if [[ $(cat {input.tax} | wc -l ) == $(find {params.outdir} -name '*.gz' | wc -l ) ]]
		then
		  rm -rf {params.tmpdir}
		  cp {input.tax} {output.tax_prot}
		fi
		"""

rule custom_gtdb_proteins_pre_derep:
	input:
		tax = config["rdir"] + "/tax_combined/gtdb_derep_taxonomy.txt",
		tax_prot = config["rdir"] + "/gtdb/metadata/gtdb_protein_taxonomy.txt"
	output:
		tax_prot_added = config["rdir"] + "/gtdb/metadata/gtdb_protein_taxonomy_added.txt"
	params:
		add = config["custom_gtdb_pre_derep"],
		outdir = config["rdir"] + "/gtdb/proteins/"
	shell:
		"""
		if [[ "{params.add}" != "" ]]
		then
		  awk -v FS="\\t" -v OFS="\\t" '$4 != "NA"' {params.add} | grep -F -f <(cut -f1 {input.tax}) > "{params.outdir}/tmp"
		  cut -f4 "{params.outdir}/tmp" | while read line
		  do
		    ln -sf "$line" {params.outdir}
		  done
		  cat {input.tax_prot} > {output.tax_prot_added}
		  cut -f1,2 "{params.outdir}/tmp" >> {output.tax_prot_added}
		else
		  cp {input.tax_prot} {output.tax_prot_added}
		fi
		"""

rule collect_gtdb_proteins:
	input:
		tax_prot = config["rdir"] + "/gtdb/metadata/gtdb_protein_taxonomy_added.txt"
	output:
		tax_prot = config["rdir"] + "/tax_combined/gtdb_protein_taxonomy.txt"
	params:
		pdir = config["rdir"] + "/gtdb/proteins",
		outdir = config["rdir"] + "/proteins_all/"
	shell:
		"""
		mkdir -p {params.outdir}
		find {params.pdir} -name '*.gz' | while read line
		do
		  ln -sf "$line" {params.outdir}
		done
		cp {input.tax_prot} {output.tax_prot}
		"""

if config["custom_pro_prot"]:
	rule add_custom_pro_prot:
		output:
			prot_added = config["rdir"] + "/tax_combined/pro_custom_protein_taxonomy.txt"
		params:
			add = config["custom_pro_prot"],
			outdir = config["rdir"] + "/proteins_all/"
		shell:
			"""
			mkdir -p {params.outdir}
			cut -f3 {params.add} | while read line
			do
			  ln -sf "$line" {params.outdir}
			done
			cut -f1,2 {params.add} > {output.prot_added}
			"""

