if config["kraken2_preset"] == "coarse":
	rule select_euk_coarse:
		input:
			tax = config["rdir"] + "/tax_combined/{library_name}_derep_taxonomy.txt",
			meta = config["rdir"] + "/{library_name}/metadata/genome_metadata.txt",
			gen2taxid = config["rdir"] + "/tax_combined/full_genome2taxid.txt"
		output:
			tax_select = config["rdir"] + "/" + config["db_name"] + "/{library_name}_kraken2_select_taxonomy.txt",
			kraken2_select = config["rdir"] + "/" + config["db_name"] + "/{library_name}_kraken2_select_accessions.txt"
		params:
			script = config["wdir"] + "/scripts/coarse_ncbi_selection.R",
			rank = lambda wildcards: config["rank_coarse"][wildcards.library_name],
			sublineage = lambda wildcards: config["coarse_sublineage"][wildcards.library_name],
			rank_sublineage = lambda wildcards: config["rank_sublineage"][wildcards.library_name],
			nmax = config["nmax_coarse"]
		conda:
			config["wdir"] + "/envs/r.yaml"
		log:
			config["rdir"] + "/logs/select_coarse_{library_name}.log"
		shell:
			"""
			{params.script} -t {input.tax} -m {input.meta} -r {params.rank} -s "{params.sublineage}" -R {params.rank_sublineage} -n {params.nmax} -o {output.tax_select} &>>{log}
			cut -f1 {output.tax_select} | grep -F -f - {input.gen2taxid} > {output.kraken2_select}
			"""

	rule collect_euk_coarse:
		input:
			kraken2_select = config["rdir"] + "/" + config["db_name"] + "/{library_name}_kraken2_select_accessions.txt"
		output:
			linked = config["rdir"] + "/" + config["db_name"] + "/genomes/{library_name}_done"
		params:
			gendir = config["rdir"] + "/derep_combined",
			outdir = config["rdir"] + "/" + config["db_name"] + "/genomes/"
		shell:
			"""
			mkdir -p {params.outdir}
			find {params.gendir} -name '*.gz' | grep -F -f <(cut -f1 {input.kraken2_select}) | while read line
			do
			  ln -sf "$line" {params.outdir}
			done
			touch {output.linked}
                        """

	rule collect_gtdb_coarse:
		input:
			gtdb_reps = config["rdir"] + "/gtdb/metadata/gtdb_reps_tax.txt",
			gen2taxid = config["rdir"] + "/tax_combined/full_genome2taxid.txt"
		output:
			kraken2_select = config["rdir"] + "/" + config["db_name"] + "/gtdb_kraken2_select_accessions.txt"
		params:
			gendir = config["rdir"] + "/gtdb/reps_genomes",
			outdir = config["rdir"] + "/" + config["db_name"] + "/genomes/"
		shell:
			"""
			mkdir -p {params.outdir}
			find {params.gendir} -type f -name '*.gz' | while read line
			do
			  ln -sf "$line" {params.outdir}
			done
			cut -f1 {input.gtdb_reps} | grep -F -f - {input.gen2taxid} > {output.kraken2_select}
			"""

	rule collect_checkv_coarse:
		input:
			reps_metadata = config["rdir"] + "/checkv/checkv_reps_metadata.txt",
			reps_fna = config["rdir"] + "/checkv/checkv_reps.fna",
			gen2taxid = config["rdir"] + "/tax_combined/full_genome2taxid.txt"
		output:
			kraken2_select = config["rdir"] + "/" + config["db_name"] + "/checkv_kraken2_select_accessions.txt",
			linked = config["rdir"] + "/" + config["db_name"] + "/genomes/checkv_done"
		params:
			gendir = config["rdir"] + "/checkv/reps_genomes",
			outdir = config["rdir"] + "/" + config["db_name"] + "/genomes/"
		conda:
			config["wdir"] + "/envs/parallel.yaml"
		threads: config["parallel_threads"]
		shell:
			"""
			mkdir -p {params.gendir}
			cd {params.gendir}
			cat {input.reps_fna} | awk '{{ if (substr($0, 1, 1)==">") {{filename=(substr($0,2) ".fa")}} print $0 > filename }}'
			find . -type f -name '*.fa' | parallel -j {threads} gzip {{}}
			cut -f1 {input.reps_metadata} | grep -F -f - {input.gen2taxid} > {output.kraken2_select}

			mkdir -p {params.outdir}
			find {params.gendir} -name '*.gz' | grep -F -f <(cut -f1 {output.kraken2_select}) | while read line
			do
			  ln -sf "$line" {params.outdir}
			done
			touch {output.linked}
			"""

	rule check_coarse:
		input:
			kraken2_select_checkv = config["rdir"] + "/" + config["db_name"] + "/checkv_kraken2_select_accessions.txt",
			kraken2_select_gtdb = config["rdir"] + "/" + config["db_name"] + "/gtdb_kraken2_select_accessions.txt",
			kraken2_select_euk = expand(config["rdir"] + "/" + config["db_name"] + "/{library_name}_kraken2_select_accessions.txt", library_name = LIBRARY_NAME),
			euk_linked = expand(config["rdir"] + "/" + config["db_name"] + "/genomes/{library_name}_done", library_name = LIBRARY_NAME)
		output:
			kraken2_select = config["rdir"] + "/" + config["db_name"] + "/kraken2_select_accessions.txt",
			checked = config["rdir"] + "/" + config["db_name"] + "/genomes/done"
		params:
			outdir = config["rdir"] + "/" + config["db_name"] + "/genomes/"
		shell:
			"""
			cat {input.kraken2_select_checkv} {input.kraken2_select_gtdb} {input.kraken2_select_euk} > {output.kraken2_select}
			if [[ $(cat {output.kraken2_select} | wc -l) == $(find {params.outdir} -name '*.gz' | wc -l) ]]
                        then
                          touch {output.checked}
                        fi
			"""

if config["kraken2_preset"] == "highres_pro":
	rule select_genomes_highres_pro:
		input:
			gtdb = config["rdir"] + "/tax_combined/gtdb_derep_taxonomy.txt",
			checkv = config["rdir"] + "/tax_combined/checkv_derep_taxonomy.txt",
			added_nuc_gtdb = config["rdir"] + "/tax_combined/pro_custom_post_derep_taxonomy.txt" if config["custom_gtdb_post_derep"] else [],
			added_nuc_checkv = config["rdir"] + "/tax_combined/vir_custom_post_derep_taxonomy.txt" if config["custom_checkv_post_derep"] else [],
			gen2taxid = config["rdir"] + "/tax_combined/full_genome2taxid.txt"
		output:
			kraken2_select = config["rdir"] + "/" + config["db_name"] + "/kraken2_select_accessions.txt"
		shell:
			"""
			cat {input.gtdb} {input.checkv} {input.added_nuc_gtdb} {input.added_nuc_checkv} | cut -f1 | grep -F -f - {input.gen2taxid} > {output.kraken2_select}
			"""

	rule collect_genomes_highres_pro:
		input:
			kraken2_select = config["rdir"] + "/" + config["db_name"] + "/kraken2_select_accessions.txt"
		output:
			linked = config["rdir"] + "/" + config["db_name"] + "/genomes/done"
		params:
			gendir = config["rdir"] + "/derep_combined",
			outdir = config["rdir"] + "/" + config["db_name"] + "/genomes/"
		shell:
			"""
			mkdir -p {params.outdir}
			find {params.gendir} -name '*.gz' | grep -F -f <(cut -f1 {input.kraken2_select}) | while read line
			do
			  ln -sf "$line" {params.outdir}
			done
			if [[ $(wc -l {input.kraken2_select}) == $(find {params.outdir} -type f -name '*.gz' | wc -l) ]]
			then
			  touch {output.linked}
			fi
			"""

if config["kraken2_preset"] == "user":
	rule collect_genomes:
		input:
			kraken2_select = config["kraken2_custom"]
		output:
			linked = config["rdir"] + "/" + config["db_name"] + "/genomes/done"
		params:
			gendir = config["rdir"] + "/derep_combined",
			outdir = config["rdir"] + "/" + config["db_name"] + "/genomes/"
		shell:
			"""
			mkdir -p {params.outdir}
			find {params.gendir} -name '*.gz' | grep -F -f <(cut -f1 {input.kraken2_select}) | while read line
			do
			  ln -sf "$line" {params.outdir}
			done
			if [[ $(wc -l {input.kraken2_select}) == $(find {params.outdir} -type f -name '*.gz' | wc -l) ]]
			then
			  touch {output.linked}
			fi
			"""

rule format_fasta_kraken2:
	input:
		kraken2_select = config["rdir"] + "/" + config["db_name"] + "/kraken2_select_accessions.txt",
		linked = config["rdir"] + "/" + config["db_name"] + "/genomes/done"
	output:
		library = config["rdir"] + "/" + config["db_name"] + "/tmp/library.fna"
	params:
		script = config["wdir"] + "/scripts/format_kraken2.sh",
		gendir = config["rdir"] + "/" + config["db_name"] + "/genomes"
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	threads:
		config["parallel_threads"]
	shell:
		"""
		cut -f1,3 {input.kraken2_select} | parallel -j{threads} '{params.script} {{}} {params.gendir}' | sed -e '/^>/!s/[a-z]/x/g' >> {output.library}
		"""

rule prelim_map:
	input:
		fasta = config["rdir"] + "/" + config["db_name"] + "/tmp/library.fna"
	output:
		map = config["rdir"] + "/" + config["db_name"] + "/tmp/prelim_map.txt"
	params:
		libdir = config["rdir"] + "/" + config["db_name"] + "/tmp"
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	shell:
		"""
		LC_ALL=C grep '^>' {input.fasta} | sed 's/^>//' > {params.libdir}/tmp.accnos
		NSEQ=$(wc -l {params.libdir}/tmp.accnos | cut -d' ' -f1)
		printf 'TAXID\\n%.0s' $(seq 1 $NSEQ) | paste - {params.libdir}/tmp.accnos | paste - <(cut -d'|' -f3 {params.libdir}/tmp.accnos) > {output.map}
		rm {params.libdir}/tmp.accnos
		"""

if config["kingdoms"]:
	rule detect_contamination:
		input:
			fasta = config["rdir"] + "/" + config["db_name"] + "/tmp/library.fna",
			map = config["rdir"] + "/" + config["db_name"] + "/tmp/prelim_map.txt",
			names = config["rdir"] + "/DB_taxonomy/names.dmp",
		output:
			kstring = config["rdir"] + "/" + config["db_name"] + "/decontamination/conterminator_string.txt",
			xstring = config["rdir"] + "/" + config["db_name"] + "/decontamination/conterminator_blacklist.txt",
			cmap = config["rdir"] + "/" + config["db_name"] + "/decontamination/cmap.txt",
			contam = config["rdir"] + "/" + config["db_name"] + "/decontamination/db_conterm_prediction"
		params:
			script = config["wdir"] + "/scripts/get_kingdoms_conterminator.R",
			tmpdir = config["rdir"] + "/" + config["db_name"] + "/decontamination/tmp",
			taxdir = config["rdir"] + "/DB_taxonomy/",
			prefix = config["rdir"] + "/" + config["db_name"] + "/decontamination/db",
			cmem = config["cmem"]
		conda:
			config["wdir"] + "/envs/r.yaml"
		threads: config["parallel_threads"]
		log:
			config["rdir"] + "/logs/run_conterminator.log"
		shell:
			"""
			# prepare fasta header mapping file for conterminator
			cut -f2,3 {input.map} > {output.cmap}
			# parse taxid string for conterminator kingdoms parameter
			{params.script} -t {params.taxdir} -k "prokaryotes,fungi,protists,plants,metazoa" -s "{params.taxdir}/accessionTaxa.sql" -o {output.kstring} -x {output.xstring}
			# run conterminator
			KSTR=$(cat {output.kstring})
			XSTR=$(cat {output.xstring})
			conterminator dna {input.fasta} {output.cmap} {params.prefix} {params.tmpdir} --mask-lower-case 1 --ncbi-tax-dump {params.taxdir} --threads {threads} --split-memory-limit {params.cmem} --blacklist $XSTR --kingdoms $KSTR &>> {log}
			"""

	rule filter_contamination:
		input:
			contam = config["rdir"] + "/" + config["db_name"] + "/decontamination/db_conterm_prediction",
			fasta = config["rdir"] + "/" + config["db_name"] + "/tmp/library.fna",
			map = config["rdir"] + "/" + config["db_name"] + "/tmp/prelim_map.txt"
		output:
			contam_filt = config["rdir"] + "/" + config["db_name"] + "/decontamination/db_conterm_prediction_filt",
			id_contam = config["rdir"] + "/" + config["db_name"] + "/decontamination/contam_id.accnos",
			fasta_noncontam = config["rdir"] + "/" + config["db_name"] + "/library/selection/library.fna",
			fasta_contam = config["rdir"] + "/" + config["db_name"] + "/decontamination/library_contam.fna",
			map_noncontam = config["rdir"] + "/" + config["db_name"] + "/library/selection/prelim_map.txt"
		conda:
			config["wdir"] + "/envs/bbmap.yaml"
		log:
			config["rdir"] + "/logs/contam_filter.log"
		shell:
			"""
			awk -v FS="\\t" -v OFS="\\t" '$5 >= 0 && $6 >= 0' {input.contam} > {output.contam_filt}
			cut -f2 {output.contam_filt} | sort | uniq > {output.id_contam}
			filterbyname.sh in={input.fasta} out={output.fasta_contam} names={output.id_contam} include=t ow=t &>> {log}
			filterbyname.sh in={input.fasta} out={output.fasta_noncontam} names={output.id_contam} include=f ow=t &>> {log}
			grep -v -F -f {output.id_contam} {input.map} > {output.map_noncontam}
			"""

	rule remove_contamination:
		input:
			contam = config["rdir"] + "/" + config["db_name"] + "/decontamination/db_conterm_prediction_filt",
			fasta_contam = config["rdir"] + "/" + config["db_name"] + "/decontamination/library_contam.fna",
			map_noncontam = config["rdir"] + "/" + config["db_name"] + "/library/selection/prelim_map.txt",
			fasta_noncontam = config["rdir"] + "/" + config["db_name"] + "/library/selection/library.fna"
		output:
			cleaned_fasta = config["rdir"] + "/" + config["db_name"] + "/decontamination/cleaned.fna",
			cleaned_map = config["rdir"] + "/" + config["db_name"] + "/decontamination/cleaned_map.txt"
		params:
			script = config["wdir"] + "/scripts/remove_contam_contigs.R",
			contam_dir = config["rdir"] + "/" + config["db_name"] + "/decontamination"
		conda:
			config["wdir"] + "/envs/r.yaml"
		log:
			config["rdir"] + "/logs/contam_cleaning.log"
		shell:
			"""
			{params.script} -i {input.fasta_contam} -c {input.contam} -o {output.cleaned_fasta} &>> {log}
			LC_ALL=C grep '^>' {output.cleaned_fasta} | sed 's/^>//' > "{params.contam_dir}/tmp.accnos"
			NSEQ=$(wc -l "{params.contam_dir}/tmp.accnos" | cut -d' ' -f1)
			printf 'TAXID\\n%.0s' $(seq 1 $NSEQ) | paste - "{params.contam_dir}/tmp.accnos" | paste - <(cut -d'|' -f3 "{params.contam_dir}/tmp.accnos") > {output.cleaned_map}
			rm "{params.contam_dir}/tmp.accnos"
			cat {output.cleaned_map} >> {input.map_noncontam}
			cat {output.cleaned_fasta} >> {input.fasta_noncontam}
			# to be implemented later: remove tmp folder in krakendb to save disk space
			"""

else:
	rule collect_db_files:
		input:
			fasta = config["rdir"] + "/" + config["db_name"] + "/tmp/library.fna",
			map = config["rdir"] + "/" + config["db_name"] + "/tmp/prelim_map.txt"
		output:
			fasta = config["rdir"] + "/" + config["db_name"] + "/library/selection/library.fna",
			map = config["rdir"] + "/" + config["db_name"] + "/library/selection/prelim_map.txt"
		shell:
			"""
			mv {input.fasta} {output.fasta}
			mv {input.map} {output.map}
			"""

rule copy_taxonomy:
	input:
		names = config["rdir"] + "/DB_taxonomy/names.dmp"
	output:
		names = config["rdir"] + "/" + config["db_name"] + "/taxonomy/names.dmp"
	params:
		taxdir = config["rdir"] + "/DB_taxonomy/",
		kraken_dir = config["rdir"] + "/" + config["db_name"],
		kraken_tax = config["rdir"] + "/" + config["db_name"] + "/taxonomy/"
	shell:
		"""
		mkdir -p {params.kraken_dir}
		cp -r {params.taxdir} {params.kraken_tax}
		"""

# to avoid ftp issue, recreate kraken2 code for adding UniVec files
# https://github.com/DerrickWood/kraken2/blob/561cc73fababe1dfd996e553e36ea1aff5642ef8/scripts/download_genomic_library.sh#L102-L117
if config["univec"]:
	rule add_univec:
		output:
			fasta = config["rdir"] + "/" + config["db_name"] + "/library/" + config["univec"] + "/library.fna",
			map = config["rdir"] + "/" + config["db_name"] + "/library/" + config["univec"] + "/prelim_map.txt"
		params:
			ncbi_server = config["ncbi_server"],
			uv_name = config["univec"],
			libdir = config["rdir"] + "/" + config["db_name"] + "/library/" + config["univec"]
		conda:
			config["wdir"] + "/envs/kraken2.yaml"
		log:
			config["rdir"] + "/logs/add_krakendb_univec.log"
		shell:
			"""
			wget -O "{params.libdir}/tmp.fna" "{params.ncbi_server}/pub/UniVec/{params.uv_name}"
			# choosing random artificial taxid (this taxid must not exist elsewhere in the database)
			sed -i 's/^>/>kraken:taxid|123456789|/' "{params.libdir}/tmp.fna"
			dustmasker -in "{params.libdir}/tmp.fna" -outfmt fasta | sed -e '/^>/!s/[a-z]/x/g' > {output.fasta}
			rm "{params.libdir}/tmp.fna"
			grep '^>' {output.fasta} | sed 's/^>//' > {params.libdir}/tmp.accnos
			NSEQ=$(wc -l {params.libdir}/tmp.accnos | cut -d' ' -f1)
			printf 'TAXID\\n%.0s' $(seq 1 $NSEQ) | paste - {params.libdir}/tmp.accnos | paste - <(cut -d'|' -f2 {params.libdir}/tmp.accnos) > {output.map}
			rm {params.libdir}/tmp.accnos
			"""

rule build_kraken_coarse:
	input:
		univec_fasta = config["rdir"] + "/" + config["db_name"] + "/library/" + config["univec"] + "/library.fna",
		univec_map = config["rdir"] + "/" + config["db_name"] + "/library/" + config["univec"] + "/prelim_map.txt",
		lib_fasta = config["rdir"] + "/" + config["db_name"] + "/library/selection/library.fna",
		libmap = config["rdir"] + "/" + config["db_name"] + "/library/selection/prelim_map.txt"
	output:
		hash = config["rdir"] + "/" + config["db_name"] + "/hash.k2d",
		opts = config["rdir"] + "/" + config["db_name"] + "/opts.k2d",
		map  = config["rdir"] + "/" + config["db_name"] + "/seqid2taxid.map",
		taxo = config["rdir"] + "/" + config["db_name"] + "/kraken2_db/taxo.k2d"
	params:
		dbdir = config["rdir"] + "/" + config["db_name"],
		kmer_len = config["kmer_len"],
		min_len = config["minimizer_len"],
		min_spaces = config["minimizer_spaces"],
		max_dbsize = config["max_dbsize"]
	threads: config["krakenbuild_threads"]
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	log:
		config["rdir"] + "/logs/build_kraken2_db.log"
	shell:
		"""
		kraken2-build --build --threads {threads} --db {params.dbdir} --kmer-len {params.kmer_len} --minimizer-len {params.min_len} --minimizer-spaces {params.min_spaces} --max-db-size {params.max_dbsize} &>> {log}
	"""



