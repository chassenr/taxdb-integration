




rule masking_ncbi:
	input:
		file_list = config["rdir"] + "/kraken2_genomes/file_names_derep_genomes.txt",
		ncbi = config["rdir"] + "/tax_combined/{library_highres}_derep_taxonomy.txt",
		nodes = config["rdir"] + "/kraken2_db_euk/taxonomy/nodes.dmp",
		names = config["rdir"] + "/kraken2_db_euk/taxonomy/names.dmp"
	output:
		fasta = config["rdir"] + "/kraken2_db_euk/library/{library_highres}/library.fna"
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	threads: config["masking_threads"]
	shell:
		"""
		cut -f1 {input.ncbi} | grep -F -f - {input.file_list} | parallel -j{threads} 'dustmasker -in {{}} -outfmt fasta' | sed -e '/^>/!s/[a-z]/x/g' >> {output.fasta}
		"""

rule prelim_map_ncbi:
	input:
		fasta = config["rdir"] + "/kraken2_db_euk/library/{library_highres}/library.fna"
	output:
		map = config["rdir"] + "/kraken2_db_euk/library/{library_highres}/prelim_map.txt"
	params:
		libdir = config["rdir"] + "/kraken2_db_euk/library/{library_highres}",
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

rule masking_gtdb:
        input:
                file_list = config["rdir"] + "/kraken2_genomes/file_names_derep_genomes.txt",
                gtdb = config["rdir"] + "/tax_combined/gtdb_derep_taxonomy.txt",
                nodes = config["rdir"] + "/kraken2_taxonomy/nodes.dmp",
                names = config["rdir"] + "/kraken2_taxonomy/names.dmp"
        output:
                fasta = config["rdir"] + "/kraken2_db_pro/library/gtdb/library.fna"
        conda:
                config["wdir"] + "/envs/kraken2.yaml"
        threads: config["masking_threads"]
        shell:
                """
                cut -f1 {input.gtdb} | grep -F -f - {input.file_list} | parallel -j{threads} 'dustmasker -in {{}} -outfmt fasta' | sed -e '/^>/!s/[a-z]/x/g' >> {output.fasta}
                """

rule prelim_map_gtdb:
        input:
                fasta = config["rdir"] + "/kraken2_db_pro/library/gtdb/library.fna"
        output:
                map = config["rdir"] + "/kraken2_db_pro/library/gtdb/prelim_map.txt"
        params:
                libdir = config["rdir"] + "/kraken2_db_pro/library/gtdb"
        conda:
                config["wdir"] + "/envs/kraken2.yaml"
        shell:
                """
                LC_ALL=C grep '^>' {input.fasta} | sed 's/^>//' > {params.libdir}/tmp_gtdb.accnos
                NSEQ=$(wc -l {params.libdir}/tmp_gtdb.accnos | cut -d' ' -f1)
                printf 'TAXID\\n%.0s' $(seq 1 $NSEQ) | paste - {params.libdir}/tmp_gtdb.accnos | paste - <(cut -d'|' -f3 {params.libdir}/tmp_gtdb.accnos) > {output.map}
                rm {params.libdir}/tmp_gtdb.accnos
                """


rule masking_checkv:
        input:
                file_list = config["rdir"] + "/kraken2_genomes/file_names_derep_genomes.txt",
                checkv = config["rdir"] + "/tax_combined/checkv_derep_taxonomy.txt",
                nodes = config["rdir"] + "/kraken2_db_pro/taxonomy/nodes.dmp",
                names = config["rdir"] + "/kraken2_db_pro/taxonomy/names.dmp"
        output:
                fasta = config["rdir"] + "/kraken2_db_pro/library/checkv/library.fna"
        conda:
                config["wdir"] + "/envs/kraken2.yaml"
        threads: config["masking_threads"]
        shell:
                """
                cut -f1 {input.checkv} | grep -F -f - {input.file_list} | parallel -j{threads} 'dustmasker -in {{}} -outfmt fasta' | sed -e '/^>/!s/[a-z]/x/g' >> {output.fasta}
                """

rule prelim_map_checkv:
        input:
                fasta = config["rdir"] + "/kraken2_db_pro/library/checkv/library.fna"
        output:
                map = config["rdir"] + "/kraken2_db_pro/library/checkv/prelim_map.txt"
        params:
                libdir = config["rdir"] + "/kraken2_db_pro/library/checkv"
        conda:
                config["wdir"] + "/envs/kraken2.yaml"
        shell:
                """
                LC_ALL=C grep '^>' {input.fasta} | sed 's/^>//' > {params.libdir}/tmp.accnos
                NSEQ=$(wc -l {params.libdir}/tmp.accnos | cut -d' ' -f1)
                printf 'TAXID\\n%.0s' $(seq 1 $NSEQ) | paste - {params.libdir}/tmp.accnos | paste - <(cut -d'|' -f3 {params.libdir}/tmp.accnos) > {output.map}
                rm {params.libdir}/tmp.accnos
                """


rule format_taxonomy:
	input:
		tax_combined = config["rdir"] + "/tax_combined/derep_taxonomy_combined.txt"
	output:
		tax_good = config["rdir"] + "/tax_combined/derep_taxonomy_good.txt",
		nodes = config["rdir"] + "/kraken2_taxonomy/nodes.dmp",
		names = config["rdir"] + "/kraken2_taxonomy/names.dmp",
		file_list = config["rdir"] + "/kraken2_genomes/file_names_derep_genomes.txt"
	params:
		krakendir = config["rdir"] + "/kraken2_genomes/genome_files",
		genomedir = config["rdir"] + "/derep_combined",
		tax_script = config["tax_script"]
	conda:
		config["wdir"] + "/envs/biopython.yaml"
	log:
		config["rdir"] + "/logs/format_taxonomy.log"
	shell:
		"""
		cut -f1,2 {input.tax_combined} > {output.tax_good}
		{params.tax_script} --gtdb {output.tax_good} --assemblies {params.genomedir} --nodes {output.nodes} --names {output.names} --kraken_dir {params.krakendir} &>> {log}
		# replace 'domain' with 'superkingdom (required for kaiju) to ensure that the same nodes.dmp and names.dmp files can be used for both databases
		sed -i -e 's/domain/superkingdom/g' {output.nodes}
		find {params.krakendir} -type f -name '*.fa' > {output.file_list}
		"""
# depending on the number and size of the genomes, it may be required to delete the zipped fna in the derep directory

rule link_taxonomy_pro:
	input:
		nodes = config["rdir"] + "/kraken2_taxonomy/nodes.dmp",
		names = config["rdir"] + "/kraken2_taxonomy/names.dmp"
	output:
		nodes_pro = config["rdir"] + "/kraken2_db_pro/taxonomy/nodes.dmp",
		names_pro = config["rdir"] + "/kraken2_db_pro/taxonomy/names.dmp"
	shell:
		"""
		ln -sf {input.nodes} {output.nodes_pro}
		ln -sf {input.names} {output.names_pro}
		"""

rule link_taxonomy_euk:
	input:
		nodes = config["rdir"] + "/kraken2_taxonomy/nodes.dmp",
		names = config["rdir"] + "/kraken2_taxonomy/names.dmp"
	output:
		nodes_euk = config["rdir"] + "/kraken2_db_euk/taxonomy/nodes.dmp",
		names_euk = config["rdir"] + "/kraken2_db_euk/taxonomy/names.dmp"
	shell:
		"""
		ln -sf {input.nodes} {output.nodes_euk}
		ln -sf {input.names} {output.names_euk}
		"""

rule build_krakendb_pro:
	input:
		gtdb_map = config["rdir"] + "/kraken2_db_pro/library/gtdb/prelim_map.txt",
		gtdb_fasta = config["rdir"] + "/kraken2_db_pro/library/gtdb/library.fna",
		checkv_map = config["rdir"] + "/kraken2_db_pro/library/checkv/prelim_map.txt",
		checkv_fasta = config["rdir"] + "/kraken2_db_pro/library/checkv/library.fna"
	output:
		hash = config["rdir"] + "/kraken2_db_pro/hash.k2d",
		opts = config["rdir"] + "/kraken2_db_pro/opts.k2d",
		map  = config["rdir"] + "/kraken2_db_pro/seqid2taxid.map",
		taxo = config["rdir"] + "/kraken2_db_pro/taxo.k2d"
	params:
		dbdir = config["rdir"] + "/kraken2_db_pro",
		kmer_len = config["kmer_len"],
		min_len = config["minimizer_len"],
		min_spaces = config["minimizer_spaces"],
		max_dbsize = config["max_dbsize"]
	threads: config["krakenbuild_threads"]
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	log:
		config["rdir"] + "/logs/build_kraken.log"
	shell:
		"""
		kraken2-build --build --threads {threads} --db {params.dbdir} --kmer-len {params.kmer_len} --minimizer-len {params.min_len} --minimizer-spaces {params.min_spaces} --max-db-size {params.max_dbsize} &>> {log}
		"""

rule build_krakendb_euk:
	input:
		ncbi_map = expand(config["rdir"] + "/kraken2_db_euk/library/{library_highres}/prelim_map.txt", library_highres = LIBRARY_HIGHRES),
		ncbi_fasta = expand(config["rdir"] + "/kraken2_db_euk/library/{library_highres}/library.fna", library_highres = LIBRARY_HIGHRES),
	output:
		hash = config["rdir"] + "/kraken2_db_euk/hash.k2d",
		opts = config["rdir"] + "/kraken2_db_euk/opts.k2d",
		map  = config["rdir"] + "/kraken2_db_euk/seqid2taxid.map",
		taxo = config["rdir"] + "/kraken2_db_euk/taxo.k2d"
	params:
		dbdir = config["rdir"] + "/kraken2_db_euk",
		kmer_len = config["kmer_len"],
		min_len = config["minimizer_len"],
		min_spaces = config["minimizer_spaces"],
		max_dbsize = config["max_dbsize"]
	threads: config["krakenbuild_threads"]
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	log:
		config["rdir"] + "/logs/build_kraken.log"
	shell:
		"""
		kraken2-build --build --threads {threads} --db {params.dbdir} --kmer-len {params.kmer_len} --minimizer-len {params.min_len} --minimizer-spaces {params.min_spaces} --max-db-size {params.max_dbsize} &>> {log}
		"""

