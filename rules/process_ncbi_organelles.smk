rule download_refseq_plastid:
	output:
		plas_fna = config["rdir"] + "/organelle/plastid.fna.gz",
		plas_faa = config["rdir"] + "/organelle/plastid.faa.gz",
		plas_gbff = config["rdir"] + "/organelle/plastid.gbff.gz",
		plas_gpff = config["rdir"] + "/organelle/plastid.gpff.gz"
	params:
		library_dir = config["rdir"] + "/organelle",
		ncbi_server = config["ncbi_server"]
	shell:
		"""
		mkdir -p {params.library_dir}/plastid
		i=1
		while wget -q --method=HEAD {params.ncbi_server}/genomes/refseq/plastid/plastid.${{i}}.1.genomic.fna.gz; [ $? -eq 0 ]
		do
		  wget -O {params.library_dir}/plastid/plastid.${{i}}.1.genomic.fna.gz {params.ncbi_server}/genomes/refseq/plastid/plastid.${{i}}.1.genomic.fna.gz
		  wget -O {params.library_dir}/plastid/plastid.${{i}}.protein.faa.gz {params.ncbi_server}/genomes/refseq/plastid/plastid.${{i}}.protein.faa.gz
		  wget -O {params.library_dir}/plastid/plastid.${{i}}.genomic.gbff.gz {params.ncbi_server}/genomes/refseq/plastid/plastid.${{i}}.genomic.gbff.gz
		  wget -O {params.library_dir}/plastid/plastid.${{i}}.protein.gpff.gz {params.ncbi_server}/genomes/refseq/plastid/plastid.${{i}}.protein.gpff.gz
		  ((i++))
		done
		cat {params.library_dir}/plastid/plastid.*.1.genomic.fna.gz > {output.plas_fna}
		cat {params.library_dir}/plastid/plastid.*.protein.faa.gz > {output.plas_faa}
		cat {params.library_dir}/plastid/plastid.*.genomic.gbff.gz > {output.plas_gbff}
		cat {params.library_dir}/plastid/plastid.*.protein.gpff.gz > {output.plas_gpff}
		rm -rf {params.library_dir}/plastid/
		"""

rule download_refseq_mitochondria:
	output:
		mito_fna = config["rdir"] + "/organelle/mitochondria.fna.gz",
		mito_faa = config["rdir"] + "/organelle/mitochondria.faa.gz",
		mito_gbff = config["rdir"] + "/organelle/mitochondria.gbff.gz",
		mito_gpff = config["rdir"] + "/organelle/mitochondria.gpff.gz"
	params:
		library_dir = config["rdir"] + "/organelle",
		ncbi_server = config["ncbi_server"]
	shell:
		"""
		mkdir -p {params.library_dir}/mito
		i=1
		while wget -q --method=HEAD {params.ncbi_server}/genomes/refseq/mitochondrion/mitochondrion.${{i}}.1.genomic.fna.gz; [ $? -eq 0 ]
                do
                  wget -O {params.library_dir}/mito/mitochondrion.${{i}}.1.genomic.fna.gz {params.ncbi_server}/genomes/refseq/mitochondrion/mitochondrion.${{i}}.1.genomic.fna.gz
                  wget -O {params.library_dir}/mito/mitochondrion.${{i}}.protein.faa.gz {params.ncbi_server}/genomes/refseq/mitochondrion/mitochondrion.${{i}}.protein.faa.gz
		  wget -O {params.library_dir}/mito/mitochondrion.${{i}}.genomic.gbff.gz {params.ncbi_server}/genomes/refseq/mitochondrion/mitochondrion.${{i}}.genomic.gbff.gz
		  wget -O {params.library_dir}/mito/mitochondrion.${{i}}.protein.gpff.gz {params.ncbi_server}/genomes/refseq/mitochondrion/mitochondrion.${{i}}.protein.gpff.gz
                  ((i++))
                done
                cat {params.library_dir}/mito/mitochondrion.*.1.genomic.fna.gz > {output.mito_fna}
                cat {params.library_dir}/mito/mitochondrion.*.protein.faa.gz > {output.mito_faa}
		cat {params.library_dir}/mito/mitochondrion.*.genomic.gbff.gz > {output.mito_gbff}
		cat {params.library_dir}/mito/mitochondrion.*.protein.gpff.gz > {output.mito_gpff}
                rm -rf {params.library_dir}/mito/
		"""

rule taxid_refseq_plastid:
	input:
		gbff = config["rdir"] + "/organelle/plastid.gbff.gz",
		gpff = config["rdir"] + "/organelle/plastid.gpff.gz",
		nodes = config["rdir"] + "/ncbi_taxdump/nodes.dmp"
	output:
		names_nucl = config["rdir"] + "/organelle/plastid_nucl_names.txt",
		ranks_nucl = config["rdir"] + "/organelle/plastid_nucl_ranks.txt",
		names_prot = config["rdir"] + "/organelle/plastid_prot_names.txt",
		ranks_prot = config["rdir"] + "/organelle/plastid_prot_ranks.txt"
	params:
		library_dir = config["rdir"] + "/organelle",
		taxdir = config["rdir"] + "/ncbi_taxdump"
	conda:
		config["wdir"] + "/envs/taxonkit.yaml"
	shell:
		"""
		paste <(zgrep '^VERSION' {input.gbff} | sed 's/^VERSION     //') <(zgrep -A1 '^  ORGANISM  ' {input.gbff} | sed 's/^  ORGANISM  //' | grep -v '^            Eukaryota;' | tr '\\n' ' ' | sed 's/ -- /\\n/g' | sed 's/ $/\\n/' | sed -r 's/ +/ /g') > {output.names_nucl}
		paste <(zgrep '^VERSION' {input.gpff} | sed 's/^VERSION     //') <(zgrep -A1 '^  ORGANISM  ' {input.gpff} | sed 's/^  ORGANISM  //' | grep -v '^            Eukaryota;' | tr '\\n' ' ' | sed 's/ -- /\\n/g' | sed 's/ $/\\n/' | sed -r 's/ +/ /g') > {output.names_prot}
		cut -f2 {output.names_nucl} | sort | uniq | taxonkit name2taxid --data-dir {params.taxdir} | taxonkit lineage --taxid-field 2 -t -R --data-dir {params.taxdir} > {params.library_dir}/tmp_plas_nucl.txt
		cut -f1 {params.library_dir}/tmp_plas_nucl.txt | sort | uniq -d | while read line
		do
		  grep "$line.*$line" {params.library_dir}/tmp_plas_nucl.txt
		done > {output.ranks_nucl}
		cut -f1 {params.library_dir}/tmp_plas_nucl.txt | sort | uniq -d | grep -w -v -F -f - {params.library_dir}/tmp_plas_nucl.txt >> {output.ranks_nucl}
		cut -f2 {output.names_prot} | sort | uniq | taxonkit name2taxid --data-dir {params.taxdir} | taxonkit lineage --taxid-field 2 -t -R --data-dir {params.taxdir} > {params.library_dir}/tmp_plas_prot.txt
		cut -f1 {params.library_dir}/tmp_plas_prot.txt | sort | uniq -d | while read line
		do
		  grep "$line.*$line" {params.library_dir}/tmp_plas_prot.txt
		done  > {output.ranks_prot}
		cut -f1 {params.library_dir}/tmp_plas_prot.txt | sort | uniq -d | grep -w -v -F -f - {params.library_dir}/tmp_plas_prot.txt >> {output.ranks_prot}
		rm {params.library_dir}/tmp_plas_nucl.txt {params.library_dir}/tmp_plas_prot.txt
		"""

rule taxpath_refseq_plastid:
	input:
		names_nucl = config["rdir"] + "/organelle/plastid_nucl_names.txt",
		ranks_nucl = config["rdir"] + "/organelle/plastid_nucl_ranks.txt",
		names_prot = config["rdir"] + "/organelle/plastid_prot_names.txt",
		ranks_prot = config["rdir"] + "/organelle/plastid_prot_ranks.txt"
	output:
		tax_nucl = config["rdir"] + "/organelle/plastid_taxonomy_nucl.txt",
		tax_prot = config["rdir"] + "/organelle/plastid_taxonomy_prot.txt"
	params:
		script = config["wdir"] + "/scripts/format_organelle_taxpath.R"
	conda:
		config["wdir"] + "/envs/r.yaml"
	shell:
		"""
		{params.script} -n {input.names_nucl} -r {input.ranks_nucl} -p {input.ranks_prot} -t "plas" -o {output.tax_nucl} -P {output.tax_prot}
		"""

rule split_plastid_nucl:
	input:
		plas_fna = config["rdir"] + "/organelle/plastid.fna.gz"
	output:
		done = config["rdir"] + "/organelle/genomes/done_plas"
	params:
		outdir = config["rdir"] + "/organelle/genomes/"
	conda:
		config["wdir"] + "/envs/parallel.yaml"
	threads: config["parallel_threads"]
	shell:
		"""
		mkdir -p {params.outdir}
		cd {params.outdir}
		zcat {input} | awk -F' ' '{{if(/^>/){{sub("^>", "");name=$1;sub("^", ">");print > name"_plas.fa"}}else{{print > name"_plas.fa"}}}}'
		find . -type f -name '*_plas.fa' | parallel -j {threads} gzip {{}}
		touch {output.done}
		"""

rule split_plastid_prot:
	input:
		plas_faa = config["rdir"] + "/organelle/plastid.faa.gz",
		names_prot = config["rdir"] + "/organelle/plastid_prot_names.txt"
	output:
		done = config["rdir"] + "/organelle/proteins/done_plas"
	params:
		outdir = config["rdir"] + "/organelle/proteins"
	conda:
		config["wdir"] + "/envs/parallel.yaml"
	threads: config["parallel_threads"]
	shell:
		"""
		mkdir -p {params.outdir}
		cd {params.outdir}
		cut -f2 {input.names_prot} | sed "s/[^[:alnum:]]/_/g" | paste <(cut -f1 {input.names_prot}) - > tmp_plas_names.txt
		# need to linearize faa for this to work
		seqkit seq -w 0 -j {threads} {input.plas_faa} > tmp_plas_lin.faa
		cut -f2 tmp_plas_names.txt | sort | uniq | parallel -j {threads} 'grep -w -F "{{}}" tmp_plas_names.txt | cut -f1 | grep -A1 -F -f - tmp_plas_lin.faa | sed "/^--$/d" > {{}}_plas.faa'
		rm tmp_plas_names.txt tmp_plas_lin.faa
		find . -type f -name '*_plas.faa' | parallel -j {threads} gzip {{}}
		touch {output.done}
		"""

rule taxid_refseq_mitochondria:
	input:
		gbff = config["rdir"] + "/organelle/mitochondria.gbff.gz",
		gpff = config["rdir"] + "/organelle/mitochondria.gpff.gz",
		nodes = config["rdir"] + "/ncbi_taxdump/nodes.dmp"
	output:
		names_nucl = config["rdir"] + "/organelle/mitochondria_nucl_names.txt",
		ranks_nucl = config["rdir"] + "/organelle/mitochondria_nucl_ranks.txt",
		names_prot = config["rdir"] + "/organelle/mitochondria_prot_names.txt",
		ranks_prot = config["rdir"] + "/organelle/mitochondria_prot_ranks.txt"
	params:
		library_dir = config["rdir"] + "/organelle",
		taxdir = config["rdir"] + "/ncbi_taxdump"
	conda:
		config["wdir"] + "/envs/taxonkit.yaml"
	shell:
		"""
		paste <(zgrep '^VERSION' {input.gbff} | sed 's/^VERSION     //') <(zgrep -A1 '^  ORGANISM  ' {input.gbff} | sed 's/^  ORGANISM  //' | grep -v '^            Eukaryota;' | tr '\\n' ' ' | sed 's/ -- /\\n/g' | sed 's/ $/\\n/' | sed -r 's/ +/ /g') > {output.names_nucl}
		paste <(zgrep '^VERSION' {input.gpff} | sed 's/^VERSION     //') <(zgrep -A1 '^  ORGANISM  ' {input.gpff} | sed 's/^  ORGANISM  //' | grep -v '^            Eukaryota;' | tr '\\n' ' ' | sed 's/ -- /\\n/g' | sed 's/ $/\\n/' | sed -r 's/ +/ /g') > {output.names_prot}
		cut -f2 {output.names_nucl} | sort | uniq | taxonkit name2taxid --data-dir {params.taxdir} | taxonkit lineage --taxid-field 2 -t -R --data-dir {params.taxdir} > {params.library_dir}/tmp_mito_nucl.txt
		cut -f1 {params.library_dir}/tmp_mito_nucl.txt | sort | uniq -d | while read line
		do
		  grep "$line.*$line" {params.library_dir}/tmp_mito_nucl.txt
		done > {output.ranks_nucl}
		cut -f1 {params.library_dir}/tmp_mito_nucl.txt | sort | uniq -d | grep -w -v -F -f - {params.library_dir}/tmp_mito_nucl.txt >> {output.ranks_nucl}
		cut -f2 {output.names_prot} | sort | uniq | taxonkit name2taxid --data-dir {params.taxdir} | taxonkit lineage --taxid-field 2 -t -R --data-dir {params.taxdir} > {params.library_dir}/tmp_mito_prot.txt
		cut -f1 {params.library_dir}/tmp_mito_prot.txt | sort | uniq -d | while read line
		do
		  grep "$line.*$line" {params.library_dir}/tmp_mito_prot.txt
		done  > {output.ranks_prot}
		cut -f1 {params.library_dir}/tmp_mito_prot.txt | sort | uniq -d | grep -w -v -F -f - {params.library_dir}/tmp_mito_prot.txt >> {output.ranks_prot}
		rm {params.library_dir}/tmp_mito_nucl.txt {params.library_dir}/tmp_mito_prot.txt
		"""

rule taxpath_refseq_mitochondria:
	input:
		names_nucl = config["rdir"] + "/organelle/mitochondria_nucl_names.txt",
		ranks_nucl = config["rdir"] + "/organelle/mitochondria_nucl_ranks.txt",
		names_prot = config["rdir"] + "/organelle/mitochondria_prot_names.txt",
		ranks_prot = config["rdir"] + "/organelle/mitochondria_prot_ranks.txt"
	output:
		tax_nucl = config["rdir"] + "/organelle/mitochondria_taxonomy_nucl.txt",
		tax_prot = config["rdir"] + "/organelle/mitochondria_taxonomy_prot.txt"
	params:
		script = config["wdir"] + "/scripts/format_organelle_taxpath.R"
	conda:
		config["wdir"] + "/envs/r.yaml"
	shell:
		"""
		{params.script} -n {input.names_nucl} -r {input.ranks_nucl} -p {input.ranks_prot} -t "mito" -o {output.tax_nucl} -P {output.tax_prot}
		"""

rule split_mitochondria_nucl:
	input:
		mito_fna = config["rdir"] + "/organelle/mitochondria.fna.gz"
	output:
		done = config["rdir"] + "/organelle/genomes/done_mito"
	params:
		outdir = config["rdir"] + "/organelle/genomes/"
	conda:
		config["wdir"] + "/envs/parallel.yaml"
	threads: config["parallel_threads"]
	shell:
		"""
		mkdir -p {params.outdir}
		cd {params.outdir}
		zcat {input} | awk -F' ' '{{if(/^>/){{sub("^>", "");name=$1;sub("^", ">");print > name"_mito.fa"}}else{{print > name"_mito.fa"}}}}'
		find . -type f -name '*_mito.fa' | parallel -j {threads} gzip {{}}
		touch {output.done}
		"""

rule split_mitochondria_prot:
	input:
		mito_faa = config["rdir"] + "/organelle/mitochondria.faa.gz",
		names_prot = config["rdir"] + "/organelle/mitochondria_prot_names.txt"
	output:
		done = config["rdir"] + "/organelle/proteins/done_mito"
	params:
		outdir = config["rdir"] + "/organelle/proteins"
	conda:
		config["wdir"] + "/envs/parallel.yaml"
	threads: config["parallel_threads"]
	shell:
		"""
		mkdir -p {params.outdir}
		cd {params.outdir}
		cut -f2 {input.names_prot} | sed "s/[^[:alnum:]]/_/g" | paste <(cut -f1 {input.names_prot}) - > tmp_mito_names.txt
		# need to linearize faa for this to work
		seqkit seq -w 0 -j {threads} {input.mito_faa} > tmp_mito_lin.faa
		cut -f2 tmp_mito_names.txt | sort | uniq | parallel -j {threads} 'grep -w -F "{{}}" tmp_mito_names.txt | cut -f1 | grep -A1 -F -f - tmp_mito_lin.faa | sed "/^--$/d" > {{}}_mito.faa'
		rm tmp_mito_names.txt tmp_mito_lin.faa
		find . -type f -name '*_mito.faa' | parallel -j {threads} gzip {{}}
		touch {output.done}
		"""

rule collect_refseq_organelles:
	input:
		split_mito_nucl = config["rdir"] + "/organelle/genomes/done_mito",
		split_mito_prot = config["rdir"] + "/organelle/proteins/done_mito",
		split_plas_nucl = config["rdir"] + "/organelle/genomes/done_plas",
		split_plas_prot = config["rdir"] + "/organelle/proteins/done_plas",
		tax_mito_nucl = config["rdir"] + "/organelle/mitochondria_taxonomy_nucl.txt",
		tax_mito_prot = config["rdir"] + "/organelle/mitochondria_taxonomy_prot.txt",
		tax_plas_nucl = config["rdir"] + "/organelle/plastid_taxonomy_nucl.txt",
		tax_plas_prot = config["rdir"] + "/organelle/plastid_taxonomy_prot.txt"
	output:
		tax = config["rdir"] + "/tax_combined/organelle_taxonomy.txt"
	params:
		library_dir = config["rdir"] + "/organelle",
		gdir = config["rdir"] + "/derep_combined/",
		pdir = config["rdir"] + "/proteins_all/"
	shell:
		"""
		mkdir -p {params.gdir}
		mkdir -p {params.pdir}
		cut -f1 {input.tax_mito_nucl} | while read line
		do
		  ln -sf {params.library_dir}/genomes/$line.fa.gz {params.gdir}
		done
		cut -f1 {input.tax_plas_nucl} | while read line
		do
		  ln -sf {params.library_dir}/genomes/$line.fa.gz {params.gdir}
		done
		cut -f1 {input.tax_mito_prot} | while read line
		do
		  ln -sf {params.library_dir}/proteins/$line.faa.gz {params.pdir}
		done
		cut -f1 {input.tax_plas_prot} | while read line
		do
		  ln -sf {params.library_dir}/proteins/$line.faa.gz {params.pdir}
		done
		cat {input.tax_mito_nucl} {input.tax_mito_prot} {input.tax_plas_nucl} {input.tax_plas_prot} > {output.tax}
		"""
		
