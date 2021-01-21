import glob
import pandas as pd

"""
Author: Christiane Hassenrueck, Antonio Fernandez-Guerra
Acknowledgements: Chiara Vanni
Affiliation: MARUM - Center for Marine Environmental Sciences University of Bremen
Aim: Create high resolution taxonomic database for all domains of life using GTDB and NCBI taxonomy
Run: snakemake -s Snakefile
"""

# from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
# min_version("5.1.2")

configfile: "config/config.yaml"

LIBRARY_NAME = config["library_name"]
LIBRARY_HIGHRES = config["library_highres"]

wildcard_constraints:
	library_name = '|'.join(LIBRARY_NAME),
	library_highres = '|'.join(LIBRARY_HIGHRES)

rule all:
	input:
		#genome_assemblies = expand(config["rdir"] + "/{library_name}/assembly_summary_combined.txt", library_name = LIBRARY_NAME),
		#download_complete_gtdb = config["rdir"] + "/gtdb/genomes/done",
		#download_complete_ncbi = expand(config["rdir"] + "/{library_name}/genomes/done", library_name = LIBRARY_NAME),
		#download_feature_counts = expand(config["rdir"] + "/{library_name}/feature_counts/done", library_name = LIBRARY_NAME),
		#download_contig_stats = expand(config["rdir"] + "/{library_name}/contig_stats/done", library_name = LIBRARY_NAME),
		#metadata = expand(config["rdir"] + "/{library_name}/metadata/genome_metadata.txt", library_name = LIBRARY_NAME),
		#gtdb_derep_meta = config["rdir"] + "/gtdb/metadata/gtdb_derep_taxonomy_meta.txt",
		checkv_taxonomy = config["rdir"] + "/checkv/checkv_taxonomy.txt",
		checkv_derep_meta = config["rdir"] + "/checkv/checkv_derep_taxonomy_meta.txt",
		#derep_meta = expand(config["rdir"] + "/{library_name}/derep_taxonomy_meta.txt", library_name = LIBRARY_NAME),
		#gtdb_derep_tax = config["rdir"] + "/tax_combined/gtdb_derep_taxonomy.txt",
		checkv_derep_tax = config["rdir"] + "/tax_combined/checkv_derep_taxonomy.txt",
		#ncbi_derep_tax = config["rdir"] + "/tax_combined/ncbi_derep_taxonomy.txt",
		#nodes = config["rdir"] + "/kraken2_db/taxonomy/nodes.dmp",
		#names = config["rdir"] + "/kraken2_db/taxonomy/names.dmp",
		#gtdb_map = config["rdir"] + "/kraken2_db/library/gtdb/prelim_map.txt",
		#gtdb_fasta = config["rdir"] + "/kraken2_db/library/gtdb/library.fna",
		#ncbi_map = expand(config["rdir"] + "/kraken2_db/library/{library_name}/prelim_map.txt", library_name = LIBRARY_NAME),
		#ncbi_fasta = expand(config["rdir"] + "/kraken2_db/library/{library_name}/library.fna", library_name = LIBRARY_NAME),
		#univec_map = config["rdir"] + "/kraken2_db/library/" + config["univec"] + "/prelim_map.txt" if config["univec"] else [],
		#univec_fasta = config["rdir"] + "/kraken2_db/library/" + config["univec"] + "/library.fna" if config["univec"] else []
		#hash = config["rdir"] + "/kraken2_db/hash.k2d",
		#opts = config["rdir"] + "/kraken2_db/opts.k2d",
		#map  = config["rdir"] + "/kraken2_db/seqid2taxid.map",
		#taxo = config["rdir"] + "/kraken2_db/taxo.k2d",
		ncbi_tax_coarse = expand(config["cdir"] + "/{library_name}/assembly_taxonomy_coarse.txt", library_name = LIBRARY_NAME),
		coarse_ncbi_download = expand(config["cdir"] + "/{library_name}/genomes/done", library_name = LIBRARY_NAME),
		gtdb_reps = config["cdir"] + "/gtdb/gtdb_reps_tax.txt",
		checkv_tax_coarse = config["cdir"] + "/checkv/checkv_reps_tax.txt",
		tax_all_coarse = config["cdir"] + "/tax_coarse_all.txt",
		file_list = config["cdir"] + "/kraken2_genomes/file_names_derep_genomes.txt"

'''
##### load rules #####
'''
#include: "rules/process_ncbi_genomes.smk"
#include: "rules/process_gtdb_genomes.smk"
#include: "rules/process_checkv_genomes.smk"
#include: "rules/build_kraken2.smk"
include: "rules/coarse_genomes.smk"
include: "rules/coarse_kraken2.smk"

