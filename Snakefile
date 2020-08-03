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

rule all:
	input:
		genome_assemblies = expand(config["rdir"] + "/{library_name}/assembly_summary_combined.txt", library_name = LIBRARY_NAME),
		gtdb_derep_tax = config["rdir"] + "/tax_combined/gtdb_derep_taxonomy.txt",
		ncbi_derep_tax = config["rdir"] + "/tax_combined/ncbi_derep_taxonomy.txt",
		nodes = config["rdir"] + "/kraken2_db/taxonomy/nodes.dmp",
		names = config["rdir"] + "/kraken2_db/taxonomy/names.dmp",
		hash = config["rdir"] + "/kraken2_db/hash.k2d",
		opts = config["rdir"] + "/kraken2_db/opts.k2d",
		map  = config["rdir"] + "/kraken2_db/seqid2taxid.map",
		taxo = config["rdir"] + "/kraken2_db/taxo.k2d"

'''
##### load rules #####
'''
include: "rules/process_ncbi_genomes.smk"
include: "rules/process_gtdb_genomes.smk"
include: "rules/build_kraken2.smk"
