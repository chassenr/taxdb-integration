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
LIBRARY_COARSE = config["library_coarse"]

wildcard_constraints:
	library_name = '|'.join(LIBRARY_NAME),
	library_highres = '|'.join(LIBRARY_HIGHRES),
	library_coarse = '|'.join(LIBRARY_COARSE)

rule all:
	input:
		taxonomy = expand(config["rdir"] + "/{library_name}/assembly_taxonomy.txt", library_name = LIBRARY_NAME),
		download_complete = expand(config["rdir"] + "/{library_name}/genomes/done", library_name = LIBRARY_NAME),
		tax_filt = expand(config["rdir"] + "/{library_name}/assembly_taxonomy_filtered.txt", library_name = LIBRARY_NAME)
		

'''
##### load rules #####
'''
include: "rules/process_ncbi_genomes.smk"
#include: "rules/process_gtdb_genomes.smk"
#include: "rules/process_checkv_genomes.smk"
#include: "rules/build_kraken2.smk"
#include: "rules/process_ncbi_proteins.smk"
#include: "rules/coarse_genomes.smk"
#include: "rules/coarse_kraken2.smk"
#include: "rules/coarse_proteins.smk"

