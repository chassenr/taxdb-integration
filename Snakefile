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
		derep_gtdb = config["rdir"] + "/gtdb/metadata/gtdb_derep_taxonomy_meta.txt",
		derep_ncbi = expand(config["rdir"] + "/{library_highres}/derep_taxonomy_meta.txt", library_highres = LIBRARY_HIGHRES),
		derep_checkv = config["rdir"] + "/checkv/checkv_derep_taxonomy_meta.txt",
		tax_coarse = config["cdir"] + "/tax_coarse_all.txt"

'''
##### load rules #####
'''
include: "rules/process_ncbi_genomes.smk"
include: "rules/process_gtdb_genomes.smk"
include: "rules/process_checkv_genomes.smk"
#include: "rules/build_kraken2.smk"
include: "rules/coarse_genomes.smk"
#include: "rules/coarse_kraken2.smk"

