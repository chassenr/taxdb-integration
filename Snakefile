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
		fasta = expand(config["rdir"] + "/kraken2_db/library/{library_highres}/library.fna", library_highres = LIBRARY_HIGHRES),
		cleaned_gtdb = config["rdir"] + "/decontamination/gtdb_cleaned.fna" if config["kingdoms_highres"] else [],
		cleaned_ncbi = expand(config["rdir"] + "/decontamination/{library_highres}_cleaned.fna", library_highres = LIBRARY_HIGHRES) if config["kingdoms_highres"] else [],
		hash = config["rdir"] + "/kraken2_db/hash.k2d"
		

'''
##### load rules #####
'''
include: "rules/process_ncbi_genomes.smk"
include: "rules/process_gtdb_genomes.smk"
include: "rules/process_checkv_genomes.smk"
include: "rules/build_kraken2.smk"
#include: "rules/coarse_genomes.smk"
#include: "rules/coarse_kraken2.smk"

