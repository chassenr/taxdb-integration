import glob
import pandas as pd
"""
Author: 
Affiliation: 
Aim: 
Run: snakemake -s Snakefile
"""
# from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
# min_version("5.1.2")

configfile: "config/config.yaml"
# report: "report/workflow.rst"

# This should be placed in the Snakefile.

'''
Working directory
'''
workdir: config["wdir"]
# message("The current working directory is " + WDIR)
'''
 The list of samples to be processed
'''
LIBRARY_NAME = ["fungi","invertebrate","viral"]
# LIBRARY_NAME = ["fungi","invertebrate","plant","protozoa","vertebrate_mammalian","vertebrate_other","viral"]

rule all:
	input:
		genome_assemblies = expand(config["rdir"] + "/" + "{library_name}" + "/assembly_summary_combined.txt", library_name = LIBRARY_NAME),
		download_complete = expand(config["rdir"] + "/" + "{library_name}" + "/genomes/done", library_name = LIBRARY_NAME),
		genome_taxonomy = expand(config["rdir"] + "/" + "{library_name}" + "/assembly_taxonomy.txt", library_name = LIBRARY_NAME),
		derep_taxonomy = expand(config["rdir"] + "/" + "{library_name}" + "/derep_assembly_taxonomy.txt", library_name = LIBRARY_NAME)
        

'''
##### load rules #####
'''
include: "rules/process_ncbi_genomes.smk"

