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

rule all:
	input:
		genome_assemblies = config["rdir"] + "/" + config["library_name"] + "/assembly_summary_combined.txt",
		download_complete = config["rdir"] + "/" + config["library_name"] + "/genomes/done",
		genome_taxonomy = config["rdir"] + "/" + config["library_name"] + "/assembly_taxonomy.txt",
		derep_taxonomy = config["rdir"] + "/" + config["library_name"] + "/derep_assembly_taxonomy.txt"
        

'''
##### load rules #####
'''
include: "rules/process_ncbi_genomes.smk"

