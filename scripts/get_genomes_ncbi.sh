#!/bin/bash

# download non-prokaryotic refseq partitions from NCBI
# adapted from: https://github.com/DerrickWood/kraken2/blob/master/scripts/download_genomic_library.sh

# Open questions:
#   take all assemblies or only complete ones (depending on refseq partition)?

# reading input options
library_dir=$1
library_name=$2 # has to be one of: fungi, invertebrate, plant, protozoa, vertebrate_mammalian, vertebrate_other, viral
assembly_level=$3 # has to be one of Contig, Scaffold, Chromosome, Complete Genome (indicates minimum assembly level, i.e. if "Contig" all assemblies will be included)

# set download options
NCBI_SERVER="http://ftp.ncbi.nlm.nih.gov"

# set up directory for each library
mkdir -p "${library_dir}"/"${library_name}"

# get list of refseq assembly locations from assembly summary file
aria2c --max-tries=20 --retry-wait=5 --dir "${library_dir}"/"${library_name}" --out assembly_summary_refseq.txt "${NCBI_SERVER}"/genomes/refseq/"${library_name}"/assembly_summary.txt

# get list of genbank assembly locations from assembly summary file
aria2c --max-tries=20 --retry-wait=5 --dir "${library_dir}"/"${library_name}" --out assembly_summary_genbank.txt "${NCBI_SERVER}"/genomes/genbank/"${library_name}"/assembly_summary.txt

# remove genbank genomes from list, which are redundant with refseq
sed '/^\#/d' "${library_dir}"/"${library_name}"/assembly_summary_refseq.txt | cut -f1 | grep -v -f - "${library_dir}"/"${library_name}"/assembly_summary_genbank.txt > "${library_dir}"/"${library_name}"/assembly_summary_genbank_nr.txt

# only use the latest assembly version of genbank assemblies with full genome representation
sed '/^\#/d' "${library_dir}"/"${library_name}"/assembly_summary_genbank_nr.txt | awk -v FS='\t' -v OFS='\t' '$11 == "latest" && $14 == "Full"' > "${library_dir}"/"${library_name}"/assembly_summary_genbank_clean.txt

# filter genbank assemblies further by assembly_level
# (I did not know how to nicely make this selection in one condition...)
if [ "${assembly_level}" = "Contig" ]
then
  cat "${library_dir}"/"${library_name}"/assembly_summary_genbank_clean.txt > "${library_dir}"/"${library_name}"/assembly_summary_genbank_filt.txt
fi
if [ "${assembly_level}" = "Scaffold" ]
then
  awk -v FS='\t' -v OFS='\t' '$12 != "Contig"' "${library_dir}"/"${library_name}"/assembly_summary_genbank_clean.txt > "${library_dir}"/"${library_name}"/assembly_summary_genbank_filt.txt
fi
if [ "${assembly_level}" = "Chromosome" ]
then
  awk -v FS='\t' -v OFS='\t' '$12 == "Chromosome" || $12 == "Complete Genome"' "${library_dir}"/"${library_name}"/assembly_summary_genbank_clean.txt > "${library_dir}"/"${library_name}"/assembly_summary_genbank_filt.txt
fi
if [ "${assembly_level}" = "Complete Genome" ]
then
  awk -v FS='\t' -v OFS='\t' '$12 == "Complete Genome"' "${library_dir}"/"${library_name}"/assembly_summary_genbank_clean.txt > "${library_dir}"/"${library_name}"/assembly_summary_genbank_filt.txt
fi

# combine assembly summary tables
sed '/^\#/d' "${library_dir}"/"${library_name}"/assembly_summary_refseq.txt > "${library_dir}"/"${library_name}"/assembly_summary_combined.txt
cat "${library_dir}"/"${library_name}"/assembly_summary_genbank_filt.txt >> "${library_dir}"/"${library_name}"/assembly_summary_combined.txt

# format input file for aria2
cut -f20 "${library_dir}"/"${library_name}"/assembly_summary_combined.txt | sed -e 's/.*\///' -e 's/$/_genomic.fna.gz/' | paste -d'/' <(cut -f20 "${library_dir}"/"${library_name}"/assembly_summary_combined.txt) - | sed 's/^ftp/http/' > "${library_dir}"/"${library_name}"/assembly_url_genomic.txt


