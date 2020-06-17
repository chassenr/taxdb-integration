# download non-prokaryotic refseq partitions from NCBI
# see: https://github.com/DerrickWood/kraken2/blob/master/scripts/download_genomic_library.sh

# Open questions:
#   take all assemblies or only complete ones (depending on refseq partition)?
#   supplement refseq with high quality genbank assemblies?

# testing with fungi section
LIBRARY_DIR="/home/chh/Documents/Projects/NCBI_taxdb_integration"
NCBI_SERVER="ftp.ncbi.nlm.nih.gov"
FTP_SERVER="ftp://$NCBI_SERVER"
library_name="fungi"

cd $LIBRARY_DIR
mkdir $library_name
cd $library_name

# get list of assembly locations
wget -q ${FTP_SERVER}"/genomes/refseq/$library_name/assembly_summary.txt"
mv assembly_summary.txt assembly_summary_refseq.txt

# also download genbank files?
# maybe restrict to assembly level Chromosome, Complete Genome, (Scaffold), i.e. exclude Contig?
# only take those not already included in refseq
wget -q ${FTP_SERVER}"/genomes/genbank/$library_name/assembly_summary.txt"
mv assembly_summary.txt assembly_summary_genbank.txt

# filter genbank assemblies by
#   assembly_level != Contig (remove don't be as stringest here...)
#   genome_rep == Full
#   version_status == latest
#   redundant information with refseq
sed '/^\#/d' assembly_summary_genbank.txt | awk -v FS='\t' -v OFS='\t' '$11 == "latest" && $12 != "Contig" && $14 == "Full"' | grep -v -f <(cut -f1 assembly_summary_refseq.txt) > assembly_summary_genbank_filtered.txt

# combine assembly summary tables
sed '/^\#/d' assembly_summary_refseq.txt > assembly_summary_combined.txt
cat assembly_summary_genbank_filtered.txt >> assembly_summary_combined.txt

# for testing only: take some unique and highly replicated genomes
cut -f7 assembly_summary_combined.txt | sort | uniq -c | sed -e 's/\s\+/\t/g' -e 's/^\t//' | awk '$1 == 1' | head -20 | cut -f2 > taxid_select.txt
cut -f7 assembly_summary_combined.txt | sort | uniq -c | sed -e 's/\s\+/\t/g' -e 's/^\t//' | awk '$1 >= 3 && $1 <= 5' | head -10 | cut -f2 >> taxid_select.txt
cut -f7 assembly_summary_combined.txt | sort | uniq -c | sed -e 's/\s\+/\t/g' -e 's/^\t//' | awk '$1 >= 20' | head -5 | cut -f2 >> taxid_select.txt

# filter assembly_summary table
while read line
do
  awk -v FS='\t' -v OFS='\t' -v taxid=${line} '$7 == taxid' assembly_summary_combined.txt >> assembly_summary_subset.txt
done < taxid_select.txt

# download *genomic.fna.gz files for each genome
cut -f20 assembly_summary_subset.txt | sed -e 's/.*\///' -e 's/$/_genomic.fna.gz/' | paste -d'/' <(cut -f20 assembly_summary_subset.txt) - > assembly_url_genomic.txt
mkdir genomes
cd genomes
cat ../assembly_url_genomic.txt | xargs -n 1 -P 6 wget -q

# get taxonomic path from species taxid
# instead of using efetch, map directly using taxdump files?
# wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz

# download files for kaiju for dereplicated assemblies?
