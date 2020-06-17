# dereplicate NCBI refseq assemblies
# see: https://github.com/rrwick/Metagenomics-Index-Correction/blob/master/dereplicate_assemblies.py
# this script can be used directly :)

/home/chh/Documents/Repos/Metagenomics-Index-Correction/dereplicate_assemblies.py --threads 4 --threshold 0.005 genomes derep_genomes fungi_taxonomy.txt

# only select dereplicated genomes from taxonomy table
ls -1 derep_genomes/ | sed 's/\([0-9]\)_.*/\1/' | grep -F -f - fungi_taxonomy.txt > fungi_taxonomy_dereplicated.txt

