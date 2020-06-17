# merge databases and format kraken2 compatible files

# copy all genome files into one directory (will be deleted again later)
cd /home/chh/Documents/Projects/NCBI_taxdb_integration/
mkdir genomes_all
cp ./fungi/derep_genomes/*.gz ./genomes_all/

# for testing: get the first 10 gtdb bac genomes
cd /home/chh/Documents/Projects/NCBI_taxdb_integration/GTDB
wget -q ${FTP_SERVER}"/genomes/refseq/bacteria/assembly_summary.txt"
mv assembly_summary.txt assembly_summary_refseq.txt
cut -f1 bac120_taxonomy_r89.tsv | shuf -n 20 | sed 's/^RS_//' | grep -F -f - assembly_summary_refseq.txt > tmp
cut -f20 tmp | sed -e 's/.*\///' -e 's/$/_genomic.fna.gz/' | paste -d'/' <(cut -f20 tmp) - | xargs -n 1 -P 6 wget -q
mv *.gz ../genomes_all/
cut -f1 tmp | grep -F -f - bac120_taxonomy_r89.tsv > gtdb_taxonomy.txt
cd ../

# combine taxonomy tables
cat ./GTDB/gtdb_taxonomy.txt ./fungi/fungi_taxonomy_dereplicated.txt > combined_taxonomy.txt

# generate kraken2 compatible files
/home/chh/Documents/Repos/Metagenomics-Index-Correction/tax_from_gtdb.py --gtdb combined_taxonomy.txt --assemblies genomes_all --nodes nodes.dmp --names names.dmp --kraken_dir kraken2_genomes
mkdir kraken2_db
mkdir kraken2_db/taxonomy
cp nodes.dmp kraken2_db/taxonomy/
cp names.dmp kraken2_db/taxonomy/

# create kraken2 database
# for snakemake, use: https://github.com/leylabmpi/Struo/blob/f8fdf3d6f04678502fb8d6b094cb4135b7c361e3/bin/kraken2/Snakefile
for fa in ./kraken2_genomes/*.fa
do
  kraken2-build --db kraken2_db --add-to-library $fa
done
kraken2-build --build --threads 1 --db kraken2_db
# error with more than 1 thread
