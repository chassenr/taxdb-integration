# database summary stats
mkdir /vol/cloud/christiane/NCBI_taxdb_integration/Testing4/Summary_stats
cd /vol/cloud/christiane/NCBI_taxdb_integration/Testing4/tax_combined

# ANI threshold 0.1 for non-microbial taxa (on genus level)
# ANI 0.05 for viruses
# ANI 0.03 for fungi, protozoa, micro-algae
# ANI 0.01 for GTDB

cat > ../Summary_stats/tmp_taxlevels.txt
phylum
class
order
family
genus
species
genomes

# number of taxa per partition
PART="fungi invertebrate plant protozoa vertebrate_mammalian vertebrate_other viral gtdb"
for part in $(echo "$PART")
do
  for i in $(seq 2 7)
  do 
    cut -f2 ${part}"_derep_taxonomy.txt" | cut -d';' -f-$i | sort | uniq | wc -l
  done > ../Summary_stats/"tmp"
  cat ${part}"_derep_taxonomy.txt" | wc -l >> ../Summary_stats/tmp
  paste ../Summary_stats/tmp_taxlevels.txt ../Summary_stats/tmp > ../Summary_stats/${part}"_tax_stats.txt"
  rm ../Summary_stats/tmp
done

# number of taxa per partition before dereplication
PART="fungi invertebrate plant protozoa vertebrate_mammalian vertebrate_other viral"
for part in $(echo "$PART")
do
  for i in $(seq 2 7)
  do 
    cut -f2 ../${part}/assembly_taxonomy.txt | cut -d';' -f-$i | sort | uniq | wc -l
  done > ../Summary_stats/"tmp"
  cat ../${part}/assembly_taxonomy.txt | wc -l >> ../Summary_stats/tmp
  paste ../Summary_stats/tmp_taxlevels.txt ../Summary_stats/tmp > ../Summary_stats/${part}"_tax_all.txt"
  rm ../Summary_stats/tmp
done
for i in $(seq 2 7)
do 
  cut -f2 ../gtdb/metadata/gtdb_taxonomy.txt | cut -d';' -f-$i | sort | uniq | wc -l
done > ../Summary_stats/"tmp"
cat ../gtdb/metadata/gtdb_taxonomy.txt | wc -l >> ../Summary_stats/tmp
paste ../Summary_stats/tmp_taxlevels.txt ../Summary_stats/tmp > ../Summary_stats/gtdb_tax_all.txt
rm ../Summary_stats/tmp

# size if only one genome (largest) per taxlevel is retained
PART="fungi invertebrate plant protozoa vertebrate_mammalian vertebrate_other viral gtdb"
for part in $(echo "$PART")
do
  for i in $(seq 2 7)
  do 
    cut -f1 ${part}"_derep_taxonomy.txt" |  sed -e 's/^/\.\.\/kraken2_genomes\/genome_files\//' -e 's/$/\.fa/'| xargs ls -l | sed 's/\s\+/\t/g' | cut -f5,9 | sed -e 's/\.\.\/kraken2_genomes\/genome_files\///' -e 's/\.fa//' | sort -k2,2 | paste - <(sort -k1,1 ${part}"_derep_taxonomy.txt" | cut -f2 | cut -d';' -f-$i) | sort -t$'\t' -k3,3 -k1,1gr | sort -t$'\t' -u -k3,3 --merge | cut -f1 | paste -sd+ | bc -l
  done > ../Summary_stats/"tmp"
  cut -f1 ${part}"_derep_taxonomy.txt" |  sed -e 's/^/\.\.\/kraken2_genomes\/genome_files\//' -e 's/$/\.fa/'| xargs ls -l | sed 's/\s\+/\t/g' | cut -f5 | paste -sd+ | bc -l >> ../Summary_stats/tmp
  paste ../Summary_stats/tmp_taxlevels.txt ../Summary_stats/tmp > ../Summary_stats/${part}"_size_stats.txt"
  rm ../Summary_stats/tmp
done


