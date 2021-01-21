# database validation

### simulate metaG reads from genomes not in the database
# select from assembly_taxonomy files accnos which are not in derep_assembly_taxonomy
# inlcude human as well (to make sure that this contamination will be recognized)
# select 2 random genomes per phylum
cd /vol/cloud/christiane/NCBI_taxdb_integration/Testing4
mkdir Simulation

PART="fungi invertebrate plant protozoa vertebrate_mammalian vertebrate_other viral"
cut -f1 tax_combined/gtdb_derep_taxonomy.txt | grep -v -F -f - gtdb/metadata/gtdb_taxonomy.txt > Simulation/non_db_genomes_taxonomy.txt
for part in $(echo "$PART")
do
  cut -f1 tax_combined/${part}"_derep_taxonomy.txt" | grep -v -F -f - ${part}/assembly_taxonomy.txt >> Simulation/non_db_genomes_taxonomy.txt
done

cd Simulation
sort -t$'\t' -k2,2 non_db_genomes_taxonomy.txt | sort -t$'\t' -u -k2,2 --merge > non_db_genomes_taxonomy_nr.txt
rm sim_genomes_taxonomy.txt
cut -f2 non_db_genomes_taxonomy_nr.txt | cut -d';' -f-2 | sort -t$'\t' | uniq > phylum_list.txt
while read line
do
  grep "${line}" non_db_genomes_taxonomy_nr.txt | shuf -n 2 >> sim_genomes_taxonomy.txt
done < phylum_list.txt
mv sim_genomes_taxonomy.txt tmp
sort -t$'\t' -k2,2 tmp | uniq > sim_genomes_taxonomy.txt # shuf may duplicate lines...
rm tmp

HUMAN=$(grep "s__Homo sapiens" sim_genomes_taxonomy.txt)
if [[ $HUMAN == "" ]]
then
  grep "s__Homo sapiens" non_db_genomes_taxonomy.txt | shuf -n 1 >> sim_genomes_taxonomy.txt
fi

cut -f2 ../gtdb/metadata/gtdb_download_info.txt | cat - ../*/assembly_url_genomic.txt | grep -F -f <(cut -f1 sim_genomes_taxonomy.txt) - > sim_genomes_download.txt
sed 's/.*\///' sim_genomes_download.txt | cut -d'_' -f2 | grep -v -F -f - sim_genomes_taxonomy.txt | cut -f1 | cut -d'_' -f2 | sed 's/\.[0-9]\+//' | grep -F -f - <(cut -f2 ../gtdb/metadata/gtdb_download_info.txt) >> sim_genomes_download.txt
mkdir genomes
aria2c -i sim_genomes_download.txt -c -l links.log --dir genomes --max-tries=20 --retry-wait=5 --max-connection-per-server=1 --max-concurrent-downloads=28
gunzip genomes/*.gz

# metaG simulator
# https://github.com/nick-youngblut/MGSIM
# install dependencies first, then run pip install command
# maybe better to use setup.py?

# format genome table
# Important: no duplicate taxa allowed in genome table
paste <(sed 's/.*\//genomes\//' sim_genomes_download.txt | sed 's/\.gz//') <(sed 's/.*\//genomes\//' sim_genomes_download.txt | cut -d'_' -f2 | sed 's/\.[0-9]\+//') | sort -t$'\t' -k2,2 > tmp1
cut -f1 sim_genomes_taxonomy.txt | cut -d'_' -f2 | sed 's/\.[0-9]\+//' | paste sim_genomes_taxonomy.txt - | sort -t$'\t' -k3,3 > tmp2
diff <(cut -f2 tmp1) <(cut -f3 tmp2)
echo -e "Accnos\tTaxon\tFasta" > sim_genomes_table.txt
cut -f2,1 tmp2 | paste - <(cut -f1 tmp1) >> sim_genomes_table.txt
rm tmp1 tmp2

# run community simulation
MGSIM communities sim_genomes_table.txt sim_out
MGSIM reads --sr-seq-depth 1e7 --art-paired -n 16 sim_genomes_table.txt sim_out_abund.txt sim_out_reads

# simulate shorter reads (80bp single-end)
MGSIM reads --sr-seq-depth 1e7 --art-len=80 --art-mflen=0 -n 16 sim_genomes_table.txt sim_out_abund.txt sim_out_reads_short_SE


### create a database subset with low resolution for euks and high resolution for microbes
cat ../tax_combined/gtdb_derep_taxonomy.txt ../tax_combined/viral_derep_taxonomy.txt ../tax_combined/fungi_derep_taxonomy.txt ../tax_combined/protozoa_derep_taxonomy.txt | cut -f1 > db_subset_micro.accnos
grep -v "p__Streptophyta" ../tax_combined/plant_derep_taxonomy.txt | cut -f1 >> db_subset_micro.accnos
grep -F -f db_subset_micro.accnos ../tax_combined/derep_taxonomy_combined.txt > db_subset_micro_taxonomy.txt

grep "p__Streptophyta" ../tax_combined/plant_derep_taxonomy.txt | cut -f1 > db_subset_macro.accnos
cut -f1 ../tax_combined/invertebrate_derep_taxonomy.txt >> db_subset_macro.accnos
cut -f1 ../tax_combined/vertebrate_mammalian_derep_taxonomy.txt >> db_subset_macro.accnos
cut -f1 ../tax_combined/vertebrate_other_derep_taxonomy.txt >> db_subset_macro.accnos
grep -F -f db_subset_macro.accnos ../tax_combined/derep_taxonomy_combined.txt > db_subset_macro_taxonomy.txt

# estimate sizes
for i in $(seq 2 7)
do 
  cut -f1 db_subset_macro_taxonomy.txt |  sed -e 's/^/\.\.\/kraken2_genomes\/genome_files\//' -e 's/$/\.fa/' | xargs ls -l | sed 's/\s\+/\t/g' | cut -f5,9 | sed -e 's/\.\.\/kraken2_genomes\/genome_files\///' -e 's/\.fa//' | sort -k2,2 | paste - <(sort -k1,1 db_subset_macro_taxonomy.txt | cut -f2 | cut -d';' -f-$i) | sort -t$'\t' -k3,3 -k1,1gr | sort -t$'\t' -u -k3,3 --merge | cut -f1 | paste -sd+ | bc -l
done > tmp
cut -f1 db_subset_macro_taxonomy.txt |  sed -e 's/^/\.\.\/kraken2_genomes\/genome_files\//' -e 's/$/\.fa/' | xargs ls -l | sed 's/\s\+/\t/g' | cut -f5 | paste -sd+ | bc -l >> tmp
paste ../Summary_stats/tmp_taxlevels.txt tmp > db_subset_macro_sizes.txt
rm tmp
cut -f1 db_subset_micro_taxonomy.txt |  sed -e 's/^/\.\.\/kraken2_genomes\/genome_files\//' -e 's/$/\.fa/' | xargs ls -l | sed 's/\s\+/\t/g' | cut -f5 | paste -sd+ | bc -l
# 321.5 GB micro + 173.5 GB macro (class level): too large without memory limitation
# only one representative per species for non-viral microbes
grep -v "d__Virus" db_subset_micro_taxonomy.txt > db_subset_micro_taxonomy_nv.txt
cut -f1 db_subset_micro_taxonomy_nv.txt | sed -e 's/^/\.\.\/kraken2_genomes\/genome_files\//' -e 's/$/\.fa/' | xargs ls -l | sed 's/\s\+/\t/g' | cut -f5,9 | sed -e 's/\.\.\/kraken2_genomes\/genome_files\///' -e 's/\.fa//' | sort -k2,2 | paste - <(sort -k1,1 db_subset_micro_taxonomy_nv.txt | cut -f2) | sort -t$'\t' -k3,3 -k1,1gr | sort -t$'\t' -u -k3,3 --merge | cut -f1 | paste -sd+ | bc -l
grep "d__Virus" db_subset_micro_taxonomy.txt > db_subset_micro_taxonomy_vi.txt
cut -f1 db_subset_micro_taxonomy_vi.txt |  sed -e 's/^/\.\.\/kraken2_genomes\/genome_files\//' -e 's/$/\.fa/' | xargs ls -l | sed 's/\s\+/\t/g' | cut -f5 | paste -sd+ | bc -l
# 232.6 GB micro: that should work with class-level macro

# select largest genome per class for macro
cut -f1 db_subset_macro_taxonomy.txt |  sed -e 's/^/\.\.\/kraken2_genomes\/genome_files\//' -e 's/$/\.fa/' | xargs ls -l | sed 's/\s\+/\t/g' | cut -f5,9 | sed -e 's/\.\.\/kraken2_genomes\/genome_files\///' -e 's/\.fa//' | sort -k2,2 | paste - <(sort -k1,1 db_subset_macro_taxonomy.txt | cut -f2 | cut -d';' -f-3) | sort -t$'\t' -k3,3 -k1,1gr | sort -t$'\t' -u -k3,3 --merge | cut -f2 | grep -F -f - db_subset_macro_taxonomy.txt > db_subset_class_macro_taxonomy.txt

# select largest genome per species for non-viral micro
cut -f1 db_subset_micro_taxonomy_nv.txt |  sed -e 's/^/\.\.\/kraken2_genomes\/genome_files\//' -e 's/$/\.fa/' | xargs ls -l | sed 's/\s\+/\t/g' | cut -f5,9 | sed -e 's/\.\.\/kraken2_genomes\/genome_files\///' -e 's/\.fa//' | sort -k2,2 | paste - <(sort -k1,1 db_subset_micro_taxonomy_nv.txt | cut -f2) | sort -t$'\t' -k3,3 -k1,1gr | sort -t$'\t' -u -k3,3 --merge | cut -f2 | grep -F -f - db_subset_micro_taxonomy_nv.txt > db_subset_species_micro_taxonomy_nv.txt

# combine final tax tables
cat db_subset_class_macro_taxonomy.txt db_subset_species_micro_taxonomy_nv.txt db_subset_micro_taxonomy_vi.txt > db_subset_v1_taxonomy.txt
cut -f1 db_subset_v1_taxonomy.txt |  sed -e 's/^/\.\.\/kraken2_genomes\/genome_files\//' -e 's/$/\.fa/' | xargs ls -l | sed 's/\s\+/\t/g' | cut -f5 | paste -sd+ | bc -l

# prepare db directory (already has UniVec_Core)
cp -r ../kraken2_db/ ./kraken2_db1

# mask low complexity and build prelim_map
mkdir kraken2_db1/library/class/ # IMPORTANT: don't call directory 'added'
cut -f1 db_subset_v1_taxonomy.txt |  sed -e 's/^/\.\.\/kraken2_genomes\/genome_files\//' -e 's/$/\.fa/' | parallel -j28 'dustmasker -in {} -outfmt fasta' | sed -e '/^>/!s/[a-z]/x/g' >> kraken2_db1/library/class/library.fna
LC_ALL=C grep '^>' kraken2_db1/library/class/library.fna | sed 's/^>//' > kraken2_db1/library/class/tmp.accnos
NSEQ=$(wc -l kraken2_db1/library/added/tmp.accnos | cut -d' ' -f1)
printf 'TAXID\n%.0s' $(seq 1 $NSEQ) | paste - kraken2_db1/library/class/tmp.accnos | paste - <(cut -d'|' -f3 kraken2_db1/library/class/tmp.accnos) > kraken2_db1/library/class/prelim_map.txt
rm kraken2_db1/library/class/tmp.accnos

# build database
# use screen
scontrol -o show nodes
salloc --cpus-per-task=28 --nodelist="bibigrid-worker1-3-jdzmrjtottq9ctr" bash
conda activate source-tracking
srun /usr/bin/time -v kraken2-build --build --threads 26 --db kraken2_db1 --kmer-len 31 --minimizer-len 25 --minimizer-spaces 5 --max-db-size 460000000000 >>db1_build.log 2>&1

# run classification
/usr/bin/time -v kraken2 --db kraken2_db1 --threads 10 --report sim_v1.kreport --output sim_v1.kraken --paired sim_out_reads/1/R1.fq sim_out_reads/1/R2.fq >>sim_v1_classify.log 2>&1
/usr/bin/time -v kraken2 --db kraken2_db1 --threads 10 --report sim_v1_short.kreport --output sim_v1_short.kraken sim_out_reads_short_SE/1/R1.fq >>sim_v1_short_classify.log 2>&1


# with memory limitation
# also build size limited database to be run on 256GB RAM nodes
mkdir kraken2_db2
mv kraken2_db1/taxonomy/ kraken2_db2/
mv kraken2_db1/library/ kraken2_db2/
scontrol -o show nodes
salloc --cpus-per-task=28 --nodelist="bibigrid-worker1-3-jdzmrjtottq9ctr" bash
conda activate source-tracking
srun /usr/bin/time -v kraken2-build --build --threads 26 --db kraken2_db2 --kmer-len 31 --minimizer-len 25 --minimizer-spaces 5 --max-db-size 200000000000 >>db2_build.log 2>&1
/usr/bin/time -v kraken2 --db kraken2_db2 --threads 10 --report sim_v2.kreport --output sim_v2.kraken --paired sim_out_reads/1/R1.fq sim_out_reads/1/R2.fq >>sim_v2_classify.log 2>&1
/usr/bin/time -v kraken2 --db kraken2_db2 --threads 10 --report sim_v2_short.kreport --output sim_v2_short.kraken sim_out_reads_short_SE/1/R1.fq >>sim_v2_short_classify.log 2>&1


### inspecting results 
# calculate individual ANI
mkdir Check_output
mash sketch -p 16 -o Check_output/Galdieria_sulphuraria -s 10000 genomes/GCA_001704855.1_ASM170485v1_genomic.fna ../kraken2_genomes/genome_files/GCA_006232545.1.fa
mash dist -p 16 Check_output/Galdieria_sulphuraria.msh Check_output/Galdieria_sulphuraria.msh > Check_output/Galdieria_sulphuraria_dist.txt
fastANI -q genomes/GCA_001704855.1_ASM170485v1_genomic.fna -r ../kraken2_genomes/genome_files/GCA_006232545.1.fa -o Check_output/Galdieria_sulphuraria_fastani.txt


### testing 2-step DB with 70bp simulated reads

# simulating reads
MGSIM reads --sr-seq-depth 1e7 --art-len=70 --art-mflen=0 -n 16 sim_genomes_table.txt sim_out_abund.txt sim_out_reads_70bp_SE

# building coarse database
mkdir kraken2_coarse_k25
mv kraken2_db2/taxonomy/ kraken2_coarse_k25/
mv kraken2_db2/library/ kraken2_coarse_k25/
scontrol -o show nodes
salloc --cpus-per-task=28 --nodelist="bibigrid-worker1-1-jdzmrjtottq9ctr" bash
conda activate source-tracking
srun /usr/bin/time -v kraken2-build --build --threads 26 --db kraken2_coarse_k25 --kmer-len 25 --minimizer-len 23 --minimizer-spaces 4 --max-db-size 460000000000 >>db_coarse_build.log 2>&1
# this taking forever
# try on worker 3 and with vmtouch and latest kraken2 version
cp -r kraken2_coarse_k25/ kraken2_coarse2_k25
salloc --cpus-per-task=28 --nodelist="bibigrid-worker1-3-jdzmrjtottq9ctr" bash
cd /vol/cloud/christiane/Repos/vmtouch
srun ./vmtouch -l -d -t -f -v /vol/cloud/christiane/NCBI_taxdb_integration/Testing4/Simulation/kraken2_coarse2_k25
srun /usr/bin/time -v /vol/cloud/christiane/Software/kraken2_v2.1.1/kraken2-build --build --threads 26 --db kraken2_coarse2_k25 --kmer-len 25 --minimizer-len 23 --minimizer-spaces 4 --max-db-size 460000000000 >>db_coarse2_build.log 2>&1
# not working...
# Antonio is building this at MPI

# also build db with default params
salloc --cpus-per-task=28 --nodelist="bibigrid-worker1-3-jdzmrjtottq9ctr" bash
export PATH="/vol/cloud/christiane/Software/kraken2_v2.1.1:$PATH"
srun /usr/bin/time -v kraken2-build --build --threads 26 --db kraken2_coarse_default --fast-build >>db_coarse_default_build.log 2>&1



# building high-res database (GTDB and virus only)
cp -r ../kraken2_db/ ./kraken2_highres_k25 # with UniVec
mkdir kraken2_highres_k25/library/highres/
cat ../tax_combined/gtdb_derep_taxonomy.txt ../tax_combined/viral_derep_taxonomy.txt | cut -f1 | grep -F -f - ../tax_combined/derep_taxonomy_combined.txt > db_highres_micro_taxonomy.txt
cut -f1 db_highres_micro_taxonomy.txt | grep -F -f - ../kraken2_genomes/file_names_derep_genomes.txt | parallel -j26 'dustmasker -in {} -outfmt fasta' | sed -e '/^>/!s/[a-z]/x/g' >> kraken2_highres_k25/library/highres/library.fna
LC_ALL=C grep '^>' kraken2_highres_k25/library/highres/library.fna | sed 's/^>//' > kraken2_highres_k25/library/highres/tmp.accnos
NSEQ=$(wc -l kraken2_highres_k25/library/highres/tmp.accnos | cut -d' ' -f1)
printf 'TAXID\n%.0s' $(seq 1 $NSEQ) | paste - kraken2_highres_k25/library/highres/tmp.accnos | paste - <(cut -d'|' -f3 kraken2_highres_k25/library/highres/tmp.accnos) > kraken2_highres_k25/library/highres/prelim_map.txt
rm kraken2_highres_k25/library/highres/tmp.accnos
scontrol -o show nodes
salloc --cpus-per-task=28 --nodelist="bibigrid-worker1-3-jdzmrjtottq9ctr" bash
conda activate source-tracking
srun /usr/bin/time -v kraken2-build --build --threads 26 --db kraken2_highres_k25 --kmer-len 25 --minimizer-len 23 --minimizer-spaces 4 --max-db-size 460000000000 >>db_highres_build.log 2>&1
# with new kraken2 version
export PATH="/vol/cloud/christiane/Software/kraken2_v2.1.1:$PATH"
srun /usr/bin/time -v kraken2-build --build --threads 26 --db kraken2_highres_k25 --kmer-len 25 --minimizer-len 23 --minimizer-spaces 4 --max-db-size 430000000000 --fast-build
# also build db with default params
salloc --cpus-per-task=28 --nodelist="bibigrid-worker1-1-jdzmrjtottq9ctr" bash
export PATH="/vol/cloud/christiane/Software/kraken2_v2.1.1:$PATH"
srun /usr/bin/time -v kraken2-build --build --threads 26 --db kraken2_highres_default --fast-build >>db_highres_default_build.log 2>&1


# classification coarse
srun /usr/bin/time -v kraken2 --db /vol/cloud/DBs/kraken2_coarse_k25 --threads 10 --report sim_70bp_SE_coarse.kreport --output sim_70bp_SE_coarse.kraken --report-minimizer-data sim_out_reads_70bp_SE/1/R1.fq >>sim_70bp_SE_coarse_classify.log 2>&1

# classification high-res
srun /usr/bin/time -v kraken2 --db kraken2_highres_k25 --threads 10 --report sim_70bp_SE_highres.kreport --output sim_70bp_SE_highres.kraken --report-minimizer-data sim_out_reads_70bp_SE/1/R1.fq >>sim_70bp_SE_highres_classify.log 2>&1

srun /usr/bin/time -v kraken2 --db kraken2_highres_default --threads 10 --report sim_70bp_SE_highres_default.kreport --output sim_70bp_SE_highres_default.kraken --report-minimizer-data sim_out_reads_70bp_SE/1/R1.fq >>sim_70bp_SE_highres_default_classify.log 2>&1


### build bracken for highres k25 and 70bp
export PATH="/vol/cloud/christiane/Repos/Bracken:/vol/cloud/christiane/Software/kraken2_v2.1.1:$PATH"
srun /usr/bin/time -v bracken-build -d kraken2_highres_k25 -t 26 -k 25 -l 70


### run full parameter sweep
# kmer length: 25 and 35
# minimizer length: 23
# minimizer spaces: 4 and 7
# max db size: inf and 220GB

# Hint: python list for parameters and for combinations use the product function from itertools
# Snakemake wildcards seems to do the trick :)
cd /vol/cloud/christiane/Repos/taxdb-integration
snakemake -s Snakefile_param_sweep --use-conda -j 100 --cluster-config config/cluster.yaml --cluster "sbatch --export=ALL -t {cluster.time} --ntasks-per-node {cluster.ntasks_per_node} --nodes {cluster.nodes} --cpus-per-task {cluster.cpus_per_task} --partition {cluster.partition} --job-name {rulename}.{jobid}"

# coarse k35
scontrol -o show nodes
salloc --cpus-per-task=28 --nodelist="bibigrid-worker1-3-jdzmrjtottq9ctr" bash
conda activate /vol/cloud/christiane/Repos/taxdb-integration/.snakemake/conda/0a1123c6
srun /usr/bin/time -v kraken2-build --build --threads 28 --db /vol/cloud/christiane/NCBI_taxdb_integration/Testing4/Simulation/db_coarse_k35_m23_s4_nolim --kmer-len 35 --minimizer-len 23 --minimizer-spaces 4 --fast-build &>> /vol/cloud/christiane/NCBI_taxdb_integration/Testing4/Simulation/ps_db_build_coarse_k35_m23_s4_nolim.log

# coarse k25 with mem220
srun /usr/bin/time -v kraken2-build --build --threads 28 --db /vol/cloud/christiane/NCBI_taxdb_integration/Testing4/Simulation/db_coarse_k25_m23_s4_mem220 --kmer-len 25 --minimizer-len 23 --minimizer-spaces 4 --max-db-size 220000000000 --fast-build &>> /vol/cloud/christiane/NCBI_taxdb_integration/Testing4/Simulation/ps_db_build_coarse_k25_m23_s4_mem220.log

# classification coarse
scontrol -o show nodes
salloc --cpus-per-task=28 --nodelist="bibigrid-worker1-3-jdzmrjtottq9ctr" bash
conda activate /vol/cloud/christiane/Repos/taxdb-integration/.snakemake/conda/0a1123c6
# k35 (nolim and mem220 will be identical)
srun /usr/bin/time -v kraken2 --db db_coarse_k35_m23_s4_nolim --threads 20 --report sim_70bp_SE_coarse_k35_m23_s4_nolim.kreport --output sim_70bp_SE_coarse_k35_m23_s4_nolim.kraken --report-minimizer-data sim_out_reads_70bp_SE/1/R1.fq >>sim_70bp_SE_coarse_k35_m23_s4_nolim_classify.log 2>&1
# k25 mem220
srun /usr/bin/time -v kraken2 --db db_coarse_k25_m23_s4_mem220 --threads 20 --report sim_70bp_SE_coarse_k25_m23_s4_mem220.kreport --output sim_70bp_SE_coarse_k25_m23_s4_mem220.kraken --report-minimizer-data sim_out_reads_70bp_SE/1/R1.fq >>sim_70bp_SE_coarse_k25_m23_s4_mem220_classify.log 2>&1


# run conifer on kraken output
# conifer does not change taxid assignment (not even with filter option)
# coarse k25 nolim
srun /vol/cloud/christiane/Repos/Conifer/conifer --all --both_scores -i sim_70bp_SE_coarse.kraken -d /vol/cloud/DBs/kraken2_coarse_k25/taxo.k2d > sim_70bp_SE_coarse.conifer_out
# coarse k25 mem220
srun /vol/cloud/christiane/Repos/Conifer/conifer --all --both_scores -i sim_70bp_SE_coarse_k25_m23_s4_mem220.kraken -d db_coarse_k25_m23_s4_mem220/taxo.k2d > sim_70bp_SE_coarse_k25_m23_s4_mem220.conifer_out
# coarse k35
srun /vol/cloud/christiane/Repos/Conifer/conifer --all --both_scores -i sim_70bp_SE_coarse_k35_m23_s4_nolim.kraken -d db_coarse_k35_m23_s4_nolim/taxo.k2d > sim_70bp_SE_coarse_k35_m23_s4_nolim.conifer_out
# highres k25
srun /vol/cloud/christiane/Repos/Conifer/conifer --all --both_scores -i sim_70bp_SE_k25_m23_s4_nolim.kraken -d db_k25_m23_s4_nolim/taxo.k2d > sim_70bp_SE_k25_m23_s4_nolim.conifer_out
# highres k35
srun /vol/cloud/christiane/Repos/Conifer/conifer --all --both_scores -i sim_70bp_SE_k35_m23_s4_nolim.kraken -d db_k35_m23_s4_nolim/taxo.k2d > sim_70bp_SE_k35_m23_s4_nolim.conifer_out

# also run kraken2 classify with confidence thresholds
# coarse k35, confidence 0.1 (rtl 0.25)
srun /usr/bin/time -v kraken2 --db db_coarse_k35_m23_s4_nolim --threads 20 --confidence 0.1 --report sim_70bp_SE_coarse_k35_m23_s4_nolim_c01.kreport --output sim_70bp_SE_coarse_k35_m23_s4_nolim_c01.kraken --report-minimizer-data sim_out_reads_70bp_SE/1/R1.fq >>sim_70bp_SE_coarse_k35_m23_s4_nolim_c01_classify.log 2>&1
/vol/cloud/christiane/Repos/Conifer/conifer --all --both_scores -i sim_70bp_SE_coarse_k35_m23_s4_nolim_c01.kraken -d db_coarse_k35_m23_s4_nolim/taxo.k2d > sim_70bp_SE_coarse_k35_m23_s4_nolim_c01.conifer_out
# coarse k25 mem220, confidence 0.1 (rtl 0.25)
srun /usr/bin/time -v kraken2 --db db_coarse_k25_m23_s4_mem220 --threads 20 --confidence 0.1 --report sim_70bp_SE_coarse_k25_m23_s4_mem220_c01.kreport --output sim_70bp_SE_coarse_k25_m23_s4_mem220_c01.kraken --report-minimizer-data sim_out_reads_70bp_SE/1/R1.fq >>sim_70bp_SE_coarse_k25_m23_s4_mem220_c01_classify.log 2>&1
/vol/cloud/christiane/Repos/Conifer/conifer --all --both_scores -i sim_70bp_SE_coarse_k25_m23_s4_mem220_c01.kraken -d db_coarse_k25_m23_s4_mem220/taxo.k2d > sim_70bp_SE_coarse_k25_m23_s4_mem220_c01.conifer_out
# coarse k25 nolim, confidence 0.1 (rtl 0.25)
srun /usr/bin/time -v kraken2 --db /vol/cloud/DBs/kraken2_coarse_k25 --threads 20 --confidence 0.1 --report sim_70bp_SE_coarse_c01.kreport --output sim_70bp_SE_coarse_c01.kraken --report-minimizer-data sim_out_reads_70bp_SE/1/R1.fq >>sim_70bp_SE_coarse_c01_classify.log 2>&1
/vol/cloud/christiane/Repos/Conifer/conifer --all --both_scores -i sim_70bp_SE_coarse_c01.kraken -d /vol/cloud/DBs/kraken2_coarse_k25/taxo.k2d > sim_70bp_SE_coarse_c01.conifer_out

salloc --cpus-per-task=28 --nodelist="bibigrid-worker1-1-jdzmrjtottq9ctr" bash
conda activate /vol/cloud/christiane/Repos/taxdb-integration/.snakemake/conda/0a1123c6
# highres k25, confidence 0.1 (rtl 0.25)
srun /usr/bin/time -v kraken2 --db db_k25_m23_s4_nolim --threads 1 --confidence 0.1 --report sim_70bp_SE_k25_m23_s4_nolim_c01.kreport --output sim_70bp_SE_k25_m23_s4_nolim_c01.kraken --report-minimizer-data sim_out_reads_70bp_SE/1/R1.fq >>ps_classify_k25_m23_s4_nolim_c01.log 2>&1
/vol/cloud/christiane/Repos/Conifer/conifer --all --both_scores -i sim_70bp_SE_k25_m23_s4_nolim_c01.kraken -d db_k25_m23_s4_nolim/taxo.k2d > sim_70bp_SE_k25_m23_s4_nolim_c01.conifer_out
# highres k35, confidence 0.1 (rtl 0.25)
srun /usr/bin/time -v kraken2 --db db_k35_m23_s4_nolim --threads 1 --confidence 0.1 --report sim_70bp_SE_k35_m23_s4_nolim_c01.kreport --output sim_70bp_SE_k35_m23_s4_nolim_c01.kraken --report-minimizer-data sim_out_reads_70bp_SE/1/R1.fq >>ps_classify_k35_m23_s4_nolim_c01.log 2>&1
/vol/cloud/christiane/Repos/Conifer/conifer --all --both_scores -i sim_70bp_SE_k35_m23_s4_nolim_c01.kraken -d db_k35_m23_s4_nolim/taxo.k2d > sim_70bp_SE_k35_m23_s4_nolim_c01.conifer_out



### repeat benchmark with 25% random sequences

# generate random genome (based on E.coli https://www.ncbi.nlm.nih.gov/assembly/GCF_000005845.2)
# needs to have file ending .fasta

# run random sequence generation
python /vol/cloud/christiane/Repos/NullSeq/nullseq.py -n 1 --seq Ecoli_GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fasta -o genomes/random_genomic.fna
# this is not working as apparently NullSeq only works with CDS and not full genomes
# work with CDS
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < GCF_000005845.2_ASM584v2_protein.faa | sed '1d' > test_protein_lin.fasta
mkdir tmp
split -d -l 2 test_protein_lin.fasta tmp/tmp
cd tmp
ls | parallel -j20 'python /vol/cloud/christiane/Repos/NullSeq/nullseq.py -n 1 --AA {} --GC 50' >> ../test_cds_random.fna
# not working with GNU parallel this way... 

# let's try randomseq
# linearize fasta
cat Ecoli_GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fasta | tr '\n' '\t' | sed -e 's/\t/__/' -e 's/\t//g' -e 's/__/\t/' -e 's/$/\t/' | tr '\t' '\n' > Ecoli_GCF_000005845.2/test.fasta
python /vol/cloud/christiane/Repos/bactome/randomseq.py shuffle --sequence="$(tail -1 Ecoli_GCF_000005845.2/test.fasta)" > Ecoli_GCF_000005845.2/random.seq
# argument list too long for bash

# let's shuffle in R (using sample without replacement)
cp Ecoli_GCF_000005845.2/random.fasta genomes/random.fna

# build sim data set with 25% random, 25% euks, 10% arch, 10% viral, 30% bac
MGSIM communities sim_genomes_table.txt sim_out_wr
# let's use the proportion suggested for teh individual genomes, but scale according to domain proportion
# in R
x <- read.table("sim_out_wr_abund.txt", sep = "\t", h = T, stringsAsFactors = F)
domains <- c(25, 25, 10, 30, 10)
names(domains) <- c("d__random", "d__Eukaryota", "d__Archaea", "d__Bacteria", "d__Viruses")
x.adj <- vector("list", length = length(domains))
for(i in 1:length(domains)) {
  tmp <- x[grepl(names(domains)[i], x$Taxon), ]
  tmp$Perc_rel_abund <- tmp$Perc_rel_abund/sum(tmp$Perc_rel_abund) * domains[i]
  x.adj[[i]] <- tmp
}
x.adj.df <- do.call("rbind", x.adj)
x.adj.df <- x.adj.df[order(x.adj.df$Perc_rel_abund, decreasing = T),]
x.adj.df$Rank <- 1:nrow(x.adj.df)
write.table(x.adj.df, "sim_out_wr_abund_adj.txt", quote = F, sep = "\t", row.names = F)

# simulate 70bp SE reads
MGSIM reads --sr-seq-depth 1e7 --art-len=70 --art-mflen=0 -n 16 sim_genomes_table.txt sim_out_wr_abund_adj.txt sim_out_wr_reads_70bp_SE

# classification
scontrol -o show nodes
salloc --cpus-per-task=24 --nodelist="bibigrid-worker1-3-jdzmrjtottq9ctr" bash
conda activate /vol/cloud/christiane/Repos/taxdb-integration/.snakemake/conda/0a1123c6

# without confidence threshold
# coarse k25 nolim
srun /usr/bin/time -v kraken2 --db /vol/cloud/DBs/kraken2_coarse_k25 --threads 20 --report sim_70bp_SE_wr_coarse_k25_nolim.kreport --output sim_70bp_SE_wr_coarse_k25_nolim.kraken --report-minimizer-data sim_out_wr_reads_70bp_SE/1/R1.fq >>sim_70bp_SE_wr_coarse_k25_nolim_classify.log 2>&1
srun /vol/cloud/christiane/Repos/Conifer/conifer --all --both_scores -i sim_70bp_SE_wr_coarse_k25_nolim.kraken -d /vol/cloud/DBs/kraken2_coarse_k25/taxo.k2d > sim_70bp_SE_wr_coarse_k25_nolim.conifer_out
# coarse k25 mem220
srun /usr/bin/time -v kraken2 --db db_coarse_k25_m23_s4_mem220 --threads 20 --report sim_70bp_SE_wr_coarse_k25_mem220.kreport --output sim_70bp_SE_wr_coarse_k25_mem220.kraken --report-minimizer-data sim_out_wr_reads_70bp_SE/1/R1.fq >>sim_70bp_SE_wr_coarse_k25_mem220_classify.log 2>&1
srun /vol/cloud/christiane/Repos/Conifer/conifer --all --both_scores -i sim_70bp_SE_wr_coarse_k25_mem220.kraken -d db_coarse_k25_m23_s4_mem220/taxo.k2d > sim_70bp_SE_wr_coarse_k25_mem220.conifer_out
# coarse k35 (nolim and mem220 will be identical)
srun /usr/bin/time -v kraken2 --db db_coarse_k35_m23_s4_nolim --threads 20 --report sim_70bp_SE_wr_coarse_k35_nolim.kreport --output sim_70bp_SE_wr_coarse_k35_nolim.kraken --report-minimizer-data sim_out_wr_reads_70bp_SE/1/R1.fq >>sim_70bp_SE_wr_coarse_k35_nolim_classify.log 2>&1
srun /vol/cloud/christiane/Repos/Conifer/conifer --all --both_scores -i sim_70bp_SE_wr_coarse_k35_nolim.kraken -d db_coarse_k35_m23_s4_nolim/taxo.k2d > sim_70bp_SE_wr_coarse_k35_nolim.conifer_out
# highres k25 (nolim and mem220 will be identical)
srun /usr/bin/time -v kraken2 --db db_k25_m23_s4_nolim --threads 20 --report sim_70bp_SE_wr_highres_k25_nolim.kreport --output sim_70bp_SE_wr_highres_k25_nolim.kraken --report-minimizer-data sim_out_wr_reads_70bp_SE/1/R1.fq >>sim_70bp_SE_wr_highres_k25_nolim_classify.log 2>&1
srun /vol/cloud/christiane/Repos/Conifer/conifer --all --both_scores -i sim_70bp_SE_wr_highres_k25_nolim.kraken -d db_k25_m23_s4_nolim/taxo.k2d > sim_70bp_SE_wr_highres_k25_nolim.conifer_out
# highres k35 (nolim and mem220 will be identical)
srun /usr/bin/time -v kraken2 --db db_k35_m23_s4_nolim --threads 20 --report sim_70bp_SE_wr_highres_k35_nolim.kreport --output sim_70bp_SE_wr_highres_k35_nolim.kraken --report-minimizer-data sim_out_wr_reads_70bp_SE/1/R1.fq >>sim_70bp_SE_wr_highres_k35_nolim_classify.log 2>&1
srun /vol/cloud/christiane/Repos/Conifer/conifer --all --both_scores -i sim_70bp_SE_wr_highres_k35_nolim.kraken -d db_k35_m23_s4_nolim/taxo.k2d > sim_70bp_SE_wr_highres_k35_nolim.conifer_out

# with confidence threshold
# coarse k25 nolim
srun /usr/bin/time -v kraken2 --db /vol/cloud/DBs/kraken2_coarse_k25 --threads 20 --confidence 0.1 --report sim_70bp_SE_wr_coarse_k25_nolim_c01.kreport --output sim_70bp_SE_wr_coarse_k25_nolim_c01.kraken --report-minimizer-data sim_out_wr_reads_70bp_SE/1/R1.fq >>sim_70bp_SE_wr_coarse_k25_nolim_c01_classify.log 2>&1
srun /vol/cloud/christiane/Repos/Conifer/conifer --all --both_scores -i sim_70bp_SE_wr_coarse_k25_nolim_c01.kraken -d /vol/cloud/DBs/kraken2_coarse_k25/taxo.k2d > sim_70bp_SE_wr_coarse_k25_nolim_c01.conifer_out
# coarse k25 mem220
srun /usr/bin/time -v kraken2 --db db_coarse_k25_m23_s4_mem220 --threads 20 --confidence 0.1 --report sim_70bp_SE_wr_coarse_k25_mem220_c01.kreport --output sim_70bp_SE_wr_coarse_k25_mem220_c01.kraken --report-minimizer-data sim_out_wr_reads_70bp_SE/1/R1.fq >>sim_70bp_SE_wr_coarse_k25_mem220_c01_classify.log 2>&1
srun /vol/cloud/christiane/Repos/Conifer/conifer --all --both_scores -i sim_70bp_SE_wr_coarse_k25_mem220_c01.kraken -d db_coarse_k25_m23_s4_mem220/taxo.k2d > sim_70bp_SE_wr_coarse_k25_mem220_c01.conifer_out
# coarse k35 (nolim and mem220 will be identical)
srun /usr/bin/time -v kraken2 --db db_coarse_k35_m23_s4_nolim --threads 20 --confidence 0.1 --report sim_70bp_SE_wr_coarse_k35_nolim_c01.kreport --output sim_70bp_SE_wr_coarse_k35_nolim_c01.kraken --report-minimizer-data sim_out_wr_reads_70bp_SE/1/R1.fq >>sim_70bp_SE_wr_coarse_k35_nolim_c01_classify.log 2>&1
srun /vol/cloud/christiane/Repos/Conifer/conifer --all --both_scores -i sim_70bp_SE_wr_coarse_k35_nolim_c01.kraken -d db_coarse_k35_m23_s4_nolim/taxo.k2d > sim_70bp_SE_wr_coarse_k35_nolim_c01.conifer_out
# highres k25 (nolim and mem220 will be identical)
srun /usr/bin/time -v kraken2 --db db_k25_m23_s4_nolim --threads 20 --confidence 0.1 --report sim_70bp_SE_wr_highres_k25_nolim_c01.kreport --output sim_70bp_SE_wr_highres_k25_nolim_c01.kraken --report-minimizer-data sim_out_wr_reads_70bp_SE/1/R1.fq >>sim_70bp_SE_wr_highres_k25_nolim_c01_classify.log 2>&1
srun /vol/cloud/christiane/Repos/Conifer/conifer --all --both_scores -i sim_70bp_SE_wr_highres_k25_nolim_c01.kraken -d db_k25_m23_s4_nolim/taxo.k2d > sim_70bp_SE_wr_highres_k25_nolim_c01.conifer_out
# highres k35 (nolim and mem220 will be identical)
srun /usr/bin/time -v kraken2 --db db_k35_m23_s4_nolim --threads 20 --confidence 0.1 --report sim_70bp_SE_wr_highres_k35_nolim_c01.kreport --output sim_70bp_SE_wr_highres_k35_nolim_c01.kraken --report-minimizer-data sim_out_wr_reads_70bp_SE/1/R1.fq >>sim_70bp_SE_wr_highres_k35_nolim_c01_classify.log 2>&1
srun /vol/cloud/christiane/Repos/Conifer/conifer --all --both_scores -i sim_70bp_SE_wr_highres_k35_nolim_c01.kraken -d db_k35_m23_s4_nolim/taxo.k2d > sim_70bp_SE_wr_highres_k35_nolim_c01.conifer_out

# with higher confidence threshold in kraken2
# not worth it...
# with lower conf
# coarse k25 mem220
srun /usr/bin/time -v kraken2 --db db_coarse_k25_m23_s4_mem220 --threads 20 --confidence 0.05 --report sim_70bp_SE_wr_coarse_k25_mem220_c005.kreport --output sim_70bp_SE_wr_coarse_k25_mem220_c005.kraken --report-minimizer-data sim_out_wr_reads_70bp_SE/1/R1.fq >>sim_70bp_SE_wr_coarse_k25_mem220_c005_classify.log 2>&1
srun /vol/cloud/christiane/Repos/Conifer/conifer --all --both_scores -i sim_70bp_SE_wr_coarse_k25_mem220_c005.kraken -d db_coarse_k25_m23_s4_mem220/taxo.k2d > sim_70bp_SE_wr_coarse_k25_mem220_c005.conifer_out
# rtl 0.2 better than confidence 0.005

