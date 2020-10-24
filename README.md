# Overview
This workflow will format NCBI genomes of non-prokaryotic RefSeq and Genbank partitions to integrate with [GTDB](https://gtdb.ecogenomic.org/). The output will be a [kraken2](https://github.com/DerrickWood/kraken2) database for all domains of life. The workflow is based on the idea of [Correcting index databases to improve metagenomic studies](https://www.biorxiv.org/content/10.1101/712166v1) and relies on some of the scripts used for that paper. Therefore, if you use this workflow, please also cite:

* Méric, G., Wick, R. R., Watts, S. C., Holt, K. E., & Inouye, M. (2019). Correcting index databases improves metagenomic studies. BioRxiv, 712166. https://doi.org/10.1101/712166

As well as [kraken2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1891-0):

* Wood, D. E., Lu, J., & Langmead, B. (2019). Improved metagenomic analysis with Kraken 2. Genome Biology, 20(257), 1–13. https://doi.org/10.1186/s13059-019-1891-0

### Workflow summary
* Step 1a: Download metadata for all available genomes from RefSeq and Genbank and filter entries based on assembly level and genome representation to ensure a high quality of input genomes for the database construction. Download the assemblies for the selected genomes. Genomes can be included from the following NBCI partitions: fungi, invertebrate, plant, protozoa, vertebrate-mammalian, vertebrate-other, viral. Plasmid, chloroplast, mitochondrion databases are not included separetaly, since according to the NCBI readme: "The genomic.fna.gz file includes all top-level sequences in the assembly (chromosomes, plasmids, organelles, unlocalized scaffolds, unplaced scaffolds, and any alternate loci or patch scaffolds)".
* Step 1b: Download genomes included in [GTDB](https://gtdb.ecogenomic.org/) before species clustering.
* Step 1c (optional): Also include sequences from [UniVec](https://www.ncbi.nlm.nih.gov/tools/vecscreen/univec/) in the database.
* Step 2: Dereplicate genomes per species (for all domains) at selected average nucleotide identity (ANI) threshold using the approach implemented in [correct index databases](https://github.com/rrwick/Metagenomics-Index-Correction/blob/master/dereplicate_assemblies.py).
* Step 3: Parse parse taxonomy using [tax_from_gtdb.py](https://github.com/rrwick/Metagenomics-Index-Correction/blob/master/tax_from_gtdb.py) and build kraken2 database with all genomes (gtdb plus ncbi).

### Next steps and ToDos:
* Include dynamic dereplication approach based on ANI network.
* Create Bracken and KrakenUniq databases.
* Configure workflow to create 2 databases: coarse for sorting reads by lineages of interest, high-resolution for specific taxonomic assignment.
* Systematic parameter sweep to fine-tune defaults and assess sensitivity and specificity of databse.
* Use dereplicated genome selection to format database for [kaiju](https://github.com/bioinformatics-centre/kaiju) relying on available annotation for eukaryotic and viral genomes, but updating prokaryotic annotation.

# Workflow description

### Installation
You will need python 3 (e.g. via [anaconda](https://docs.anaconda.com/anaconda/install/linux/) or its light-weight version [miniconda](https://docs.conda.io/en/latest/miniconda.html)) and [snakemake](https://snakemake.readthedocs.io/en/stable/) installed on your system to run this workflow. All further dependencies will be installed into individual conda environments. To get started, clone this repository:
```
git clone https://github.com/chassenr/taxdb-integration.git
```

### Configuration
To customize the workflow to your user environment, you will need to modify the [config.yaml](https://github.com/chassenr/taxdb-integration/blob/master/config/config.yaml) and [cluster.yaml](https://github.com/chassenr/taxdb-integration/blob/master/config/cluster.yaml) files. 

The [cluster.yaml](https://github.com/chassenr/taxdb-integration/blob/master/config/cluster.yaml) file is specifying how the individual steps of the workflow will be distributed on your compute cluster. It defines, e.g., maximum run time and maximum number of cpus to use per node. At the moment, the custer configuration is customized for [slurm](https://slurm.schedmd.com/documentation.html) as workload manager. If you are not running the workflow on a cluster, this file will be ignored.

The [config.yaml](https://github.com/chassenr/taxdb-integration/blob/master/config/config.yaml) sets the parameters for the workflow. Here you will need to modify the following values:
* Working directory (wdir): Absolute path to this repository.
* Results directory (rdir): Absolute path to the location where the database and all associated files will be created.
* library_name: NCBI partitions to be included in the database (has to be one or more of ```fungi```, ```invertebrate```, ```plant```, ```protozoa```, ```vertebrate-mammalian```, ```vertebrate-other```, ```viral```).
* download_threads: Number of parallel processes for downloading genomes.
* assembly_level: Minimum required assembly level for each NCBI partition according to [NCBI](https://www.ncbi.nlm.nih.gov/assembly/help/) (has to be one of ```Contig```, ```Scaffold```, ```Chromosome```, ```Complete Genome```, ```variable```). With this parameter the quality of the assemblies considered for the database creation can be controlled. It may therefore strongly affect the size of your database. ```Contig``` is the most relaxed setting, with which all available genomes will be considered. ```Complete Genome``` is the most stringent setting. If ```variable```, for species with ```Chromosome``` and/or ```Complete Genome```-level assemblies available, ```Contig``` and ```Scaffold```-level assemblies will not be included. This constitutes a compromise between excluding incomplete and low-quality assemblies and including species with incomplete assemblies. Apart from the setting for assembly level, this workflow is only considering the latest assembly of NCBI genomes with full genome representation.
* derep_script: Location of script to perform genome dereplication available in the [correct index databases](https://github.com/rrwick/Metagenomics-Index-Correction) repository.
* derep_threads: Number of parallel processes for dereplicating genomes.
* derep_lineage: With this option it is possible to specify separate dereplication parameters for one (or more, '|' separated) sub-lineages within each NCBI parition.
* derep_threshold_gtdb and derep_threshold_ncbi_main/derep_threshold_ncbi_sub: ANI dereplication threshold for GTDB and each NCBI partition (separated into main and sub-lineages as specified in ```derep_lineage```. For each species (see ```derep_taxlevel_ncbi_main``` and ```derep_taxlevel_ncbi_sub```), only one genome per cluster defined at this ANI threshold will be included in the database. The lower the treshold value, the more likely it is that more genomes per species will be included in the database, which may increase your classification success. This parameter will strongly affect the size of your database. We suggest values between 0.005 and 0.05 as compromise between sensitivity and database size for microbial taxa. For eukaryotic genomes, we selected an ANI distance of 0.1, but it is also possible to go as lower. For viral genomes, we select an ANI distance threshold of 0.05 (only considering complete genomes) based on the recommendations in [Roux et al. 2018](https://www.nature.com/articles/nbt.4306). For GTDB we suggest to follow the recommendations in [Méric et al. 2019](https://www.biorxiv.org/content/10.1101/712166v1) or slightly higher (0.01) to avoid too large a database.
* derep_taxlevel_ncbi_main and derep_taxlevel_ncbi_sub: Specify the taxonomic level (1: domain, 7: species) at which the dereplication should be performed for the main and sub-lineages (see ```derep_lineage```). For higher ANI distance thresholds lower taxonomic levels are recommended, e.g. use ANI distance of 0.1 on genus level.
* tax_script: Location of script to prepare the required taxonomy files for kraken2 and format the genomes to be compatible with ```kraken2-build```, available in [correct index databases](https://github.com/rrwick/Metagenomics-Index-Correction) repository. 
* univec: Include UniVec (either ```UniVec``` or ```UniVec_Core```) in the database.
* krakenbuild_threads: Number of parallel processes for building kraken2 database.
* kmer_len, minimizer_len, minimizer_spaces: Settings for ```kraken2-build --build``` command. The current settings are slightly lower than the kraken2 defaults, which may considerably improve classification success without noticably increasing the rate of false positives (according the the benchmarks in the kraken2 paper by [Wood et al. 2019](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1891-0)).

### Running the workflow
After adapting the [config.yaml](https://github.com/chassenr/taxdb-integration/blob/master/config/config.yaml) and [cluster.yaml](https://github.com/chassenr/taxdb-integration/blob/master/config/cluster.yaml) files, the workflow can be executed with the following command from within this repository:
```
snakemake --use-conda -j 100 --cluster-config config/cluster.yaml --cluster "sbatch --export=ALL -t {cluster.time} --ntasks-per-node {cluster.ntasks_per_node} --nodes {cluster.nodes} --cpus-per-task {cluster.cpus_per_task} --partition {cluster.partition} --job-name {rulename}.{jobid}"
```
If you are not working on a cluster, only run: ```snakemake --use-conda -j 100```. You will need to adjust the maximum number of available threads to your system.

### Trouble-shooting
If the download of the genomes throws an error (without any obvious reason), thereby interrupting the snakemake execution, the following piece of code will add the missing genomes so that snakemake can continue with the next rule (example: downloading GTDB genomes):
```
# if aria2 error occurs (for unknown reasons)
# navigate to the genome directory, here: yourResultsDirectory/gtdb/genome
cut -f4 ../metadata/gtdb_download_info.txt | sort > tmp1
find . -type f -name '*.gz' | xargs -n 1 basename | sort > tmp2
grep -v -F -f tmp2 tmp1 > tmp3
if [ -s tmp3 ]
then
  grep -B1 -F -f tmp3 links | sed '/^--$/d' > retry_links
  aria2c -i retry_links -c -l retry_links.log --dir ./ --max-tries=20 --retry-wait=5 --max-connection-per-server=1 --max-concurrent-downloads=4 &>> ../../logs/download_gtdb_ncbi_retry.log
fi
# if download successful, generate input file for next snakemake rule and remove tmp files (keep logs for later debugging)
touch done
rm links retry_links tmp1 tmp2 tmp3
```

If you encounter any other bugs, please don't hesitate to open an issue.

