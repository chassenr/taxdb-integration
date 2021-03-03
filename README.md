# Overview
Current metagenomic reference databases are often either limited in their taxonomic representation, restricted to a sepcific domain of life, biased towards model taxa, prone to contamination, or too large to handle computationally. This workflow aims to address these limitations by providing a flexible, user-configurable pipeline to generate cross-domain databases that allow for the classification of metagenomic reads at a high specificity and sensitivity. Its main feature include:
* Integration of reference genomes from various sources ([GTDB](https://gtdb.ecogenomic.org/), [CheckV](https://bitbucket.org/berkeleylab/checkv/src/master/), NCBI refseq and genbank) and the option to include user-specified custom assemblies. 
* 2-step classification process where reads are first sorted by domain or phylum using a cross-domain database, followed by the classification at a high taxonomic resolution of reads within specific target groups (e.g. microbes: archaea, bacteria, viruses, fungi, protists), thereby reducing computational demands without lowering performance
* Even representation of known biodiversity by a [network-guided dereplication of genomes](https://github.com/genomewalker/derep-genomes) based on ANI
* Detection and removal of cross-domain contamination in the reference genomes using [conterminator](https://github.com/martin-steinegger/conterminator)
* Highly customizable workflow depending on individual user needs
* Portable: The workflow is implemented in snakemake and relies heavily on conda environments
The output of the workflow will be 2 [kraken2](https://github.com/DerrickWood/kraken2) and 2 [kaiju](https://github.com/bioinformatics-centre/kaiju) databases: (1) coarse taxonomic resolution for sorting by domain or phylum, (2) high resolution for target group. The workflow further uses scripts from [Correcting index databases to improve metagenomic studies](https://www.biorxiv.org/content/10.1101/712166v1). Therefore, if you use this workflow, please cite:

* Méric, G., Wick, R. R., Watts, S. C., Holt, K. E., & Inouye, M. (2019). Correcting index databases improves metagenomic studies. BioRxiv, 712166. https://doi.org/10.1101/712166
* Wood, D. E., Lu, J., & Langmead, B. (2019). Improved metagenomic analysis with Kraken 2. Genome Biology, 20(257), 1–13. https://doi.org/10.1186/s13059-019-1891-0
* Steinegger, M., Salzberg, S. L. (2020). Terminating contamination: large-scale search identifies more than 2,000,000 contaminated entries in GenBank. Genome Biol 21,115. https://doi.org/10.1186/s13059-020-02023-1
* Menzel, P., Krogh, A. (2016) Kaiju : Fast and sensitive taxonomic classification for metagenomics. Nat Commun 7: 11257. http://dx.doi.org/10.1038/ncomms11257

### Installation and dependencies

You will need python 3 (e.g. via [anaconda](https://docs.anaconda.com/anaconda/install/linux/) or its light-weight version [miniconda](https://docs.conda.io/en/latest/miniconda.html)) and [snakemake](https://snakemake.readthedocs.io/en/stable/) installed on your system to run this workflow. All further dependencies will be installed into individual conda environments. To get started, clone this repository:
```
git clone https://github.com/chassenr/taxdb-integration.git
```

You will also need to compile the latest version of conterminator as described [here](https://github.com/martin-steinegger/conterminator#optional-install-by-compilation). 

### Workflow description and configuration

![Illustration of the main steps in the workflow]()

The first step in the workflow is to screen the available collection of genomes on NCBI refseq and genbank based on their assembly quality. This first filter is based on the metadata parameter [assembly_level](https://www.ncbi.nlm.nih.gov/assembly/help/). Possible choices are: ```Contig```, ```Scaffold```, ```Chromosome```, ```Complete Genome```, or ```variable```. ```Contig``` is the most relaxed setting, with which all available genomes will be considered. ```Complete Genome``` is the most stringent setting. If ```variable```, for species with ```Chromosome``` and/or ```Complete Genome```-level assemblies available ```Contig``` and ```Scaffold```-level assemblies will not be included. This constitutes a compromise between excluding incomplete and low-quality assemblies and including species with incomplete assemblies. Apart from the setting for assembly level, this workflow is only considering the latest assembly of NCBI genomes with full genome representation.

While it is possible to restrict the taxonomic groups to be represented in the coarse database, it is generally recommended to build this database covering all domains of life. As a database with genomes from all available species would not be computationally feasible, the coase database will consist of a subset of the known diversity for each larger taxonomic group. Bacteria and Archaea will be represented by the GTDB representative genomes, viruses by the CehckV cluster representatives, and eukaryotes by genomes at a user-defined taxonomic resolution. For instance for higher eukaryotes (metazoa and streptophytes) each class will be represented by n genomes depending on the number of orders within that class, whereas for fungi, protists and algae this selection may be done at family level. To change these setting, it is possible to specify the taxonomic level(```rank_coarse```) in the workflow configuration for each of the main NCBI partitions (fungi, invertebrates, plants, protozoa, vertebrate_other, vertebrate_mammalian). To distinguish between lineages within each of these partitions, a seperate values can be specified for ```sublineage``` (e.g. to separate higher plants from algae). 

The coarse database generation also automatically runs conterminator to detect and remove contigs in assemblies with cross-domain contamination. The default here is to compare prokaryotes, fungi, protists, highler plants, and metazoa. Viruses are not included in the decontamination. Furthermore, common vector sequences ([UniVec](https://www.ncbi.nlm.nih.gov/tools/vecscreen/univec/)) are included in the database.

To combine all these genomes within a common taxonomic framework, we use the script tax_from_gtdb.py available [here](https://github.com/rrwick/Metagenomics-Index-Correction/blob/master/tax_from_gtdb.py).

For the high resolution database, available genomes can further be filtered based on assembly statistics that are calculated by [statswrapper.sh](https://github.com/BioInfoTools/BBMap/blob/master/sh/statswrapper.sh) of the [bbmap suite](https://sourceforge.net/projects/bbmap/), e.g. genome size, GC content, contig N50. It is possible to set absolute as well as relative cut-offs based on the data distribution (quantiles).

Before the netword guided dereplication of genomes within each species, it is possible to add genomes that are not included in the aboce data sources, e.g. high quality MAGs previously assembled from the studies environment, to increase the diversity coverage of the database. Also here, it is possible to exclude lineages that should not be included in the high resolution database (e.g. Streptophyta if your target group is restricted to microbes).

The remaining step in the generation of the high resolution database are similar to those for the coarse database: create a common taxonomic framework, remove contamination, and build the kraken database. As it is possible to restrict the database to only one domain or larger taxonomic group, the decontamination step can also be disabled. 

For building the kraken2 database, the parameters kmer length, minimizer length, and minimizer spaces




### Next steps and ToDos:
* Kaiju implementation (update prokaryotic annotation)
* Extend framework to Bracken and KrakenUniq databases
* Systematic parameter sweep to fine-tune defaults and assess sensitivity and specificity of databse

### Configuration
To customize the workflow to your user environment, you will need to modify the [config.yaml](https://github.com/chassenr/taxdb-integration/blob/master/config/config.yaml) and [cluster.yaml](https://github.com/chassenr/taxdb-integration/blob/master/config/cluster.yaml) files. 

The [cluster.yaml](https://github.com/chassenr/taxdb-integration/blob/master/config/cluster.yaml) file is specifying how the individual steps of the workflow will be distributed on your compute cluster. It defines, e.g., maximum run time and maximum number of cpus to use per node. At the moment, the custer configuration is customized for [slurm](https://slurm.schedmd.com/documentation.html) as workload manager. If you are not running the workflow on a cluster, this file will be ignored.

The [config.yaml](https://github.com/chassenr/taxdb-integration/blob/master/config/config.yaml) sets the parameters for the workflow. Here you will need to modify the following values:
* Working directory (wdir): Absolute path to this repository.
* Results directory (rdir): Absolute path to the location where the database and all associated files will be created.
* library_name: NCBI partitions to be included in the database (has to be one or more of ```fungi```, ```invertebrate```, ```plant```, ```protozoa```, ```vertebrate-mammalian```, ```vertebrate-other```, ```viral```).
* download_threads: Number of parallel processes for downloading genomes.
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

