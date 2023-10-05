## :warning: This repository is not maintained anymore. Because of github enforcing 2FA, the repository moved to a less restrictive platform at: https://git.io-warnemuende.de/hassenru/taxdb-integration


# Overview
Current metagenomic reference databases are often either limited in their taxonomic representation, restricted to a sepcific domain of life, biased towards model taxa, prone to contamination, or too large to handle computationally. This workflow aims to address these limitations by providing a flexible, user-configurable pipeline to generate cross-domain databases that allow for the classification of metagenomic reads at a high specificity and sensitivity. 
To allow a flexible selection of genomes (nucleotide and protein sequences) for a variety of classification tools, the workflow is structured into modules.
Its main feature include:
* The basis of the workflow is the **COLLECT** module, which retrieves nucleotide and protein sequences from various sources, such as [GTDB](https://gtdb.ecogenomic.org/), [CheckV](https://bitbucket.org/berkeleylab/checkv/src/master/), NCBI refseq and genbank, and includes the option to provide additional user-specified custom assemblies as input. At the moment, already formatted taxonomy tables for the [TARA MAGs and SMAGs](https://www.genoscope.cns.fr/tara/) and [giant virus collection](https://figshare.com/ndownloader/files/32385995) are included, as well as resources from [EukProt](https://ndownloader.figshare.com/files/23580944) and [MMETSP](https://zenodo.org/record/3247846/files/mmetsp_dib_trinity_zenodo.fasta.tar.gz?download=1). Suitable genomes for the reference database are selected by a [network-guided dereplication](https://github.com/genomewalker/derep-genomes) based on ANI to achieve an even representation of the known biodiversity. The **COLLECT** module further generates a common taxonomic framework for all input genomes and provides NCBI-style taxdump files and accession2taxid mapping files as output.
* The output of the **COLLECT** module can then be used as input for several sequence classification tools. We provide here the workflow to generate [kraken2](https://github.com/DerrickWood/kraken2),  [kaiju](https://github.com/bioinformatics-centre/kaiju), [CAT/BAT](https://github.com/dutilh/CAT) databases offering several presets to determine the taxonomic resolution (and size) of the database. Specifically, we recommend a 2-step classification process for nucleotide level where reads are first sorted by domain using a cross-domain database, followed by the classification at a high taxonomic resolution of reads within each domain or larger taxonomic target group (e.g. prokaryotes and viruses, lower eukaryotes and plants, metazoa), thereby reducing computational demands without lowering performance.
* The kraken and kaiju database generation modules include the option to remove cross-domain contamination in the reference genomes using [conterminator](https://github.com/martin-steinegger/conterminator). This feature is experimental at the moment and needs further testing.
* Portable: The workflow is implemented in snakemake and relies heavily on conda environments


This workflow integrates several open source bioinformatic tools. Therefore, if you use this workflow, please cite:

* Wood, D. E., Lu, J., & Langmead, B. (2019). Improved metagenomic analysis with Kraken 2. Genome Biology, 20(257), 1–13. https://doi.org/10.1186/s13059-019-1891-0
* Steinegger, M., Salzberg, S. L. (2020). Terminating contamination: large-scale search identifies more than 2,000,000 contaminated entries in GenBank. Genome Biol 21,115. https://doi.org/10.1186/s13059-020-02023-1
* Menzel, P., Krogh, A. (2016) Kaiju : Fast and sensitive taxonomic classification for metagenomics. Nat Commun 7: 11257. http://dx.doi.org/10.1038/ncomms11257
* Steinegger, M. & Söding, J. (2017): MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nat Biotechnol 35:1026–1028. https://doi.org/10.1038/nbt.3988
* Hyatt, D., Chen, GL., LoCascio, P.F. et al. (2010) Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC Bioinformatics 11, 119. https://doi.org/10.1186/1471-2105-11-119


### Repository contents

* **Main directory**: Snakefiles for each workflow module
* **rules**: Snakemake rules for each module of the workflow
* **scripts**: Additional scripts called by the rules
* **envs**: Conda environment configuration files (.yaml) 
* **config**: Configuration files (.yaml) for each module of the workflow and for [SLURM](https://slurm.schedmd.com/overview.html) usage
* **assets**: Taxonomy files for additional user-supplied genomes and protein sequences
* **benchmark**: Results of the kraken2 and kaiju database performance benchmark (parameter sweep) with simulated metagenomic reads
* **stats**: Overview statistics of the DB contents (number of genomes/sequences per taxonomic level)


### Installation and dependencies

You will need python 3 (e.g. via [anaconda](https://docs.anaconda.com/anaconda/install/linux/) or its light-weight version [miniconda](https://docs.conda.io/en/latest/miniconda.html)) and [snakemake](https://snakemake.readthedocs.io/en/stable/) installed on your system to run this workflow. All further dependencies will be installed into individual conda environments. To get started, clone this repository:
```
git clone https://github.com/chassenr/taxdb-integration.git
```

You may also need to compile the latest version of conterminator as described [here](https://github.com/martin-steinegger/conterminator#optional-install-by-compilation). However, conterminator is disabled at the moment by default.


### Workflow description and configuration

![Illustration of the main steps in the workflow](https://github.com/chassenr/taxdb-integration/blob/master/images/taxdb_workflow.jpg)

Further details about the configurable parameters are included in the comments of the config files for each workflow module.

#### Collect module

Eukaryotic reference sequences are retrieved from NCBI. The first step in the workflow is to screen the available collection of genomes on NCBI refseq and genbank based on their assembly quality. This first filter is based on the metadata parameter [assembly_level](https://www.ncbi.nlm.nih.gov/assembly/help/). Possible choices are: ```Contig```, ```Scaffold```, ```Chromosome```, ```Complete Genome```, or ```variable```. ```Contig``` is the most relaxed setting, with which all available genomes will be considered. ```Complete Genome``` is the most stringent setting. If ```variable```, for species with ```Chromosome``` and/or ```Complete Genome```-level assemblies available ```Contig``` and ```Scaffold```-level assemblies will not be included. This constitutes a compromise between excluding incomplete and low-quality assemblies and including species with incomplete assemblies. Apart from the setting for assembly level, this workflow is only considering the latest assembly of NCBI genomes with full genome representation. Furthermore, if genomes with predicted proteins are available per species, assemblies without predicted proteins are excluded from downstream processing. User-supplied genomes are not subject to any filtering. After genome download, assemblies per species are further filtered by genome size, GC-content, and N50 to remove possible outliers. The strignency of this filter is user-configurable.

Prokaryotic and viral genomes are taken from the full GTDB and checkV databases. If necessary, proteins will be predicted with [prodigal](https://github.com/hyattpd/Prodigal).

To avoid redundancy, all genomes per species are dereplicated using [derepG](https://github.com/genomewalker/derep-genomes). Please refer to the program repository for further information about parameter setting.

User-supplied genomes and protein sequences can be added pre- or post-dereplication. The preformatted data sets for the TARA, EukProt, and MMETSP resources are included post-derep.

Lastly, the predicted protein sequences for the dereplicated genomes (if avaiable) are downloaded. Also here, the user has the option to provide custom files. Protein sequences are then subject to a second dereplication/clustering using [mmseqs2](https://github.com/soedinglab/MMseqs2).

To combine all these genomes within a common taxonomic framework, we use the script gtdb_to_taxdump.py available [here](https://github.com/nick-youngblut/gtdb_to_taxdump).
We modified this script to include 2 additional taxonomic ranks (lineage and kingdom) to better represent the high-level taxonomy of eukaryotes available via the ambiguous rank clade and kingdom on NCBI. 


#### Kraken2 module

As the file size of all genomes after dereplication is still too large to be used as input for building a kraken2 database on most systems (>> 1TB), it is necessary to further subset the available genomes to a more manageable number. This can be done manually or by using the presets offered in the kraken2 module of the workflow. Thes presets include:
* coarse: GTDB and checkV representative genomes, up to 3 (user-configurable) genomes per phylum and family for higher and microbial eukaryotes, respectively. NCBI RefSeq organelle (mitochondria and plastids) are also included. The purpose of this database is to sort metagenomic reads by domain for a second classification step with a higher resolved database. The coarse preset includes the option to remove cross-domain contamination (experimental).
* highres_pro: Dereplicated GTDB and checkV genomes for high resolution classification of prokaryotic reads. Not recommended without prior sorting by domain.
* microeuk: Dereplicated (high resolution) protozoa and fungi with representatives of higher plants and metazoa (3 per phylum), NCBI RefSeq organelles. All user-supplied taxa are included, which are not represented in the afore-mentioned macro-eukaryotic lineages.
* highres_plant :warning: not available yet :warning:
* highres_metazoa :warning: not available yet :warning:
* onestep: Dereplicated GTDB and checkV genomes, plus up to 3 (user-configurable) genomes per phylum and family for higher and microbial eukaryotes, respectively, plus NCBI RefSeq organelles, for high resolution classification of prokaryotic reads in one go (CAUTION: DB size exceeding 350GB).

For building the kraken2 database, the parameters kmer length, minimizer length, and minimizer spaces can be adjusted depending on the type of reads that you want to classify. We recommend a kmer of ~31-35 and ~25 for modern and ancient metagenomic reads, respectively.


#### Kaiju module

All clustered protein sequences are converted into a kaiju database index. No further subsetting is required trim the size of the input data (~130GB).


#### CAT/BAT module

All clustered protein sequences are converted into a diamond database. The taxonomy is formatted to be compatible with the use of CAT/BAT.


### DB composition

![Number of genomes per domain and taxonomic level represented in the database](https://github.com/chassenr/taxdb-integration/blob/master/images/db_stats.jpg)


### Next steps and ToDos:
* Improve conterminator implementation: add maximum runtime limit to avoid endless runs if no/low contamination
* Extend framework to Bracken and KrakenUniq databases


### Running the workflow
After adapting the workflow and cluster files, the workflow can be executed with the following command from within this repository:
```
snakemake --use-conda -j 100 --cluster-config config/cluster.yaml --cluster "sbatch --export=ALL -t {cluster.time} --ntasks-per-node {cluster.ntasks_per_node} --nodes {cluster.nodes} --cpus-per-task {cluster.cpus_per_task} --partition {cluster.partition} --job-name {rulename}.{jobid}"
```
:warning: Snakemake has changed the way the cluster configuration is provided in the command. This workflow has not yet been updated to the latest version of snakemake. It has been tested for snakemake 6.13.1.

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

