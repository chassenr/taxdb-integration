# parse taxonomic path for NCBI genomes

# conda install -c bioconda r-taxonomizr
require(taxonomizr)
require(data.table)
require(tidyverse)
require(purrr)

# format NCBI taxdump
read.names.sql(
  '/home/chh/Documents/Projects/NCBI_taxdb_integration/taxdump/names.dmp',
  '/home/chh/Documents/Projects/NCBI_taxdb_integration/taxdump/accessionTaxa.sql'
)
read.nodes.sql(
  '/home/chh/Documents/Projects/NCBI_taxdb_integration/taxdump/nodes.dmp',
  '/home/chh/Documents/Projects/NCBI_taxdb_integration/taxdump/accessionTaxa.sql'
)

# map taxid
assembly_summary <- fread(
  "/home/chh/Documents/Projects/NCBI_taxdb_integration/fungi/assembly_summary_subset.txt",
  h = F,
  sep = "\t"
)
taxpath <- getTaxonomy(assembly_summary$V7, '/home/chh/Documents/Projects/NCBI_taxdb_integration/taxdump/accessionTaxa.sql') %>% 
  as_tibble()

# parse taxonomic path
# qiime format:
# d__Archaea;p__Halobacterota;c__Methanosarcinia;o__Methanosarcinales;f__Methanosarcinaceae;g__Methanosarcina;s__Methanosarcina mazei
# repeat taxon name if na
taxpath_parsed <- taxpath %>% 
  mutate(
    superkingdom = paste0("k__", superkingdom),
    phylum = paste0("p__", ifelse(is.na(phylum), gsub("k__", "", superkingdom), phylum)),
    class = paste0("c__", ifelse(is.na(class), gsub("p__", "", phylum), class)),
    order = paste0("o__", ifelse(is.na(order), gsub("c__", "", class), order)),
    family = paste0("f__", ifelse(is.na(family), gsub("o__", "", order), family)),
    genus = paste0("g__", genus),
    species = paste0("s__", species)
  ) %>% 
  unite("path", sep = ";") %>% 
  mutate(
    accnos = assembly_summary$V1,
    .before = path
  )

# write output table
write_delim(
  taxpath_parsed,
  "/home/chh/Documents/Projects/NCBI_taxdb_integration/fungi/fungi_taxonomy.txt",
  delim = "\t",
  col_names = F
)
