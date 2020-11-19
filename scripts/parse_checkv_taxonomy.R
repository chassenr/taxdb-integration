#!/usr/bin/env Rscript

# parse taxonomic path for NCBI genomes

### setting up environment ####

# Check if packages are installed
package.list <- c(
  "crayon",
  "optparse",
  "config",
  "taxonomizr",
  "tidyverse",
  "data.table",
  "purrr",
  "taxonomizr",
  "rentrez"
)

# Function to check if packages are installed
is.installed <- function(pkg){
  is.element(pkg, installed.packages()[,1])
}

# If not all packages are available
if(any(!is.installed(package.list))) {
  cat("Not all required packages are available. They will now be installed.\n")
  
  # give the user the chance to abort manually
  Sys.sleep(20)
  
  # then install packages
  for(i in which(!is.installed(package.list))) {
    suppressMessages(install.packages(package.list[i], repos = "http://cran.us.r-project.org"))
  }
}

# Break the script if the package installation was unsuccessful
if(any(!is.installed(package.list))) {
  cat(
    paste0(
      "Unable to install all required packages.\nPlease install ",
      paste0(package.list[!is.installed(package.list)], collapse = ", "),
      " manually."
    )
  )
  break
}

# Load packages
cat("Loading libraries...")
silent <- suppressMessages(lapply(package.list, function(X) {require(X, character.only = TRUE)}))
rm(silent)
cat(" done\n")

# Some functions for message output
msg <- function(X){
  cat(crayon::white(paste0("[",format(Sys.time(), "%T"), "]")), X)
}
msg_sub <- function(X){
  cat(crayon::white(paste0("  [",format(Sys.time(), "%T"), "]")), X)
}


### Reading command line options ####

# define command line options
option_list <- list(
  make_option(
    c("-g", "--genbank"), 
    type = "character", 
    default = NULL,
    help = "genbank metadata", 
    metavar = "character"
  ),
  make_option(
    c("-i", "--img"), 
    type = "character", 
    default = NULL,
    help = "metadata for IMG viral contigs", 
    metavar = "character"
  ),
  make_option(
    c("-c", "--clusters"), 
    type = "character", 
    default = NULL,
    help = "checkV cluster mapping file", 
    metavar = "character"
  ),
  make_option(
    c("-t", "--taxdump"), 
    type = "character", 
    default = NULL,
    help = "directory containing NCBI nodes.dmp and names.dmp", 
    metavar = "character"
  ),
  make_option(
    c("-s", "--sql"), 
    type = "character", 
    default = NULL,
    help = "location and name of sql database that will be generated from the taxdump", 
    metavar = "character"
  ),
  make_option(
    c("-o", "--output"), 
    type = "character", 
    default = NULL,
    help = "gtdb-style taxonomy table",
    metavar = "character"
  ),
  make_option(
    c("-m", "--metadata"), 
    type = "character", 
    default = NULL,
    help = "additional metadata for checkV cluster representatives",
    metavar = "character"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$img) | is.null(opt$genbank) | is.null(opt$clusters) | is.null(opt$output) | 
    is.null(opt$meta) | is.null(opt$sql) | is.null(opt$taxdump)) {
  print_help(opt_parser)
  stop(
    "All parameters are mandatory.\n", 
    call. = FALSE
  )
}


### read NCBI taxonomy info ####

# format NCBI taxdump database
read.names.sql(
  paste0(opt$taxdump, "/names.dmp"),
  opt$sql
)
read.nodes.sql(
  paste0(opt$taxdump, "/nodes.dmp"),
  opt$sql
)


### read checkv metadata tables ####

# IMG contigs
img <- fread(
  opt$img,
  h = T,
  sep = "\t",
  quote = ""
) %>% 
  mutate(path = "d__Viruses;p__Viruses;c__Viruses;o__Viruses;f__Viruses;g__Viruses;s__Viruses")

# genbank contigs
# as lineage provided by checkV does not contain rank information
# retrieve info using taxonomizr
genbank <- fread(
  opt$genbank,
  h = T,
  sep = "\t",
  quote = ""
)

# map taxid
taxpath <- getTaxonomy(genbank$ncbi_id, opt$sql) %>% 
  as_tibble() %>% 
  mutate(accnos = genbank$checkv_id)

# if no taxonomic path found (i.e. deleted or merged taxids that were not updated in the assembly summary file),
# retrieve correct taxid, and repeat getTaxonomy command
if(anyNA(taxpath$superkingdom)) {
  taxpath_new <- genbank %>% 
    filter(is.na(taxpath$superkingdom)) %>% 
    pull(1) %>% 
    map_dfr(., function(X) {
      out_search <- entrez_search(db = "assembly", term = X)
      out_links <- entrez_link(dbfrom = "assembly", id = out_search$ids[1], db = "all")
      taxid_new <- out_links$links$assembly_taxonomy[1]
      getTaxonomy(taxid_new, opt$sql) %>% 
        as_tibble() %>% 
        mutate(accnos = X)
    }) %>% 
    rows_update(taxpath, ., by = "accnos", copy = T)
} else {
  taxpath_new <- taxpath
}

# parse taxonomic path
# qiime format:
# e.g. d__Archaea;p__Halobacterota;c__Methanosarcinia;o__Methanosarcinales;f__Methanosarcinaceae;g__Methanosarcina;s__Methanosarcina mazei
# repeat taxon name if na for intermediate ranks
genbank_parsed <- taxpath_new %>% 
  mutate(
    superkingdom = paste0("d__", superkingdom),
    phylum = paste0("p__", ifelse(is.na(phylum), gsub("d__", "", superkingdom), phylum)),
    class = paste0("c__", ifelse(is.na(class), gsub("p__", "", phylum), class)),
    order = paste0("o__", ifelse(is.na(order), gsub("c__", "", class), order)),
    family = paste0("f__", ifelse(is.na(family), gsub("o__", "", order), family)),
    genus = paste0("g__", ifelse(is.na(genus), gsub("f__", "", family), genus)),
    species = ifelse(
      is.na(species), 
      paste0("s__", gsub("g__", "", genus)),
      paste0("s__", species)
    )
  ) %>% 
  # mutate( root = "root", .before = superkingdom) %>% 
  unite("path", -accnos, sep = ";") %>% 
  relocate(accnos, .before = path) %>% 
  filter(!is.na(taxpath_new$superkingdom)) %>% 
  right_join(genbank, ., by = c("checkv_id" = "accnos"))

# combine taxonomy for img and genbank genomes
full_tax <- bind_rows(
  img %>% select(checkv_id, path),
  genbank_parsed %>% select(checkv_id, path)
)

# checkV clusters
clusters <- fread(
  opt$clusters,
  h = T,
  sep = "\t",
  quote = ""
) %>% 
  mutate(tmp_id = ifelse(genbank_rep == "NULL", circular_rep, genbank_rep)) %>% 
  left_join(., full_tax, by = c("tmp_id" = "checkv_id")) %>% 
  mutate(path_new = paste(path, rep_genome, sep = "_")) %>% 
  select(-path, -tmp_id)
# choose better naming strategy for clusters?


### parse taxonomy table ####

full_tax_parsed <- data.frame(
  accnos = unlist(strsplit(clusters$genome_ids, ",", fixed = T)),
  path = rep(clusters$path_new, sapply(strsplit(clusters$genome_ids, ",", fixed = T), length)),
  stringsAsFactors = F
)

# write output tables
write_delim(
  clusters,
  opt$meta,
  delim = "\t",
  col_names = T
)

write_delim(
  full_tax_parsed,
  opt$output,
  delim = "\t",
  col_names = F
)