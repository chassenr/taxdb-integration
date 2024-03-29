#!/usr/bin/env Rscript

# parse taxonomic path for NCBI genomes

### setting up environment ####

# Check if packages are installed
package.list <- c(
  "crayon",
  "optparse",
  "config",
  "taxonomizr",
  "rentrez",
  "tidyverse",
  "data.table",
  "purrr"
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
    c("-i", "--input"), 
    type = "character", 
    default = NULL,
    help = "NCBI-style assembly summary table", 
    metavar = "character"
  ),
  make_option(
    c("-s", "--sql"), 
    type = "character", 
    default = NULL,
    help = "location and name of sql database  for the NCBI taxdump", 
    metavar = "character"
  ),
  make_option(
    c("-o", "--output"), 
    type = "character", 
    default = NULL,
    help = "GTDB-style formatted taxid to taxonomy map",
    metavar = "character"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input) | is.null(opt$output) | is.null(opt$sql)) {
  print_help(opt_parser)
  stop(
    "You need to provide the assembly summary table, the taxdump sql database, and the name of the output file.\n", 
    call. = FALSE
  )
}


### retrieve taxonomy ####

# read assembly summary
assembly_summary <- fread(
  opt$input,
  h = F,
  sep = "\t",
  quote = ""
)

# map taxid
taxpath <- getTaxonomy(assembly_summary$V7, opt$sql) %>% 
  as_tibble() %>% 
  mutate(accnos = assembly_summary$V1)

# if no taxonomic path found (i.e. deleted or merged taxids that were not updated in the assembly summary file),
# retrieve correct taxid, and repeat getTaxonomy command
if(anyNA(taxpath$superkingdom)) {
  taxpath_new <- assembly_summary %>% 
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
taxpath_parsed <- taxpath_new %>% 
  mutate(
    superkingdom = paste0("d__", superkingdom),
    phylum = paste0("p__", ifelse(is.na(phylum), gsub("d__", "", superkingdom), phylum)),
    class = paste0("c__", ifelse(is.na(class), gsub("p__", "", phylum), class)),
    order = paste0("o__", ifelse(is.na(order), gsub("c__", "", class), order)),
    family = paste0("f__", ifelse(is.na(family), gsub("o__", "", order), family)),
    genus = paste0("g__", ifelse(is.na(genus), gsub("f__", "", family), genus)),
    species = ifelse(
      is.na(species), 
      paste0("s__", gsub("g__", "", genus), " unclassified"),
      paste0("s__", species)
    ),
    species = ifelse(
      !grepl(" ", species),
      paste0(species, " unknown species"),
      species
    )
  ) %>% 
  # mutate( root = "root", .before = superkingdom) %>% 
  unite("path", -accnos, sep = ";") %>% 
  relocate(accnos, .before = path) %>% 
  filter(!is.na(taxpath_new$superkingdom))

# write output table
write_delim(
  taxpath_parsed,
  opt$output,
  delim = "\t",
  col_names = F
)
