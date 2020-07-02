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
    help = "GTDB-style formatted taxid to taxonomy map",
    metavar = "character"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
if (is.null(opt$input) | is.null(opt$output) | is.null(opt$taxdump) | is.null(opt$sql)) {
  print_help(opt_parser)
  stop(
    "You need to provide the assembly summary table, the location of the taxdump files and sql database, and the name of the output file.\n", 
    call. = FALSE
  )
}


### retrieve taxonomy ####

# format NCBI taxdump database
read.names.sql(
  paste0(opt$taxdump, "/names.dmp"),
  opt$sql
)
read.nodes.sql(
  paste0(opt$taxdump, "/nodes.dmp"),
  opt$sql
)

# read assembly summary
assembly_summary <- fread(
  opt$input,
  h = F,
  sep = "\t",
  quote = ""
)

# map taxid
taxpath <- getTaxonomy(assembly_summary$V7, opt$sql) %>% 
  as_tibble()

# parse taxonomic path
# qiime format:
# e.g. d__Archaea;p__Halobacterota;c__Methanosarcinia;o__Methanosarcinales;f__Methanosarcinaceae;g__Methanosarcina;s__Methanosarcina mazei
# repeat taxon name if na for intermediate ranks
taxpath_parsed <- taxpath %>% 
  mutate(
    superkingdom = paste0("d__", superkingdom),
    phylum = paste0("p__", ifelse(is.na(phylum), gsub("d__", "", superkingdom), phylum)),
    class = paste0("c__", ifelse(is.na(class), gsub("p__", "", phylum), class)),
    order = paste0("o__", ifelse(is.na(order), gsub("c__", "", class), order)),
    family = paste0("f__", ifelse(is.na(family), gsub("o__", "", order), family)),
    genus = paste0("g__", ifelse(is.na(genus), gsub("f__", "", family), genus)),
    species = paste0("s__", species)
  )

# fix duplicate tax levels
# adapt for tidyverse later
taxpath_df <- as.data.frame(taxpath_parsed)
taxpath_fixed <- taxpath_df
for(j in 2:ncol(taxpath_df)) {
  for(i in unique(taxpath_df[, j])) {
    tmp1 <- taxpath_df[taxpath_df[, j] == i, j-1]
    if(length(unique(tmp1)) > 1) {
      tmp2 <- as.factor(tmp1)
      levels(tmp2) <- LETTERS[1:length(unique(tmp1))]
      taxpath_fixed[taxpath_df[, j] == i, j] <- paste0(taxpath_df[taxpath_df[, j] == i, j], "_", as.character(tmp2))
    }
  }
} 

# remove unknown domain and format path
taxpath_out <- taxpath_fixed %>% 
  as_tibble() %>% 
  # mutate( root = "root", .before = superkingdom) %>% 
  unite("path", sep = ";") %>% 
  mutate(
    accnos = assembly_summary$V1,
    .before = path
  ) %>% 
  filter(!is.na(taxpath$superkingdom))

# write output table
write_delim(
  taxpath_out,
  opt$output,
  delim = "\t",
  col_names = F
)
