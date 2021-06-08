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
    c("-t", "--taxid"), 
    type = "character", 
    default = NULL,
    help = "list of accession numbers and taxids in DB (if 3 columns: 1st genome, 2nd sequence, 3rd taxid)", 
    metavar = "character"
  ),
  make_option(
    c("-s", "--sql"), 
    type = "character", 
    default = NULL,
    help = "location and name of sql database for the DB taxdump", 
    metavar = "character"
  ),
  make_option(
    c("-o", "--output"), 
    type = "character", 
    default = NULL,
    help = "Summary stats of DB content per taxonomic level",
    metavar = "character"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$taxid) | is.null(opt$output) | is.null(opt$sql)) {
  print_help(opt_parser)
  stop(
    "You need to provide the assembly summary table, the taxdump sql database, and the name of the output file.\n", 
    call. = FALSE
  )
}


### Summarize counts of DB entries ####

# read accessions and taxids
acc2taxid <- fread(
  opt$taxid,
  h = F,
  sep = "\t",
  quote = ""
)
if(ncol(acc2taxid) == 2) {
  colnames(acc2taxid) <- c("accnos", "taxid")
} else {
  colnames(acc2taxid) <- c("accnos", "seq", "taxid")
}

# map taxid
taxpath <- getTaxonomy(acc2taxid$taxid, opt$sql) %>% 
  as_tibble() %>%
  mutate(genome = acc2taxid$accnos) %>%
  group_by(superkingdom)

# calculate genome stats
db_stats <- taxpath %>%
  summarise(
    n_phylum = n_distinct(phylum),
    n_class = n_distinct(class),
    n_order = n_distinct(order),
    n_family = n_distinct(family),
    n_genus = n_distinct(genus),
    n_species = n_distinct(species),
    n_genomes = n_distinct(genome),
    n_seqs = n()
  )

# write output table
write_delim(
  db_stats,
  opt$output,
  delim = "\t",
  col_names = T
)
