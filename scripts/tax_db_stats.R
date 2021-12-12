#!/usr/bin/env Rscript

# parse taxonomic path for NCBI genomes

### setting up environment ####

# Check if packages are installed
package.list <- c(
  "crayon",
  "optparse",
  "config",
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
    c("-o", "--output"), 
    type = "character", 
    default = NULL,
    help = "Summary stats of DB content per taxonomic level",
    metavar = "character"
  ),
  make_option(
    c("-c", "--cpus"),
    type = "integer",
    default = 1,
    help = "number of cpus to use [default: 1]",
    metavar = "number"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$taxid) | is.null(opt$output)) {
  print_help(opt_parser)
  stop(
    "All parameters are mandatory.\n", 
    call. = FALSE
  )
}


### Summarize counts of DB entries ####

# read accessions and taxids
acc2taxid <- fread(
  opt$taxid,
  h = T,
  sep = "\t",
  quote = "",
  nThread = opt$cpus
) %>%
  separate(
    col = "path",
    into = c("domain", "lineage", "kingdom", "phylum", "class", "order", "family", "genus", "species"),
    sep = ";",
    remove = FALSE
  ) %>%
  group_by(domain)

if("genome" %in% colnames(acc2taxid)) {
  if("sequence" %in% colnames(acc2taxid)) {
    db_stats <- acc2taxid %>%
      summarise(
        n_lineage = n_distinct(lineage),
	n_kingdom = n_distinct(kingdom),
	n_phylum = n_distinct(phylum),
        n_class = n_distinct(class),
        n_order = n_distinct(order),
        n_family = n_distinct(family),
        n_genus = n_distinct(genus),
        n_species = n_distinct(species),
        n_genomes = n_distinct(genome),
        n_seqs = n()
      )
  } else {
    db_stats <- acc2taxid %>%
      summarise(
        n_lineage = n_distinct(lineage),
        n_kingdom = n_distinct(kingdom),
        n_phylum = n_distinct(phylum),
        n_class = n_distinct(class),
        n_order = n_distinct(order),
        n_family = n_distinct(family),
        n_genus = n_distinct(genus),
        n_species = n_distinct(species),
        n_genomes = n()
      )
  }
} else {
  db_stats <- acc2taxid %>%
    summarise(
      n_lineage = n_distinct(lineage),
      n_kingdom = n_distinct(kingdom),
      n_phylum = n_distinct(phylum),
      n_class = n_distinct(class),
      n_order = n_distinct(order),
      n_family = n_distinct(family),
      n_genus = n_distinct(genus),
      n_species = n_distinct(species),
      n_seqs = n()
    )
}	

# write output table
write_delim(
  db_stats,
  opt$output,
  delim = "\t",
  col_names = T
)
