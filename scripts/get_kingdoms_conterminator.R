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
  "data.table"
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
    help = "taxonomy table", 
    metavar = "character"
  ),
  make_option(
    c("-t", "--taxdump"), 
    type = "character", 
    default = NULL,
    help = "directory containing the kraken2 nodes.dmp and names.dmp", 
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
    help = "taxonomy table with taxid appended",
    metavar = "character"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input) | is.null(opt$output) | is.null(opt$taxdump) | is.null(opt$sql)) {
  print_help(opt_parser)
  stop(
    "All parameters are mandatory.\n", 
    call. = FALSE
  )
}

### retrieve taxid for each path in dereplicated database ####

# format NCBI taxdump database
read.names.sql(
  paste0(opt$taxdump, "/names.dmp"),
  opt$sql
)
read.nodes.sql(
  paste0(opt$taxdump, "/nodes.dmp"),
  opt$sql
)

# read taxonomy table
tax_table <- fread(
  opt$input,
  h = F,
  sep = "\t",
  quote = ""
)

# map taxid
phylum_table <- separate(tax_table, V2, into = c("domain", "phylum", "class", "order", "family", "genus", "species"), sep = ";") %>% 
  filter(domain == "d__Eukaryota") %>% 
  select(2, 3) %>% 
  unique()

# select the phylum groupings for:
phylum_list <- list(
  # higher plants
  c("p__Streptophyta"),
  # metazoa
  c(
    "p__Acanthocephala",
    "p__Annelida",
    "p__Arthropoda",
    "p__Brachiopoda",
    "p__Bryozoa",
    "p__Chordata",
    "p__Cnidaria",
    "p__Ctenophora",
    "p__Dicyemida",
    "p__Echinodermata",
    "p__Hemichordata",
    "p__Mollusca",
    "p__Nematoda",
    "p__Nemertea",
    "p__Onychophora",
    "p__Orthonectida",
    "p__Phoronida",
    "p__Placozoa",
    "p__Platyhelminthes",
    "p__Porifera",
    "p__Priapulida",
    "p__Rotifera",
    "p__Tardigrada",
    "p__Xenacoelomorpha"
  ),
  # fungi
  c(
    "p__Ascomycota",
    "p__Basidiomycota",
    "p__Blastocladiomycota",
    "p__Chytridiomycota",
    "p__Cryptomycota",
    "p__Microsporidia",
    "p__Mucoromycota",
    "p__Zoopagomycota"
  )
)
names(phylum_list) <- c("plants", "metazoa", "fungi")
#   protist and algae: to be selected based on exclusion

# retrieve taxid
taxid_list <- map(
  1:length(phylum_list),
  function(X) {
    tmp <- getId(gsub("^p__", "", phylum_list[[X]]), opt$sql) %>% 
      strsplit(., split = ",") %>% 
      unlist() %>% 
      getTaxonomy(., opt$sql)
    gsub(" ", "", rownames(tmp)[is.na(tmp[, "class"])])
  }
)
names(taxid_list) <- names(phylum_list)
all.equal(sapply(phylum_list, length), sapply(taxid_list, length))

# parse string for conterminator
# taxid 2 and 3 will always be archaea and bacteria
# euks are taxid 4
conterminator_string <- paste0(
  "'(2||3),(",
  paste(taxid_list$fungi, collapse = "||"),
  "),(",
  taxid_list$plants,
  "),(",
  paste(taxid_list$metazoa, collapse = "||"),
  "),(4&&!",
  paste(unlist(taxid_list), collapse = "&&!"),
  ")'"
)

# write output table
write.table(
  conterminator_string,
  opt$output,
  sep = "\t",
  col.names = F,
  row.names = F,
  quote = F
)
