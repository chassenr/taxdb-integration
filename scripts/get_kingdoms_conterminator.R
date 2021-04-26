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
    c("-t", "--taxdump"), 
    type = "character", 
    default = NULL,
    help = "directory containing the kraken2 nodes.dmp and names.dmp", 
    metavar = "character"
  ),
  make_option(
    c("-k", "--kingdoms"), 
    type = "character", 
    default = "prokaryotes,fungi,protists,plants,metazoa",
    help = "space-separated list of kingdoms to compare", 
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
    help = "conterminator string for kingdoms",
    metavar = "character"
  ),
  make_option(
    c("-x", "--blacklist"), 
    type = "character", 
    default = NULL,
    help = "conterminator string for blacklist (always viruses)", 
    metavar = "character"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$output) | is.null(opt$taxdump) | is.null(opt$sql)) {
  print_help(opt_parser)
  stop(
    "All parameters are mandatory.\n", 
    call. = FALSE
  )
}

### retrieve taxid for each path in dereplicated database ####

# format NCBI taxdump database
if(!file.exists(opt$sql)) {
  read.names.sql(
    paste0(opt$taxdump, "/names.dmp"),
    opt$sql
  )
  read.nodes.sql(
    paste0(opt$taxdump, "/nodes.dmp"),
    opt$sql
  )
}

# select the phylum groupings for:
phylum_list <- list(
  # bacteria, archaea (taken from domain list)
  c(),
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
  ),
  # protist: to be selected based on exclusion
  c(),
  # Viridiplantae
  c(
    "p__Streptophyta",
    "p__Chlorophyta",
    "p__Rhodophyta"
  ),
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
  )
)
names(phylum_list) <- c("prokaryotes", "fungi", "protists", "plants", "metazoa")

# domain list
domains <- c("d__Archaea", "d__Bacteria", "d__Eukaryota", "d__Viruses")
tmp <- getId(domains, opt$sql) %>% 
  strsplit(., split = ",") %>% 
  unlist() %>% 
  getTaxonomy(., opt$sql)
domains_taxid <- gsub(" ", "", rownames(tmp)[is.na(tmp[, "phylum"])])
names(domains_taxid) <- domains

# subset to kingdoms to compare
kingdoms <- strsplit(opt$kingdoms, ",")[[1]]
phylum_list <- phylum_list[kingdoms]

# retrieve taxid
taxid_list <- map(
  1:length(phylum_list),
  function(X) {
    if(!is.null(phylum_list[[X]])) {
      tmp <- getId(phylum_list[[X]], opt$sql) %>% 
        strsplit(., split = ",") %>% 
        unlist() %>% 
        getTaxonomy(., opt$sql)
      gsub(" ", "", rownames(tmp)[is.na(tmp[, "class"]) & !is.na(tmp[, "phylum"])]) 
    } else {
      c()
    }
  }
)
names(taxid_list) <- names(phylum_list)

# parse string for conterminator
conterminator_string <- gsub(
  ")(", 
  "),(",
  paste0(
    "'",
    ifelse(c("prokaryotes") %in% kingdoms, paste0("(",paste(domains_taxid[1:2], collapse = "||"), ")"), ""),
    ifelse(c("fungi") %in% kingdoms, paste0("(",paste(taxid_list$fungi, collapse = "||"), ")"), ""),
    ifelse(c("plants") %in% kingdoms, paste0("(", paste(taxid_list$plants, collapse = "||"), ")"), ""),
    ifelse(c("metazoa") %in% kingdoms, paste0("(", paste(taxid_list$metazoa, collapse = "||"), ")"), ""),
    ifelse(c("protists") %in% kingdoms, paste0("(", domains_taxid[3], "&&!", paste(unlist(taxid_list), collapse = "&&!"), ")"), ""),
    "'"
  ),
  fixed = T
)

# blacklist viruses and p__Eukaryota
tmp <- getId("p__Eukaryota", opt$sql) %>%
  strsplit(., split = ",") %>%
  unlist() %>%
  getTaxonomy(., opt$sql)
blacklist_string <- paste0(
  "'", 
  domains_taxid[4],
  ",",
  gsub(" ", "", rownames(tmp)[is.na(tmp[, "class"]) & !is.na(tmp[, "phylum"])]),
  "'"
)
rm(tmp)

# write output
write.table(
  conterminator_string,
  opt$output,
  sep = "\t",
  col.names = F,
  row.names = F,
  quote = F
)
write.table(
  blacklist_string,
  opt$blacklist,
  sep = "\t",
  col.names = F,
  row.names = F,
  quote = F
)
