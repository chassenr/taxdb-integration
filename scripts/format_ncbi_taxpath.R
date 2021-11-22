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
    c("-i", "--input"), 
    type = "character", 
    default = NULL,
    help = "NCBI-style assembly summary table", 
    metavar = "character"
  ),
  make_option(
    c("-t", "--taxon"),
    type = "character",
    default = NULL,
    help = "ncbi taxonomic group (library_name)",
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

if (is.null(opt$input) | is.null(opt$taxon) | is.null(opt$output)) {
  print_help(opt_parser)
  stop(
    "You need to provide the assembly summary table and the name of the output file.\n", 
    call. = FALSE
  )
}


### retrieve taxonomy ####

# read taxonkit paths and ranks
tax_all_ranks <- fread(
  opt$input,
  h = F,
  sep = "\t",
  quote = "",
  col.names = c("accnos", "taxid", "path", "ranks")
) %>%
  as_tibble() %>%
  # remove genomes without taxonomic paths
  filter(!(is.na(path) | path == ""))

# select ranks to keep
keep_ranks <- c("superkingdom", "lineage", "kingdom", "phylum", "class", "order", "family", "genus", "species")

# retain taxonomic names between superkingdom and phylum
tax_keep_ranks <- map_dfr(
  1:nrow(tax_all_ranks),
  function(x) {
    tmp_tax <- strsplit(tax_all_ranks$path[x], ";", fixed = T)[[1]]
    tmp_rank <- strsplit(tax_all_ranks$ranks[x], ";", fixed = T)[[1]]
    # we can assume that first 2 levels are known (cellular organis, eukaryote)
    # protozoa: no kingdom, for SAR clade is used twice in a row
    if(opt$taxon == "protozoa") {
      if(tmp_rank[3] == "clade") {
        if(tmp_rank[4] == "clade") {
          tmp_rank[3:4] <- c("lineage", "kingdom")
        } else {
          tmp_rank[3] <- "lineage"
        }
      }
    }
    # in plants as of 19.11.2021 only Rhodophyta without kingdom rank, also duplicate kingdom for lineage
    if(opt$taxon == "plant") {
      if(!c("kingdom") %in% tmp_rank) {
        tmp_rank <- c(tmp_rank, "kingdom")
        tmp_tax <- c(tmp_tax, "Rhodophyta")
      }
      tmp_rank <- c(tmp_rank, "lineage")
      tmp_tax <- c(tmp_tax, tmp_tax[tmp_rank == "kingdom"])
    }
    # all invertebrates, fungi and vertebrates are opithokonts
    if(opt$taxon %in% c("invertebrate", "fungi", "vertebrate_other", "vertebrate_mammalian")) {
      tmp_rank <- c(tmp_rank, "lineage")
      tmp_tax <- c(tmp_tax, "Opisthokonta")
    }
    tmp <- matrix(c(tax_all_ranks$accnos[x], tmp_tax[match(keep_ranks, tmp_rank)]), nrow = 1, ncol = length(keep_ranks) + 1)
    colnames(tmp) <- c("accnos", keep_ranks)
    tmp <- as_tibble(tmp)
  }
)

# parse taxonomic path
# qiime format:
# e.g. d__;k__;p__;c__;o__;f__;g__;s__
# repeat taxon name if NA for intermediate ranks
taxpath_parsed <- tax_keep_ranks %>%
  mutate(
    superkingdom = paste0("d__", superkingdom),
    lineage = paste0("l__", ifelse(is.na(lineage), gsub("d__", "", superkingdom), lineage)),
    kingdom = paste0("k__", ifelse(is.na(kingdom), gsub("l__", "", lineage), kingdom)),
    phylum = paste0("p__", ifelse(is.na(phylum), gsub("k__", "", kingdom), phylum)),
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
  unite("path", -accnos, sep = ";") %>%
  relocate(accnos, .before = path) %>%
  filter(!is.na(tax_keep_ranks$superkingdom))

# write output table
write_delim(
  taxpath_parsed,
  opt$output,
  delim = "\t",
  col_names = F
)
