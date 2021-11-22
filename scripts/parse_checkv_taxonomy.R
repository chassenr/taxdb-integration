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
    c("-t", "--tax_genbank"), 
    type = "character", 
    default = NULL,
    help = "ncbi taxonomic paths for genbank genomes", 
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
    is.null(opt$meta) | is.null(opt$tax_genbank)) {
  print_help(opt_parser)
  stop(
    "All parameters are mandatory.\n", 
    call. = FALSE
  )
}


### read checkv metadata tables ####

# IMG contigs
img <- fread(
  opt$img,
  h = T,
  sep = "\t",
  quote = ""
) %>% 
  mutate(path = "d__Viruses;l__Viruses;k__Viruses;p__Viruses;c__Viruses;o__Viruses;f__Viruses;g__Viruses;s__Viruses unclassified")

# genbank contigs
genbank <- fread(
  opt$genbank,
  h = T,
  sep = "\t",
  quote = ""
)

# taxonomy of genbank viruses
genbank_tax <- fread(
  opt$tax_genbank,
  h = F,
  sep = "\t",
  quote = "",
  col.names = c("checkv_id", "ncbi_id", "path", "ranks")
) %>%
  as_tibble()

# select ranks to keep
keep_ranks <- c("superkingdom", "lineage", "kingdom", "phylum", "class", "order", "family", "genus", "species")

# retain taxonomic names between superkingdom and phylum
tax_keep_ranks <- map_dfr(
  1:nrow(genbank_tax),
  function(x) {
    tmp_tax <- strsplit(genbank_tax$path[x], ";", fixed = T)[[1]]
    tmp_rank <- strsplit(genbank_tax$ranks[x], ";", fixed = T)[[1]]
    # we can assume that the path always starts with superkingdom
    if(tmp_rank[2] == "clade") {
      tmp_rank[2] <- "lineage"
    }
    tmp <- matrix(c(genbank_tax$checkv_id[x], tmp_tax[match(keep_ranks, tmp_rank)]), nrow = 1, ncol = length(keep_ranks) + 1)
    colnames(tmp) <- c("checkv_id", keep_ranks)
    tmp <- as_tibble(tmp)
  }
)

# parse taxonomic path
# qiime format:
# e.g. d__;k__;p__;c__;o__;f__;g__;s__
# repeat taxon name if na for intermediate ranks
genbank_parsed <- tax_keep_ranks %>% 
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
  unite("path", -checkv_id, sep = ";") %>% 
  relocate(checkv_id, .before = path) %>% 
  filter(!is.na(tax_keep_ranks$superkingdom)) %>% 
  right_join(genbank, ., by = c("checkv_id" = "checkv_id"))

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
  mutate(path_new = paste0(path, "___", rep_genome))

# map taxonomy of cluster reps to all genomes
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
