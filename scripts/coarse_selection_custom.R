#!/usr/bin/env Rscript

# parse assembly metadata for NCBI genomes

### setting up environment ####

# Check if packages are installed
package.list <- c(
  "crayon",
  "optparse",
  "tidyverse",
  "scales"
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
    c("-t", "--tax"),
    type = "character",
    default = NULL,
    help = "taxonomy table (previously selected taxa)",
    metavar = "character"
  ),
  make_option(
    c("-c", "--custom"),
    type = "character",
    default = NULL,
    help = "custom genome selection (post derep)",
    metavar = "character"
  ),
  make_option(
    c("-r", "--rank"),
    type = "character",
    default = NULL,
    help = "taxonomic rank at which to perform the genome selection",
    metavar = "character"
  ),
  make_option(
    c("-n", "--n_max"),
    type = "integer",
    default = 5,
    help = "maximum number of genomes to select for each level of rank",
    metavar = "number"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = NULL,
    help = "output file name",
    metavar = "character"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$tax) | is.null(opt$custom) |
    is.null(opt$rank) | is.null(opt$output)) {
  print_help(opt_parser)
  stop("All parameters are mandatory.\n", call. = FALSE)
}


### select additional genomes from custom genomes, if their classes (specified rank) aren't already represented in previous selection ####

# possible taxonomic ranks
tax_ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")

# maximum number of genomes per level of selected rank
n_max = opt$n_max # set scaling factor (user defined)

# read previously selected taxonomic paths
# extract and collapse taxonomic paths until selected rank
tax_included <- read.table(opt$tax, h = F, sep = "\t", stringsAsFactors = F) %>%
  separate("V2", into = tax_ranks, sep = ";", remove = F) %>%
  dplyr::rename(accession = V1, path = V2) %>%
  unite("path_trunc", 3:(which(tax_ranks == opt$rank) + 2), sep = ";") %>%
  pull(path_trunc) %>%
  unique()

# if rank is not species, append ; to ensure uniqueness of path
if(opt$rank != "species") {
  tax_included <- gsub("$", ";", tax_included)
}

# read custom genomes metadata table
# extract taxa not yet represented
# calculate genome size (approximated by zipped file size)
custom_genomes <- read.table(opt$custom, h = T, sep = "\t", stringsAsFactors = F) %>%
  filter(!grepl(paste(tax_included, collapse = "|"), taxon)) %>%
  mutate(file_size = file.size(src)) %>%
  # remove strain names as they are not relevant for selection here
  mutate(taxon2 = gsub(";t__.*", "", taxon)) %>%
  separate("taxon2", into = tax_ranks, sep = ";", remove = F)
rownames(custom_genomes) <- custom_genomes$accession

# select genomes
# get largest genome per subtax (refseq if available)
# clean-up this code later...
get_largest <- function(genome_data, rank_select, rank_subtax) {
  max_genome_subtax <- data.frame(unique(genome_data[, c(rank_select, rank_subtax)]))
  max_genome_subtax <- max_genome_subtax[order(max_genome_subtax[, rank_select]), ]
  rownames(max_genome_subtax) <- NULL
  max_genome_subtax$accession <- c()
  max_genome_subtax$genome_size <- c()
  for(i in unique(max_genome_subtax[, rank_select])) {
    tmp <- genome_data[genome_data[, rank_select] == i, ]
    for(j in unique(tmp[, rank_subtax])) {
      tmp.sub <- tmp[tmp[, rank_subtax] == j, ]
      max_genome_subtax[max_genome_subtax[, rank_select] == i & max_genome_subtax[, rank_subtax] == j, "file_size"] <- max(tmp.sub$file_size)
      max_genome_subtax[max_genome_subtax[, rank_select] == i & max_genome_subtax[, rank_subtax] == j, "accession"] <- tmp.sub[which.max(tmp.sub$file_size), "accession"]
    }
  }
  rm(tmp, tmp.sub)
  return(max_genome_subtax)
}

# processing custom genomes

# scenario 1: rank is not species
if(opt$rank != "species") {
  
  rank_select <- opt$rank
  rank_subtax <- tax_ranks[which(tax_ranks == rank_select) + 1]

  # get number of taxa in rank below selected
  div_subtax <- c(by(custom_genomes[, rank_subtax], custom_genomes[, rank_select], function(x) length(unique(x))))
  msg(paste0("Maximum number of ", rank_subtax, " per ", rank_select, ":", max(div_subtax), ".\n"))

  # for now the purpose is to evenly cover diversity at next lower level
  # assuming that any taxonomy is reflecting phylogeny
  # it is more important to capture 2 orders (higher divergence) with few families (less divergence),
  # than multiple families from 1 order
  n_select <- floor(rescale(div_subtax, to = c(1, ifelse(max(div_subtax) > n_max, n_max, max(div_subtax)))))

  # number of taxa at rank
  msg(paste0("Number of taxa at ", rank_select, " level: ", length(div_subtax), ".\n"))
  # number of taxa at subrank
  msg(paste0("Number of taxa at ", rank_subtax, " level: ", sum(div_subtax), ".\n"))
  # number selected
  msg(paste0("Total number of genomes selected: ", sum(n_select), ".\n"))
  # seems to be a reasonable compromise

  # get largest genome per subtax
  max_genome_subtax <- get_largest(custom_genomes, rank_select, rank_subtax)

  # select n subtax at random (?) or from min to max genome size (?)
  # assuming that more closely related taxa should have similar genome length
  # choose by other metrics: number of genes, N50?
  # rewrite metadata rules to get N50 from assembly_stats.txt file (NCBI) rather then with bbmap from genome
  # (avoid genome download this way)
  genome_select <- c()
  tmp <- max_genome_subtax
  for(j in names(n_select)) {
    tmp.sub <- tmp[tmp[, rank_select] == j, ]
    genome_select <- c(genome_select, tmp.sub$accession[order(tmp.sub$file_size, decreasing = T)[seq(1, nrow(tmp.sub), length.out = n_select[j])]])
  }
  genome_select <- custom_genomes[genome_select, c("accession", "taxon")]
  rm(tmp, tmp.sub)
} 

# scenario 2: rank is species
if(opt$rank == "species") {

  rank_select <- opt$rank

  # get number of genomes per species
  div_subtax <- c(by(custom_genomes[, "accession"], custom_genomes[, rank_select], function(x) length(unique(x))))
  msg(paste0("Maximum number of genomes per ", rank_select, ": ", max(div_subtax), ".\n"))
  n_select <- floor(rescale(div_subtax, to = c(1, ifelse(max(div_subtax) > n_max, n_max, max(div_subtax)))))
  msg(paste0("In total ", sum(n_select), " genomes selected.\n"))
  
  # select based on file_size
  genome_select <- c()
  tmp <- custom_genomes
  for(j in names(n_select)) {
    tmp.sub <- tmp[tmp[, rank_select] == j, ]
    genome_select <- c(genome_select, tmp.sub$accession[order(tmp.sub$file_size, decreasing = T)[seq(1, nrow(tmp.sub), length.out = n_select[j])]])
  }
  genome_select <- custom_genomes[genome_select, c("accession", "taxon")]
  rm(tmp, tmp.sub)
}

# write output
write.table(
  genome_select,
  opt$output,
  quote = F,
  row.names = F,
  col.names = F,
  sep = "\t"
)


