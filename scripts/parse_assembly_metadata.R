#!/usr/bin/env Rscript

# parse assembly metadata for NCBI genomes

### setting up environment ####

# Check if packages are installed
package.list <- c(
  "crayon",
  "optparse",
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
    c("-g", "--gene_count"),
    type = "character", 
    default = NULL, 
    help = "gene counts extracted from feature count table",
    metavar = "character"
  ),
  make_option(
    c("-c", "--contig_stats"),
    type = "character",
    default = NULL, 
    help = "genome length extracted from assembly stats table", 
    metavar = "character"
  ),
  make_option(
    c("-s", "--assembly_summary"),
    type = "character", 
    default = NULL, 
    help = "assembly summary table",
    metavar = "character"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = NULL, 
    help = "Name of output file with combined metadata", 
    metavar = "character"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$gene_count) | is.null(opt$contig_stats) |
    is.null(opt$assembly_summary) | is.null(opt$output)) {
  print_help(opt_parser)
  stop("All parameters are mandatory.\n", call. = FALSE)
}


### get assembly metadata ####
msg("Reading assembly summary...\n")
# only keep accnos and info about refseq genome (refseq_category)
asm_summary <- fread(
  opt$assembly_summary, 
  h = F, 
  sep = "\t"
) %>%
  as_tibble() %>% 
  select(1, 5) %>% 
  setNames(c("accession", "refseq_category")) 

msg("Reading assembly stats...\n")
# only using genome length right now (N50 not provided for all genomes)
contig_stats <- fread(
  opt$contig_stats, 
  h = F, 
  sep = "\t"
) %>%
  as_tibble() %>% 
  select(1, 7) %>% 
  rename(accession = V1, genome_size = V7)
  
msg("Reading gene counts...\n")
# keep gene counts for all (incl. non-nuclear if present), otherwise Primary accession
# this selection is not considering taxon-specific terminology (e.g. macronucleus) or misleading descriptions in the feature_counts.txt file
gene_counts <- fread(
  opt$gene_count, 
  h = F, 
  sep = "\t"
) %>%
  as_tibble() %>% 
  select(3, 5, 7) %>% 
  setNames(c("accession", "asm_unit_name", "protein_count"))
gene_counts_all <- gene_counts %>% 
  filter(asm_unit_name == "all")
gene_counts_primary <- gene_counts %>% 
  filter(
    !accession %in% gene_counts_all$accession,
    asm_unit_name == "Primary Assembly"
    )
gene_counts_nr <- bind_rows(gene_counts_all, gene_counts_primary)
msg(paste0("Found gene counts for ", nrow(gene_counts_nr), " genomes\n"))


### merge assembly metadata ####
asm_metadata <- left_join(asm_summary, contig_stats, by = "accession") %>% 
  left_join(., gene_counts_nr, by = "accession")


### write final output ####
msg("writing output...")
write_tsv(
  asm_metadata,
  opt$output,
  col_names = TRUE
)
cat(" done\n")

