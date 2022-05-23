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
    c("-n", "--names"), 
    type = "character", 
    default = NULL,
    help = "nucl accession to scientific name map", 
    metavar = "character"
  ),
  make_option(
    c("-r", "--ranks"),
    type = "character",
    default = NULL,
    help = "nucl full taxonomic path and ranks",
    metavar = "character"
  ),
  make_option(
    c("-p", "--protein"),
    type = "character",
    default = NULL,
    help = "protein full taxonomic path and ranks",
    metavar = "character"
  ),
  make_option(
    c("-t", "--organelle"),
    type = "character",
    default = NULL,
    help = "mito or plas",
    metavar = "character"
  ),
  make_option(
    c("-o", "--output_nucl"), 
    type = "character", 
    default = NULL,
    help = "GTDB-style formatted accession to taxonomy map",
    metavar = "character"
  ),
  make_option(
    c("-P", "--output_prot"), 
    type = "character", 
    default = NULL,
    help = "GTDB-style formatted accession to taxonomy map",
    metavar = "character"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$names) | is.null(opt$ranks) | is.null(opt$protein) |
    is.null(opt$organelle) | is.null(opt$output_nucl) | is.null(opt$output_prot)) {
  print_help(opt_parser)
  stop(
    "You need to provide the assembly summary table and the name of the output file.\n", 
    call. = FALSE
  )
}


### parse taxonomic path for organelles (ncbi refseq) ####

# read scientific name to taxpath and ranks map (nucl)
# prune taxonomic path to species
# format taxonomic path to be compatible with GTDB-style framework
nucl_ranks <- read.table(
  opt$ranks,
  h = F,
  sep = "\t",
  quote = "",
  comment.char = "",
  col.names = c("scientific_name", "taxid_name", "path", "taxid_path", "ranks")
) %>% 
  filter(!is.na(path) & !path == "")

# select ranks to keep
keep_ranks <- c("superkingdom", "lineage", "kingdom", "phylum", "class", "order", "family", "genus", "species")

# retain taxonomic names between superkingdom and phylum
tax_keep_ranks <- map_dfr(
  1:nrow(nucl_ranks),
  function(x) {
    tmp_tax <- strsplit(nucl_ranks$path[x], ";", fixed = T)[[1]]
    tmp_rank <- strsplit(nucl_ranks$ranks[x], ";", fixed = T)[[1]]
    # we can assume that first 2 levels are known (cellular organis, eukaryote)
    if(!"Rhodophyta" %in% tmp_tax & !"Viridiplantae" %in% tmp_tax & !"Opisthokonta" %in% tmp_tax) {
      # protozoa: no kingdom, for SAR clade is used twice in a row
      if(tmp_rank[3] == "clade") {
        if(tmp_rank[4] == "clade") {
          tmp_rank[3:4] <- c("lineage", "kingdom")
        } else {
          tmp_rank[3] <- "lineage"
        }
      }
    } else {
      if("Rhodophyta" %in% tmp_tax | "Viridiplantae" %in% tmp_tax) {
        # in plants as of 19.11.2021 only Rhodophyta without kingdom rank, also duplicate kingdom for lineage
        if(!c("kingdom") %in% tmp_rank) {
          tmp_rank <- c(tmp_rank, "kingdom")
          tmp_tax <- c(tmp_tax, "Rhodophyta")
        }
        tmp_rank <- c(tmp_rank, "lineage")
        tmp_tax <- c(tmp_tax, tmp_tax[tmp_rank == "kingdom"])
      }
      if("Opisthokonta" %in% tmp_tax) {
        tmp_rank <- c(tmp_rank, "lineage")
        tmp_tax <- c(tmp_tax, "Opisthokonta")
      }
    }
    tmp <- matrix(c(nucl_ranks$scientific_name[x], tmp_tax[match(keep_ranks, tmp_rank)]), nrow = 1, ncol = length(keep_ranks) + 1)
    colnames(tmp) <- c("scientific_name", keep_ranks)
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
  unite("path", -scientific_name, sep = ";") %>%
  relocate(scientific_name, .before = path) %>%
  filter(!is.na(tax_keep_ranks$superkingdom)) %>% 
  column_to_rownames("scientific_name")

# read accession to scientific name map (nucl)
# append formatted taxonomic path
# retain accession and path
nucl_names <- read.table(
  opt$names,
  h = F,
  sep = "\t",
  col.names = c("accnos", "scientific_name")
)
nucl_names$path <- taxpath_parsed[nucl_names$scientific_name, "path"]
nucl_names <- nucl_names[!is.na(nucl_names$path), ]
nucl_names$accnos_mod <- paste(nucl_names$accnos, opt$organelle, sep = "_")

# read  scientific name to taxpath and ranks map (protein)
# accession for faa files is based on scientific name with non-alphanum replaced by '_'
# prune taxonomic path to species
# format taxonomic path to be compatible with GTDB-style framework
# retain accession and path
prot_ranks <- read.table(
  opt$protein,
  h = F,
  sep = "\t",
  quote = "",
  comment.char = "",
  col.names = c("scientific_name", "taxid_name", "path", "taxid_path", "ranks")
) %>% 
  filter(!is.na(path) & !path == "")

# retain taxonomic names between superkingdom and phylum
tax_keep_ranks <- map_dfr(
  1:nrow(prot_ranks),
  function(x) {
    tmp_tax <- strsplit(prot_ranks$path[x], ";", fixed = T)[[1]]
    tmp_rank <- strsplit(prot_ranks$ranks[x], ";", fixed = T)[[1]]
    # we can assume that first 2 levels are known (cellular organis, eukaryote)
    if(!"Rhodophyta" %in% tmp_tax & !"Viridiplantae" %in% tmp_tax & !"Opisthokonta" %in% tmp_tax) {
      # protozoa: no kingdom, for SAR clade is used twice in a row
      if(tmp_rank[3] == "clade") {
        if(tmp_rank[4] == "clade") {
          tmp_rank[3:4] <- c("lineage", "kingdom")
        } else {
          tmp_rank[3] <- "lineage"
        }
      }
    } else {
      if("Rhodophyta" %in% tmp_tax | "Viridiplantae" %in% tmp_tax) {
        # in plants as of 19.11.2021 only Rhodophyta without kingdom rank, also duplicate kingdom for lineage
        if(!c("kingdom") %in% tmp_rank) {
          tmp_rank <- c(tmp_rank, "kingdom")
          tmp_tax <- c(tmp_tax, "Rhodophyta")
        }
        tmp_rank <- c(tmp_rank, "lineage")
        tmp_tax <- c(tmp_tax, tmp_tax[tmp_rank == "kingdom"])
      }
      if("Opisthokonta" %in% tmp_tax) {
        tmp_rank <- c(tmp_rank, "lineage")
        tmp_tax <- c(tmp_tax, "Opisthokonta")
      }
    }
    tmp <- matrix(c(prot_ranks$scientific_name[x], tmp_tax[match(keep_ranks, tmp_rank)]), nrow = 1, ncol = length(keep_ranks) + 1)
    colnames(tmp) <- c("scientific_name", keep_ranks)
    tmp <- as_tibble(tmp)
  }
)

# parse taxonomic path
# qiime format:
# e.g. d__;k__;p__;c__;o__;f__;g__;s__
# repeat taxon name if NA for intermediate ranks
protein_parsed <- tax_keep_ranks %>%
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
  unite("path", -scientific_name, sep = ";") %>%
  relocate(scientific_name, .before = path) %>%
  filter(!is.na(tax_keep_ranks$superkingdom)) %>% 
  mutate(accnos = paste(gsub("[^[:alnum:]]", "_", scientific_name), opt$organelle, sep = "_"))

# write output
write.table(
  nucl_names[, c(4, 3)],
  opt$output_nucl,
  sep = "\t",
  col.names = F,
  row.names = F,
  quote = F
)

write.table(
  protein_parsed[, c(3, 2)],
  opt$output_prot,
  sep = "\t",
  col.names = F,
  row.names = F,
  quote = F
)

