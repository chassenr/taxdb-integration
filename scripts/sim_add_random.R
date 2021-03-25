#!/usr/bin/env Rscript

# parse taxonomic path for NCBI genomes
# To do: improve parallelization for speed!

### setting up environment ####

# Check if packages are installed
package.list <- c(
  "crayon",
  "optparse",
  "config",
  "seqinr"
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
    c("-a", "--abundance"), 
    type = "character", 
    default = NULL,
    help = "simulated community proportions by MGSIM", 
    metavar = "character"
  ),
  make_option(
    c("-g", "--genome"), 
    type = "character", 
    default = NULL,
    help = "genome table for simulation", 
    metavar = "character"
  ),
  make_option(
    c("-f", "--fasta"), 
    type = "character", 
    default = NULL,
    help = "fasta file of genome to be used as template to generate random sequence", 
    metavar = "character"
  ),
  make_option(
    c("-A", "--output_abund"), 
    type = "character", 
    default = NULL,
    help = "adjusted community proportions for simulation",
    metavar = "character"
  ),
  make_option(
    c("-G", "--output_genome"), 
    type = "character", 
    default = NULL,
    help = "adjusted genome table for simulation",
    metavar = "character"
  ),
  make_option(
    c("-r", "--random"), 
    type = "character", 
    default = NULL,
    help = "fasta output file for random genome",
    metavar = "character"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$abundance) | is.null(opt$genome) | is.null(opt$fasta) | 
    is.null(opt$output_abund) | is.null(opt$output_genome) | is.null(opt$random)) {
  print_help(opt_parser)
  stop(
    "All parameters are mandatory.\n", 
    call. = FALSE
  )
}


### generate random genome ####

fasta <- read.fasta(opt$fasta, forceDNAtolower = F)
fasta.random <- sapply(
  fasta,
  function(x) {
    sample(x, replace = F)
  }
)
write.fasta(
  fasta.random,
  paste("random", 1:length(fasta.random), sep = "_"),
  file.out = opt$random
)


### adjust community proportions and genome table ####

genome_table <- read.table(
  opt$genome,
  sep = "\t",
  h = T, 
  quote = "",
  stringsAsFactors = F
)

abund_table <- read.table(
  opt$abundance, 
  sep = "\t",
  h = T, 
  quote = "",
  stringsAsFactors = F
)

# append random
genome_table_adj <- rbind(
  genome_table,
  data.frame(
    Accnos = "random",
    Taxon = "d__random;p__random;c__random;o__random;f__random;g__random;s__random",
    Fasta = opt$random
  )
)

# here I am just entering dummy data for Perc_rel_abund and Rank
abund_table_wr <- rbind(
  abund_table,
  data.frame(
    Community = 1,
    Taxon = "d__random;p__random;c__random;o__random;f__random;g__random;s__random",
    Perc_rel_abund = 1,
    Rank = nrow(abund_table) + 1
  )
)

# re-define percentages for each domain
domains <- c(10, 30, 10, 40, 10)
names(domains) <- c("d__random", "d__Eukaryota", "d__Archaea", "d__Bacteria", "d__Viruses")

phylum_list <- list(
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
  )
)
names(phylum_list) <- c("fungi", "plants", "metazoa")
phylum_groups <- c(30, 20, 20, 30)
names(phylum_groups) <- c("fungi", "plants", "metazoa", "protists")

# recalculate percentages for each domain (and phylum group for euks)
abund_table_adj <- vector("list", length = length(domains))
for(i in 1:length(domains)) {
  tmp <- abund_table_wr[grepl(names(domains)[i], abund_table_wr$Taxon), ]
  if(names(domains)[i] == "d__Eukaryota") {
    tmp.list <- vector("list", length = length(phylum_groups))
    for(j in 1:length(phylum_groups)) {
      if(j != length(phylum_groups)) {
        tmp.sub <- tmp[grepl(paste(phylum_list[[j]], collapse = "_|"), tmp$Taxon), ]
      } else {
        tmp.sub <- tmp[!grepl(paste(unlist(phylum_list), collapse = "_|"), tmp$Taxon), ]
      }
      tmp.sub$Perc_rel_abund <- tmp.sub$Perc_rel_abund/sum(tmp.sub$Perc_rel_abund) * phylum_groups[j] * domains[i]/100
      tmp.list[[j]] <- tmp.sub
    }
    tmp <- do.call("rbind", tmp.list)
  } else {
    tmp$Perc_rel_abund <- tmp$Perc_rel_abund/sum(tmp$Perc_rel_abund) * domains[i]
  }
  abund_table_adj[[i]] <- tmp
}
abund_table_adj_df <- do.call("rbind", abund_table_adj)

# update ranks
abund_table_adj_df <- abund_table_adj_df[order(abund_table_adj_df$Perc_rel_abund, decreasing = T), ]
abund_table_adj_df$Rank <- 1:nrow(abund_table_adj_df)

# write output tables
write.table(
  abund_table_adj_df, 
  opt$output_abund,
  quote = F, 
  sep = "\t", 
  row.names = F
)

write.table(
  genome_table_adj, 
  opt$output_genome,
  quote = F, 
  sep = "\t", 
  row.names = F
)
