#!/usr/bin/env Rscript

# parse taxonomic path for NCBI genomes

### setting up environment ####

# Check if packages are installed
package.list <- c(
  "crayon",
  "optparse",
  "config",
  "taxonomizr",
  "tidytable",
  "tidyverse",
  "data.table",
  "caret",
  "vegan"
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
    c("-i", "--kraken"), 
    type = "character", 
    default = NULL,
    help = "kraken output after conifer", 
    metavar = "character"
  ),
  make_option(
    c("-t", "--taxdump"), 
    type = "character", 
    default = NULL,
    help = "directory containing nodes.dmp and names.dmp of the taxonomy used for the classification", 
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
    c("-r", "--rtl"), 
    type = "double", 
    default = 0,
    help = "RTL cut-off to be applied [default: 0]", 
    metavar = "number"
  ),
  make_option(
    c("-m", "--mode"), 
    type = "character", 
    default = NULL,
    help = "Are reads simulated as paired or single end (options: PE or SE)", 
    metavar = "character"
  ),
  make_option(
    c("-o", "--output"), 
    type = "character", 
    default = NULL,
    help = "base name for output files",
    metavar = "character"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$kraken) | is.null(opt$output) | is.null(opt$taxdump) | is.null(opt$sql) | is.null(opt$mode)) {
  print_help(opt_parser)
  stop(
    "You need to provide the assembly summary table, the location of the taxdump files and sql database, and the name of the output file.\n", 
    call. = FALSE
  )
}


### retrieve taxonomy ####

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


### read kraken conifer output ####
select_cols <- c(2:3, if(opt$mode == "PE") c(8, 11) else c(6, 7))
kraken_out <- fread(
  opt$kraken, 
  header = FALSE,
  fill = TRUE, 
  sep = "\t"
)[, ..select_cols]
colnames(kraken_out) <- c("tax_path_in", "taxid_out", "conf", "rtl")

# assign taxonomy to output
# convert any non alpha-numeric character to underscore
# (this is done for input taxonomic path by simulation program)
tax_out_class <- getTaxonomy(
  kraken_out$taxid_out,
  opt$sql
) %>% 
  gsub("[ ;]", "_", .)

# recode NA as character
tax_out_class[is.na(tax_out_class)] <- "NA"

# splitting input taxonomic path
tax_in_class <- kraken_out$tax_path_in %>% 
  gsub("__SEQ[0-9]*$", "", .) %>%
  gsub("^d__", "", .) %>%
  as_tibble() %>%
  separate(
    col = "value",
    into = c("superkingdom", "phylum", "class", "order", "family", "genus", "species"),
    sep = "_[pcofgs]__"
  ) %>% 
  as.matrix()
# set random to NA
tax_in_class[tax_in_class == "random"] <- "NA"

# apply RTL threshold
tax_out_class[kraken_out$rtl < opt$rtl, ] <- "NA"

# directly compare tax assignments
tax_comparison <- tax_in_class == tax_out_class


### calculate confusion matrix ####

# comparison between domains
cm_domain <- confusionMatrix(
  data = factor(data.frame(tax_out_class)[, 1], levels = c("Eukaryota", "Bacteria", "Archaea", "Viruses", "NA")),
  reference = factor(data.frame(tax_in_class)[, 1], levels = c("Eukaryota", "Bacteria", "Archaea", "Viruses", "NA")),
  mode = "everything"
)

# comparison between larger taxonomic groups
# archaea, bacteria, viruses, fungi, protists (incl lower plants), higher plants, metazoa
phylum_list <- list(
  # fungi
  c(
    "Ascomycota",
    "Basidiomycota",
    "Blastocladiomycota",
    "Chytridiomycota",
    "Cryptomycota",
    "Microsporidia",
    "Mucoromycota",
    "Zoopagomycota"
  ),
  # higher plants
  c("Streptophyta"),
  # metazoa
  c(
    "Acanthocephala",
    "Annelida",
    "Arthropoda",
    "Brachiopoda",
    "Bryozoa",
    "Chordata",
    "Cnidaria",
    "Ctenophora",
    "Dicyemida",
    "Echinodermata",
    "Hemichordata",
    "Mollusca",
    "Nematoda",
    "Nemertea",
    "Onychophora",
    "Orthonectida",
    "Phoronida",
    "Placozoa",
    "Platyhelminthes",
    "Porifera",
    "Priapulida",
    "Rotifera",
    "Tardigrada",
    "Xenacoelomorpha"
  )
)
names(phylum_list) <- c("fungi", "plants", "metazoa")

taxon_groups_in <- vector(mode = "character", length = nrow(tax_in_class))
taxon_groups_in[tax_in_class[, 1] == "Bacteria"] <- "Bacteria"
taxon_groups_in[tax_in_class[, 1] == "Archaea"] <- "Archaea"
taxon_groups_in[tax_in_class[, 1] == "Viruses"] <- "Viruses"
taxon_groups_in[tax_in_class[, 1] == "NA"] <- "NA"
taxon_groups_in[tax_in_class[, 1] == "Eukaryota" & tax_in_class[, 2] %in% phylum_list$fungi] <- "Fungi"
taxon_groups_in[tax_in_class[, 1] == "Eukaryota" & tax_in_class[, 2] %in% phylum_list$plants] <- "Plants"
taxon_groups_in[tax_in_class[, 1] == "Eukaryota" & tax_in_class[, 2] %in% phylum_list$metazoa] <- "Metazoa"
taxon_groups_in[tax_in_class[, 1] == "Eukaryota" & !tax_in_class[, 2] %in% unlist(phylum_list)] <- "Protists"

taxon_groups_out <- vector(mode = "character", length = nrow(tax_out_class))
taxon_groups_out[tax_out_class[, 1] == "Bacteria"] <- "Bacteria"
taxon_groups_out[tax_out_class[, 1] == "Archaea"] <- "Archaea"
taxon_groups_out[tax_out_class[, 1] == "Viruses"] <- "Viruses"
taxon_groups_out[tax_out_class[, 1] == "NA"] <- "NA"
taxon_groups_out[tax_out_class[, 1] == "Eukaryota" & tax_out_class[, 2] %in% phylum_list$fungi] <- "Fungi"
taxon_groups_out[tax_out_class[, 1] == "Eukaryota" & tax_out_class[, 2] %in% phylum_list$plants] <- "Plants"
taxon_groups_out[tax_out_class[, 1] == "Eukaryota" & tax_out_class[, 2] %in% phylum_list$metazoa] <- "Metazoa"
taxon_groups_out[tax_out_class[, 1] == "Eukaryota" & !tax_out_class[, 2] %in% unlist(phylum_list)] <- "Protists"

cm_taxon_groups <- confusionMatrix(
  data = factor(taxon_groups_out, levels = c("Bacteria", "Archaea", "Viruses", "NA", "Fungi", "Plants", "Metazoa", "Protists")),
  reference = factor(taxon_groups_in, levels = c("Bacteria", "Archaea", "Viruses", "NA", "Fungi", "Plants", "Metazoa", "Protists")),
  mode = "everything"
)


### write output ####

write.table(
  data.frame(cm_domain$table),
  paste0(opt$output, "_cm_domain.txt"),
  sep = "\t",
  row.names = F,
  quote = F
)

write.table(
  data.frame(cm_taxon_groups$table),
  paste0(opt$output, "_cm_taxongroups.txt"),
  sep = "\t",
  row.names = F,
  quote = F
)

write.table(
  data.frame(cm_domain$byClass),
  paste0(opt$output, "_stats_domain.txt"),
  sep = "\t",
  quote = F
)

write.table(
  data.frame(cm_taxon_groups$byClass),
  paste0(opt$output, "_stats_taxongroups.txt"),
  sep = "\t",
  quote = F
)