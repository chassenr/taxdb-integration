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
    c("-i", "--kaiju"), 
    type = "character", 
    default = NULL,
    help = "kaiju output", 
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
    c("-o", "--output"), 
    type = "character", 
    default = NULL,
    help = "base name for output files",
    metavar = "character"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$kaiju) | is.null(opt$output) | is.null(opt$taxdump) | is.null(opt$sql)) {
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


### read kaiju output ####
kaiju_out <- fread(
  opt$kaiju, 
  header = FALSE,
  fill = TRUE, 
  sep = "\t"
)
colnames(kaiju_out) <- c("classified", "tax_path_in", "taxid_out")

# assign taxonomy to output
# convert any non alpha-numeric character to underscore
# (this is done for input taxonomic path by simulation program)
tax_out_class <- getTaxonomy(
  kaiju_out$taxid_out,
  opt$sql
) %>% 
  gsub("[ ;]", "_", .) %>%
  gsub("[dpcofgs]__", "", .)

# recode NA as character
tax_out_class[is.na(tax_out_class)] <- "NA"

# splitting input taxonomic path
tax_in_class <- kaiju_out$tax_path_in %>% 
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


### calculate confusion matrix ####

# comparison between domains
cm_domain <- confusionMatrix(
  data = factor(data.frame(tax_out_class)[, 1], levels = c("Eukaryota", "Bacteria", "Archaea", "Viruses", "NA")),
  reference = factor(data.frame(tax_in_class)[, 1], levels = c("Eukaryota", "Bacteria", "Archaea", "Viruses", "NA")),
  mode = "everything"
)


# comparison on each taxonomic level:
#   was correct taxon assigned?
#   was incorrect taxon assigned?
#   was no taxon assigned?

tax_out_cat <- data.frame(tax_out_class, stringsAsFactors = F)
for(i in 1:ncol(tax_out_cat)) {
  tax_out_cat[, i] <- factor(
    base::ifelse(
      tax_out_class[, i] == "NA",
      "unassigned",
      ifelse(
        tax_out_class[, i] == tax_in_class[, i],
        "correct",
        "incorrect"
      )
    ),
    levels = c("correct", "incorrect", "unassigned")
  )
}
summary_taxlevels <- tax_out_cat %>%
  mutate(domain = tax_in_class[, 1]) %>%
  pivot_longer(., 1:7) %>%
  group_by(domain, name, value) %>% 
  count() %>%
  pivot_wider(., names_from = "name", values_from = "n") %>%
  relocate("superkingdom", "phylum", "class", "order", "family", "genus", "species", .after = "value")


### write output ####

write.table(
  data.frame(cm_domain$table),
  paste0(opt$output, "_cm_domain.txt"),
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
  data.frame(summary_taxlevels),
  paste0(opt$output, "_taxlevels.txt"),
  sep = "\t",
  quote = F
)


