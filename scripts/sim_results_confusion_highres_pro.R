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
    help = "directory containing nodes.dmp and names.dmp of the taxonomy used for the highres classification", 
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
    c("-c", "--coarse"), 
    type = "character", 
    default = NULL,
    help = "file name of conifer (kraken) output of best scenario for coarse classification", 
    metavar = "character"
  ),
  make_option(
    c("-T", "--taxdump_coarse"), 
    type = "character", 
    default = NULL,
    help = "directory containing nodes.dmp and names.dmp of the taxonomy used for the highres classification", 
    metavar = "character"
  ),
  make_option(
    c("-S", "--sql_coarse"), 
    type = "character", 
    default = NULL,
    help = "location and name of sql database that will be generated from the taxdump", 
    metavar = "character"
  ),
  make_option(
    c("-R", "--rtl_coarse"), 
    type = "double", 
    default = 0,
    help = "RTL cut-off to be applied to coarse classification [default: 0]", 
    metavar = "number"
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

if (is.null(opt$kraken) | is.null(opt$output) | is.null(opt$taxdump) | 
    is.null(opt$sql) | is.null(opt$mode) | is.null(opt$coarse)) {
  print_help(opt_parser)
  stop(
    "You need to provide the assembly summary table, the location of the taxdump files and sql database, and the name of the output file.\n", 
    call. = FALSE
  )
}


### retrieve taxonomy ####

# format highres taxdump database
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

# format coarse taxdump database
if(!file.exists(opt$sql)) {
  read.names.sql(
    paste0(opt$taxdump_coarse, "/names.dmp"),
    opt$sql_coarse
  )
  read.nodes.sql(
    paste0(opt$taxdump_coarse, "/nodes.dmp"),
    opt$sql_coarse
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


### read coarse classification ####
kraken_coarse <- fread(
  opt$coarse, 
  header = FALSE,
  fill = TRUE, 
  sep = "\t"
)[, ..select_cols]
colnames(kraken_coarse) <- c("tax_path_in", "taxid_out", "conf", "rtl")

# assign taxonomy to output
# convert any non alpha-numeric character to underscore
# (this is done for input taxonomic path by simulation program)
tax_out_coarse <- getTaxonomy(
  kraken_coarse$taxid_out,
  opt$sql_coarse
) %>% 
  gsub("[ ;]", "_", .)

# recode NA as character
tax_out_coarse[is.na(tax_out_coarse)] <- "NA"

# apply RTL threshold
tax_out_coarse[kraken_coarse$rtl < opt$rtl_coarse, ] <- "NA"


### apply filters ####

# only keep sequences affiliated with Bacteria, Archaea, Viruses in both the coarse and highres output
tax_out_class_filt <- tax_out_class[tax_out_coarse[, 1] %in% c("Bacteria", "Archaea", "Viruses") & tax_out_class[, 1] %in% c("Bacteria", "Archaea", "Viruses"), ]
tax_in_class_filt <- tax_in_class[tax_out_coarse[, 1] %in% c("Bacteria", "Archaea", "Viruses") & tax_out_class[, 1] %in% c("Bacteria", "Archaea", "Viruses"), ]

# only keep sequences affiliated with Bacteria, Archaea, Viruses in the highres output (no coarse filter)
tax_out_class_nofilt <- tax_out_class[tax_out_class[, 1] %in% c("Bacteria", "Archaea", "Viruses"), ]
tax_in_class_nofilt <- tax_in_class[tax_out_class[, 1] %in% c("Bacteria", "Archaea", "Viruses"), ]


### calculate confusion matrix ####

# comparison between domains
cm_domain_filt <- confusionMatrix(
  data = factor(data.frame(tax_out_class_filt)[, 1], levels = c("Eukaryota", "Bacteria", "Archaea", "Viruses", "NA")),
  reference = factor(data.frame(tax_in_class_filt)[, 1], levels = c("Eukaryota", "Bacteria", "Archaea", "Viruses", "NA")),
  mode = "everything"
)
cm_domain_nofilt <- confusionMatrix(
  data = factor(data.frame(tax_out_class_nofilt)[, 1], levels = c("Eukaryota", "Bacteria", "Archaea", "Viruses", "NA")),
  reference = factor(data.frame(tax_in_class_nofilt)[, 1], levels = c("Eukaryota", "Bacteria", "Archaea", "Viruses", "NA")),
  mode = "everything"
)

# comparison on each taxonomic level:
#   was correct taxon assigned?
#   was incorrect taxon assigned?
#   was no taxon assigned?

tax_out_cat_filt <- data.frame(tax_out_class_filt, stringsAsFactors = F)
for(i in 1:ncol(tax_out_cat_filt)) {
  tax_out_cat_filt[, i] <- factor(
    ifelse(
      tax_in_class_filt[, 1] %in% c("Bacteria", "Archaea", "Viruses"),
      ifelse(
        tax_out_class_filt[, i] == tax_in_class_filt[, i],
        "correct",
        "incorrect_pro"
      ),
      ifelse(
        tax_out_class_filt[, i] == "NA",
        "unassigned",
        "incorrect_euk"
      )
    ),
    levels = c("correct", "incorrect_pro", "incorrect_euk", "unassigned")
  )
}
summary_taxlevels_filt <- apply(tax_out_cat_filt, 2, table) %>% 
  bind_rows() %>% 
  mutate(taxlevel = colnames(tax_out_cat_filt))

tax_out_cat_nofilt <- data.frame(tax_out_class_nofilt, stringsAsFactors = F)
for(i in 1:ncol(tax_out_cat_nofilt)) {
  tax_out_cat_nofilt[, i] <- factor(
    ifelse(
      tax_in_class_nofilt[, 1] %in% c("Bacteria", "Archaea", "Viruses"),
      ifelse(
        tax_out_class_nofilt[, i] == tax_in_class_nofilt[, i],
        "correct",
        "incorrect_pro"
      ),
      ifelse(
        tax_out_class_nofilt[, i] == "NA",
        "unassigned",
        "incorrect_euk"
      )
    ),
    levels = c("correct", "incorrect_pro", "incorrect_euk", "unassigned")
  )
}
summary_taxlevels_nofilt <- apply(tax_out_cat_nofilt, 2, table) %>% 
  bind_rows() %>% 
  mutate(taxlevel = colnames(tax_out_cat_nofilt))


### write output ####

write.table(
  data.frame(cm_domain_filt$table),
  paste0(opt$output, "_cm_domain_filt.txt"),
  sep = "\t",
  row.names = F,
  quote = F
)

write.table(
  data.frame(cm_domain_filt$byClass),
  paste0(opt$output, "_stats_domain_filt.txt"),
  sep = "\t",
  quote = F
)

write.table(
  data.frame(cm_domain_nofilt$table),
  paste0(opt$output, "_cm_domain_nofilt.txt"),
  sep = "\t",
  row.names = F,
  quote = F
)

write.table(
  data.frame(cm_domain_nofilt$byClass),
  paste0(opt$output, "_stats_domain_nofilt.txt"),
  sep = "\t",
  quote = F
)

write.table(
  data.frame(summary_taxlevels_filt),
  paste0(opt$output, "_summary_taxlevels_filt.txt"),
  sep = "\t",
  quote = F
)

write.table(
  data.frame(summary_taxlevels_nofilt),
  paste0(opt$output, "_summary_taxlevels_nofilt.txt"),
  sep = "\t",
  quote = F
)