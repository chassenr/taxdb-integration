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
  "Biostrings"
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
    help = "faa file to be formatted",
    metavar = "character"
  ),
#  make_option(
#    c("-t", "--taxid"),
#    type = "character",
#    default = NULL,
#    help = "accession2taxid for input faa",
#    metavar = "character"
#  ),
  make_option(
    c("-c", "--cpus"),
    type = "integer",
    default = 1,
    help = "number of cpus to use [default: 1]",
    metavar = "number"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = NULL,
    help = "formatted faa output",
    metavar = "character"
  ),
  make_option(
    c("-m", "--map"),
    type = "character",
    default = NULL,
    help = "mapping file of old to new fasta headers",
    metavar = "character"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$output) | is.null(opt$map) | is.null(opt$input)) {
  print_help(opt_parser)
  stop(
    "You need to provide the names of the input and output faa and the name of the faa accesion mapping file.\n",
    call. = FALSE
  )
}


### format fasta headers for kaiju ####
# prepend taxid by running number (for each taxid)

faa <- readAAStringSet(opt$input)
# acc2taxid <- fread(opt$taxid, h = F, sep = "\t", col.names = c("accnos", "taxid"), nThread = opt$cpus) %>%
#   group_by(taxid) %>%
#   mutate(id = row_number(), accnos_new = paste(id, taxid, sep = "_"))
acc2taxid <- names(faa) %>% 
  as_tibble() %>% 
  separate(
    value, 
    into = c("taxid", "accnos", "add_info"), 
    sep = " ", 
    extra = "merge",
    fill = "right"
  ) %>%
  group_by(taxid) %>%
  mutate(id = row_number(), accnos_new = paste(id, taxid, sep = "_"))
names(faa) <- acc2taxid$accnos_new
writeXStringSet(faa, opt$output)
acc2taxid %>%
  ungroup() %>%
  select(taxid, accnos, accnos_new) %>%
  fwrite(., opt$map, quote = F, sep = "\t", nThread = opt$cpus, col.names = F)


