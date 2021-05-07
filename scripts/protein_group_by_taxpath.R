#!/usr/bin/env Rscript

# parse taxonomic path for NCBI genomes

### setting up environment ####

# Check if packages are installed
package.list <- c(
  "crayon",
  "optparse",
  "config",
  "data.table",
  "tidyverse"
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
    help = "protein taxonomy file",
    metavar = "character"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = NULL,
    help = "list of taxonomic paths that will be used to groups assemblies for clustering",
    metavar = "character"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
if (is.null(opt$input) | is.null(opt$output)) {
  print_help(opt_parser)
  stop(
    "You need to provide the input and output files.\n",
    call. = FALSE
  )
}


### Selecting taxonomic paths ####

taxpath <- fread(
  opt$input,
  h = F,
  sep = "\t",
  quote = "",
  col.names = c("accnos", "path")
) %>%
  separate(
    path,
    sep = ";",
    into = c("superkingdom", "phylum", "class", "order", "family", "genus", "species"),
    remove = F
  ) %>%
  # choose genus level unless genus is unknown (then take species)
  unite("path_genus", superkingdom:genus, sep = ";", remove = F) %>%
  mutate(
    taxpath_group = if_else(
      str_replace(genus, "g__", "") == str_replace(family, "f__", ""),
      path,
      str_replace(path_genus, "$", ";")
    )
  )

taxpath_grouped <- unique(taxpath$taxpath_group)

msg(paste0("Read ", nrow(taxpath), " accessions for ", length(unique(taxpath$path)), " species.\n"))
msg(paste0("Identified ", length(taxpath_grouped), " taxa for clustering.\n"))

# write output
write_delim(
  data.frame(taxpath_grouped),
  opt$output,
  delim = "\t",
  col_names = F
)



