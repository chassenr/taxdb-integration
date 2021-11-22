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
    c("-i", "--input"), 
    type = "character", 
    default = NULL,
    help = "taxonomy file", 
    metavar = "character"
  ),
  make_option(
    c("-o", "--output"), 
    type = "character", 
    default = NULL,
    help = "taxonomy with duplicates resolved",
    metavar = "character"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
if (is.null(opt$input) | is.null(opt$output)) {
  print_help(opt_parser)
  stop(
    "You need to provide the input and output taxonomy files.\n", 
    call. = FALSE
  )
}


### fix NCBI duplicate tax levels ####

# this includes all concatenated taxonomy
taxpath <- fread(
  opt$input,
  h = F,
  sep = "\t",
  quote = "",
  col.names = c("accnos", "path")
) %>%
  arrange(accnos) %>%
  filter(duplicated(accnos) == FALSE) %>%
  separate(
    path,
    sep = ";",
    into = c("superkingdom", "lineage", "kingdom", "phylum", "class", "order", "family", "genus", "species")
  )

# add index for duplicated tax names within NCBI taxonomy
# adapt for tidyverse later
taxpath_df <- as.data.frame(taxpath)
taxpath_fixed <- taxpath_df
for(j in 3:ncol(taxpath_df)) {
  for(i in unique(taxpath_df[, j])) {
    tmp1 <- taxpath_df[taxpath_df[, j] == i, j-1]
    if(length(unique(tmp1)) > 1) {
      tmp2 <- as.factor(tmp1)
      levels(tmp2) <- LETTERS[1:length(unique(tmp1))]
      taxpath_fixed[taxpath_df[, j] == i, j] <- paste0(taxpath_df[taxpath_df[, j] == i, j], "_", as.character(tmp2))
    }
  }
} 

# add mapping between old and new taxpaths
taxpath_out <- taxpath_df %>%
  select(accnos) %>% 
  mutate(
    taxpath_fixed %>% 
      select(-accnos) %>% 
      unite("path_fixed", sep = ";"),
    taxpath_df %>% 
      select(-accnos) %>% 
      unite("path_old", sep = ";")
  ) %>% 
  as_tibble()
  
# write output
write_delim(
  taxpath_out,
  opt$output,
  delim = "\t",
  col_names = F
)
