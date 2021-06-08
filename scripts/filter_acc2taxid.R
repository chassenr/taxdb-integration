#!/usr/bin/env Rscript

# parse taxonomic path for NCBI genomes

### setting up environment ####

# Check if packages are installed
package.list <- c(
  "crayon",
  "optparse",
  "config",
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
    c("-i", "--input"),
    type = "character",
    default = NULL,
    help = "accesions of clustered protein representatives",
    metavar = "character"
  ),
  make_option(
    c("-a", "--acc2taxid"),
    type = "character",
    default = NULL,
    help = "protein accession2taxid file",
    metavar = "character"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = NULL,
    help = "accession2taxid file for clustered proteins (CAUTION: genome and sequence accession are already pasted together with _",
    metavar = "character"
  ),
    make_option(
    c("-c", "--cpus"),
    type = "integer",
    default = 1,
    help = "number of cpus to use [default: 1]",
    metavar = "number"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input) | is.null(opt$output) | is.null(opt$acc2taxid)) {
  print_help(opt_parser)
  stop(
    "You need to provide the assembly summary table, the taxdump sql database, and the name of the output file.\n",
    call. = FALSE
  )
}


### filter by clustered protein accessions ####

acc <- fread(opt$acc2taxid, h = F, sep = "\t", nThread = opt$cpus)
# columns: genome accession, sequence accession, taxid
acc$V4 <- paste(acc$V1, acc$V2, sep = "_")
setkey(acc, V4)
prot <- fread(opt$input,  h = F, sep = "\t", nThread = opt$cpus)
acc_prot <- acc[prot$V1]
fwrite(acc_prot[, c(4,3)], opt$output, quote = F, sep = "\t", nThread = opt$cpus, col.names = F)

msg(paste0("There are ", nrow(prot), " sequences in the input, and ", nrow(acc_prot), " in the output.\n"))


