#!/usr/bin/env Rscript

# parse taxonomic path for NCBI genomes

### setting up environment ####

# Check if packages are installed
package.list <- c(
  "crayon",
  "optparse",
  "config",
  "tidyverse",
  "seqinr",
  "IRanges"
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
    help = "fasta file with contaminated sequence", 
    metavar = "character"
  ),
  make_option(
    c("-c", "--contam_pred"), 
    type = "character", 
    default = NULL,
    help = "predicted contamination (conterminator output table)", 
    metavar = "character"
  ),
  make_option(
    c("-o", "--output"), 
    type = "character", 
    default = NULL,
    help = "file name for cleaned output sequences",
    metavar = "character"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input) | is.null(opt$output) | is.null(opt$contam_pred)) {
  print_help(opt_parser)
  stop(
    "All parameters are mandatory.\n", 
    call. = FALSE
  )
}

### reading fasta and trimming based on alignment positions ####
# for now not just the contamination but the whole contig are removed

fna <- read.fasta(opt$input, forceDNAtolower = F)
msg(paste0(length(fna), " sequences contain contamination\n"))
contam <- read.table(opt$contam_pred, h = F, sep = "\t", stringsAsFactors = F)
msg(paste0(nrow(contam), " contaminations have been predicted\n"))

# collect clean seqs into a list
cleaned <- list()

# for each contaminated sequence (scaffold)
for(i in names(fna)) {
  fna.sub <- fna[[i]]
  contam.sub <- contam[contam$V2 == i, ]
  contam.ranges <- data.frame(reduce(IRanges(contam.sub$V5,contam.sub$V6)))
  # define the ranges of clean sequences
  clean.ranges <- data.frame(
    start = c(1, contam.ranges[, 2] + 1),
    end = c(contam.ranges[, 1] - 1, length(fna.sub))
  )
  # get positions of Ns to identify contigs in scaffolds
  pos.n <- which(fna.sub %in% c("n", "N"))
  # if no Ns (i.e. no scaffolding), remove whole contig (i.e. don't write any output)
  if(length(pos.n) == 0) {
    cleaned <- cleaned
  } else {
    # split scaffold into clean sections
    fna.split <- apply(
      clean.ranges,
      1,
      function(x) {
        fna.sub[x[1]:x[2]]
      }
    )
    # get position of Ns in each clean section
    pos.n.split <- lapply(fna.split, function(x) which(x %in% c("n", "N")))
    fna.cleaned <- vector("list", length = length(fna.split))
    # for each clean section
    for(k in 1:length(fna.split)) {
      # only process if Ns are present, otherwise remove, since contamination from same contig
      if(length(pos.n.split[[k]]) > 0) {
        # trim leading and trailing Ns and masked low-complexity regions
        fna.cleaned[[k]] <- fna.split[[k]][
          ifelse(k == 1, 1, pos.n.split[[k]][1]):ifelse(k == length(fna.split), length(fna.split[[k]]), pos.n.split[[k]][length(pos.n.split[[k]])])
        ] %>% 
          paste(., collapse = "") %>% 
          gsub("[nNx]*$", "", .) %>% 
          gsub("^[nNx]*", "", .) %>% 
          strsplit(., "") %>% 
          unlist()
      }
    }
    # filter out sequences where trimming removed whole sequence
    fna.cleaned <- fna.cleaned[!sapply(fna.cleaned, is.null) & sapply(fna.cleaned, length) != 0]
    # if none are left, don't write output
    if(length(fna.cleaned) == 0) {
      cleaned <- cleaned
    } else {
      # if only 1 sequence is left, keep name of original sequence
      if(length(fna.cleaned) == 1) {
        names(fna.cleaned) <- i
      } else {
        # if more than 1 split sequence survived, add additional index to sequence ID
        names(fna.cleaned) <- paste0(
          gsub("\\|kraken:taxid\\|[0-9]*$", "_", i),
          1:length(fna.cleaned),
          gsub("^.*\\|kraken:taxid\\|", "\\|kraken:taxid\\|", i)
        )
      }
      cleaned <- c(cleaned, fna.cleaned)
    }
  }
}
msg(paste0(length(cleaned), " sequences were obtained after cleaning\n"))

# write output
write.fasta(
  cleaned,
  names(cleaned),
  file.out = opt$output
)
