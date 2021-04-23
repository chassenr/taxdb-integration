#!/usr/bin/env Rscript

# parse assembly metadata for NCBI genomes

### setting up environment ####

# Check if packages are installed
package.list <- c(
  "crayon",
  "optparse"
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
    c("-t", "--tax"),
    type = "character", 
    default = NULL, 
    help = "taxonomy table for assemblies",
    metavar = "character"
  ),
  make_option(
    c("-m", "--meta"),
    type = "character",
    default = NULL, 
    help = "assembly metadata (column names have to match selected parameters)", 
    metavar = "character"
  ),
  make_option(
    c("-n", "--ncbi"),
    type = "character",
    default = NULL,
    help = "assembly summary table from NCBI",
    metavar = "character"
  ),
  make_option(
    c("-p", "--parameters"),
    type = "character", 
    default = NULL, 
    help = "comma-separated list of parameter names to use for filtering",
    metavar = "character"
  ),
  make_option(
    c("-a", "--absolute_cutoff"),
    type = "character",
    default = NULL, 
    help = "comma-separated list of absolute cut-off to use for each parameter", 
    metavar = "character"
  ),
  make_option(
    c("-r", "--relative_cutoff"),
    type = "character",
    default = NULL, 
    help = "comma-separated list of absolute cut-off to use for each parameter", 
    metavar = "character"
  ),
  make_option(
    c("-d", "--difference"),
    type = "character",
    default = NULL, 
    help = "comma-separated list of minimum difference for relative cutoff to be implemented", 
    metavar = "character"
  ),
  make_option(
    c("-s", "--sign"),
    type = "character",
    default = NULL, 
    help = "comma-separated list of sign of filter (greater, less, both)", 
    metavar = "character"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = NULL, 
    help = "Name of output file with combined metadata", 
    metavar = "character"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$tax) | is.null(opt$meta) | is.null(opt$ncbi) |
    is.null(opt$parameters) | is.null(opt$absolute_cutoff) |
    is.null(opt$relative_cutoff) | is.null(opt$difference) |
    is.null(opt$sign) | is.null(opt$output)) {
  print_help(opt_parser)
  stop("All parameters are mandatory.\n", call. = FALSE)
}


### define functions ####

filter_assemblies <- function(accnos, meta, abs.cut, rel.cut, cut.dif, sign) {
  if(sign == "greater") {
    if(abs.cut == -1 | sum(meta >= abs.cut) == 0) {
      threshold <- quantile(meta, rel.cut)
      if(median(meta) - threshold >= cut.dif) {
        accnos.select <- accnos[meta >= threshold]
      } else {
        accnos.select <- accnos
      }
    } else {
      accnos.select <- accnos[meta >= abs.cut]
    }
  }
  if(sign == "less") {
    if(abs.cut == -1 | sum(meta <= abs.cut) == 0) {
      threshold <- quantile(meta, rel.cut)
      if(threshold - median(meta) >= cut.dif) {
        accnos.select <- accnos[meta <= threshold]
      } else {
        accnos.select <- accnos
      }
    } else {
      accnos.select <- accnos[meta <= abs.cut]
    }
  }
  if(sign == "both") {
    threshold <- quantile(meta, c(rel.cut, 1 - rel.cut))
    if(threshold[2] - threshold[1] >= 2 * cut.dif & sum(meta >= threshold[1] & meta <= threshold[2]) > 0) {
      accnos.select <- accnos[meta >= threshold[1] & meta <= threshold[2]]
    } else {
      accnos.select <- accnos
    }
  }
  return(accnos.select)
}


### read input data ####

metadata <- read.table(
  opt$meta, 
  h = T, 
  sep = "\t",
  stringsAsFactors = F
)

tax <- read.table(
  opt$tax, 
  h = F, 
  sep = "\t",
  stringsAsFactors = F
)

ncbi <- read.table(
  opt$ncbi,
  h = F,
  sep = "\t",
  quote = "",
  comment.char = "",
  stringsAsFactors = F
)

metadata <- metadata[match(tax$V1, metadata$accession), ]
ncbi <- ncbi[match(tax$V1, ncbi$V1), ]

if(!all.equal(metadata$accession, tax$V1) | !all.equal(ncbi$V1, tax$V1)) {
  stop("taxonomy and metadata files are not in the correct order", call. = FALSE)
}

meta <- data.frame(
  metadata,
  ncbi,
  stringsAsFactors = F
)

### parse command line options ####

params <- strsplit(opt$parameters, ",")[[1]]
abs_cut <- as.numeric(strsplit(opt$absolute_cutoff, ",")[[1]])
rel_cut <- as.numeric(strsplit(opt$relative_cutoff, ",")[[1]])
cut_dif <- as.numeric(strsplit(opt$difference, ",")[[1]])
cut_sign <- strsplit(opt$sign, ",")[[1]]

if(any(!params %in% colnames(metadata))) {
  stop("Parameters for filtering not found in metadata", call. = FALSE)
}

### some pre-filtering to remove inconsistently small assemblies
# if species has representative genome, remove assemblies that are smaller than 5% of genome length of representative

prefilter_accnos <- c()
for(i in unique(tax$V2)) {
  tmp <- meta[tax$V2 == i, ]
  if(any(tmp$V5 != "na")) {
    prefilter_accnos <- c(
      prefilter_accnos,
      tmp$accession[tmp$genome_size >= 0.05 * min(tmp$genome_size[tmp$V5 != "na"])]
    )
  } else {
    prefilter_accnos <- c(
      prefilter_accnos,
      tmp$accession
    )
  }
}

meta_pre <- meta[meta$accession %in% prefilter_accnos, ]
tax_pre <- tax[meta$accession %in% prefilter_accnos, ]
msg(paste0(nrow(tax_pre), " assemblies passed the genome size pre-filter\n"))


### loop through individual filters ####

metadata_filt <- meta_pre
tax_filt <- tax_pre
for(k in 1:length(params)) {
  filter_accnos <- c()
  for(i in unique(tax_filt$V2)) {
    tmp <- metadata_filt[tax_filt$V2 == i, ]
    filter_accnos <- c(
      filter_accnos,
      filter_assemblies(
        tmp$accession, 
        tmp[, params[k]], 
        abs.cut = abs_cut[k], 
        rel.cut = rel_cut[k], 
        cut.dif = ifelse(cut_dif[k] >= 1, cut_dif[k], cut_dif[k] * median(tmp[, params[k]])), 
        sign = cut_sign[k]
      )
    )
  }
  metadata_next <- metadata_filt[metadata_filt$accession %in% filter_accnos, ]
  tax_next <- tax_filt[metadata_filt$accession %in% filter_accnos, ]
  msg(paste0(nrow(tax_next), " assemblies passed the ", params[k], " filter\n"))
  metadata_filt <- metadata_next
  tax_filt <- tax_next
}


### write outputfile (only taxonomy) ####
msg("writing output...")
write.table(
  tax_filt,
  opt$output,
  quote = F,
  sep = "\t",
  row.names = F,
  col.names = F
)
cat(" done\n")
