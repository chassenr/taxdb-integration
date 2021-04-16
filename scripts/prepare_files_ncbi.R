#!/usr/bin/env Rscript

# parse NCBI assembly metadata to select genomes for download

### setting up environment ####

# Check if packages are installed
package.list <- c(
  "crayon",
  "optparse",
  "tidyverse",
  "data.table",
  "httr",
  "furrr"
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
    c("-g", "--genbank"),
    type = "character", 
    default = NULL, 
    help = "genbank assembly summary file",
    metavar = "character"
  ),
  make_option(
    c("-r", "--refseq"),
    type = "character",
    default = NULL, 
    help = "refseq assembly summary file", 
    metavar = "character"
  ),
  make_option(
    c("-a", "--assembly-level"),
    type = "character", 
    default = NULL,
    help = "minimum assembly level",
    metavar = "character"
  ),
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
    help = "Name of output file with download links", 
    metavar = "character"
  ),
  make_option(
    c("-s", "--summary"),
    type = "character",
    default = NULL, 
    help = "Name of output file combined assembly summary", 
    metavar = "character"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$refseq) | is.null(opt$genbank) | 
    is.null(opt$'assembly-level') |
    is.null(opt$output) | is.null(opt$summary)) {
  print_help(opt_parser)
  stop("All parameters are mandatory.\n", call. = FALSE)
}


### Filter NCBI genomes ####

# refseq summary
# filter out any genomes without download link
refseq <- read.table(
  file = opt$refseq, 
  sep = "\t",
  h = F,
  skip = 2,
  quote = "",
  comment.char = "",
  stringsAsFactors = F
) %>% 
  separate(
    V1,
    into = c("group", "acc_short"), 
    sep = "_", 
    remove = FALSE, 
    extra = "drop"
  ) %>% 
  filter(V20 != "na")


# genbank summary
# remove genbank genomes from list, which are redundant with refseq
# only use the latest assembly version of genbank assemblies with full genome representation
# also remove genomes with missing links
genbank <- read.table(
  file = opt$genbank, 
  sep = "\t",
  h = F,
  skip = 2,
  quote = "",
  comment.char = "",
  stringsAsFactors = F
) %>% 
  separate(
    V1,
    into = c("group", "acc_short"), 
    sep = "_", 
    remove = FALSE, 
    extra = "drop"
  ) %>% 
  filter(
    !acc_short %in% refseq$acc_short,
    V11 == "latest",
    V14 == "Full",
    V20 != "na"
  )

# filter genbank assemblies further by assembly_level
if(opt$'assembly-level' == "Contig") {
  genbank_clean <- genbank
}
if(opt$'assembly-level' == "Scaffold") {
  genbank_clean <- genbank %>% 
    filter(V12 != "Contig")
}
if(opt$'assembly-level' == "Chromosome") {
  genbank_clean <- genbank %>% 
    filter(V12 %in% c("Chromosome", "Complete Genome"))
}
if(opt$'assembly-level' == "Complete Genome") {
  genbank_clean <- genbank %>% 
    filter(V12 == "Complete Genome")
}
if(opt$'assembly-level' == "variable") {
  genbank_clean <- genbank %>% 
    group_by(V7) %>% 
    group_modify(
      ~ {
        if(any(.x$V12 %in% c("Chromosome", "Complete Genome"))) {
          .x %>% filter(.x$V12 %in% c("Chromosome", "Complete Genome"))
        } else {
          .x
        }
      }
    )
}

# combine assembly summary tables (and check for missing links in refseq file
summary_combined <- rbind(
  refseq,
  genbank_clean %>% ungroup()
) %>% 
  select(-group, -acc_short)

# run url check in parallel
plan(multisession, workers = opt$cpus)

# format input file for aria2 (nucleotide and protein)
summary_combined_links <- summary_combined %>% 
  mutate(
    link = gsub("^ftp", "http", V20),
    file_fna = paste0(basename(V20), "_genomic.fna.gz"),
    file_faa = paste0(basename(V20), "_protein.faa.gz"),
    link_fna = paste(link, file_fna, sep = "/"),
    link_faa = paste(link, file_faa, sep = "/"),
    # check if protein file is available
    faa_available = !future_map_lgl(
      link_faa,
      http_error
    )
  ) %>% 
  select(-link, -file_fna, -file_faa)

# only keep genomes with proteins if available for each taxid
summary_combined_with_faa <- summary_combined_links %>% 
  filter(faa_available)
summary_combined_no_faa <- summary_combined_links %>% 
  filter(!V7 %in% summary_combined_with_faa$V7)
summary_combined_good <- bind_rows(summary_combined_with_faa, summary_combined_no_faa)

# print some summary stats
msg(
  paste0(
    "There are ",
    nrow(refseq),
    " RefSeq assemblies.\n"
  )
)
msg(
  paste0(
    "There are ",
    nrow(genbank), 
    " Genbank assemblies, ", 
    nrow(genbank_clean),
    " were retained at the selected assembly level.\n"
  )
)
msg(
  paste0(
    nrow(summary_combined_with_faa),
    " genomes of ",
    length(table(summary_combined_with_faa$V7)),
    " species have predicted proteins, ", 
    length(table(summary_combined_no_faa$V7)),
    " species (",
    nrow(summary_combined_no_faa),
    " genomes) do not have any genomes with predicted proteins.\n"
  )
)


### write output ###

write.table(
  summary_combined_good[, 1:22],
  opt$summary,
  quote = F,
  sep = "\t",
  row.names = F,
  col.names = F
)

write.table(
  summary_combined_good[, 23:25],
  opt$output,
  quote = F,
  sep = "\t",
  row.names = F,
  col.names = F
)
