#!/usr/bin/env Rscript

# parse gtdb taxonomy and metadata
# prepared by Antonio

### setting up environment ####

# Check if packages are installed
package.list <- c(
  "crayon",
  "optparse",
  "tidyverse",
  "data.table",
  "scales",
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
    c("-a", "--ar_taxonomy"),
    type = "character", 
    default = NULL, 
    help = "GTDB Archaea taxonomy",
    metavar = "character"
  ),
  make_option(
    c("-b", "--bac_taxonomy"),
    type = "character",
    default = NULL, 
    help = "GTDB Bacteria taxonomy", 
    metavar = "character"
  ),
  make_option(
    c("-A", "--ar_metadata"),
    type = "character", 
    default = NULL, 
    help = "GTDB Archaea metadata",
    metavar = "character"
  ),
  make_option(
    c("-B", "--bac_metadata"),
    type = "character",
    default = NULL, 
    help = "GTDB Bacteria metadata", 
    metavar = "character"
  ),
  make_option(
    c("-r", "--refseq"),
    type = "character",
    default = NULL, 
    help = "RefSeq summary", 
    metavar = "character"
  ),
  make_option(
    c("-g", "--genbank"),
    type = "character", 
    default = NULL,
    help = "GenBank summary",
    metavar = "character"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = NULL, 
    help = "Name of output file with download links and taxonomy", 
    metavar = "character"
  ),
  make_option(
    c("-m", "--metadata-output"),
    type = "character",
    default = NULL, 
    help = "Name of output file with assembly metadata", 
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

if (is.null(opt$ar_taxonomy) | is.null(opt$bac_taxonomy) |
    is.null(opt$ar_metadata) | is.null(opt$bac_metadata) |
    is.null(opt$refseq) | is.null(opt$genbank) | 
    is.null(opt$output) | is.null(opt$'metadata-output')) {
  print_help(opt_parser)
  stop("All parameters are mandatory.\n", call. = FALSE)
}


### Get necessary files ####

# Get GTDB taxonomy
msg("Reading GTDB taxonomy...")
bac_taxonomy <- read_tsv(
  file = opt$bac_taxonomy, 
  col_names = FALSE, 
  na = "N/A",
  col_types = cols()
) %>%
  setNames(c("gtdb_genome", "tax_string")) %>%
  separate(
    col = gtdb_genome,
    into = c("class", "acc"), 
    extra = "merge",
    remove = FALSE, 
    fill = "left"
  ) %>%
  mutate(class = ifelse(grepl("UBA", gtdb_genome), "UBA", class)) %>%
  separate(
    acc,
    into = "acc_short", 
    sep = "\\.", 
    remove = FALSE, 
    extra = "drop"
  )
arc_taxonomy <- read_tsv(
  file = opt$ar_taxonomy,
  col_names = FALSE, 
  na = "N/A",
  col_types = cols()
) %>%
  setNames(c("gtdb_genome", "tax_string")) %>%
  separate(
    col = gtdb_genome,
    into = c("class", "acc"), 
    extra = "merge", 
    remove = FALSE, 
    fill = "left"
  ) %>%
  mutate(class = ifelse(grepl("UBA", gtdb_genome), "UBA", class)) %>%
  separate(
    acc,
    into = "acc_short",
    sep = "\\.", 
    remove = FALSE, 
    extra = "drop"
  )
# gtdb_taxonomy_cat <- bind_rows(bac_taxonomy, arc_taxonomy)
gtdb_taxonomy <- bind_rows(bac_taxonomy, arc_taxonomy)
cat(" done\n")

# CAUTION: there may be duplicated accessions in gtdb
# pick most recent assembly based on version number
# This bug seems to have been fixed in release 95
# gtdb_taxonomy <- gtdb_taxonomy_cat %>%
#   mutate(
#     accnos = gsub("GCA_|GCF_", "", acc_short),
#     version = gsub(".*\\.", "", acc)
#   ) %>% 
#   group_by(accnos) %>% 
#   arrange(desc(version), .by_group = T) %>% 
#   filter(!duplicated(accnos)) %>% 
#   ungroup() %>% 
#   select(-accnos, -version)

# Report some numbers
msg(paste("Entries in GTDB taxonomy:", comma(nrow(gtdb_taxonomy)), "\n"))
# no UBA anymore in version 95

# Get RefSeq data
msg("Reading RefSeq assembly summaries...")
refseq_metadata <- read_tsv(
  file = opt$refseq, 
  col_names = FALSE, 
  na = "na",
  col_types = cols()
) %>% 
  setNames(c("acc1", "acc2", "link")) %>% 
  separate(
    acc1,
    into = "acc1_short", 
    sep = "\\.", 
    remove = FALSE, 
    extra = "drop"
  ) %>%
  separate(
    acc2,
    into = "acc2_short", 
    sep = "\\.", 
    remove = FALSE, 
    extra = "drop"
  )
cat(" done\n")

# Get GenBank data
msg("Reading GenBank assembly summaries...")
genbank_metadata <- read_tsv(
  file = opt$genbank,
  col_names = FALSE, 
  na = "na",
  col_types = cols()
) %>% 
  setNames(c("acc1", "acc2", "link")) %>% 
  separate(
    acc1,
    into = "acc1_short", 
    sep = "\\.", 
    remove = FALSE, 
    extra = "drop"
  ) %>%
  separate(
    acc2,
    into = "acc2_short", 
    sep = "\\.", 
    remove = FALSE, 
    extra = "drop"
  )
cat(" done\n")

msg(paste("Entries in RefSeq file:", comma(nrow(refseq_metadata)), "\n"))
msg(paste("Entries in GenBank file:", comma(nrow(genbank_metadata)), "\n"))


### Find overlaps between GTDB and NCBI ####

msg("Identifying common genomes between GTDB and RefSeq/GenBank data...")
# search first based on primary accession
gtdb_ncbi_acc1 <- bind_rows(
  inner_join(gtdb_taxonomy, refseq_metadata, by = c("acc_short" = "acc1_short")),
  inner_join(gtdb_taxonomy, genbank_metadata, by = c("acc_short" = "acc1_short"))
)
# if not found, also consider secondary accession
gtdb_ncbi_acc2 <- bind_rows(
  gtdb_taxonomy %>% 
    filter(!acc_short %in% gtdb_ncbi_acc1$acc_short) %>% 
    inner_join(refseq_metadata, by = c("acc_short" = "acc2_short")),
  gtdb_taxonomy %>% 
    filter(!acc_short %in% gtdb_ncbi_acc1$acc_short) %>% 
    inner_join(genbank_metadata, by = c("acc_short" = "acc2_short"))
)
msg(paste("We found", comma(sum(nrow(gtdb_ncbi_acc1), nrow(gtdb_ncbi_acc2))), "GTDB entries in NCBI\n"))

# Which GTDB genomes are not in NCBI
gtdb_ncbi_missing <- gtdb_taxonomy %>%
  filter(class != "UBA") %>%
  filter(!(gtdb_genome %in% c(gtdb_ncbi_acc1$gtdb_genome, gtdb_ncbi_acc2$gtdb_genome))) %>%
  mutate(
    acc_short_mod = case_when(
      grepl("GCF", acc_short) ~ gsub("GCF", "GCA", acc_short),
      grepl("GCA", acc_short) ~ gsub("GCA", "GCF", acc_short)
    )
  )
msg(paste("We are missing", comma(nrow(gtdb_ncbi_missing)), "entries in NCBI (might be updated, removed...)\n"))

# rescue some missing genomes by taking their modified accession number
gtdb_ncbi_rescued1 <- bind_rows(
  inner_join(gtdb_ncbi_missing, refseq_metadata, by = c("acc_short_mod" = "acc1_short")),
  inner_join(gtdb_ncbi_missing, genbank_metadata, by = c("acc_short_mod" = "acc1_short"))
)
gtdb_ncbi_rescued2 <- bind_rows(
  gtdb_ncbi_missing %>% 
    filter(!acc_short_mod %in% gtdb_ncbi_rescued1$acc_short_mod) %>% 
    inner_join(refseq_metadata, by = c("acc_short_mod" = "acc2_short")),
  gtdb_ncbi_missing %>% 
    filter(!acc_short_mod %in% gtdb_ncbi_rescued1$acc_short_mod) %>% 
    inner_join(genbank_metadata, by = c("acc_short_mod" = "acc2_short"))
)

# concatenate download links
gtdb_links <- bind_rows(
  gtdb_ncbi_acc1 %>% select(gtdb_genome, link, tax_string),
  gtdb_ncbi_acc2 %>% select(gtdb_genome, link, tax_string),
  gtdb_ncbi_rescued1 %>% select(gtdb_genome, link, tax_string),
  gtdb_ncbi_rescued2 %>% select(gtdb_genome, link, tax_string)
) %>% 
  filter(!is.na(link))

msg(paste("Generating", comma(nrow(gtdb_links %>% filter(!grepl("UBA", gtdb_genome)))), "links for download...\n"))


### Create necessary files For downloading fasta files ####
# and check for availability of protein data

# run url check in parallel
plan(multisession, workers = opt$cpus)

out_file_links <- gtdb_links %>%
  mutate(
    link = gsub("^ftp", "http", link),
    filename = paste0(basename(link), "_genomic.fna.gz"),
    download_link = paste(link, filename, sep = "/"),
    outfile = paste0(gsub("GB_|RS_", "", gtdb_genome), "_genomic.fna.gz"),
    acc = gsub("GB_|RS_", "", gtdb_genome),
    md5file = paste(link, "md5checksums.txt", sep = "/"),
    filename_faa = paste0(basename(link), "_protein.faa.gz"),
    download_link_faa = paste(link, filename_faa, sep = "/"),
    outfile_faa = paste0(gsub("GB_|RS_", "", gtdb_genome), "_protein.faa.gz"),
    # adding artificial taxlevels for kingdom and lineage to make taxonomy compatible with euks
    tax_string = gsub("d__Bacteria;", "d__Bacteria;l__Bacteria;k__Bacteria;", tax_string),
    tax_string = gsub("d__Archaea;", "d__Archaea;l__Archaea;k__Archaea;", tax_string),
    faa_available = !future_map_lgl(
      download_link_faa,
      http_error
    )
  ) %>%
  select(gtdb_genome, download_link, filename, outfile, md5file, acc, tax_string, download_link_faa, filename_faa, outfile_faa, faa_available)


### get assembly metadata ####
msg("Reading GTDB metadata...")
bac_metadata <- fread(
  opt$bac_metadata, 
  h = T, 
  sep = "\t"
) %>%
  as_tibble() %>% 
  select(
    accession,
    gtdb_representative,
    checkm_completeness,
    checkm_contamination,
    checkm_strain_heterogeneity,
    contig_count,
    genome_size,
    longest_contig,
    n50_contigs,
    protein_count,
    gc_percentage
  )
ar_metadata <- fread(
  opt$ar_metadata, 
  h = T, 
  sep = "\t"
) %>%
  as_tibble() %>% 
  select(
    accession,
    gtdb_representative,
    checkm_completeness,
    checkm_contamination,
    checkm_strain_heterogeneity,
    contig_count,
    genome_size,
    longest_contig,
    n50_contigs,
    protein_count,
    gc_percentage
  )
gtdb_metadata <- bind_rows(bac_metadata, ar_metadata)
gtdb_metadata <- gtdb_metadata[match(out_file_links$gtdb_genome, gtdb_metadata$accession), ]


### write final output ####
msg("writing output...")
write_tsv(
  out_file_links,
  opt$output,
  col_names = FALSE
)
write_tsv(
  gtdb_metadata,
  opt$'metadata-output',
  col_names = TRUE
)
cat(" done\n")

# summarizing stats for missing genomes
found <- nrow(gtdb_links)
total <- nrow(gtdb_taxonomy)
msg(
  paste0("Found ", comma(found), "/", comma(total), " genomes\n\n")
)
