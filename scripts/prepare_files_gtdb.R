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
  "scales"
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
    c("-a", "--ar_metadata"),
    type = "character", default = NULL,
    help = "GTDB Archaea metadata", 
    metavar = "character"
  ),
  make_option(
    c("-b", "--bac_metadata"),
    type = "character", 
    default = NULL, 
    help = "GTDB Bacteria metadata", 
    metavar = "character"
  ),
  make_option(
    c("-A", "--ar_taxonomy"),
    type = "character", 
    default = NULL, 
    help = "GTDB Archaea taxonomy",
    metavar = "character"
  ),
  make_option(
    c("-B", "--bac_taxonomy"),
    type = "character",
    default = NULL, 
    help = "GTDB Bacteria taxonomy", 
    metavar = "character"
  ),
  make_option(
    c("-c", "--gtdb_clusters"),
    type = "character", 
    default = NULL, 
    help = "GTDB clusters",
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
    c("-o", "--outdir"),
    type = "character",
    default = NULL, 
    help = "Output folder", 
    metavar = "character"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$ar_metadata) | is.null(opt$bac_metadata) |
    is.null(opt$ar_taxonomy) | is.null(opt$bac_taxonomy) |
    is.null(opt$gtdb_clusters) | is.null(opt$refseq) |
    is.null(opt$genbank) | is.null(opt$outdir)) {
  print_help(opt_parser)
  stop("All parameters are mandatory.\n", call. = FALSE)
}


### Get necessary files ####

# Get GTDB clusters
msg("Reading GTDB clusters...")
sp_clusters <- read_tsv(
  file = opt$gtdb_clusters, 
  col_names = TRUE, 
  na = "N/A", 
  col_types = cols()
) %>% 
  select(1, 10) %>%
  setNames(c("rep_genome", "cl_genomes")) %>%
  separate_rows(sep = ",", cl_genomes) %>%
  select(cl_genomes) %>%
  separate(
    col = cl_genomes,
    into = c("class", "acc"), 
    extra = "merge", 
    remove = FALSE, 
    fill = "left"
  ) %>%
  mutate(class = ifelse(grepl("UBA", cl_genomes), "UBA", class)) %>%
  separate(
    acc,
    into = "acc_short",
    sep = "\\.", 
    remove = FALSE, 
    extra = "drop"
  )
cat(" done\n")

# Get GTDB metadata
msg("Reading GTDB metadata...")
bac_metadata <- read_tsv(file = opt$bac_metadata, col_names = TRUE, na = "N/A", col_types = cols())
arc_metadata <- read_tsv(file = opt$ar_metadata, col_names = TRUE, na = "N/A", col_types = cols())
gtdb_metadata <- bind_rows(bac_metadata, arc_metadata)
cat(" done\n")

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
gtdb_taxonomy <- bind_rows(bac_taxonomy, arc_taxonomy)
cat(" done\n")

# Report some numbers
msg(paste("Entries in GTDB taxonomy:", comma(nrow(gtdb_taxonomy)), "\n"))
msg(paste("Entries in GTDB metadata:", comma(nrow(gtdb_taxonomy)), "\n"))
msg(paste(comma(nrow(gtdb_taxonomy %>% filter(class == "UBA"))), "are UBA entries\n"))

# Get RefSeq data
msg("Reading RefSeq assembly summaries...")
refseq_metadata <- read_tsv(
  file = opt$refseq, 
  col_names = FALSE, 
  na = "na",
  col_types = cols()
) %>% 
  setNames(c("acc1", "acc2", "link"))
cat(" done\n")

# Get GenBank data
msg("Reading GenBank assembly summaries...")
genbank_metadata <- read_tsv(
  file = opt$genbank,
  col_names = FALSE, 
  na = "na",
  col_types = cols()
) %>% 
  setNames(c("acc1", "acc2", "link"))
cat(" done\n")

# parse NCBI metadata
ncbi_metadata <- bind_rows(genbank_metadata, refseq_metadata) %>%
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
msg(paste("Entries in RefSeq file:", comma(nrow(refseq_metadata)), "\n"))
msg(paste("Entries in GenBank file:", comma(nrow(genbank_metadata)), "\n"))


### Find overlaps between GTDB and NCBI ####
msg("Identifying common genomes between GTDB and RefSeq/GenBank data...")
gtdb_ncbi_tmp <- bind_rows(
  ncbi_metadata %>% inner_join(gtdb_taxonomy, by = c(acc1_short = "acc_short")),
  ncbi_metadata %>% inner_join(gtdb_taxonomy, by = c(acc2_short = "acc_short"))
) %>%
  arrange(gtdb_genome) %>%
  separate(
    acc,
    into = "acc_short", 
    sep = "\\.", 
    remove = FALSE,
    extra = "drop"
  )
cat(" done\n")

# Find duplicated entries, because mix between genbank and refseq, we will keep the ones from GTDB
gtdb_ncbi_tmp_dup <- gtdb_ncbi_tmp %>%
  group_by(gtdb_genome) %>%
  count() %>%
  filter(n > 1)
gtdb_ncbi_tmp_sng <- gtdb_ncbi_tmp %>%
  group_by(gtdb_genome) %>%
  count() %>%
  filter(n == 1)
msg(paste("We found", comma(nrow(gtdb_ncbi_tmp_dup)), "duplicated entries between GTDB/NCBI files\n"))
msg(paste("We found", comma(nrow(gtdb_ncbi_tmp_sng)), "singleton entries between GTDB/NCBI files\n"))

msg("Filtering out duplicates...")
# Find those that the GTDB accession is the same as the primary accession
gtdb_ncbi_tmp_1 <- gtdb_ncbi_tmp %>%
  filter(gtdb_genome %in% gtdb_ncbi_tmp_dup$gtdb_genome) %>%
  filter(acc_short == acc1_short)

# Find the singletons
gtdb_ncbi_tmp_2 <- gtdb_ncbi_tmp %>%
  filter(gtdb_genome %in% gtdb_ncbi_tmp_sng$gtdb_genome)

# check if each duplicate has 1 representative (empty character expected)
diff_acc <- setdiff(gtdb_ncbi_tmp_dup$gtdb_genome, gtdb_ncbi_tmp_1$gtdb_genome) %>%
  unique()

cat(" done\n")

# Which GTDB genomes are not in NCBI
gtdb_ncbi_missing <- gtdb_taxonomy %>%
  filter(class != "UBA") %>%
  filter(!(gtdb_genome %in% c(gtdb_ncbi_tmp_dup$gtdb_genome, gtdb_ncbi_tmp_sng$gtdb_genome))) %>%
  mutate(
    acc_short_mod = case_when(
      grepl("GCF", acc_short) ~ gsub("GCF", "GCA", acc_short),
      grepl("GCA", acc_short) ~ gsub("GCA", "GCF", acc_short)
    )
  )

msg(paste("We are missing", comma(nrow(gtdb_ncbi_missing)), "entries in NCBI (might be updated, removed...)\n"))

# rescue some missing genomes by taking their modified accession number
gtdb_ncbi_tmp_3 <- bind_rows(
  ncbi_metadata %>% inner_join(gtdb_ncbi_missing, by = c(acc1_short = "acc_short_mod")),
  ncbi_metadata %>% inner_join(gtdb_ncbi_missing, by = c(acc2_short = "acc_short_mod"))
) %>%
  arrange(gtdb_genome) %>%
  filter(!is.na(link))

# concatenate download links
gtdb_links <- bind_rows(gtdb_ncbi_tmp_1, gtdb_ncbi_tmp_2, gtdb_ncbi_tmp_3) %>%
  select(gtdb_genome, link, tax_string) %>%
  bind_rows(
    gtdb_taxonomy %>% 
      filter(class == "UBA") %>%
      mutate(link = NA) %>%
      select(gtdb_genome, link, tax_string)
    )

msg(paste("Generating", comma(nrow(gtdb_links %>% filter(!grepl("UBA", gtdb_genome)))), "links for download..."))

### Create necessary files For downloading fasta files ####
out_file_links <- gtdb_links %>%
  # Is there any particular reason for changing from GCF to GCA here?
  # mutate(
  #   link = ifelse(
  #     grepl("GCF_001865635", gtdb_genome), 
  #     "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/865/635/GCA_001865635.3_ASM186563v3",
  #     link
  #   )
  # ) %>%
  filter(!grepl("UBA", gtdb_genome)) %>%
  select(gtdb_genome, link, tax_string) %>%
  mutate(
    link = gsub("^ftp", "http", link),
    filename = paste0(gsub("^.*/", "", link), "_genomic.fna.gz"),
    download_link = paste(link, filename, sep = "/"),
    outfile = paste0(gsub("GB_|RS_", "", gtdb_genome), "_genomic.fna.gz"),
    acc = gsub("GB_|RS_", "", gtdb_genome),
    md5file = paste(link, "md5checksums.txt", sep = "/")
  ) %>%
  select(gtdb_genome, download_link, filename, outfile, md5file, acc, tax_string)

write_tsv(
  out_file_links,
  path = file.path(opt$outdir, "gtdb_download_info.txt"),
  col_names = FALSE
)
cat(" done\n")

msg(paste("Exporting accessions for", comma(nrow(gtdb_links %>% filter(grepl("UBA", gtdb_genome)))), "UBA genomes..."))

out_file_uba <- gtdb_links %>%
  filter(grepl("UBA", gtdb_genome)) %>%
  select(gtdb_genome, tax_string)

write_tsv(
  out_file_uba,
  path = file.path(opt$outdir, "gtdb_uba.txt"),
  col_names = FALSE
)
cat(" done\n\n")

# summarizing stats for missing genomes
found <- nrow(gtdb_links)
total <- nrow(gtdb_taxonomy)
msg(
  paste0("Found ", comma(found), "/", comma(total), " genomes\n\n")
)