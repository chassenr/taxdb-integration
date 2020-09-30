# calculating sensitivity and specificity stats for taxdb validation
# setwd("/home/chh/Documents/Projects/NCBI_taxdb_integration/Simulation")
setwd("C:/Users/chassenrueck/Documents/Bioinf_projects/NCBI_taxdb_integration/Simulation")
# save.image("sim_db_stats.Rdata")
# load("sim_db_stats.Rdata")

require(data.table)
require(tidyverse)
require(tidytable)
require(taxonomizr)
require(caret)

# # some file prep in bash (faster than R)
# # taxonomy of input sequences (same for v1 and v2 since identical input)
# cut -f2 sim_v1.kraken | sed 's/__SEQ[0-9]\+$//' | sort -t$'\t' | uniq -c | sed -e 's/^\s*//' -e 's/ /\t/' > sim_tax_in_freq.txt
# cut -f2 sim_v1_short.kraken | sed 's/__SEQ[0-9]\+$//' | sort -t$'\t' | uniq -c | sed -e 's/^\s*//' -e 's/ /\t/' > sim_short_tax_in_freq.txt
# 
# # db1 stats 2x150bp
# grep '^C' sim_v1.kraken > sim_v1_kraken_classified.txt
# grep '^U' sim_v1.kraken > sim_v1_kraken_unclassified.txt
# cut -f2 sim_v1_kraken_unclassified.txt | sed 's/__SEQ[0-9]\+$//' | sort -t$'\t' | uniq -c | sed -e 's/^\s*//' -e 's/ /\t/' > sim_v1_unclass_freq.txt
# 
# # db1 stats 1x80bp
# grep '^C' sim_v1_short.kraken > sim_v1_short_kraken_classified.txt
# grep '^U' sim_v1_short.kraken > sim_v1_short_kraken_unclassified.txt
# cut -f2 sim_v1_short_kraken_unclassified.txt | sed 's/__SEQ[0-9]\+$//' | sort -t$'\t' | uniq -c | sed -e 's/^\s*//' -e 's/ /\t/' > sim_v1_short_unclass_freq.txt
# 
# # db2 stats 2x150bp
# grep '^C' sim_v2.kraken > sim_v2_kraken_classified.txt
# grep '^U' sim_v2.kraken > sim_v2_kraken_unclassified.txt
# cut -f2 sim_v2_kraken_unclassified.txt | sed 's/__SEQ[0-9]\+$//' | sort -t$'\t' | uniq -c | sed -e 's/^\s*//' -e 's/ /\t/' > sim_v2_unclass_freq.txt
# 
# # db2 stats 1x80bp
# grep '^C' sim_v2_short.kraken > sim_v2_short_kraken_classified.txt
# grep '^U' sim_v2_short.kraken > sim_v2_short_kraken_unclassified.txt
# cut -f2 sim_v2_short_kraken_unclassified.txt | sed 's/__SEQ[0-9]\+$//' | sort -t$'\t' | uniq -c | sed -e 's/^\s*//' -e 's/ /\t/' > sim_v2_short_unclass_freq.txt


# read tax files
taxaNodes <- read.nodes.sql("nodes.dmp", sqlFile = "taxdump.sqlite", overwrite=TRUE)
taxaNames <- read.names.sql("names.dmp", sqlFile = "taxdump.sqlite", overwrite=TRUE)
tax_coarse_res_db <- read.table("db_subset_macro_taxonomy.txt", h = F, sep = "\t")[, 2] %>% 
  gsub(";s__.*$", "", .) %>% 
  gsub(";", "_", .) %>% 
  unique(.)

### results for 2x150bp simulation ####

file_tax_freq_in <- "sim_tax_in_freq.txt"
file_list_class_out <- c("sim_v1_kraken_classified.txt", "sim_v2_kraken_classified.txt")
file_list_unclass_freq <- c("sim_V1_unclass_freq.txt", "sim_V2_unclass_freq.txt")
names_scenarios <- c("no_mem_lim", "mem_lim_200G")

### results for 2x150bp simulation ####
# for 80bp sim metaG
file_tax_freq_in <- "sim_short_tax_in_freq.txt"
file_list_class_out <- c("sim_v1_short_kraken_classified.txt", "sim_v2_short_kraken_classified.txt")
file_list_unclass_freq <- c("sim_v1_short_unclass_freq.txt", "sim_v2_short_unclass_freq.txt")
names_scenarios <- c("no_mem_lim", "mem_lim_200G")

### to be converted to a function later on ####

# read input taxon frequencies
tax_in_freq <- read.table(file_tax_freq_in, h = F) %>% 
  rename("in_reads" = "V1", "tax_path" = "V2") %>% 
  # add information about tax domain and db resolution
  mutate(
    domain = gsub("d__", "", gsub("_p__.*", "", tax_path)),
    db_res = ifelse(
      grepl(paste(tax_coarse_res_db, collapse = "|"), tax_path),
      "coarse", # one genome per class
      "high" # one genome per species (95% ANI for viruses)
    )
  ) %>% 
  relocate(in_reads, .after = db_res)

# proportion of reads per input taxon that were not classified
unclass_stats <- map(
  1:length(file_list_unclass_freq),
  function(X) {
    
    # read number of unclassified reads per taxon
    unclass_freq <- read.table(file_list_unclass_freq[X], h = F, stringsAsFactors = F)
    
    # calculate percentage of unclassified reads per taxon compared to number of input reads
    test <- unclass_freq %>% 
      as_tibble() %>% 
      rename("unclass_reads" = "V1", "tax_path" = "V2") %>% 
      full_join(tax_in_freq, by = "tax_path") %>% 
      relocate(unclass_reads, .after = in_reads) %>% 
      mutate(
        unclass_reads = replace_na(unclass_reads, replace = 0),
        unclass_perc = unclass_reads/in_reads * 100
      ) %>% 
      arrange(tax_path) 
    
  }
)
names(unclass_stats) <- names_scenarios
View(unclass_stats$no_mem_lim)

# quick and dirty first look at completely unclassified sequences (no kmer hits at all)
par(mfrow = c(1, length(file_list_unclass_freq)))
map(1:length(file_list_unclass_freq), function(X) {
  boxplot(
    unclass_stats[[X]]$unclass_perc ~ unclass_stats[[X]]$domain,
    las = 2,
    ylab = "% unclassified reads",
    xlab = "",
    main = names_scenarios[X]
  )
})
map(1:length(file_list_unclass_freq), function(X) {
  boxplot(
    unclass_stats[[X]]$unclass_perc ~ unclass_stats[[X]]$db_res,
    las = 2,
    ylab = "% unclassified reads",
    xlab = "",
    main = names_scenarios[X]
  )
})
unclass_stats_summary_domain <- map_dfr(unclass_stats, function(X) {
  group_by(X, domain) %>% 
    summarise(sum_unclass = sum(unclass_reads)/sum(in_reads) * 100)
}) %>% 
  mutate(
    scenario = rep(names_scenarios, each = 4),
    .before = domain
  )
unclass_stats_summary_db_res <- map_dfr(unclass_stats, function(X) {
  group_by(X, db_res) %>% 
    summarise(sum_unclass = sum(unclass_reads)/sum(in_reads) * 100)
}) %>% 
  mutate(
    scenario = rep(names_scenarios, each = 2),
    .before = db_res
  )
# less than 0.35% unclassified for any taxon
# in total 0.5 of euk reads unclassified (0.003% for bacteria, and 0.159 for viruses)

# for 80bp: 2-3 times more reads per input genome not classified in mem limited verion
# in general 1 order of magnitude more unclassified reads
# considerable more unclassified reads for non-microbial genomes (coarse)

# correct assignment for classified reads 
# (i.e. those with at least 1 kmer match, may still include those unclassified according to LCA)

# reading classified reads
class_reads <- map(
  1:length(file_list_class_out),
  function(X) {
    fread(file_list_class_out[X], header = FALSE, fill = TRUE, sep = "\t")[, 2:3] %>% 
      rename("tax_path_in" = "V2", "taxid_out" = "V3")
  }
)
names(class_reads) <- names_scenarios

# adding taxonomy to classified reads
# convert any non alpha-numeric character to underscore
# (this is done for input taxonomic path by simulation program)
tax_out_class <-  map(
  class_reads,
  function(X) {
    getTaxonomy(X$taxid_out, taxaNodes, taxaNames) %>% 
      gsub("[^[:alnum:].-]", "_", .)
  }
)
names(tax_out_class) <- names_scenarios

# splitting input taxonomic path
tax_in_class <- map(
  class_reads,
  function(X) {
    X$tax_path_in %>% 
      gsub("__SEQ[0-9]*$", "", .) %>%
      gsub("^d__", "", .) %>%
      as_tibble() %>%
      separate(
        col = "value",
        into = c("superkingdom", "phylum", "class", "order", "family", "genus", "species"),
        sep = "_[pcofgs]__"
      )
  }
)
names(tax_in_class) <- names_scenarios

# generate summary for input taxonomic paths
tax_in_class_summary <- map(
  class_reads,
  function(X) {
    X$tax_path_in %>% 
      gsub("__SEQ[0-9]*$", "", .) %>%
      as_tibble() %>%
      rename("tax_path_in" = "value") %>% 
      mutate(
        db_res = ifelse(tax_path_in %in% tax_in_freq$tax_path[tax_in_freq$db_res == "coarse"], "coarse", "high"),
        domain = gsub("d__", "", gsub("_p__.*", "", tax_path_in))
      )
  }
)
names(tax_in_class_summary) <- names_scenarios

# for each taxonomic level, how many sequences were assigned to correct taxon?
tax_comparison <- map(
  1:length(names_scenarios),
  function(X) {
    tax_in_class[[X]] == tax_out_class[[X]]
  }
)
names(tax_comparison) <- names_scenarios

# overall stats
class_stats_summary_domain <- map(
  1:length(names_scenarios),
  function(X) {
    apply(tax_comparison[[X]], 2, function(x) {
      c(by(x, tax_in_class_summary[[X]]$domain, function(y) sum(y, na.rm = T)/length(y) * 100))
    }) %>% 
      as.data.frame() %>% 
      rownames_to_column(var = "domain")
  }
) %>% 
  do.call("rbind", .) %>% 
  mutate(
    scenario = rep(names_scenarios, each = 4),
    .before = domain
  )
class_stats_summary_db_res <- map(
  1:length(names_scenarios),
  function(X) {
    apply(tax_comparison[[X]], 2, function(x) {
      c(by(x, tax_in_class_summary[[X]]$db_res, function(y) sum(y, na.rm = T)/length(y) * 100))
    }) %>% 
      as.data.frame() %>% 
      rownames_to_column(var = "db_res")
  }
) %>% 
  do.call("rbind", .) %>% 
  mutate(
    scenario = rep(names_scenarios, each = 2),
    .before = db_res
  )

# per taxonomic path summary
class_stats_summary_path <- map(
  1:length(names_scenarios),
  function(X) {
    apply(tax_comparison[[X]], 2, function(x) {
      c(by(x, tax_in_class_summary[[X]]$tax_path_in, function(y) sum(y, na.rm = T)/length(y) * 100))
    }) %>% 
      as.data.frame() %>% 
      rownames_to_column(var = "tax_path_in") %>% 
      left_join(., tax_in_freq[, -4], by = c("tax_path_in" = "tax_path"))
  }
)
names(class_stats_summary_path) <- names_scenarios
View(class_stats_summary_path$no_mem_lim)

# quick and dirty first look
par(mfrow = c(1, length(names_scenarios)))
for(i in 2:8) {
  map(1:length(names_scenarios), function(X) {
    boxplot(
      class_stats_summary_path[[X]][, i] ~ class_stats_summary_path[[X]]$db_res,
      las = 2,
      ylab = paste0("% correctly assigned reads per genome, tax level: ", colnames(class_stats_summary_path[[X]])[i]),
      xlab = "",
      main = names_scenarios[X],
      ylim = c(0, 100)
    )
  })
}

# which input path (i.e. genome) has high specificity in coarse DB?
tax_coarse_res_db_reps <- read.table("db_subset_class_macro_taxonomy.txt", h = F, sep = "\t")[, 2]
class_stats_summary_path$no_mem_lim %>% 
  filter(
    db_res == "coarse",
    species > 60
  )
grep("Priapulus", tax_coarse_res_db_reps, value = T) # species present (but different genome)
grep("Hypsibius", tax_coarse_res_db_reps, value = T) # species present (but different genome)
grep("Hofstenia", tax_coarse_res_db_reps, value = T) # species present (but different genome)

# which input path (i.e. genome) has low specificity in high res DB?
tax_high_res_db_reps <- c(
  read.table("db_subset_species_micro_taxonomy_nv.txt", h = F, sep = "\t")[, 2],
  read.table("db_subset_micro_taxonomy_vi.txt", h = F, sep = "\t")[, 2]
) %>% unique()
class_stats_summary_path$no_mem_lim %>% 
  filter(
    db_res == "high",
    phylum < 40
  )
grep("g__Bacillariophyta", tax_high_res_db_reps, value = T) 
# species not defined (s__Bacillariophyta sp.)
# since only 1 representative picked per species (as defined based on path), the db rep is too distantly related to classify seqs
grep("Galdieria", tax_high_res_db_reps, value = T)
# species present (but different genome)
table(
  tax_out_class$no_mem_lim[tax_in_class_summary$no_mem_lim$tax_path_in == "d__Eukaryota_p__Rhodophyta_c__Bangiophyceae_o__Cyanidiales_f__Cyanidiaceae_g__Galdieria_s__Galdieria_sulphuraria", "phylum"],
  useNA = "always"
)
# very weird assignment...
# mash distance between sim and db genome of same species is 0.233808
# ANI (fastani) is 78.1384
# none of these genome is refseq
# I suspect that the one of them (or both) may be of bad quality
class_stats_summary_path$no_mem_lim %>% 
  filter(
    db_res == "high",
    species < 70,
    phylum > 70
  )
# missclassification on species level only
# may be better with more extensive db

# confusion matrix
conf_m_domain_high_res <- map(
  1:length(names_scenarios),
  function(X) {
    tmp_out <- as.matrix(tax_out_class[[X]][tax_in_class_summary[[X]]$db_res == "high", ])
    tmp_in <- as.matrix(tax_in_class[[X]][tax_in_class_summary[[X]]$db_res == "high", ])
    
    confusionMatrix(
      factor(tmp_out[, 1], levels = unique(tmp_in[, 1])),
      factor(tmp_in[, 1], levels = unique(tmp_in[, 1]))
    )
  }
)
lapply(conf_m_domain_high_res, function(x) round(prop.table(as.matrix(x$table), 2) * 100, 2))

conf_m_domain_coarse <- map(
  1:length(names_scenarios),
  function(X) {
    tmp_out <- as.matrix(tax_out_class[[X]][tax_in_class_summary[[X]]$db_res == "coarse", ])
    data.frame(table(tmp_out[, 1]), row.names = 1)
  }
)
lapply(conf_m_domain_coarse, function(x) round(prop.table(as.matrix(x), 2) * 100, 2))

# per genome % classified between the 2 DB version
par(mfrow = c(4, 2))
for(i in 2:8) {
  plot(
    class_stats_summary_path$no_mem_lim[, i],
    class_stats_summary_path$mem_lim_200G[, i],
    col = as.numeric(as.factor(class_stats_summary_path$no_mem_lim$domain))
  )
  abline(0, 1)
}
