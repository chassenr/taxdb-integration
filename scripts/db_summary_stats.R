# reading database summary files
setwd("/home/chh/Documents/Projects/NCBI_taxdb_integration/Summary_stats")
require(tidyverse)
require(R.utils)

# number of taxa at each level
files <- list.files(pattern = "_tax_stats.txt")
tax_stats <- map_dfr(files, function(X) {
  read.table(X, h = F, sep = "\t", stringsAsFactors = F) %>% 
    rename(taxlevel = V1, N = V2) %>% 
    mutate(partition = gsub("_tax_stats\\.txt", "", X), .before = "taxlevel")
})
# number of taxa without dereplication
files <- list.files(pattern = "_tax_all.txt")
tax_all <- map_dfr(files, function(X) {
  read.table(X, h = F, sep = "\t", stringsAsFactors = F) %>% 
    rename(taxlevel = V1, N = V2) %>% 
    mutate(partition = gsub("_tax_all.txt", "", X), .before = "taxlevel")
})

plot(
  0, 0,
  xlim = c(0.5, 7.5),
  ylim = c(0, log10(max(c(tax_stats$N, tax_all$N)))),
  type = "n",
  xlab = "",
  ylab = "Number of taxa [log10]",
  axes = F
)
axis(1, at = 1:7, labels = tax_stats$taxlevel[1:7], las = 2)
axis(2)
box("plot")
for(i in 1:length(unique(tax_stats$partition))) {
  lines(
    1:7,
    log10(tax_stats$N[tax_stats$partition == unique(tax_stats$partition)[i]]),
    col = rainbow(length(unique(tax_stats$partition)))[i]
  )
  lines(
    1:7,
    log10(tax_all$N[tax_all$partition == unique(tax_stats$partition)[i]]),
    col = rainbow(length(unique(tax_stats$partition)))[i],
    lty = 3
  )
}
legend(
  "bottomright",
  legend = unique(tax_stats$partition),
  col = rainbow(length(unique(tax_stats$partition))),
  pch = 15,
  cex = 0.7,
  pt.cex = 1
)
# withou GTDB and viruses
tax_stats_sub <- tax_stats[!tax_stats$partition %in% c("gtdb", "viral", "fungi"), ]
tax_all_sub <- tax_all[!tax_all$partition %in% c("gtdb", "viral", "fungi"), ]
plot(
  0, 0,
  xlim = c(0.5, 7.5),
  ylim = c(0, max(c(tax_stats_sub$N, tax_all_sub$N))),
  type = "n",
  xlab = "",
  ylab = "Number of taxa",
  axes = F
)
axis(1, at = 1:7, labels = tax_stats_sub$taxlevel[1:7], las = 2)
axis(2)
box("plot")
for(i in 1:length(unique(tax_stats_sub$partition))) {
  lines(
    1:7,
    tax_stats_sub$N[tax_stats_sub$partition == unique(tax_stats_sub$partition)[i]],
    col = rainbow(length(unique(tax_stats_sub$partition)))[i]
  )
  lines(
    1:7,
    tax_all_sub$N[tax_all_sub$partition == unique(tax_stats_sub$partition)[i]],
    col = rainbow(length(unique(tax_stats_sub$partition)))[i],
    lty = 3
  )
}
legend(
  "topleft",
  legend = unique(tax_stats_sub$partition),
  col = rainbow(length(unique(tax_stats_sub$partition))),
  pch = 15,
  cex = 0.7,
  pt.cex = 1
)


# database size if only 1 representative is picked
files <- list.files(pattern = "_size_stats.txt")
size_stats <- map_dfr(files, function(X) {
  read.table(X, h = F, sep = "\t", stringsAsFactors = F) %>% 
    rename(taxlevel = V1, size = V2) %>% 
    mutate(partition = gsub("_size_stats\\.txt", "", X), .before = "taxlevel")
})

plot(
  0, 0,
  xlim = c(0.5, 7.5),
  ylim = c(0, max(size_stats$size)),
  type = "n",
  xlab = "",
  ylab = "File size",
  axes = F
)
axis(1, at = 1:7, labels = size_stats$taxlevel[1:7], las = 2)
axis(2, at = seq(0, max(size_stats$size), length.out = 5), labels = hsize(seq(0, max(size_stats$size), length.out = 5)))
box("plot")
for(i in 1:length(unique(size_stats$partition))) {
  lines(
    1:7,
    (size_stats$size[size_stats$partition == unique(size_stats$partition)[i]]),
    col = rainbow(length(unique(size_stats$partition)))[i]
  )
}
legend(
  "topleft",
  legend = unique(size_stats$partition),
  col = rainbow(length(unique(size_stats$partition))),
  pch = 15,
  cex = 0.7,
  pt.cex = 1
)

hsize(sum(size_stats$size[size_stats$taxlevel == "class" & !size_stats$partition %in% c("gtdb", "viral", "fungi")]))
# 1.2TB family, 511GB order, 183GB class
hsize(sum(size_stats$size[size_stats$taxlevel == "genomes" & size_stats$partition %in% c("gtdb", "viral", "fungi")]))
# 285GB
# even at family level, size of database is too large
# on order level, should be ok with 1TB RAM
# but I fear that too many reads would end up as unclassified 
# (or worse, be wrongly assigned to different partition)
# a large number of taxa would be excluded
# 
