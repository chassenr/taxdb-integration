# inspect coarse benchmark results
require(tidyverse)
require(reshape)

sim_summarize <- function(simpath) {
  file_list_cm_domain <- list.files(path = simpath, pattern = "cm_domain.txt")
  file_list_stats_domain <- list.files(path = simpath, pattern = "stats_domain.txt")
  file_list_cm_taxongroups <- list.files(path = simpath, pattern = "cm_taxongroups.txt")
  file_list_stats_taxongroups <- list.files(path = simpath, pattern = "stats_taxongroups.txt")
  
  
  ### by domain ####
  cm_domain <- map_dfr(
    file_list_cm_domain,
    function(X) {
      tmp <- read.table(
        file.path(simpath, X),
        h = T,
        stringsAsFactors = F,
        sep = "\t",
        na.strings = ""
      )
      scenario <- strsplit(X, "_")[[1]][1:7]
      tmp$db <- scenario[1]
      tmp$kmer <- as.numeric(gsub("k", "", scenario[2]))
      tmp$minimizer <- as.numeric(gsub("m", "", scenario[3]))
      tmp$mspaces <- as.numeric(gsub("s", "", scenario[4]))
      tmp$memory <- scenario[5]
      tmp$confidence <- as.numeric(gsub("c", "", scenario[6]))
      tmp$rtl <- as.numeric(gsub("rtl", "", scenario[7]))
      return(tmp)
    }
  )
  cm_domain$Prediction <- factor(cm_domain$Prediction, levels = c("Archaea", "Bacteria", "Eukaryota", "Viruses", "NA"))
  cm_domain$Reference <- factor(cm_domain$Reference, levels = c("Archaea", "Bacteria", "Eukaryota", "Viruses", "NA"))
  cm_domain$mem_kmer <- paste(cm_domain$memory, cm_domain$kmer, sep = "_")
  
  pdf(file.path(simpath, "cm_domain.pdf"), width = 10, height = 9)
  par(
    ann = F,
    mfcol = c(
      length(unique(sort(cm_domain$rtl))), 
      length(unique(sort(paste(cm_domain$mem_kmer, cm_domain$confidence))))
    ), 
    mar = rep(0.5, 4), 
    oma = c(3, 5, 7, 3), 
    xpd = NA
  )
  for(i in unique(sort(cm_domain$mem_kmer))) {
    for(j in unique(sort(cm_domain$confidence))) {
      for(k in unique(sort(cm_domain$rtl))) {
        tmp <- cm_domain[cm_domain$mem_kmer == i & cm_domain$confidence == j & cm_domain$rtl == k, ]
        tmp_table <- cast(tmp, "Prediction ~ Reference", value = "Freq", fun.aggregate = sum)
        rownames(tmp_table) <- tmp_table$Prediction
        tmp_table <- data.matrix(tmp_table[, -1])
        tmp_perc <- prop.table(tmp_table, 2) * 100
        image(
          x = 1:ncol(tmp_perc),
          y = 1:nrow(tmp_perc),
          z = sqrt(t(tmp_perc[nrow(tmp_perc):1, ])),
          col = colorRampPalette(c("grey95", "yellow", "orange", "red", "darkred"))(50),
          zlim = c(0, 10),
          axes = F
        )
        box("plot")
        text(
          1:ncol(tmp_perc),
          nrow(tmp_perc):1,
          labels = round(diag(tmp_perc)),
          font = 2,
          col = "white"
        )
        if(k == min(unique(sort(cm_domain$rtl)))) {
          axis(3, at = 1:ncol(tmp_perc), labels = levels(cm_domain$Reference), las = 2, tcl = -0.3, mgp = c(2, 0.5, 0), cex.axis = 0.8)
          title(main = unique(tmp$memory), font.main = 1, line = 6)
          title(main = paste("kmer length", unique(tmp$kmer)), font.main = 1, line = 5)
          title(main = paste("confidence", j), font.main = 1, line = 4)
        }
        if(k == max(unique(sort(cm_domain$rtl)))) {
          axis(1, at = 1:ncol(tmp_perc), labels = levels(cm_domain$Reference), las = 2, tcl = -0.3, mgp = c(2, 0.5, 0), cex.axis = 0.8)
        }
        if(i == unique(sort(cm_domain$mem_kmer))[1] & j == min(unique(sort(cm_domain$confidence)))) {
          axis(2, at = nrow(tmp_perc):1, labels = levels(cm_domain$Reference), las = 2, tcl = -0.3, mgp = c(2, 0.5, 0), cex.axis = 0.8)
          title(ylab = paste("RTL", k), mgp = c(4, 0.5, 0))
        }
        if(i == unique(sort(cm_domain$mem_kmer))[length(unique(sort(cm_domain$mem_kmer)))] & j == max(unique(sort(cm_domain$confidence)))) {
          axis(4, at = nrow(tmp_perc):1, labels = levels(cm_domain$Reference), las = 2, tcl = -0.3, mgp = c(2, 0.5, 0), cex.axis = 0.8)
        }
      }
    }
  }
  dev.off()
  
  stats_domain <- map_dfr(
    file_list_stats_domain,
    function(X) {
      tmp <- read.table(
        file.path(simpath, X),
        h = T,
        stringsAsFactors = F,
        sep = "\t",
        na.strings = ""
      )
      tmp$class <- gsub("Class: ", "", rownames(tmp))
      rownames(tmp) <- NULL
      scenario <- strsplit(X, "_")[[1]][1:7]
      tmp$db <- scenario[1]
      tmp$kmer <- as.numeric(gsub("k", "", scenario[2]))
      tmp$minimizer <- as.numeric(gsub("m", "", scenario[3]))
      tmp$mspaces <- as.numeric(gsub("s", "", scenario[4]))
      tmp$memory <- scenario[5]
      tmp$confidence <- as.numeric(gsub("c", "", scenario[6]))
      tmp$rtl <- as.numeric(gsub("rtl", "", scenario[7]))
      return(tmp)
    }
  )
  stats_domain$class <- factor(stats_domain$class, levels = c("Archaea", "Bacteria", "Eukaryota", "Viruses", "NA"))
  stats_domain$mem_kmer <- paste(stats_domain$memory, stats_domain$kmer, sep = "_")
  
  pdf(file.path(simpath, "stats_domain.pdf"), width = 10, height = 9)
  par(
    ann = F,
    mfcol = c(
      length(unique(sort(stats_domain$rtl))), 
      length(unique(sort(paste(stats_domain$mem_kmer, stats_domain$confidence))))
    ), 
    mar = rep(0.5, 4),
    oma = c(3, 3, 7, 3),
    xpd = NA
  )
  for(i in unique(sort(stats_domain$mem_kmer))) {
    for(j in unique(sort(stats_domain$confidence))) {
      for(k in unique(sort(stats_domain$rtl))) {
        tmp <- stats_domain[stats_domain$mem_kmer == i & stats_domain$confidence == j & stats_domain$rtl == k, ]
        plot(
          0, 0, 
          type = "n",
          xlim = c(0.7, length(levels(stats_domain$class)) + 0.3),
          ylim = c(0, 1),
          axes = F
        )
        box("plot")
        points(
          as.numeric(tmp$class),
          tmp$Specificity,
          pch = 3,
          col = "blue",
          lwd = 2
        )
        points(
          as.numeric(tmp$class),
          tmp$Sensitivity,
          pch = 1,
          col = "red",
          lwd = 2
        )
        points(
          as.numeric(tmp$class),
          tmp$F1,
          pch = 5,
          col = "green",
          lwd = 2
        )
        if(k == min(unique(sort(stats_domain$rtl)))) {
          axis(3, at = 1:length(levels(stats_domain$class)), labels = levels(stats_domain$class), las = 2, tcl = -0.3, mgp = c(2, 0.5, 0), cex.axis = 0.8)
          title(main = unique(tmp$memory), font.main = 1, line = 6)
          title(main = paste("kmer length", unique(tmp$kmer)), font.main = 1, line = 5)
          title(main = paste("confidence", j), font.main = 1, line = 4)
        }
        if(k == max(unique(sort(stats_domain$rtl)))) {
          axis(1, at = 1:length(levels(stats_domain$class)), labels = levels(stats_domain$class), las = 2, tcl = -0.3, mgp = c(2, 0.5, 0), cex.axis = 0.8)
        }
        if(i == unique(sort(stats_domain$mem_kmer))[1] & j == min(unique(sort(stats_domain$confidence)))) {
          axis(2, las = 2, tcl = -0.3, mgp = c(2, 0.5, 0), cex.axis = 0.8)
          title(ylab = paste("RTL", k), mgp = c(2, 0.5, 0))
        }
        if(i == unique(sort(stats_domain$mem_kmer))[length(unique(sort(stats_domain$mem_kmer)))] & j == max(unique(sort(stats_domain$confidence)))) {
          axis(4, las = 2, tcl = -0.3, mgp = c(2, 0.5, 0), cex.axis = 0.8)
        }
      }
    }
  }
  dev.off()
  
  #####
  
  
  ### by larger taxonomic groups ####
  
  cm_taxongroups <- map_dfr(
    file_list_cm_taxongroups,
    function(X) {
      tmp <- read.table(
        file.path(simpath, X),
        h = T,
        stringsAsFactors = F,
        sep = "\t",
        na.strings = ""
      )
      scenario <- strsplit(X, "_")[[1]][1:7]
      tmp$db <- scenario[1]
      tmp$kmer <- as.numeric(gsub("k", "", scenario[2]))
      tmp$minimizer <- as.numeric(gsub("m", "", scenario[3]))
      tmp$mspaces <- as.numeric(gsub("s", "", scenario[4]))
      tmp$memory <- scenario[5]
      tmp$confidence <- as.numeric(gsub("c", "", scenario[6]))
      tmp$rtl <- as.numeric(gsub("rtl", "", scenario[7]))
      return(tmp)
    }
  )
  cm_taxongroups$Prediction <- factor(cm_taxongroups$Prediction, levels = c("Archaea", "Bacteria", "Viruses", "Fungi", "Protists", "Plants", "Metazoa", "NA"))
  cm_taxongroups$Reference <- factor(cm_taxongroups$Reference, levels = c("Archaea", "Bacteria", "Viruses", "Fungi", "Protists", "Plants", "Metazoa", "NA"))
  cm_taxongroups$mem_kmer <- paste(cm_taxongroups$memory, cm_taxongroups$kmer, sep = "_")
  
  pdf(file.path(simpath, "cm_taxongroups.pdf"), width = 10, height = 9)
  par(
    ann = F, 
    mfcol = c(
      length(unique(sort(cm_taxongroups$rtl))), 
      length(unique(sort(paste(cm_taxongroups$mem_kmer, cm_taxongroups$confidence))))
    ), 
    mar = rep(0.5, 4), 
    oma = c(3, 5, 7, 3),
    xpd = NA
  )
  for(i in unique(sort(cm_taxongroups$mem_kmer))) {
    for(j in unique(sort(cm_taxongroups$confidence))) {
      for(k in unique(sort(cm_taxongroups$rtl))) {
        tmp <- cm_taxongroups[cm_taxongroups$mem_kmer == i & cm_taxongroups$confidence == j & cm_taxongroups$rtl == k, ]
        tmp_table <- cast(tmp, "Prediction ~ Reference", value = "Freq", fun.aggregate = sum)
        rownames(tmp_table) <- tmp_table$Prediction
        tmp_table <- data.matrix(tmp_table[, -1])
        tmp_perc <- prop.table(tmp_table, 2) * 100
        image(
          x = 1:ncol(tmp_perc),
          y = 1:nrow(tmp_perc),
          z = sqrt(t(tmp_perc[nrow(tmp_perc):1, ])),
          col = colorRampPalette(c("grey95", "yellow", "orange", "red", "darkred"))(50),
          zlim = c(0, 10),
          axes = F
        )
        box("plot")
        text(
          1:ncol(tmp_perc),
          nrow(tmp_perc):1,
          labels = round(diag(tmp_perc)),
          font = 2,
          col = "white"
        )
        if(k == min(unique(sort(cm_taxongroups$rtl)))) {
          axis(3, at = 1:ncol(tmp_perc), labels = levels(cm_taxongroups$Reference), las = 2, tcl = -0.3, mgp = c(2, 0.5, 0), cex.axis = 0.8)
          title(main = unique(tmp$memory), font.main = 1, line = 6)
          title(main = paste("kmer length", unique(tmp$kmer)), font.main = 1, line = 5)
          title(main = paste("confidence", j), font.main = 1, line = 4)
        }
        if(k == max(unique(sort(cm_taxongroups$rtl)))) {
          axis(1, at = 1:ncol(tmp_perc), labels = levels(cm_taxongroups$Reference), las = 2, tcl = -0.3, mgp = c(2, 0.5, 0), cex.axis = 0.8)
        }
        if(i == unique(sort(cm_taxongroups$mem_kmer))[1] & j == min(unique(sort(cm_taxongroups$confidence)))) {
          axis(2, at = nrow(tmp_perc):1, labels = levels(cm_taxongroups$Reference), las = 2, tcl = -0.3, mgp = c(2, 0.5, 0), cex.axis = 0.8)
          title(ylab = paste("RTL", k), mgp = c(4, 0.5, 0))
        }
        if(i == unique(sort(cm_taxongroups$mem_kmer))[length(unique(sort(cm_taxongroups$mem_kmer)))] & j == max(unique(sort(cm_taxongroups$confidence)))) {
          axis(4, at = nrow(tmp_perc):1, labels = levels(cm_taxongroups$Reference), las = 2, tcl = -0.3, mgp = c(2, 0.5, 0), cex.axis = 0.8)
        }
      }
    }
  }
  dev.off()
  
  stats_taxongroups <- map_dfr(
    file_list_stats_taxongroups,
    function(X) {
      tmp <- read.table(
        file.path(simpath, X),
        h = T,
        stringsAsFactors = F,
        sep = "\t",
        na.strings = ""
      )
      tmp$class <- gsub("Class: ", "", rownames(tmp))
      rownames(tmp) <- NULL
      scenario <- strsplit(X, "_")[[1]][1:7]
      tmp$db <- scenario[1]
      tmp$kmer <- as.numeric(gsub("k", "", scenario[2]))
      tmp$minimizer <- as.numeric(gsub("m", "", scenario[3]))
      tmp$mspaces <- as.numeric(gsub("s", "", scenario[4]))
      tmp$memory <- scenario[5]
      tmp$confidence <- as.numeric(gsub("c", "", scenario[6]))
      tmp$rtl <- as.numeric(gsub("rtl", "", scenario[7]))
      return(tmp)
    }
  )
  stats_taxongroups$class <- factor(stats_taxongroups$class, levels = c("Archaea", "Bacteria", "Viruses", "Fungi", "Protists", "Plants", "Metazoa", "NA"))
  stats_taxongroups$mem_kmer <- paste(stats_taxongroups$memory, stats_taxongroups$kmer, sep = "_")
  
  pdf(file.path(simpath, "stats_taxongroups.pdf"), width = 10, height = 9)
  par(
    ann = F, 
    mfcol = c(
      length(unique(sort(stats_taxongroups$rtl))), 
      length(unique(sort(paste(stats_taxongroups$mem_kmer, stats_taxongroups$confidence))))
    ),
    mar = rep(0.5, 4), 
    oma = c(3, 3, 7, 3), 
    xpd = NA
  )
  for(i in unique(sort(stats_taxongroups$mem_kmer))) {
    for(j in unique(sort(stats_taxongroups$confidence))) {
      for(k in unique(sort(stats_taxongroups$rtl))) {
        tmp <- stats_taxongroups[stats_taxongroups$mem_kmer == i & stats_taxongroups$confidence == j & stats_taxongroups$rtl == k, ]
        plot(
          0, 0, 
          type = "n",
          xlim = c(0.7, length(levels(stats_taxongroups$class)) + 0.3),
          ylim = c(0, 1),
          axes = F
        )
        box("plot")
        points(
          as.numeric(tmp$class),
          tmp$Specificity,
          pch = 3,
          col = "blue",
          lwd = 2
        )
        points(
          as.numeric(tmp$class),
          tmp$Sensitivity,
          pch = 1,
          col = "red",
          lwd = 2
        )
        points(
          as.numeric(tmp$class),
          tmp$F1,
          pch = 5,
          col = "green",
          lwd = 2
        )
        if(k == min(unique(sort(stats_taxongroups$rtl)))) {
          axis(3, at = 1:length(levels(stats_taxongroups$class)), labels = levels(stats_taxongroups$class), las = 2, tcl = -0.3, mgp = c(2, 0.5, 0), cex.axis = 0.8)
          title(main = unique(tmp$memory), font.main = 1, line = 6)
          title(main = paste("kmer length", unique(tmp$kmer)), font.main = 1, line = 5)
          title(main = paste("confidence", j), font.main = 1, line = 4)
        }
        if(k == max(unique(sort(stats_taxongroups$rtl)))) {
          axis(1, at = 1:length(levels(stats_taxongroups$class)), labels = levels(stats_taxongroups$class), las = 2, tcl = -0.3, mgp = c(2, 0.5, 0), cex.axis = 0.8)
        }
        if(i == unique(sort(stats_taxongroups$mem_kmer))[1] & j == min(unique(sort(stats_taxongroups$confidence)))) {
          axis(2, las = 2, tcl = -0.3, mgp = c(2, 0.5, 0), cex.axis = 0.8)
          title(ylab = paste("RTL", k), mgp = c(2, 0.5, 0))
        }
        if(i == unique(sort(stats_taxongroups$mem_kmer))[length(unique(sort(stats_taxongroups$mem_kmer)))] & j == max(unique(sort(stats_taxongroups$confidence)))) {
          axis(4, las = 2, tcl = -0.3, mgp = c(2, 0.5, 0), cex.axis = 0.8)
        }
      }
    }
  }
  dev.off()
  
  #####
  
  return(
    list(
      cm_domain,
      stats_domain,
      cm_taxongroups,
      stats_taxongroups
    )
  )
}

SE_40 <- sim_summarize("C:/Users/chassenrueck/Documents/Bioinf_projects/NCBI_taxdb_integration/Sim_olorin/SE_40")
SE_80 <- sim_summarize("C:/Users/chassenrueck/Documents/Bioinf_projects/NCBI_taxdb_integration/Sim_olorin/SE_80")
PE_150 <- sim_summarize("C:/Users/chassenrueck/Documents/Bioinf_projects/NCBI_taxdb_integration/Sim_olorin/PE_150")

pdf("C:/Users/chassenrueck/Documents/Bioinf_projects/NCBI_taxdb_integration/Sim_olorin/sim_legend.pdf", width = 5, height = 5)
par(ann = F, mfrow = c(2, 1), mar = c(5, 3, 5, 3), xpd = NA)
image(
  x = seq(0, 10, length.out = 50),
  y = 1,
  z = matrix(seq(0, 10, length.out = 50), ncol = 1, nrow = 50),
  col = colorRampPalette(c("grey95", "yellow", "orange", "red", "darkred"))(50),
  axes = F
)
box("plot")
axis(
  1, 
  at = sqrt(c(0, 1, 5, 10, 20, 30, 50, 70, 100)),
  labels = c(0, 1, 5, 10, 20, 30, 50, 70, 100),
  mgp = c(2, 0.5, 0),
  tcl = -0.3,
  cex.axis = 0.8
)
title(main = "Legend confusion matrix heatmap", font.main = 1)
title(xlab = "Percentage of input reads per domain", mgp = c(2, 0.5, 0))
par(mar = c(5, 3, 7, 3))
plot.new()
legend(
  "top",
  legend = c("Specificity", "Sensitivity", "F1 score"),
  col = c("blue", "red", "green"),
  pch = c(3, 1, 5),
  pt.cex = 1.3,
  pt.lwd = 2,
  bty = "n"
)
title(main = "Legend classification statistics", font.main = 1)
dev.off()
