# inspect coarse benchmark results
require(tidyverse)
require(reshape)
require(dlfUtils)

format_axis <- function(param, label, xmar = 0, yoffset = 0.1) {
  tmp_axis <- c(1, which(diff(param) != 0) + 1, length(param) + 1)
  axis(
    1,
    at = tmp_axis,
    labels = F,
    tcl = -0.5,
    mgp = c(3, xmar + 0.5, xmar)
  )
  text(
    apply(data.frame(tmp_axis[-length(tmp_axis)], tmp_axis[-1]), 1, mean),
    rep(line2user(line = xmar, side = 1, outer = F), length(tmp_axis) - 1),
    labels = param[tmp_axis][-length(tmp_axis)],
    cex = 0.5,
    pos = 1,
    offset = 0.2,
    adj = 0.5
  )
  if(!is.na(label)) {
    text(
      par("usr")[2] + (par("usr")[2] - par("usr")[1]) * yoffset,
      line2user(line = xmar, side = 1, outer = F),
      labels = label,
      cex = 0.8,
      pos = 1,
      offset = 0.2,
      adj = 0
    )
  }
}

sim_summarize <- function(simpath, ground_truth) {
  file_list_taxlevels_filt <- list.files(path = simpath, pattern = "_summary_taxlevels_filt.txt")
  file_list_taxlevels_nofilt <- list.files(path = simpath, pattern = "_summary_taxlevels_nofilt.txt")

  # read filtered data
  taxlevels_filt <- map_dfr(
    file_list_taxlevels_filt,
    function(X) {
      tmp <- read.table(
        file.path(simpath, X),
        h = T,
        stringsAsFactors = F,
        sep = "\t"
      )
      scenario <- strsplit(X, "_")[[1]][c(1:7, 10, 11)]
      tmp$db <- scenario[1]
      tmp$kmer <- as.numeric(gsub("k", "", scenario[2]))
      tmp$minimizer <- as.numeric(gsub("m", "", scenario[3]))
      tmp$mspaces <- as.numeric(gsub("s", "", scenario[4]))
      tmp$memory <- scenario[5]
      tmp$conf_highres <- as.numeric(gsub("c", "", scenario[6]))
      tmp$rtl_highres <- as.numeric(gsub("rtl", "", scenario[7]))
      tmp$conf_coarse <- as.numeric(gsub("c", "", scenario[8]))
      tmp$rtl_coarse <- as.numeric(gsub("rtl", "", scenario[9]))
      return(tmp)
    }
  )
  taxlevels_filt$taxlevel <- factor(taxlevels_filt$taxlevel, levels = c("superkingdom", "phylum", "class", "order", "family", "genus", "species"))
  taxlevels_filt$unassigned[is.na(taxlevels_filt$unassigned)] <- 0
  
  # read data without coarse filter
  taxlevels_nofilt <- map_dfr(
    file_list_taxlevels_nofilt,
    function(X) {
      tmp <- read.table(
        file.path(simpath, X),
        h = T,
        stringsAsFactors = F,
        sep = "\t"
      )
      scenario <- strsplit(X, "_")[[1]][c(1:7, 10, 11)]
      tmp$db <- scenario[1]
      tmp$kmer <- as.numeric(gsub("k", "", scenario[2]))
      tmp$minimizer <- as.numeric(gsub("m", "", scenario[3]))
      tmp$mspaces <- as.numeric(gsub("s", "", scenario[4]))
      tmp$memory <- scenario[5]
      tmp$conf_highres <- as.numeric(gsub("c", "", scenario[6]))
      tmp$rtl_highres <- as.numeric(gsub("rtl", "", scenario[7]))
      tmp$conf_coarse <- as.numeric(gsub("c", "", scenario[8]))
      tmp$rtl_coarse <- as.numeric(gsub("rtl", "", scenario[9]))
      return(tmp)
    }
  )
  taxlevels_nofilt$taxlevel <- factor(taxlevels_nofilt$taxlevel, levels = c("superkingdom", "phylum", "class", "order", "family", "genus", "species"))
  taxlevels_nofilt$unassigned[is.na(taxlevels_nofilt$unassigned)] <- 0
  
  # show absolute read numbers
  pdf(file.path(simpath, paste0(basename(simpath), "_taxlevels_abs.pdf")), width = 10, height = 12, onefile = T)
  par(ann = F, mfrow = c(length(levels(taxlevels_filt$taxlevel)), 2), mar = rep(0.5, 4), oma = c(4, 5, 4, 6), xpd = NA)
  for(k in unique(taxlevels_filt$kmer)) {
    filt_sub <- taxlevels_filt[taxlevels_filt$kmer == k, ]
    nofilt_sub <- taxlevels_nofilt[taxlevels_nofilt$kmer == k, ]
    for(m in unique(filt_sub$memory)) {
      filt_sub_sub <- filt_sub[filt_sub$memory == m, ]
      nofilt_sub_sub <- nofilt_sub[nofilt_sub$memory == m, ]
      for(j in levels(filt_sub_sub$taxlevel)) {
        filt_sub_taxlevel <- filt_sub_sub[filt_sub_sub$taxlevel == j, ] %>% 
          arrange(conf_coarse, rtl_coarse, conf_highres, rtl_highres) %>% 
          unite(col = plot_group, conf_coarse, rtl_coarse, conf_highres, remove = F)
        nofilt_sub_taxlevel <- nofilt_sub_sub[nofilt_sub_sub$taxlevel == j, ] %>% 
          arrange(conf_coarse, rtl_coarse, conf_highres, rtl_highres) %>% 
          unite(col = plot_group, conf_coarse, rtl_coarse, conf_highres, remove = F)
        
        # plot absolute seq counts (filtered)
        plot(
          0, 0, 
          type = "n",
          xlim = c(1, nrow(filt_sub_taxlevel)),
          # ylim = c(
          #   min(c(taxlevels_filt$correct, taxlevels_nofilt$correct)), 
          #   max(c(apply(taxlevels_filt[, 1:4], 1, sum), apply(taxlevels_nofilt[, 1:4], 1, sum)))
          # ),
          ylim = c(
            0, 
            max(c(apply(taxlevels_filt[, 1:4], 1, sum), apply(taxlevels_nofilt[, 1:4], 1, sum)))
          ),
          axes = F
        )
        par(xpd = F)
        abline(h = ground_truth, col = "green")
        par(xpd = NA)
        for(i in unique(filt_sub_taxlevel$plot_group)) {
          lines(
            which(filt_sub_taxlevel$plot_group == i),
            (filt_sub_taxlevel$correct + filt_sub_taxlevel$incorrect_pro + filt_sub_taxlevel$incorrect_euk + filt_sub_taxlevel$unassigned)[filt_sub_taxlevel$plot_group == i],
            col = "grey"
          )
          lines(
            which(filt_sub_taxlevel$plot_group == i),
            (filt_sub_taxlevel$correct + filt_sub_taxlevel$incorrect_pro + filt_sub_taxlevel$incorrect_euk)[filt_sub_taxlevel$plot_group == i],
            col = "red"
          )
          lines(
            which(filt_sub_taxlevel$plot_group == i),
            (filt_sub_taxlevel$correct + filt_sub_taxlevel$incorrect_pro)[filt_sub_taxlevel$plot_group == i],
            col = "blue"
          )
          lines(
            which(filt_sub_taxlevel$plot_group == i),
            filt_sub_taxlevel$correct[filt_sub_taxlevel$plot_group == i],
            col = "black"
          )
        }
        box("plot")
        axis(2, las = 2, mgp = c(2, 0.5, 0), cex.axis = 0.8)
        title(ylab = j, mgp = c(2.5, 0.5, 0))
        if(j == "superkingdom") {
          title(main = "With coarse filter", line = 1, font.main = 1)
        }
        if(j == "species") {
          text(
            c(which.min(filt_sub_taxlevel$rtl_highres), which.max(filt_sub_taxlevel$rtl_highres)),
            rep(line2user(line = 0, side = 1, outer = F), 2),
            labels = c(min(filt_sub_taxlevel$rtl_highres), max(filt_sub_taxlevel$rtl_highres)),
            cex = 0.5,
            pos = 3,
            offset = 0.2
          )
          format_axis(param = filt_sub_taxlevel$conf_highres, label = NA, xmar = 0)
          format_axis(param = filt_sub_taxlevel$rtl_coarse, label = NA, xmar = 1)
          format_axis(param = filt_sub_taxlevel$conf_coarse, label = NA, xmar = 2)
        }
        
        # plot absolute seq counts (without coarse filtering)
        plot(
          0, 0, 
          type = "n",
          xlim = c(1, nrow(nofilt_sub_taxlevel)),
          # ylim = c(
          #   min(c(taxlevels_filt$correct, taxlevels_nofilt$correct)), 
          #   max(c(apply(taxlevels_filt[, 1:4], 1, sum), apply(taxlevels_nofilt[, 1:4], 1, sum)))
          # ),
          ylim = c(
            0, 
            max(c(apply(taxlevels_filt[, 1:4], 1, sum), apply(taxlevels_nofilt[, 1:4], 1, sum)))
          ),
          axes = F
        )
        par(xpd = F)
        abline(h = ground_truth, col = "green")
        par(xpd = NA)
        for(i in unique(nofilt_sub_taxlevel$plot_group)) {
          lines(
            which(nofilt_sub_taxlevel$plot_group == i),
            (nofilt_sub_taxlevel$correct + nofilt_sub_taxlevel$incorrect_pro + nofilt_sub_taxlevel$incorrect_euk + nofilt_sub_taxlevel$unassigned)[nofilt_sub_taxlevel$plot_group == i],
            col = "grey"
          )
          lines(
            which(nofilt_sub_taxlevel$plot_group == i),
            (nofilt_sub_taxlevel$correct + nofilt_sub_taxlevel$incorrect_pro + nofilt_sub_taxlevel$incorrect_euk)[nofilt_sub_taxlevel$plot_group == i],
            col = "red"
          )
          lines(
            which(nofilt_sub_taxlevel$plot_group == i),
            (nofilt_sub_taxlevel$correct + nofilt_sub_taxlevel$incorrect_pro)[nofilt_sub_taxlevel$plot_group == i],
            col = "blue"
          )
          lines(
            which(nofilt_sub_taxlevel$plot_group == i),
            nofilt_sub_taxlevel$correct[nofilt_sub_taxlevel$plot_group == i],
            col = "black"
          )
        }
        box("plot")
        axis(2, labels = NA, mgp = c(2, 0.5, 0))
        if(j == "superkingdom") {
          title(main = "No coarse filter", line = 1, font.main = 1)
        }
        if(j == "species") {
          text(
            c(which.min(filt_sub_taxlevel$rtl_highres), which.max(filt_sub_taxlevel$rtl_highres)),
            rep(line2user(line = 0, side = 1, outer = F), 2),
            labels = c(min(filt_sub_taxlevel$rtl_highres), max(filt_sub_taxlevel$rtl_highres)),
            cex = 0.5,
            pos = 3,
            offset = 0.2
          )
          text(
            par("usr")[2] + (par("usr")[2] - par("usr")[1]) * 0.1,
            line2user(line = 0, side = 1, outer = F),
            labels = "RTL highres",
            cex = 0.8,
            pos = 3,
            offset = 0.2,
            adj = 0
          )
          format_axis(param = filt_sub_taxlevel$conf_highres, label = "Confidence highres", xmar = 0)
          format_axis(param = filt_sub_taxlevel$rtl_coarse, label = "RTL coarse", xmar = 1)
          format_axis(param = filt_sub_taxlevel$conf_coarse, label = "Confidence coarse", xmar = 2)
          legend(
            x = par("usr")[2],
            y = par("usr")[4],
            legend = c("expected", "correct", "incorrect (pro)", "incorrect (euk)", "unassigned"),
            lty = 1,
            col = c("green", "black", "blue", "red", "grey"),
            cex = 0.8,
            bty = "n"
          )
        }
      }
      mtext("Number of reads", side = 2, line = 3.5, outer = T)
      mtext(paste(m, "kmer", k), side = 3, line = 2, outer = T)
    }
  }
  dev.off()
  
  # show percentages
  taxlevels_filt_rel <- taxlevels_filt
  taxlevels_filt_rel[, 1:4] <- taxlevels_filt[, 1:4]/rowSums(taxlevels_filt[, 1:4]) * 100
  taxlevels_nofilt_rel <- taxlevels_nofilt
  taxlevels_nofilt_rel[, 1:4] <- taxlevels_nofilt[, 1:4]/rowSums(taxlevels_nofilt[, 1:4]) * 100
  
  pdf(file.path(simpath, paste0(basename(simpath), "_taxlevels_rel.pdf")), width = 10, height = 12, onefile = T)
  par(ann = F, mfrow = c(length(levels(taxlevels_filt$taxlevel)), 2), mar = rep(0.5, 4), oma = c(4, 5, 4, 6), xpd = NA)
  for(k in unique(taxlevels_filt$kmer)) {
    filt_sub <- taxlevels_filt_rel[taxlevels_filt$kmer == k, ]
    nofilt_sub <- taxlevels_nofilt_rel[taxlevels_nofilt$kmer == k, ]
    for(m in unique(filt_sub$memory)) {
      filt_sub_sub <- filt_sub[filt_sub$memory == m, ]
      nofilt_sub_sub <- nofilt_sub[nofilt_sub$memory == m, ]
      for(j in levels(filt_sub_sub$taxlevel)) {
        filt_sub_taxlevel <- filt_sub_sub[filt_sub_sub$taxlevel == j, ] %>% 
          arrange(conf_coarse, rtl_coarse, conf_highres, rtl_highres) %>% 
          unite(col = plot_group, conf_coarse, rtl_coarse, conf_highres, remove = F)
        nofilt_sub_taxlevel <- nofilt_sub_sub[nofilt_sub_sub$taxlevel == j, ] %>% 
          arrange(conf_coarse, rtl_coarse, conf_highres, rtl_highres) %>% 
          unite(col = plot_group, conf_coarse, rtl_coarse, conf_highres, remove = F)
        
        # plot absolute seq counts (filtered)
        plot(
          0, 0, 
          type = "n",
          xlim = c(1, nrow(filt_sub_taxlevel)),
          ylim = c(0, 100),
          axes = F
        )
        for(i in unique(filt_sub_taxlevel$plot_group)) {
          lines(
            which(nofilt_sub_taxlevel$plot_group == i),
            (nofilt_sub_taxlevel$correct + nofilt_sub_taxlevel$incorrect_pro + nofilt_sub_taxlevel$incorrect_euk + nofilt_sub_taxlevel$unassigned)[nofilt_sub_taxlevel$plot_group == i],
            col = "grey"
          )
          lines(
            which(filt_sub_taxlevel$plot_group == i),
            (filt_sub_taxlevel$correct + filt_sub_taxlevel$incorrect_pro + filt_sub_taxlevel$incorrect_euk)[filt_sub_taxlevel$plot_group == i],
            col = "red"
          )
          lines(
            which(filt_sub_taxlevel$plot_group == i),
            (filt_sub_taxlevel$correct + filt_sub_taxlevel$incorrect_pro)[filt_sub_taxlevel$plot_group == i],
            col = "blue"
          )
          lines(
            which(filt_sub_taxlevel$plot_group == i),
            filt_sub_taxlevel$correct[filt_sub_taxlevel$plot_group == i],
            col = "black"
          )
        }
        box("plot")
        axis(2, las = 2, mgp = c(2, 0.5, 0), cex.axis = 0.8)
        title(ylab = j, mgp = c(2.5, 0.5, 0))
        if(j == "superkingdom") {
          title(main = "With coarse filter", line = 1, font.main = 1)
        }
        if(j == "species") {
          text(
            c(which.min(filt_sub_taxlevel$rtl_highres), which.max(filt_sub_taxlevel$rtl_highres)),
            rep(line2user(line = 0, side = 1, outer = F), 2),
            labels = c(min(filt_sub_taxlevel$rtl_highres), max(filt_sub_taxlevel$rtl_highres)),
            cex = 0.5,
            pos = 3,
            offset = 0.2
          )
          format_axis(param = filt_sub_taxlevel$conf_highres, label = NA, xmar = 0)
          format_axis(param = filt_sub_taxlevel$rtl_coarse, label = NA, xmar = 1)
          format_axis(param = filt_sub_taxlevel$conf_coarse, label = NA, xmar = 2)
        }
        
        # plot absolute seq counts (without coarse filtering)
        plot(
          0, 0, 
          type = "n",
          xlim = c(1, nrow(nofilt_sub_taxlevel)),
          ylim = c(0, 100),
          axes = F
        )
        for(i in unique(nofilt_sub_taxlevel$plot_group)) {
          lines(
            which(nofilt_sub_taxlevel$plot_group == i),
            (nofilt_sub_taxlevel$correct + nofilt_sub_taxlevel$incorrect_pro + nofilt_sub_taxlevel$incorrect_euk + nofilt_sub_taxlevel$unassigned)[nofilt_sub_taxlevel$plot_group == i],
            col = "grey"
          )
          lines(
            which(nofilt_sub_taxlevel$plot_group == i),
            (nofilt_sub_taxlevel$correct + nofilt_sub_taxlevel$incorrect_pro + nofilt_sub_taxlevel$incorrect_euk)[nofilt_sub_taxlevel$plot_group == i],
            col = "red"
          )
          lines(
            which(nofilt_sub_taxlevel$plot_group == i),
            (nofilt_sub_taxlevel$correct + nofilt_sub_taxlevel$incorrect_pro)[nofilt_sub_taxlevel$plot_group == i],
            col = "blue"
          )
          lines(
            which(nofilt_sub_taxlevel$plot_group == i),
            nofilt_sub_taxlevel$correct[nofilt_sub_taxlevel$plot_group == i],
            col = "black"
          )
        }
        box("plot")
        axis(2, labels = NA, mgp = c(2, 0.5, 0))
        if(j == "superkingdom") {
          title(main = "No coarse filter", line = 1, font.main = 1)
        }
        if(j == "species") {
          text(
            c(which.min(filt_sub_taxlevel$rtl_highres), which.max(filt_sub_taxlevel$rtl_highres)),
            rep(line2user(line = 0, side = 1, outer = F), 2),
            labels = c(min(filt_sub_taxlevel$rtl_highres), max(filt_sub_taxlevel$rtl_highres)),
            cex = 0.5,
            pos = 3,
            offset = 0.2
          )
          text(
            par("usr")[2] + (par("usr")[2] - par("usr")[1]) * 0.1,
            line2user(line = 0, side = 1, outer = F),
            labels = "RTL highres",
            cex = 0.8,
            pos = 3,
            offset = 0.2,
            adj = 0
          )
          format_axis(param = filt_sub_taxlevel$conf_highres, label = "Confidence highres", xmar = 0)
          format_axis(param = filt_sub_taxlevel$rtl_coarse, label = "RTL coarse", xmar = 1)
          format_axis(param = filt_sub_taxlevel$conf_coarse, label = "Confidence coarse", xmar = 2)
          legend(
            x = par("usr")[2],
            y = par("usr")[4],
            legend = c("correct", "incorrect (pro)", "incorrect (euk)", "unassigned"),
            lty = 1,
            col = c("black", "blue", "red", "grey"),
            cex = 0.8,
            bty = "n"
          )
        }
      }
      mtext("Percentage of reads", side = 2, line = 3.5, outer = T)
      mtext(paste(m, "kmer", k), side = 3, line = 2, outer = T)
    }
  }
  dev.off()
  
  return(
    list(
      taxlevels_filt,
      taxlevels_nofilt
    )
  )
}

ground_truth_SE_40 <- 5937901
ground_truth_SE_80 <- 5894762
ground_truth_PE_150 <- 5812590

SE_40_highres_pro <- sim_summarize("/storage/hdd1/chh/TaxDB/Test1/simulation/sim_out/sim_SE_40_out", ground_truth_SE_40)
SE_80_highres_pro <- sim_summarize("/storage/hdd1/chh/TaxDB/Test1/simulation/sim_out/sim_SE_80_out", ground_truth_SE_80)
PE_150_highres_pro <- sim_summarize("/storage/hdd1/chh/TaxDB/Test1/simulation/sim_out/sim_PE_150_out", ground_truth_PE_150)
