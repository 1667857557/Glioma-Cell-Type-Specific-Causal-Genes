suppressMessages({
  pacman::p_load(ggplot2, ggtext, data.table, ggrepel, cowplot, ggsci, ggnewscale)
})

TWAS_plot <- function(twas, sig_z = 4.73, sig_p = NULL, output, width = 8.5, height = 6, res = 600,
                      title = NULL) {
  suppressMessages(library(data.table))
  suppressMessages(library(ggrepel))
  suppressMessages(library(ggplot2))
  suppressMessages(library(cowplot))
  suppressMessages(library(ggsci))
  suppressMessages(library(ggnewscale))
  
  old_scipen <- getOption("scipen")
  options(scipen = 999)
  
  TWAS_manhattan = function(dataframe, title = title,
                            ylimit = max(abs(dataframe$TWAS.Z), na.rm = TRUE) + 1,
                            Sig_Z_Thresh = qnorm(1 - (0.05 / length(dataframe$TWAS.Z)) / 2)) {
    
    d = dataframe[order(CHR, P0), ]
    d = d[!is.na(TWAS.P), ]
    d[, PP.H4.abf := as.numeric(PP.H4.abf)]
    d[, pos := NA_real_]
    ticks = NULL
    lastbase = 0
    sorted_chrs <- sort(unique(d$CHR))
    num_chrs <- length(sorted_chrs)
    base_palette <- pal_jama("default")(6)
    colors <- rep(base_palette, length.out = num_chrs)
    names(colors) <- as.character(sorted_chrs)
    
    for (i in sorted_chrs) {
      if (i == min(sorted_chrs)) {
        d[CHR == i, pos := P0]
      } else {
        prev_max_pos <- max(d[CHR == (i - 1), pos], na.rm = TRUE)
        lastbase = prev_max_pos
        d[CHR == i, pos := P0 + lastbase]
      }
      chrom_pos = d[CHR == i, pos]
      mid_index = floor(length(chrom_pos) / 2) + 1
      ticks = c(ticks, chrom_pos[mid_index])
    }
    
    ticklim = c(min(d$pos), max(d$pos))
    d[, Sig_Z_Thresh := Sig_Z_Thresh]
    d_nonsig <- d[abs(TWAS.Z) <= Sig_Z_Thresh]
    d_sig <- d[abs(TWAS.Z) > Sig_Z_Thresh]
    
    d_sig_min <- d_sig[order(TWAS.P), .SD[1], by = ID]
    d_sig_non_min <- d_sig[!ID %in% d_sig_min$ID]
    d_sig <- rbind(d_sig_min, d_sig_non_min)
    d_sig[, is_min_p := ifelse(ID %in% d_sig_min$ID & TWAS.P == min(TWAS.P), TRUE, FALSE), by = ID]
    
    d_sig[is_min_p == TRUE, status := 'Not duplicate']
    d_sig[duplicate == 'Y' & PP.H4.abf >= 0.7 & Causal == 'Y' & is_min_p == TRUE, 
          status := 'Duplicate, PP.H4.abf ≥ 0.7 & Causal == Y']
    d_sig[duplicate == 'Y' & (PP.H4.abf < 0.7 | Causal != 'Y') & is_min_p == TRUE, 
          status := 'Duplicate, PP.H4.abf < 0.7 or Causal != Y']
    d_sig[is_min_p == FALSE, status := 'No color']
    d_sig[, status := factor(status,
                             levels = c('No color',
                                        'Not duplicate',
                                        'Duplicate, PP.H4.abf < 0.7 or Causal != Y',
                                        'Duplicate, PP.H4.abf ≥ 0.7 & Causal == Y'))]
    
    chr_labs <- as.character(sorted_chrs)
    
    p <- ggplot() +
      geom_point(data = d_nonsig, aes(x = pos, y = TWAS.Z, color = factor(CHR)), shape = 16, size = 0.5, show.legend = FALSE) +
      scale_color_manual(values = colors, guide = "none") +
      new_scale_color() +
      
      geom_point(data = d_sig[status == 'Not duplicate'], 
                 aes(x = pos, y = TWAS.Z), 
                 color = 'gray', 
                 shape = 16, 
                 size = 1,
                 show.legend = FALSE) +
      
      geom_point(data = d_sig[status == 'Duplicate, PP.H4.abf < 0.7 or Causal != Y'], 
                 aes(x = pos, y = TWAS.Z), 
                 color = 'black', 
                 shape = 17, 
                 size = 1,
                 show.legend = FALSE) +
      
      geom_point(data = d_sig[status == 'Duplicate, PP.H4.abf ≥ 0.7 & Causal == Y'], 
                 aes(x = pos, y = TWAS.Z), 
                 color = 'red', 
                 shape = 17, 
                 size = 1,
                 show.legend = FALSE) +
      
      scale_x_continuous(name = "Chromosome", breaks = ticks, labels = chr_labs) +
      scale_y_continuous(name = "Z score", limits = c(-ylimit, ylimit)) +
      geom_hline(yintercept = 0, colour = "black") +
      geom_hline(yintercept = Sig_Z_Thresh, colour = "#00A1D5FF") +
      geom_hline(yintercept = -Sig_Z_Thresh, colour = "#00A1D5FF") +
      geom_text_repel(data = d_sig[duplicate == 'Y' & TWAS.Z > 0 & is_min_p == TRUE],
                      aes(x = pos, y = TWAS.Z, label = ID),
                      colour = 'black',
                      nudge_y = 1, size = 2.5, force = 5,
                      segment.alpha = 0.25,
                      segment.size = 0.2,
                      ylim = c(Sig_Z_Thresh + 0.1, NA)) +
      
      geom_text_repel(data = d_sig[duplicate == 'Y' & TWAS.Z < 0 & is_min_p == TRUE],
                      aes(x = pos, y = TWAS.Z, label = ID),
                      colour = 'black',
                      nudge_y = -1, size = 2.5, force = 5,
                      segment.alpha = 0.25,
                      segment.size = 0.2,
                      ylim = c(NA, -Sig_Z_Thresh - 0.1)) +
      theme_cowplot() +
      theme(
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA),
        axis.text.x = element_text(angle = 45, size = 8, hjust = 1),
        plot.title = element_text(hjust = 0.5)
      ) +
      labs(title = title)
    
    options(scipen = old_scipen)
    p
  }
  
  if (!is.data.table(twas)) {
    twas <- as.data.table(twas)
  }
  
  
  if (!is.na(sig_z)) {
    pdf(file = paste0(output, '.pdf'), width = width / 600 * res, height = height / 600 * res, bg = "white")
    print(TWAS_manhattan(dataframe = twas, Sig_Z_Thresh = sig_z, title = title))
    dev.off()
  }
}
setwd("D:/Work/Ten/TWAS图")
A<-fread("Inhibitory.neurons_GBM.txt")
str(A)
setwd("G:/")
TWAS_plot(twas = A, output = "Inhibitory.neurons_GBM_TWAS"
          ,title = "Inhibitory Neuron-Specific Glioblastomagenesis Causal Genes")
setwd("D:/Work/Ten/TWAS图")
A<-fread("Astrocyte_GBM.txt")
str(A)
setwd("G:/")
TWAS_plot(twas = A, output = "Astrocyte_GBM_TWAS"
          ,title = "Astrocyte-Specific Glioblastomagenesis Causal Genes")
setwd("D:/Work/Ten/TWAS图")
A<-fread("Oligodendrocytes_GBM.txt")
str(A)
setwd("G:/")
TWAS_plot(twas = A, output = "Oligodendrocyte_GBM_TWAS"
          ,title = "Oligodendrocyte-Specific Glioblastomagenesis Causal Genes")
setwd("D:/Work/Ten/TWAS图")
A<-fread("Excitatory.neurons_GBM.txt")
str(A)
setwd("G:/")
TWAS_plot(twas = A, output = "Excitatory.neurons_GBM_TWAS"
          ,title = "Excitatory Neuron-Specific Glioblastomagenesis Causal Genes")
setwd("D:/Work/Ten/TWAS图")
A<-fread("Endothelial.cells_GBM.txt")
str(A)
setwd("G:/")
TWAS_plot(twas = A, output = "Endothelial.cell_GBM_TWAS"
          ,title = "Endothelial-cell-Specific Glioblastomagenesis Causal Genes")

setwd("D:/Work/Ten/TWAS图")
A<-fread("Microglia_GBM.txt")
str(A)
setwd("G:/")
TWAS_plot(twas = A, output = "Microglia_GBM_TWAS"
          ,title = "Microglia-Specific Glioblastomagenesis Causal Genes")
setwd("D:/Work/Ten/TWAS图")
A<-fread("OPC_GBM.txt")
str(A)
setwd("G:/")

TWAS_plot(twas = A, output = "OPC_GBM_TWAS"
          ,title = "OPC-Specific Glioblastomagenesis Causal Genes")

setwd("D:/Work/Ten/TWAS图")
A<-fread("Pericytes_GBM.txt")
str(A)
setwd("G:/")

TWAS_plot(twas = A, output = "Pericytes_GBM_TWAS"
          ,title = "Pericytes-Specific Glioblastomagenesis Causal Genes")

#####figure legend-----

pdf("图例.pdf", width = 11, height = 8)
df <- data.frame(
  x = rnorm(30),
  y = rnorm(30),
  category = factor(rep(c("Cell-type-specific Causal Gene",
                          "Potential Cell-type-specific Causal Gene",
                          "Unduplicated in QTL-base Association Study"), each = 10))
)

legend_labels <- c(
  "Cell-type-specific Causal Gene\n(Duplicate in QTL-base Association Study; PP.H4.abf ≥ 70%; Causal Significant)",
  "Potential Cell-type-specific Causal Gene\n(Duplicate in QTL-base Association Study; PP.H4.abf < 70%)",
  "Unduplicated in QTL-base Association Study"
)
legend_colors <- c("red", "black", "grey")
legend_shapes <- c(17, 17, 16) 
p<-ggplot(df, aes(x = x, y = y, color = category, shape = category)) +
  geom_point(size = 3) +
  scale_color_manual(values = legend_colors, labels = legend_labels) +
  scale_shape_manual(values = legend_shapes, labels = legend_labels) +
  theme_bw() +
  theme(
    legend.position = "right",                      
    legend.title = element_blank(),                 
    legend.text = element_text(size = 10),        
    legend.spacing.y = unit(1, 'cm'),             
    legend.key.height = unit(1.2, "cm")           
  )
ggsave("my_plot_highres.pdf", plot = p, width = 8, height = 6, units = "in", dpi = 300)
