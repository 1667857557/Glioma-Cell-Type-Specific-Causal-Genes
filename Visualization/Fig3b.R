Bootstrap_P_Barplot <- function(p_results, p_names, title = "Test scPagwas", figurenames = NULL, 
                                width = 5, height = 7, do_plot = TRUE) {
  logp <- -log10(p_results)
  sig <- rep("b", length(p_results))
  sig[which(p_results < 0.05/16)] <- "a"
  gg <- data.frame(logp, sig, p_names)
  if (sum(p_results < 0.05/16) > 0) {
    p1 <- ggplot2::ggplot(gg, aes(x = stats::reorder(p_names, 
                                                     logp), y = logp, fill = sig)) + geom_bar(position = "dodge", 
                                                                                              stat = "identity") + theme_classic() + labs(x = "", 
                                                                                                                                          y = "-log10(P)", title = title) + coord_flip() + geom_hline(aes(yintercept = 2.505), 
                                                                                                                                                                                                      colour = "#990000", linetype = "dashed") + scale_fill_manual(values = c("#BB6464", 
                                                                                                                                                                                                                                                                              "#C3DBD9")) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
  } else {
    p1 <- ggplot2::ggplot(gg, aes(x = stats::reorder(p_names, 
                                                     logp), y = logp)) + geom_bar(position = "dodge", 
                                                                                  stat = "identity", color = "#C3DBD9") + theme_classic() + 
      labs(x = "", y = "-log10(P bonferroni)", title = title) + coord_flip() + 
      theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
  }
  if (do_plot) 
    print(p1)
  if (!is.null(figurenames)) {
    grDevices::pdf(figurenames, width = width, height = height)
    print(p1)
    grDevices::dev.off()
  }
}


A<-vroom("G:/A/cell_specify/scPagwas/GBM_cell_type_GBmap/GBmap_cell_type.txt")
Bootstrap_P_Barplot(p_results=A$pvalue,
                    p_names=A$celltype,
                    title = "GBM Risk Cell",
                    figurenames = "GBM_related_celltypes.pdf",
                    width = 5,
                    height = 3,
                    do_plot=TRUE)
