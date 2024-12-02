query_snp <- function (snp, input_type = "rsID", eqtl_tissue = "Whole_Blood", verbose = FALSE) {
  if (!input_type %in% c("rsID", "hg19", "hg38")) {
    return(message("Please select 1 input type from rsID, hg19 or hg38!"))
  }
  if (verbose) {
    print("Get information from 3D-SNP v2 ...")
  }
  if (input_type == "hg38") {
    snp <- ViSNP::get_snp_rsID(snp, assembly = "hg38")
  } else if (input_type == "hg19") {
    snp <- ViSNP::get_snp_rsID(snp, assembly = "hg19")
  }
  url <- paste("https://www.omic.tech/3dsnpv2/api.do?id=", snp, "&format=json&type=basic,motif,ccre,eqtl,physcores,3dsnp,3dgene,tfbs", sep = "")
  res = httr::GET(url, config = httr::config(ssl_verifypeer = FALSE))
  res_data = jsonlite::fromJSON(rawToChar(res$content))
  rsID <- snp
  location <- res_data$disp_chr_pos[!is.na(res_data$disp_chr_pos)]
  chrom <- stringr::str_split_fixed(location, ":", 2)[1]
  position_hg19 <- as.numeric(stringr::str_split_fixed(location, ":", 2)[2]) + 1
  location_hg19 <- paste(chrom, position_hg19, sep = ":")
  location_hg38 <- ViSNP::get_snp_loc(snp, assembly = "hg38")
  ref_allele <- res_data$alleles_Ref[!is.na(res_data$alleles_Ref)]
  alt_allele <- res_data$alleles_Alt[!is.na(res_data$alleles_Alt)]
  EAS_AF <- as.numeric(res_data$alleles_EAS[1])
  AMR_AF <- as.numeric(res_data$alleles_AMR[1])
  AFR_AF <- as.numeric(res_data$alleles_AFR[1])
  EUR_AF <- as.numeric(res_data$alleles_EUR[1])
  SAS_AF <- as.numeric(res_data$alleles_SAS[1])
  AF_info <- paste(EAS_AF, "(EAS),", AMR_AF, "(AMR),", AFR_AF, "(AFR),", EUR_AF, "(EUR),", SAS_AF, "(SAS)", sep = "")
  linear_closest_gene <- res_data$data_gene[[1]]$name
  linear_closest_gene_id <- res_data$data_gene[[1]]$id
  linear_closest_gene_type <- res_data$data_gene[[1]]$loc
  linear_closest_gene_merged <- c()
  if (is.null(linear_closest_gene)) {
    linear_closest_gene_merged <- "-"
  } else {
    for (i in 1:length(linear_closest_gene)) {
      new <- paste(linear_closest_gene[i], "(", linear_closest_gene_id[i], ")", sep = "")
      linear_closest_gene_merged <- c(linear_closest_gene_merged, new)
    }
  }
  linear_closest_gene_type <- paste(linear_closest_gene_type, collapse = ",")
  if (is.null(linear_closest_gene)) {
    linear_closest_gene <- "-"
    linear_closest_gene_id <- "-"
    linear_closest_gene_type <- "-"
    linear_closest_gene_description <- "-"
  }
  loop_gene_col <- res_data$data_loop_gene[!is.null(res_data$data_loop_gene)]
  for (i in 1:length(loop_gene_col)) {
    if (!is.null(loop_gene_col[[i]])) {
      loop_gene_info <- loop_gene_col[[i]]
      break
    }
  }
  loop_gene <- paste(unique(loop_gene_info$gene)[unique(loop_gene_info$gene) != ""], collapse = ",")
  if (loop_gene == "") {
    loop_gene <- "-"
  }
  loop_snp_col <- res_data$data_loop_snp
  for (i in 1:length(loop_snp_col)) {
    if (!is.null(loop_snp_col[[i]])) {
      loop_snp_info <- loop_snp_col[[i]]
      break
    }
  }
  loop_snps <- unique(loop_snp_info$SNP_B)[unique(loop_snp_info$SNP_B) != ""]
  if (length(loop_snps) != 0) {
    loop_snps <- loop_snps[order(loop_snps)]
    loop_snps <- paste(loop_snps, collapse = ",")
  } else {
    loop_snps <- "-"
  }
  if (loop_snps == "") {
    loop_snps <- "-"
  }
  tfbs_col <- res_data$data_tfbs
  tfbs <- "-"
  if (!is.null(tfbs_col)) {
    for (i in 1:length(tfbs_col)) {
      if (!is.null(tfbs_col[[i]])) {
        tfbs_info <- tfbs_col[[i]]
        break
      }
    }
    tfbs <- unique(tfbs_info$gene)
    tfbs <- paste(tfbs, collapse = ",")
  }
  if (verbose) {
    print("Get information from SCREEN and GWAS catalog...")
  }
  cCRE_table <- ViSNP::get_snp_ccre(snp)
  if (nrow(cCRE_table) == 0) {
    cCRE_info <- "-"
  } else {
    cCRE_info <- paste(cCRE_table$type, " (", cCRE_table$cCRE_accession, ")", sep = "")
  }
  gwas_table <- ViSNP::get_snp_gwas(snp)
  if (nrow(gwas_table) == 0) {
    gwas_info <- "-"
  } else {
    gwas_info <- unique(gwas_table$DISEASE.TRAIT)
    gwas_info <- paste(gwas_info, collapse = "; ")
  }
  if (verbose) {
    cat(paste("rsID:", rsID, "\n"))
    cat(paste("chrom:", chrom, "\n"))
    cat(paste("Location(hg19):", location_hg19, "\n"))
    cat(paste("Location(hg38):", location_hg38, "\n"))
    cat(paste("Ref allele:", ref_allele, "\n"))
    cat(paste("Alt allele:", alt_allele, "\n"))
    cat(paste("AF:", AF_info, "\n"))
    if (length(linear_closest_gene) <= 1) {
      if (linear_closest_gene != "-") {
        cat(paste("Linear closest gene: ", linear_closest_gene_merged, "\n", sep = ""))
      } else {
        cat(paste("Linear closest gene:", linear_closest_gene, "\n", sep = ""))
      }
    } else {
      cat(paste("Linear closest gene:", paste(linear_closest_gene_merged, collapse = ","), "\n", sep = ""))
    }
    cat(paste("Linear closest gene type:", linear_closest_gene_type, "\n"))
    cat(paste("3D Loop Genes:", loop_gene, "\n"))
    cat(paste("3D-Interacting SNPs:", loop_snps, "\n"))
    cat(paste("TFBS:", tfbs, "\n"))
    cat(paste("cCRE:", cCRE_info, "\n"))
    cat(paste("GWAS:", gwas_info, "\n"))
  }
  info_table <- data.frame(
    Info = c("rsID", "Chromosome", "Location(hg19)", "Location(hg38)", "Ref allele", "Alt allele", "AF", "Linear closest gene", "Linear closest gene type", "3D-Loop Genes", "3D-Interacting SNPs", "TFBS", "cCRE", "GWAS"),
    Value = c(rsID, chrom, location_hg19, location_hg38, ref_allele, alt_allele, AF_info, paste(linear_closest_gene_merged, collapse = ","), linear_closest_gene_type, loop_gene, loop_snps, tfbs, cCRE_info, gwas_info)
  )
  return(info_table)
}
plot_snp_circos <- function(
    snp_info_table, 
    user_snp_list = NULL,     
    output_assembly = "hg19", 
    window_size = 2e5, 
    circos_pos = "Relative", 
    savefile = "circos.pdf"
){
  
  library(circlize)
  library(stringr)
  library(ComplexHeatmap)
  library(ViSNP)
  
  snp <- snp_info_table[snp_info_table$Info == "rsID", 2]
  chrom <- snp_info_table[snp_info_table$Info == "Chromosome", 2]
  window_size <- as.numeric(window_size)
  
  if (output_assembly == "hg19"){
    location <- snp_info_table[snp_info_table$Info == "Location(hg19)", 2]
    chrom_info <- chrom_info_hg19
  } else {
    location <- snp_info_table[snp_info_table$Info == "Location(hg38)", 2]
    chrom_info <- chrom_info_hg38
  }
  
  chrom_max_len <- chrom_info[chrom_info$V1 == chrom, 2]
  pos <- as.numeric(str_split_fixed(location, ":", 2)[2])
  
  loop_genes <- snp_info_table[snp_info_table$Info == "3D-Loop Genes", 2]
  n_genes <- str_count(loop_genes, ",")
  loop_genes <- as.vector(str_split_fixed(loop_genes, ",", n_genes + 1))
  
  input_type <- ifelse(output_assembly == "hg38", "hg38", "hg19")
  gene_loc_df <- get_genes_loc(loop_genes, assembly = input_type)
  
  gene_loc_input <- data.frame(
    chrom = rep(chrom, nrow(gene_loc_df)), 
    start_pos = gene_loc_df$start_position, 
    end_pos = gene_loc_df$end_position, 
    gene = gene_loc_df$gene_symbol
  )
  
  loop_snps <- snp_info_table[snp_info_table$Info == "3D-Interacting SNPs", 2]
  n_snps <- str_count(loop_snps, ",")
  loop_snps <- as.vector(str_split_fixed(loop_snps, ",", n_snps + 1))
  
  if (!is.null(user_snp_list)) {
    loop_snps <- intersect(loop_snps, user_snp_list)
  }
  
  loop_snps_filtered <- c()
  snps_loc <- c()
  for (s in loop_snps){
    s_loc <- get_snp_loc(s, assembly = input_type)
    if (s_loc != "NA:NA"){
      snps_loc <- c(snps_loc, s_loc)
      loop_snps_filtered <- c(loop_snps_filtered, s)
    }
  }
    curr_loc_input <- data.frame(
    chrom = chrom, 
    start_pos = pos, 
    end_pos = pos, 
    rsID = snp
  )
  
  if(length(snps_loc) > 0){
    snps_loc_input <- data.frame(
      chrom = str_split_fixed(snps_loc, ":", 2)[,1], 
      start_pos = as.numeric(str_split_fixed(snps_loc, ":", 2)[,2]), 
      end_pos = as.numeric(str_split_fixed(snps_loc, ":", 2)[,2]), 
      rsID = loop_snps_filtered
    )
  } else {
    snps_loc_input <- data.frame(
      chrom = character(0), 
      start_pos = numeric(0), 
      end_pos = numeric(0), 
      rsID = character(0)
    )
  }
  
  snps_loc_input <- rbind(snps_loc_input, curr_loc_input)
  
  start_pos <- max(0, pos - window_size)
  end_pos <- min(pos + window_size, chrom_max_len)
  chrom_input <- data.frame(
    chrom = chrom, 
    start_pos = start_pos, 
    end_pos = end_pos
  )
  
  if (output_assembly == "hg19"){
    cCRE_data <- cCRE_data_hg19
  } else {
    cCRE_data <- cCRE_data_hg38
  }
  
  colnames(cCRE_data) <- c("chrom", "start", "end", "accession", "cCRE_accession", "type")
  cCRE_data_intersect <- cCRE_data[
    (cCRE_data$chrom == chrom) & 
      (cCRE_data$start >= start_pos) & 
      (cCRE_data$end <= end_pos), 
  ]
  
  cCRE_data_group <- cCRE_data_intersect
  cCRE_data_group$type <- gsub(pattern = ",CTCF-bound", "", cCRE_data_intersect$type)
  cCRE_data_input <- cCRE_data_group[, c("chrom", "start", "end", "type")]
  
  cCRE_all_colors <- data.frame(
    type = c("dELS", "pELS", "PLS", "DNase-H3K4me3", "CTCF-only"), 
    cols = c("#FFCD00", "#FFA700", "#FF0000", "#ffaaaa", "#00B0F0")
  )
  
  cCRE_colors <- sapply(cCRE_data_input$type, function(t){
    color <- cCRE_all_colors$cols[cCRE_all_colors$type == t]
    if(length(color) == 0){
      color <- "#000000" 
    }
    return(color)
  })
  
  pdf(file = savefile, width = 10, height = 7)
  
  circos.par("start.degree" = 90)
  
  if (circos_pos == "Absolute"){
    circos.genomicInitialize(
      chrom_input, 
      plotType = 'axis', 
      axis.labels.cex = 0.45, 
      tickLabelsStartFromZero = FALSE
    )
  } else {
    circos.genomicInitialize(
      chrom_input, 
      plotType = 'axis', 
      axis.labels.cex = 0.7, 
      tickLabelsStartFromZero = TRUE
    )
  }
  
  if(nrow(cCRE_data_input) > 0){
    circos.genomicTrackPlotRegion(
      cCRE_data_input, 
      track.height = 0.1, 
      stack = TRUE, 
      bg.border = "gray",
      panel.fun = function(region, value, ...) {
        circos.genomicRect(region, value, col = cCRE_colors, border = NA, ...)
      }
    )
  }
  
  if(nrow(gene_loc_input) > 0){
    gene_cols <- pal_igv(alpha = 0.6)(nrow(gene_loc_input))
    
    circos.genomicTrackPlotRegion(
      gene_loc_input, 
      track.height = 0.1, 
      stack = TRUE, 
      bg.border = "gray",
      panel.fun = function(region, value, ...) {
        circos.genomicRect(region, value, col = gene_cols, border = NA, ...)
      }
    )
  }
  
  if(nrow(snps_loc_input) > 0){
    circos.genomicTrackPlotRegion(
      snps_loc_input, 
      track.height = 0.1, 
      stack = TRUE, 
      bg.border = NA,
      panel.fun = function(region, value, ...) {
        circos.genomicPoints(
          region, value, 
          pch = 18, 
          cex = 1.2, 
          col = ifelse(snps_loc_input$rsID == snp, "#FF1300", "#00C322"),
          border = NA, ...
        )
      }
    )
  }
  
  if (nrow(gene_loc_input) != 0){
    curr_snp <- data.frame(
      chrom = chrom, 
      start_pos = pos, 
      end_pos = pos, 
      rsID = snp, 
      freq = 1:nrow(gene_loc_input)
    )
    rcols <- scales::alpha(gene_cols, alpha = 0.4)
    circos.genomicLink(curr_snp, gene_loc_input, directional = 1, col = rcols)
  }
  
  if(nrow(snps_loc_input) > 0){
    circos.genomicLabels(
      snps_loc_input, 
      labels.column = 4, 
      side = "inside", 
      connection_height = 0.2, 
      padding = 3,
      col = ifelse(snps_loc_input$rsID == snp, "#FF1300", "#009393"), 
      cex = 0.8, 
      labels_height = 0.25, 
      niceFacing = TRUE
    )
  }
  
  if(nrow(cCRE_data_input) > 0){
    legend_labels <- unique(cCRE_data_input$type)
    at_colors <- sapply(legend_labels, function(t){
      color <- cCRE_all_colors$cols[cCRE_all_colors$type == t]
      if(length(color) == 0){
        color <- "#000000" 
      }
      return(color)
    })
    
    cCRE_legend <- Legend(
      labels = legend_labels, 
      legend_gp = gpar(fill = at_colors, fontsize = 6),
      grid_height = unit(0.5, 'cm'), 
      grid_width = unit(0.5, 'cm'), 
      title = "cCREs"
    )
    pushViewport(viewport(x = 0.9, y = 0.7))
    grid.draw(cCRE_legend)
    upViewport()
  }
  
  if(nrow(gene_loc_input) > 0){
    gene_legend <- Legend(
      labels = gene_loc_input$gene, 
      legend_gp = gpar(fill = gene_cols, fontsize = 6),
      grid_height = unit(0.5, 'cm'), 
      grid_width = unit(0.5, 'cm'), 
      title = "3D-interacting genes"
    )
    pushViewport(viewport(x = 0.9, y = 0.4))
    grid.draw(gene_legend)
    upViewport()
  }
  
  title(chrom)
  dev.off()
  circos.clear()
  
}

info_table <- query_snp("rs6960438")
cCRE_intersect <- get_snp_ccre("rs6960438")
head(info_table)
cCRE_intersect
setwd("G:/")
snp_ids <- c("rs759166", "rs1024749", "rs1024750", "rs6964933", "rs723526", "rs1344307", "rs1861007", "rs6960438", "rs759168")

plot_snp_circos(
  snp_info_table = info_table, 
  user_snp_list = snp_ids, 
  output_assembly = "hg19", 
  window_size = 2e5, 
  circos_pos = "Relative", 
  savefile = "circos5.pdf"
)