MDA.cellchat<-readRDS("G:/A/cell_specify/pseudotime/CellChatDB_cellphonedb_cellchat_GBM_maligant_cell_type.rds.gz")
ks_cellchatV2_sig_interCount <- function(cellchat, 
                                         count = F, 
                                         source_celltype, 
                                         showInter = FALSE,
                                         color_set = NULL, 
                                         celltype_order = NULL, 
                                         celltype.size = FALSE, 
                                         flipped = TRUE 
){
  
  require("reshape2")
  require("ggplot2")
  require("ggraph")
  require("tidygraph")
  require("dplyr")
  require("igraph")
  require("CellChat")
  require("tools") 
  
  if(count){
    mat <- as.data.frame(cellchat@net$count)
  } else {
    mat <- as.data.frame(cellchat@net$weight)
  }
  
  sourecell <- which(colnames(mat) == source_celltype)
  
  mat <- mat[order(mat[, sourecell], decreasing = TRUE), ] 
  
  group.size <- as.data.frame(table(as.character(cellchat@idents))) 
  
  df <- data.frame(from = rep(source_celltype, nrow(mat)),
                   to = rownames(mat),
                   inter = mat[, source_celltype],
                   stringsAsFactors = FALSE)
  
  nodes <- data.frame(name = as.character(df$to), stringsAsFactors = FALSE)
  nodes$inter <- df$inter
  
  nodes$name <- tools::toTitleCase(nodes$name)
  
  if(!celltype.size){
    size <- rep(5, nrow(mat))
    nodes$size <- size
    nodes <- nodes[order(nodes[, "inter"], decreasing = TRUE), ]
  } else {
    colnames(group.size) <- c('name', "size")
    group.size$name <- tools::toTitleCase(as.character(group.size$name))
    nodes <- merge(nodes, group.size, by = 'name', all = FALSE)
    nodes <- nodes[order(nodes[, "size"], decreasing = TRUE), ]
  }
  
  edges <- df[c("from", "to", "inter")]
  edges$from <- tools::toTitleCase(as.character(edges$from))
  edges$to <- tools::toTitleCase(as.character(edges$to))
  
  if(is.null(color_set)){
    color.use <- c("#D51F26","#272E6A","#208A42","#89288F","#F47D2B", 
                   "#FEE500","#8A9FD1","#C06CAB", "#D8A767",
                   "#90D5E4", "#89C75F","#F37B7D","#9983BD",
                   "#D24B27","#3BBCA8", "#6E4B9E","#0C727C", "#7E1416")
    
    color.use <- color.use[1:nrow(mat)]
    names(color.use) <- nodes$name
  } else {
    color.use <- color_set
    if(is.null(celltype_order)){
      names(color.use) <- nodes$name
    } else {
      names(color.use) <- celltype_order
    }
  }
  
  net <- tbl_graph(nodes = nodes, edges = edges)
  
  p <- ggraph(net, layout = 'igraph', algorithm = 'circle') +
    geom_edge_bend(aes(edge_width = inter),
                   strength = 0.2, alpha = 0.8,
                   flipped = flipped, edge_color = "#A9AAAA",
                   n = 50, show.legend = FALSE,
                   check_overlap = TRUE) +
    geom_edge_loop(aes(edge_width = inter,
                       direction = (from - 1) * 360 / length(net)),
                   colour = "#A9AAAA",
                   alpha = 0.5, show.legend = FALSE) +
    scale_edge_width_continuous(range = c(0, 5))
  
  if(!showInter){
    p <- p + 
      geom_node_point(aes(size = size, colour = name), show.legend = FALSE) +
      geom_node_point(aes(size = size), show.legend = FALSE,
                      shape = 21, colour = 'black', stroke = 1.5) +
      geom_node_text(aes(x = x * 1.06, y = y * 1.06, label = name),
                     angle = 0, hjust = 0, size = 3, family = "sans") + 
      scale_size_continuous(range = c(1, 15)) +
      scale_color_manual(values = color.use) +
      theme_graph() +
      theme(plot.margin = unit(c(2, 2, 2, 2), "cm"))
  } else {
    if(celltype.size){
      p <- p + 
        geom_node_point(aes(size = size, colour = inter)) +
        geom_node_point(aes(size = size), show.legend = FALSE,
                        shape = 21, colour = 'black', stroke = 1.5) +
        geom_node_text(aes(x = x * 1.06, y = y * 1.06, label = name),
                       angle = 0, hjust = 0, fontface = "bold", size = 3, family = "sans") + 
        scale_color_gradientn(colors = colorRampPalette(c("#2166AC",'#478ABF','#90C0DC', 
                                                          "white",'#EF8C65','#CF4F45',"#B2182B"))(100)) +
        scale_size_continuous(range = c(1, 10)) +
        theme_graph() +
        theme(plot.margin = unit(c(2, 2, 2, 2), "cm"))
    } else {
      p <- p + 
        geom_node_point(aes(size = inter, colour = inter)) +
        geom_node_point(aes(size = inter), show.legend = FALSE,
                        shape = 21, colour = 'black', stroke = 1.5) +
        geom_node_text(aes(x = x * 1.06, y = y * 1.06, label = name),
                       angle = 0, hjust = 0, fontface = "bold", size = 3, family = "sans") + 
        scale_color_gradientn(colors = colorRampPalette(c("#2166AC",'#478ABF','#90C0DC', 
                                                          "white",'#EF8C65','#CF4F45',"#B2182B"))(100)) +
        scale_size_continuous(range = c(1, 10)) +
        theme_graph() +
        theme(plot.margin = unit(c(2, 2, 2, 2), "cm"))
    }
  }
  
  return(p)
}

setwd("G:/")

pdf("cellchat_malignant_weight.pdf", width = 11, height = 8)
tiff('cellchat_malignant_weight.tiff',units = "in",width = 11, height = 8,res= 600,compression = "lzw")

p <- ks_cellchatV2_sig_interCount(cellchat = MDA.cellchat,
                                  source_celltype = "malignant cell",
                                  count = F,
                                  celltype.size = TRUE,
                                  showInter = TRUE,
                                  flipped = TRUE)
print(p)
dev.off()
gc()
