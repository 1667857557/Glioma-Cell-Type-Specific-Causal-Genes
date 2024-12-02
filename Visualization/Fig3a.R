suppressPackageStartupMessages({
  pacman::p_load(ggplot2, Seurat, dittoSeq, patchwork, plot1cell, preprocessCore, vroom)
})

setwd("D:/single_eqtl")
reference <- readRDS('C:/Users/huang/Downloads/643b2784-22c1-43e9-876f-e0fa0f5d75f4.rds')
meta<-reference@meta.data
unique(meta$annotation_level_3)
gc()
meta <- meta %>%
  mutate(annotation_level_3 = if_else(
    annotation_level_3 %in% c("MES-like", "OPC-like", "AC-like", "NPC-like"),
    annotation_level_4,  
    annotation_level_3   
  ))
reference<-AddMetaData(reference,meta)
reference$Cluster <- reference$annotation_level_3
unique(reference$Cluster)
circ_data <- prepare_circlize_data(reference, scale = 0.75)
gc()
circ_data$Cluster<-circ_data$annotation_level_3
set.seed(1234)
unique(circ_data$annotation_level_4)
circ_data$Cluster <- factor(x = circ_data$Cluster, 
                            levels = c(
                              'Mural cell','Endothelial',
                              rev(c('Mast', 'Mono','TAM-BDM','TAM-MG','DC')),
                              rev(c('NK','CD4/CD8','B cell','Plasma B')),
                              'AC-like','AC-like Prolif','MES-like hypoxia/MHC','MES-like hypoxia independent','OPC-like','OPC-like Prolif','NPC-like OPC','NPC-like Prolif','NPC-like neural',
                              'Astrocyte','Oligodendrocyte','RG', 'OPC','Neuron'
                              
                              
                            ))

cluster_colors<- c(
  '#756bb1','#17becf',
  rev(c('#8c564b','#393b79','#637939','#8c6d31','#636363')),
  rev(c('#7fc97f','#7b4173','#beaed4','#fdc086')),
  '#d62728',"#A2005699", '#ff7f0e',"#DF8F44FF", '#3182bd',"#374E55FF",'#31a354',"#008B4599","#79AF97FF",
  '#e7298a','#e6ab02', '#bcbd22','#bf5b17','#00F6B3'
)
pdf('circlize_plot.pdf', width = 10, height = 10)
png(filename =  'circlize_plot.png', width = 8, height = 8, units = 'in', res = 300)

plot_circlize(circ_data,do.label = T, pt.size = 0.5, col.use = cluster_colors,
              bg.color = 'white', kde2d.n = 1000, repel = T, label.cex = 0.7)

add_track(circ_data, group = "author", colors = c(RColorBrewer::brewer.pal(8, "Dark2"),
                                                  RColorBrewer::brewer.pal(8, "Accent")), 
          track_num = 2) ## can change it to one of the columns in the meta data of your seurat object
add_track(circ_data, group = "assay",colors = RColorBrewer::brewer.pal(9, "Set1"), 
          track_num = 3) ## can change it to one of the columns in the meta data of your seurat object
add_track(circ_data, group = "donor_id",colors = dittoColors(), 
          track_num = 4)
add_track(circ_data, group = "tissue", colors = c(RColorBrewer::brewer.pal(8, "Dark2"),
                                                  RColorBrewer::brewer.pal(8, "Accent")), 
          track_num = 5) 
dev.off()
gc()