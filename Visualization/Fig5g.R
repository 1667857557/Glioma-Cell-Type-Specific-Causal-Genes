pacman::p_load(monocle3, dplyr, ggbreak, ggplot2)
cds_data<-readRDS("GBmap_core_astrocte_EGFR_monocle3_without_partition.rds.gz")
pseudotime <- pseudotime(cds_data) %>% as.data.frame()
pseudotime$cell <- rownames(pseudotime)
colnames(pseudotime)[1] <- "peu"
A <- readRDS("C:/Users/huang/Downloads/643b2784-22c1-43e9-876f-e0fa0f5d75f4.rds")
A<-A@meta.data
A$cell<-rownames(A)
unique(A$annotation_level_4)
celltype<-select(A,c(cell,annotation_level_4))
gc()

merge <-merge(pseudotime, celltype, by = 'cell')
merge <- merge[order(merge$peu), ]
unique(merge$annotation_level_4)
colnames(merge)[colnames(merge) == "annotation_level_4"] <- "Cell_type"
p=ggplot(merge,aes(peu,fill=Cell_type, color=Cell_type)) +  
  geom_density(alpha = 0.5,size=1.2)+
  scale_fill_manual(values=c("#374E55FF","#00A1D5FF","#B24745FF","#79AF97FF","#6A6599FF","#80796BFF","#DF8F44FF",
                             "#3B499299","#EE000099","#008B4599"))+  #设置颜色 
  scale_color_manual(values=c("#374E55FF","#00A1D5FF","#B24745FF","#79AF97FF","#6A6599FF","#80796BFF","#DF8F44FF",
                              "#3B499299","#EE000099","#008B4599"))+  
  theme_classic()+ 
  scale_y_continuous(expand = c(0,0))+
  theme(axis.line = element_blank(),
        axis.ticks.y=element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x=element_text(colour='black',size=10))+
  labs(title = "Pseudotime")
p
pdf("moncle_cell_differentation.pdf", width = 10, height = 4)
p+scale_y_break(c(1, 2),
                scales = c(0.5,2),
                expand=expansion(add = c(0, 0)),
                space = 0.2)
dev.off()

head(merge)
gc()

