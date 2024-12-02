setwd("G:/A/target/enrichment")
pacman::p_load("vroom","DOSE","org.Hs.eg.db","topGO",
               "tidyr","dplyr","clusterProfiler","pathview",
               "ggplot2","forcats","aplot","ggtree","reshape2","tidyverse",
               "ggsci","scales","data.table")
genes <- fread("G:/A/target/enrichment/noGBM.txt",header=T,sep='\t')
name_ID = bitr(genes$gene, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
kegg <- enrichKEGG(
  gene = name_ID$ENTREZID,  
  keyType = 'kegg', 
  organism = 'hsa', 
  pAdjustMethod = 'bonferroni',
  pvalueCutoff = 0.05)  
barplot(kegg,showCategory = 20,title="The KEGG enrichment analysis of all identified genes ") 
goplot(go_bp)
write.table(kegg@result, 'noGBM_kegg.txt', sep = '\t', quote = FALSE, row.names = FALSE)


data <- vroom('all_glioma_kegg.txt')
colnames(data)
data$Log10P<-(-log10(data$pvalue))
data<- subset(data, p.adjust<0.05)
data<-select(data,c(Description,Count,Log10P))
colnames(data)<-c("KEGG_terms","Gene_count","-Log10(P value)")
p1 <- ggplot(data) +
  geom_col(aes(x = reorder(KEGG_terms, -`-Log10(P value)`),
               y = Gene_count),
           color = 'black',   
           width = 0.6,       
           fill = '#A59ACA') +
  labs(x = NULL, y = 'Gene count',
       title = 'Enriched KEGG Terms in All-Glioma (P.adjust<0.05)') +
  theme_test(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(color = 'black', face = 'bold', size = 13),
    plot.margin = margin(1, 0.5, 0.5, 2.5, 'cm'),
    panel.border = element_rect(size = 1),
    axis.title = element_text(face = 'bold', size = 18), 
    axis.title.y = element_text(face = 'bold', size = 14), 
    plot.title = element_text(face = 'bold', size = 14, hjust = 0.5)
  )
p1
p2 <- p1+
  scale_y_continuous(expand = c(0,0),limits = c(0,12),
                     sec.axis = sec_axis(~./1,
                                         name = '-Log10(P value)',
                                         breaks = seq(0,10,5)))+
  geom_line(aes(x= reorder(KEGG_terms,-`-Log10(P value)`),
                y=`-Log10(P value)`*1,
                group=1),
            linetype=3,cex=0.6)+
  geom_point(aes(x= reorder(KEGG_terms,-`-Log10(P value)`),
                 y=`-Log10(P value)`*1),
             color = "black", fill = '#88D277', shape = 21, size=3.5)+
  theme(axis.title.y.right = element_text(size = 14))
p2
p3 <- p2+
  geom_text(aes(x= reorder(KEGG_terms,-`-Log10(P value)`),
                y=Gene_count,
                label=Gene_count),
            vjust=-0.5,size=3.5,fontface='bold')
p3
df <- data.frame(a=c(1.5,1.5,11.5,11.5), 
                 b=c(12,11,11,12))
p4 <- p3+
  geom_line(data = df,aes(a,b),cex=0.5)+
  geom_rect(aes(xmin=2,xmax=3,ymin=11.35,ymax=11.65),
            fill='#A59ACA',color='black')+ 
  annotate(geom='text',x=4.5,y=11.5,label='Gene Count',
           fontface='bold',size=4)+
  annotate('segment',x=7.1,xend = 8.1,y=11.5,yend = 11.5,
           linetype=3,cex=0.5)+
  annotate(geom='point', x=7.6,y=11.5,
           color = "black", fill = '#88D277', shape = 21, size = 3)+
  annotate('text',x=9.8,y=11.5,label='-Log10(P value)',
           fontface='bold',size=4)
p4
ggsave("All-Glioma_KEGG_doubleY-enrih-bar.png", width = 12, height = 7, dpi = 600)
ggsave("All-Glioma_KEGG_doubleY-enrih-bar.pdf", width = 12, height = 7)
