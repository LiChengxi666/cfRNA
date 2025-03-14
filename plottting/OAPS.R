library(dplyr)
library(ComplexHeatmap) 
library(ggplot2)
library(tidyr)
library(edgeR)
library(circlize)
library(org.Hs.eg.db)
library(clusterProfiler)

file_path <- "D:/OAPS/"
oaps <- read.table(paste0(file_path,"gencode.txt"),header = T,row.names = 1, 
                   check.names = F)
oaps <- oaps[,!grepl("1123",colnames(oaps))]

feature <- data.frame(geneid=rownames(oaps))
rownames(feature) <- feature$geneid
feature <- separate(data = feature,col=geneid, into=c("geneid","length","genename","type"),sep="\\|")

feature[grepl("protein",feature$type),]$type="mRNA"
feature[grepl("vault_RNA",feature$type),]$type="sncRNA"
feature[grepl("scaRNA",feature$type),]$type="snoRNA"
feature[grepl("TEC",feature$type),]$type="misc_RNA"
feature[grepl("scRNA",feature$type),]$type="sncRNA"
feature[grepl("snoRNA",feature$type),]$type="sncRNA"
feature[grepl("miRNA",feature$type),]$type="sncRNA"
feature[grepl("snRNA",feature$type),]$type="sncRNA"
feature[grepl("sRNA",feature$type),]$type="sncRNA"
feature[grepl("ribozyme",feature$type),]$type="sncRNA"
feature[grepl("Mt",feature$type),]$type="sncRNA"


rownames(oaps) <- feature$geneid
rownames(feature) <- feature$geneid

feature <- feature[!grepl("pseudogene",feature$type),]
feature <- feature[!grepl("misc_RNA",feature$type),]

anno <- data.frame(sample=colnames(oaps))
anno$type <- c(1,0,1,0,1,0,0,0,0,0)
anno$type <- factor(ifelse(anno$type==1,"OAPS","Non-OAPS"))
rownames(anno) <- anno$sample
design <- model.matrix(~anno$type)
y <- DGEList(counts = oaps[feature$geneid,],group = anno$type)
y <- calcNormFactors(y, method="TMM")
y <- estimateDisp(y, design = design)
fit <- glmFit(y,coef=2)
lrt <- glmLRT(fit,coef=2)
diff.table <- topTags(lrt, n = nrow(y))$table
diff.table.filtered <- diff.table[abs(diff.table$logFC) > 0 & diff.table$PValue<0.01,]

CPM <- cpm(oaps,log = F)
CPM <- CPM[rowSums(CPM)>0,]
CPM <- log(CPM+1)
CPM <- (CPM-rowMeans(CPM))/apply(CPM,1,sd)

selected <- intersect(rownames(CPM),rownames(diff.table.filtered))
scount <- CPM[selected,]
rownames(scount) <- feature[rownames(scount),]$genename
dft <- diff.table.filtered[selected,]
dft$no <- c(1:nrow(dft))

maincol <- colorRamp2(c(min(scount),max(scount)),c("#2B0E4C","#E9E45C"))
annocol <- data.frame(row.names = colnames(CPM),type=anno[colnames(CPM),]$type)
annosd <- data.frame(row.names = feature[selected,]$genename,FDR=-log(dft$FDR), 
                     type=feature[selected,]$type) 
topanno <- HeatmapAnnotation(df=annocol,show_annotation_name = F,show_legend = T,
                             col=list(type=c('Non-OAPS'="#00D65C","OAPS"="#FF9189")))
rightanno <- HeatmapAnnotation(df=annosd,which='row',annotation_name_align = T,
                               show_annotation_name = F)
labelfeature <- rownames(dft)[1:20]
leftanno <- rowAnnotation(selected=anno_mark(at=c(dft[labelfeature,]$no),side = "left", 
                                             labels = feature[labelfeature,]$genename))
dcgnew <- Heatmap(scount,name="z-score", 
                  show_row_dend = F,show_row_names = F,show_column_dend = F,
                  show_column_names = T,column_split = annocol$type, top_annotation = topanno,
                  right_annotation = rightanno,col = maincol,cluster_rows = T, left_annotation = leftanno,
                  cluster_columns =T,
                  heatmap_width = unit(10, "cm"), heatmap_height = unit(18, "cm")) 
print(dcgnew)



mfeature <- feature[selected,]
mfeature <- mfeature %>% filter(mfeature$type=="mRNA" | mfeature$type=="IG_C_gene"| mfeature$type == "IG_V_gene" | mfeature$type=="TR_J_gene"| mfeature$type == "TR_V_gene")
enrichgo <- enrichGO(mfeature$genename, 
                     OrgDb = org.Hs.eg.db,keyType = "SYMBOL",ont="ALL",
                     pAdjustMethod = 'fdr',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
enrichtb <- enrichgo@result
go <-data.frame(term = enrichtb$Description,
                ONTOLOGY=enrichtb$ONTOLOGY,p.adjust=enrichtb$p.adjust)

go <- go[order(go$ONTOLOGY, go$p.adjust), ]
go$term <- factor(go$term, levels = go$term)

ggplot(go, aes(term, -log10(p.adjust))) +
  geom_col(aes(fill = p.adjust), width = 0.5, show.legend = FALSE) +
  #scale_fill_manual(values = c('#D06660', '#5AAD36', '#6C85F5')) +
  facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y') +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
  coord_flip() +
  labs(x = '', y = '-Log10 P-Value\n')+
  geom_hline(yintercept = -log10(0.05),colour='black',linetype=2,show.legend = T)


entrez_ids <- mapIds(
  org.Hs.eg.db,
  keys = mfeature$genename,
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)

enrich <- enrichKEGG(gene=entrez_ids, keyType = "kegg", organism = "hsa",
                     pAdjustMethod = 'fdr',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
enrichtb <- enrich@result
suppressMessages(enrichtb$geneID <- unlist(lapply(enrichtb$geneID,KEGGresidtransform)))
KEGGls[[paste0(name,"_",type[i])]] <- enrichtb

ggplot(enrichtb,aes(y=Description,x=Count,fill=qvalue))+
  
  geom_bar(stat = "identity",width=0.8)+ #柱状图宽度设置
  
  scale_fill_gradient(low = "red",high ="blue" )+
  
  labs(title = "KEGG Pathways Enrichment",  #设置标题、x轴和Y轴名称
       
       x = "Gene number",
       
       y = "Pathway")+
  
  theme(axis.title.x = element_text(face = "bold",size = 16),
        
        axis.title.y = element_text(face = "bold",size = 16),
        
        legend.title = element_text(face = "bold",size = 16))+
  
  theme_bw()



barplot(enrich,showCategory = 40,title = 'KEGG Pathway')