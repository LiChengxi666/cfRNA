# 加载必要的库
library(pheatmap)

# 读取计数矩阵
count_data <- read.table("D:/cfRNA/selected_data_short.txt",sep='\t', header = TRUE, row.names = 1)
metadata <- read.table("D:/cfRNA/metadata.short.txt", sep='\t', header = TRUE, row.names = 1)
count_matrix <- t(count_data)
selected_columns <- rownames(metadata)[metadata$label %in% c("HD", "CRC")]
count_matrix <- count_matrix[, selected_columns]
metadata <- metadata[selected_columns, ]


# 设置注释
annotation <- metadata$label
annotation_col = data.frame(Type = factor(annotation))
rownames(annotation_col) = colnames(count_matrix)
ann_colors = list(Type = c("CRC"="brown3", "HD"="darkgreen"))
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))

# 绘图并保存
pheatmap(count_matrix,
         clustering_method = "ward.D2",
         cutree_row = 2, cutree_col = 2, 
         cluster_rows=TRUE, show_rownames = TRUE, 
         cluster_cols=FALSE, show_colnames = FALSE,
         annotation_col = annotation_col, annotation_colors = ann_colors,
         scale = "none",
         color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),
                   colorRampPalette(colors = c("white","red"))(length(bk)/2)),
         breaks = bk,
)
