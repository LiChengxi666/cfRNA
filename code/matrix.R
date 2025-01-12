# 1.过滤低表达基因
library(edgeR)

# 读取计数矩阵
count_data <- read.table("D:/cfRNA/count.matrix.long.txt", header = TRUE, row.names = 1)
matadta <- read.table("D:/cfRNA/metadata.long.txt", header = TRUE, row.names = 1)

# 创建DGEList对象
dge <- DGEList(counts = count_data)

# 过滤低表达基因
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]

# 2.标准化
dge <- calcNormFactors(dge, method="TMM")
cpm.matrix <- cpm(dge)

library(EDASeq)
plotRLE(cpm.matrix)

# 3.去除批次效应
suppressPackageStartupMessages(library(RUVSeq))
# calculate coefficient of variation
cv.index <- apply(cpm.matrix,1,function(x){sd(x)/mean(x)})
cv.index.sorted <- sort(cv.index)
# using 20% RNA with smallest cv as empirical control gene as there is little knowledge for stably expressed miRNA
empirical.control <- head(names(cv.index.sorted),as.integer(length(cv.index.sorted)*0.2))
# perform batch  correction with RUVg
res <- RUVg(cpm.matrix, empirical.control, k=2,isLog=TRUE)
log.tmm.matrix.ruvg <- res$normalizedCounts
log.tmm.matrix.ruvg <- as.data.frame(log.tmm.matrix.ruvg)
write.table(log.tmm.matrix.ruvg,"D:/cfRNA/count.matrix.long.tmm.txt",sep="\t",quote = F)
# saving results
# write.table(log.tmm.matrix.ruvg,"small/miRNA.TMM.logged.ruvg.txt",sep="\t",quote = F)