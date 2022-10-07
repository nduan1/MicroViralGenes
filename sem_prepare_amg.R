colnames(microbe_gene_table)
microbe_gene_cate <- microbe_gene_table[,c(2:14)] %>% 
  group_by(genes) %>% 
  summarise(across(everything(),sum))
dim(microbe_gene_cate)
colnames(microbe_gene_cate)
colnames(amg_gene_tab.stack)[1] <- "genes"
colnames(amg_gene_tab.stack)
amg.gene <- rbind(microbe_gene_cate,amg_gene_tab.stack)
dim(amg.gene)
library(DESeq2)
rownames(amg.gene)
amg.gene.new <- amg.gene[,-1]
colnames(amg.gene.new)
rownames(amg.gene.new) <-(amg.gene)$genes
rownames(amg.gene.new)
dds.new3 <- DESeqDataSetFromMatrix(countData = amg.gene.new,
                                   colData = meta_data.1,
                                   design = ~ Nitrogen*Cover)
nrow(dds.new3)
dds.new3 <- DESeq(dds.new3)
nrow(dds.new3)
sizeFactors(dds.new3)
plotDispEsts(dds.new3)
resultsNames(dds.new3)
amg.gene.nor <- as.data.frame(counts(dds.new3,normalized=T))
head(amg.gene.nor)
view(amg.gene.nor)
amg.gene.nor.t <- as.data.frame(t(amg.gene.nor))
view(amg.gene.nor.t)
write.csv(amg.gene.nor.t,"amg.gene.nor.t.csv")


