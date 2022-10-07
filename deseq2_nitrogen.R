#ref:http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#countmat
#DESeq want their normalization to handle:1.differences in library sizes
#2. differences in library composition
#log won't influence 
head(nitrogen_gene_tab.new)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DESeq2")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("pasilla")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("apeglm")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("IHW")
library(DESeq2)
library(pasilla)
library(apeglm)
library(IHW)
####raw dataset input####
# The DESeq2 model internally corrects for library size, so transformed or normalized values such as counts scaled by library size should not be used as input.
colnames(nitrogen_gene_tab.new)
#unique the gene with the same gene annotation
nitrogen.deseq <- read.csv("dataset/nitrogen_deseq.csv")
colnames(nitrogen.deseq);dim(nitrogen.deseq)
nitrogen.deseq.up <- nitrogen.deseq[,c(6,9:20)]
dim(nitrogen.deseq.up)
#nitrogen.deseq.up.1 <- nitrogen.deseq.up[rowSums(nitrogen.deseq.up==0)<9,]#remove the 0 happened more than 3 times
nitrogen.deseq.up.1 <- nitrogen.deseq.up
dim(nitrogen.deseq.up.1)#300 clean dataset
nitrogen.deseq.tab <- nitrogen.deseq.up.1 %>% group_by(comp_gene) %>% summarise(across(everything(), sum))
dim(nitrogen.deseq.tab)#19 13
head(nitrogen.deseq.tab)
#keep Nitrification=hao+amoA+amoB+amoC,
# Denitrification=norB+nirK+nirS+nosZ,
# Assimilatory=nirA+nasA,
# Dissimilatory=nirB+nrfA+nrfH,
# Transport=`NRT,narK,nrtP,nasA`
# nitrogen.deseq.tab <- subset(nitrogen.deseq.tab,comp_gene!="napA"&comp_gene!="napB"&
#                                comp_gene!="narG,narZ,nxrA"&comp_gene!="narH,narY,nxrB"&
#                                comp_gene!="narI,narV")

nitrocts.tab <- nitrogen.deseq.tab
colnames(nitrocts.tab)[1] <- "comp_gene"
# rownames(nitrocts) <- nitrocts$comp_gene#counts dataset
# head(nitrocts) 
nitroanno.tab <- unique(nitrogen.deseq[,5:7]) %>% inner_join(nitrocts.tab,by="comp_gene")
head(nitroanno.tab)
nitroanno <- nitroanno.tab[,1:3]#annotation dataset
rownames(as.data.frame(nitrocts.tab)) <- as.data.frame(nitrocts.tab)[,1]
nitrocts <- (nitrocts.tab[,-1])#counts dataset
rownames(nitrocts) <- nitrocts.tab$comp_gene
metadata.1 <- read.csv("dataset/meta_data.csv",row.names = 1)#metadata
head(metadata.1)
metadata <- cbind(as.data.frame(metadata.1[,1:4]),as.data.frame(scale(metadata.1[,5:23])))
all(rownames(metadata) %in% colnames(nitrocts))#true 
all(rownames(metadata)==colnames(nitrocts))#true
#####take a look at the data###
head(nitrocts)
barplot(colSums(nitrocts),las=3)#check the total counts for each sample
hist(nitrocts$NCNTN60_2,br=100)#check the skew
nitro_logcount <- log2(1+nitrocts)
hist(nitro_logcount$NCNTN60_2,br=100)#check the skew
plot(nitro_logcount[,1],nitro_logcount[,2])#if they are identical among the replicates, 
#it should be the 45 degree angle if replicate are ideally identical
plot(nitro_logcount[,1],nitro_logcount[,6])#different treatment
coldata <- metadata
head(coldata)
####create dds object####
dds <- DESeqDataSetFromMatrix(countData = nitrocts,
                              colData = metadata,
                              design = ~ Nitrogen + Cover+Nitrogen:Cover)
nrow(dds)#23
####pre-filter#
keep <- rowSums(counts(dds)) >= 5 #elimited low counting genes
dds <- dds[keep,]
nrow(dds)#22
###Differential expression analysis####
colData(dds)
#levels(ddsMF$Cover)#
dds <- DESeq(dds)
nrow(dds)
sizeFactors(dds)#check the reads number in each treatment, it will go with your size factor

####PCA###
library(ggplot2)
rld <- rlog(dds)
plotPCA(rld,intgroup="Nitrogen")
plotPCA(rld,intgroup="Cover")
plotPCA(rld,intgroup=c("Nitrogen","Cover"))
####heatmap#####
#counts heatmap
library(pheatmap)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                   decreasing=TRUE)
df <- as.data.frame(colData(dds)[,c("Cover","Nitrogen")])
pheatmap(assay(rld), cluster_rows=T, show_rownames=T,
         cluster_cols=T, annotation_col=df)

#overall check if there any gene expression difference
head(assay(rld))
sampleDists <- dist(t(assay(rld)))#compute the distances between the rows of a data matrix., so t()
library("RColorBrewer")
library(pheatmap)
sampleDistMatrix <- as.matrix(sampleDists)
#rownames(sampleDistMatrix) <- paste(rld$Cover, rld$Nitrogen, sep="-")
#colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)# the results should be the same treatment clustered together
#I found NCNTNO_1 is different from 2 and 4, and it seems like fertilization have effect on the cluster
#So I am planning to remove the NCNTN0_1 and then do the downstream analysis
##remove outlier samples#
# nitrocts.up <- nitrocts[,-1]#remove NCNTN0_1 
# nitroanno.up <- nitroanno
# metadata.up <- metadata[-1,]
# head(metadata.up)
# dds.up <- DESeqDataSetFromMatrix(countData = nitrocts.up,
#                                  colData = metadata.up,
#                                  design = ~ Nitrogen + Cover+Nitrogen*Cover)
# keep <- rowSums(counts(dds.up)) >= 5 #elimited low counting genes
# dds.up <- dds.up[keep,]
# nrow(dds)
# dds.up <- DESeq(dds.up)
# nrow(dds.up)
# sizeFactors(dds.up)
# ##PCA_clean#
# library(ggplot2)
# rld.up <- rlog(dds.up)
# plotPCA(rld.up,intgroup="Nitrogen")
# plotPCA(rld.up,intgroup="Cover")
# plotPCA(rld.up,intgroup=c("Nitrogen","Cover"))
# #heatmap_clean#
# head(assay(rld.up))
# sampleDists.up <- dist(t(assay(rld.up)))#compute the distances between the rows of a data matrix., so t()
# library("RColorBrewer")
# library(pheatmap)
# sampleDistMatrix.up <- as.matrix(sampleDists.up)
# colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
# #distance heatmap
# pheatmap(sampleDistMatrix.up,
#          clustering_distance_rows=sampleDists.up,
#          clustering_distance_cols=sampleDists.up,
#          col=colors)
# #counts heatmap
# select.up <- order(rowMeans(counts(dds.up,normalized=TRUE)),
#                 decreasing=TRUE)
# df.up <- as.data.frame(colData(dds.up)[,c("Cover","Nitrogen")])
# pheatmap(assay(rld.up), cluster_rows=FALSE, show_rownames=T,
#          cluster_cols=FALSE, annotation_col=df.up)

#PCA for paper
library(phyloseq)
library(vegan)
#install.packages("ggforce")
library(ggforce)
rownames(nitroanno.up) <- nitroanno.up$comp_gene
gene.phyloseq <- merge_phyloseq(sample_data(metadata.up), otu_table(as.data.frame(assay(rld.up)),taxa_are_rows = TRUE))
gene.phyloseq
df_pca <- prcomp(t(otu_table(gene.phyloseq)))
pca_tab <- as.data.frame(df_pca$x)
pca_tab$Nitrogen <- metadata.up$Nitrogen
pca_tab$Cover <- metadata.up$Cover

#RDA

# formularda <- rda(t(otu_table(gene.phyloseq)) ~ Cover+Nitrogen+Cover*Nitrogen+pH+GWC+NH4+NO3+
#                     AG+BG+CB+LAP+NAG+PHOS+XYL,data = metadata.up)
# summary(formularda)
rda_data <- ordinate(gene.phyloseq,method = "RDA",formula = ~Cover+Nitrogen)
env.rda <- rda_data$CCA
arrow.dat <- as.data.frame(env.rda$biplot)
set.seed(2021)
step.backward <-
  ordistep(rda_data,
           permutations = how(nperm = 999)
  )
Anova(rda_data)
anova.perm <- anova(rda_data,permutations = how(nperm = 9999), by="margin");anova.perm


adj <- RsquareAdj(rda_data);adj
adonis(t(otu_table(gene.phyloseq))~Nitrogen+Cover+pH+GWC+NH4+NO3,
       data=data.frame(sample_data(gene.phyloseq)),
     permutations=9999,by="margin")
# #only nitrogen significant
# tiff('figures/PCA_nitrogen.tiff', units="in", width=7, height=5, res=300)
# rda.plot <- plot_ordination(gene.phyloseq,rda_data,
#                 type="sample",color="Nitrogen")+
#   geom_point(size=3)+
#   geom_text(aes(label=rownames(pca_tab),vjust=-1),size=2.5)+
#   ggtitle("N-cycling genes")+
#   scale_color_manual(values = c("#106494", "#bf6a0a"))+
# #  xlab(label = "PC1 [53.1%]")+
# #  ylab(label = "PC2 [20%]")+
#   geom_mark_ellipse(aes(color = Nitrogen, group=Nitrogen),expand=0)+
#   geom_vline(xintercept = c(0), size=0.2,linetype=2)+
#   geom_hline(yintercept = c(0), size=0.2,linetype=2)+
#   theme(legend.position="bottom",panel.background = element_rect(fill = "white",color="black"),
#         panel.grid.major = element_line(colour="grey",size = 0.1),panel.grid.minor = element_blank(),
#         plot.title=(element_text(size=13,family = "Arial",face="plain",vjust = 3,hjust=0.5)),
#         text=element_text(family = "Arial",face="plain",size = 15),
#         panel.border = element_rect(colour = "black", fill=NA, size=0.5),
#         axis.text.y=element_text(size=13,family="Arial",face = "plain"),
#         axis.title=element_text(size =13,face="plain",family = "Arial",vjust = 1),
#         legend.title = element_text(size=13,face="plain",family = "Arial"),
#         legend.text = (element_text(size=13,family = "Arial")))
# rda_env_plot <-rda.plot+ 
#   geom_segment(data=arrow.dat, aes(x=0, xend=RDA1, y=0, yend=RDA2), 
#                          color="black", arrow=arrow(length=unit(0.01,"npc"))) +
#   geom_text(data=arrow.dat, 
#             aes(x=RDA1,y=RDA2,label=rownames(arrow.dat),
#                 hjust=0.5*(1-sign(RDA1)),vjust=0.5*(1-sign(RDA1))), 
#             color="red", size=4)
tiff('figures/pca.tiff', units="in", width=10, height=10, res=300)
ggplot(pca_tab,aes(x=PC1,y=PC2,type="samples",color=Nitrogen, shape = Cover))+
  geom_point(size=3)+
  scale_shape_manual(values = c(17,8))+
  geom_text(aes(label=rownames(pca_tab),vjust=-1),size=2.5)+
  ggtitle("N-cycling genes")+
  scale_color_manual(values = c("#106494", "#bf6a0a"))+
  xlab(label = "PC1 [53.1%]")+
  ylab(label = "PC2 [20%]")+
  geom_mark_ellipse(aes(color = Nitrogen, group=Nitrogen),expand=0)+
  geom_vline(xintercept = c(0), size=0.2,linetype=2)+
  geom_hline(yintercept = c(0), size=0.2,linetype=2)+
  theme(legend.position="bottom",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_line(colour="grey",size = 0.1),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=15,family = "Arial",face="plain",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="plain",size = 15),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(size=15,family="Arial",face = "plain"),
        axis.title=element_text(size =15,face="plain",family = "Arial",vjust = 1),
        legend.title = element_text(size=15,face="plain",family = "Arial"),
        legend.text = (element_text(size=15,family = "Arial")))
  # annotate(geom="text", x=1.0, y=1.5, size=5,label="Adj.R^{2}",parse=T,color="black",fontface="italic")+
  # annotate(geom="text", x=1.3, y=1.5, size=5,label=paste(" =",round(adj$adj.r.squared,2)),color="black", fontface="plain")

dev.off()


#######get results of deseq####
plotDispEsts(dds)
resultsNames(dds)
res_nitrogen_nocover<- results(dds,contrast=c("Nitrogen","N60","N0"),lfcThreshold = 0.01)#2
res_cover<- results(dds,list(c("Cover_Vetch_vs_No_cover","NitrogenN60.CoverVetch")),lfcThreshold = 0.01)#not sig
res_nitrogen<- results(dds,list(c("Nitrogen_N60_vs_N0","NitrogenN60.CoverVetch")),lfcThreshold = 0.01)#3 and 2
res_all<- results(dds,lfcThreshold = 0.01)#1 and 1, no
res_cover_n0<- results(dds,contrast=c("Cover","Vetch","No_cover"),lfcThreshold = 0.01)#no

head(res)
summary(res)# the interaction term  

#####ma plot####
plotMA(res, ylim=c(-5,5))
#blue points are the gene significant and grey one is not significant based on deseq calculation
#above the horziontal line, the gene is upregulated
#x is the number of normalized reads of particluar gene has
#nomalized reads and log fold change relationship, 
#if we have logfold change is really high,then we can get a more lower normalized read counts 
####identify significant##
library(dplyr)
res1 <- as.data.frame(res)
res1 <- mutate(res1,sig=ifelse(res1$padj<0.1,"FDR<0.1","Not Sig"))
res1[which(abs(res1$log2FoldChange)<0.1),"sig"]="Not Sig" #greater than 1.0 keep it as significant, otherwise not significant
res1##




















##################################################################################################################################
#"sig" here is significant column
library(ggplot2)
ggplot(res1, aes(log2FoldChange, -log10(padj)))+geom_point(aes(col=sig))+
  scale_color_manual(values = c("red","black"))#sig red, not sig is black
#broadly check the change

####check individual genes####
resMF.nitrogen=resMF.nitrogen[order(abs(resMF.nitrogen$log2FoldChange), decreasing = T),] #order the rows from the largest to the lowest
topgene <- rownames(resMF.nitrogen)[1]
topgene
plotCounts(ddsMF,gene = topgene,intgroup = c("Nitrogen","Cover"))
####fit the gene names###
library(AnnotationDbi)#convert refseqid to gene symbol
library(org.Mm.eg.db)#how to convert between geneid and the symbol
#if we want to replace the refseq geneid to gene symbol, we need to remove the refseq id ".1"
resMF.nitrogen$refseq=gsub("\\..*","",row.names(resMF.nitrogen))

#####gene ontology analysis####
library(Go.db)
resMF.nitrogen <- results(ddsMF, contrast = c("Nitrogen","N0","N60"))
#null is that the logarithmic fold change (LFC) between treatment
#and control for a gene’s expression is exactly zero,
#even if statistically highly significant,
#might not be the most interesting candidates for further investigation
head(resMF.nitrogen)
resMF.nitrogen.ordered <- resMF.nitrogen[order(resMF.nitrogen$pvalue),]
summary(resMF.nitrogen.ordered)
sum(resMF.nitrogen.ordered$padj < 0.1, na.rm=TRUE)#58 genes padj<0.1
resMF.nitrogen05 <- results(ddsMF, contrast = c("Nitrogen","N0","N60"),alpha = 0.05)
summary(resMF.nitrogen05)
sum(resMF.nitrogen05$padj < 0.05, na.rm=TRUE)#44

resMF.cover <- results(ddsMF, contrast = c("Cover","No_cover","Vetch"))
head(resMF.cover)
resMF.cover.ordered <- resMF.cover[order(resMF.cover$pvalue),]
summary(resMF.cover.ordered)
sum(resMF.cover.ordered$padj < 0.1, na.rm=TRUE)#5 genes padj<0.1
resMF.cover05 <- results(ddsMF, contrast = c("Cover","No_cover","Vetch"),alpha = 0.05)
summary(resMF.cover05)
names(resMF.cover05)
sum(resMF.cover05$padj < 0.05, na.rm=TRUE)#5


plotDispEsts(ddsMF)
plotMA(resMF.cover)

featureData <- data.frame((nitroanno))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)
dds
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dim(dds)
dds$Treat <- relevel(dds$Treat, ref = "NCNTN0")
dds.deseq <- DESeq(dds)
dds.res <- results(dds.deseq,contrast = c("Treat","NCNTN0","NCNTN60"))#use contrast argument to extract the treate you want to compare
dds.res  
resultsNames(dds.deseq) 
resLFC <- lfcShrink(dds.deseq,coef = "Treat_NCNTN60_vs_NCNTN0",type="apeglm")
resLFC  
resOrdered <- dds.res[order(dds.res$pvalue),]  
summary(dds.res)
sum(dds.res$padj < 0.1, na.rm=TRUE)
res05 <- results(dds.deseq, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)

#######Independent hypothesis weighting####
library(IHW)
resIHW <- results(dds.deseq, filterFun=ihw)
summary(resIHW)
sum(resIHW$padj < 0.1, na.rm=TRUE)
metadata(resIHW)$ihwResult
plotMA(dds.res, ylim=c(-2,2))
plotMA(resLFC, ylim=c(-2,2))
idx <- identify(dds.res$baseMean, dds.res$log2FoldChange)
rownames(dds.res)[idx]
