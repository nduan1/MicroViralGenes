library(tidyverse)
library(DESeq2)
nitrogen.deseq <- read.csv("dataset/nitrogen_deseq.csv")
colnames(nitrogen.deseq);dim(nitrogen.deseq)

head(nitrogen.deseq)
dim(nitrogen.deseq)
colnames(nitrogen.deseq)
dim(nitrogen.deseq)
nitrogen.deseq.1 <- nitrogen.deseq[rowSums(nitrogen.deseq[,9:20]==0)<8,] #%>% 
  #filter(NCNTN60_4!=0)#60
nitro_all_gene.tab <- nitrogen.deseq.1[,-c(1,3,4,5,7,8)] %>% mutate(gene_contig=gsub("nonvirus_contigs_k127_","",contig))
colnames(nitro_all_gene.tab)
nitro_all_gene.dat <- nitro_all_gene.tab %>% 
  mutate(newname=paste(gene_contig,comp_gene,sep = "_"))
colnames(nitro_all_gene.dat)
nitro_new.dat <- nitro_all_gene.dat[,c(16,3:14)] 
##_----change
nitro_all_gene.dat$newname <- gsub(" ","",nitro_all_gene.dat$newname)
nitro_all_gene.dat$newname <- gsub(",","_",nitro_all_gene.dat$newname)

rownames(nitro_new.dat) <- nitro_all_gene.dat$newname
# nitro_new.2 <- nitro_new.dat %>% 
#   left_join(unique(nitro_all_gene.dat[,c(2,16)])) 
# nitro_new.1 <- nitro_new.2[,-1] %>% 
#   group_by(comp_gene) %>%
#   summarise(across(everything(), sum))
# rownames(nitro_new.1) <- nitro_new.1$comp_gene
# colnames(nitro_new.1)
nitro_new <- nitro_new.dat[,-1]
rownames(nitro_new)
# rownames(nitro_new) <- nitro_new.dat$newname
# colnames(nitro_new)
# rownames(nitro_new)
nitro_dds.new <- DESeqDataSetFromMatrix(countData = nitro_new,
                                         colData = metadata.1,
                                         design = ~ Nitrogen+Cover+Nitrogen:Cover)
nrow(nitro_dds.new)#312

colData(nitro_dds.new)
#levels(ddsMF$Cover)#
nitro_dds.new <- DESeq(nitro_dds.new)
nrow(nitro_dds.new)
sizeFactors(nitro_dds.new)
plotDispEsts(nitro_dds.new)
resultsNames(nitro_dds.new)
#under no cover, the differences between n60 and n0
res_nitrogen_nocover<- results(nitro_dds.new,contrast=c("Nitrogen","N60","N0"),lfcThreshold = 0.01)
summary(res_nitrogen_nocover)#sig 0 0
#undervetch, the differences between n60 and n0
res_nitrogen<- results(nitro_dds.new,list(c("Nitrogen_N60_vs_N0","NitrogenN60.CoverVetch")),lfcThreshold = 0.01)#sig
summary(res_nitrogen)# sig 1 1
#under n60, the differences between vetch and no cover
res_cover<- results(nitro_dds.new,list(c("Cover_Vetch_vs_No_cover","NitrogenN60.CoverVetch")),lfcThreshold = 0.01)#not sig
summary(res_cover)#sig0  1
#under n0, the differences between vetch and no cover
res_cover_n0<- results(nitro_dds.new,contrast=c("Cover","Vetch","No_cover"),lfcThreshold = 0.01)#no
summary(res_cover_n0)#sig 0 0
# the differences between vetch with fertilzation and no cover with no fertlization, we don't consider it here
res_all<- results(nitro_dds.new,lfcThreshold = 0.01)#sig
summary(res_all)

res_nitrogen1 <- as.data.frame(res_nitrogen)
res_nitrogen2 <- mutate(res_nitrogen1,sig=ifelse(res_nitrogen1$padj<0.1,"FDR<0.1","Not Sig"))
res_nitrogen2[which(abs(res_nitrogen2$log2FoldChange)<4.5),"sig"]="Not Sig" #greater than 0.1 keep it as significant, otherwise not significant
res_nitrogen2
res_nitrogen.sig <- subset(res_nitrogen2,sig=="FDR<0.1")
res_nitrogen.sig.new <- mutate(res_nitrogen.sig,gene_name=rownames(res_nitrogen.sig))
dim(res_nitrogen.sig.new)

res_nitrogen_nocover1 <- as.data.frame(res_nitrogen_nocover)
res_nitrogen_nocover2 <- mutate(res_nitrogen_nocover1,sig=ifelse(res_nitrogen_nocover1$padj<0.1,"FDR<0.1","Not Sig"))
res_nitrogen_nocover2[which(abs(res_nitrogen_nocover2$log2FoldChange)<4.5),"sig"]="Not Sig" #greater than 0.1 keep it as significant, otherwise not significant
res_nitrogen_nocover2
res_nitrogen_nocover.sig <- subset(res_nitrogen_nocover2,sig=="FDR<0.1")
res_nitrogen_nocover.sig.new <- mutate(res_nitrogen_nocover.sig,gene_name=rownames(res_nitrogen_nocover.sig))
dim(res_nitrogen_nocover.sig.new)

res_cover1 <- as.data.frame(res_cover)
res_cover2 <- mutate(res_cover1,sig=ifelse(res_cover1$padj<0.1,"FDR<0.1","Not Sig"))
res_cover2[which(abs(res_cover2$log2FoldChange)<4.5),"sig"]="Not Sig" #greater than 0.1 keep it as significant, otherwise not significant
res_cover2
res_cover.sig <- subset(res_cover2,sig=="FDR<0.1")
res_cover.sig.new <- mutate(res_cover.sig,gene_name=rownames(res_cover.sig))
dim(res_cover.sig.new)

res_cover_n01 <- as.data.frame(res_cover_n0)
res_cover_n02 <- mutate(res_cover_n01,sig=ifelse(res_cover_n01$padj<0.1,"FDR<0.1","Not Sig"))
res_cover_n02[which(abs(res_cover_n02$log2FoldChange)<4.5),"sig"]="Not Sig" #greater than 0.1 keep it as significant, otherwise not significant
res_cover_n02
res_cover_n0.sig <- subset(res_cover_n02,sig=="FDR<0.1")
res_cover_n0.sig.new <- mutate(res_cover_n0.sig,gene_name=rownames(res_cover_n0.sig))
dim(res_cover.sig.new)

res_all1 <- as.data.frame(res_all)
res_all2 <- mutate(res_all1,sig=ifelse(res_all1$padj<0.1,"FDR<0.1","Not Sig"))
res_all2[which(abs(res_all2$log2FoldChange)<4.5),"sig"]="Not Sig" #greater than 0.1 keep it as significant, otherwise not significant
res_all2
res_all.sig <- subset(res_all2,sig=="FDR<0.1")
res_all.sig.new <- mutate(res_all.sig,gene_name=rownames(res_all.sig))
dim(res_all.sig.new)
######extract the gene I need#######3
all_results <- as.data.frame(rbind(res_nitrogen.sig.new,res_nitrogen_nocover.sig.new,
                                   res_cover.sig.new, res_cover_n0.sig.new,res_all.sig.new))
genes <- as.data.frame(unique(all_results$gene_name))
colnames(genes) <- "newname"
dim(genes)
genes.data <- left_join(genes,nitro_all_gene.dat,by="newname")
genes.data.all <-genes.data #left_join(genes.data,nitrogen.deseq[,c(2,6)],by="contig")
dim(genes.data.all)
colnames(genes.data.all)
rownames(genes.data.all)
glimpse(genes.data.all)
#genes.data.all$NCNTN60_2 <- as.integer(genes.data.all$NCNTN60_2)
genes.data.all.group <- data.frame(genes.data.all[,c(3:15)])%>% 
  group_by(comp_gene) %>% 
  summarise(across(everything(), sum))
dim(genes.data.all.group)
unique(genes.data.all.group$comp_gene)
# new.cate.dat.amo <- apply(genes.data.all.group.1[1:3,2:13],2,sum)
# new.cate.dat.nap <-apply(genes.data.all.group.1[c(4),2:13],2,sum)
# new.cate.dat.nar <- apply(genes.data.all.group.1[c(5:7),2:13],2,sum)
# new.cate.dat.nas <- apply(genes.data.all.group.1[c(8),2:13],2,sum)
# new.cate.dat.nira <- apply(genes.data.all.group.1[c(9),2:13],2,sum)
# new.cate.dat.nirb <- apply(genes.data.all.group.1[c(10),2:13],2,sum)
# new.cate.dat.nirk <- apply(genes.data.all.group.1[c(11),2:13],2,sum)
# new.cate.dat.norb <- apply(genes.data.all.group.1[c(12),2:13],2,sum)
# new.cate.dat.nrf <- apply(genes.data.all.group.1[c(13:14),2:13],2,sum)
# new.cate.dat.nrt <- apply(genes.data.all.group.1[c(15),2:13],2,sum)
# genes.data.all.group <- rbind(new.cate.dat.amo,new.cate.dat.nap,new.cate.dat.nar,new.cate.dat.nas,
#                               new.cate.dat.nira,new.cate.dat.nirb,new.cate.dat.nirk,new.cate.dat.norb,
#                               new.cate.dat.nrf,new.cate.dat.nrt)
# rownames(genes.data.all.group) <- c("amoABC","nap","narGHI","nasA","nirA","nirB","nirK","norB","nrfAH","NRT")
colnames(genes.data.all.group)
#head(new.cate.dat.amo)
colnames(genes.data.all.group)
genes.data.all.group$comp_gene
#genes.data.all.group.1 <- genes.data.all.group %>% mutate(process=c(rep("ammonia_hydroxylamine",3),rep("nitrate_nitrite",3),))

write.csv(genes.data.all.group,"dataset/genedata.csv")#could be use as the supplyment material
library(DESeq2)
genes.data.group <- as.data.frame(genes.data.all.group)[,-c(11)]
rownames(genes.data.group) <- (genes.data.group)$comp_gene
rownames(genes.data.group)
new.dat2 <- genes.data.group[,-1]
#new.dat2 <- genes.data.all.group[,-c(1,10)]
# rownames(new.dat2)[c(6,7,8,16)] <- c("narG","narH","narI","NRT")
# rownames(new.dat2)
dds.new2 <- DESeqDataSetFromMatrix(countData = new.dat2[-15,],
                                   colData = metadata.1[-c(10),],
                                   design = ~ Nitrogen*Cover)
nrow(dds.new2)#15

#####deseq2###3
colData(dds.new2)
#levels(ddsMF$Cover)#
dds.new2 <- DESeq(dds.new2)
nrow(dds.new2)
sizeFactors(dds.new2)
plotDispEsts(dds.new2)
resultsNames(dds.new2)
# library(pheatmap)
# select <- order(rowMeans(counts(dds.new2,normalized=TRUE)),
#                 decreasing=TRUE)
# df <- as.data.frame(colData(dds.new2)[,c("Cover","Nitrogen")])
# pheatmap(assay(rld2), cluster_rows=T, show_rownames=T,
#          cluster_cols=T, annotation_col=df)

res_nitrogen_nocover2<- results(dds.new2,contrast=c("Nitrogen","N60","N0"),lfcThreshold = 0.01)
summary(res_nitrogen_nocover2)#sig 1 1
#undervetch, the differences between n60 and n0
res_nitrogen2<- results(dds.new2,list(c("Nitrogen_N60_vs_N0","NitrogenN60.CoverVetch")),lfcThreshold = 0.01)#sig
summary(res_nitrogen2)# sig 3 3
#under n60, the differences between vetch and no cover
res_cover2<- results(dds.new2,list(c("Cover_Vetch_vs_No_cover","NitrogenN60.CoverVetch")),lfcThreshold = 0.01)#not sig
summary(res_cover2)#sig0  1
#under n0, the differences between vetch and no cover
res_cover_n02<- results(dds.new2,contrast=c("Cover","Vetch","No_cover"),lfcThreshold = 0.01)#no
summary(res_cover_n02)#sig 0 0
# # # the differences between vetch with fertilzation and no cover with no fertlization, we don't consider it here
# res_all2<- results(dds.new2,lfcThreshold = 0.01)#sig
# summary(res_all2)

res2_nitrogen_nocover1 <- as.data.frame(res_nitrogen_nocover2)
res2_nitrogen_nocover2 <- mutate(res2_nitrogen_nocover1,sig=ifelse(res2_nitrogen_nocover1$padj<0.1,"FDR<0.1","Not Sig"))
res2_nitrogen_nocover2[which(abs(res2_nitrogen_nocover2$log2FoldChange)<0.1),"sig"]="Not Sig" #greater than 0.1 keep it as significant, otherwise not significant
res2_nitrogen_nocover2
res2_nitrogen_nocover.sig <- subset(res2_nitrogen_nocover2,sig=="FDR<0.1")
res2_nitrogen_nocover.sig.new <- mutate(res2_nitrogen_nocover.sig,gene_name=rownames(res2_nitrogen_nocover.sig))
write.csv(res2_nitrogen_nocover.sig.new,"dataset/res2_nitrogen_nocover.sig.new.csv")
res2.nocover.plot <- 
  ggplot(data = res2_nitrogen_nocover.sig.new,
         aes(x = reorder(gene_name,log2FoldChange), y = log2FoldChange, fill=log2FoldChange>0))+
  scale_fill_manual(name="Treatment",values=c("#e3c232","grey"),label=c("NCN60","NCN0"))+
  geom_bar(stat = "identity",width = 0.7)+
  #ylim(-1.5,1)+
  coord_flip()+
  labs(title = "NCN60 vs. NCN0",x = "Genes", y = "log2FoldChange")+
  #guides(fill = FALSE)+
  theme(legend.position="right",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_line(colour = "#f0f0f0",size=0.75),
        # panel.grid.minor = element_blank(),
        plot.title=(element_text(size=15,family = "Arial",face="plain",hjust = 0,vjust = 3)),
        text=element_text(family = "Arial",face="plain"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text=element_text(size=15,family="Arial"),
        axis.title=element_text(size = 15,face="plain",family = "Arial"),
        legend.title = element_text(size=15,face="plain",family = "Arial"),
        legend.text = (element_text(size=15,family = "Arial")))

res2_nitrogen1 <- as.data.frame(res_nitrogen2)
res2_nitrogen2 <- mutate(res2_nitrogen1,sig=ifelse(res2_nitrogen1$padj<0.1,"FDR<0.1","Not Sig"))
res2_nitrogen2[which(abs(res2_nitrogen2$log2FoldChange)<0.1),"sig"]="Not Sig" #greater than 0.1 keep it as significant, otherwise not significant
res2_nitrogen2
res2_nitrogen.sig <- subset(res2_nitrogen2,sig=="FDR<0.1")
res2_nitrogen.sig.new <- mutate(res2_nitrogen.sig,gene_name=rownames(res2_nitrogen.sig))
write.csv(res2_nitrogen.sig.new,"dataset/res2_nitrogen.sig.new.csv")

res2.cover.plot <-
  ggplot(data = res2_nitrogen.sig.new,
         aes(x = reorder(gene_name,log2FoldChange), y = log2FoldChange, fill=log2FoldChange>0))+
  scale_fill_manual(name="Treatment",values=c("grey","#e3c232"),label=c("VN0","VN60"))+
  geom_bar(stat = "identity",width = 0.7)+
  #ylim(-1.5,1)+
  coord_flip()+
  labs(title = "VN60 vs. VN0",x = "", y = "log2FoldChange")+
  #guides(fill = FALSE)+
  theme(legend.position="right",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_line(colour = "#f0f0f0",size=0.75),
        # panel.grid.minor = element_blank(),
        plot.title=(element_text(size=15,family = "Arial",face="plain",hjust = 0,vjust = 3)),
        text=element_text(family = "Arial",face="plain"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text=element_text(size=15,family="Arial"),
        axis.title=element_text(size = 15,face="plain",family = "Arial"),
        legend.title = element_text(size=15,face="plain",family = "Arial"),
        legend.text = (element_text(size=15,family = "Arial")))

library(dplyr)
res2_cover1 <- as.data.frame(res_cover2)
res2_cover2 <- mutate(res2_cover1,sig=ifelse(res2_cover1$padj<0.1,"FDR<0.1","Not Sig"))
res2_cover2[which(abs(res2_cover2$log2FoldChange)<0.1),"sig"]="Not Sig" #greater than 0.1 keep it as significant, otherwise not significant
res2_cover2
res2_cover.sig <- subset(res2_cover2,sig=="FDR<0.1")
res2_cover.sig.new <- mutate(res2_cover.sig,gene_name=rownames(res2_cover.sig))
dim(res2_cover.sig.new)
write.csv(res2_cover.sig.new,"dataset/res2_cover.sig.new.csv")

res2_n60.sig.plot <-
  ggplot(data = res2_cover.sig.new,
         aes(x = reorder(gene_name,log2FoldChange), y = log2FoldChange, fill=log2FoldChange>0))+
  scale_fill_manual(name="Regulation",values=c("grey","#91b04a"),label=c("Down","Up"))+
  #ylim(-1.5,1)+
  geom_bar(stat = "identity",width = 0.7)+
  coord_flip()+
  labs(title = "VN60 (Treatment) vs. NCN60",x = "", y = "log2FoldChange")+
  #guides(fill = FALSE)+
  theme(legend.position="right",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_line(colour = "#f0f0f0",size=0.75),
        # panel.grid.minor = element_blank(),
        plot.title=(element_text(size=15,family = "Arial",face="plain",hjust = 0,vjust = 3)),
        text=element_text(family = "Arial",face="plain"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text=element_text(size=15,family="Arial"),
        axis.title=element_text(size = 15,face="plain",family = "Arial"),
        legend.title = element_text(size=15,face="plain",family = "Arial"),
        legend.text = (element_text(size=15,family = "Arial")))


res2_cover_n01 <- as.data.frame(res_cover_n02)
res2_cover_n02 <- mutate(res2_cover_n01,sig=ifelse(res2_cover_n01$padj<0.1,"FDR<0.1","Not Sig"))
res2_cover_n02[which(abs(res2_cover_n02$log2FoldChange)<0.1),"sig"]="Not Sig" #greater than 0.1 keep it as significant, otherwise not significant
res2_cover_n02
res2_cover_n0.sig <- subset(res2_cover_n02,sig=="FDR<0.1")
res2_cover_n0.sig.new <- mutate(res2_cover_n0.sig,gene_name=rownames(res2_cover_n0.sig))
dim(res2_cover_n0.sig.new)
write.csv(res2_cover_n0.sig.new,"dataset/res2_cover_n0.sig.new.csv")

res2_n0.sig.plot <-
  ggplot(data = res2_cover_n0.sig.new,
         aes(x = reorder(gene_name,log2FoldChange), y = log2FoldChange, fill=log2FoldChange>0))+
  scale_fill_manual(name="Treatment",values=c("grey","#91b04a"),label=c("NCN0","VN0"))+
  #ylim(-1.5,1)+
  geom_bar(stat = "identity",width = 0.7)+
  coord_flip()+
  labs(title = "VN0 vs. NCN0",x = "Genes", y = "log2FoldChange")+
  #guides(fill = FALSE)+
  theme(legend.position="right",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_line(colour = "#f0f0f0",size=0.75),
        # panel.grid.minor = element_blank(),
        plot.title=(element_text(size=15,family = "Arial",face="plain",hjust = 0,vjust = 3)),
        text=element_text(family = "Arial",face="plain"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text=element_text(size=15,family="Arial"),
        axis.title=element_text(size = 15,face="plain",family = "Arial"),
        legend.title = element_text(size=15,face="plain",family = "Arial"),
        legend.text = (element_text(size=15,family = "Arial")))

#RDA
library(ggplot2)
library(phyloseq)
library(vegan)
library(ggforce)
library(MASS)
library(car)
library(lme4)
library(tidyverse)
library(fBasics)
rld2 <- rlog(dds.new2)
gene.phyloseq <- merge_phyloseq(sample_data(metadata.1), otu_table(as.data.frame(assay(rld2)),taxa_are_rows = TRUE))
gene.phyloseq
# formularda <- rda(t(otu_table(gene.phyloseq)) ~ Cover+Nitrogen+Cover*Nitrogen+pH+GWC+NH4+NO3+
#                     AG+BG+CB+LAP+NAG+PHOS+XYL,data = metadata.up)
# summary(formularda)
rda_data <- ordinate(gene.phyloseq,method = "RDA",formula = ~Nitrogen+Cover)
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
adonis(t(otu_table(gene.phyloseq))~Nitrogen+Cover,
       data=data.frame(sample_data(gene.phyloseq)),
       permutations=9999,by="margin")#only nitrogen significant

df_pca <- prcomp(t(otu_table(gene.phyloseq)))
pca_tab <- as.data.frame(df_pca$x)
pca_tab$Nitrogen <- metadata[-10,]$Nitrogen
pca_tab$Cover <- metadata[-10,]$Cover

rda.plot <- 
  plot_ordination(gene.phyloseq,rda_data,type="sample",color="Nitrogen",shape = "Cover")+
  geom_point(size=3)+
  geom_text(aes(label=rownames(pca_tab),vjust=-2,hjust=0.5),size=2)+
  ggtitle("N genes")+
  scale_color_manual(values = c("#106494", "#bf6a0a"))+
  scale_shape_manual(values = c(17,8))+
  #xlab(label = "PC1 [53.1%]")+
  #ylab(label = "PC2 [20%]")+
  geom_mark_ellipse(aes(color = Nitrogen, group=Nitrogen),expand=0)+
  geom_vline(xintercept = c(0), size=0.2,linetype=2)+
  geom_hline(yintercept = c(0), size=0.2,linetype=2)+
  theme(legend.position="bottom",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_line(colour="grey",size = 0.1),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=15,family = "Arial",face="plain",vjust = 3,hjust=0)),
        text=element_text(family = "Arial",face="plain",size = 15),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(size=15,family="Arial",face = "plain"),
        axis.title=element_text(size =15,face="plain",family = "Arial",vjust = 1),
        legend.title = element_text(size=15,face="plain",family = "Arial"),
        legend.text = (element_text(size=15,family = "Arial")))#+
#  annotate(geom="text", x=1.32, y=-0.42, size=8,label="*",color="red", fontface="bold")
library(ggpubr)
tiff('figures/allgene_cate_cazy_divergent.tiff', units="in", width=15, height=15, res=300)
ggarrange(res2.nocover.plot,res2.cover.plot,res2_n0.sig.plot,rda.plot,
          labels = c("A","B","C","D"),ncol = 2, nrow = 2,common.legend = F ,legend = "bottom", widths=c(1,1))
dev.off()

library(vegan)
library(ggplot2)
arrowmat <- vegan::scores(rda_data, display = "bp")
head(arrowmat)
#arrowmat.sub <-arrowmat[-c(1,),] 
# Add labels, make a data.frame
arrowdf <- data.frame(labels = c("N60","GWC","NH4","NO3"), arrowmat)
# Define the arrow aesthetic mapping
arrow_map <- aes(xend = RDA1, yend = RDA2, x = 0,  y = 0,  shape = NULL,color = NULL, label = labels)
label_map <- aes(x = 1.3*RDA1, y = 1.3*RDA2, shape = NULL, color = NULL,label = labels)
arrowhead = arrow(length = unit(0.01, "npc"))
rda.witharrow <- rda.plot+
  geom_segment(mapping = arrow_map, size = .5,data = arrowdf,color = "black",arrow = arrowhead )+
  geom_text(mapping = label_map, size = 5,
            data = arrowdf,  show.legend = FALSE,
            color="black",
            #nudge_x = -0.01,nudge_y = 0.15,
            hjust=0.5*(1-sign(arrowdf$RDA1)),vjust=0.5*(1-sign(arrowdf$RDA2)))
rda.witharrow

# library(ggpubr)
# tiff('figures/rda.tiff', units="in", width=18, height=10, res=300)
# ggarrange(res2.nocover.plot,rda.witharrow,
#           labels = c("A","B"),ncol = 2, nrow = 1,common.legend = F ,legend = "bottom", widths=c(1,1)) 
# dev.off()

# library(corrplot)#correlation with env
# rownames(nitrogen.deseq.cate)
# colnames(nitrogen.deseq.cate)
# library(DESeq2)
# library(tidyverse)
# gene_tab <- as.data.frame(assay(rld2))
# gene_tab.t <- as.data.frame(t(gene_tab))
# gene_tab.t.tab <- mutate(gene_tab.t,sample=rownames(gene_tab.t))
# head(gene_tab.t.tab)
# metadata.new <- mutate(metadata,sample=rownames(metadata))
# tab <- left_join(metadata.new,gene_tab.t.tab)
# normalTest(log(tab$norB))
# nitrogen.deseq.cate.n60 <- subset(tab,Nitrogen=="N60")
# colnames(nitrogen.deseq.cate)
# colnames(tab)
# #N_gene <- as.matrix((nitrogen.deseq.cate[,c(40,36,37,25,32,28,34)]))
# N_gene <- as.data.frame((nitrogen.deseq.cate.n60[,c(35,32,33,26,30,27,31)]))[-c(10),]
# glimpse(N_gene)
# normalTest((N_gene$norB))
# N_gene.cor <- cor((N_gene))
# N_gene.res1 <- cor.mtest((N_gene), conf.level = .95)
# pAdj <- p.adjust(c(N_gene.res1[[1]]), method = "BH")
# resAdj <- matrix(pAdj, ncol = dim(N_gene.res1[[1]])[1])
# corrplot(N_gene.cor,type = "lower",na.label = "NA",p.mat = resAdj, sig.level = .05,insig = "blank",tl.col = "black",tl.srt = 45)
# N_gene.cor.plot

# normalTest(N_gene$amoC[-c(1,2,6,7,8)])
# amoc <- N_gene$amoC[-c(1,2,6,7,8)]
# 
# cor.test(N_gene$amoC[-c(1,2,6,7,8)],N_gene$norB[-c(1,2,6,7,8)])
# cor.test(nitrogen.deseq.cate$nirA,(nitrogen.deseq.cate$nirB))

#####denitrification and dissimlatory nitrate reduction#####



