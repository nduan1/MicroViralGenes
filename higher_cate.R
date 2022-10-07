# colnames(carb_all_gene.dat)
# starch_new.dat <- carb_all_gene.dat[,c(16,2:13)]
# rownames(starch_new.dat) <- starch_new.dat$newname
# starch_new <- starch_new.dat[,-1]
# 
# starch_dds.newhigher <- DESeqDataSetFromMatrix(countData = starch_new,
#                                          colData = metadata,
#                                          design = ~ Nitrogen*Cover)
# nrow(starch_dds.newhigher)#39
# 
# colData(starch_dds.new)
# #levels(ddsMF$Cover)#
# starch_dds.newhigher <- DESeq(starch_dds.newhigher)
# nrow(starch_dds.newhigher)
# sizeFactors(starch_dds.newhigher)
# plotDispEsts(starch_dds.newhigher)
# resultsNames(starch_dds.newhigher)
# #under no cover, the differences between n60 and n0
# res_nitrogen_nocoverhigher<- results(starch_dds.newhigher,contrast=c("Nitrogen","N60","N0"),lfcThreshold = 0.01)
# summary(res_nitrogen_nocoverhigher)#sig 0 0
# #undervetch, the differences between n60 and n0
# res_nitrogenhigher<- results(starch_dds.newhigher,list(c("Nitrogen_N60_vs_N0","NitrogenN60.CoverVetch")),lfcThreshold = 0.01)#sig
# summary(res_nitrogenhigher)# sig 1 1
# #under n60, the differences between vetch and no cover
# res_coverhigher<- results(starch_dds.newhigher,list(c("Cover_Vetch_vs_No_cover","NitrogenN60.CoverVetch")),lfcThreshold = 0.01)#not sig
# summary(res_cover)#sig0  1
# #under n0, the differences between vetch and no cover
# res_cover_n0higher<- results(starch_dds.newhigher,contrast=c("Cover","Vetch","No_cover"),lfcThreshold = 0.01)#no
# summary(res_cover_n0higher)#sig 0 0
# # the differences between vetch with fertilzation and no cover with no fertlization, we don't consider it here
# res_allhigher<- results(starch_dds.newhigher,lfcThreshold = 0.01)#sig
# summary(res_allhigher)
# 
# 
# library(tidyverse)
# res_nitrogen1higher <- as.data.frame(res_nitrogenhigher)
# res_nitrogen2higher <- mutate(res_nitrogen1higher,sig=ifelse(res_nitrogen1higher$padj<0.1,"FDR<0.01","Not Sig"))
# res_nitrogen2higher[which(abs(res_nitrogen2higher$log2FoldChange)<2),"sig"]="Not Sig" #greater than 0.1 keep it as significant, otherwise not significant
# res_nitrogen2higher
# res_nitrogen.sighigher <- subset(res_nitrogen2higher,sig=="FDR<0.01")
# res_nitrogen.sig.newhigher <- mutate(res_nitrogen.sighigher,gene_name=rownames(res_nitrogen.sighigher))
# dim(res_nitrogen.sig.newhigher)
# 
# res_nitrogen_nocover1higher <- as.data.frame(res_nitrogen_nocoverhigher)
# res_nitrogen_nocover2higher <- mutate(res_nitrogen_nocover1higher,sig=ifelse(res_nitrogen_nocover1higher$padj<0.1,"FDR<0.01","Not Sig"))
# res_nitrogen_nocover2higher[which(abs(res_nitrogen_nocover2higher$log2FoldChange)<2),"sig"]="Not Sig" #greater than 0.1 keep it as significant, otherwise not significant
# res_nitrogen_nocover2higher
# res_nitrogen_nocover.sighigher <- subset(res_nitrogen_nocover2higher,sig=="FDR<0.01")
# res_nitrogen_nocover.sig.newhigher <- mutate(res_nitrogen_nocover.sighigher,gene_name=rownames(res_nitrogen_nocover.sighigher))
# dim(res_nitrogen_nocover.sig.newhigher)
# 
# res_cover1higher <- as.data.frame(res_coverhigher)
# res_cover2higher <- mutate(res_cover1higher,sig=ifelse(res_cover1higher$padj<0.1,"FDR<0.01","Not Sig"))
# res_cover2higher[which(abs(res_cover2higher$log2FoldChange)<2),"sig"]="Not Sig" #greater than 0.1 keep it as significant, otherwise not significant
# res_cover2higher
# res_cover.sighigher <- subset(res_cover2higher,sig=="FDR<0.01")
# res_cover.sig.newhigher <- mutate(res_cover.sighigher,gene_name=rownames(res_cover.sighigher))
# dim(res_cover.sig.newhigher)
# 
# res_cover_n01higher <- as.data.frame(res_cover_n0higher)
# res_cover_n02higher <- mutate(res_cover_n01higher,sig=ifelse(res_cover_n01higher$padj<0.1,"FDR<0.01","Not Sig"))
# res_cover_n02higher[which(abs(res_cover2higher$log2FoldChange)<2),"sig"]="Not Sig" #greater than 0.1 keep it as significant, otherwise not significant
# res_cover_n02higher
# res_cover_n0.sighigher <- subset(res_cover_n02higher,sig=="FDR<0.01")
# res_cover_n0.sig.newhigher <- mutate(res_cover_n0.sighigher,gene_name=rownames(res_cover_n0.sighigher))
# dim(res_cover_n0.sig.newhigher)

# res_all_n01 <- as.data.frame(res_all)
# res_all_n02higher <- mutate(res_all_n01,sig=ifelse(res_all_n01$padj<0.1,"FDR<0.01","Not Sig"))
# res_all_n02higher[which(abs(res_all2higher$log2FoldChange)<2),"sig"]="Not Sig" #greater than 0.1 keep it as significant, otherwise not significant
# res_all_n02higher
# res_all_n0.sighigher <- subset(res_all_n02higher,sig=="FDR<0.01")
# res_all_n0.sig.newhigher <- mutate(res_all_n0.sighigher,gene_name=rownames(res_all_n0.sighigher))
# dim(res_all_n0.sig.newhigher)
# ###################
# all_resultshigher <- as.data.frame(rbind(res_nitrogen.sig.newhigher,
#                                          res_nitrogen_nocover.sig.newhigher,
#                                   res_cover.sig.newhigher, 
#                                   res_cover_n0.sig.newhigher))
# geneshigher <- as.data.frame(unique(all_resultshigher$gene_name))
# colnames(geneshigher) <- "newname"
# dim(geneshigher)
# genes.datahigher <- left_join(geneshigher,carb_all_gene.dat,by="newname")
# genes.data.allhigher <- left_join(genes.datahigher,carb_all_gene[,c(1,17)],by="contig")
# colnames(carb_all_gene)
# dim(genes.data.allhigher)
# colnames(genes.data.allhigher)
# rownames(genes.data.allhigher)
# genes.data.allhigher$NCNTN60_2 <- as.integer(genes.data.allhigher$NCNTN60_2)
# ######________________________________________
# genes.data.all.group.new <- data.frame(genes.data.all[,c(3:14,17)])%>% 
#   group_by(group) %>% 
#   summarise(across(everything(), sum))
# dim(genes.data.all.group.new)
# colnames(genes.data.all.group.new)
# library(DESeq2)
# genes.data.group.new<- as.data.frame(genes.data.all.group.new)
# rownames(genes.data.group.new) <- (genes.data.group.new)$group
# rownames(genes.data.group.new)
# new.dat2.new <- genes.data.group.new[,-1]
# dds.new2.new <- DESeqDataSetFromMatrix(countData = new.dat2.new,#[,-c(5)],
#                                    colData = metadata,#[-c(5),],
#                                    design = ~ Nitrogen*Cover)
# # nrow(dds.new2.new)#24
# # #####deseq2###3
# # colData(dds.new2.new)
# # #levels(ddsMF$Cover)#
# dds.new2.new <- DESeq(dds.new2.new)
# nrow(dds.new2.new)
# sizeFactors(dds.new2.new)
# plotDispEsts(dds.new2.new)
# resultsNames(dds.new2.new)
# 
# res_nitrogen_nocover2.new<- results(dds.new2.new,contrast=c("Nitrogen","N60","N0"),lfcThreshold = 0.01)
# summary(res_nitrogen_nocover2.new)#sig 1 1
# res_nitrogen2.new<- results(dds.new2.new,list(c("Nitrogen_N60_vs_N0","NitrogenN60.CoverVetch")),lfcThreshold = 0.01)#sig
# summary(res_nitrogen2.new)# sig 3 3
# res_cover2.new<- results(dds.new2.new,list(c("Cover_Vetch_vs_No_cover","NitrogenN60.CoverVetch")),lfcThreshold = 0.01)#not sig
# summary(res_cover2.new)#sig0  1
# res_cover_n02.new<- results(dds.new2.new,contrast=c("Cover","Vetch","No_cover"),lfcThreshold = 0.01)#no
# summary(res_cover_n02.new)#sig 0 0
# 
# res2_nitrogen_nocover1.new <- as.data.frame(res_nitrogen_nocover2.new)
# res2_nitrogen_nocover2.new<- mutate(res2_nitrogen_nocover1.new,sig=ifelse(res2_nitrogen_nocover1.new$padj<0.1,"FDR<0.1","Not Sig"))
# res2_nitrogen_nocover2.new[which(abs(res2_nitrogen_nocover2.new$log2FoldChange)<0.1),"sig"]="Not Sig" #greater than 0.1 keep it as significant, otherwise not significant
# res2_nitrogen_nocover2.new
# res2_nitrogen_nocover.sig.new <- subset(res2_nitrogen_nocover2.new,sig=="FDR<0.1")
# res2_nitrogen_nocover.sig.new.new <- mutate(res2_nitrogen_nocover.sig.new,gene_name=rownames(res2_nitrogen_nocover.sig.new))
# dim(res2_nitrogen_nocover.sig.new.new)
# res2.nocover.plot.new <- 
#   ggplot(data = res2_nitrogen_nocover.sig.new.new,
#          aes(x = reorder(gene_name,log2FoldChange), y = log2FoldChange, fill=log2FoldChange>0))+
#   scale_fill_manual(name="Treatment",values=c("grey","#e3c232"),label=c("NCN0","NCN60"))+
#   geom_bar(stat = "identity",width = 0.7)+
#   #ylim(-1.5,1)+
#   coord_flip()+
#   labs(title = "",x = "Genes", y = "log2FoldChange")+
#   #guides(fill = FALSE)+
#   theme(legend.position="right",panel.background = element_rect(fill = "white",color="black"),
#         panel.grid.major = element_line(colour = "#f0f0f0",size=0.75),
#         # panel.grid.minor = element_blank(),
#         plot.title=(element_text(size=15,family = "Arial",face="plain",hjust = 0,vjust = 3)),
#         text=element_text(family = "Arial",face="plain"),
#         panel.border = element_rect(colour = "black", fill=NA, size=0.5),
#         axis.text=element_text(size=15,family="Arial"),
#         axis.title=element_text(size = 15,face="plain",family = "Arial"),
#         legend.title = element_text(size=15,face="plain",family = "Arial"),
#         legend.text = (element_text(size=15,family = "Arial")))
# 
# res2_nitrogen1.new <- as.data.frame(res_nitrogen2.new)
# res2_nitrogen2.new <- mutate(res2_nitrogen1.new,sig=ifelse(res2_nitrogen1.new$padj<0.1,"FDR<0.1","Not Sig"))
# res2_nitrogen2.new[which(abs(res2_nitrogen2.new$log2FoldChange)<0.1),"sig"]="Not Sig" #greater than 0.1 keep it as significant, otherwise not significant
# res2_nitrogen2.new
# res2_nitrogen.sig.new <- subset(res2_nitrogen2.new,sig=="FDR<0.1")
# res2_nitrogen.sig.new.new<- mutate(res2_nitrogen.sig.new,gene_name=rownames(res2_nitrogen.sig.new))
# 
# res2.cover.plot.new <- 
#   ggplot(data = res2_nitrogen.sig.new.new,
#          aes(x = reorder(gene_name,log2FoldChange), y = log2FoldChange, fill=log2FoldChange>0))+
#   scale_fill_manual(name="Treatment",values=c("grey","#e3c232"),label=c("VN0","VN60"))+
#   geom_bar(stat = "identity",width = 0.7)+
#   #ylim(-1.5,1)+
#   coord_flip()+
#   labs(title = "",x = "", y = "log2FoldChange")+
#   #guides(fill = FALSE)+
#   theme(legend.position="right",panel.background = element_rect(fill = "white",color="black"),
#         panel.grid.major = element_line(colour = "#f0f0f0",size=0.75),
#         # panel.grid.minor = element_blank(),
#         plot.title=(element_text(size=15,family = "Arial",face="plain",hjust = 0,vjust = 3)),
#         text=element_text(family = "Arial",face="plain"),
#         panel.border = element_rect(colour = "black", fill=NA, size=0.5),
#         axis.text=element_text(size=15,family="Arial"),
#         axis.title=element_text(size = 15,face="plain",family = "Arial"),
#         legend.title = element_text(size=15,face="plain",family = "Arial"),
#         legend.text = (element_text(size=15,family = "Arial")))
# 
# library(dplyr)
# res2_cover1.new <- as.data.frame(res_cover2.new)
# res2_cover2.new <- mutate(res2_cover1.new,sig=ifelse(res2_cover1.new$padj<0.1,"FDR<0.1","Not Sig"))
# res2_cover2.new[which(abs(res2_cover2.new$log2FoldChange)<0.1),"sig"]="Not Sig" #greater than 0.1 keep it as significant, otherwise not significant
# res2_cover2.new
# res2_cover.sig.new <- subset(res2_cover2.new,sig=="FDR<0.1")
# res2_cover.sig.new.new <- mutate(res2_cover.sig.new,gene_name=rownames(res2_cover.sig.new))
# dim(res2_cover.sig.new.new)
# 
# res2_n60.sig.plot.new <-
#   ggplot(data = res2_cover.sig.new.new,
#          aes(x = reorder(gene_name,log2FoldChange), y = log2FoldChange, fill=log2FoldChange>0))+
#   scale_fill_manual(name="Treatment",values=c("grey","#91b04a"),label=c("NCN60","VN60"))+
#   #ylim(-1.5,1)+
#   geom_bar(stat = "identity",width = 0.7)+
#   coord_flip()+
#   labs(title = "",x = "Genes", y = "log2FoldChange")+
#   #guides(fill = FALSE)+
#   theme(legend.position="right",panel.background = element_rect(fill = "white",color="black"),
#         panel.grid.major = element_line(colour = "#f0f0f0",size=0.75),
#         # panel.grid.minor = element_blank(),
#         plot.title=(element_text(size=15,family = "Arial",face="plain",hjust = 0,vjust = 3)),
#         text=element_text(family = "Arial",face="plain"),
#         panel.border = element_rect(colour = "black", fill=NA, size=0.5),
#         axis.text=element_text(size=15,family="Arial"),
#         axis.title=element_text(size = 15,face="plain",family = "Arial"),
#         legend.title = element_text(size=15,face="plain",family = "Arial"),
#         legend.text = (element_text(size=15,family = "Arial")))
# 
# 
# res2_cover_n01.new <- as.data.frame(res_cover_n02.new)
# res2_cover_n02.new <- mutate(res2_cover_n01.new,sig=ifelse(res2_cover_n01.new$padj<0.1,"FDR<0.1","Not Sig"))
# res2_cover_n02.new[which(abs(res2_cover_n02.new$log2FoldChange)<0.1),"sig"]="Not Sig" #greater than 0.1 keep it as significant, otherwise not significant
# res2_cover_n02.new
# res2_cover_n0.sig.new <- subset(res2_cover_n02.new,sig=="FDR<0.1")
# res2_cover_n0.sig.new.new <- mutate(res2_cover_n0.sig.new,gene_name=rownames(res2_cover_n0.sig.new))
# dim(res2_cover_n0.sig.new.new)
# res2_n0.sig.plot.new <-
#   ggplot(data = res2_cover_n0.sig.new.new,
#          aes(x = reorder(gene_name,log2FoldChange), y = log2FoldChange, fill=log2FoldChange>0))+
#   scale_fill_manual(name="Treatment",values=c("grey","#91b04a"),label=c("NCN0","VN0"))+
#   #ylim(-1.5,1)+
#   geom_bar(stat = "identity",width = 0.7)+
#   coord_flip()+
#   labs(title = "",x = "Genes", y = "log2FoldChange")+
#   #guides(fill = FALSE)+
#   theme(legend.position="right",panel.background = element_rect(fill = "white",color="black"),
#         panel.grid.major = element_line(colour = "#f0f0f0",size=0.75),
#         # panel.grid.minor = element_blank(),
#         plot.title=(element_text(size=15,family = "Arial",face="plain",hjust = 0,vjust = 3)),
#         text=element_text(family = "Arial",face="plain"),
#         panel.border = element_rect(colour = "black", fill=NA, size=0.5),
#         axis.text=element_text(size=15,family="Arial"),
#         axis.title=element_text(size = 15,face="plain",family = "Arial"),
#         legend.title = element_text(size=15,face="plain",family = "Arial"),
#         legend.text = (element_text(size=15,family = "Arial")))
# 
# tiff('figures/supplygene_cate_cazy_divergent.tiff', units="in", width=17, height=17, res=300)
# ggarrange(res2.nocover.plot.new,res2.cover.plot.new,res2_n0.sig.plot.new,
#           labels = c("A","B","C"),ncol = 2, nrow = 2,common.legend = F,legend = "bottom", widths=c(1,1),heights = c(1,1)) 
# dev.off()
# 
# 
# ######bar plot
# library(DESeq2)
# bar.higher.tab <- counts(dds.new2.new,normalized=T)
# head(bar.higher.tab)
# bar.higher.tab.t <- as.data.frame(t(bar.higher.tab))
# head(bar.higher.tab.t)
# head(metadata)
# bar.higher.tab.t.meta <- cbind(metadata[-c(5,10),],bar.higher.tab.t)
# head(bar.higher.tab.t.meta)
# colnames(bar.higher.tab.t.meta)
# library(corrplot)
# N_gene <- as.data.frame((bar.higher.tab.t.meta[,c(6:7,13:17,25:30)]))
# N_gene.cor <- cor((N_gene))
# N_gene.res1 <- cor.mtest((N_gene), conf.level = .95)
# pAdj <- p.adjust(c(N_gene.res1[[1]]), method = "BH")
# resAdj <- matrix(pAdj, ncol = dim(N_gene.res1[[1]])[1])
# quartz()
# corrplot(N_gene.cor,type = "lower",na.label = "NA",p.mat = resAdj, sig.level = .05,insig = "blank",tl.col = "black",tl.srt = 45)
# N_gene.cor.plot






#######using the raw counts of each genes
org.counts <- as.data.frame(counts(dds.new2,normalized=T))
head(org.counts)
org.counts.1 <- mutate(org.counts, gene_name=rownames(org.counts))
org.counts.tab <- left_join(org.counts.1,unique(genes.data.all[,c(15,17)]),by=c("gene_name"="gene_cate"))
genes.data.all.group.new <- data.frame(org.counts.tab[,c(1:11,13)])%>% 
  group_by(group) %>% 
  summarise(across(everything(), sum))
dim(genes.data.all.group.new)
colnames(genes.data.all.group.new)
rownames(genes.data.all.group.new) <- genes.data.all.group.new$group
genes.data.all.group.new.1 <- genes.data.all.group.new[,-1]
rownames(genes.data.all.group.new.1) <- genes.data.all.group.new$group
colnames(genes.data.all.group.new.1)
library(tidyverse)
# bar.higher.tab.t <- as.data.frame(t(genes.data.all.group.new.1))
# head(bar.higher.tab.t)
# head(metadata)
# # bar.higher.tab.t.meta <- cbind(metadata,bar.higher.tab.t)
# # 
# # #bar.higher.tab.t.meta <- cbind(metadata[-c(5,10),],bar.higher.tab.t)
# # head(bar.higher.tab.t.meta)
# # colnames(bar.higher.tab.t.meta)
# # bar.higher.tab.t.meta.long <- pivot_longer(bar.higher.tab.t.meta,
# #                                            cols = colnames(bar.higher.tab.t.meta)[25:30])
# # head(bar.higher.tab.t.meta.long)

# colnames(nitrogen.deseq.cate.long)
# df_cate=bar.higher.tab.t.meta.long %>% group_by(name, Treat) %>% mutate(sem_cate=sd(value)/sqrt(n())) #stdev
# mean_cate = df_cate %>% group_by(name, Treat) %>% mutate(avg_cate=sum(value)/(n()))#combine mean and stdev in "mean_shannon"
# cate.order <- c("Starch","Hemicellulose","Cellulose","Chitin","Pectin","Lignin")#make sure the order is the same to the bac.mag.tree
# mean_cate$name <- factor(mean_cate$name, levels = as.factor(cate.order))

library(beyonce)
library(ggplot2)
######relative abundance calculation step#
#cazy_gene_tab.nor.len <- cazy_gene_tab.nor[,c(4:15)]*cazy_gene_tab.nor$Length
colnames(seq.depth)
seq.depth.trim <- seq.depth[-c(5),]
gene_tab.nor.r <- function(file){
  gene.r.matrix <- matrix(nrow=dim(file)[1],ncol = 0)#if you don't add 0 here, you will get an extra useless columns
  for(a in 1:11){
    newcol <- file[,a]/seq.depth.trim$depth[a]
    gene.r.matrix <- cbind(gene.r.matrix,newcol)
  }
  return(gene.r.matrix)
}

r.nor.cazy.group <- gene_tab.nor.r(genes.data.all.group.new.1)
rownames(r.nor.cazy.group) <- rownames(genes.data.all.group.new.1)
rownames(r.nor.cazy.group) 
head(r.nor.cazy.group)
library(tidyverse)
r.bar.higher.tab.t <- as.data.frame(t(r.nor.cazy.group))
head(r.bar.higher.tab.t)
head(metadata)
#r.bar.higher.tab.t.meta <- cbind(metadata,r.bar.higher.tab.t)
r.bar.higher.tab.t.meta <- cbind(metadata[-c(5),],r.bar.higher.tab.t) %>% 
  mutate(Treat2=c(rep("NCN0",3),rep("NCN60",2),rep("VN0",3),rep("VN60",3)))

#r.bar.higher.tab.t.meta <- cbind(metadata[-c(5,10),],r.bar.higher.tab.t)
head(r.bar.higher.tab.t.meta)
colnames(r.bar.higher.tab.t.meta)
r.bar.higher.tab.t.meta.long <- pivot_longer(r.bar.higher.tab.t.meta,
                                           cols = colnames(r.bar.higher.tab.t.meta)[25:30])
head(r.bar.higher.tab.t.meta.long)
normalTest((r.bar.higher.tab.t.meta$Cellulose))
normalTest((r.bar.higher.tab.t.meta$Chitin))
normalTest((r.bar.higher.tab.t.meta$Lignin))
normalTest((r.bar.higher.tab.t.meta$Starch))
normalTest((r.bar.higher.tab.t.meta$Pectin))
normalTest((r.bar.higher.tab.t.meta$Hemicellulose))
Cellulose.mod <- lmer((Cellulose)~Nitrogen+Cover+(1|block),data =r.bar.higher.tab.t.meta)
Anova(Cellulose.mod)
anova(Cellulose.mod)#Nitrogen
mod1 <-lmer((Cellulose)~Cover+(1|block),data =data_N60)
Anova(mod1)#sig
mod2 <-lmer((Cellulose)~Cover+(1|block),data =data_N0)
Anova(mod2)#sig
post2 <- glht(mod2, linfct=mcp(Cover="Tukey"))
summary(post)
mod3 <-lmer((Cellulose)~Nitrogen+(1|block),data =data_vetch)
Anova(mod3)
mod4 <-lmer((Cellulose)~Nitrogen+(1|block),data =data_nocover)
Anova(mod4)
head(data_N60)
post4 <- glht(mod4, linfct=mcp(Nitrogen="Tukey"))
summary(post4)


Hemicellulose.mod <- lmer((Hemicellulose)~Nitrogen+Cover+(1|block),data =r.bar.higher.tab.t.meta)
Anova(Hemicellulose.mod)#Nitrogen
anova(Hemicellulose.mod)
mod1 <-lmer((Hemicellulose)~Cover+(1|block),data =data_N60)
Anova(mod1)#sig
mod2 <-lmer((Hemicellulose)~Cover+(1|block),data =data_N0)
Anova(mod2)#sig
post2 <- glht(mod2, linfct=mcp(Cover="Tukey"))
summary(post2)
mod3 <-lmer((Hemicellulose)~Nitrogen+(1|block),data =data_vetch)
Anova(mod3)
mod4 <-lmer((Hemicellulose)~Nitrogen+(1|block),data =data_nocover)
Anova(mod4)
head(data_N60)
post4 <- glht(mod4, linfct=mcp(Nitrogen="Tukey"))
summary(post4)


Lignin.mod <- lmer((Lignin)~Nitrogen+Cover+(1|block),data =r.bar.higher.tab.t.meta)
Anova(Lignin.mod)#Nitrogen
anova(Lignin.mod)
mod1 <-lmer((Lignin)~Cover+(1|block),data =data_N60)
Anova(mod1)#sig
mod2 <-lmer((Lignin)~Cover+(1|block),data =data_N0)
Anova(mod2)#sig
post2 <- glht(mod2, linfct=mcp(Cover="Tukey"))
summary(post2)
mod3 <-lmer((Lignin)~Nitrogen+(1|block),data =data_vetch)
Anova(mod3)
mod4 <-lmer((Lignin)~Nitrogen+(1|block),data =data_nocover)
Anova(mod4)
head(data_N60)
post4 <- glht(mod4, linfct=mcp(Nitrogen="Tukey"))
summary(post4)

Pectin.mod <- lmer((Pectin)~Nitrogen+Cover+(1|block),data =r.bar.higher.tab.t.meta)
Anova(Pectin.mod)
anova(Pectin.mod)
mod1 <-lmer((Pectin)~Cover+(1|block),data =data_N60)
Anova(mod1)#sig
mod2 <-lmer((Pectin)~Cover+(1|block),data =data_N0)
Anova(mod1)#sig
post2 <- glht(mod2, linfct=mcp(Cover="Tukey"))
summary(post)
mod3 <-lmer((Pectin)~Nitrogen+(1|block),data =data_vetch)
Anova(mod1)
mod4 <-lmer((Pectin)~Nitrogen+(1|block),data =data_nocover)
Anova(mod4)
head(data_N60)
post4 <- glht(mod4, linfct=mcp(Nitrogen="Tukey"))
summary(post4)
Starch.mod <- lmer((Starch)~Nitrogen+Cover+(1|block),data =r.bar.higher.tab.t.meta)
Anova(Starch.mod)#nothing
anova(Starch.mod)

Chitin.mod <- lmer((Chitin)~Nitrogen*Cover+(1|block),data =r.bar.higher.tab.t.meta)
Anova(Chitin.mod)
anova(Chitin.mod)
data_N60 <- subset(r.bar.higher.tab.t.meta, Nitrogen=="N60")
data_N0 <- subset(r.bar.higher.tab.t.meta, Nitrogen=="N0")
data_vetch <- subset(r.bar.higher.tab.t.meta, Cover=="Vetch")
data_nocover <- subset(r.bar.higher.tab.t.meta, Cover=="No_cover")
library(multcomp)
dim(data_N60)
mod1 <-lmer((Chitin)~Cover+(1|block),data =data_N60)
Anova(mod1)#not sig
mod2 <-lmer((Chitin)~Cover+(1|block),data =data_N0)
Anova(mod2)#not sig
post2 <- glht(mod2, linfct=mcp(Cover="Tukey"))
summary(post)
mod3 <-lmer((Chitin)~Nitrogen+(1|block),data =data_vetch)
Anova(mod3)
mod4 <-lmer((Chitin)~Nitrogen+(1|block),data =data_nocover)
Anova(mod4)
head(data_N60)
post4 <- glht(mod4, linfct=mcp(Nitrogen="Tukey"))
summary(post4)

#cate_gene_bar barplot
colnames(nitrogen.deseq.cate.long)
df_cate.treat=r.bar.higher.tab.t.meta.long %>% group_by(name, Treat) %>% mutate(sem_cate=sd(value*100)/sqrt(n())) #stdev
mean_cate.treat = df_cate.treat %>% group_by(name, Treat) %>% mutate(avg_cate=sum(value*100)/(n()))#combine mean and stdev in "mean_shannon"
cate.order <- c("Starch","Hemicellulose","Cellulose","Chitin","Pectin","Lignin")#make sure the order is the same to the bac.mag.tree
mean_cate.treat$name <- factor(mean_cate.treat$name, levels = as.factor(cate.order))
fig1.all <- subset(mean_cate.treat,name=="Starch"|name=="Hemicellulose"|name=="Cellulose")
plot1 <-
  ggplot(fig1.all, aes(x = name, y = avg_cate, fill = Treat2))+
  geom_col(alpha = 1, position = 'dodge')+
  #  ylim(0,8)+
  scale_fill_manual(values = beyonce_palette(129,4,type = "discrete"),name="Treatment")+
  labs(title="",x="Genes",y="Relative abundance (%)")+
  geom_errorbar(aes(ymin =fig1.all$avg_cate - fig1.all$sem_cate,
                    ymax = fig1.all$avg_cate + fig1.all$sem_cate),width = 0.2, colour = "black", position = position_dodge(.9))+
  #facet_wrap(~Tillage)+
  #theme_bw()+
  theme(legend.position="right",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_line(colour = "#f0f0f0",size=0.75),
        # panel.grid.minor = element_blank(),
        plot.title=(element_text(size=5,family = "Arial",face="plain",hjust = 0.5,vjust = 3)),
        text=element_text(family = "Arial",face="plain"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text=element_text(size=15,family="Arial"),
        axis.title=element_text(size = 15,face="plain",family = "Arial"),
        legend.title = element_text(size=15,face="plain",family = "Arial"),
        legend.text = (element_text(size=15,family = "Arial")))#+
    # annotate(geom="text", x=1.23, y=0.00005, size=8,label="*",color="black", fontface="bold")+
    # annotate(geom="text", x=0.73, y=0.00005, size=8,label="*",color="black", fontface="bold")#+

fig2.all <- subset(mean_cate.treat,name=="Chitin"|name=="Pectin"|name=="Lignin")
plot2 <-
  ggplot(fig2.all, aes(x = name, y = avg_cate, fill = Treat2))+
  geom_col(alpha = 1, position = 'dodge')+
  #  ylim(0,8)+
  scale_fill_manual(values = beyonce_palette(129,4,type = "discrete"),name="Treatment")+
  labs(title="",x="Genes",y="")+
  geom_errorbar(aes(ymin =fig2.all$avg_cate - fig2.all$sem_cate,
                    ymax = fig2.all$avg_cate + fig2.all$sem_cate),width = 0.2, colour = "black", position = position_dodge(.9))+
  #facet_wrap(~Tillage)+
  #theme_bw()+
  theme(legend.position="right",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_line(colour = "#f0f0f0",size=0.75),
        # panel.grid.minor = element_blank(),
        plot.title=(element_text(size=5,family = "Arial",face="plain",hjust = 0.5,vjust = 3)),
        text=element_text(family = "Arial",face="plain"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text=element_text(size=15,family="Arial"),
        axis.title=element_text(size = 15,face="plain",family = "Arial"),
        legend.title = element_text(size=15,face="plain",family = "Arial"),
        legend.text = (element_text(size=15,family = "Arial")))
# #   annotate(geom="text", x=1.23, y=0.00005, size=8,label="*",color="black", fontface="bold")+
# #   annotate(geom="text", x=0.73, y=0.00005, size=8,label="*",color="black", fontface="bold")#+
# #   annotate(geom="text", x=2, y=0.00027, size=8,label="*",color="black", fontface="bold")+
# #   annotate(geom="text", x=4.11, y=0.00024, size=8,label="*",color="black", fontface="bold")+
# #   annotate(geom="text", x=4.23, y=0.00023, size=8,label="*",color="black", fontface="bold")+
# #   annotate(geom="text", x=5.11, y=0.00019, size=8,label="*",color="black", fontface="bold")+
# #   annotate(geom="text", x=5.23, y=0.00018, size=8,label="*",color="black", fontface="bold")+
# #   annotate(geom="segment", x=5.13, xend=5.33, y=0.00018,yend = 0.00018, size=0.5)+
# #   annotate(geom="segment", x=4.89, xend=5.33, y=0.00019,yend = 0.00019, size=0.5)+
# #   annotate(geom="segment", x=4.13, xend=4.33, y=0.00023,yend = 0.00023, size=0.5)+
# #   annotate(geom="segment", x=3.89, xend=4.33, y=0.00024,yend = 0.00024, size=0.5)+
# #   annotate(geom="segment", x=3.89, xend=4.33, y=0.00024,yend = 0.00024, size=0.5)+
# #   annotate(geom="segment", x=1.13, xend=1.33, y=0.00005,yend = 0.00005, size=0.5)+
# #   annotate(geom="segment", x=0.63, xend=0.83, y=0.00005,yend = 0.00005, size=0.5)
tiff('figures/highergene_treatbar.tiff', units="in", width=16, height=8, res=300)
ggarrange(plot1,plot2,
          labels = c("A","B"),ncol = 2, nrow = 1,common.legend =T ,legend = "bottom", widths=c(1,1),heights = c(1,1))
dev.off()
# colnames(nitrogen.deseq.cate.long)
# df_cate=r.bar.higher.tab.t.meta.long %>% group_by(name, Cover) %>% mutate(sem_cate=sd(value*100)/sqrt(n())) #stdev
# mean_cate = df_cate %>% group_by(name, Cover) %>% mutate(avg_cate=sum(value*100)/(n()))#combine mean and stdev in "mean_shannon"
# cate.order <- c("Starch","Hemicellulose","Cellulose","Chitin","Pectin","Lignin")#make sure the order is the same to the bac.mag.tree
# mean_cate$name <- factor(mean_cate$name, levels = as.factor(cate.order))
# library(beyonce)
# library(ggplot2)
# #cate_gene_bar barplot
# fig1 <- subset(mean_cate,name=="Starch"|name=="Hemicellulose"|name=="Cellulose")
# fig2 <- subset(mean_cate,name=="Chitin"|name=="Pectin"|name=="Lignin")
# 
# plot1 <-
#   ggplot(fig1, aes(x = name, y = avg_cate, fill = Cover))+
#   geom_col(alpha = 1, position = 'dodge')+
#   #  ylim(0,8)+
#   scale_fill_manual(values = c("#99b872","#4c6e21"),name="Treatment")+
#   labs(title="",x="Genes",y="Relative abundance (%)")+
#   geom_errorbar(aes(ymin =fig1$avg_cate - fig1$sem_cate, ymax = fig1$avg_cate + fig1$sem_cate),width = 0.2, colour = "black", position = position_dodge(.9))+
#   #facet_wrap(~Tillage)+
#   #theme_bw()+
#   theme(legend.position="right",panel.background = element_rect(fill = "white",color="black"),
#         panel.grid.major = element_line(colour = "#f0f0f0",size=0.75),
#         # panel.grid.minor = element_blank(),
#         plot.title=(element_text(size=5,family = "Arial",face="plain",hjust = 0.5,vjust = 3)),
#         text=element_text(family = "Arial",face="plain"),
#         panel.border = element_rect(colour = "black", fill=NA, size=0.5),
#         axis.text=element_text(size=15,family="Arial"),
#         axis.title=element_text(size = 15,face="plain",family = "Arial"),
#         legend.title = element_text(size=15,face="plain",family = "Arial"),
#         legend.text = (element_text(size=15,family = "Arial")))
# plot2 <-
#   ggplot(fig2, aes(x = name, y = avg_cate, fill = Cover))+
#   geom_col(alpha = 1, position = 'dodge')+
#   #  ylim(0,8)+
#   scale_fill_manual(values = c("#99b872","#4c6e21"),name="Treatment")+
#   labs(title="",x="Genes",y="")+
#   geom_errorbar(aes(ymin =fig2$avg_cate - fig2$sem_cate, ymax = fig2$avg_cate + fig2$sem_cate),width = 0.2, colour = "black", position = position_dodge(.9))+
#   #facet_wrap(~Tillage)+
#   #theme_bw()+
#   theme(legend.position="right",panel.background = element_rect(fill = "white",color="black"),
#         panel.grid.major = element_line(colour = "#f0f0f0",size=0.75),
#         # panel.grid.minor = element_blank(),
#         plot.title=(element_text(size=5,family = "Arial",face="plain",hjust = 0.5,vjust = 3)),
#         text=element_text(family = "Arial",face="plain"),
#         panel.border = element_rect(colour = "black", fill=NA, size=0.5),
#         axis.text=element_text(size=15,family="Arial"),
#         axis.title=element_text(size = 15,face="plain",family = "Arial"),
#         legend.title = element_text(size=15,face="plain",family = "Arial"),
#         legend.text = (element_text(size=15,family = "Arial")))
# # +
# #     annotate(geom="text", x=1, y=0.00045, size=12,label="**",color="red", fontface="bold")+
# #     annotate(geom="segment", x=0.8, xend=1.2, y=0.00045,yend = 0.00045, size=0.5)+
# #     annotate(geom="text", x=3, y=0.0006, size=12,label="**",color="red", fontface="bold")+
# #     annotate(geom="segment", x=2.8, xend=3.2, y=0.0006,yend = 0.0006, size=0.5)
# 
# 
# 
# colnames(nitrogen.deseq.cate.long)
# df_cate.2=r.bar.higher.tab.t.meta.long %>% group_by(name, Nitrogen) %>% mutate(sem_cate=sd(value*100)/sqrt(n())) #stdev
# mean_cate.2 = df_cate.2%>% group_by(name, Nitrogen) %>% mutate(avg_cate=sum(value*100)/(n()))#combine mean and stdev in "mean_shannon"
# cate.order.2 <- c("Starch","Hemicellulose","Cellulose","Chitin","Pectin","Lignin")#make sure the order is the same to the bac.mag.tree
# mean_cate.2$name <- factor(mean_cate$name, levels = as.factor(cate.order))
# 
# 
# fig1.2 <- subset(mean_cate.2,name=="Starch"|name=="Hemicellulose"|name=="Cellulose")
# fig2.2 <- subset(mean_cate.2,name=="Chitin"|name=="Pectin"|name=="Lignin")
# 
# 
# plot1_N <-
# ggplot(fig1.2, aes(x = name, y = avg_cate, fill = Nitrogen))+
#     geom_col(alpha = 1, position = 'dodge')+
#     #  ylim(0,8)+
#     scale_fill_manual(values = c("#b8a272","#a1720b"),name="Treatment")+
#     labs(title="",x="Genes",y="Relative abundance (%)")+
#     geom_errorbar(aes(ymin =fig1.2$avg_cate - fig1.2$sem_cate, ymax = fig1.2$avg_cate + fig1.2$sem_cate),width = 0.2, colour = "black", position = position_dodge(.9))+
#     #facet_wrap(~Tillage)+
#     #theme_bw()+
#     theme(legend.position="right",panel.background = element_rect(fill = "white",color="black"),
#           panel.grid.major = element_line(colour = "#f0f0f0",size=0.75),
#           # panel.grid.minor = element_blank(),
#           plot.title=(element_text(size=5,family = "Arial",face="plain",hjust = 0.5,vjust = 3)),
#           text=element_text(family = "Arial",face="plain"),
#           panel.border = element_rect(colour = "black", fill=NA, size=0.5),
#           axis.text=element_text(size=15,family="Arial"),
#           axis.title=element_text(size = 15,face="plain",family = "Arial"),
#           legend.title = element_text(size=15,face="plain",family = "Arial"),
#           legend.text = (element_text(size=15,family = "Arial")))+
#   annotate(geom="text", x=3, y=0.01, size=12,label="*",color="red", fontface="bold")+
#   annotate(geom="segment", x=2.8, xend=3.2, y=0.01,yend = 0.01, size=0.5)+#cellulose
#   annotate(geom="text", x=2, y=0.0075, size=12,label="*",color="red", fontface="bold")+
#   annotate(geom="segment", x=1.8, xend=2.2, y=0.0075,yend = 0.0075, size=0.5)#hemicellulose
# #hemicellulose
# plot2_N <-
# ggplot(fig2.2, aes(x = name, y = avg_cate, fill = Nitrogen))+
#     geom_col(alpha = 1, position = 'dodge')+
#     #  ylim(0,8)+
#     scale_fill_manual(values = c("#b8a272","#a1720b"),name="Treatment")+
#     labs(title="",x="Genes",y="")+
#     geom_errorbar(aes(ymin =fig2.2$avg_cate - fig2.2$sem_cate, ymax = fig2.2$avg_cate + fig2.2$sem_cate),width = 0.2, colour = "black", position = position_dodge(.9))+
#     #facet_wrap(~Tillage)+
#     #theme_bw()+
#     theme(legend.position="right",panel.background = element_rect(fill = "white",color="black"),
#           panel.grid.major = element_line(colour = "#f0f0f0",size=0.75),
#           # panel.grid.minor = element_blank(),
#           plot.title=(element_text(size=5,family = "Arial",face="plain",hjust = 0.5,vjust = 3)),
#           text=element_text(family = "Arial",face="plain"),
#           panel.border = element_rect(colour = "black", fill=NA, size=0.5),
#           axis.text=element_text(size=15,family="Arial"),
#           axis.title=element_text(size = 15,face="plain",family = "Arial"),
#           legend.title = element_text(size=15,face="plain",family = "Arial"),
#           legend.text = (element_text(size=15,family = "Arial")))+
#   annotate(geom="text", x=3, y=0.0007, size=12,label="*",color="red", fontface="bold")+
#   annotate(geom="segment", x=2.8, xend=3.2, y=0.0007,yend = 0.0007, size=0.5)
# library(ggpubr)
# tiff('figures/highergene_cate.tiff', units="in", width=12, height=14, res=300)
# ggarrange(plot1_N,plot2_N,plot1,plot2,
#           labels = c("A","B","C","D"),ncol = 2, nrow = 2,common.legend = F,legend = "bottom", widths=c(1,1),heights = c(1,1))
# dev.off()

# ####higher level rda####3
# library(phyloseq)
# library(vegan)
# library(ggforce)
# library(MASS)
# library(car)
# library(lme4)
# library(tidyverse)
# library(fBasics)
#rld2.new <- rlog(dds.new2.new)
gene.phyloseq.new <- merge_phyloseq(sample_data(metadata), otu_table(as.data.frame(r.nor.cazy.group),taxa_are_rows = TRUE))
gene.phyloseq.new
formularda <- rda(t(otu_table(gene.phyloseq)) ~ Cover+Nitrogen+Cover*Nitrogen+pH+GWC+NH4+NO3+
                    AG+BG+CB+LAP+NAG+PHOS+XYL,data = metadata.up)
summary(formularda)
rda_data <- ordinate(gene.phyloseq,method = "RDA",formula = ~Cover+Nitrogen+Cover*Nitrogen+pH+GWC+NH4+NO3)
env.rda <- rda_data$CCA
arrow.dat <- as.data.frame(env.rda$biplot)
set.seed(2021)
step.backward <-
  ordistep(rda_data,
           permutations = how(nperm = 999)
  )
Anova(rda_data)
anova(rda_data)
anova.perm <- anova(rda_data,permutations = how(nperm = 9999), by="margin");anova.perm


adj <- RsquareAdj(rda_data);adj
adonis(t(otu_table(gene.phyloseq))~Nitrogen*Cover+pH+GWC+NH4+NO3,
       data=data.frame(sample_data(gene.phyloseq)),
       permutations=9999,by="margin")#only nitrogen significant

df_pca <- prcomp(t(otu_table(gene.phyloseq)))
pca_tab <- as.data.frame(df_pca$x)
pca_tab$Nitrogen <- metadata$Nitrogen
pca_tab$Cover <- metadata$Cover

higherrda.plot <-
  plot_ordination(gene.phyloseq,rda_data,type="sample",color="Nitrogen",shape = "Cover")+
  geom_point(size=3)+
  geom_text(aes(label=rownames(pca_tab),vjust=-2,hjust=0.5),size=2)+
  ggtitle("C genes")+
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
        legend.text = (element_text(size=15,family = "Arial")))

library(ggpubr)
#divergent <- ggarrange(res2.cover.plot,res2.nocover.plot,res2_n60.sig.plot,
#                      labels = c("A","","","B"),ncol = 2, nrow = 2,common.legend = ,legend = "bottom", widths=c(1,1))

tiff('figures/highergene_cate_cazy_divergent.tiff', units="in", width=18, height=18, res=300)
ggarrange(higherrda.plot,
          labels = c("A","B","C","D"),ncol = 2, nrow = 2,common.legend =F ,legend = "bottom", widths=c(1,1),heights = c(1,1))
dev.off()

# #####correlation with environmental factor
# library(DESeq2)
# library(corrplot)
# rld.data <- as.data.frame(r.nor.cazy.group)
# rld.data.t <- as.data.frame(t(rld.data))
# colnames(rld.data.t)
# #view(rld.data.t)
# cor.rld.data <- cbind(metadata,rld.data.t)
# colnames(cor.rld.data)
# cor.rld.data.n60 <-subset(cor.rld.data, Nitrogen=="N0") 

# N_gene <- as.data.frame((cor.rld.data.n60[,c(6:7,13:15,25:30)]))
# N_gene.cor <- cor((N_gene))
# N_gene.res1 <- cor.mtest((N_gene), conf.level = .95)
# pAdj <- p.adjust(c(N_gene.res1[[1]]), method = "BH")
# resAdj <- matrix(pAdj, ncol = dim(N_gene.res1[[1]])[1])
# quartz()
# corrplot(N_gene.cor,type = "lower",na.label = "NA",p.mat = resAdj, sig.level = .05,insig = "blank",tl.col = "black",tl.srt = 45)
# N_gene.cor.plot

