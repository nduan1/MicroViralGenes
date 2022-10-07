library(tidyverse)
library(phyloseq)
library(vegan)
library(ggplot2)
library(tidyr)
library(readxl)
library(base)
#gene mapping data
seq.depth <- read.csv("dataset/seq.depth.csv", header = T, row.names = 1);head(seq.depth)
nitrogen_gene <- read.csv("dataset/filtered_counts_coverage_table.csv", header = T, sep = "|")
dim(nitrogen_gene)
head(nitrogen_gene)
nitrogen_gene_tab <- separate(nitrogen_gene[,c(2:15)],col = 1,sep = "rank:",into = c("contig","gene_name"))
head(nitrogen_gene_tab)
colnames(nitrogen_gene_tab)
colnames(nitrogen_gene_tab)[4:15] <- c("NCNTN0_4","VNTN0_4","VNTN60_2","NCNTN0_1","VNTN0_1",
                                   "VNTN60_1","NCNTN60_1","VNTN0_2","NCNTN0_2","NCNTN60_4","VNTN60_4","NCNTN60_2")
colnames(nitrogen_gene_tab)
seq.depth$sample
nitrogen_gene_tab.new <- nitrogen_gene_tab[,c(1:3,7,12,4,10,15,13,8,11,5,9,6,14)]
colnames(nitrogen_gene_tab.new)
head(nitrogen_gene_tab.new)
#identical(colnames(nitrogen_gene_tab.new)[4:15],seq.depth$sample)
#normalization step
votu.nor <- function(file){
  nor.matrix <- matrix(nrow=dim(file)[1],ncol = 0)#if you don't add 0 here, you will get an extra useless columns
  for(a in 1:12){
    newcol <- round(file[,a]*(seq.depth$mutiplier[a]),3)
    nor.matrix <- cbind(nor.matrix,newcol)
  }
  return(nor.matrix)
}

nitrogen_gene_tab.nor <- cbind(nitrogen_gene_tab.new[,1:3],votu.nor(nitrogen_gene_tab.new[,c(4:15)]))
colnames(nitrogen_gene_tab.nor) <- colnames(nitrogen_gene_tab.new)
rownames(nitrogen_gene_tab.nor) <- rownames(nitrogen_gene_tab.new)
head(nitrogen_gene_tab.nor)
dim(nitrogen_gene_tab.nor)#[1] 363  15
colnames(nitrogen_gene_tab.nor)
#######relative abundance calculation step#
# nitrogen_gene_tab.nor.len <- nitrogen_gene_tab.nor[,c(4:15)]*nitrogen_gene_tab.nor$Length
# gene_tab.nor.r <- function(file){
#   gene.r.matrix <- matrix(nrow=dim(file)[1],ncol = 0)#if you don't add 0 here, you will get an extra useless columns
#   for(a in 1:12){
#     newcol <- file[,a]/seq.depth$depth[a]
#     gene.r.matrix <- cbind(gene.r.matrix,newcol)
#   }
#   return(gene.r.matrix)
# }
# nitrogen_gene_tab.nor.r <- cbind(nitrogen_gene_tab.nor[,1:3],gene_tab.nor.r(nitrogen_gene_tab.nor.len))
# colnames(nitrogen_gene_tab.nor.r) <- colnames(nitrogen_gene_tab.nor)
# 
write.csv(nitrogen_gene_tab.nor,"dataset/nitrogen_gene_tab.nor.csv")###manually organize files
# write.csv(nitrogen_gene_tab.nor.r,"dataset/nitrogen_gene_tab.nor.r.csv")
###manually organize files and I don't need the reference file to categorize it
#use vim to remove \t\n
nitrogen_man <- read.csv("dataset/nitrogen_gene_man.nor.csv", row.names = 1)#
head(nitrogen_man);dim(nitrogen_man)
###relative abundance calculation###
nitrogen_man_sum <- nitrogen_man[,c(5,8:19)]%>% group_by(comp_gene) %>% summarise(across(everything(),sum))
head(nitrogen_man_sum)
nitrogen_man_sum.r <- apply(nitrogen_man_sum[,c(2:13)], 2, function(x) x/sum(x))
head(nitrogen_man_sum.r)
#write.csv(nitrogen_man_sum,"nitrogen_man_sum.csv")
nitrogen_man_sum_long <- pivot_longer(cbind(nitrogen_man_sum[1],nitrogen_man_sum.r),cols = colnames(nitrogen_man_sum)[2:13])
head(nitrogen_man_sum_long)
colnames(nitrogen_man_sum_long)[2:3] <- c("sample_ID","relative_abundance")
head(nitrogen_man_sum_long)
meta_data <- read.csv("dataset/meta_data.csv")
head(meta_data)
colnames(meta_data)[1] <-"sample_ID"
head(meta_data)
nitrogen_man_tab <- left_join(nitrogen_man_sum_long,meta_data,by=c("sample_ID"="sample_ID"))
head(nitrogen_man_tab)
nitrogen_tab <- inner_join(nitrogen_man[,3:6],nitrogen_man_tab,by=c("comp_gene"="comp_gene"))
head(nitrogen_tab)#for plot
#######statistical analysis####
library(car)
library(fBasics)
library(MASS)
library(lme4)
library(multcomp)
nitrogen_stats_tab <- (cbind(nitrogen_man_sum[1],nitrogen_man_sum.r))
rownames(nitrogen_stats_tab) <- nitrogen_stats_tab$comp_gene
head(nitrogen_stats_tab)
nitrogen_stats <- apply(t(nitrogen_stats_tab)[-c(1),],2,as.numeric)
nitrogen_stats_tb <- mutate(as.data.frame(nitrogen_stats),sample_ID =colnames(nitrogen_stats_tab)[2:13])
dim(nitrogen_stats_tb)
nitrogen_stats_dat <- left_join(nitrogen_stats_tb,meta_data,by=c("sample_ID"="sample_ID"))
normalTest(nitrogen_stats_dat$hao)
normalTest(nitrogen_stats_dat$napA)
normalTest(nitrogen_stats_dat$napB)
normalTest(nitrogen_stats_dat$narB)#not good
normalTest(nitrogen_stats_dat$`narG,narZ,nxrA`)
normalTest(nitrogen_stats_dat$`narH,narY,nxrB`)
normalTest(nitrogen_stats_dat$`narI,narV`)
normalTest(nitrogen_stats_dat$`narJ,narW`)#not good
normalTest(sqrt(nitrogen_stats_dat$nasA))#not good
normalTest((nitrogen_stats_dat$nifK))#not good
normalTest(nitrogen_stats_dat$nirA)
normalTest(nitrogen_stats_dat$nirB)#not good
normalTest(nitrogen_stats_dat$nirK)
normalTest(nitrogen_stats_dat$nirS)
normalTest(nitrogen_stats_dat$norB)
normalTest(nitrogen_stats_dat$norC)#not good
normalTest(nitrogen_stats_dat$NosZ)
normalTest(nitrogen_stats_dat$nrfA)
normalTest(nitrogen_stats_dat$nrfH)
normalTest(nitrogen_stats_dat$`NRT,narK,nrtP,nasA`)
normalTest(nitrogen_stats_dat$`pmoA-amoA`)
normalTest(nitrogen_stats_dat$`pmoB-amoB`)
normalTest(nitrogen_stats_dat$`pmoC-amoC`)
hao.mod <- lmer(hao~Nitrogen*Cover+(1|block),data = nitrogen_stats_dat)
Anova(hao.mod)#nitrogen
anova(hao.mod)
napA.mod <- lmer(napA~Nitrogen*Cover+(1|block),data = nitrogen_stats_dat)
Anova(napA.mod)#nitrogen
anova(napA.mod)
napB.mod <- lmer(napB~Nitrogen*Cover+(1|block),data = nitrogen_stats_dat)
Anova(napB.mod)#nitrogen and cover
anova(napB.mod)
hao.mod <- lmer(hao~Nitrogen*Cover+(1|block),data = nitrogen_stats_dat)
Anova(hao.mod)#nitrogen
anova(hao.mod)



########bar plot with error bar###
df_nitrogen=nitrogen_tab %>% group_by(Treat, comp_gene) %>% mutate(sem_nitrogen=sd(relative_abundance)/sqrt(n())) #stdev
mean_nitrogen = df_nitrogen %>% group_by(Treat, comp_gene) %>% mutate(avg_nitrogen=sum(relative_abundance)/(n()))#combine mean and stdev in "mean_shannon"
#group.order <- c("Starch","Hemicellulose","Cellulose","Chitin","Pectin","Aromatic","Lignin")#make sure the order is the same to the bac.mag.tree
#mean_carbon$group <- factor(mean_carbon$group, levels = as.factor(group.order))
#pectin_data <- subset(mean_carbon_rev,group=="Pectin")
nitrogen_bar <- ggplot(mean_nitrogen, aes(x = comp_gene, y = avg_nitrogen, fill = Treat))+
  geom_col(alpha = 0.5,width = 0.7, position = 'dodge')+
  #  ylim(0,8)+
  scale_fill_manual(values = c("#8a0101","#254870","#a0a3a2","#7a4c0f"),name="Treatment")+
  labs(title="Pectin",x="Treatment",y="")+
  geom_errorbar(aes(ymin =mean_nitrogen$avg_nitrogen - mean_nitrogen$sem_nitrogen, 
                    ymax = mean_nitrogen$avg_nitrogen + mean_nitrogen$sem_nitrogen),
                width = 0.2, colour = "black", position = position_dodge(.9))+
  #facet_wrap(~Tillage)+
  #  theme_bw()+
  theme(legend.position="right",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title=(element_text(size=15,family = "Arial",face="plain",hjust = 0.5,vjust = 3)),
        text=element_text(family = "Arial",face="plain"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text=element_text(size=15,family="Arial"),
        axis.title=element_text(size = 15,face="plain",family = "Arial"),
        legend.title = element_text(size=15,face="plain",family = "Arial"),
        legend.text = (element_text(size=15,family = "Arial")))+
  annotate(geom="text", x=1, y=0.0025, size=8,label="*",color="black", fontface="bold")


