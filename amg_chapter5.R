library(tidyverse)
library(phyloseq)
library(vegan)
library(ggplot2)
library(tidyr)
library(base)
library(tidyverse)
library(ggpubr)
seq.depth <- read.csv("dataset/seq.depth.csv", header = T, row.names = 1);head(seq.depth)
amg_gene <- read.csv("dataset/counts_coverage_table.tsv", header = T, sep = "\t")
dim(amg_gene)
head(amg_gene)
colnames(amg_gene)
amg_gene_tab <- separate(amg_gene,col = 1,sep = "#",into = c("contig","others"))
head(amg_gene_tab)
amg_gene_dat <- amg_gene_tab[,-c(2:3)]
head(amg_gene_dat)
colnames(amg_gene_dat)
colnames(amg_gene_dat)[4:15] <- c("NCN0_4","VN0_4","VN60_2","NCN0_1","VN0_1",
                                       "VN60_1","NCN60_1","VN0_2","NCN0_2","NCN60_4","VN60_4","NCN60_2")
colnames(amg_gene_dat)
amg_gene_tab.new <- amg_gene_dat[,c(1:3,7,12,4,10,15,13,8,11,5,9,6,14)]
colnames(amg_gene_tab.new)
head(amg_gene_tab.new)

meta_data <- read.csv("dataset/meta_data.csv")
colnames(meta_data)

meta_data.1 <- read.csv("dataset/meta_data.csv", row.names = 1)
colnames(meta_data.1)
###microbes###
microbe_gene <- read.csv("dataset/microbes_filtered_counts_coverage_table.csv",sep = "|", row.names = 1)
head(microbe_gene)
microbe_gene_tab.1 <- separate(microbe_gene,col = 1,sep = "rank:*",into = c("contig","genes"))
microbe_gene_tab.2 <- separate(microbe_gene_tab.1,col = 2,sep = ";",into = c("rank","genes"))
microbe_gene_tab <- microbe_gene_tab.2[,-c(2,4)]
head(microbe_gene_tab)
colnames(microbe_gene_tab)
colnames(microbe_gene_tab)[3:14] <- c("NCN0_4","VN0_4","VN60_2","NCN0_1","VN0_1",
                                  "VN60_1","NCN60_1","VN0_2","NCN0_2","NCN60_4","VN60_4","NCN60_2")
colnames(microbe_gene_tab)
microbe_gene_tab.new<- microbe_gene_tab[,c(1:2,6,11,3,9,14,12,7,10,4,8,5,13)]
colnames(microbe_gene_tab.new)
dim(microbe_gene_tab.new)
microbe_gene_table.1 <- microbe_gene_tab.new %>% 
  filter(!grepl("malate",genes)&!grepl("protoporphyrinogen",genes)) 
microbe_gene_table.1$genes <- gsub(" (db=kegg)","",microbe_gene_table.1$genes,fixed = T)
unique(microbe_gene_table.1$genes)
colnames(microbe_gene_table.1)
microbe_gene_table <- microbe_gene_table.1[rowSums(microbe_gene_table.1==0)<12,]
dim(microbe_gene_table)
dim(microbe_gene_table)#[1] 340
UDP_glucose <- subset(microbe_gene_table, genes==" UDP-glucose 4-epimerase [EC:5.1.3.2]")
dim(UDP_glucose)#[1] 206  14
#View(UDP_glucose)
write.csv(UDP_glucose$contig,"dataset/microbe_UDP_glucose.txt")

UDP_galactopyranose <- subset(microbe_gene_table, genes==" UDP-galactopyranose mutase [EC:5.4.99.9]")
dim(UDP_galactopyranose)#[1] 54 14
write.csv(UDP_galactopyranose$contig,"dataset/microbe_UDP_galactopyranose.txt")
GDPmannose <- subset(microbe_gene_table, genes==" GDPmannose 4,6-dehydratase [EC:4.2.1.47]")
dim(GDPmannose)#[1] 94 14
write.csv(GDPmannose$contig,"dataset/microbe_GDPmannose.txt")
glycogen_synthase <- subset(microbe_gene_table, genes==" glycogen synthase [EC:2.4.1.11]")
dim(glycogen_synthase)#[1] 5 14
write.csv(glycogen_synthase$contig,"dataset/microbe_glycogen_synthase.txt")



###amg stacked bar###
colnames(amg_gene_tab.new)
rownames(amg_gene_tab.new)
amg_gene_tab.stack <- amg_gene_tab.new[,c(3:15)] %>% 
  group_by(gene) %>% 
  summarise(across(everything(),sum))
#view(amg_gene_tab.stack)
write.csv(amg_gene_tab.stack,"amg_gene_tab.stack.csv")
gene_tab.nor.r <- function(file){
  gene.r.matrix <- matrix(nrow=dim(file)[1],ncol = 0)#if you don't add 0 here, you will get an extra useless columns
  for(a in 1:12){
    newcol <- file[,a]/seq.depth$depth[a]
    gene.r.matrix <- cbind(gene.r.matrix,newcol)
  }
  return(gene.r.matrix)
}
amg_gene_tab.stack.r <- cbind(amg_gene_tab.stack[,1],gene_tab.nor.r(amg_gene_tab.stack[2:13]))
colnames(amg_gene_tab.stack.r) 

#amg_gene_r <-right_join(amg_gene_tab.new[2:3],amg_gene_tab.stack.r)
amg_gene_r_long <- pivot_longer(as.data.frame(amg_gene_tab.stack.r),cols=colnames(amg_gene_tab.stack.r)[2:13])
amg_gene_short.stack <- right_join(unique(amg_gene_tab.new[2:3]),amg_gene_r_long,by=c("gene"="gene"))

amg_gene_short.stack.dat <- right_join(meta_data,amg_gene_short.stack,by=c("X"="name.y"))
colnames(amg_gene_short.stack.dat)
df_amg=amg_gene_short.stack.dat %>% group_by(Treat, name.x) %>% mutate(sum_treat=sum(value)) #stdev
head(df_amg);dim(df_amg)
cols_phylum <- c("#a195ab", "#9B110E", "#688787" ,"#d49783", "#550307", "#61a1c9", "#FDD262", "#D3DDDC", "#C7B19C","#899DA4", "#FAEFD1", "#76aac4","#798E87", "#C27D38","#CCC591","#29211F")
colnames(df_amg)
df_amg.unique <- unique(df_amg[,c(3,27,29)])
dim(df_amg.unique)
library(beyonce)
unique((df_amg.unique)$gene)
dim(df_amg.unique)
df_amg.unique.new <- subset(df_amg.unique,gene!="neuA"&
                              gene!="RGP_UTM"&
                              gene!="RHM")
dim(df_amg.unique.new)
tiff('figures/stacked_amg.treat.tiff', units="in", width=10, height=5, res=300)
#stacked_bar <- 
  ggplot(df_amg.unique.new, aes(x =X, y = value*100, fill=name.x))+
  geom_bar(stat='identity',colour="grey",size=0.1, width = 0.6,alpha=1)+
  labs(title="",x="Treatment",y="Relative abundance (%)")+
  scale_fill_manual(name="AMGs",values=cols_phylum)+
    guides(fill=guide_legend(nrow=7,ncol=1,bycol=TRUE))+
  theme(legend.position="right",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_line(colour="grey",size = 0.1),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=10,family = "Arial",face="plain",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="plain",size = 15),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(size=10,family="Arial",face = "plain"),
        axis.title=element_text(size =10,face="plain",family = "Arial",vjust = 1),
        legend.title = element_text(size=10,face="plain",family = "Arial"),
        legend.text = (element_text(size=10,family = "Arial")))
dev.off()

####heatmap
colnames(df_amg.unique.new)
dim(df_amg.unique.new)#84 28
amg_hp <- df_amg.unique.new#subset(amg_gene_short.stack.dat,name.x!="N-acylneuraminate cytidylyltransferase [EC:2.7.7.43]"&
                   #name.x!="reversibly glycosylated polypeptide / UDP-arabinopyranose mutase [EC:2.4.1.- 5.4.99.30]"&
                   #name.x!="UDP-glucose 4,6-dehydratase [EC:4.2.1.76]")
dim(amg_hp)
colnames(amg_hp)
glimpse(amg_hp)
unique(amg_hp$X)
hp_dat <- amg_hp#[,c(1,27:28)] 
glimpse(hp_dat)
hp_dat$sum_treat[hp_dat$sum_treat==0] <- NA
head(hp_dat)
tiff('figures/hp_amg.treat.tiff', units="in", width=8, height=5, res=300)
ggplot(hp_dat,aes(x = Treat, y = gene, fill = (sum_treat)*100))+
  guides(fill=guide_legend(title="Relative abundnace (%)"))+
  geom_tile(colour="white",stat = "identity",size=0.25) +
  labs(x="",y="",title="AMGs")+
  scale_fill_distiller(palette = 'Greens',direction = -1,na.value =NA)+
  theme(legend.position="right",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_line(colour="grey",size = 0.1),panel.grid.minor = element_line(colour="grey",size = 0.1),
        plot.title=(element_text(size=10,family = "Arial",face="plain",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="plain",size = 15),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(size=10,family="Arial",face = "plain"),
        axis.title=element_text(size =10,face="plain",family = "Arial",vjust = 1),
        legend.title = element_text(size=10,face="plain",family = "Arial"),
        legend.text = (element_text(size=10,family = "Arial")))
dev.off()
#export figure
#ggsave(p,filename="measles-mod3.png",height=5.5,width=8.8,units="in",dpi=300)
  

dev.off()
####amg deseq####
head(amg_gene_tab.stack)
rownames(amg_gene_tab.stack) <- (amg_gene_tab.stack)$contig
rownames(amg_gene_tab.stack)
amg.new.dat2 <- amg_gene_tab.stack[,-1]
rownames(amg.new.dat2) <- (amg_gene_tab.stack)$contig
rownames(amg.new.dat2)
colnames(amg.new.dat2)
colnames(meta_data.1)
amg.dds.new2 <- DESeqDataSetFromMatrix(countData = amg.new.dat2,
                                   colData = meta_data.1,
                                   design = ~ Treat)
nrow(amg.dds.new2)
amg.dds.new2 <- DESeq(amg.dds.new2)#doesn't work
nrow(amg.dds.new2)
sizeFactors(amg.dds.new2)
plotDispEsts(amg.dds.new2)
resultsNames(amg.dds.new2)


#####only GDPmannose 4,6-dehydratase [EC:4.2.1.47] in virus have enough data to do correlation and also can be found in the microbial contigs###
GDPmannose <- subset(microbe_gene_table, genes==" GDPmannose 4,6-dehydratase [EC:4.2.1.47]")
dim(GDPmannose)
head(GDPmannose)
GDPmannose.tab <- GDPmannose[,-2]
head(GDPmannose.tab)
head(amg_gene_tab.stack)
colnames(amg_gene_tab.stack)[1]<-"contig" 
GDPmannose.amg.microbe <- rbind(GDPmannose.tab,amg_gene_tab.stack[4,])
head(GDPmannose.amg.microbe)
dim(GDPmannose.amg.microbe)
####deseq2 normalization
library(DESeq2)
rownames(GDPmannose.amg.microbe) <- (GDPmannose.amg.microbe)$contig
rownames(GDPmannose.amg.microbe)
new.dat2 <- rbind(GDPmannose.amg.microbe[,-1],amg.new.dat2[4,])[-96,]
rownames(new.dat2)
colnames(new.dat2)
colnames(meta_data)
dds.new2 <- DESeqDataSetFromMatrix(countData = new.dat2,
                                   colData = meta_data.1,
                                   design = ~ Nitrogen*Cover)
nrow(dds.new2)
dds.new2 <- DESeq(dds.new2)
nrow(dds.new2)
sizeFactors(dds.new2)
plotDispEsts(dds.new2)
resultsNames(dds.new2)
GDPmannose.amg.microbe <- as.data.frame(counts(dds.new2,normalized=T))
dim(GDPmannose.amg.microbe)
GDPmannose.amg.microbe.t <- as.data.frame(t(GDPmannose.amg.microbe))[,c(14,40,41,51,53,95)]

microbe <- apply(GDPmannose.amg.microbe[1:94,], 2,sum)
amg <- GDPmannose.amg.microbe[95,]
microbesum.amg <- rbind(microbe,amg)
microbesum.amg.t <- as.data.frame(t(microbesum.amg))
#cor.test(microbe,amg)
library(corrplot)
N_gene.cor <- cor(microbesum.amg.t)#[,c(21:50,95)])
N_gene.res1 <- cor.mtest((GDPmannose.amg.microbe.t), conf.level = .95)
pAdj <- p.adjust(c(N_gene.res1[[1]]), method = "BH")
resAdj <- matrix(pAdj, ncol = dim(N_gene.res1[[1]])[1])
corrplot(N_gene.cor,type = "lower",na.label = "NA",p.mat = resAdj, sig.level = .05,insig = "blank",tl.col = "black",tl.srt = 45)
N_gene.cor.plot
#####nothing significant######## Done




