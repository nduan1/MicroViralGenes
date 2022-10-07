library(ggplot2)
library(gggenes)
lysogen.df <- read.csv("dataset/lysogenetic_gene.csv", header = T)###this file originated from the clean.filterednov4_cover_table.phages_lysogenic.ffn.tsv
lytic.df <- read.csv("dataset/lytic_genes.csv", header = T)#the similar above
head(lysogen.df)
dim(lysogen.df)
head(lytic.df)
cols_phylum <- c("#4D897C", "#C6B7EC", "#ba892f","#5e5959","#c76d70","#5e5959","#804b0f", "#2d5180", "#0072B2" ,"#00AFBB", "#66A61E", "#E6AB02", "#A6761D","#666666", "#56B4E9","#ABB065", "#999999","#F0E442")
cols_phylum <- cog.62276.df$color_manual


tiff('figures/lysogenetic_amg.tiff', units="in", width=16, height=8, res=300)
ggplot(lysogen.df, aes(xmin = start, xmax = end, y =vOTUs, fill =color_manual, color=protein,forward = strand)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm"),color="black")+
#  geom_gene_label(align = "left")+
#  geom_blank(data = dummies) +
   facet_wrap(~ viral.genome, scales = "free", ncol = 1) +
#  scale_fill_brewer(palette = "Set3") +
  scale_fill_identity(guide = "legend",name="Gene type",labels=c("Recombinase","Functional genes", "Hallmark,signature,connectors","Unknow","AMGs")) +
  ylab("vOTUs")+
  labs(fill="Gene type")+
  theme_genes()
dev.off()
unique(lysogen.df$vOTUs)

tiff('figures/lytic_amg1.tiff', units="in", width=16, height=12, res=300)
ggplot(lytic.df[1:1036,], aes(xmin = start, xmax = end, y =vOTUs, fill =color_manual, color=protein,forward = strand)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm"),color="black")+
  #  geom_gene_label(align = "left")+
  #  geom_blank(data = dummies) +
  facet_wrap(~ viral.genome, scales = "free", ncol = 1) +
  #  scale_fill_brewer(palette = "Set3") +
  scale_fill_identity(guide = "legend",name="Gene type",labels=c("Functional genes", "Hallmark,signature,connectors","Unknow","AMGs")) +
  ylab("vOTUs")+
  labs(fill="Gene type")+
  theme_genes()
dev.off()

tiff('figures/lytic_amg2.tiff', units="in", width=16, height=12, res=300)
ggplot(lytic.df[1037:1987,], aes(xmin = start, xmax = end, y =vOTUs, fill =color_manual, color=protein,forward = strand)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm"),color="black")+
  #  geom_gene_label(align = "left")+
  #  geom_blank(data = dummies) +
  facet_wrap(~ viral.genome, scales = "free", ncol = 1) +
  #  scale_fill_brewer(palette = "Set3") +
  scale_fill_identity(guide = "legend",name="Gene type",labels=c("Functional genes", "Hallmark,signature,connectors","Unknow","AMGs")) +
  ylab("vOTUs")+
  labs(fill="Gene type")+
  theme_genes()
dev.off()

tiff('figures/lytic_amg3.tiff', units="in", width=16, height=12, res=300)
ggplot(lytic.df[1988:2939,], aes(xmin = start, xmax = end, y =vOTUs, fill =color_manual, color=protein,forward = strand)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm"),color="black")+
  #  geom_gene_label(align = "left")+
  #  geom_blank(data = dummies) +
  facet_wrap(~ viral.genome, scales = "free", ncol = 1) +
  #  scale_fill_brewer(palette = "Set3") +
  scale_fill_identity(guide = "legend",name="Gene type",labels=c("Functional genes", "Hallmark,signature,connectors","Unknow","AMGs")) +
  ylab("vOTUs")+
  labs(fill="Gene type")+
  theme_genes()
dev.off()

tiff('figures/lytic_amg4.tiff', units="in", width=16, height=12, res=300)
ggplot(lytic.df[2939:3895,], aes(xmin = start, xmax = end, y =vOTUs, fill =color_manual, color=protein,forward = strand)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm"),color="black")+
  #  geom_gene_label(align = "left")+
  #  geom_blank(data = dummies) +
  facet_wrap(~ viral.genome, scales = "free", ncol = 1) +
  #  scale_fill_brewer(palette = "Set3") +
  scale_fill_identity(guide = "legend",name="Gene type",labels=c("Functional genes", "Hallmark,signature,connectors","Unknow","AMGs")) +
  ylab("vOTUs")+
  labs(fill="Gene type")+
  theme_genes()
dev.off()

tiff('figures/lytic_amg5.tiff', units="in", width=16, height=12, res=300)
ggplot(lytic.df[3896:4859,], aes(xmin = start, xmax = end, y =vOTUs, fill =color_manual, color=protein,forward = strand)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm"),color="black")+
  #  geom_gene_label(align = "left")+
  #  geom_blank(data = dummies) +
  facet_wrap(~ viral.genome, scales = "free", ncol = 1) +
  #  scale_fill_brewer(palette = "Set3") +
  scale_fill_identity(guide = "legend",name="Gene type",labels=c("Functional genes", "Hallmark,signature,connectors","Unknow","AMGs")) +
  ylab("vOTUs")+
  labs(fill="Gene type")+
  theme_genes()
dev.off()

tiff('figures/lytic_amg6.tiff', units="in", width=16, height=12, res=300)
ggplot(lytic.df[4860:5438,], aes(xmin = start, xmax = end, y =vOTUs, fill =color_manual, color=protein,forward = strand)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm"),color="black")+
  #  geom_gene_label(align = "left")+
  #  geom_blank(data = dummies) +
  facet_wrap(~ viral.genome, scales = "free", ncol = 1) +
  #  scale_fill_brewer(palette = "Set3") +
  scale_fill_identity(guide = "legend",name="Gene type",labels=c("Functional genes", "Hallmark,signature,connectors","Unknow","AMGs")) +
  ylab("vOTUs")+
  labs(fill="Gene type")+
  theme_genes()
dev.off()



#too big to draw, so I am planing to subset the data and extract the vOTUs which only contain the AMG gene
#lytic
sub.amg.votu <- function(file){
  nor.matrix <- matrix(nrow=0,ncol =dim(file)[2])
  for (i in 1:dim(file)[1]){
    a <- file[i,9]
    b <- file[i,6]
    if (a=="AMG"){
      df.sub <- subset(file,viral.genome==b)
      nor.matrix <- rbind(nor.matrix,df.sub)
    }
  }
  return(nor.matrix)
}

sub.lytic.df <- sub.amg.votu(lytic.df)
head(sub.lytic.df)
tiff('figures/lytic_amg.tiff', units="in", width=16, height=8, res=300)
ggplot(sub.lytic.df, aes(xmin = start, xmax = end, y =vOTUs, fill =color_manual, color=protein,forward = strand)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm"),color="black")+
  #  geom_gene_label(align = "left")+
  #  geom_blank(data = dummies) +
  facet_wrap(~ viral.genome, scales = "free", ncol = 1) +
  #  scale_fill_brewer(palette = "Set3") +
  scale_fill_identity(guide = "legend",name="Gene type",labels=c("Functional genes", "Hallmark,signature,connectors","Unknow","AMGs")) +
  ylab("vOTUs")+
  labs(fill="Gene type")+
  theme_genes()
dev.off()
########I need to extract the votus with confidence and explore the AMG genes#######
#votu_124
head(sub.lytic.df)
votu_124 <- subset(sub.lytic.df, vOTUs=="vOTU_124")
head(votu_124)
tiff('figures/votu_124.tiff', units="in", width=16, height=8, res=300)
votu_124_plot <- ggplot(votu_124, aes(xmin = start, xmax = end, y =vOTUs, fill =color_manual, color=protein,forward = strand)) +
#  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm"),color="black")+
  geom_gene_arrow(color="black")+
  #geom_gene_label(align = "left")+
  # geom_blank(data = dummies) +
#  facet_wrap(~ viral.genome, scales = "free", ncol = 1) +
  #  scale_fill_brewer(palette = "Set3") +
  scale_fill_identity(guide = "legend",name="Gene type",labels=c("Hallmark,signature,connectors","hypothetical protein","AMGs")) +
  xlab("Length (bp)")+
  ylab("vOTU_124 (17,102 bp)")+
  ggtitle("Starch and sucrose metabolism")+
#  labs(fill="Gene type")+
  theme_genes()+ #%+replace% theme(panel.grid.major.y = element_line(colour = "red")) 
  geom_text(data=votu_124[c(3,11),] %>% mutate(start = (start + end)/2),
            aes(x=start, label = protein), nudge_y = 0.1,size=4,color="#751e1e")+
  geom_text(data=votu_124[c(5),] %>% mutate(start = (start + end)/1.5),
            aes(x=start, label = protein), nudge_y = -0.1,size=4,color="#664c03")+
  geom_segment(data=votu_124,aes(x=18000,y=1.5,xend=18000,yend=0.5), arrow = arrow(length =unit(0.2, "cm")),color="black")+
  annotate("text", x = 18000, y = 1.55, label = "Start")+
  annotate("text", x = 18000, y = 0.45, label = "End")
dev.off()

library(ggpubr)
tiff('figures/amg_gene_3.tiff', units="in", width=17, height=12, res=300)
ggarrange(votu_41_plot,votu_46_plot,votu_124_plot,labels = c("A", "B","C"),ncol = 1, nrow = 3,common.legend = T,legend = "right") 
dev.off()

#lysogenetic

sub.lyso.df <- sub.amg.votu(lysogen.df)
tiff('figures/lyso_amg.tiff', units="in", width=16, height=8, res=300)
ggplot(sub.lyso.df, aes(xmin = start, xmax = end, y =vOTUs, fill =color_manual, color=protein,forward = strand)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm"),color="black")+
  #  geom_gene_label(align = "left")+
  #  geom_blank(data = dummies) +
  facet_wrap(~ viral.genome, scales = "free", ncol = 1) +
  #  scale_fill_brewer(palette = "Set3") +
  scale_fill_identity(guide = "legend",name="Gene type",labels=c("Functional genes", "Hallmark,signature,connectors","Unknow","AMGs")) +
  ylab("vOTUs")+
  labs(fill="Gene type")+
  theme_genes()
dev.off()

#votu_46
head(sub.lyso.df)
votu_46 <- unique(subset(sub.lyso.df, vOTUs=="vOTU_46"))[-75,]#there is overlapping on gene names, so I removed it
head(votu_46)
tiff('figures/votu_46.tiff', units="in", width=16, height=8, res=300)
votu_46_plot <- 
  ggplot(votu_46, aes(xmin = start, xmax = end, y =vOTUs, fill =color_manual, color=protein,forward = strand)) +
  #  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm"),color="black")+
  geom_gene_arrow(color="black")+
  #geom_gene_label(align = "left")+
  # geom_blank(data = dummies) +
  #  facet_wrap(~ viral.genome, scales = "free", ncol = 1) +
  #  scale_fill_brewer(palette = "Set3") +
  scale_fill_identity(guide = "legend",name="Gene type",labels=c("Genes","Hallmark,signature,connectors","hypothetical protein","AMGs")) +
  xlab("Length (bp)")+ #put x before y 
  ylab("vOTU_46 (47,770 bp)")+
  ggtitle("Galactose metabolism")+
  theme_genes()+ #%+replace% theme(panel.grid.major.y = element_line(colour = "red")) 
  geom_text(data=subset(votu_46, gene_type=="hallmark")[1,] %>% mutate(start = (start + end)/0.15),
            aes(x=start, label = protein), nudge_x=1,nudge_y = 0.5,size=4,color="#751e1e")+
  geom_text(data=subset(votu_46, gene_type=="hallmark")[2,] %>% mutate(start = (start + end)/2.5),
            aes(x=start, label = protein), nudge_x=1,nudge_y = 0.45,size=4,color="#751e1e")+
  geom_text(data=subset(votu_46, gene_type=="hallmark")[3,] %>% mutate(start = (start + end)/0.8),
            aes(x=start, label = protein), nudge_x=1,nudge_y = 0.4,size=4,color="#751e1e")+
  geom_text(data=subset(votu_46, gene_type=="hallmark")[4,] %>% mutate(start = (start + end)/2),
            aes(x=start, label = protein), nudge_x=1,nudge_y = 0.35,size=4,color="#751e1e")+
  geom_text(data=subset(votu_46, gene_type=="hallmark")[5,] %>% mutate(start = (start + end)/2),
            aes(x=start, label = protein), nudge_x=1,nudge_y = 0.3,size=4,color="#751e1e")+
  geom_text(data=subset(votu_46, gene_type=="hallmark")[6,] %>% mutate(start = (start + end)/2),
            aes(x=start, label = protein), nudge_x=1,nudge_y = 0.25,size=4,color="#751e1e")+
  geom_text(data=subset(votu_46, gene_type=="hallmark")[7,] %>% mutate(start = (start + end)/2),
            aes(x=start, label = protein), nudge_x=1,nudge_y = 0.2,size=4,color="#751e1e")+
  geom_text(data=subset(votu_46, gene_type=="hallmark")[8,] %>% mutate(start = (start + end)/2),
            aes(x=start, label = protein), nudge_x=1,nudge_y = 0.15,size=4,color="#751e1e")+
  geom_text(data=subset(votu_46, gene_type=="hallmark")[9,] %>% mutate(start = (start + end)/2.25),
            aes(x=start, label = protein), nudge_x=1,nudge_y = 0.1,size=4,color="#751e1e")+
  geom_text(data=subset(votu_46, gene_type=="AMG")[6,] %>% mutate(start = (start + end)/2.5),
            aes(x=start, label = protein), nudge_y = -0.35,size=4,color="#332500")+
  geom_text(data=subset(votu_46, gene_type=="AMG")[5,] %>% mutate(start = (start + end)/3),
            aes(x=start, label = protein), nudge_y = -0.3,size=4,color="#332500")+
  geom_text(data=subset(votu_46, gene_type=="AMG")[4,] %>% mutate(start = (start + end)/3),
            aes(x=start, label = protein), nudge_y = -0.25,size=4,color="#332500")+
  geom_text(data=subset(votu_46, gene_type=="AMG")[3,] %>% mutate(start = (start + end)/4),
            aes(x=start, label = protein), nudge_y = -0.2,size=4,color="#332500")+
  geom_text(data=subset(votu_46, gene_type=="AMG")[2,] %>% mutate(start = (start + end)/3),
            aes(x=start, label = protein), nudge_y = -0.15,size=4,color="#332500")+
  geom_text(data=subset(votu_46, gene_type=="AMG")[1,] %>% mutate(start = (start + end)/3),
            aes(x=start, label = protein), nudge_y = -0.1,size=4,color="#332500")+
  geom_segment(data=votu_46,aes(x=55000,y=1.5,xend=55000,yend=0.5), arrow = arrow(length =unit(0.2, "cm")),color="black")+
  annotate("text", x = 55000, y = 1.55, label = "Start")+
  annotate("text", x = 55000, y = 0.45, label = "End")
dev.off()

#votu_41
head(sub.lytic.df)
votu_41 <- unique(subset(sub.lytic.df, vOTUs=="vOTU_41"))
head(votu_41);dim(votu_41)
tiff('figures/votu_41.tiff', units="in", width=16, height=8, res=300)
votu_41_plot <- ggplot(votu_41, aes(xmin = start, xmax = end, y =vOTUs, fill =color_manual, color=protein,forward = strand)) +
  geom_gene_arrow(color="black")+
  scale_fill_identity(guide = "legend",name="Gene type",labels=c("Genes","Hallmark,signature,connectors","hypothetical protein","AMGs")) +
  xlab("Length (bp)")+
  ylab("vOTU_41 (51,796 bp)")+
  ggtitle("Fructose and mannose metabolism")+
#  labs(fill="Gene type")+
  theme_genes()+
  geom_text(data=subset(votu_41, gene_type=="hallmark")[1,] %>% mutate(start = (start + end)/0.5),
            aes(x=start, label = protein), nudge_x=1,nudge_y = 0.35,size=4,color="#751e1e")+
  geom_text(data=subset(votu_41, gene_type=="hallmark")[2,] %>% mutate(start = (start + end)/1.5),
            aes(x=start, label = protein), nudge_x=1,nudge_y = 0.3,size=4,color="#751e1e")+
  geom_text(data=subset(votu_41, gene_type=="hallmark")[3,] %>% mutate(start = (start + end)/2),
            aes(x=start, label = protein), nudge_x=1,nudge_y = 0.25,size=4,color="#751e1e")+
  geom_text(data=subset(votu_41, gene_type=="hallmark")[4,] %>% mutate(start = (start + end)/2),
            aes(x=start, label = protein), nudge_x=1,nudge_y = 0.2,size=4,color="#751e1e")+
  geom_text(data=subset(votu_41, gene_type=="hallmark")[5,] %>% mutate(start = (start + end)/2),
            aes(x=start, label = protein), nudge_x=1,nudge_y = 0.15,size=4,color="#751e1e")+
  geom_text(data=subset(votu_41, gene_type=="hallmark")[6,] %>% mutate(start = (start + end)/2),
            aes(x=start, label = protein), nudge_x=1,nudge_y = 0.1,size=4,color="#751e1e")+
  geom_text(data=subset(votu_41, gene_type=="AMG")[2,] %>% mutate(start = (start + end)/2),
            aes(x=start, label = protein), nudge_x=1,nudge_y = -0.15,size=4,color="#332500")+
  geom_text(data=subset(votu_41, gene_type=="AMG")[1,] %>% mutate(start = (start + end)/2),
            aes(x=start, label = protein), nudge_x=1,nudge_y = -0.1,size=4,color="#332500")+
  geom_segment(data=votu_46,aes(x=53000,y=1.5,xend=53000,yend=0.5), arrow = arrow(length =unit(0.2, "cm")),color="black")+
  annotate("text", x = 53000, y = 1.55, label = "Start")+
  annotate("text", x = 53000, y = 0.45, label = "End")
dev.off()

###########I need to clearify when writing papers
#possible family_Myoviridae
votu_30 <- subset(lytic.df,vOTUs=="vOTU_30")
votu_30 
votu_63 <- subset(lysogen.df ,vOTUs=="vOTU_63")
votu_93 <- subset(lysogen.df ,vOTUs=="vOTU_93")
votu_149<- subset(lytic.df,vOTUs=="vOTU_149")
votu_149<- subset(lytic.df,vOTUs=="vOTU_149")
votu_229 <-  subset(lytic.df,vOTUs=="vOTU_229");dim(votu_229)#have problem to assign
votu_69 <-  subset(lytic.df,vOTUs=="vOTU_69");dim(votu_69)#have problem to assign
votu_145 <-  subset(lysogen.df,vOTUs=="vOTU_145");dim(votu_145)
votu_100 <-  subset(lytic.df,vOTUs=="vOTU_100");dim(votu_100)
votu_67 <-  subset(lytic.df,vOTUs=="vOTU_67");dim(votu_67)
votu_9 <-  subset(lytic.df,vOTUs=="vOTU_9");dim(votu_9)
View(votu_9)

#################
unique(lysogen.df$vOTUs)
votu_29 <-  subset(lysogen.df,vOTUs=="vOTU_29");dim(votu_29)#didn't find integrase bu recobination protein
votu_91 <-  subset(lysogen.df,vOTUs=="vOTU_91");dim(votu_91)# find integrase
votu_116 <-  subset(lysogen.df,vOTUs=="vOTU_116");dim(votu_116) #find integrase
votu_31 <-  subset(lysogen.df,vOTUs=="vOTU_31");dim(votu_31)# find integrase
#https://www.technologynetworks.com/immunology/articles/lytic-vs-lysogenic-understanding-bacteriophage-life-cycles-308094
votu_52 <-  subset(lysogen.df,vOTUs=="vOTU_52");dim(votu_52)# find integrase
votu_87 <-  subset(lysogen.df,vOTUs=="vOTU_87");dim(votu_87)# find integrase
votu_87 <-  subset(lysogen.df,vOTUs=="vOTU_87");dim(votu_87)# find integrase
votu_49 <-  subset(lysogen.df,vOTUs=="vOTU_49");dim(votu_49)# find integrase
votu_214 <-  subset(lysogen.df,vOTUs=="vOTU_214");dim(votu_214)# find integrase
votu_58 <-  subset(lysogen.df,vOTUs=="vOTU_58");dim(votu_58)# find integrase
votu_141 <-  subset(lysogen.df,vOTUs=="vOTU_141");dim(votu_141)# find integrase
votu_53 <-  subset(lysogen.df,vOTUs=="vOTU_53");dim(votu_53)# find integrase
votu_46 <-  subset(lysogen.df,vOTUs=="vOTU_46");dim(votu_46)# find integrase
votu_74 <-  subset(lysogen.df,vOTUs=="vOTU_74");dim(votu_74)# find integrase
votu_103 <-  subset(lysogen.df,vOTUs=="vOTU_103");dim(votu_103)#find integrase
lyso_list <- data.frame(unique(lysogen.df$vOTUs))
write.csv(lyso_list,"lyso_list.csv")


###########archaea votu
#votu_159
votu_159 <-  subset(lytic.df,vOTUs=="vOTU_159");dim(votu_159)# cannot be classfied as lytic or lysogenic therefore no info in those dataset
votu_128 <-  subset(lytic.df,vOTUs=="vOTU_128");dim(votu_128)
votu_251 <-  subset(lysogen.df,vOTUs=="vOTU_251");dim(votu_251)# cannot be classfied as lytic or lysogenic therefore no info in those dataset
votu_120 <-  subset(lysogen.df,vOTUs=="vOTU_120");dim(votu_120)# cannot be classfied as lytic or lysogenic therefore no info in those dataset

##########
unique(lysogen.df$vOTUs)#number of lysogeny phages 22
votu_60 <- subset(lysogen.df,vOTUs=="vOTU_60")
votu_135 <- subset(lysogen.df,vOTUs=="vOTU_135")
votu_207 <- subset(lysogen.df,vOTUs=="vOTU_207")
votu_39 <- subset(lysogen.df,vOTUs=="vOTU_39")






########votu 119 85 78 62#####
head(sub.lytic.df)
votu_119 <- unique(subset(sub.lytic.df, vOTUs=="vOTU_119"))
head(votu_119);dim(votu_119)
tiff('figures/votu_119.tiff', units="in", width=16, height=8, res=300)
votu_119_plot <- 
  ggplot(votu_119, aes(xmin = start, xmax = end, y =vOTUs, fill =color_manual, color=protein,forward = strand)) +
  geom_gene_arrow(color="black")+
  scale_fill_identity(guide = "legend",name="Gene type",labels=c("Genes","Hallmark,signature,connectors","hypothetical protein","AMGs")) +
  xlab("Length (bp)")+
  ylab("vOTU_119 (17,521 bp)")+
  ggtitle("Amino sugar and nucleotide sugar metabolism")+
  #  labs(fill="Gene type")+
  theme_genes()+
  geom_text(data=subset(votu_119, gene_type=="hallmark")[1,] %>% mutate(start = (start + end)/2),
            aes(x=start, label = protein), nudge_x=1,nudge_y = 0.35,size=4,color="#751e1e")+
  geom_text(data=subset(votu_119, gene_type=="hallmark")[2,] %>% mutate(start = (start + end)/1.5),
            aes(x=start, label = protein), nudge_x=1,nudge_y = 0.3,size=4,color="#751e1e")+
  geom_text(data=subset(votu_119, gene_type=="hallmark")[3,] %>% mutate(start = (start + end)/2),
            aes(x=start, label = protein), nudge_x=1,nudge_y = 0.25,size=4,color="#751e1e")+
  geom_text(data=subset(votu_119, gene_type=="hallmark")[4,] %>% mutate(start = (start + end)/2),
            aes(x=start, label = protein), nudge_x=1,nudge_y = 0.2,size=4,color="#751e1e")+
  geom_text(data=subset(votu_119, gene_type=="hallmark")[5,] %>% mutate(start = (start + end)/2),
            aes(x=start, label = protein), nudge_x=1,nudge_y = 0.15,size=4,color="#751e1e")+
  geom_text(data=subset(votu_119, gene_type=="AMG")[3,] %>% mutate(start = (start + end)/2),
            aes(x=start, label = protein), nudge_x=1,nudge_y = -0.2,size=4,color="#332500")+
  geom_text(data=subset(votu_119, gene_type=="AMG")[2,] %>% mutate(start = (start + end)/2),
            aes(x=start, label = protein), nudge_x=1,nudge_y = -0.15,size=4,color="#332500")+
  geom_text(data=subset(votu_119, gene_type=="AMG")[1,] %>% mutate(start = (start + end)/2),
            aes(x=start, label = protein), nudge_x=1,nudge_y = -0.1,size=4,color="#332500")+
  geom_segment(data=votu_91,aes(x=25000,y=1.5,xend=25000,yend=0.5), arrow = arrow(length =unit(0.2, "cm")),color="black")+
  annotate("text", x = 25000, y = 1.55, label = "Start")+
  annotate("text", x = 25000, y = 0.45, label = "End")
dev.off()

head(sub.lytic.df)
votu_85 <- unique(subset(sub.lytic.df, vOTUs=="vOTU_85"))
head(votu_85);dim(votu_85)
tiff('figures/votu_85.tiff', units="in", width=16, height=8, res=300)

votu_85_plot <- 
  ggplot(votu_85, aes(xmin = start, xmax = end, y =vOTUs, fill =color_manual, color=protein,forward = strand)) +
  geom_gene_arrow(color="black")+
  scale_fill_identity(guide = "legend",name="Gene type",labels=c("Genes","Hallmark,signature,connectors","hypothetical protein","AMGs")) +
  xlab("Length (bp)")+
  ylab("vOTU_85 (23,282 bp)")+
  ggtitle("Galactose metabolism")+
  #  labs(fill="Gene type")+
  theme_genes()+
  geom_text(data=subset(votu_85, gene_type=="hallmark")[1,] %>% mutate(start = (start + end)/0.5),
            aes(x=start, label = protein), nudge_x=1,nudge_y = 0.35,size=4,color="#751e1e")+
  geom_text(data=subset(votu_85, gene_type=="hallmark")[2,] %>% mutate(start = (start + end)/1),
            aes(x=start, label = protein), nudge_x=1,nudge_y = 0.3,size=4,color="#751e1e")+
  geom_text(data=subset(votu_85, gene_type=="hallmark")[3,] %>% mutate(start = (start + end)/1.5),
            aes(x=start, label = protein), nudge_x=1,nudge_y = 0.25,size=4,color="#751e1e")+
  geom_text(data=subset(votu_85, gene_type=="AMG")[1,] %>% mutate(start = (start + end)/2),
            aes(x=start, label = protein), nudge_x=1,nudge_y = -0.1,size=4,color="#332500")+
  geom_segment(data=votu_91,aes(x=25000,y=1.5,xend=25000,yend=0.5), arrow = arrow(length =unit(0.2, "cm")),color="black")+
  annotate("text", x = 25000, y = 1.55, label = "Start")+
  annotate("text", x = 25000, y = 0.45, label = "End")
dev.off()
  
  
head(sub.lytic.df)
votu_78 <- unique(subset(sub.lytic.df, vOTUs=="vOTU_78"))
head(votu_78);dim(votu_78)
tiff('figures/votu_78.tiff', units="in", width=16, height=8, res=300)
votu_78_plot <- 
  ggplot(votu_78, aes(xmin = start, xmax = end, y =vOTUs, fill =color_manual, color=protein,forward = strand)) +
  geom_gene_arrow(color="black")+
  scale_fill_identity(guide = "legend",name="Gene type",labels=c("Genes","Hallmark,signature,connectors","hypothetical protein","AMGs")) +
  xlab("Length (bp)")+
  ylab("vOTU_78 (26,849 bp)")+
  ggtitle("Amino sugar and nucleotide sugar metabolism")+
  #  labs(fill="Gene type")+
  theme_genes()+
  geom_text(data=subset(votu_78, gene_type=="hallmark")[1,] %>% mutate(start = (start + end)/0.3),
            aes(x=start, label = protein), nudge_x=1,nudge_y = 0.35,size=4,color="#751e1e")+
  geom_text(data=subset(votu_78, gene_type=="hallmark")[2,] %>% mutate(start = (start + end)/0.8),
            aes(x=start, label = protein), nudge_x=1,nudge_y = 0.3,size=4,color="#751e1e")+
    geom_text(data=subset(votu_78, gene_type=="AMG")[1,] %>% mutate(start = (start + end)/2),
            aes(x=start, label = protein), nudge_x=1,nudge_y = -0.1,size=4,color="#332500")+
  geom_segment(data=votu_91,aes(x=29000,y=1.5,xend=29000,yend=0.5), arrow = arrow(length =unit(0.2, "cm")),color="black")+
  annotate("text", x = 29000, y = 1.55, label = "Start")+
  annotate("text", x = 29000, y = 0.45, label = "End")
dev.off()


votu_62 <- unique(subset(sub.lytic.df, vOTUs=="vOTU_62"))
head(votu_62);dim(votu_62)
tiff('figures/votu_62.tiff', units="in", width=16, height=8, res=300)
votu_62_plot <- 
  ggplot(votu_62, aes(xmin = start, xmax = end, y =vOTUs, fill =color_manual, color=protein,forward = strand)) +
  geom_gene_arrow(color="black")+
  scale_fill_identity(guide = "legend",name="Gene type",labels=c("Genes","Hallmark,signature,connectors","hypothetical protein","AMGs")) +
  xlab("Length (bp)")+
  ylab("vOTU_62 (41,319 bp)")+
  ggtitle("Amino sugar and nucleotide sugar metabolism")+
  #  labs(fill="Gene type")+
  theme_genes()+
  geom_text(data=subset(votu_62, gene_type=="hallmark")[1,] %>% mutate(start = (start + end)/0.6),
            aes(x=start, label = protein), nudge_x=1,nudge_y = 0.35,size=4,color="#751e1e")+
  geom_text(data=subset(votu_62, gene_type=="hallmark")[2,] %>% mutate(start = (start + end)/1.8),
            aes(x=start, label = protein), nudge_x=1,nudge_y = 0.3,size=4,color="#751e1e")+
  geom_text(data=subset(votu_62, gene_type=="hallmark")[3,] %>% mutate(start = (start + end)/1.9),
            aes(x=start, label = protein), nudge_x=1,nudge_y = 0.25,size=4,color="#751e1e")+
  geom_text(data=subset(votu_62, gene_type=="hallmark")[4,] %>% mutate(start = (start + end)/3),
            aes(x=start, label = protein), nudge_x=1,nudge_y = 0.2,size=4,color="#751e1e")+
  geom_text(data=subset(votu_62, gene_type=="AMG")[2,] %>% mutate(start = (start + end)/0.5),
            aes(x=start, label = protein), nudge_x=1,nudge_y = -0.15,size=4,color="#332500")+
  geom_text(data=subset(votu_62, gene_type=="AMG")[1,] %>% mutate(start = (start + end)/0.3),
            aes(x=start, label = protein), nudge_x=1,nudge_y = -0.1,size=4,color="#332500")+
  geom_segment(data=votu_91,aes(x=52000,y=1.5,xend=52000,yend=0.5), arrow = arrow(length =unit(0.2, "cm")),color="black")+
  annotate("text", x = 52000, y = 1.55, label = "Start")+
  annotate("text", x = 52000, y = 0.45, label = "End")
dev.off()


######amino acid metabolism
library(ggpubr)
tiff('figures/amg_gene_78_119_62.tiff', units="in", width=17, height=12, res=300)
ggarrange(votu_119_plot,votu_78_plot,votu_62_plot,labels = c("A", "B","C"),ncol = 1, nrow = 3,common.legend = T,legend = "right") 
dev.off()

