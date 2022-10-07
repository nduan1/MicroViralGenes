library(readxl)
library(tidyverse)
library(fs)
M107_carbon <- read_excel("raw_dataset/M107_metabolism_summary.xlsx", sheet = "carbon utilization")
head(M107_carbon)

M112_carbon <- read_excel("raw_dataset/M112_metabolism_summary.xlsx", sheet = "carbon utilization")
head(M112_carbon)

M16_carbon <- read_excel("raw_dataset/M16_metabolism_summary.xlsx", sheet = "carbon utilization")
head(M16_carbon)

M35_carbon <- read_excel("raw_dataset/M35_metabolism_summary.xlsx", sheet = "carbon utilization")
head(M35_carbon)

M37_carbon <- read_excel("raw_dataset/M37_metabolism_summary.xlsx", sheet = "carbon utilization")
head(M37_carbon)

M51_carbon <- read_excel("raw_dataset/M51_metabolism_summary.xlsx", sheet = "carbon utilization")
head(M51_carbon)

M56_carbon <- read_excel("raw_dataset/M56_metabolism_summary.xlsx", sheet = "carbon utilization")
head(M56_carbon)

M62_carbon <- read_excel("raw_dataset/M62_metabolism_summary.xlsx", sheet = "carbon utilization")
head(M62_carbon)

M63_carbon <- read_excel("raw_dataset/M63_metabolism_summary.xlsx", sheet = "carbon utilization")
head(M63_carbon)

M77_carbon <- read_excel("raw_dataset/M77_metabolism_summary.xlsx", sheet = "carbon utilization")
head(M77_carbon)

M79_carbon <- read_excel("raw_dataset/M79_metabolism_summary.xlsx", sheet = "carbon utilization")
head(M79_carbon)

M9_carbon <- read_excel("raw_dataset/M9_metabolism_summary.xlsx", sheet = "carbon utilization")
head(M9_carbon)

#######combine all the dataset because they have the same rows#######
C_gene <- cbind(M107_carbon,M112_carbon,M16_carbon,M35_carbon,
                M37_carbon,M51_carbon,M56_carbon,M62_carbon,
                M63_carbon,M77_carbon,M79_carbon,M9_carbon)
colnames(C_gene)
#double check it
for (i in 1:11){
  a <- identical(C_gene[,1],C_gene[,1+6*i])
  print(a)
}

for (i in 1:11){
  a <- identical(C_gene[,2],C_gene[,2+6*i])
  print(a)
}

for (i in 1:11){
  a <- identical(C_gene[,3],C_gene[,3+6*i])
  print(a)
}

for (i in 1:11){
  a <- identical(C_gene[,4],C_gene[,4+6*i])
  print(a)
}

for (i in 1:11){
  a <- identical(C_gene[,5],C_gene[,5+6*i])
  print(a)
}

for (i in 1:11){
  a <- identical(C_gene[,6],C_gene[,6+6*i])
  print(a)
}

#####therefore, we could delete some of the datacols########

colnames(C_gene)
C_gene.df <- C_gene[,c(1:6,12,18,24,30,36,42,48,54,60,66,72)];dim(C_gene.df)#subheader is useless, so removed it
head(C_gene.df)
colnames(C_gene.df)
colnames(C_gene.df)[6:17] <- c("NCNTN0_4","VNTN0_4","VNTN60_2","NCNTN0_1","VNTN0_1","VNTN60_1",
                                  "NCNTN60_1","VNTN0_2","NCNTN0_2","NCNTN60_4","VNTN60_4","NCNTN60_2")
head(C_gene.df)
#I realize there is "Deleted"data in our dataset which is useless, and I didn't found any data in "Deleted" rows, therefore we could just delete it
C_gene.df.1 <- C_gene.df %>% 
  filter(!grepl("Deleted",gene_description))
C_gene.df.clean <- C_gene.df.1[rowSums(C_gene.df.1[6:17])>0,]#remove the rows that contains only 0
dim(C_gene.df.clean)#869  17
#####relative abundance based upon all the gene in carbon utilize########
r_c_gene <- apply(C_gene.df.clean[6:17],2,function(x) x/sum(x));head(r_c_gene)
dim(r_c_gene)
sum_gene <- apply(C_gene.df.clean[6:17],2,sum);sum_gene#give you all the number of gene
r_c_gene.df<- cbind(C_gene.df.clean[1:5],r_c_gene);dim(r_c_gene.df)
####subset dataset based on the different header##########
unique(r_c_gene.df$header)
r_cazy_gene <- subset(r_c_gene.df,header=="CAZY");dim(r_cazy_gene)#228  17#Carbohydrate-Active Enzymes (CAZymes)
r_central_gene <- subset(r_c_gene.df,header=="central carbon");dim(r_central_gene)#195  17
r_hydrocarbon_gene <- subset(r_c_gene.df,header=="hydrocarbon degradation");dim(r_hydrocarbon_gene)#68 17
r_pyruvate_gene <- subset(r_c_gene.df,header=="pyruvate metabolism");dim(r_pyruvate_gene)#30 17
r_na_gene <- r_c_gene.df[c(520:522),]

#1. now I am going to extract the starch from cazy dataset#
r_starch_cazy_gene <- r_cazy_gene %>% 
  filter(grepl("starch",gene_description)|grepl("Starch",subheader)|grepl("Starch",gene_description)|grepl("starch",subheader))
r_starch_cazy_gene.1 <-mutate(r_starch_cazy_gene,group=c(rep("Starch",dim(r_starch_cazy_gene)[1])))
head(r_starch_cazy_gene)
dim(r_starch_cazy_gene)#22 18

#2. cellulose and hemicellulose
r_allcellose_cazy_gene <- r_cazy_gene %>% 
  filter(grepl("cellulose",gene_description)|grepl("cellulose",subheader)|grepl("Cellulose",subheader)|grepl("Cellulose",gene_description))
dim(r_allcellose_cazy_gene)#50
#2.1
r_hemicellose_cazy_gene <- r_allcellose_cazy_gene %>% 
  filter(grepl("hemicellulose",gene_description)|grepl("hemicellulose",subheader)|grepl("Hemicellulose",subheader)|grepl("Hemicellulose",gene_description)) 
r_hemicellose_cazy_gene.1 <-mutate(r_hemicellose_cazy_gene,group=c(rep("Hemicellulose",dim(r_hemicellose_cazy_gene)[1])))
dim(r_hemicellose_cazy_gene)#36
r_cellose_cazy_gene <- r_allcellose_cazy_gene %>% 
  filter(!(grepl("hemicellulose",gene_description)|grepl("hemicellulose",subheader)|grepl("Hemicellulose",subheader)|grepl("Hemicellulose",gene_description))) 
r_cellose_cazy_gene.1 <-mutate(r_cellose_cazy_gene,group=c(rep("Cellulose",dim(r_cellose_cazy_gene)[1])))
dim(r_cellose_cazy_gene)#14

#3. Chitin
r_chitin_cazy_gene <- r_cazy_gene %>% 
  filter(grepl("chitin",gene_description)|grepl("Chitin",subheader)|grepl("Chitin",gene_description)|grepl("chitin",subheader))
r_chitin_cazy_gene.1 <- mutate(r_chitin_cazy_gene,group=c(rep("Chitin",dim(r_chitin_cazy_gene)[1])))
dim(r_chitin_cazy_gene)#20 17

#4. lignin in cazy group and other lignin in the hydrocarbon degradation groups
r_lignin_cazy_gene <- r_cazy_gene %>% 
  filter(grepl("lignin",gene_description)|grepl("Lignin",subheader)|grepl("Lignin",gene_description)|grepl("lignin",subheader))
dim(r_lignin_cazy_gene)#1 17

r_lignin_hd_gene <- r_hydrocarbon_gene %>% 
  filter(grepl("lignin",gene_description)|grepl("Lignin",module)|grepl("Lignin",gene_description)|grepl("lignin",module))
dim(r_lignin_hd_gene)##5 17 

r_lignin_cazy_gene.1 <- mutate(rbind(r_lignin_cazy_gene,r_lignin_hd_gene),group=c(rep("Lignin",dim(rbind(r_lignin_cazy_gene,r_lignin_hd_gene))[1])))
#6. Pectin
r_pectin_cazy_gene <- r_cazy_gene %>% 
  filter(grepl("pectin",gene_description)|grepl("Pectin",module)|grepl("Pectin",gene_description)|grepl("pectin",module))
dim(r_pectin_cazy_gene)#4 17
r_pectin_cazy_gene.1 <- mutate(r_pectin_cazy_gene,group=c(rep("Pectin",dim(r_pectin_cazy_gene)[1])))
#7. aromatic, go to the hydrocarbon cate.#I realize that all the function in this category are all the aromatic with benzene rings, which also include the lignin
r_aromatic_hd_exlignin <- r_hydrocarbon_gene %>% filter(!grepl("lignin",module))
r_aromatic_hd_exlignin.1 <- mutate(r_aromatic_hd_exlignin,group=c(rep("Aromatic",dim(r_aromatic_hd_exlignin)[1])))
dim(r_aromatic_hd_exlignin)#63 17

######I will start making figures for the different carbon degradation####
carbon.cat.dat <- rbind(r_starch_cazy_gene.1,r_hemicellose_cazy_gene.1,r_cellose_cazy_gene.1,r_chitin_cazy_gene.1,r_pectin_cazy_gene.1,r_aromatic_hd_exlignin.1,r_lignin_cazy_gene.1 )
dim(carbon.cat.dat)
sample.data.longer <- pivot_longer(carbon.cat.dat,cols = c("NCNTN0_4","VNTN0_4","VNTN60_2","NCNTN0_1","VNTN0_1","VNTN60_1",
                                                           "NCNTN60_1","VNTN0_2","NCNTN0_2","NCNTN60_4","VNTN60_4","NCNTN60_2"))
head(sample.data.longer)
dim(carbon.cat.dat)
#view(carbon.cat.dat)
carbon.cat.dat.df <- carbon.cat.dat[,c(1:5,18,6:17)]
colnames(carbon.cat.dat.df)
sample.data <- read.csv("dataset/meta_data.csv")
head(sample.data)
carbon.cat.boxplot <- inner_join(sample.data,sample.data.longer,by=c("Sample"="name")) 
carbon.cat.boxplot$id_description_module <- apply(carbon.cat.boxplot[,c(9,10)],1,paste,collapse="_")
colnames(carbon.cat.boxplot)
##all together too messy

#1.starch
tiff('figures/box_c_degre_gene_treat.tiff', units="in", width=10, height=6, res=300)
starch.carbon.cat.box <- subset(carbon.cat.boxplot,group=="Starch") %>% 
  mutate(gene_name=c(rep(unique(starch.carbon.cat.box$gene_id),12)))
ggplot(starch.carbon.cat.box, aes(x =gene_id, y = value, fill=Treat))+
  geom_boxplot(size=0.3,width=0.5,alpha=1,notch = FALSE,na.rm=TRUE,outlier.shape=NA,position=position_dodge(preserve = "single"), color="black")+
  geom_jitter(aes(color=Treat),shape=16,size=0.5,position=position_jitterdodge(0.2),na.rm = TRUE)+
  labs(title="Starch",x="Gene ID",y="Relative abundance",size=12)+
  scale_fill_manual(values=c("aquamarine4","royalblue4","grey","#e6b325"))+
  scale_color_manual(values = c("aquamarine4","royalblue4","grey","#e6b325"))+
  theme(legend.position="bottom",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_line(colour="grey",size = 0.1),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=12,family = "Arial",face="plain",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="plain",size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(size=12,family="Arial",face = "bold"),
        axis.title=element_text(size = 12,face="plain",family = "Arial",vjust = 1),
        legend.title = element_text(size=12,face="plain",family = "Arial"),
        legend.text = (element_text(size=12,family = "Arial")))
dev.off()
#####Starche stats#####
#CBM20
CBM20 <- subset(starch.carbon.cat.box,gene_name=="CBM20");dim(CBM20)
CBM20.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = CBM20)
Anova(CBM20.mod)#
summary(CBM20.mod)
coef(CBM20.mod)

#CBM25
CBM25 <- subset(starch.carbon.cat.box,gene_name=="CBM25");dim(CBM25)
CBM25.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = CBM25)
Anova(CBM25.mod)#
summary(CBM25.mod)
coef(CBM25.mod)

#CBM34
CBM34 <- subset(starch.carbon.cat.box,gene_name=="CBM34");dim(CBM34)
CBM34.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = CBM34)
Anova(CBM34.mod)#
summary(CBM34.mod)
coef(CBM34.mod)

#GH126
GH126 <- subset(starch.carbon.cat.box,gene_name=="GH126");dim(GH126)
GH126.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = GH126)
Anova(GH126.mod)#
summary(GH126.mod)
coef(GH126.mod)

#GH13
normalTest(GH13$value,"sw")#good
GH13 <- subset(starch.carbon.cat.box,gene_name=="GH13");dim(GH13)
GH13.mod <- lmer(value~ Nitrogen +(1|block), data = GH13)
Anova(GH13.mod)#Nitrogen<no
summary(GH13.mod)
coef(GH13.mod)

#GH133
normalTest(GH133$value,"sw")#good
GH133 <- subset(starch.carbon.cat.box,gene_name=="GH133");dim(GH133)
GH133.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = GH133)
Anova(GH133.mod)#Cover:Nitrogen
summary(GH133.mod)
coef(GH133.mod)
nocover.GH133<- subset(GH133, Cover=="No_cover")
GH133.mod.1 <- lmer(value~ Nitrogen +(1|block), data = nocover.GH133)
nitrogen.GH133.1<-glht(GH133.mod.1,linfct=mcp(Nitrogen="Tukey"))
summary(nitrogen.GH133.1)#fer>no
cover.GH133<- subset(GH133, Cover=="Vetch")
GH133.mod.2 <- lmer(value~ Nitrogen +(1|block), data = cover.GH133)
nitrogen.GH133.2<-glht(GH133.mod.2,linfct=mcp(Nitrogen="Tukey"))
summary(nitrogen.GH133.2)# fer < no

n0.GH133<- subset(GH133, Nitrogen=="N0")
GH133.mod.3 <- lmer(value~ Cover +(1|block), data = n0.GH133)
cover.GH133.3<-glht(GH133.mod.3,linfct=mcp(Cover="Tukey"))
summary(cover.GH133.3)#no sig
n60.GH133<- subset(GH133,Nitrogen=="N60")
GH133.mod.4 <- lmer(value~ Cover +(1|block), data = n60.GH133)
cover.GH133.4<-glht(GH133.mod.4,linfct=mcp(Cover="Tukey"))
summary(cover.GH133.4)#no > vetch
#GH15
GH15 <- subset(starch.carbon.cat.box,gene_name=="GH15");dim(GH15)
GH15.mod <- lmer(value~ Nitrogen +(1|block), data = GH15)
Anova(GH15.mod)#Nitrogen>no
summary(GH15.mod)
coef(GH15.mod)

#GH57
GH57 <- subset(starch.carbon.cat.box,gene_name=="GH57");dim(GH57)
GH57.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = GH57)
Anova(GH57.mod)#
summary(GH57.mod)
coef(GH57.mod)

#GH97
GH97 <- subset(starch.carbon.cat.box,gene_name=="GH97");dim(GH97)
GH97.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = GH97)
Anova(GH97.mod)#
summary(GH97.mod)
coef(GH97.mod)

#GT35
normalTest(GT35$value,'sw')#not good enough
GT35 <- subset(starch.carbon.cat.box,gene_name=="GT35");dim(GT35)
GT35.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = GT35)
Anova(GT35.mod)#Cover:Nitrogen
summary(GT35.mod)
coef(GT35.mod)
nocover.GT35<- subset(GT35, Cover=="No_cover")
GT35.mod.1 <- lmer(value~ Nitrogen +(1|block), data = nocover.GT35)
nitrogen.GT35.1<-glht(GT35.mod.1,linfct=mcp(Nitrogen="Tukey"))
summary(nitrogen.GT35.1)#no sig
cover.GT35<- subset(GT35, Cover=="Vetch")
GT35.mod.2 <- lmer(value~ Nitrogen +(1|block), data = cover.GT35)
nitrogen.GT35.2<-glht(GT35.mod.2,linfct=mcp(Nitrogen="Tukey"))
summary(nitrogen.GT35.2)# NO SIG

n0.GT35<- subset(GT35, Nitrogen=="N0")
GT35.mod.3 <- lmer(value~ Cover +(1|block), data = n0.GT35)
cover.GT35.3<-glht(GT35.mod.3,linfct=mcp(Cover="Tukey"))
summary(cover.GT35.3)#no sig
n60.GT35<- subset(GT35,Nitrogen=="N60")
GT35.mod.4 <- lmer(value~ Cover +(1|block), data = n60.GT35)
cover.GT35.4<-glht(GT35.mod.4,linfct=mcp(Cover="Tukey"))
summary(cover.GT35.4)# no sig

#GT5
normalTest(GT5$value,'sw')#good
GT5 <- subset(starch.carbon.cat.box,gene_name=="GT5");dim(GT5)
GT5.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = GT5)
Anova(GT5.mod)#Cover:Nitrogen
summary(GT5.mod)
coef(GT5.mod)

nocover.GT5<- subset(GT5, Cover=="No_cover")
GT5.mod.1 <- lmer(value~ Nitrogen +(1|block), data = nocover.GT5)
nitrogen.GT5.1<-glht(GT5.mod.1,linfct=mcp(Nitrogen="Tukey"))
summary(nitrogen.GT5.1)#no >fer
cover.GT5<- subset(GT5, Cover=="Vetch")
GT5.mod.2 <- lmer(value~ Nitrogen +(1|block), data = cover.GT5)
nitrogen.GT5.2<-glht(GT5.mod.2,linfct=mcp(Nitrogen="Tukey"))
summary(nitrogen.GT5.2)# NO SIG

n0.GT5<- subset(GT5, Nitrogen=="N0")
GT5.mod.3 <- lmer(value~ Cover +(1|block), data = n0.GT5)
cover.GT5.3<-glht(GT5.mod.3,linfct=mcp(Cover="Tukey"))
summary(cover.GT5.3)#vetch<no
n60.GT5<- subset(GT5,Nitrogen=="N60")
GT5.mod.4 <- lmer(value~ Cover +(1|block), data = n60.GT5)
cover.GT5.4<-glht(GT5.mod.4,linfct=mcp(Cover="Tukey"))
summary(cover.GT5.4)# no sig


#2.cellulose
cellulose.carbon.cat.box <- subset(carbon.cat.boxplot,group=="Cellulose") %>% 
  mutate(gene_name=c(rep(unique(cellulose.carbon.cat.box$gene_id),12)))
tiff('figures/box_c_degre_gene_treat.tiff', units="in", width=10, height=6, res=300)
ggplot(cellulose.carbon.cat.box, aes(x =gene_id, y = value, fill=Treat))+
  geom_boxplot(size=0.3,width=0.5,alpha=1,notch = FALSE,na.rm=TRUE,outlier.shape=NA,position=position_dodge(preserve = "single"), color="black")+
  geom_jitter(aes(color=Treat),shape=16,size=0.5,position=position_jitterdodge(0.2),na.rm = TRUE)+
  labs(title="Cellulose",x="Gene ID",y="Relative abundance",size=12)+
  scale_fill_manual(values=c("aquamarine4","royalblue4","grey","#e6b325"))+
  scale_color_manual(values = c("aquamarine4","royalblue4","grey","#e6b325"))+
  theme(legend.position="bottom",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_line(colour="grey",size = 0.1),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=12,family = "Arial",face="plain",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="plain",size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(size=12,family="Arial",face = "bold"),
        axis.title=element_text(size = 12,face="plain",family = "Arial",vjust = 1),
        legend.title = element_text(size=12,face="plain",family = "Arial"),
        legend.text = (element_text(size=12,family = "Arial")))
dev.off()
#####Cellulose stats#####
#AA10
AA10 <- subset(cellulose.carbon.cat.box,gene_name=="AA10");dim(AA10)
AA10.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = AA10)
Anova(AA10.mod)#
summary(AA10.mod)
coef(AA10.mod)

#CBM2
normalTest(CBM2$value,"sw")#not good
CBM2 <- subset(cellulose.carbon.cat.box,gene_name=="CBM2");dim(CBM2)
CBM2.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = CBM2)
Anova(CBM2.mod)#Cover:Nitrogen
summary(CBM2.mod)
coef(CBM2.mod)

#CBM3
CBM3 <- subset(cellulose.carbon.cat.box,gene_name=="CBM3");dim(CBM3)
CBM3.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = CBM3)
Anova(CBM3.mod)#
summary(CBM3.mod)
coef(CBM3.mod)

#CBM4
CBM4 <- subset(cellulose.carbon.cat.box,gene_name=="CBM4");dim(CBM4)
CBM4.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = CBM4)
Anova(CBM4.mod)#
summary(CBM4.mod)
coef(CBM4.mod)

#CBM44
CBM44 <- subset(cellulose.carbon.cat.box,gene_name=="CBM44");dim(CBM44)
CBM44.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = CBM44)
Anova(CBM44.mod)#
summary(CBM44.mod)
coef(CBM44.mod)

#CBM46
CBM46 <- subset(cellulose.carbon.cat.box,gene_name=="CBM46");dim(CBM46)
CBM46.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = CBM46)
Anova(CBM46.mod)#
summary(CBM46.mod)
coef(CBM46.mod)

#CBM6
CBM6 <- subset(cellulose.carbon.cat.box,gene_name=="CBM6");dim(CBM6)
CBM6.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = CBM6)
Anova(CBM6.mod)#
summary(CBM6.mod)
coef(CBM6.mod)

#CBM63
CBM63 <- subset(cellulose.carbon.cat.box,gene_name=="CBM63");dim(CBM63)
CBM63.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = CBM63)
Anova(CBM63.mod)#
summary(CBM63.mod)
coef(CBM63.mod)

#CBM8
CBM8 <- subset(cellulose.carbon.cat.box,gene_name=="CBM8");dim(CBM8)
CBM8.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = CBM8)
Anova(CBM8.mod)#
summary(CBM8.mod)
coef(CBM8.mod)

#CBM9
normalTest(CBM9$value,"sw")#good
CBM9 <- subset(cellulose.carbon.cat.box,gene_name=="CBM9");dim(CBM9)
CBM9.mod <- lmer(value~ Cover+(1|block), data = CBM9)
Anova(CBM9.mod)#no Cover> vetch
summary(CBM9.mod)
coef(CBM9.mod)

#GH116
GH116 <- subset(cellulose.carbon.cat.box,gene_name=="GH116");dim(GH116)
GH116.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = GH116)
Anova(GH116.mod)#
summary(GH116.mod)
coef(GH116.mod)

#GH6
GH6 <- subset(cellulose.carbon.cat.box,gene_name=="GH6");dim(GH6)
GH6.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = GH6)
Anova(GH6.mod)#
summary(GH6.mod)
coef(GH6.mod)

#GH94
GH94 <- subset(cellulose.carbon.cat.box,gene_name=="GH94");dim(GH94)
GH94.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = GH94)
Anova(GH94.mod)#
summary(GH94.mod)
coef(GH94.mod)

#GT2
GT2 <- subset(cellulose.carbon.cat.box,gene_name=="GT2");dim(GT2)
GT2.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = GT2)
Anova(GT2.mod)#
summary(GT2.mod)
coef(GT2.mod)

#3.hemicellulose
hemicellulose.carbon.cat.box <- subset(carbon.cat.boxplot,group=="Hemicellulose") %>% 
  mutate(gene_name=c(rep(unique(hemicellulose.carbon.cat.box$gene_id),12)))
tiff('figures/box_c_degre_gene_treat.tiff', units="in", width=10, height=6, res=300)
ggplot(hemicellulose.carbon.cat.box, aes(x =gene_name, y = value, fill=Treat))+
  geom_boxplot(size=0.3,width=0.5,alpha=1,notch = FALSE,na.rm=TRUE,outlier.shape=NA,position=position_dodge(preserve = "single"), color="black")+
  geom_jitter(aes(color=Treat),shape=16,size=0.5,position=position_jitterdodge(0.2),na.rm = TRUE)+
  labs(title="Hemicellulose",x="Gene ID",y="Relative abundance",size=12)+
  scale_fill_manual(values=c("aquamarine4","royalblue4","grey","#e6b325"))+
  scale_color_manual(values = c("aquamarine4","royalblue4","grey","#e6b325"))+
  theme(legend.position="bottom",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_line(colour="grey",size = 0.1),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=12,family = "Arial",face="plain",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="plain",size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(size=12,family="Arial",face = "bold"),
        axis.title=element_text(size = 12,face="plain",family = "Arial",vjust = 1),
        legend.title = element_text(size=12,face="plain",family = "Arial"),
        legend.text = (element_text(size=12,family = "Arial")))
dev.off()
#####Hemicellulose stats#####
#GH1
GH1 <- subset(hemicellulose.carbon.cat.box,gene_name=="GH1");dim(GH1)
GH1.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = GH1)
Anova(GH1.mod)#
summary(GH1.mod)
coef(GH1.mod)

#GH10
GH10 <- subset(hemicellulose.carbon.cat.box,gene_name=="GH10");dim(GH10)
GH10.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = GH10)
Anova(GH10.mod)#
summary(GH10.mod)
coef(GH10.mod)

#GH11
GH11 <- subset(hemicellulose.carbon.cat.box,gene_name=="GH11");dim(GH11)
GH11.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = GH11)
Anova(GH11.mod)#
summary(GH11.mod)
coef(GH11.mod)

#GH113
normalTest(GH113$value,"sw")#not good
GH113 <- subset(hemicellulose.carbon.cat.box,gene_name=="GH113");dim(GH113)
GH113.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = GH113)
Anova(GH113.mod)#Cover:Nitrogen
summary(GH113.mod)
coef(GH113.mod)
nocover.GH113<- subset(GH113, Cover=="No_cover")
GH113.mod.1 <- lmer(value~ Nitrogen +(1|block), data = nocover.GH113)
nitrogen.GH113.1<-glht(GH113.mod.1,linfct=mcp(Nitrogen="Tukey"))
summary(nitrogen.GH113.1)#no >fer
cover.GH113<- subset(GH113, Cover=="Vetch")
GH113.mod.2 <- lmer(value~ Nitrogen +(1|block), data = cover.GH113)
nitrogen.GH113.2<-glht(GH113.mod.2,linfct=mcp(Nitrogen="Tukey"))
summary(nitrogen.GH113.2)# fer > no

n0.GH113<- subset(GH113, Nitrogen=="N0")
GH113.mod.3 <- lmer(value~ Cover +(1|block), data = n0.GH113)
cover.GH113.3<-glht(GH113.mod.3,linfct=mcp(Cover="Tukey"))
summary(cover.GH113.3)#vetch<no
n60.GH113<- subset(GH113,Nitrogen=="N60")
GH113.mod.4 <- lmer(value~ Cover +(1|block), data = n60.GH113)
cover.GH113.4<-glht(GH113.mod.4,linfct=mcp(Cover="Tukey"))
summary(cover.GH113.4)#vetch>no
#GH115
GH115 <- subset(hemicellulose.carbon.cat.box,gene_name=="GH115");dim(GH115)
GH115.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = GH115)
Anova(GH115.mod)#
summary(GH115.mod)
coef(GH115.mod)

#GH12
GH12 <- subset(hemicellulose.carbon.cat.box,gene_name=="GH12");dim(GH12)
GH12.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = GH12)
Anova(GH12.mod)#
summary(GH12.mod)
coef(GH12.mod)

#GH120
GH120 <- subset(hemicellulose.carbon.cat.box,gene_name=="GH120");dim(GH120)
GH120.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = GH120)
Anova(GH120.mod)#
summary(GH120.mod)
coef(GH120.mod)

#GH125
GH125 <- subset(hemicellulose.carbon.cat.box,gene_name=="GH125");dim(GH125)
GH125.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = GH125)
Anova(GH125.mod)#
summary(GH125.mod)
coef(GH125.mod)

#GH130
GH130 <- subset(hemicellulose.carbon.cat.box,gene_name=="GH130");dim(GH130)
GH130.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = GH130)
Anova(GH130.mod)#
summary(GH130.mod)
coef(GH130.mod)

#GH141
normalTest(GH141$value,"sw")#not good so you need to double check it
GH141 <- subset(hemicellulose.carbon.cat.box,gene_name=="GH141");dim(GH141)
GH141.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = GH141)
Anova(GH141.mod)#Cover:Nitrogen 
summary(GH141.mod)
coef(GH141.mod)
nocover.GH141<- subset(GH141, Cover=="No_cover")
GH141.mod.1 <- lmer(value~ Nitrogen +(1|block), data = nocover.GH141)
nitrogen.GH141.1<-glht(GH141.mod.1,linfct=mcp(Nitrogen="Tukey"))
summary(nitrogen.GH141.1)#no sig
cover.GH141<- subset(GH141, Cover=="Vetch")
GH141.mod.2 <- lmer(value~ Nitrogen +(1|block), data = cover.GH141)
nitrogen.GH141.2<-glht(GH141.mod.2,linfct=mcp(Nitrogen="Tukey"))
summary(nitrogen.GH141.2)# no sig

n0.GH141<- subset(GH141, Nitrogen=="N0")
GH141.mod.3 <- lmer(value~ Cover +(1|block), data = n0.GH141)
cover.GH141.3<-glht(GH141.mod.3,linfct=mcp(Cover="Tukey"))
summary(cover.GH141.3)#no sig
n60.GH141<- subset(GH141,Nitrogen=="N60")
GH141.mod.4 <- lmer(value~ Cover +(1|block), data = n60.GH141)
cover.GH141.4<-glht(GH141.mod.4,linfct=mcp(Cover="Tukey"))
summary(cover.GH141.4)#no sig
#GH158
GH158 <- subset(hemicellulose.carbon.cat.box,gene_name=="GH158");dim(GH158)
GH158.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = GH158)
Anova(GH158.mod)#
summary(GH158.mod)
coef(GH158.mod)

#GH16
normalTest(GH16$value,"sw")
GH16 <- subset(hemicellulose.carbon.cat.box,gene_name=="GH16");dim(GH16)
GH16.mod <- lmer(value~ Nitrogen +(1|block), data = GH16)
Anova(GH16.mod)#Nitrogen>no 
summary(GH16.mod)
coef(GH16.mod)

#GH17
normalTest(GH17$value,"sw")
GH17 <- subset(hemicellulose.carbon.cat.box,gene_name=="GH17");dim(GH17)
GH17.mod <- lmer(value~ Nitrogen +(1|block), data = GH17)
Anova(GH17.mod)#Nitrogen> no
summary(GH17.mod)
coef(GH17.mod)

#GH2
normalTest(GH2$value,"sw")
GH2 <- subset(hemicellulose.carbon.cat.box,gene_name=="GH2");dim(GH2)
GH2.mod <- lmer(value~ Nitrogen +(1|block), data = GH2)
Anova(GH2.mod)#Nitrogen>no 
summary(GH2.mod)
coef(GH2.mod)

#GH26
GH26 <- subset(hemicellulose.carbon.cat.box,gene_name=="GH26");dim(GH26)
GH26.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = GH26)
Anova(GH26.mod)#
summary(GH26.mod)
coef(GH26.mod)

#GH3
normalTest(GH3$value,"sw")
GH3 <- subset(hemicellulose.carbon.cat.box,gene_name=="GH3");dim(GH3)
GH3.mod <- lmer(value~ Nitrogen +(1|block), data = GH3)
Anova(GH3.mod)#Nitrogen>no
summary(GH3.mod)
coef(GH3.mod)

#GH30
GH30 <- subset(hemicellulose.carbon.cat.box,gene_name=="GH30");dim(GH30)
GH30.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = GH30)
Anova(GH30.mod)#
summary(GH30.mod)
coef(GH30.mod)

#GH31
GH31 <- subset(hemicellulose.carbon.cat.box,gene_name=="GH31");dim(GH31)
GH31.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = GH31)
Anova(GH31.mod)#
summary(GH31.mod)
coef(GH31.mod)

#GH38
normalTest(GH38$value,"sw")
GH38 <- subset(hemicellulose.carbon.cat.box,gene_name=="GH38");dim(GH38)
GH38.mod <- lmer(value~ Nitrogen +(1|block), data = GH38)
Anova(GH38.mod)#Nitrogen> no
summary(GH38.mod)
coef(GH38.mod)

#GH39
GH39 <- subset(hemicellulose.carbon.cat.box,gene_name=="GH39");dim(GH39)
GH39.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = GH39)
Anova(GH39.mod)#
summary(GH39.mod)
coef(GH39.mod)

#GH43
normalTest(GH43$value,"sw")#good
GH43 <- subset(hemicellulose.carbon.cat.box,gene_name=="GH43");dim(GH43)
GH43.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = GH43)
Anova(GH43.mod)#Cover:Nitrogen
summary(GH43.mod)
coef(GH43.mod)
nocover.GH43<- subset(GH43, Cover=="No_cover")
GH43.mod.1 <- lmer(value~ Nitrogen +(1|block), data = nocover.GH43)
nitrogen.GH43.1<-glht(GH43.mod.1,linfct=mcp(Nitrogen="Tukey"))
summary(nitrogen.GH43.1)#no sig
cover.GH43<- subset(GH43, Cover=="Vetch")
GH43.mod.2 <- lmer(value~ Nitrogen +(1|block), data = cover.GH43)
nitrogen.GH43.2<-glht(GH43.mod.2,linfct=mcp(Nitrogen="Tukey"))
summary(nitrogen.GH43.2)# no fer<fer

n0.GH43<- subset(GH43, Nitrogen=="N0")
GH43.mod.3 <- lmer(value~ Cover +(1|block), data = n0.GH43)
cover.GH43.3<-glht(GH43.mod.3,linfct=mcp(Cover="Tukey"))
summary(cover.GH43.3)#vetch<no 
n60.GH43<- subset(GH43,Nitrogen=="N60")
GH43.mod.4 <- lmer(value~ Cover +(1|block), data = n60.GH43)
cover.GH43.4<-glht(GH43.mod.4,linfct=mcp(Cover="Tukey"))
summary(cover.GH43.4)#vetch>no 

#GH44
GH44 <- subset(hemicellulose.carbon.cat.box,gene_name=="GH44");dim(GH44)
GH44.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = GH44)
Anova(GH44.mod)#
summary(GH44.mod)
coef(GH44.mod)

#GH47
normalTest(GH47$value,"sw")
GH47 <- subset(hemicellulose.carbon.cat.box,gene_name=="GH47");dim(GH47)
GH47.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = GH47)
Anova(GH47.mod)#Cover:Nitrogen 
summary(GH47.mod)
coef(GH47.mod)
nocover.GH47<- subset(GH47, Cover=="No_cover")
GH47.mod.1 <- lmer(value~ Nitrogen +(1|block), data = nocover.GH47)
nitrogen.GH47.1<-glht(GH47.mod.1,linfct=mcp(Nitrogen="Tukey"))
summary(nitrogen.GH47.1)#no fer < fer under no cover
cover.GH47<- subset(GH47, Cover=="Vetch")
GH47.mod.2 <- lmer(value~ Nitrogen +(1|block), data = cover.GH47)
nitrogen.GH47.2<-glht(GH47.mod.2,linfct=mcp(Nitrogen="Tukey"))
summary(nitrogen.GH47.2)# no fer< fer

n0.GH47<- subset(GH47, Nitrogen=="N0")
GH47.mod.3 <- lmer(value~ Cover +(1|block), data = n0.GH47)
cover.GH47.3<-glht(GH47.mod.3,linfct=mcp(Cover="Tukey"))
summary(cover.GH47.3)#no sig 
n60.GH47<- subset(GH47,Nitrogen=="N60")
GH47.mod.4 <- lmer(value~ Cover +(1|block), data = n60.GH47)
cover.GH47.4<-glht(GH47.mod.4,linfct=mcp(Cover="Tukey"))
summary(cover.GH47.4)#no sig
#GH5
normalTest(GH5$value,"sw")#good
GH5 <- subset(hemicellulose.carbon.cat.box,gene_name=="GH5");dim(GH5)
GH5.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = GH5)
Anova(GH5.mod)#Cover:Nitrogen
summary(GH5.mod)
coef(GH5.mod)
nocover.GH5<- subset(GH5, Cover=="No_cover")
GH5.mod.1 <- lmer(value~ Nitrogen +(1|block), data = nocover.GH5)
nitrogen.GH5.1<-glht(GH5.mod.1,linfct=mcp(Nitrogen="Tukey"))
summary(nitrogen.GH5.1)#no fer < fer under no cover
cover.GH5<- subset(GH5, Cover=="Vetch")
GH5.mod.2 <- lmer(value~ Nitrogen +(1|block), data = cover.GH5)
nitrogen.GH5.2<-glht(GH5.mod.2,linfct=mcp(Nitrogen="Tukey"))
summary(nitrogen.GH5.2)# no fer> fer

n0.GH5<- subset(GH5, Nitrogen=="N0")
GH5.mod.3 <- lmer(value~ Cover +(1|block), data = n0.GH5)
cover.GH5.3<-glht(GH5.mod.3,linfct=mcp(Cover="Tukey"))
summary(cover.GH5.3)#vetch>no 
n60.GH5<- subset(GH5,Nitrogen=="N60")
GH5.mod.4 <- lmer(value~ Cover +(1|block), data = n60.GH5)
cover.GH5.4<-glht(GH5.mod.4,linfct=mcp(Cover="Tukey"))
summary(cover.GH5.4)#vetch<no

#GH52
GH52 <- subset(hemicellulose.carbon.cat.box,gene_name=="GH52");dim(GH52)
GH52.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = GH52)
Anova(GH52.mod)#
summary(GH52.mod)
coef(GH52.mod)

#GH55
GH55 <- subset(hemicellulose.carbon.cat.box,gene_name=="GH55");dim(GH55)
GH55.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = GH55)
Anova(GH55.mod)#
summary(GH55.mod)
coef(GH55.mod)

#GH63
normalTest(GH63$value,"sw")#good
GH63 <- subset(hemicellulose.carbon.cat.box,gene_name=="GH63");dim(GH63)
GH63.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = GH63)
Anova(GH63.mod)#Cover:Nitrogen
summary(GH63.mod)
coef(GH63.mod)
nocover.GH63<- subset(GH63, Cover=="No_cover")
GH63.mod.1 <- lmer(value~ Nitrogen +(1|block), data = nocover.GH63)
nitrogen.GH63.1<-glht(GH63.mod.1,linfct=mcp(Nitrogen="Tukey"))
summary(nitrogen.GH63.1)#no fer < fer under no cover
cover.GH63<- subset(GH63, Cover=="Vetch")
GH63.mod.2 <- lmer(value~ Nitrogen +(1|block), data = cover.GH63)
nitrogen.GH63.2<-glht(GH63.mod.2,linfct=mcp(Nitrogen="Tukey"))
summary(nitrogen.GH63.2)# no fer> fer

n0.GH63<- subset(GH63, Nitrogen=="N0")
GH63.mod.3 <- lmer(value~ Cover +(1|block), data = n0.GH63)
cover.GH63.3<-glht(GH63.mod.3,linfct=mcp(Cover="Tukey"))
summary(cover.GH63.3)#vetch>no 
n60.GH63<- subset(GH63,Nitrogen=="N60")
GH63.mod.4 <- lmer(value~ Cover +(1|block), data = n60.GH63)
cover.GH63.4<-glht(GH63.mod.4,linfct=mcp(Cover="Tukey"))
summary(cover.GH63.4)#vetch<no
#GH64
GH64 <- subset(hemicellulose.carbon.cat.box,gene_name=="GH64");dim(GH64)
GH64.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = GH64)
Anova(GH64.mod)#
summary(GH64.mod)
coef(GH64.mod)

#GH67
GH67 <- subset(hemicellulose.carbon.cat.box,gene_name=="GH67");dim(GH67)
GH67.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = GH67)
Anova(GH67.mod)#
summary(GH67.mod)
coef(GH67.mod)

#GH74
GH74 <- subset(hemicellulose.carbon.cat.box,gene_name=="GH74");dim(GH74)
GH74.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = GH74)
Anova(GH74.mod)#
summary(GH74.mod)
coef(GH74.mod)

#GH76
GH76 <- subset(hemicellulose.carbon.cat.box,gene_name=="GH76");dim(GH76)
GH76.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = GH76)
Anova(GH76.mod)#
summary(GH76.mod)
coef(GH76.mod)

#GH8
normalTest(GH8$value,"sw")
GH8 <- subset(hemicellulose.carbon.cat.box,gene_name=="GH8");dim(GH8)
GH8.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = GH8)
Anova(GH8.mod)
summary(GH8.mod)
coef(GH8.mod)

#GH81
GH81 <- subset(hemicellulose.carbon.cat.box,gene_name=="GH81");dim(GH81)
GH81.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = GH81)
Anova(GH81.mod)#
summary(GH81.mod)
coef(GH81.mod)

#GH9
GH9 <- subset(hemicellulose.carbon.cat.box,gene_name=="GH9");dim(GH9)
GH9.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = GH9)
Anova(GH9.mod)#
summary(GH9.mod)
coef(GH9.mod)

#GH92
GH92 <- subset(hemicellulose.carbon.cat.box,gene_name=="GH92");dim(GH92)
GH92.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = GH92)
Anova(GH92.mod)#
summary(GH92.mod)
coef(GH92.mod)

#GH99
normalTest(GH99$value,"sw")#good
GH99 <- subset(hemicellulose.carbon.cat.box,gene_name=="GH99");dim(GH99)
GH99.mod <- lmer(value~ Nitrogen +(1|block), data = GH99)
Anova(GH99.mod)#Nitrogen > no
summary(GH99.mod)
coef(GH99.mod)



#4. Chitin
chitin.carbon.cat.box <- subset(carbon.cat.boxplot,group=="Chitin") %>% 
  mutate(gene_name=c(rep(c("GH5","GH116", "GH3" ,  "GH16"  ,"GH18" , "GH19" , "GH23" , "GH73" , "CE4" ,"GH20",
                           "GH84","AA10" ,"CBM12", "CBM2" , "CBM3"  ,"CBM5" , "CBM50" ,"CBM54","GT2" ,"GT23"),12)))
tiff('figures/box_c_degre_gene_treat.tiff', units="in", width=10, height=6, res=300)
ggplot(chitin.carbon.cat.box, aes(x =gene_name, y = value, fill=Treat))+
  geom_boxplot(size=0.3,width=0.5,alpha=1,notch = FALSE,na.rm=TRUE,outlier.shape=NA,position=position_dodge(preserve = "single"), color="black")+
  geom_jitter(aes(color=Treat),shape=16,size=0.5,position=position_jitterdodge(0.2),na.rm = TRUE)+
  labs(title="Chitin",x="Genes",y="Relative abundance",size=12)+
  scale_fill_manual(values=c("aquamarine4","royalblue4","grey","#e6b325"))+
  scale_color_manual(values = c("aquamarine4","royalblue4","grey","#e6b325"))+
  theme(legend.position="bottom",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_line(colour="grey",size = 0.1),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=12,family = "Arial",face="plain",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="plain",size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(size=12,family="Arial",face = "bold"),
        axis.title=element_text(size = 12,face="plain",family = "Arial",vjust = 1),
        legend.title = element_text(size=12,face="plain",family = "Arial"),
        legend.text = (element_text(size=12,family = "Arial")))
dev.off()
#####Chitin stats#####
#AA10
AA10 <- subset(chitin.carbon.cat.box,gene_name=="AA10");dim(AA10)
AA10.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = AA10)
Anova(AA10.mod)#
summary(AA10.mod)
coef(AA10.mod)
#CBM12
CBM12 <- subset(chitin.carbon.cat.box,gene_name=="CBM12");dim(CBM12)
CBM12.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = CBM12)
Anova(CBM12.mod)#
summary(CBM12.mod)
coef(CBM12.mod)
#CBM3
CBM3 <- subset(chitin.carbon.cat.box,gene_name=="CBM3");dim(CBM3)
CBM3.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = CBM3)
Anova(CBM3.mod)#
summary(CBM3.mod)
coef(CBM3.mod)
#CBM5
CBM5 <- subset(chitin.carbon.cat.box,gene_name=="CBM5");dim(CBM5)
CBM5.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = CBM5)
Anova(CBM5.mod)#Cover:Nitrogen
summary(CBM5.mod)
coef(CBM5.mod)
nocover.CBM5<- subset(CBM5, Cover=="No_cover")
CBM5.mod.1 <- lmer(value~ Nitrogen +(1|block), data = nocover.CBM5)
nitrogen.CBM5.1<-glht(CBM5.mod.1,linfct=mcp(Nitrogen="Tukey"))
summary(nitrogen.CBM5.1)#no fer < fer under no cover
cover.CBM5<- subset(CBM5, Cover=="Vetch")
CBM5.mod.2 <- lmer(value~ Nitrogen +(1|block), data = cover.CBM5)
nitrogen.CBM5.2<-glht(CBM5.mod.2,linfct=mcp(Nitrogen="Tukey"))
summary(nitrogen.CBM5.2)# no fer> fer

n0.CBM5<- subset(CBM5, Nitrogen=="N0")
CBM5.mod.3 <- lmer(value~ Cover +(1|block), data = n0.CBM5)
cover.CBM5.3<-glht(CBM5.mod.3,linfct=mcp(Cover="Tukey"))
summary(cover.CBM5.3)#vetch>no 
n60.CBM5<- subset(CBM5,Nitrogen=="N60")
CBM5.mod.4 <- lmer(value~ Cover +(1|block), data = n60.CBM5)
cover.CBM5.4<-glht(CBM5.mod.4,linfct=mcp(Cover="Tukey"))
summary(cover.CBM5.4)#vetch<no
#CBM50
CBM50 <- subset(chitin.carbon.cat.box,gene_name=="CBM50");dim(CBM50)
CBM50.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = CBM50)
Anova(CBM50.mod)#
summary(CBM50.mod)
coef(CBM50.mod)
#CBM54
CBM54 <- subset(chitin.carbon.cat.box,gene_name=="CBM54");dim(CBM54)
CBM54.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = CBM54)
Anova(CBM54.mod)#
summary(CBM54.mod)
coef(CBM54.mod)
#CE4
CE4 <- subset(chitin.carbon.cat.box,gene_name=="CE4");dim(CE4)
CE4.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = CE4)
Anova(CE4.mod)#Cover:Nitrogen
summary(CE4.mod)
coef(CE4.mod)
nocover.CE4<- subset(CE4, Cover=="No_cover")
CE4.mod.1 <- lmer(value~ Nitrogen +(1|block), data = nocover.CE4)
nitrogen.CE4.1<-glht(CE4.mod.1,linfct=mcp(Nitrogen="Tukey"))
summary(nitrogen.CE4.1)#no fer < fer under no cover
cover.CE4<- subset(CE4, Cover=="Vetch")
CE4.mod.2 <- lmer(value~ Nitrogen +(1|block), data = cover.CE4)
nitrogen.CE4.2<-glht(CE4.mod.2,linfct=mcp(Nitrogen="Tukey"))
summary(nitrogen.CE4.2)# no sig

n0.CE4<- subset(CE4, Nitrogen=="N0")
CE4.mod.3 <- lmer(value~ Cover +(1|block), data = n0.CE4)
cover.CE4.3<-glht(CE4.mod.3,linfct=mcp(Cover="Tukey"))
summary(cover.CE4.3)#vetch>no 
n60.CE4<- subset(CE4,Nitrogen=="N60")
CE4.mod.4 <- lmer(value~ Cover +(1|block), data = n60.CE4)
cover.CE4.4<-glht(CE4.mod.4,linfct=mcp(Cover="Tukey"))
summary(cover.CE4.4)#no sig


#GH116
GH116 <- subset(chitin.carbon.cat.box,gene_name=="GH116");dim(GH116)
GH116.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = GH116)
Anova(GH116.mod)#
summary(GH116.mod)
coef(GH116.mod)
#GH16
GH16 <- subset(chitin.carbon.cat.box,gene_name=="GH16");dim(GH16)
GH16.mod <- lmer(value~ Nitrogen +(1|block), data = GH16)
Anova(GH16.mod)#Nitrogen>no 
summary(GH16.mod)
coef(GH16.mod)
#GH18
GH18 <- subset(chitin.carbon.cat.box,gene_name=="GH18");dim(GH18)
GH18.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = GH18)
Anova(GH18.mod)#
summary(GH18.mod)
coef(GH18.mod)
#GH19
GH19 <- subset(chitin.carbon.cat.box,gene_name=="GH19");dim(GH19)
GH19.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = GH19)
Anova(GH19.mod)#
summary(GH19.mod)
coef(GH19.mod)
#GH20
GH20 <- subset(chitin.carbon.cat.box,gene_name=="GH20");dim(GH20)
GH20.mod <- lmer(value~ Nitrogen +(1|block), data = GH20)
Anova(GH20.mod)#Nitrogen>no
summary(GH20.mod)
coef(GH20.mod)
#GH23
GH23 <- subset(chitin.carbon.cat.box,gene_name=="GH23");dim(GH23)
GH23.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = GH23)
Anova(GH23.mod)#
summary(GH23.mod)
coef(GH23.mod)
#GH3
GH3 <- subset(chitin.carbon.cat.box,gene_name=="GH3");dim(GH3)
GH3.mod <- lmer(value~ Nitrogen +(1|block), data = GH3)
Anova(GH3.mod)#Nitrogen>no 
summary(GH3.mod)
coef(GH3.mod)
#GH5
GH5 <- subset(chitin.carbon.cat.box,gene_name=="GH5");dim(GH5)
GH5.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = GH5)
Anova(GH5.mod)#Cover:Nitrogen 
summary(GH5.mod)
coef(GH5.mod)
nocover.GH5<- subset(GH5, Cover=="No_cover")
GH5.mod.1 <- lmer(value~ Nitrogen +(1|block), data = nocover.GH5)
nitrogen.GH5.1<-glht(GH5.mod.1,linfct=mcp(Nitrogen="Tukey"))
summary(nitrogen.GH5.1)#no fer < fer under no cover
cover.GH5<- subset(GH5, Cover=="Vetch")
GH5.mod.2 <- lmer(value~ Nitrogen +(1|block), data = cover.GH5)
nitrogen.GH5.2<-glht(GH5.mod.2,linfct=mcp(Nitrogen="Tukey"))
summary(nitrogen.GH5.2)# no fer > fer under cover

n0.GH5<- subset(GH5, Nitrogen=="N0")
GH5.mod.3 <- lmer(value~ Cover +(1|block), data = n0.GH5)
cover.GH5.3<-glht(GH5.mod.3,linfct=mcp(Cover="Tukey"))
summary(cover.GH5.3)#vetch>no 
n60.GH5<- subset(GH5,Nitrogen=="N60")
GH5.mod.4 <- lmer(value~ Cover +(1|block), data = n60.GH5)
cover.GH5.4<-glht(GH5.mod.4,linfct=mcp(Cover="Tukey"))
summary(cover.GH5.4)#vetch<no 


#GH73
GH73 <- subset(chitin.carbon.cat.box,gene_name=="GH73");dim(GH73)
GH73.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = GH73)
Anova(GH73.mod)#Nitrogen
summary(GH73.mod)
coef(GH73.mod)
nocover.GH73<- subset(GH73, Cover=="No_cover")
GH73.mod.1 <- lmer(value~ Nitrogen +(1|block), data = nocover.GH73)
nitrogen.GH73.1<-glht(GH73.mod.1,linfct=mcp(Nitrogen="Tukey"))
summary(nitrogen.GH73.1)#no fer > fer under no cover
cover.GH73<- subset(GH73, Cover=="Vetch")
GH73.mod.2 <- lmer(value~ Nitrogen +(1|block), data = cover.GH73)
nitrogen.GH73.2<-glht(GH73.mod.2,linfct=mcp(Nitrogen="Tukey"))
summary(nitrogen.GH73.2)# not sig

n0.GH73<- subset(GH73, Nitrogen=="N0")
GH73.mod.3 <- lmer(value~ Cover +(1|block), data = n0.GH73)
cover.GH73.3<-glht(GH73.mod.3,linfct=mcp(Cover="Tukey"))
summary(cover.GH73.3)#no sig
n60.GH73<- subset(GH73,Nitrogen=="N60")
GH73.mod.4 <- lmer(value~ Cover +(1|block), data = n60.GH73)
cover.GH73.4<-glht(GH73.mod.4,linfct=mcp(Cover="Tukey"))
summary(cover.GH73.4)#no sig ALL 0


#GH84
GH84 <- subset(chitin.carbon.cat.box,gene_name=="GH84");dim(GH84)
GH84.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = GH84)
Anova(GH84.mod)
summary(GH84.mod)
coef(GH84.mod)
#GT2
GT2 <- subset(chitin.carbon.cat.box,gene_name=="GT2");dim(GT2)
GT2.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = GT2)
Anova(GT2.mod)
summary(GT2.mod)
coef(GT2.mod)
#GT23
GT23 <- subset(chitin.carbon.cat.box,gene_name=="GT23");dim(GT23)
GT23.mod <- lmer(value~ Nitrogen +(1|block), data = GT23)
Anova(GT23.mod)#no fer>fer 
summary(GT23.mod)
coef(GT23.mod)

#6. Pectin
pectin.carbon.cat.box <- subset(carbon.cat.boxplot,group=="Pectin") %>% 
  mutate(gene_name=c(rep(c("PL1","CE12","CE8","CE13"),12)))
tiff('figures/box_c_degre_gene_treat.tiff', units="in", width=10, height=6, res=300)
ggplot(pectin.carbon.cat.box, aes(x =gene_name, y = value, fill=Treat))+
  geom_boxplot(size=0.3,width=0.5,alpha=1,notch = FALSE,na.rm=TRUE,outlier.shape=NA,position=position_dodge(preserve = "single"), color="black")+
  geom_jitter(aes(color=Treat),shape=16,size=0.5,position=position_jitterdodge(0.2),na.rm = TRUE)+
  labs(title="Pectin",x="Gene ID",y="Relative abundance",size=12)+
  scale_fill_manual(values=c("aquamarine4","royalblue4","grey","#e6b325"))+
  scale_color_manual(values = c("aquamarine4","royalblue4","grey","#e6b325"))+
  theme(legend.position="bottom",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_line(colour="grey",size = 0.1),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=12,family = "Arial",face="plain",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="plain",size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(size=12,family="Arial",face = "bold"),
        axis.title=element_text(size = 12,face="plain",family = "Arial",vjust = 1),
        legend.title = element_text(size=12,face="plain",family = "Arial"),
        legend.text = (element_text(size=12,family = "Arial")))
  
dev.off()
#####Pectin stats#####
#CE12
CE12 <- subset(pectin.carbon.cat.box,gene_name=="CE12");dim(CE12)
CE12.mod <- lmer(value~ Cover+Nitrogen +(1|block), data = CE12)
Anova(CE12.mod)#cover and nitrogen
summary(CE12.mod)
coef(CE12.mod)
nocover.CE12<- subset(CE12, Cover=="No_cover")
CE12.mod.1 <- lmer(value~ Nitrogen +(1|block), data = nocover.CE12)
nitrogen.CE12.1<-glht(CE12.mod.1,linfct=mcp(Nitrogen="Tukey"))
summary(nitrogen.CE12.1)#no fer > fer under no cover
cover.CE12<- subset(CE12, Cover=="Vetch")
CE12.mod.2 <- lmer(value~ Nitrogen +(1|block), data = cover.CE12)
nitrogen.CE12.2<-glht(CE12.mod.2,linfct=mcp(Nitrogen="Tukey"))
summary(nitrogen.CE12.2)# not sig

n0.CE12<- subset(CE12, Nitrogen=="N0")
CE12.mod.3 <- lmer(value~ Cover +(1|block), data = n0.CE12)
cover.CE12.3<-glht(CE12.mod.3,linfct=mcp(Cover="Tukey"))
summary(cover.CE12.3)#no sig
n60.CE12<- subset(CE12,Nitrogen=="N60")
CE12.mod.4 <- lmer(value~ Cover +(1|block), data = n60.CE12)
cover.CE12.4<-glht(CE12.mod.4,linfct=mcp(Cover="Tukey"))
summary(cover.CE12.4)#no cover > vetch under N60

#CE13
CE13 <- subset(pectin.carbon.cat.box,gene_name=="CE13");dim(CE13)
CE13.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = CE13)
Anova(CE13.mod)#cover and nitrogen
summary(CE13.mod)
coef(CE13.mod)

#CE8
CE8 <- subset(pectin.carbon.cat.box,gene_name=="CE8");dim(CE8)
CE8.mod <- lmer(value~ Cover*Nitrogen +(1|block), data = CE8)
Anova(CE8.mod)#cover and nitrogen
summary(CE8.mod)
coef(CE8.mod)
nocover.CE8<- subset(CE8, Cover=="No_cover")
CE8.mod.1 <- lmer(value~ Nitrogen +(1|block), data = nocover.CE8)
nitrogen.CE8.1<-glht(CE8.mod.1,linfct=mcp(Nitrogen="Tukey"))
summary(nitrogen.CE8.1)#no fer > fer under no cover
cover.CE8<- subset(CE8, Cover=="Vetch")
CE8.mod.2 <- lmer(value~ Nitrogen +(1|block), data = cover.CE8)
nitrogen.CE8.2<-glht(CE8.mod.2,linfct=mcp(Nitrogen="Tukey"))
summary(nitrogen.CE8.2)# not sig

n0.CE8<- subset(CE8, Nitrogen=="N0")
CE8.mod.3 <- lmer(value~ Cover +(1|block), data = n0.CE8)
cover.CE8.3<-glht(CE8.mod.3,linfct=mcp(Cover="Tukey"))
summary(cover.CE8.3)#no cover>vetch under N0
n60.CE8<- subset(CE8,Nitrogen=="N60")
CE8.mod.4 <- lmer(value~ Cover +(1|block), data = n60.CE8)
cover.CE8.4<-glht(CE8.mod.4,linfct=mcp(Cover="Tukey"))
summary(cover.CE8.4)# no sig

#PL1
PL1 <- subset(pectin.carbon.cat.box,gene_name=="PL1");dim(PL1)
PL1.mod <- lmer(value~ Cover +(1|block), data = PL1)
Anova(PL1.mod)#no cover> cover
summary(PL1.mod)
coef(PL1.mod)


#7. aromatic
aromatic.carbon.cat.box.df <- subset(carbon.cat.boxplot,group=="Aromatic") %>% 
  mutate(gene_name=c(rep(c("badH","badI","bcrC","bcrB","bcrA","bcrD","dch","had","oah","dmpB","mhpE_Catechol","praC","mhpD_Catechol","catE",
                           "dmpC","pcaD","catB","catA","catC","pcaL_Catechol","cmtAb","cmtAc","hpaE","hpaD","hpaF","hpaG","ligC_Homoprotocatechuate","pht5",
                           "pht4","pht3","pht2","bbsE","bbsF","bbsH","xylB","tmoA","tmoB","tmoC","tmoD","tmoE","tmoF","hcaD","mhpE_cinnamate","mhpD_cinnamate",
                           "hcaC","mhpA","mhpB","mhpC","E1.1.1.90","lpdC","lpdD","pcaG","pcaH","pcaC","pcaB","pcaL_protocatechuate","ligA","ligB","ligK","ligC_protocatechuate","ligJ","ligI","galD"),12)))
degradation <- separate(as.data.frame(aromatic.carbon.cat.box$module),1,sep=",",into = c("Category","other_info"))# seperate by ,
#view(degradation)
aromatic.carbon.cat.box <- cbind(aromatic.carbon.cat.box.df,degradation) 
ggplot(aromatic.carbon.cat.box, aes(x =gene_name, y = value, fill=Treat))+
  geom_boxplot(size=0.3,width=0.5,alpha=1,notch = FALSE,na.rm=TRUE,outlier.shape=NA,position=position_dodge(preserve = "single"), color="black")+
  geom_jitter(aes(color=Treat),shape=16,size=0.5,position=position_jitterdodge(0.2),na.rm = TRUE)+
  labs(title="Aromatic",x="Gene ID",y="Relative abundance",size=12)+
  scale_fill_manual(values=c("aquamarine4","royalblue4","grey","#e6b325"))+
  scale_color_manual(values = c("aquamarine4","royalblue4","grey","#e6b325"))+
  theme(legend.position="bottom",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_line(colour="grey",size = 0.1),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=12,family = "Arial",face="plain",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="plain",size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(size=12,family="Arial",face = "bold"),
        axis.title=element_text(size = 12,face="plain",family = "Arial",vjust = 1),
        legend.title = element_text(size=12,face="plain",family = "Arial"),
        legend.text = (element_text(size=12,family = "Arial")))+
  facet_grid(~Category)
dev.off()

unique(aromatic.carbon.cat.box$Category)
Benzoate <- subset(aromatic.carbon.cat.box,Category=="Benzoate degradation" )
ggplot(Benzoate, aes(x =gene_name, y = value, fill=Treat))+
  geom_boxplot(size=0.3,width=0.5,alpha=1,notch = FALSE,na.rm=TRUE,outlier.shape=NA,position=position_dodge(preserve = "single"), color="black")+
  geom_jitter(aes(color=Treat),shape=16,size=0.5,position=position_jitterdodge(0.2),na.rm = TRUE)+
  labs(title="Aromatic",x="Genes",y="Relative abundance",size=12)+
  scale_fill_manual(values=c("aquamarine4","royalblue4","grey","#e6b325"))+
  scale_color_manual(values = c("aquamarine4","royalblue4","grey","#e6b325"))+
  theme(legend.position="bottom",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_line(colour="grey",size = 0.1),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=12,family = "Arial",face="plain",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="plain",size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(size=12,family="Arial",face = "bold"),
        axis.title=element_text(size = 12,face="plain",family = "Arial",vjust = 1),
        legend.title = element_text(size=12,face="plain",family = "Arial"),
        legend.text = (element_text(size=12,family = "Arial")))

#badH
badH <- subset(Benzoate,gene_name=="badH");dim(badH)
badH.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = badH)
Anova(badH.mod)
summary(badH.mod)
coef(badH.mod)
#badI
badI <- subset(Benzoate,gene_name=="badI");dim(badI)
badI.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = badI)
Anova(badI.mod)
summary(badI.mod)
coef(badI.mod)

Benzoyl <- subset(aromatic.carbon.cat.box,Category=="Benzoyl-CoA degradation" )
ggplot(Benzoyl, aes(x =gene_name, y = value, fill=Treat))+
  geom_boxplot(size=0.3,width=0.5,alpha=1,notch = FALSE,na.rm=TRUE,outlier.shape=NA,position=position_dodge(preserve = "single"), color="black")+
  geom_jitter(aes(color=Treat),shape=16,size=0.5,position=position_jitterdodge(0.2),na.rm = TRUE)+
  labs(title="Aromatic",x="Genes",y="Relative abundance",size=12)+
  scale_fill_manual(values=c("aquamarine4","royalblue4","grey","#e6b325"))+
  scale_color_manual(values = c("aquamarine4","royalblue4","grey","#e6b325"))+
  theme(legend.position="bottom",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_line(colour="grey",size = 0.1),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=12,family = "Arial",face="plain",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="plain",size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(size=12,family="Arial",face = "bold"),
        axis.title=element_text(size = 12,face="plain",family = "Arial",vjust = 1),
        legend.title = element_text(size=12,face="plain",family = "Arial"),
        legend.text = (element_text(size=12,family = "Arial")))
#bcrA
bcrA <- subset(Benzoyl,gene_name=="bcrA");dim(bcrA)
bcrA.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = bcrA)
Anova(bcrA.mod)
summary(bcrA.mod)
coef(bcrA.mod)
#bcrB
bcrB <- subset(Benzoyl,gene_name=="bcrB");dim(bcrB)
bcrB.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = bcrB)
Anova(bcrB.mod)
summary(bcrB.mod)
coef(bcrB.mod)
#bcrC
bcrC <- subset(Benzoyl,gene_name=="bcrC");dim(bcrC)
bcrC.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = bcrC)
Anova(bcrC.mod)
summary(bcrC.mod)
coef(bcrC.mod)
#bcrD
bcrD <- subset(Benzoyl,gene_name=="bcrD");dim(bcrD)
bcrD.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = bcrD)
Anova(bcrD.mod)
summary(bcrD.mod)
coef(bcrD.mod)
#dch
dch <- subset(Benzoyl,gene_name=="dch");dim(dch)
dch.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = dch)
Anova(dch.mod)
summary(dch.mod)
coef(dch.mod)
#had
had <- subset(Benzoyl,gene_name=="had");dim(had)
had.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = had)
Anova(had.mod)
summary(had.mod)
coef(had.mod)
#oah
oah <- subset(Benzoyl,gene_name=="oah");dim(oah)
oah.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = oah)
Anova(oah.mod)
summary(oah.mod)
coef(oah.mod)

unique(aromatic.carbon.cat.box$Category)
Catechol.meta <- subset(aromatic.carbon.cat.box,Category=="Catechol meta-cleavage" )
ggplot(Catechol.meta, aes(x =gene_name, y = value, fill=Treat))+
  geom_boxplot(size=0.3,width=0.5,alpha=1,notch = FALSE,na.rm=TRUE,outlier.shape=NA,position=position_dodge(preserve = "single"), color="black")+
  geom_jitter(aes(color=Treat),shape=16,size=0.5,position=position_jitterdodge(0.2),na.rm = TRUE)+
  labs(title="Aromatic",x="Genes",y="Relative abundance",size=12)+
  scale_fill_manual(values=c("aquamarine4","royalblue4","grey","#e6b325"))+
  scale_color_manual(values = c("aquamarine4","royalblue4","grey","#e6b325"))+
  theme(legend.position="bottom",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_line(colour="grey",size = 0.1),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=12,family = "Arial",face="plain",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="plain",size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(size=12,family="Arial",face = "bold"),
        axis.title=element_text(size = 12,face="plain",family = "Arial",vjust = 1),
        legend.title = element_text(size=12,face="plain",family = "Arial"),
        legend.text = (element_text(size=12,family = "Arial")))
#catE
catE <- subset(Catechol.meta,gene_name=="catE");dim(catE)
catE.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = catE)
Anova(catE.mod)#Cover:Nitrogen
summary(catE.mod)
coef(catE.mod)
#dmpB
dmpB <- subset(Catechol.meta,gene_name=="dmpB");dim(dmpB)
dmpB.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = dmpB)
Anova(dmpB.mod)
summary(dmpB.mod)
coef(dmpB.mod)
#dmpC
dmpC <- subset(Catechol.meta,gene_name=="dmpC");dim(dmpC)
dmpC.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = dmpC)
Anova(dmpC.mod)
summary(dmpC.mod)
coef(dmpC.mod)
#mnpD(Catechol)
mhpD.Catechol <- subset(Catechol.meta,gene_name=="mhpD_Catechol");dim(mhpD.Catechol)
mhpD.Catechol.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = mhpD.Catechol)
Anova(mhpD.Catechol.mod)#Cover 
summary(mhpD.Catechol.mod)
coef(mhpD.Catechol.mod)
#mnpE(Catechol)
mhpE.Catechol <- subset(Catechol.meta,gene_name=="mhpE_Catechol");dim(mhpE.Catechol)
mhpE.Catechol.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = mhpE.Catechol)
Anova(mhpE.Catechol.mod)
summary(mhpE.Catechol.mod)
coef(mhpE.Catechol.mod)
#praC
praC <- subset(Catechol.meta,gene_name=="praC");dim(praC)
praC.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = praC)
Anova(praC.mod)
summary(praC.mod)
coef(praC.mod)

unique(aromatic.carbon.cat.box$Category)
Catechol.ortho <- subset(aromatic.carbon.cat.box,Category=="Catechol ortho-cleavage" )
ggplot(Catechol.ortho, aes(x =gene_name, y = value, fill=Treat))+
  geom_boxplot(size=0.3,width=0.5,alpha=1,notch = FALSE,na.rm=TRUE,outlier.shape=NA,position=position_dodge(preserve = "single"), color="black")+
  geom_jitter(aes(color=Treat),shape=16,size=0.5,position=position_jitterdodge(0.2),na.rm = TRUE)+
  labs(title="Aromatic",x="Genes",y="Relative abundance",size=12)+
  scale_fill_manual(values=c("aquamarine4","royalblue4","grey","#e6b325"))+
  scale_color_manual(values = c("aquamarine4","royalblue4","grey","#e6b325"))+
  theme(legend.position="bottom",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_line(colour="grey",size = 0.1),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=12,family = "Arial",face="plain",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="plain",size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(size=12,family="Arial",face = "bold"),
        axis.title=element_text(size = 12,face="plain",family = "Arial",vjust = 1),
        legend.title = element_text(size=12,face="plain",family = "Arial"),
        legend.text = (element_text(size=12,family = "Arial")))

#catA
catA <- subset(Catechol.ortho,gene_name=="catA");dim(catA)
catA.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = catA)
Anova(catA.mod)
summary(catA.mod)
coef(catA.mod)

#catB
catB <- subset(Catechol.ortho,gene_name=="catB");dim(catB)
catB.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = catB)
Anova(catB.mod)
summary(catB.mod)
coef(catB.mod)

#catC
catC <- subset(Catechol.ortho,gene_name=="catC");dim(catC)
catC.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = catC)
Anova(catC.mod)#Cover
summary(catC.mod)
coef(catC.mod)

#pcaD
pcaD <- subset(Catechol.ortho,gene_name=="pcaD");dim(pcaD)
pcaD.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = pcaD)
Anova(pcaD.mod)#Nitrogen 
summary(pcaD.mod)
coef(pcaD.mod)

#pcaL_Catechol
pcaL_Catechol <- subset(Catechol.ortho,gene_name=="pcaL_Catechol");dim(pcaL_Catechol)
pcaL_Catechol.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = pcaL_Catechol)
Anova(pcaL_Catechol.mod)#Cover
summary(pcaL_Catechol.mod)
coef(pcaL_Catechol.mod)

unique(aromatic.carbon.cat.box$Category)
Cumate<- subset(aromatic.carbon.cat.box,Category=="Cumate degradation" )
ggplot(Cumate, aes(x =gene_name, y = value, fill=Treat))+
  geom_boxplot(size=0.3,width=0.5,alpha=1,notch = FALSE,na.rm=TRUE,outlier.shape=NA,position=position_dodge(preserve = "single"), color="black")+
  geom_jitter(aes(color=Treat),shape=16,size=0.5,position=position_jitterdodge(0.2),na.rm = TRUE)+
  labs(title="Aromatic",x="Genes",y="Relative abundance",size=12)+
  scale_fill_manual(values=c("aquamarine4","royalblue4","grey","#e6b325"))+
  scale_color_manual(values = c("aquamarine4","royalblue4","grey","#e6b325"))+
  theme(legend.position="bottom",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_line(colour="grey",size = 0.1),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=12,family = "Arial",face="plain",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="plain",size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(size=12,family="Arial",face = "bold"),
        axis.title=element_text(size = 12,face="plain",family = "Arial",vjust = 1),
        legend.title = element_text(size=12,face="plain",family = "Arial"),
        legend.text = (element_text(size=12,family = "Arial")))
#cmtAb
cmtAb <- subset(Cumate,gene_name=="cmtAb");dim(cmtAb)
cmtAb.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = cmtAb)
Anova(cmtAb.mod)#Nitrogen
summary(cmtAb.mod)
coef(cmtAb.mod)

#cmtAc
cmtAc <- subset(Cumate,gene_name=="cmtAc");dim(cmtAc)
cmtAc.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = cmtAc)
Anova(cmtAc.mod)#
summary(cmtAc.mod)
coef(cmtAc.mod)

unique(aromatic.carbon.cat.box$Category)
Homoprotocatechuate<- subset(aromatic.carbon.cat.box,Category=="Homoprotocatechuate degradation" )
ggplot(Homoprotocatechuate, aes(x =gene_name, y = value, fill=Treat))+
  geom_boxplot(size=0.3,width=0.5,alpha=1,notch = FALSE,na.rm=TRUE,outlier.shape=NA,position=position_dodge(preserve = "single"), color="black")+
  geom_jitter(aes(color=Treat),shape=16,size=0.5,position=position_jitterdodge(0.2),na.rm = TRUE)+
  labs(title="Aromatic",x="Genes",y="Relative abundance",size=12)+
  scale_fill_manual(values=c("aquamarine4","royalblue4","grey","#e6b325"))+
  scale_color_manual(values = c("aquamarine4","royalblue4","grey","#e6b325"))+
  theme(legend.position="bottom",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_line(colour="grey",size = 0.1),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=12,family = "Arial",face="plain",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="plain",size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(size=12,family="Arial",face = "bold"),
        axis.title=element_text(size = 12,face="plain",family = "Arial",vjust = 1),
        legend.title = element_text(size=12,face="plain",family = "Arial"),
        legend.text = (element_text(size=12,family = "Arial")))
#hpaD
hpaD <- subset(Homoprotocatechuate,gene_name=="hpaD");dim(hpaD)
hpaD.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = hpaD)
Anova(hpaD.mod)#
summary(hpaD.mod)
coef(hpaD.mod)
#hpaE
hpaE <- subset(Homoprotocatechuate,gene_name=="hpaE");dim(hpaE)
hpaE.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = hpaE)
Anova(hpaE.mod)#
summary(hpaE.mod)
coef(hpaE.mod)
#hpaF
hpaF <- subset(Homoprotocatechuate,gene_name=="hpaF");dim(hpaF)
hpaF.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = hpaF)
Anova(hpaF.mod)#
summary(hpaF.mod)
coef(hpaF.mod)
#hpaG
hpaG <- subset(Homoprotocatechuate,gene_name=="hpaG");dim(hpaG)
hpaG.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = hpaG)
Anova(hpaG.mod)#Cover:Nitrogen
summary(hpaG.mod)
coef(hpaG.mod)
#ligC_Homoprotocatechuate
ligC_Homoprotocatechuate <- subset(Homoprotocatechuate,gene_name=="ligC_Homoprotocatechuate");dim(ligC_Homoprotocatechuate)
ligC_Homoprotocatechuate.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = ligC_Homoprotocatechuate)
Anova(ligC_Homoprotocatechuate.mod)#
summary(ligC_Homoprotocatechuate.mod)
coef(ligC_Homoprotocatechuate.mod)

unique(aromatic.carbon.cat.box$Category)
Phthalate<- subset(aromatic.carbon.cat.box,Category=="Phthalate degradation" )
ggplot(Phthalate, aes(x =gene_name, y = value, fill=Treat))+
  geom_boxplot(size=0.3,width=0.5,alpha=1,notch = FALSE,na.rm=TRUE,outlier.shape=NA,position=position_dodge(preserve = "single"), color="black")+
  geom_jitter(aes(color=Treat),shape=16,size=0.5,position=position_jitterdodge(0.2),na.rm = TRUE)+
  labs(title="Aromatic",x="Genes",y="Relative abundance",size=12)+
  scale_fill_manual(values=c("aquamarine4","royalblue4","grey","#e6b325"))+
  scale_color_manual(values = c("aquamarine4","royalblue4","grey","#e6b325"))+
  theme(legend.position="bottom",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_line(colour="grey",size = 0.1),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=12,family = "Arial",face="plain",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="plain",size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(size=12,family="Arial",face = "bold"),
        axis.title=element_text(size = 12,face="plain",family = "Arial",vjust = 1),
        legend.title = element_text(size=12,face="plain",family = "Arial"),
        legend.text = (element_text(size=12,family = "Arial")))
#pht2
pht2 <- subset(Phthalate,gene_name=="pht2");dim(pht2)
pht2.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = pht2)
Anova(pht2.mod)
summary(pht2.mod)
coef(pht2.mod)

#pht3
pht3 <- subset(Phthalate,gene_name=="pht3");dim(pht3)
pht3.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = pht3)
Anova(pht3.mod)
summary(pht3.mod)
coef(pht3.mod)

#pht4
pht4 <- subset(Phthalate,gene_name=="pht4");dim(pht4)
pht4.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = pht4)
Anova(pht4.mod)#Nitrogen
summary(pht4.mod)
coef(pht4.mod)

#pht5
pht5 <- subset(Phthalate,gene_name=="pht5");dim(pht5)
pht5.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = pht5)
Anova(pht5.mod)#Nitrogen 
summary(pht5.mod)
coef(pht5.mod)

unique(aromatic.carbon.cat.box$Category)
Toluene<- subset(aromatic.carbon.cat.box,Category=="Toluene degradation" )
ggplot(Toluene, aes(x =gene_name, y = value, fill=Treat))+
  geom_boxplot(size=0.3,width=0.5,alpha=1,notch = FALSE,na.rm=TRUE,outlier.shape=NA,position=position_dodge(preserve = "single"), color="black")+
  geom_jitter(aes(color=Treat),shape=16,size=0.5,position=position_jitterdodge(0.2),na.rm = TRUE)+
  labs(title="Aromatic",x="Genes",y="Relative abundance",size=12)+
  scale_fill_manual(values=c("aquamarine4","royalblue4","grey","#e6b325"))+
  scale_color_manual(values = c("aquamarine4","royalblue4","grey","#e6b325"))+
  theme(legend.position="bottom",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_line(colour="grey",size = 0.1),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=12,family = "Arial",face="plain",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="plain",size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(size=12,family="Arial",face = "bold"),
        axis.title=element_text(size = 12,face="plain",family = "Arial",vjust = 1),
        legend.title = element_text(size=12,face="plain",family = "Arial"),
        legend.text = (element_text(size=12,family = "Arial")))
#bbsE
bbsE <- subset(Toluene,gene_name=="bbsE");dim(bbsE)
bbsE.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = bbsE)
Anova(bbsE.mod) 
summary(bbsE.mod)
coef(bbsE.mod)
#bbsF
bbsF <- subset(Toluene,gene_name=="bbsF");dim(bbsF)
bbsF.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = bbsF)
Anova(bbsF.mod) #Nitrogen
summary(bbsF.mod)
coef(bbsF.mod)
#bbsH
bbsH <- subset(Toluene,gene_name=="bbsH");dim(bbsH)
bbsH.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = bbsH)
Anova(bbsH.mod) 
summary(bbsH.mod)
coef(bbsH.mod)
#tmoA
tmoA <- subset(Toluene,gene_name=="tmoA");dim(tmoA)
tmoA.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = tmoA)
Anova(tmoA.mod) #Cover 
summary(tmoA.mod)
coef(tmoA.mod)
#tmoB
tmoB <- subset(Toluene,gene_name=="tmoB");dim(tmoB)
tmoB.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = tmoB)
Anova(tmoB.mod) #Cover
summary(tmoB.mod)
coef(tmoB.mod)
#tmoC
tmoC <- subset(Toluene,gene_name=="tmoC");dim(tmoC)
tmoC.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = tmoC)
Anova(tmoC.mod) 
summary(tmoC.mod)
coef(tmoC.mod)
#tmoD
tmoD <- subset(Toluene,gene_name=="tmoD");dim(tmoD)
tmoD.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = tmoD)
Anova(tmoD.mod) 
summary(tmoD.mod)
coef(tmoD.mod)
#tmoE
tmoE <- subset(Toluene,gene_name=="tmoE");dim(tmoE)
tmoE.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = tmoE)
Anova(tmoE.mod) 
summary(tmoE.mod)
coef(tmoE.mod)
#tmoF
tmoF <- subset(Toluene,gene_name=="tmoF");dim(tmoF)
tmoF.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = tmoF)
Anova(tmoF.mod) 
summary(tmoF.mod)
coef(tmoF.mod)
#xylB
xylB <- subset(Toluene,gene_name=="xylB");dim(xylB)
xylB.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = xylB)
Anova(xylB.mod) #Cover 
summary(xylB.mod)
coef(xylB.mod)

unique(aromatic.carbon.cat.box$Category)
cinnama<- subset(aromatic.carbon.cat.box,Category=="Trans-cinnamate degradation" )
ggplot(cinnama, aes(x =gene_name, y = value, fill=Treat))+
  geom_boxplot(size=0.3,width=0.5,alpha=1,notch = FALSE,na.rm=TRUE,outlier.shape=NA,position=position_dodge(preserve = "single"), color="black")+
  geom_jitter(aes(color=Treat),shape=16,size=0.5,position=position_jitterdodge(0.2),na.rm = TRUE)+
  labs(title="Aromatic",x="Genes",y="Relative abundance",size=12)+
  scale_fill_manual(values=c("aquamarine4","royalblue4","grey","#e6b325"))+
  scale_color_manual(values = c("aquamarine4","royalblue4","grey","#e6b325"))+
  theme(legend.position="bottom",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_line(colour="grey",size = 0.1),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=12,family = "Arial",face="plain",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="plain",size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(size=12,family="Arial",face = "bold"),
        axis.title=element_text(size = 12,face="plain",family = "Arial",vjust = 1),
        legend.title = element_text(size=12,face="plain",family = "Arial"),
        legend.text = (element_text(size=12,family = "Arial")))
#hcaC
hcaC <- subset(cinnama,gene_name=="hcaC");dim(hcaC)
hcaC.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = hcaC)
Anova(hcaC.mod)  
summary(hcaC.mod)
coef(hcaC.mod)
#hcaD
hcaD <- subset(cinnama,gene_name=="hcaD");dim(hcaD)
hcaD.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = hcaD)
Anova(hcaD.mod)  
summary(hcaD.mod)
coef(hcaD.mod)
#mhpA
mhpA <- subset(cinnama,gene_name=="mhpA");dim(mhpA)
mhpA.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = mhpA)
Anova(mhpA.mod)
summary(mhpA.mod)
coef(mhpA.mod)
#mhpB
mhpB <- subset(cinnama,gene_name=="mhpB");dim(mhpB)
mhpB.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = mhpB)
Anova(mhpB.mod)  
summary(mhpB.mod)
coef(mhpB.mod)
#mhpC
mhpC <- subset(cinnama,gene_name=="mhpC");dim(mhpC)
mhpC.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = mhpC)
Anova(mhpC.mod) 
summary(mhpC.mod)
coef(mhpC.mod)
#mhpD_Trans-cinnamate
mhpD_cinnamate <- subset(cinnama,gene_name=="mhpD_cinnamate");dim(mhpD_cinnamate)
mhpD_cinnamate.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = mhpD_cinnamate)
Anova(mhpD_cinnamate.mod)#Cover  
summary(mhpD_cinnamate.mod)
coef(mhpD_cinnamate.mod)
#mhpEcinnamate
mhpE_cinnamate <- subset(cinnama,gene_name=="mhpE_cinnamate");dim(mhpE_cinnamate)
mhpE_cinnamate.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = mhpE_cinnamate)
Anova(mhpE_cinnamate.mod)  
summary(mhpE_cinnamate.mod)
coef(mhpE_cinnamate.mod)

unique(aromatic.carbon.cat.box$Category)
Xylene<- subset(aromatic.carbon.cat.box,Category=="Xylene degradation" )
ggplot(Xylene, aes(x =gene_name, y = value, fill=Treat))+
  geom_boxplot(size=0.3,width=0.5,alpha=1,notch = FALSE,na.rm=TRUE,outlier.shape=NA,position=position_dodge(preserve = "single"), color="black")+
  geom_jitter(aes(color=Treat),shape=16,size=0.5,position=position_jitterdodge(0.2),na.rm = TRUE)+
  labs(title="Aromatic",x="Genes",y="Relative abundance",size=12)+
  scale_fill_manual(values=c("aquamarine4","royalblue4","grey","#e6b325"))+
  scale_color_manual(values = c("aquamarine4","royalblue4","grey","#e6b325"))+
  theme(legend.position="bottom",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_line(colour="grey",size = 0.1),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=12,family = "Arial",face="plain",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="plain",size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(size=12,family="Arial",face = "bold"),
        axis.title=element_text(size = 12,face="plain",family = "Arial",vjust = 1),
        legend.title = element_text(size=12,face="plain",family = "Arial"),
        legend.text = (element_text(size=12,family = "Arial")))
#E1.1.1.90
E1.1.1.90 <- subset(Xylene,gene_name=="E1.1.1.90");dim(E1.1.1.90)
E1.1.1.90.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = E1.1.1.90)
Anova(E1.1.1.90.mod)#Cover  
summary(E1.1.1.90.mod)
coef(E1.1.1.90.mod)

unique(aromatic.carbon.cat.box$Category)
gallate<- subset(aromatic.carbon.cat.box,Category=="anaerobic gallate degradation" )
ggplot(gallate, aes(x =gene_name, y = value, fill=Treat))+
  geom_boxplot(size=0.3,width=0.5,alpha=1,notch = FALSE,na.rm=TRUE,outlier.shape=NA,position=position_dodge(preserve = "single"), color="black")+
  geom_jitter(aes(color=Treat),shape=16,size=0.5,position=position_jitterdodge(0.2),na.rm = TRUE)+
  labs(title="Aromatic",x="Genes",y="Relative abundance",size=12)+
  scale_fill_manual(values=c("aquamarine4","royalblue4","grey","#e6b325"))+
  scale_color_manual(values = c("aquamarine4","royalblue4","grey","#e6b325"))+
  theme(legend.position="bottom",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_line(colour="grey",size = 0.1),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=12,family = "Arial",face="plain",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="plain",size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(size=12,family="Arial",face = "bold"),
        axis.title=element_text(size = 12,face="plain",family = "Arial",vjust = 1),
        legend.title = element_text(size=12,face="plain",family = "Arial"),
        legend.text = (element_text(size=12,family = "Arial")))

#lpdC
lpdC <- subset(gallate,gene_name=="lpdC");dim(lpdC)
lpdC.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = lpdC)
Anova(lpdC.mod)  
summary(lpdC.mod)
coef(lpdC.mod)
#lpdD
lpdD <- subset(gallate,gene_name=="lpdD");dim(lpdD)
lpdD.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = lpdD)
Anova(lpdD.mod)  
summary(lpdD.mod)
coef(lpdD.mod)

unique(aromatic.carbon.cat.box$Category)
protocatechuate<- subset(aromatic.carbon.cat.box,Category=="protocatechuate degradation" )
ggplot(protocatechuate, aes(x =gene_name, y = value, fill=Treat))+
  geom_boxplot(size=0.3,width=0.5,alpha=1,notch = FALSE,na.rm=TRUE,outlier.shape=NA,position=position_dodge(preserve = "single"), color="black")+
  geom_jitter(aes(color=Treat),shape=16,size=0.5,position=position_jitterdodge(0.2),na.rm = TRUE)+
  labs(title="Aromatic",x="Genes",y="Relative abundance",size=12)+
  scale_fill_manual(values=c("aquamarine4","royalblue4","grey","#e6b325"))+
  scale_color_manual(values = c("aquamarine4","royalblue4","grey","#e6b325"))+
  theme(legend.position="bottom",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_line(colour="grey",size = 0.1),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=12,family = "Arial",face="plain",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="plain",size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(size=12,family="Arial",face = "bold"),
        axis.title=element_text(size = 12,face="plain",family = "Arial",vjust = 1),
        legend.title = element_text(size=12,face="plain",family = "Arial"),
        legend.text = (element_text(size=12,family = "Arial")))

#galD
galD <- subset(protocatechuate,gene_name=="galD");dim(galD)
galD.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = galD)
Anova(galD.mod)  
summary(galD.mod)
coef(galD.mod)
#ligA
ligA <- subset(protocatechuate,gene_name=="ligA");dim(ligA)
ligA.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = ligA)
Anova(ligA.mod)#Nitrogen  
summary(ligA.mod)
coef(ligA.mod)
#ligB
ligB <- subset(protocatechuate,gene_name=="ligB");dim(ligB)
ligB.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = ligB)
Anova(ligB.mod)  
summary(ligB.mod)
coef(ligB.mod)
#ligC_protocatechuate
ligC_protocatechuate <- subset(protocatechuate,gene_name=="ligC_protocatechuate");dim(ligC_protocatechuate)
ligC_protocatechuate.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = ligC_protocatechuate)
Anova(ligC_protocatechuate.mod)  
summary(ligC_protocatechuate.mod)
coef(ligC_protocatechuate.mod)
#ligI
ligI <- subset(protocatechuate,gene_name=="ligI");dim(ligI)
ligI.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = ligI)
Anova(ligI.mod) #Cover:Nitrogen 
summary(ligI.mod)
coef(ligI.mod)
#ligJ
ligJ <- subset(protocatechuate,gene_name=="ligJ");dim(ligJ)
ligJ.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = ligJ)
Anova(ligJ.mod) #Nitrogen 
summary(ligJ.mod)
coef(ligJ.mod)
#ligK
ligK <- subset(protocatechuate,gene_name=="ligK");dim(ligK)
ligK.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = ligK)
Anova(ligK.mod) #Nitrogen  
summary(ligK.mod)
coef(ligK.mod)
#pcaB
pcaB <- subset(protocatechuate,gene_name=="pcaB");dim(pcaB)
pcaB.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = pcaB)
Anova(pcaB.mod) #Cover 
summary(pcaB.mod)
coef(pcaB.mod)
#pcaC
pcaC <- subset(protocatechuate,gene_name=="pcaC");dim(pcaC)
pcaC.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = pcaC)
Anova(pcaC.mod)#Cover  Nitrogen
summary(pcaC.mod)
coef(pcaC.mod)
#pcaG
pcaG <- subset(protocatechuate,gene_name=="pcaG");dim(pcaG)
pcaG.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = pcaG)
Anova(pcaG.mod)  
summary(pcaG.mod)
coef(pcaG.mod)
#pcaH
pcaH <- subset(protocatechuate,gene_name=="pcaH");dim(pcaH)
pcaH.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = pcaH)
Anova(pcaH.mod)  
summary(pcaH.mod)
coef(pcaH.mod)
#pcaL_protocatechuate
pcaL_protocatechuate <- subset(protocatechuate,gene_name=="pcaL_protocatechuate");dim(pcaL_protocatechuate)
pcaL_protocatechuate.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = pcaL_protocatechuate)
Anova(pcaL_protocatechuate.mod) #Cover 
summary(pcaL_protocatechuate.mod)
coef(pcaL_protocatechuate.mod)


#8.lignin
lignin.carbon.cat.box <- subset(carbon.cat.boxplot,group=="Lignin") %>% 
  mutate(gene_name=c(rep(c("Mnp","vanA","ligXa","ligY","ligW","ligM"),12)))
ggplot(lignin.carbon.cat.box, aes(x =gene_name, y = value, fill=Treat))+
  geom_boxplot(size=0.3,width=0.5,alpha=1,notch = FALSE,na.rm=TRUE,outlier.shape=NA,position=position_dodge(preserve = "single"), color="black")+
  geom_jitter(aes(color=Treat),shape=16,size=0.5,position=position_jitterdodge(0.2),na.rm = TRUE)+
  labs(title="Lignin",x="Genes",y="Relative abundance",size=12)+
  scale_fill_manual(values=c("aquamarine4","royalblue4","grey","#e6b325"))+
  scale_color_manual(values = c("aquamarine4","royalblue4","grey","#e6b325"))+
  theme(legend.position="bottom",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_line(colour="grey",size = 0.1),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=12,family = "Arial",face="plain",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="plain",size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(size=12,family="Arial",face = "bold"),
        axis.title=element_text(size = 12,face="plain",family = "Arial",vjust = 1),
        legend.title = element_text(size=12,face="plain",family = "Arial"),
        legend.text = (element_text(size=12,family = "Arial")))
dev.off()
#####lignin stats#####
#ligM
ligM <- subset(lignin.carbon.cat.box,gene_name=="ligM");dim(ligM)
ligM.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = ligM)
Anova(ligM.mod)
summary(ligM.mod)
coef(ligM.mod)
#ligW
ligW <- subset(lignin.carbon.cat.box,gene_name=="ligW");dim(ligW)
ligW.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = ligW)
Anova(ligW.mod)
summary(ligW.mod)
coef(ligW.mod)
#ligXa
ligXa <- subset(lignin.carbon.cat.box,gene_name=="ligXa");dim(ligXa)
ligXa.mod <- lmer(value~ Nitrogen +(1|block), data = ligXa)
Anova(ligXa.mod)#nitrogen
summary(ligXa.mod)
coef(ligXa.mod)
nitrogen.ligXa<-glht(ligXa.mod,linfct=mcp(Nitrogen="Tukey"))
summary(nitrogen.ligXa)#no fer>fer
#ligY
ligY <- subset(lignin.carbon.cat.box,gene_name=="ligY");dim(ligY)
ligY.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = ligY)
Anova(ligY.mod)
summary(ligY.mod)
coef(ligY.mod)
#Mnp
Mnp <- subset(lignin.carbon.cat.box,gene_name=="Mnp");dim(Mnp)
Mnp.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = Mnp)
Anova(Mnp.mod)
summary(Mnp.mod)
coef(Mnp.mod)
#vanA
vanA <- subset(lignin.carbon.cat.box,gene_name=="vanA");dim(vanA)
vanA.mod <- lmer(value~ Cover* Nitrogen +(1|block), data = vanA)
Anova(vanA.mod)
summary(vanA.mod)
coef(vanA.mod)

#9 put all category together and see if there have a overall variation 
colnames(carbon.cat.boxplot)
carbon.cat.boxplot.group <- carbon.cat.boxplot%>% group_by(group,Sample) %>% summarise(sum.value=sum(value))
carbon.cat.boxplot.group.df <- inner_join(sample.data,carbon.cat.boxplot.group,by=c("Sample"="Sample"))
view(carbon.cat.boxplot.group.df)
ggplot(carbon.cat.boxplot.group.df, aes(x =group, y = sum.value, fill=Treat))+
  geom_boxplot(size=0.3,width=0.5,alpha=1,notch = FALSE,na.rm=TRUE,outlier.shape=NA,position=position_dodge(preserve = "single"), color="black")+
  geom_jitter(aes(color=Treat),shape=16,size=0.5,position=position_jitterdodge(0.2),na.rm = TRUE)+
  labs(title="",x="Gene ID",y="Relative abundance",size=12)+
  scale_fill_manual(values=c("aquamarine4","royalblue4","grey","#e6b325"))+
  scale_color_manual(values = c("aquamarine4","royalblue4","grey","#e6b325"))+
  theme(legend.position="bottom",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_line(colour="grey",size = 0.1),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=12,family = "Arial",face="plain",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="plain",size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(size=12,family="Arial",face = "bold"),
        axis.title=element_text(size = 12,face="plain",family = "Arial",vjust = 1),
        legend.title = element_text(size=12,face="plain",family = "Arial"),
        legend.text = (element_text(size=12,family = "Arial")))


df_carbon=carbon.cat.boxplot.group.df %>% group_by(Treat, group) %>% mutate(sem_carbon=sd(sum.value)/sqrt(n())) #stdev
mean_carbon = df_carbon %>% group_by(Treat, group) %>% mutate(avg_carbon=sum(sum.value)/(n()))#combine mean and stdev in "mean_shannon"
group.order <- c("Starch","Hemicellulose","Cellulose","Chitin","Pectin","Aromatic","Lignin")#make sure the order is the same to the bac.mag.tree
mean_carbon$group <- factor(mean_carbon$group, levels = as.factor(group.order))

##shannon_tillage
#mutiple barplot
tiff('figures/carbon_gene_bar.tiff', units="in", width=12, height=5, res=300)
carbon_gene_bar <- ggplot(mean_carbon, aes(x = group, y = avg_carbon, fill = Treat))+
  geom_col(alpha = 0.5, position = 'dodge')+
#  ylim(0,8)+
  scale_fill_manual(values = c("#8a0101","#254870","#a0a3a2","#7a4c0f"),name="Treatment")+
  labs(title="",x="Treatment",y="Relative abundance")+
  geom_errorbar(aes(ymin =mean_carbon$avg_carbon - mean_carbon$sem_carbon, ymax = mean_carbon$avg_carbon + mean_carbon$sem_carbon),width = 0.2, colour = "black", position = position_dodge(.9))+
#facet_wrap(~Tillage)+
#  theme_bw()+
  theme(legend.position="right",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title=(element_text(size=5,family = "Arial",face="plain",hjust = 0.5,vjust = 3)),
        text=element_text(family = "Arial",face="plain"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text=element_text(size=15,family="Arial"),
        axis.title=element_text(size = 15,face="plain",family = "Arial"),
        legend.title = element_text(size=15,face="plain",family = "Arial"),
        legend.text = (element_text(size=15,family = "Arial")))+
  annotate(geom="text", x=1.89, y=0.08, size=8,label="*",color="black", fontface="bold")+
  annotate(geom="text", x=2.33, y=0.08, size=8,label="*",color="black", fontface="bold")+
  annotate(geom="text", x=3.89, y=0.125, size=8,label="*",color="black", fontface="bold")+
  annotate(geom="text", x=4.33, y=0.125, size=8,label="*",color="black", fontface="bold")+
  annotate(geom="text", x=4.67, y=0.01, size=8,label="*",color="black", fontface="bold")+
  annotate(geom="text", x=6.67, y=0.01, size=8,label="*",color="black", fontface="bold")+
  annotate(geom="text", x=7.1, y=0.01, size=8,label="*",color="black", fontface="bold")

tiff('figures/carbon_gene_bar.tiff', units="in", width=12, height=5, res=300)
carbon_gene_bar <- ggplot(mean_carbon, aes(x = group, y = avg_carbon, fill = Treat))+
  geom_col(alpha = 0.5, position = 'dodge')+
  #  ylim(0,8)+
  scale_fill_manual(values = c("#8a0101","#254870","#a0a3a2","#7a4c0f"),name="Treatment")+
  labs(title="",x="Treatment",y="Relative abundance")+
  geom_errorbar(aes(ymin =mean_carbon$avg_carbon - mean_carbon$sem_carbon, ymax = mean_carbon$avg_carbon + mean_carbon$sem_carbon),width = 0.2, colour = "black", position = position_dodge(.9))+
  #facet_wrap(~Tillage)+
  #  theme_bw()+
  theme(legend.position="right",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title=(element_text(size=5,family = "Arial",face="plain",hjust = 0.5,vjust = 3)),
        text=element_text(family = "Arial",face="plain"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text=element_text(size=15,family="Arial"),
        axis.title=element_text(size = 15,face="plain",family = "Arial"),
        legend.title = element_text(size=15,face="plain",family = "Arial"),
        legend.text = (element_text(size=15,family = "Arial")))+
  annotate(geom="text", x=1.89, y=0.08, size=8,label="*",color="black", fontface="bold")+
  annotate(geom="text", x=2.33, y=0.08, size=8,label="*",color="black", fontface="bold")+
  annotate(geom="text", x=3.89, y=0.125, size=8,label="*",color="black", fontface="bold")+
  annotate(geom="text", x=4.33, y=0.125, size=8,label="*",color="black", fontface="bold")+
  annotate(geom="text", x=4.67, y=0.01, size=8,label="*",color="black", fontface="bold")+
  annotate(geom="text", x=6.67, y=0.01, size=8,label="*",color="black", fontface="bold")+
  annotate(geom="text", x=7.1, y=0.01, size=8,label="*",color="black", fontface="bold")  
  
dev.off()


tiff('figures/hemi_gene_bar.tiff', units="in", width=12, height=5, res=300)
mean_carbon_rev <- mutate(as.data.frame(mean_carbon), t2=c(rep("NCN0",21),rep("NCN60",21),rep("VN0",21),rep("VN60",21)))
hemi_data <- subset(mean_carbon_rev,group=="Hemicellulose")
hemi_gene_bar <- ggplot(hemi_data, aes(x = t2, y = avg_carbon, fill = Treat))+
  geom_col(alpha = 0.5,width = 0.7, position = 'dodge')+
  #  ylim(0,8)+
  scale_fill_manual(values = c("#8a0101","#254870","#a0a3a2","#7a4c0f"),name="Treatment")+
  labs(title="Hemicellulose",x="Treatment",y="Relative abundance")+
  geom_errorbar(aes(ymin =hemi_data$avg_carbon - hemi_data$sem_carbon, ymax = hemi_data$avg_carbon + hemi_data$sem_carbon),width = 0.2, colour = "black", position = position_dodge(.9))+
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
  annotate(geom="text", x=2, y=0.08, size=8,label="*",color="black", fontface="bold")+
  annotate(geom="text", x=4, y=0.08, size=8,label="*",color="black", fontface="bold")

dev.off()

chitin_data <- subset(mean_carbon_rev,group=="Chitin")
chitin_gene_bar <- ggplot(chitin_data, aes(x = t2, y = avg_carbon, fill = Treat))+
  geom_col(alpha = 0.5,width = 0.7, position = 'dodge')+
  #  ylim(0,8)+
  scale_fill_manual(values = c("#8a0101","#254870","#a0a3a2","#7a4c0f"),name="Treatment")+
  labs(title="Chitin",x="Treatment",y="")+
  geom_errorbar(aes(ymin =chitin_data$avg_carbon - chitin_data$sem_carbon, ymax = chitin_data$avg_carbon + chitin_data$sem_carbon),width = 0.2, colour = "black", position = position_dodge(.9))+
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
  annotate(geom="text", x=2, y=0.125, size=8,label="*",color="black", fontface="bold")+
  annotate(geom="text", x=4, y=0.125, size=8,label="*",color="black", fontface="bold")

pectin_data <- subset(mean_carbon_rev,group=="Pectin")
pectin_gene_bar <- ggplot(pectin_data , aes(x = t2, y = avg_carbon, fill = Treat))+
  geom_col(alpha = 0.5,width = 0.7, position = 'dodge')+
  #  ylim(0,8)+
  scale_fill_manual(values = c("#8a0101","#254870","#a0a3a2","#7a4c0f"),name="Treatment")+
  labs(title="Pectin",x="Treatment",y="")+
  geom_errorbar(aes(ymin =pectin_data $avg_carbon - pectin_data $sem_carbon, ymax = pectin_data $avg_carbon + pectin_data $sem_carbon),width = 0.2, colour = "black", position = position_dodge(.9))+
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

lignin_data <- subset(mean_carbon_rev,group=="Lignin")
lignin_gene_bar <- ggplot(lignin_data, aes(x = t2, y = avg_carbon, fill = Treat))+
  geom_col(alpha = 0.5,width = 0.7, position = 'dodge')+
  #  ylim(0,8)+
  scale_fill_manual(values = c("#8a0101","#254870","#a0a3a2","#7a4c0f"),name="Treatment")+
  labs(title="Lignin",x="Treatment",y="")+
  geom_errorbar(aes(ymin =lignin_data $avg_carbon - lignin_data$sem_carbon, ymax = lignin_data$avg_carbon + lignin_data$sem_carbon),width = 0.2, colour = "black", position = position_dodge(.9))+
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
  annotate(geom="text", x=1, y=0.0050, size=8,label="*",color="black", fontface="bold")+
  annotate(geom="text", x=3, y=0.0050, size=8,label="*",color="black", fontface="bold")

library("ggpubr")
tiff('figures/semianr_gene.tiff', units="in", width=15, height=5, res=300)
ggarrange(hemi_gene_bar,chitin_gene_bar ,pectin_gene_bar,lignin_gene_bar,labels = c("A", "B","C","D"),ncol = 4, nrow = 1,common.legend = T,legend = "bottom") 
dev.off()

##########stats for category########
library(lme4)
library(fBasics)
library(MASS)
library(car)
library(multcomp)
aromatic.all <- subset(carbon.cat.boxplot.group.df,group=="Aromatic");dim(aromatic.all)
aromatic.all.mod <- lmer(sum.value~ Cover* Nitrogen +(1|block), data = aromatic.all)
Anova(aromatic.all.mod)
summary(aromatic.all.mod)
coef(aromatic.all.mod)#no sig

cellulose.all <- subset(carbon.cat.boxplot.group.df,group=="Cellulose");dim(cellulose.all)
cellulose.all.mod <- lmer(sum.value~ Cover* Nitrogen +(1|block), data = cellulose.all)
Anova(cellulose.all.mod)
summary(cellulose.all.mod)
coef(cellulose.all.mod)#no sig

hemicellulose.all <- subset(carbon.cat.boxplot.group.df,group=="Hemicellulose");dim(hemicellulose.all)
hemicellulose.all.mod <- lmer(sum.value~ Cover* Nitrogen +(1|block), data = hemicellulose.all)
Anova(hemicellulose.all.mod)#sig on nitrogen
summary(hemicellulose.all.mod)
coef(hemicellulose.all.mod)#
nitrogen.hemicellulose<-glht(hemicellulose.all.mod ,linfct=mcp(Nitrogen="Tukey"))
summary(nitrogen.hemicellulose)#fertilization > no fer

chitin.all <- subset(carbon.cat.boxplot.group.df,group=="Chitin");dim(chitin.all)
chitin.all.mod <- lmer(sum.value~ Cover* Nitrogen +(1|block), data = chitin.all)
Anova(chitin.all.mod)
summary(chitin.all.mod)
coef(chitin.all.mod)
nitrogen.chitin<-glht(chitin.all.mod ,linfct=mcp(Nitrogen="Tukey"))
summary(nitrogen.chitin)#fertilization > no fer

lignin.all <- subset(carbon.cat.boxplot.group.df,group=="Lignin");dim(lignin.all)
lignin.all.mod <- lmer(sum.value~ Cover* Nitrogen +(1|block), data = lignin.all)
Anova(lignin.all.mod)
summary(lignin.all.mod)
coef(lignin.all.mod)#
nitrogen.lignin<-glht(lignin.all.mod ,linfct=mcp(Nitrogen="Tukey"))
summary(nitrogen.lignin)#fertilization < no fer

pectin.all <- subset(carbon.cat.boxplot.group.df,group=="Pectin");dim(pectin.all)
pectin.all.mod <- lmer(sum.value~ Cover* Nitrogen +(1|block), data = pectin.all)
Anova(pectin.all.mod)
summary(pectin.all.mod)
coef(pectin.all.mod)
nocover.pectin.all <- subset(pectin.all, Cover=="No_cover")
pectin.all.mod.1 <- lmer(sum.value~ Nitrogen +(1|block), data = nocover.pectin.all)
nitrogen.pectin.1<-glht(pectin.all.mod.1,linfct=mcp(Nitrogen="Tukey"))
summary(nitrogen.pectin.1)#no fer > fer under no cover
cover.pectin.all <- subset(pectin.all, Cover=="Vetch")
pectin.all.mod.2 <- lmer(sum.value~ Nitrogen +(1|block), data = cover.pectin.all)
nitrogen.pectin.2<-glht(pectin.all.mod.2,linfct=mcp(Nitrogen="Tukey"))
summary(nitrogen.pectin.2)#no sig

starch.all <- subset(carbon.cat.boxplot.group.df,group=="Starch");dim(starch.all)
starch.all.mod <- lmer(sum.value~ Cover* Nitrogen +(1|block), data = starch.all)
Anova(starch.all.mod)
summary(starch.all.mod)
coef(starch.all.mod) #no sig







##########This time I am going to subset the pathway of the carbon related gene in central carbon dataset ########
#1. TCA cycling
r_tca_central <- r_central_gene %>% subset(subheader=="TCA")
dim(r_tca_central)#79 17
#2. glycolysis
r_glycolysis_central <- r_central_gene %>% subset(subheader=="glycolysis")
dim(r_glycolysis_central)#88 17
#3.pentose pathway
r_pentose_central <- r_central_gene %>% subset(subheader=="pentose pathway")
dim(r_pentose_central)#28 17
############pyruvate metabolism
head(r_pyruvate_gene);dim(r_pyruvate_gene)



