# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("ggtree")

library("treeio")
library("ggtree")
library(ggplot2)
library(cowplot)
library(magick)
#install.packages("png")             # Install png package
library("png") 
#install.packages("patchwork")       # Install patchwork package
library("patchwork") 
library(ggforce)
library(here)
#https://yulab-smu.top/treedata-book/chapter4.html
set.seed(2017-02-16)
GDPmannose.tree <- read.tree("dataset/GDPmannose_best.tre")
GDPmannose.tree$tip.label <- c("Planctomycetes hermogutta sp.(NMC18978.1)","Planctomycetes (RMG00240.1) ",
                                "Acidobacteria (PYX52591.1)","Acidobacteria (TAM84186.1)","AMG_64073_40",
                               "Actinobacteria Streptomyces sp. BK340 (TVZ92930.1)",
                               "Actinobacteria Streptomyces geysiriensis (GGZ10435.1)")

# GDPmannose.tree$tip.label <- c("Proteobacteria Pseudomonas aeruginosa (1RPN)","M51_1683822_7",
#                                "M51_1276914_19","M51_1312024_5","Acidobacteria (PYX52591.1)",
#                                "M51_235266_2","Acidobacteria (TAM84186.1)","Actinobacteria Acidimicrobiales (MXV86448.1)",
#                                "Actinobacteria Acidimicrobiaceae (HCB34168.1)","M35_1469610_2","Actinobacteria (MBI4261762.1)",
#                                "AMG_64073_40","Actinobacteria Streptomyces sp. BK340 (TVZ92930.1)","Actinobacteria Streptomyces geysiriensis (GGZ10435.1)")
GDPmannose.tree$tip.label
length(GDPmannose.tree$tip.label)#7
GDPmannose.tree$node.label <- as.character(round(as.numeric(GDPmannose.tree$node.label),2)*100)
GDPmannose.tree$node.label
GDPmannose.node.type <- c("Ref","Ref","Ref","Ref","AMGs","Ref","Ref")
length(GDPmannose.node.type)
GDPmannose.info <-data.frame(GDPmannose.tree$tip.label,GDPmannose.node.type) 
colnames(GDPmannose.info) <- c("tip.label","Genes")
head(GDPmannose.info)
GDPmannose.p <- ggtree(GDPmannose.tree, layout="rectangular")+ 
  geom_treescale(x=2,y=0,fontsize=4, linesize=0.5, offset=0.2)+
 # geom_point(aes(shape=isTip, color=isTip), size=3)+
  geom_nodepoint(color="#b5e521", alpha=0, size=1)+
#  geom_tiplab(size=4, color="black")+
 # geom_hilight(node=12, fill="gold")+
  geom_nodelab(nudge_x = -0.023,nudge_y = 0.28)
GDPmannose.p2 <- GDPmannose.p %<+% GDPmannose.info
####
GDPmannose.p3 <- GDPmannose.p2 + geom_tiplab(offset =0.05, hjust = 0,size=5) +
  geom_tippoint(aes(color = Genes),size=4) +
  ggtitle("GDP-mannose 4,6-dehydratase [EC:4.2.1.47]")+
  theme(legend.position = "right",
        legend.title = element_text(size=15,face="plain",family = "Arial"),
        legend.text = element_text(size=15,family = "Arial"),
        plot.title=(element_text(size=15,family = "Arial",face="bold",vjust = 3,hjust=0.5)))+
scale_color_manual(values = c("#ad3232","#327aad", "black"))
#how to insert image into the tree plot
#https://stackoverflow.com/questions/9917049/inserting-an-image-to-ggplot2
  # ggsave(filename = paste0(here("/"), last_plot()$labels$title, ".png"),
  #        width = 5, height = 4, dpi = 300)
#annotate(geom="text", x=1.8, y=10, size=12,label="3D structre of AMG",color="black", fontface="plain")+#bold.italic
#annotation_raster(GDPmannose.structure, xmin = -1, xmax = 1, ymin = -1, ymax = 1)

#GDPmannose.structure <- readPNG("/Users/Ning/Desktop/R/New project/amg_chapter5/itasser/GDPmannose.png", native = TRUE)
#GDPmannose.structure <- image_read("/Users/Ning/Desktop/R/New project/amg_chapter5/itasser/GDPmannose.png")
# GDPmannose.structure %>%
# #  image_scale("100") %>%
# #  image_background("grey", flatten = TRUE) %>%
# # image_border("grey", "600x10") %>%
#   image_annotate("GDP-mannose 4,6-dehydratase (AMGs)", color = "white", size = 30,
#                  location = "+10+50", gravity = "northeast")
# 
#final_plot <- image_append(image_scale(c(GDPmannose.p3, GDPmannose.structure), "500"), stack = TRUE)
# GDPmannose.plot <- GDPmannose.p3+
#   inset_element(p=GDPmannose.structure,0.8,0.5,1,1)
# #  image_annotate("a", color = "white", size = 30, location = "+10+50", gravity = "northeast")
# GDPmannose.plot
##############################glycogen_synthase#####
set.seed(1010)
glycogen_synthase.tree <- read.tree("dataset/glycogen_synthase_best.tre")
glycogen_synthase.tree$tip.label
glycogen_synthase.tree$tip.label <- c("Bacteroidetes Chryseolinea sp. (WP_202010767.1) ","Bacteroidetes Cyclobacteriaceae (MBL7877674.1)",
                                      "Cyanobacteria Synechococcus elongatus (5ZE7)","Cyanobacteria Synechococcaceae WB9_3_282 (NDE23345.1)",
                                      "AMG_60918_5","Firmicutes Bacillus anthracis (3MBO)",
                                       "Firmicutes Staphylococcus aureus (6D9T)","Actinobacteria (MBN1288709.1)",
                                      "Actinobacteria Nitriliruptor alkaliphilus (WP_052667785.1)","Actinobacteria Nitriliruptoraceae (MTV26987.1)")

# glycogen_synthase.tree$tip.label <- c("M16_125857_18","Cyanobacteria Gloeobacter kilaueensis (WP_041244128.1)","M35_990922_1",
#                                       "Actinobacteria (MBN1288709.1)","M9_8016_3","Acidobacteria (MBI4729783.1)",
#                                       "M107_1193588_3","Bacteroidetes Cyclobacteriaceae (MBL7872191.1)","Bacteroidetes Cyclobacteriaceae (MBL7877674.1)",
#                                       "AMG_60918_5","Cyanobacteria Synechococcaceae WB9_3_282 (NDE23345.1)")
glycogen_synthase.tree$tip.label
length(glycogen_synthase.tree$tip.label)#10
glycogen_synthase.tree$node.label <- as.character(round(as.numeric(glycogen_synthase.tree$node.label),2)*100)
glycogen_synthase.tree$node.label
glycogen_synthase.node.type <- c("Ref","Ref",
                                 "Ref","Ref",
                                 "AMGs","Ref",
                                 "Ref","Ref",
                                 "Ref","Ref")
length(glycogen_synthase.node.type)
glycogen_synthase.info <-data.frame(glycogen_synthase.tree$tip.label,glycogen_synthase.node.type) 
colnames(glycogen_synthase.info) <- c("tip.label","Genes")
head(glycogen_synthase.info)
glycogen_synthase.p <- ggtree(glycogen_synthase.tree, layout="rectangular")+ 
  geom_treescale(x=6.5,y=0,fontsize=4, linesize=0.5, offset=0.2)+
  geom_nodepoint(color="#b5e521", alpha=0, size=1)+
  geom_nodelab(nudge_x = -0.08,nudge_y = 0.28)
glycogen_synthase.p2 <- glycogen_synthase.p %<+% glycogen_synthase.info
####
glycogen_synthase.p3 <- glycogen_synthase.p2 + geom_tiplab(offset =0.05, hjust = 0,size=5) +
  geom_tippoint(aes(color = Genes),size=4) +
  ggtitle("Glycogen synthase [EC:2.4.1.11]")+
  theme(legend.position = "right",
        legend.title = element_text(size=15,face="plain",family = "Arial"),
        legend.text = element_text(size=15,family = "Arial"),
        plot.title=(element_text(size=15,family = "Arial",face="bold",vjust = 3,hjust=0.5)))+
  scale_color_manual(values = c("#ad3232","#327aad", "black"))
#glycogen_synthase.structure <- readPNG("/Users/Ning/Desktop/R/New project/amg_chapter5/itasser/glycogen_synthase.png", native = TRUE)
#glycogen_synthase.plot <- glycogen_synthase.p3+
#  inset_element(p=glycogen_synthase.structure,0.8,0.5,1,1)
#glycogen_synthase.plot

##############################UDP_glucose#####
set.seed(1020)
UDP_glucose.tree <- read.tree("dataset/UDP_gluco_best.tre")
UDP_glucose.tree$tip.label
UDP_glucose.tree$tip.label <- c("Actinobacteria (TMM23743.1)","Actinobacteria (TMK27275.1)",
                                "Proteobacteria Escherichia coli (5GY7)","AMG_63920_73",
                                "Proteobacteria Betaproteobacteria (OQA34708.1)","Proteobacteria Alphaproteobacteria (MBI1252513.1)",
                                "Actinobacteria Streptomyces sp. (WP_186319230.1)","Actinobacteria Streptomyces luteocolor (WP_069883537.1)",
                                "Actinobacteria Streptomyces sp.(WP_121700703.1)","AMG_63951_18",
                                "Actinobacteria Streptomyces marianii (WP_171053074.1)","Actinobacteria Desertiactinospora gelatinilytica (WP_111166099.1)",
                                "Actinobacteria Salinispora arenicola (TQL39087.1)","Actinobacteria Streptomyces europaeiscabiei (WP_205576239.1)",
                                "Actinobacteria Streptomyces sp.(WP_151787052.1)","AMG_63948_14",
                                "Actinobacteria Saccharothrix ecbatanensis (MBB5803140.1)")
                                
                                # UDP_glucose.tree$tip.label <- c("Proteobacteria Escherichia coli K12 (5GY7)","AMG_63920_73","Proteobacteria Betaproteobacteria ADurb.Bin341 (OQA34708.1)",
#                                 "M16_508859_3","M51_1973616_3","M112_990769_4",
#                                 "Actinobacteria (TMM23743.1)","Actinobacteria (TMK27275.1)","AMG_63951_18",
#                                 "Actinobacteria Salinispora arenicola (TQL39087.1)","Actinobacteria Saccharothrix ecbatanensis (MBB5803140.1)","AMG_63948_14")

UDP_glucose.tree$tip.label
length(UDP_glucose.tree$tip.label)#19
UDP_glucose.tree$node.label <- as.character(round(as.numeric(UDP_glucose.tree$node.label),2)*100)
UDP_glucose.tree$node.label
UDP_glucose.node.type <- c("Ref","Ref",
                           "Ref","AMGs",
                           "Ref","Ref",
                           "Ref","Ref",
                           "Ref","AMGs",
                           "Ref","Ref",
                           "Ref","Ref",
                           "Ref","AMGs",
                           "Ref")
length(UDP_glucose.node.type)
UDP_glucose.info <-data.frame(UDP_glucose.tree$tip.label,UDP_glucose.node.type) 
colnames(UDP_glucose.info) <- c("tip.label","Genes")
head(UDP_glucose.info)
UDP_glucose.p <- ggtree(UDP_glucose.tree, layout="rectangular")+ 
  geom_treescale(x=4,y=0,fontsize=4, linesize=0.5, offset=0.2)+
  geom_nodepoint(color="#b5e521", alpha=0, size=1)+
  geom_nodelab(nudge_x = -0.05,nudge_y = 0.28)
UDP_glucose.p2 <- UDP_glucose.p %<+% UDP_glucose.info
####
UDP_glucose.p3 <- UDP_glucose.p2 + geom_tiplab(offset =0.05, hjust = 0,size=5) +
  geom_tippoint(aes(color = Genes),size=4) +
  ggtitle("UDP-glucose 4-epimerase [EC:5.1.3.2]")+
  theme(legend.position = "right",
        legend.title = element_text(size=15,face="plain",family = "Arial"),
        legend.text = element_text(size=15,family = "Arial"),
        plot.title=(element_text(size=15,family = "Arial",face="bold",vjust = 3,hjust=0.5)))+
  scale_color_manual(values = c("#ad3232","#327aad", "black"))
# UDP_glucose.structure <- readPNG("/Users/Ning/Desktop/R/New project/amg_chapter5/itasser/UDP_glucose_73.png", native = TRUE)
# UDP_glucose.plot <- UDP_glucose.p3+
#   inset_element(p=UDP_glucose.structure,0.8,0.5,1,1)
# UDP_glucose.plot


##############################UDP_galactopyranose#####
set.seed(1030)
UDP_galactopyranose.tree <- read.tree("dataset/UDP_galact_best.tre")
UDP_galactopyranose.tree$tip.label
UDP_galactopyranose.tree$tip.label <- c("Proteobacteria Klebsiella pneumoniae (2BI8)","Proteobacteria Escherichia coli (1I8T)",
                                        "Proteobacteria Campylobacter (4MO2)","Acidobacteriia (PWU02268.1)",
                                        "Acidobacteria Candidatus Solibacter usitatus (ABJ81361.1)","Planctomycetes Blastopirellula sp. (MAR08810.1)",
                                        "AMG_63920_84"
                                        )

# UDP_galactopyranose.tree$tip.label <- c("Fibrobacteres Chitinispirillaceae (MBN1600069.1)","M112_706273_2","Elusimicrobia (MBI4425081.1)",
#                                         "Proteobacteria Deltaproteobacteria (MBI5827104.1)","Proteobacteria Klebsiella pneumoniae (2BI8)","Proteobacteria Escherichia coli (1I8T)",
#                                         "AMG_63920_84","Proteobacteria Nisaea sp. (MAJ00350.1)")
UDP_galactopyranose.tree$tip.label
length(UDP_galactopyranose.tree$tip.label)#8
UDP_galactopyranose.tree$node.label <- as.character(round(as.numeric(UDP_galactopyranose.tree$node.label),2)*100)
UDP_galactopyranose.tree$node.label
UDP_galactopyranose.node.type <- c("Ref","Ref",
                                 "Ref","Ref",
                                 "Ref","Ref",
                                 "AMGs"
                                 )
length(UDP_galactopyranose.node.type)
UDP_galactopyranose.info <-data.frame(UDP_galactopyranose.tree$tip.label,UDP_galactopyranose.node.type) 
colnames(UDP_galactopyranose.info) <- c("tip.label","Genes")
head(UDP_galactopyranose.info)
UDP_galactopyranose.p <- ggtree(UDP_galactopyranose.tree, layout="rectangular")+ 
  geom_treescale(x=5,y=0,fontsize=4, linesize=0.5, offset=0.2)+
  geom_nodepoint(color="#b5e521", alpha=0, size=1)+
  geom_nodelab(nudge_x = -0.07,nudge_y = 0.36)
UDP_galactopyranose.p
UDP_galactopyranose.p2 <- UDP_galactopyranose.p %<+% UDP_galactopyranose.info
####
UDP_galactopyranose.p3 <- UDP_galactopyranose.p2 + geom_tiplab(offset =0.05, hjust = 0,size=5) +
  geom_tippoint(aes(color = Genes),size=4) +
  ggtitle("UDP-galactopyranose mutase [EC:5.4.99.9]")+
  theme(legend.position = "right",
        legend.title = element_text(size=15,face="plain",family = "Arial"),
        legend.text = element_text(size=15,family = "Arial"),
        plot.title=(element_text(size=15,family = "Arial",face="bold",vjust = 3,hjust=0.5)))+
  scale_color_manual(values = c("#ad3232","#327aad", "black"))
# UDP_galactopyranose.structure <- readPNG("/Users/Ning/Desktop/R/New project/amg_chapter5/itasser/UDP_galactopyranose.png", native = TRUE)
# UDP_galactopyranose.plot <- UDP_galactopyranose.p3+
#   inset_element(p=UDP_galactopyranose.structure,0.8,0.5,1,1)
 
  
#UDP_galactopyranose.plot

library("ggpubr")
tiff('figures/tree.png', units="in", width=17, height=17, res=300)
ggarrange(glycogen_synthase.p3,GDPmannose.p3,UDP_galactopyranose.p3,UDP_glucose.p3,
          labels = c("A","B","C","D"),ncol = 1, nrow = 4,common.legend = T,legend = "bottom",
          label.x =-0.001, label.y = 1.002)
dev.off()


####combine figures and image
#GDPmannose.plot <- ggsave(GDPmannose.p3,device = "png",path = "/Users/Ning/Desktop/R/New project/amg_chapter5/figures/",dpi = 300,)
GDPmannose.structure <- image_read("/Users/Ning/Desktop/R/New project/amg_chapter5/itasser/GDPmannose.png") %>%
  image_scale("550") %>%
  image_annotate("AMG_60743_40", color = "white", size = 60,
                 location = "+70+0", gravity = "northeast")
glycogen_synthase.structure <- image_read("/Users/Ning/Desktop/R/New project/amg_chapter5/itasser/glycogen_synthase.png")%>%
  image_scale("550") %>%
  image_annotate("AMG_60918_5", color = "white", size = 60,
                 location = "+70+0", gravity = "northeast")
UDP_glucose_73.structure <- image_read("/Users/Ning/Desktop/R/New project/amg_chapter5/itasser/UDP_glucose_73.png")%>%
  image_scale("550") %>%
  image_annotate("AMG_63951_18", color = "white", size = 60,
                 location = "+70+0", gravity = "northeast")
UDP_galactopyranose.structure <- image_read("/Users/Ning/Desktop/R/New project/amg_chapter5/itasser/UDP_galactopyranose.png")%>%
  image_scale("550") %>%
  image_annotate("AMG_63920_84", color = "white", size = 60,
                 location = "+70+0", gravity = "northeast")

four.tree <- image_read("/Users/Ning/Desktop/R/New project/amg_chapter5/figures/tree.png")

#GDPmannose_plot <- image_append(image_scale(c(GDPmannose.p3, four.tree), "500"))
four.tree_glycogen <- image_composite(image_scale(four.tree), glycogen_synthase.structure, offset = "+4500+100")
four.tree_glycogen_mannose <- image_composite(image_scale(four.tree_glycogen),GDPmannose.structure,offset = "+4500+1300")
four.tree_glycogen_mannose_galacto <- image_composite(image_scale(four.tree_glycogen_mannose),UDP_galactopyranose.structure,offset = "+4500+2600")
four.tree_glycogen_mannose_galacto_gluco73 <- image_composite(image_scale(four.tree_glycogen_mannose_galacto),UDP_glucose_73.structure,offset = "+4500+3900")
image_write(four.tree_glycogen_mannose_galacto_gluco73,path = "/Users/Ning/Desktop/R/New project/amg_chapter5/figures/advtree.tiff",density="300",format = "tiff")

