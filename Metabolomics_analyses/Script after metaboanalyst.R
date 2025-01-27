## Metabolomics analysis with R

library(ggplot2)
library(DT)
library(reshape2)
library(ggrepel)
library(vegan)
library(pairwiseAdonis) 
library(metR) 
library(cowplot)
library(ggpubr)
library(ape)
library(dplyr)
library(agricolae) #for statistical tests
library(RVAideMemoire) #for statistical tests
setwd("your path here")

#dir.create("1_Class_difference_results")
#dir.create("2_Photic_zone_difference_all_results")
#dir.create("3_Photic_zone_difference_demospongiae_results")
#dir.create("4_Photic_zone_difference_hexactinellida_results")


# 1. Analysis Q1 : difference between Demospongiae and Hexactinellida ####

color_vector_class = c("Demospongiae" = "chartreuse", 
                       "Hexactinellida" = "purple3")

color_vector_zone = c( "1.Mesophotic" = "yellow1",
                       "2.Upper rariphotic" = "aquamarine",
                       "3.Lower rariphotic" = "royalblue")

color_vector_genus = c("Aciculites" = "forestgreen",
                       "Calthropella" = "chartreuse4",
                       "Cinachyrella" = "lightgreen", 
                       "Gastrophanella" = "lawngreen", 
                       "Geodia" = "springgreen",
                       "Penares" = "seagreen", 
                       "Neopetrosia" = "orangered",
                       "Petrosia" = "darkorange",
                       "Svenzea" = "deeppink3",
                       "Topsentia" = "hotpink",
                       "Biemna" = "firebrick",
                       "Conorete" = "navyblue",
                       "Dictyoplax" = "dodgerblue3",
                       "Lefroyella" = "dodgerblue2",
                       "Dactylocalyx" = "skyblue",
                       "Heterotella" = "lightslateblue",
                       "Hexactinella" = "purple",
                       "Myliusia" = "darkslategray", 
                       "Regadrella" = "turquoise",
                       "Verrucocoeloidea" = "darkcyan", 
                       "Hexactinellida" = "skyblue4")


## 1.1. PCA analysis all samples ####

PCA_score = read.csv(file = "./pca_score_pareto.csv" , header = TRUE , sep = ";" , row.names = 1)
dim(PCA_score)
datatable(PCA_score)

plot_pca = ggplot(PCA_score, aes(x = PC1, y = PC2))
plot_pca = plot_pca + geom_point(aes(fill = Photic_zone, pch = Class_host) , stroke = 1.5 , size = 7, alpha = 0.7, color = "black" , show.legend = T)
plot_pca = plot_pca + scale_shape_manual(values = c(22, 23))
plot_pca = plot_pca + theme_bw(base_size = 20) 
plot_pca = plot_pca + scale_fill_manual(values = color_vector_zone)
plot_pca = plot_pca +  xlab("PC1(18 %)") + ylab("PC2(11 %)")
plot_pca = plot_pca + guides(fill = guide_legend(override.aes = list(shape = 21)))
plot_pca = plot_pca +  labs(shape = "Sponge class", fill = "Photic zone")
plot_pca

plot_pca_all_class = plot_pca

ggsave(filename = "PCA_plot_all_class.pdf", 
       plot = plot_pca_all_class, 
       device = "pdf" , 
       width = 35 , height = 25, units = "cm", 
       path = "./1_Class_difference_results")


plot_pca = ggplot(PCA_score, aes(x = PC1, y = PC2))
plot_pca = plot_pca + geom_point(aes(fill = Genus_host, pch = Class_host) , stroke = 1.5 , size = 7, alpha = 1, color = "black" , show.legend = T)
plot_pca = plot_pca + scale_shape_manual(values = c(22, 23))
plot_pca = plot_pca + scale_fill_manual(values = color_vector_genus)
plot_pca = plot_pca + theme_bw(base_size = 20) 
plot_pca = plot_pca + geom_text_repel(aes(label = Code, color = Genus_host), size = 6, box.padding = 0.8)
plot_pca = plot_pca + scale_color_manual(values = color_vector_genus)
plot_pca = plot_pca +  xlab("PC1(18 %)") + ylab("PC2(11 %)")
plot_pca = plot_pca + labs(shape = "Sponge class", fill = "Sponge genus", color = "Sponge genus")
plot_pca = plot_pca + guides(fill = guide_legend(override.aes = list(shape = 21)))
plot_pca
plot_pca_all_genus = plot_pca

ggsave(filename = "PCA_plot_all_genus.pdf", 
       plot = plot_pca_all_genus, 
       device = "pdf" , 
       width = 35 , height = 25, units = "cm", 
       path = "./1_Class_difference_results")



### to do ordisurf line with the depth, the PCA need to be done with R directly with the normalized dataset

Data_table_normalized = read.csv(file = "./data_normalized_adonis_filtered.csv" , header = TRUE , sep = ";" , row.names = 1)
datatable(Data_table_normalized)

data_distbeta = vegdist(Data_table_normalized, method="euclidean")
data_distbeta

Metadata = read.csv(file = "./Metadata_table.csv" , header = TRUE , sep = ";" , row.names = 1)
datatable(Metadata)


pcoa <- wcmdscale(data_distbeta, eig = TRUE)
pcoa$points
pcoa

data_pca = pcoa$points
data_pca
#colnames(data_pca) <- c("PC1", "PC2")

data_pca
ordiplot(pcoa, display = 'sites', type = 'text')

ordi = ordisurf(pcoa ~ Depth_m, Metadata)
ordi$grid
ordi.grid <- ordi$grid #extracts the ordisurf object
ordi.grid
str(ordi.grid) #it's a list though - cannot be plotted as is
ordi.mite <- expand.grid(x = ordi.grid$x, y = ordi.grid$y) #get x and ys
ordi.mite
ordi.mite$z <- as.vector(ordi.grid$z) #unravel the matrix for the z scores
ordi.mite.na <- data.frame(na.omit(ordi.mite)) #gets rid of the nas
ordi.mite.na

ordi.mite.na$z


data_pca = cbind(data_pca , Metadata)

plot_pca = ggplot(data_pca, aes(x = Dim1  , y = Dim2))
plot_pca = plot_pca + geom_point(aes(fill = Photic_zone, pch = Class_host), size = 7, alpha = 0.8)
plot_pca = plot_pca + scale_shape_manual(values = c(22,23))
plot_pca = plot_pca + scale_fill_manual(values = color_vector_zone, labels = c("Mesophotic", "Upper rariphotic", "Lower rariphotic"))
plot_pca = plot_pca + theme(plot.title = element_text(size=22, face="bold"), axis.text=element_text(size=12),axis.title=element_text(size=14))
plot_pca = plot_pca + guides(fill = guide_legend(override.aes = list(shape = 22, size = 7)))
plot_pca = plot_pca + stat_contour(data = ordi.mite.na, aes(x = x, y = y, z = z, colour = after_stat(level)))
plot_pca =  plot_pca + scale_colour_gradient2(low = "gold", mid = "aquamarine", high = "darkblue", midpoint = 200) 
plot_pca = plot_pca +  xlab("PC1(25 %)") + ylab("PC2(13 %)")+ theme_bw(base_size = 25)
plot_pca = plot_pca + labs(shape = "Sponge class", fill = "Photic zone", colour = "Depth")
plot_pca
plot_pca_classzone_ordi = plot_pca


plot_pca = ggplot(ordi.mite.na, aes(x = x  , y = y))
plot_pca = plot_pca + theme(plot.title = element_text(size=22, face="bold"), axis.text=element_text(size=12),axis.title=element_text(size=14))
plot_pca = plot_pca + geom_contour2(aes(z = z, colour = stat(level))) + metR::geom_text_contour(aes(z = z), stroke = 0.15) 
plot_pca = plot_pca + scale_colour_gradient2(low = "gold", mid = "aquamarine", high = "darkblue", midpoint = 200, guide=guide_colourbar(reverse = TRUE)) 
plot_pca = plot_pca + geom_point(data = data_pca , aes(x = Dim1  , y = Dim2, fill = Photic_zone, pch = Class_host), size = 7, alpha = 0.8)
plot_pca = plot_pca + scale_shape_manual(values = c(22,23))
plot_pca = plot_pca + guides(fill = guide_legend(override.aes = list(shape = 22, size = 7)))
plot_pca = plot_pca + scale_fill_manual(values = color_vector_zone, labels = c("Mesophotic", "Upper rariphotic", "Lower rariphotic")) 
plot_pca = plot_pca + xlab("PC1(25 %)") + ylab("PC2(13 %)") + theme_bw(base_size = 25)
plot_pca = plot_pca + labs(shape = "Sponge class", fill = "Photic zone", colour = "Depth")  
plot_pca_classzone_ordi2 = plot_pca 
plot_pca_classzone_ordi2
  
ggsave(filename = "PCA_plot_all_class_ordi2_filtered.pdf", 
       plot = plot_pca_classzone_ordi2, 
       device = "pdf" , 
       width = 35 , height = 25, units = "cm", 
       path = "./1_Class_difference_results")

plot_pca = ggplot(ordi.mite.na, aes(x = x  , y = y))
plot_pca = plot_pca + theme(plot.title = element_text(size=22, face="bold"), axis.text=element_text(size=12),axis.title=element_text(size=14))
plot_pca = plot_pca + geom_contour2(aes(z = z, colour = stat(level))) + metR::geom_text_contour(aes(z = z), stroke = 0.15) 
plot_pca = plot_pca + scale_colour_gradient2(low = "gold", mid = "aquamarine", high = "darkblue", midpoint = 200) 
plot_pca = plot_pca + geom_point(data = data_pca , aes(x = Dim1  , y = Dim2, fill = Genus_host, pch = Class_host), size = 7, alpha = 0.8)
plot_pca = plot_pca + scale_shape_manual(values = c(22,23))
plot_pca = plot_pca + guides(fill = guide_legend(override.aes = list(shape = 22, size = 7)))
plot_pca = plot_pca + scale_fill_manual(values = color_vector_genus) 
plot_pca = plot_pca + xlab("PC1(25 %)") + ylab("PC2(13 %)") + theme_bw(base_size = 25)
plot_pca = plot_pca + labs(shape = "Sponge class", fill = "Sponge genus", colour = "Depth")  
plot_pca = plot_pca + guides(shape = guide_none(), colour = guide_none())
plot_pca_genus_ordi2 = plot_pca 
plot_pca_genus_ordi2

ggsave(filename = "PCA_plot_all_genus_ordi2_filtered.pdf", 
       plot = plot_pca_genus_ordi2, 
       device = "pdf" , 
       width = 35 , height = 25, units = "cm", 
       path = "./1_Class_difference_results")



## 1.2. PLSDA plot analysis #### 

PLSDA_score = read.csv(file = "./plsda_score_pareto.csv", header = TRUE , sep = ";" , row.names = 1)
dim(PLSDA_score)
datatable(PLSDA_score)


plot_plsda = ggplot(PLSDA_score, aes(x = Comp.1, y = Comp.2))
plot_plsda = plot_plsda + geom_point(aes(fill = Photic_zone, pch = Class_host) , stroke = 1.5 , size = 7, alpha = 0.7, color = "black" , show.legend = T)
plot_plsda = plot_plsda + scale_shape_manual(values = c(21, 23))
plot_plsda = plot_plsda + theme_bw(base_size = 20) 
plot_plsda = plot_plsda + scale_fill_manual(values = color_vector_zone)
plot_plsda = plot_plsda +  xlab("Component 1 (8%)") + ylab("Component 2 (5%)")
plot_plsda = plot_plsda + guides(fill = guide_legend(override.aes = list(shape = 21)))
plot_plsda = plot_plsda +  labs(shape = "Sponge class", fill = "Photic zone")
plot_plsda

plot_plsda_all_class = plot_plsda

ggsave(filename = "PLSDA_plot_all_class.pdf", 
       plot = plot_plsda_all_photic, 
       device = "pdf" , 
       width = 35 , height = 25, units = "cm", 
       path = "./1_Class_difference_results")




## 1.3. Two-way permanova class*zone ####


Data_table_normalized = read.csv(file = "./data_normalized_adonis.csv" , header = TRUE , sep = ";" , row.names = 1)
datatable(Data_table_normalized)


Metadata = read.csv(file = "./Metadata_table.csv" , header = TRUE , sep = ";" , row.names = 1)
datatable(Metadata)


data_distbeta = dist(Data_table_normalized, method="euclidean")
data_distbeta

data_permanova = adonis2(data_distbeta ~ Class_host*Photic_zone, data = Metadata)
data_permanova
datatable(data_permanova)

write.csv(as.data.frame(data_permanova), 
          file.path("./1_Class_difference_results" , "Data_twoway_permnova_class_zone.csv"))


## 1.4. Nested permanova Class / genus ####

data_permanova = adonis2(data_distbeta ~ Class_host/Genus_host, data = Metadata)
data_permanova
datatable(data_permanova)

write.csv(as.data.frame(data_permanova), 
          file.path("./1_Class_difference_results" , "Data_twoway_permnova_class_genus.csv"))


data_permanova = adonis2(data_distbeta ~ Class_host/Genus_host/Species_host, data = Metadata)
data_permanova
datatable(data_permanova)

write.csv(as.data.frame(data_permanova), 
          file.path("./1_Class_difference_results" , "Data_twoway_permnova_class_genus_species.csv"))



## 1.5.Pairwise adonis ####

Metadata["Class_zone"] <- paste(Metadata$Class_host, Metadata$Photic_zone)
datatable(Metadata)

data_pairwiseadonis_classzone= pairwise.adonis(data_distbeta, Metadata$Class_zone)
data_pairwiseadonis_classzone
datatable(data_pairwiseadonis_classzone)

write.csv(data_pairwiseadonis_classzone, 
          file.path("./1_Class_difference_results" , "Pairwise_adonis_classzone.csv"))


data_pairwiseadonis_genus = pairwise.adonis(data_distbeta, Metadata$Genus_host)
data_pairwiseadonis_genus
datatable(data_pairwiseadonis_genus)

write.csv(data_pairwiseadonis_genus, 
          file.path("./1_Class_difference_results" , "Pairwise_adonis_genus.csv"))



## 1.6. Boxplots all VIPS Demospongiae vs Hexactinellida ####

Data_table_normalized = read.csv(file = "./data_normalized_class.csv" , header = TRUE , sep = ";" , row.names = 1)

datatable(Data_table_normalized)
dim(Data_table_normalized)

VIP_table_normalized = Data_table_normalized[, c("Sample_names","Label", 
                                                 "X2064_510.35064_957.99", #Lyso-PC(C17:0)*
                                                 "X2022_544.33574_897.56", #Lyso-PC(C20:4)**
                                                 #"X2126_524.36706_1070.53", #Lyso-PC(C18:0)
                                                 #"X2106_496.37124_1030.26", #Lyso-PC(O-C17:0) NO6P
                                                 "X2015_544.33523_883.64", #Lyso-PC(C20:4)**
                                                 #"X1991_542.31997_846.88", #Lyso-PC(C20:5)
                                                 "X2073_510.35004_975.58", #Lyso-PC(C17:0)*
                                                 #"X2030_508.3351_903.47", #Lyso-PC(C17:1)
                                                 #"X2127_510.3875_1069.92", #Lyso-PC(O-C18:0) NO6P  
                                                 #"X1117_538.38277_1168.36", #Lyso-PC(C19:0)
                                                 #"X1989_542.3203_833.94", #Lyso-PC(C20:5)
                                                 #"X2148_550.38314_1095.4", #Lyso-PC(C20:1)
                                                 #"X2054_496.33453_941.01",  #Lyso-PC(C16:0)
                                                 #"X2192_597.33806_1141.14", #Lyso-PI(C19:0)
                                                 #"X1008_583.32127_1075.61", #Lyso-PI(C18:0)
                                                 "X2007_502.28847_877.11", #Lyso-PE(C20:4)_1
                                                 "X1999_502.28802_864.01", #Lyso-PE(C20:4)_2
                                                 "X1888_539.0388_522.08", #C20H24Br2N6O2 Aphrocallistin
                                                 "X1884_569.04746_504.2", #C21H26Br2N6O3 Aphrocallistin der
                                                 #"X1626_187.10489_69.52", #Alanylproline
                                                 "X1839_206.04277_275.93" #Xanthurenic acid
                                                 #"X465_246.16695_320.81" #Valerylcarnitine
                                                 )]


VIP_table_normalized
datatable(VIP_table_normalized)
dim(VIP_table_normalized)



#to group together the duplicated VIPs
colnames(VIP_table_normalized) <- c("Sample_names","Sponge_class", 
                                    "Lyso-PC(C17:0)", 
                                    "Lyso-PC(C20:4)", 
                                    #"Lyso-PC(C18:0)", 
                                    #"Lyso-PC(O-C17:0)", 
                                    "Lyso-PC(C20:4)", 
                                    #"Lyso-PC(C20:5)", 
                                    "Lyso-PC(C17:0)", 
                                    #"Lyso-PC(C17:1)", 
                                    #"Lyso-PC(O-C18:0)", 
                                    #"Lyso-PC(19:0)", 
                                    #"Lyso-PC(C20:5)", 
                                    #"Lyso-PC(C20:1)", 
                                    #"Lyso-PC(C16:0)", 
                                    #"Lyso-PI(C19:0)", 
                                    #"Lyso-PI(C18:0)", 
                                    "Lyso-PE(C20:4)", 
                                    "Lyso-PE(C20:4)", 
                                    "Aphrocallistin", 
                                    "Aphrocallistin_derivative",
                                    #"Alanylproline", 
                                    "Xanthurenic acid") 
                                   # "Valerylcarnitine")

datatable(VIP_table_normalized)
dim(VIP_table_normalized)

data_long = melt(VIP_table_normalized, id.var=c("Sample_names","Sponge_class"))
datatable(data_long)

boxplot_vip = ggplot(data_long, aes(x = Sponge_class, y = value)) 
boxplot_vip = boxplot_vip + geom_point(aes(fill = Sponge_class), size = 2, alpha = 0.8, pch = 21)
boxplot_vip = boxplot_vip + geom_boxplot(aes(fill = Sponge_class), alpha = 0.4, size = 0.8)
boxplot_vip = boxplot_vip + scale_fill_manual(values = color_vector_class)
boxplot_vip = boxplot_vip + theme_bw(base_size = 10)  + theme(legend.position="left")
boxplot_vip = boxplot_vip + facet_wrap(~variable)
#boxplot_vip = boxplot_vip + facet_grid(variable~Site, scales = "free")
boxplot_vip = boxplot_vip + theme(axis.title.x = element_blank())
boxplot_vip = boxplot_vip + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
boxplot_vip = boxplot_vip + labs(fill = "Sponge class") + ylab("Normalized concentration")
boxplot_vip

boxplot_vip_class_main = boxplot_vip

ggsave(filename = "Boxplot_VIP_class.pdf", 
       plot = boxplot_vip_class, 
       device = "pdf" , 
       width = 30 , height = 25, units = "cm", 
       path = "./1_Class_difference_results")

#### specific analysis for Aphrocallistin X1888_539.0388_522.08 C19H29Br2N2O6



VIP_table_normalized_aphro_xa = Data_table_normalized[, c("Sample_names","Label", 
                                                 "X1888_539.0388_522.08", #C19H29Br2N2O6 Aphrocallistin
                                                 "X1884_569.04746_504.2" , #Aphro derivative 
                                                 "X638_525.02068_488.92", # Aphrocallistin C ?
                                                 "X1839_206.04277_275.93" #Xanthurenic acid

)]

VIP_table_normalized_aphro_xa

dim(VIP_table_normalized_aphro_xa)

VIP_table_normalized_aphro_xa = VIP_table_normalized_aphro_xa[order(VIP_table_normalized_aphro_xa$Sample_names   ),]


Metadata = read.csv(file = "./Metadata_table.csv" , header = TRUE , sep = ";" , row.names = 1)
dim(Metadata)
Metadata


VIP_table_normalized_aphro_xa_meta = cbind(VIP_table_normalized_aphro_xa, Metadata$Genus_host)
colnames(VIP_table_normalized_aphro_xa_meta) <- c("Sample_names","Sponge_class", 
                                           "Aphrocallistin", 
                                           "Aphrocallistin_derivative", 
                                           "Aphrocallistin_c",
                                           "Xanthurenic_acid", "Genus_host")
VIP_table_normalized_aphro_xa_meta

data_long_VIP_table_normalized_aphro_xa_meta = melt(VIP_table_normalized_aphro_xa_meta, id.var=c("Sample_names","Sponge_class", "Genus_host"))
data_long_VIP_table_normalized_aphro_xa_meta
colnames(data_long_VIP_table_normalized_aphro_xa_meta) <- c("Sample_names","Sponge_class","Genus_host", "Compound", "Normalized_concentration")
data_long_VIP_table_normalized_aphro_xa_meta

boxplot_vip = ggplot(data_long_VIP_table_normalized_aphro_xa_meta, aes(x = Genus_host, y = Normalized_concentration))
boxplot_vip = boxplot_vip + facet_grid( Compound ~ Sponge_class , scales="free", space = "free_x") 
boxplot_vip = boxplot_vip + geom_point(aes(fill = Genus_host), size = 2, alpha = 0.8, pch = 21)
boxplot_vip = boxplot_vip + geom_boxplot(aes(fill = Genus_host), alpha = 0.4, size = 0.8)
boxplot_vip = boxplot_vip + scale_fill_manual(values = color_vector_genus)
boxplot_vip = boxplot_vip + theme_bw(base_size = 15)  + theme(legend.position="left")
boxplot_vip = boxplot_vip + theme(axis.title.x = element_blank())
boxplot_vip = boxplot_vip + theme(axis.text.x = element_text(angle = 45, hjust = 1))
boxplot_vip = boxplot_vip + labs(fill = "Sponge genus") + ylab("Normalized concentration")
boxplot_vip
boxplot_vip_aphro_xa = boxplot_vip

ggsave(filename = "Boxplot_VIP_aphro_xa.pdf", 
       plot = boxplot_vip_aphro_xa, 
       device = "pdf" , 
       width = 35 , height = 25, units = "cm", 
       path = "./1_Class_difference_results")

VIP_table_normalized_aphro_xa_meta$Xanthurenic_acid

aov_Aphrocallistin = aov(Aphrocallistin  ~ Genus_host, VIP_table_normalized_aphro_xa_meta)
aov_Aphrocallistin 
summary(aov_Aphrocallistin)

hsd_Aphrocallistin = HSD.test(aov(Aphrocallistin  ~ Genus_host, VIP_table_normalized_aphro_xa_meta), "Genus_host", group=T)
hsd_Aphrocallistin



aov_Aphrocallistin_der = aov(Aphrocallistin_derivative  ~ Genus_host, VIP_table_normalized_aphro_xa_meta)
aov_Aphrocallistin_der 
summary(aov_Aphrocallistin_der)

hsd_Aphrocallistin_der = HSD.test(aov(Aphrocallistin_derivative  ~ Genus_host, VIP_table_normalized_aphro_xa_meta), "Genus_host", group=T)
hsd_Aphrocallistin_der



aov_XA= aov(Xanthurenic_acid  ~ Genus_host, VIP_table_normalized_aphro_xa_meta)
aov_XA 
summary(aov_XA)

hsd_XA = HSD.test(aov(Xanthurenic_acid  ~ Genus_host, VIP_table_normalized_aphro_xa_meta), "Genus_host", group=T)
hsd_XA


# 2. Analysis Q2 : difference between photic zone ####


## 2.1. PLSDA plot analysis with all samples #### 

PLSDA_score = read.csv(file = "./plsda_score_pareto.csv", header = TRUE , sep = ";" , row.names = 1)
dim(PLSDA_score)
datatable(PLSDA_score)


plot_plsda = ggplot(PLSDA_score, aes(x = Comp.1, y = Comp.2))
plot_plsda = plot_plsda + geom_point(aes(fill = Photic_zone, pch = Class_host) , stroke = 1.5 , size = 7, alpha = 0.7, color = "black" , show.legend = T)
plot_plsda = plot_plsda + scale_shape_manual(values = c(22, 23))
plot_plsda = plot_plsda + theme_bw(base_size = 25) 
plot_plsda = plot_plsda + scale_fill_manual(values = color_vector_zone)
plot_plsda = plot_plsda +  xlab("Component 1 (7%)") + ylab("Component 2 (4%)")
plot_plsda = plot_plsda + guides(fill = guide_legend(override.aes = list(shape = 22)))
plot_plsda = plot_plsda +  labs(shape = "Sponge class", fill = "Photic zone")
plot_plsda

plot_plsda_all_photic = plot_plsda

ggsave(filename = "PLSDA_plot_all_photic.pdf", 
       plot = plot_plsda_all_photic, 
       device = "pdf" , 
       width = 35 , height = 25, units = "cm", 
       path = "./2_Photic_zone_difference_all_results")



## Triplot PCA zone, PCA genus and PLSDA zone all

legend_pca_zone = as_ggplot(get_legend(plot_pca_classzone_ordi2+ theme_bw(base_size = 20)))
legend_pca_zone

legend_pca_genus = as_ggplot(get_legend(plot_pca_genus_ordi2+ theme_bw(base_size = 20)))
legend_pca_genus

plot_pca_classzone_ordi2_f = plot_pca_classzone_ordi2 + theme(legend.position="none")
plot_pca_classzone_ordi2_f

plot_pca_genus_ordi2_f = plot_pca_genus_ordi2 + theme(legend.position="none")
plot_pca_genus_ordi2_f

plot_plsda_all_photic_f = plot_plsda_all_photic + theme(legend.position="none")
plot_plsda_all_photic_f


legend_all = plot_grid(legend_pca_zone, legend_pca_genus  , ncol = 2 )
legend_all

triplot_pca_plsda_ordi = plot_grid( plot_pca_classzone_ordi2_f , plot_pca_genus_ordi2_f , plot_plsda_all_photic_f ,  legend_all , labels=c("A", "B" ,"C", ""),rel_widths = c(1, 1 , 1, 2 )  ,  ncol = 2, nrow = 2 ,label_size = 20)
triplot_pca_plsda_ordi


ggsave(filename = "Triplot_pca_plsda_ordi_filtered.pdf", 
       plot = triplot_pca_plsda_ordi, 
       device = "pdf" , 
       width = 40 , height = 35, units = "cm", 
       path = "./2_Photic_zone_difference_all_results")



## 2.2. PLSDA plot analysis with Demospongiae #### 


PLSDA_score = read.csv(file = "./plsda_score_zone_demospongiae.csv", header = TRUE , sep = ";" , row.names = 1)
dim(PLSDA_score)
datatable(PLSDA_score)

plot_plsda = ggplot(PLSDA_score, aes(x = Comp.1, y = Comp.2))
plot_plsda = plot_plsda + geom_point(aes(fill = Photic_zone, pch = Class_host) , pch = 21, stroke = 1.5 , size = 7, alpha = 0.7, color = "black" , show.legend = T)
plot_plsda = plot_plsda + theme_bw(base_size = 20) 
plot_plsda = plot_plsda + scale_fill_manual(values = color_vector_zone, labels = c("Mesophotic", "Upper rariphotic", "Lower rariphotic"))
plot_plsda = plot_plsda +  xlab("Component 1 (7%)") + ylab("Component 2 (4%)")
plot_plsda = plot_plsda + guides(fill = guide_legend(override.aes = list(shape = 21)))
plot_plsda = plot_plsda + geom_text_repel(aes(label = Code), size = 6, box.padding = 0.8)
plot_plsda = plot_plsda +  labs(shape = "Sponge class", fill = "Photic zone")
plot_plsda

plot_plsda_demospongiae_photic = plot_plsda

ggsave(filename = "PLSDA_plot_demospongiae_photic.pdf", 
       plot = plot_plsda_demospongiae_photic, 
       device = "pdf" , 
       width = 25 , height = 15, units = "cm", 
       path = "./3_Photic_zone_difference_demospongiae_results")





## 2.3. Boxplots all VIPS zone for Demospongiae ####

Data_table_normalized = read.csv(file = "./data_normalized_demospongiae.csv" , header = TRUE , sep = ";" , row.names = 1)

datatable(Data_table_normalized)
dim(Data_table_normalized)


VIP_table_normalized = Data_table_normalized[, c("Sample_names","Label", 
                                                 "X2127_510.3875_1069.92", #Lyso-PC(O-C18:0)
                                                  "X1717_176.12565_89.24", #Rhodosamine
                                                 #"X1638_155.07952_69.96", #C7H11N2O2  Cyclo-(Pro-Gly)?
                                                 #"X406_174.12135_91.21", #C7H15N3O2 Amidinonorleucine ?
                                                  "X1186_522.38681_1217.11", #PC(O-C18:0/O-1:0)
                                                 #"X1694_166.08386_71.91", #C9H12NO2 Phenylalanine ?
                                                 #"X423_367.14706_131.15" , #C19H18N4O4 alkaloid
                                                  "X1724_180.09894_89.64",  #Salsolinol
                                                 #"X1645_146.11527_69.93", #C7H16NO2 Methyl-isoleucine
                                                 #"X1750_137.04436_94.87", #C5H5N4O Hypoxanthine
                                                 #"X454_144.07941_287.22" #C5H10N3O2 
                                                 "X402_204.13031_89.27", #Gly-Lys
                                                 "X394_226.11945_71.84")] #Gln-Pro


datatable(VIP_table_normalized)
dim(VIP_table_normalized)

colnames(VIP_table_normalized) <- c("Sample_names","Photic_zone", 
                                    "Lyso-PC(O-C18:0)", 
                                    "Rhodosamine", 
                                    #"C7H11N2O2 Cyclo-(Pro-Gly)?",
                                    #"C7H15N3O2 Amidinonorleucine?",
                                    "PC(O-C18:0/O-1:0)", 
                                    #"C9H12NO2 Phenylalanine?",
                                    #"C19H18N4O4 alkaloid",
                                    "Salsolinol",
                                    #"C7H16NO2 Methyl-isoleucine ?",
                                    #"C5H5N4O Hypoxanthine ?",
                                    #"C5H10N3O2", 
                                    "Gly-Lys",
                                    "Gln-Pro")


datatable(VIP_table_normalized)
dim(VIP_table_normalized)

data_long = melt(VIP_table_normalized, id.var=c("Sample_names","Photic_zone"))
datatable(data_long)

boxplot_vip_class = ggplot(data_long, aes(x = Photic_zone, y = value))
boxplot_vip_class = boxplot_vip_class + geom_point(aes(fill = Photic_zone), size = 2, alpha = 0.8, pch = 21)
boxplot_vip_class = boxplot_vip_class + geom_boxplot(aes(fill = Photic_zone), alpha = 0.4, size = 0.8)
boxplot_vip_class = boxplot_vip_class + scale_fill_manual(name = "Photic zone", values = color_vector_zone, labels = c("Mesophotic", "Upper rariphotic", "Lower rariphotic"))
boxplot_vip_class = boxplot_vip_class + theme_bw(base_size = 10)  + theme(legend.position="left")
boxplot_vip_class = boxplot_vip_class + facet_wrap(~variable)
#boxplot_vip_class = boxplot_vip_class + facet_grid(variable~Site, scales = "free")
boxplot_vip_class = boxplot_vip_class + theme(axis.title.x = element_blank())
boxplot_vip_class = boxplot_vip_class + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
boxplot_vip_class = boxplot_vip_class + ylab("Normalized concentration")
boxplot_vip_class
boxplot_vip_zone_demospongiae_main = boxplot_vip_class


ggsave(filename = "Boxplot_VIP_zone_demospongiae.pdf", 
       plot = boxplot_vip_zone_demospongiae, 
       device = "pdf" , 
       width = 35 , height = 25, units = "cm", 
       path = "./3_Photic_zone_difference_demospongiae_results")





## 2.4. PLSDA plot analysis with Hexactinellida #### 


PLSDA_score = read.csv(file = "./plsda_score_zone_Hexactinellida.csv", header = TRUE , sep = ";" , row.names = 1)
dim(PLSDA_score)
datatable(PLSDA_score)


plot_plsda = ggplot(PLSDA_score, aes(x = Comp.1, y = Comp.2))
plot_plsda = plot_plsda + geom_point(aes(fill = Photic_zone, pch = Class_host) , pch = 23, stroke = 1.5 , size = 7, alpha = 0.7, color = "black" , show.legend = T)
plot_plsda = plot_plsda + theme_bw(base_size = 20) 
plot_plsda = plot_plsda + scale_fill_manual(values = color_vector_zone, labels = c("Mesophotic", "Upper rariphotic", "Lower rariphotic"))
plot_plsda = plot_plsda +  xlab("Component 1 (5%)") + ylab("Component 2 (5%)")
plot_plsda = plot_plsda + guides(fill = guide_legend(override.aes = list(shape = 21)))
plot_plsda = plot_plsda + geom_text_repel(aes(label = Code), size = 6, box.padding = 0.8)
plot_plsda = plot_plsda +  labs(shape = "Sponge class", fill = "Photic zone")
plot_plsda

plot_plsda_hexactinellida_photic = plot_plsda

ggsave(filename = "PLSDA_plot_hexactinellida_photic.pdf", 
       plot = plot_plsda_hexactinellida_photic, 
       device = "pdf" , 
       width = 30 , height = 20, units = "cm", 
       path = "./4_Photic_zone_difference_hexactinellida_results")



## 2.5. Boxplots all VIPS zone for Hexactinellida ####

Data_table_normalized = read.csv(file = "./data_normalized_Hexactinellida.csv" , header = TRUE , sep = ";" , row.names = 1)

datatable(Data_table_normalized)
dim(Data_table_normalized)


VIP_table_normalized = Data_table_normalized[, c("Sample_names","Label", 
                                                 "X1983_482.24702_809.77", #Lyso-PS(C15:1)
                                                 "X1839_206.04277_275.93", #Xanthurenic acid
                                                 "X2039_570.35206_922.83", #Lyso-PC(C22:5)
                                                 "X1689_231.1668_71.99", #Leucylvaline
                                                 "X391_245.18216_71.64", #Leucylleucine
                                                 "X2129_534.35216_1068.35", #Lyso-PE(C22:2)
                                                 #"X1033_794.59263_1085.76", #PC(O-C38:5) ?
                                                 #"X1849_188.06838_303.97", #Imidazole derivative
                                                 #"X2060_530.31959_948.86", #Lyso-PE(C22:4)
                                                 #"X402_204.13031_89.27", #Glycyllysine
                                                 "X2031_570.35324_908.5")] #Lyso-PC(C22:5)
                                                 #"X2007_502.28847_877.11")]#Lyso-PE(20:4)****
                                                 #"X2222_810.59565_1160.99")] #PC(C38:4) 



datatable(VIP_table_normalized)
dim(VIP_table_normalized)

colnames(VIP_table_normalized) <- c("Sample_names","Photic_zone", 
                                    "Lyso-PS(C15:1)",
                                    "Xanthurenic acid",
                                    "Lyso-PC(C22:5)",
                                    "Leu-Val",
                                    "Leu-Leu",
                                    "Lyso-PE(C22:2)",
                                    #"PC(O-C38:5)",
                                    #"Imidazole derivative",
                                    #"Lyso-PE(C22:4)",
                                    #"Gly-Lys",
                                    "Lyso-PC(C22:5)")
                                    #"Lyso-PE(C20:4)")
                                    #"PC(C38:4)")


datatable(VIP_table_normalized)
dim(VIP_table_normalized)

data_long = melt(VIP_table_normalized, id.var=c("Sample_names","Photic_zone"))
datatable(data_long)

boxplot_vip_class = ggplot(data_long, aes(x = Photic_zone, y = value))
boxplot_vip_class = boxplot_vip_class + geom_point(aes(fill = Photic_zone), size = 2, alpha = 0.8, pch = 21)
boxplot_vip_class = boxplot_vip_class + geom_boxplot(aes(fill = Photic_zone), alpha = 0.4, size = 0.8)
boxplot_vip_class = boxplot_vip_class + scale_fill_manual(name = "Photic zone", values = color_vector_zone, labels = c("Mesophotic", "Upper rariphotic", "Lower rariphotic"))
boxplot_vip_class = boxplot_vip_class + theme_bw(base_size = 10)  + theme(legend.position="left")
boxplot_vip_class = boxplot_vip_class + facet_wrap(~variable)
#boxplot_vip_class = boxplot_vip_class + facet_grid(variable~Site, scales = "free")
boxplot_vip_class = boxplot_vip_class + theme(axis.title.x = element_blank())
boxplot_vip_class = boxplot_vip_class + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
boxplot_vip_class = boxplot_vip_class + ylab("Normalized concentration")
boxplot_vip_class
boxplot_vip_zone_hexactinellida_main = boxplot_vip_class


ggsave(filename = "Boxplot_VIP_zone_hexactinellida.pdf", 
       plot = boxplot_vip_zone_hexactinellida, 
       device = "pdf" , 
       width = 35 , height = 25, units = "cm", 
       path = "./4_Photic_zone_difference_hexactinellida_results")





# FIgure main MS with all boxplots

legend_boxplot_class = as_ggplot(get_legend(boxplot_vip_class_main+ theme_bw(base_size = 17)))
legend_boxplot_class

legend_boxplot_zone = as_ggplot(get_legend(boxplot_vip_zone_demospongiae_main + theme_bw(base_size = 17)))
legend_boxplot_zone


legend_all = plot_grid( legend_boxplot_class , legend_boxplot_zone  , labels=c("", "" ),rel_widths = c(1, 1 )  ,  ncol = 1, nrow = 2 ,label_size = 20)
legend_all


boxplot_vip_class_main_f = boxplot_vip_class_main + theme(legend.position="none")
boxplot_vip_class_main_f

boxplot_vip_zone_demospongiae_main_f = boxplot_vip_zone_demospongiae_main + theme(legend.position="none") 
boxplot_vip_zone_demospongiae_main_f

boxplot_vip_zone_hexactinellida_main_f = boxplot_vip_zone_hexactinellida_main + theme(legend.position="none") 
boxplot_vip_zone_hexactinellida_main_f

boxplot_vip_all_main = plot_grid( boxplot_vip_class_main_f , boxplot_vip_zone_demospongiae_main_f , boxplot_vip_zone_hexactinellida_main_f  , legend_all, labels=c("A", "B" ,"C", ""),rel_widths = c(1, 1 , 1, 2)  ,  ncol = 2, nrow = 2 ,label_size = 20)
boxplot_vip_all_main


ggsave(filename = "Boxplot_VIP_all.pdf", 
       plot = boxplot_vip_all_main, 
       device = "pdf" , 
       width = 19 , height = 18, units = "cm", 
       path = "./1_Class_difference_results")



