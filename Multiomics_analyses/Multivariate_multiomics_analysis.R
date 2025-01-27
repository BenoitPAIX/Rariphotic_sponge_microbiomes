# Multiomics data analyses  CAO project


#BiocManager::install('mixOmics')
library(ggplot2)
library(DT)
library(reshape2)
library(ggrepel)
library(vegan)
library(gridGraphics) 
library(mixOmics)
library(igraph)
library(DT)
library(cowplot)

setwd("your path file here")

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
                       "Dactylocalyx" = "dodgerblue2",
                       "Dictyoplax" = "skyblue",
                       "Heterotella" = "lightslateblue",
                       "Verrucocoeloidea" = "darkcyan")

# 1. Read the data tables ####

ASV_table = read.csv(file = "ASV_table_compo.csv" , header = TRUE , sep = ";" , row.names = 1)
ASV_table
dim(ASV_table)

LCMS_table = read.csv(file = "LCMS_table_normalized.csv" , sep = ";" , header = T , row.names = 1)
dim(LCMS_table)
datatable(LCMS_table)


METADATA_table = read.csv(file = "Metadata_table.csv" , header = TRUE , sep = ";" , row.names = 1)
dim(METADATA_table)
datatable(METADATA_table)
length(unique(list(row.names(ASV_table), 
                   row.names(LCMS_table), 
                   row.names(METADATA_table))))==1


# 2. Mantel tests ####

dist_matrix_asv = as.matrix(dist(ASV_table), method="bray")
dist_matrix_asv

dist_matrix_lcms = as.matrix(dist(LCMS_table), method="euclidean")
dist_matrix_lcms

data_mantel = mantel(dist_matrix_asv, dist_matrix_lcms, method = "spearman", permutations = 9999, na.rm = TRUE)
data_mantel



## 2.4 Procrust analysis ####
dist_matrix_asv
nmds_asv <- metaMDS(ASV_table, distance = "bray")
nmds_asv
plot(nmds_asv)

data.scores = as.data.frame(scores(nmds_asv)$sites)
data.scores.meta  = cbind(data.scores, METADATA_table)
data.scores.meta

plot_nmds = ggplot(data.scores.meta, aes(x = NMDS1, y = NMDS2))
plot_nmds = plot_nmds + geom_point(aes(fill = Photic_zone, pch = Class_host), size = 7, alpha = 0.8)
plot_nmds = plot_nmds + scale_shape_manual(values = c(22,23))
plot_nmds = plot_nmds + scale_fill_manual(values = color_vector_zone, labels = c("Mesophotic", "Upper rariphotic", "Lower rariphotic"))
plot_nmds = plot_nmds + theme_bw(base_size = 25)+ theme(plot.title = element_text(size=22, face="bold"), axis.text=element_text(size=12),axis.title=element_text(size=14))
plot_nmds = plot_nmds + guides(fill = guide_legend(override.aes = list(shape = 21, size = 7)))
plot_nmds = plot_nmds + annotate("text", label="2D Stress = 0.16")
plot_nmds = plot_nmds + labs(shape = "Sponge class", fill = "Photic zone")
plot_nmds
plot_nmds_classzone = plot_nmds



pca_lcms <- cmdscale(dist_matrix_lcms)
data.scores.lcms = as.data.frame(scores(pca_lcms))
data.scores.lcms.meta = cbind(data.scores.lcms, METADATA_table)
data.scores.lcms.meta

plot_pca = ggplot(data.scores.lcms.meta, aes(x = Dim1, y = Dim2))
plot_pca = plot_pca + geom_point(aes(fill = Photic_zone, pch = Class_host), size = 7, alpha = 0.8)
plot_pca = plot_pca + scale_shape_manual(values = c(22,23))
plot_pca = plot_pca + scale_fill_manual(values = color_vector_zone, labels = c("Mesophotic", "Upper rariphotic", "Lower rariphotic"))
plot_pca = plot_pca + theme_bw(base_size = 25)+ theme(plot.title = element_text(size=22, face="bold"), axis.text=element_text(size=12),axis.title=element_text(size=14))
plot_pca = plot_pca + guides(fill = guide_legend(override.aes = list(shape = 21, size = 7)))
plot_pca = plot_pca + annotate("text", label="2D Stress = 0.16")
plot_pca = plot_pca + labs(shape = "Sponge class", fill = "Photic zone")
plot_pca
plot_pca_classzone = plot_pca


pro <- procrustes( X = pca_lcms,Y = nmds_asv, symmetric = TRUE)
pro_x = as.data.frame(pro$X)
pro_x
pro_y = as.data.frame(pro$Yrot)
pro_y

data.scores.pro.meta = cbind(pro_x, pro_y, METADATA_table)
data.scores.pro.meta$Genus_host

plot_pro = ggplot(data.scores.pro.meta)
plot_pro = plot_pro + geom_segment(aes(x=Dim1,y=Dim2,xend=V1,yend=V2),colour = "black", arrow=arrow(length=unit(0.3,"cm")), alpha = 0.6, size = 1)
plot_pro = plot_pro + geom_point(aes(x = Dim1, y = Dim2, fill = Photic_zone, pch = Class_host), size = 4, alpha = 0.7, color = "#950606", stroke  = 2)
plot_pro = plot_pro + geom_point(aes(x = V1, y = V2, fill = Photic_zone, pch = Class_host), size = 4, alpha = 0.7, color = "#06402B", stroke  = 2)
plot_pro = plot_pro + scale_shape_manual(values = c(22,23))
plot_pro = plot_pro + scale_fill_manual(values = color_vector_zone, labels = c("Mesophotic", "Upper rariphotic", "Lower rariphotic"))
#plot_pro = plot_pro + scale_colour_manual(values = color_vector_zone, labels = c("Mesophotic", "Upper rariphotic", "Lower rariphotic"))
plot_pro = plot_pro + theme_bw(base_size = 15)+ theme(plot.title = element_text(size=22, face="bold"), axis.text=element_text(size=12),axis.title=element_text(size=14))
plot_pro = plot_pro + guides(fill = guide_legend(override.aes = list(shape = 21, size = 7)))
plot_pro = plot_pro + labs(shape = "Sponge class", fill = "Photic zone")
plot_pro
plot_pro_classzone = plot_pro

plot_pro = ggplot(data.scores.pro.meta)
plot_pro = plot_pro + geom_segment(aes(x=Dim1,y=Dim2,xend=V1,yend=V2),colour = "black", arrow=arrow(length=unit(0.3,"cm")), alpha = 0.6, size = 1)
plot_pro = plot_pro + geom_point(aes(x = Dim1, y = Dim2, fill = Genus_host, pch = Class_host), size = 4, alpha = 0.7, color = "#950606", stroke  = 2)
plot_pro = plot_pro + geom_point(aes(x = V1, y = V2, fill = Genus_host, pch = Class_host), size = 4, alpha = 0.7, color = "#06402B", stroke  = 2)
plot_pro = plot_pro + scale_shape_manual(values = c(22,23))
plot_pro = plot_pro + scale_fill_manual(values = color_vector_genus)
#plot_pro = plot_pro + scale_colour_manual(values = color_vector_zone, labels = c("Mesophotic", "Upper rariphotic", "Lower rariphotic"))
plot_pro = plot_pro + theme_bw(base_size = 15)+ theme(plot.title = element_text(size=22, face="bold"), axis.text=element_text(size=12),axis.title=element_text(size=14))
plot_pro = plot_pro + guides(fill = guide_legend(override.aes = list(shape = 21, size = 7)))
plot_pro = plot_pro + labs(shape = "Sponge class", fill = "Sponge genus") 
plot_pro = plot_pro + theme(axis.title.y = element_blank())
plot_pro_classgenus = plot_pro
plot_pro_classgenus

biplot_pro = plot_grid( plot_pro_classzone+ theme(legend.position="none") , plot_pro_classgenus+ theme(legend.position="none") , labels=c("A", ""),rel_widths = c(1, 1 )  ,  ncol = 2, nrow = 1 ,label_size = 20)
biplot_pro

triplot_pro_alpha = plot_grid( plot_alpha_multiomic_all , plot_pro_classzone+ theme(legend.position="none") , plot_pro_classgenus+ theme(legend.position="none") , labels=c("A", "B", ""),rel_widths = c(1, 1 , 1)  ,  ncol = 3, nrow = 1 ,label_size = 20)
triplot_pro_alpha

ggsave(filename = "Triplot_pro_alpha_multiomics_all.pdf", 
       plot = triplot_pro_alpha, 
       device = "pdf" , 
       width = 40 , height = 15, units = "cm")


# 3. mixOmics analyses ####


data = list( datalcms = LCMS_table , dataASV = ASV_table)

datagroup = METADATA_table

datagroup
datagroup["Class_zone"] <- paste(datagroup$Class_host, datagroup$Photic_zone)
datatable(datagroup)

design = matrix(1, ncol = length(data), nrow = length(data),dimnames = list(names(data), names(data)))
diag(design) =  0

ncomp = c(8)
# kept ncomp = c(4)


#block.plsda = block.splsda(X = data, Y = datagroup$Class_zone, ncomp = ncomp, design = design, near.zero.var = TRUE)
#block.plsda

# run component number tuning with repeated CV
#perf.diablo = perf(block.plsda, validation = 'Mfold',  folds = 10, nrepeat = 10) 

#plot(perf.diablo) # plot output of tuning

#ncomp = perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"] 
#ncomp

#perf.diablo$choice.ncomp$WeightedVote 

ncomp = 8


#tune.TCGA = tune.block.splsda(X = data, Y = datagroup$Class_zone, ncomp = ncomp, 
#                              test.keepX = test.keepX, design = design,
#                              validation = 'Mfold', folds = 10, nrepeat = 1,
#                              dist = "centroids.dist", near.zero.var = TRUE)


#list.keepX = tune.TCGA$choice.keepX # set the optimal values of features to retain
#list.keepX

## results of the tuning : 
#$datalcms
#[1] 25  5  5  5  5  5

#$dataASV
#[1]  6 20 30  5 30  5

list.keepX = list(datalcms = c(25,5,5,5,5,5), dataASV = c(6,20,30,5,30,5))
list.keepX


datagroup$Photic_zone
final.diablo.model = block.splsda(X = data, Y = datagroup$Class_zone, ncomp = ncomp, 
                                  keepX = list.keepX, design = design, near.zero.var = TRUE)

final.diablo.model

length(list.keepX)
#block.plsda = block.splsda(X = data, Y = datagroup$Class_zone, ncomp = ncomp, keepX = list.keepX, design = design, near.zero.var = TRUE)
#block.plsda
#kept block.plsda = block.splsda(X = data, Y = datagroup$Class_zone, ncomp = ncomp, keepX = list.keepX, design = design, near.zero.var = TRUE)



plotArrow(final.diablo.model, ind.names = FALSE, legend = TRUE)


plotDiablo(final.diablo.model , ncomp = 1)

plotIndiv(final.diablo.model, ind.names = FALSE, legend = TRUE, 
          title = 'DIABLO Sample Plots')


# circosPlot(block.plsda, cutoff = 0.8, line = TRUE, color.blocks= color_node, color.cor = c("chocolate3","grey20"), size.labels = 1.5)

#parameters saved for paper : network.result = network(block.plsda , cutoff = 0.75 , blocks = c(1,2,3) , color.node = color_node , cex.node.name = 0.3 , keysize = c(1,1) , interactive = T)

#network.result = network(block.plsda , cutoff = 0.75 , blocks = c(1,2) , color.node = color_node , cex.node.name = 0.3 , keysize = c(1,1) , interactive = F)
#network.result

#write.graph(network.result$gR, file = "Correlation_network.graphml", format = "graphml")
group_colors <- c("Demospongiae 1.Mesophotic" = "yellow1", 
                  "Demospongiae 2.Upper rariphotic" = "aquamarine", 
                  "Demospongiae 2.Upper rariphotic" = "royalblue", 
                  "Hexactinellida 2.Upper rariphotic" = "aquamarine", 
                  "Hexactinellida 3.Lower rariphotic" = "royalblue")  

# Définir un vecteur de couleurs pour chaque type de données
color_code <- c("dataLCMS" = "#2ca02c",   # Vert pour les métabolites en LCMS
                "dataASV" = "#ff7f0e")  # Orange pour les séquences ASV 16S


heatmap.block.plsda= cimDiablo(final.diablo.model, comp = 1:5, 
                                 margins = c(10, 25), 
                                 color.blocks = color_code,  
                                 color.Y = group_colors,#group_colors
                                 trim = 3, 
                                 transpose = T, 
                                 legend.position = "none",
                                 keysize = c(10, 1), 
                                 keysize.label = 1,
                                 row.cex = 0.7, 
                                 col.cex = 1, 
                                 size.legend = 1, 
                                 row.names = datagroup$Sample_names)

heatmap.block.plsda$col.names <- gsub("^X", "", heatmap.block.plsda$col.names)
heatmap.block.plsda$col.names

Annotation_table = read.csv(file = "Annotation_table.csv" , header = TRUE , sep = ";" )
selected.variables <- heatmap.block.plsda$col.names

real_names <- Annotation_table$Final_name[match(selected.variables, Annotation_table$X)]
real_names


heatmap.block.plsda$col.names



heatmap.block.plsda= cimDiablo(final.diablo.model, comp = 1:5, 
                               margins = c(10, 25), 
                               color.blocks = color_code,  
                               color.Y = group_colors,#group_colors
                               trim = 3, 
                               transpose = F, 
                               legend.position = "right",
                               keysize = c(10, 1), 
                               keysize.label = 1,
                               row.cex = 0.7, 
                               col.cex = 1, 
                               size.legend = 1, 
                               row.names = datagroup$Sample_names,
                               col.names = real_names)
