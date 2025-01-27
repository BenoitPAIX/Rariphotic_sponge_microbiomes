# Metabarcoding analysis - rariphotic sponge microbiome project #

# 1. Getting ready to start with R ####

## 1.2 Package installation ####
#install.packages("randomcoloR")
#install.packages("ggplot2")
#install.packages("stringi")
#install.packages("rlang")

#if (!requireNamespace("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}
#BiocManager::install("phyloseq")
#BiocManager::install("microbiome")
#BiocManager::install("decontam")
#install.packages("ggordiplots")
#install.packages("remotes")
#remotes::install_github("vmikk/metagMisc")
#install.packages("iNEXT")

## 1.3 Packages loading ####

library(stringi)
library(vegan)
library(Rcpp)
library(ggplot2)
library(randomcoloR)
library(rlang)
library(phyloseq)
library(ape)
library(dplyr)
library(agricolae) #for statistical tests
library(RVAideMemoire) #for statistical tests
library(microbiome)
library(hrbrthemes)
library(cowplot) #for figures
#library(ggthemes)
library(RColorBrewer) #for figures
library(decontam) #for decontamination of the dataset
library(ggrepel) #for figures
#library(microbiomeMarker) #for lefse analysis
library(multcompView)
library(rcompanion)
library(pairwiseAdonis) 
library(DT)
library(reshape2)
library(ggh4x)
library(tidyverse)
library(ggordiplots)
library(ggpubr)
library(cowplot)
library(metacoder)
library(microbiomeutilities)
#library(metagMisc)## check this package for later, many miscellaneous function that can be usefull
# library(iNEXT) not used yet
library(cowplot)
library(scales)

## 1.4 Setting and preparing your working directory ####

setwd("your path file here")
#dir.create("1_Data_prep_results")
#dir.create("2_Alpha_div_results")
#dir.create("3_Beta_div_results")
#dir.create("4_Compositional_results")
#dir.create("5_Differential_results")


# 2. Preparation of the 16S metabarcoding dataset for a ready-to-go analysis ####


## 2.1. Import the dataset ####


ASV_table = read.csv(file = "ASV_table.csv" , header = TRUE , sep = ";" , row.names = 1)
dim(ASV_table)
#datatable(ASV_table)

TAX_table = read.csv(file = "Taxonomy_table.csv" , sep = ";" , header = T , row.names = 1)
dim(TAX_table)
#datatable(TAX_table)

META_table = read.csv(file = "Metadata_table.csv" , header = TRUE , sep = ";" , row.names = 1)
dim(META_table)
#datatable(META_table)

TAX_table = as.data.frame(TAX_table)
TAX_table

TAX_table.v2 <- TAX_table

TAX_table.v2$Kingdom <- TAX_table.v2$Kingdom
TAX_table.v2$Phylum <- paste(TAX_table.v2$Kingdom,"|" ,TAX_table$Phylum)
TAX_table.v2$Class <- paste(TAX_table.v2$Phylum,"|" ,TAX_table$Class)
TAX_table.v2$Order <- paste(TAX_table.v2$Class ,"|" ,TAX_table$Order)
TAX_table.v2$Family <- paste(TAX_table.v2$Order,"|" ,TAX_table$Family)
TAX_table.v2$Genus <- paste(TAX_table.v2$Family,"|" ,TAX_table$Genus)


head(TAX_table.v2)

## 2.2. Create your main phyloseq object (output: physeq_raw) ####

ASV_phylo = otu_table(ASV_table, taxa_are_rows = TRUE)
dim(ASV_phylo)


TAX_table = as.matrix(TAX_table.v2)

TAX_phylo = tax_table(TAX_table)
dim(TAX_phylo)

META_phylo = sample_data(META_table)
dim(META_phylo)

physeq_raw = phyloseq(ASV_phylo, TAX_phylo, META_phylo)
physeq_raw


## 2.2. Decontaminate your dataset using the decontam package (output: physeq_decontam)####


#inspection of the library sizes
df <- as.data.frame(sample_data(physeq_raw)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(physeq_raw)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))

plot_lib_size = ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control))
plot_lib_size = plot_lib_size + geom_point()
plot_lib_size

ggsave(filename = "Plot_lib_size.pdf", 
       plot = plot_lib_size, 
       device = "pdf" , 
       width = 15 , height = 10, units = "cm", 
       path = "./1_Data_prep_results")


#identification of the contaminants (prevalence method)

sample_data(physeq_raw)$is.neg <- sample_data(physeq_raw)$Sample_or_Control == "Control"

contamdf.prev05 <- isContaminant(physeq_raw, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant)


contamdf.prev05 = cbind(as.data.frame(tax_table(physeq_raw)) , contamdf.prev05)
contamdf.prev05

write.csv(contamdf.prev05, file.path("./1_Data_prep_results" , "Contamination_table_prev05.csv"))


# Make phyloseq object of presence-absence in negative controls and true samples

ps.pa <- transform_sample_counts(physeq_raw, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "True sample", ps.pa)

# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev05$contaminant)

plot_prevalence = ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) 
plot_prevalence = plot_prevalence + geom_point() 
plot_prevalence = plot_prevalence + xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
plot_prevalence


ggsave(filename = "Plot_prevalence.pdf", 
       plot = plot_prevalence, 
       device = "pdf" , 
       width = 15 , height = 10, units = "cm", 
       path = "./1_Data_prep_results")



#create a phyloseq_decontam object without the contaminants
physeq_decontam = prune_taxa(!contamdf.prev05$contaminant, physeq_raw)
physeq_decontam


## 2.3. Filter additional non-prokaryotic taxa by removing ASVs with 16S from Eukaryotes, chloroplasts, and mitochondria ####


physeq_filtered = subset_taxa(physeq_decontam, 
                              (Kingdom != "Eukaryota" | is.na(Kingdom)) &
                              (Order!= "Bacteria | Cyanobacteria | Cyanobacteriia | Chloroplast"| is.na(Order)) &
                              (Family != "Bacteria | Proteobacteria | Alphaproteobacteria | Rickettsiales | Mitochondria"| is.na(Family)))

physeq_filtered


## 2.4.Checking the rarefaction curves and removing samples with not enough reads ####


as.data.frame(t(otu_table(physeq_filtered)))

rarecurve(as.data.frame(t(otu_table(physeq_filtered))), step = 20, cex = 0.5)
#save manually the plot to ./1_Data_prep_results


## 2.5. Removing unwanted samples for the analysis (control samples and those below the threshold) ####

#Set a minimum reads threshold based on the rarefaction curves results

min_reads_threshold = 10000


sample_sums(physeq_filtered)


physeq_presubsampled = prune_samples(sample_sums(physeq_filtered)>=min_reads_threshold , physeq_filtered)
physeq_presubsampled

physeq_subsampled = subset_samples(physeq_presubsampled, ID_Samples != "CAO_420")
physeq_subsampled

saveRDS(physeq_subsampled, "Physeq_subsampled.RDS")


physeq_subsampled = readRDS("Physeq_subsampled.RDS") ################################## <--- start again here !
physeq_subsampled

### re-doing the rarefaction curves with the physeq_subsampled object

rarecurve(as.data.frame(t(otu_table(physeq_subsampled))), step = 20, cex = 0.5)

data_rarecurve = rarecurve(as.data.frame(t(otu_table(physeq_subsampled))), step=50,  tidy = TRUE)
data_rarecurve

write.csv(data_rarecurve, "Data_rarecurve.csv")


plot_rarecurve = ggplot(data_rarecurve, aes(x = Sample   , y = Species, colour = Site))
plot_rarecurve = plot_rarecurve + geom_line(linewidth = 0.5, alpha = 0.7)
plot_rarecurve = plot_rarecurve + theme_bw(15)
plot_rarecurve = plot_rarecurve + scale_color_manual(values = rep(c("black"),times=58))
plot_rarecurve = plot_rarecurve +  scale_x_continuous(labels = comma_format(big.mark = ".",
                                                                            decimal.mark = ","))
plot_rarecurve = plot_rarecurve + theme(legend.position = "none")
plot_rarecurve = plot_rarecurve + xlab("Sample size (reads number)") + ylab("ASVs number")
plot_rarecurve 

ggsave(filename = "Plot_rarefaction_curves_v2.pdf", 
       plot = plot_rarecurve, 
       device = "pdf" , 
       width = 20 , height = 20, units = "cm", 
       path = "./1_Data_prep_results")


#write.csv(as.data.frame(sample_data(physeq_subsampled)), "Metadata_subsampled.csv")

# 4. Alpha diversity analyses ####

## 4.1. Create a phyloseq rarefied objects specifically for alpha-div analyses ####

# Create a phyloseq rarefied objects for alpha div analyses

physeq_rarefied = rarefy_even_depth(physeq_subsampled)
physeq_rarefied

sample_sums(physeq_rarefied)


## 4.2. Consider the factors of comparison for this analysis and prepare your color vectors ####
metadata_subsampled = sample_data(physeq_subsampled)
datatable(metadata_subsampled)


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


## 4.3. Create a data frame with the results of the alpha-diversity indices and the metadata ####
data_alpha = estimate_richness(physeq_rarefied , measures = c("Observed", "Shannon", "Chao1"))
data_alpha

Pielou = data_alpha$Shannon / log(data_alpha$Observed)
Pielou
data_alpha = cbind(sample_data(physeq_rarefied), data_alpha , Pielou)

datatable(data_alpha)
data_alpha
#write.csv(data_alpha, "Alpha_table.csv")

## 4.4. Generate a single figure for the 3 alpha-diversity indices together ####

data_alpha_subset = data_alpha[, c("Class_host", "Photic_zone", "Shannon" ,"Chao1","Pielou")]
data_alpha_subset
data_alpha_subset_long = melt(data_alpha_subset, id.var=c("Class_host", "Photic_zone"))
colnames(data_alpha_subset_long) <- c("Class_host", "Photic_zone", "Index", "Values")
data_alpha_subset_long

datatable(data_alpha_subset_long)

plot_alpha = ggplot(data_alpha_subset_long, aes(Photic_zone ,Values, fill = Photic_zone, pch = Class_host)) + guides(fill=guide_legend(title="Photic zone"))
plot_alpha = plot_alpha + geom_boxplot(alpha = 0.8, size = 1) + facet_grid( Index ~ Class_host , scales="free", space = "free_x") 
plot_alpha = plot_alpha + geom_point(size = 2, alpha = 0.8, stroke = 1)
plot_alpha = plot_alpha + theme_bw(base_size = 18) 
plot_alpha = plot_alpha + theme(legend.position="bottom",
                                      legend.box = "vertical",
                                      strip.text.x = element_text(size = 18),
                                      strip.text.y = element_text(size = 18),
                                      legend.title = element_text(size=18), 
                                      legend.text = element_text(size=18))
plot_alpha = plot_alpha + theme(axis.title.x = element_blank(),axis.title.y = element_blank())
plot_alpha = plot_alpha + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
plot_alpha = plot_alpha + scale_fill_manual(values = color_vector_zone, labels = c("Mesophotic", "Upper rariphotic", "Lower rariphotic"))
plot_alpha = plot_alpha + scale_shape_manual(values = c(21, 23))
plot_alpha = plot_alpha + labs(shape = "Sponge class", fill = "Photic zone")
plot_alpha 


ggsave(filename = "Plot_alpha_class_photic.pdf", 
       plot = plot_alpha, 
       device = "pdf" , 
       width = 20 , height = 25, units = "cm", 
       path = "./2_Alpha_div_results")



## 4.5. Test the normality of your distribution for each index ####

shapiro_data_shannon = shapiro.test(data_alpha$Shannon)
shapiro_data_shannon

shapiro_data_chao1 = shapiro.test(data_alpha$Chao1)
shapiro_data_chao1

shapiro_data_pielou = shapiro.test(data_alpha$Pielou)
shapiro_data_pielou



shapiro_data_alpha <- matrix(nrow = 3 ,  ncol=2, byrow=TRUE)
colnames(shapiro_data_alpha) = c("W","p-value")
rownames(shapiro_data_alpha) = c("Shannon","Chao1","Pielou")

shapiro_data_alpha[,1] <- c(shapiro_data_shannon$statistic,
                            shapiro_data_chao1$statistic,
                            shapiro_data_pielou$statistic)

shapiro_data_alpha[,2]<- c(shapiro_data_shannon$p.value,
                           shapiro_data_chao1$p.value,
                           shapiro_data_pielou$p.value)

shapiro_data_alpha
datatable(shapiro_data_alpha)

write.csv(shapiro_data_alpha, file.path("./2_Alpha_div_results" , "Shapiro_data_alpha.csv"))


## 4.6. Test the differences of alpha-diversity according to your factors with analyses of variances ####

data_alpha["Class_zone"] <- paste(data_alpha$Class_host, data_alpha$Photic_zone)
datatable(data_alpha)


data_kruskal_Shannon_classzone = kruskal.test(Shannon ~ Class_zone, data_alpha)
data_kruskal_Shannon_classzone

data_kruskal_Chao1_classzone = kruskal.test(Chao1 ~ Class_zone, data_alpha)
data_kruskal_Chao1_classzone

data_kruskal_Pielou_classzone = kruskal.test(Pielou ~ Class_zone, data_alpha)
data_kruskal_Pielou_classzone 


data_kruskal_alpha_classzone  <- matrix(nrow = 3 ,  ncol=3, byrow=TRUE)
colnames(data_kruskal_alpha_classzone) = c("Chi-square","Df","p-value")
rownames(data_kruskal_alpha_classzone) = c("Shannon","Chao1","Pielou")

data_kruskal_alpha_classzone

data_kruskal_alpha_classzone[,1] <- c(data_kruskal_Shannon_classzone$statistic,
                                  data_kruskal_Chao1_classzone$statistic,
                                  data_kruskal_Pielou_classzone$statistic)

data_kruskal_alpha_classzone[,2]<- c(data_kruskal_Shannon_classzone$parameter,
                                 data_kruskal_Chao1_classzone$parameter,
                                 data_kruskal_Pielou_classzone$parameter)

data_kruskal_alpha_classzone[,3]<- c(data_kruskal_Shannon_classzone$p.value,
                                 data_kruskal_Chao1_classzone$p.value,
                                 data_kruskal_Pielou_classzone$p.value)

data_kruskal_alpha_classzone
datatable(data_kruskal_alpha_classzone)

write.csv(data_kruskal_alpha_classzone, file.path("./2_Alpha_div_results" , "Data_kruskal_alpha_classzone.csv"))

## 4.7. Test the differences between groups with pairwise comparisons ####

### 4.7.1. Wilcoxon test Shannon ####

data_wilcox_Shannon_classzone = pairwise.wilcox.test(data_alpha$Shannon, data_alpha$Class_zone,
                                             p.adjust.method = "bonf")
data_wilcox_Shannon_classzone

data_wilcox_shannon_classzone_full = fullPTable(data_wilcox_Shannon_classzone$p.value)
data_wilcox_shannon_classzone_full

indices_wilcox_shannon_classzone  = multcompLetters(data_wilcox_shannon_classzone_full,
                                               compare="<",
                                               threshold=0.05,
                                               Letters=letters,
                                               reversed = FALSE)
indices_wilcox_shannon_classzone$Letters

data_wilcox_shannon_classzone_full2 = cbind(data_wilcox_shannon_classzone_full, indices_wilcox_shannon_classzone$Letters)
data_wilcox_shannon_classzone_full2


datatable(data_wilcox_shannon_classzone_full2)


write.csv(data_wilcox_shannon_classzone_full2, file.path("./2_Alpha_div_results" , "Data_wilcox_shannon_classzone.csv"))


### 4.7.2. Wilcoxon test Chao1 ####

data_wilcox_Chao1_classzone = pairwise.wilcox.test(data_alpha$Chao1, data_alpha$Class_zone,
                                                     p.adjust.method = "bonf")
data_wilcox_Chao1_classzone

data_wilcox_Chao1_classzone_full = fullPTable(data_wilcox_Chao1_classzone$p.value)
data_wilcox_Chao1_classzone_full

indices_wilcox_Chao1_classzone  = multcompLetters(data_wilcox_Chao1_classzone_full,
                                                    compare="<",
                                                    threshold=0.05,
                                                    Letters=letters,
                                                    reversed = FALSE)
indices_wilcox_Chao1_classzone$Letters

data_wilcox_Chao1_classzone_full2 = cbind(data_wilcox_Chao1_classzone_full, indices_wilcox_Chao1_classzone$Letters)
data_wilcox_Chao1_classzone_full2


datatable(data_wilcox_Chao1_classzone_full2)


write.csv(data_wilcox_Chao1_classzone_full2, file.path("./2_Alpha_div_results" , "Data_wilcox_Chao1_classzone.csv"))


### 4.7.3. Wilcoxon test Pielou ####

data_wilcox_Pielou_classzone = pairwise.wilcox.test(data_alpha$Pielou, data_alpha$Class_zone,
                                                   p.adjust.method = "bonf")
data_wilcox_Pielou_classzone

data_wilcox_Pielou_classzone_full = fullPTable(data_wilcox_Pielou_classzone$p.value)
data_wilcox_Pielou_classzone_full

indices_wilcox_Pielou_classzone  = multcompLetters(data_wilcox_Pielou_classzone_full,
                                                  compare="<",
                                                  threshold=0.05,
                                                  Letters=letters,
                                                  reversed = FALSE)
indices_wilcox_Pielou_classzone$Letters

data_wilcox_Pielou_classzone_full2 = cbind(data_wilcox_Pielou_classzone_full, indices_wilcox_Pielou_classzone$Letters)
data_wilcox_Pielou_classzone_full2


datatable(data_wilcox_Pielou_classzone_full2)


write.csv(data_wilcox_Pielou_classzone_full2, file.path("./2_Alpha_div_results" , "Data_wilcox_Pielou_classzone.csv"))




# 5. Beta diversity analyses ####


## 5.1. Create a phyloseq compositional objects specifically for beta-div and compositional analyses ####

physeq_compo = transform(physeq_subsampled, "compositional")
physeq_compo


#write.csv(as.data.frame(otu_table(physeq_compo)), "ASV_table_subsampled.csv")

## 5.2. Consider the factors of comparison for this analysis and prepare your color and shape vectors ####


metadata_subsampled = sample_data(physeq_subsampled)
datatable(metadata_subsampled)

color_vector_class 
color_vector_zone 
color_vector_genus 


## 5.3. NMDS analysis ####


nmds = ordinate(physeq_compo, "NMDS", "bray")
nmds$points
nmds$stress


### 5.3.1 NMDS plot with photic zone and class ####

data_nmds_samples <- plot_ordination(physeq_compo, nmds, type = "Samples", justDF = TRUE)
data_nmds_samples


plot_nmds = ggplot(data_nmds_samples, aes(x = NMDS1, y = NMDS2))
plot_nmds = plot_nmds + geom_point(aes(fill = Photic_zone, pch = Class_host), size = 7, alpha = 0.8)
plot_nmds = plot_nmds + scale_shape_manual(values = c(22,23))
plot_nmds = plot_nmds + scale_fill_manual(values = color_vector_zone, labels = c("Mesophotic", "Upper rariphotic", "Lower rariphotic"))
plot_nmds = plot_nmds + theme_bw(base_size = 25)+ theme(plot.title = element_text(size=22, face="bold"), axis.text=element_text(size=12),axis.title=element_text(size=14))
plot_nmds = plot_nmds + guides(fill = guide_legend(override.aes = list(shape = 21, size = 7)))
plot_nmds = plot_nmds + annotate("text", label="2D Stress = 0.16", x = -0.2, y = 0.15, size = 5)
plot_nmds = plot_nmds + labs(shape = "Sponge class", fill = "Photic zone")
plot_nmds
plot_nmds_classzone = plot_nmds


ggsave(filename = "Plot_NMDS_class_zone.pdf", 
       plot = plot_nmds_classzone, 
       device = "pdf" , 
       width = 30 , height = 20, units = "cm", 
       path = "./3_Beta_div_results")

### 5.3.2 NMDS plot with genus and class ####

plot_nmds = ggplot(data_nmds_samples, aes(x = NMDS1, y = NMDS2))
plot_nmds = plot_nmds + geom_point(aes(fill = Genus_nicole, pch = Class_host), size = 7, alpha = 0.8)
plot_nmds = plot_nmds + scale_shape_manual(values = c(22,23))
plot_nmds = plot_nmds + scale_fill_manual(values = color_vector_genus)
plot_nmds = plot_nmds + theme_bw(base_size = 20)+ theme(plot.title = element_text(size=22, face="bold"), axis.text=element_text(size=12),axis.title=element_text(size=14))
plot_nmds = plot_nmds + guides(fill = guide_legend(override.aes = list(shape = 21, size = 7)))
plot_nmds = plot_nmds + annotate("text", label="2D Stress = 0.16", x = -0.2, y = 0.15, size = 5)
plot_nmds = plot_nmds + labs(shape = "Sponge class", fill = "Sponge genus")
plot_nmds = plot_nmds + geom_text_repel(aes(label = Sample_names_nicole, color = Genus_nicole), size = 6, box.padding = 0.8)
plot_nmds = plot_nmds + scale_color_manual(values = color_vector_genus, guide = 'none')
plot_nmds = plot_nmds + guides(fill = guide_legend(ncol=2, override.aes = list(shape = 21)))
plot_nmds_nicole = plot_nmds
plot_nmds_nicole

plot_nmds = ggplot(data_nmds_samples, aes(x = NMDS1, y = NMDS2))
plot_nmds = plot_nmds + geom_point(aes(fill = Genus__celso, pch = Class_host), size = 7, alpha = 0.8)
plot_nmds = plot_nmds + scale_shape_manual(values = c(22,23))
plot_nmds = plot_nmds + scale_fill_manual(values = color_vector_genus)
plot_nmds = plot_nmds + theme_bw(base_size = 20)+ theme(plot.title = element_text(size=22, face="bold"), axis.text=element_text(size=12),axis.title=element_text(size=14))
plot_nmds = plot_nmds + guides(fill = guide_legend(override.aes = list(shape = 21, size = 7)))
plot_nmds = plot_nmds + annotate("text", label="2D Stress = 0.16", x = -0.2, y = 0.15, size = 5)
plot_nmds = plot_nmds + labs(shape = "Sponge class", fill = "Sponge genus")
plot_nmds = plot_nmds + geom_text_repel(aes(label = Sample_names_celso, color = Genus__celso), size = 6, box.padding = 0.8)
plot_nmds = plot_nmds + scale_color_manual(values = color_vector_genus, guide = 'none')
plot_nmds = plot_nmds + guides(fill = guide_legend(ncol=2, override.aes = list(shape = 21)))
plot_nmds_celso = plot_nmds 
plot_nmds_celso

plot_nmds_classgenus = plot_nmds_celso
plot_nmds_classgenus

ggsave(filename = "Plot_NMDS_classgenus.pdf", 
       plot = plot_nmds_classgenus, 
       device = "pdf" , 
       width = 30 , height = 20, units = "cm", 
       path = "./3_Beta_div_results")


### 5.3.3 same NMDS analysis with ordisurf for the depth ####


ordi =  ordisurf(nmds,data_nmds_samples$Depth_m, plot = FALSE, bs = "ds")
ordi$grid
ordi.grid <- ordi$grid #extracts the ordisurf object
ordi.grid
str(ordi.grid) #it's a list though - cannot be plotted as is
ordi.mite <- expand.grid(x = ordi.grid$x, y = ordi.grid$y) #get x and ys
ordi.mite$z <- as.vector(ordi.grid$z) #unravel the matrix for the z scores
ordi.mite.na <- data.frame(na.omit(ordi.mite)) #gets rid of the nas
ordi.mite.na



plot_nmds = ggplot(data_nmds_samples, aes(x = NMDS1, y = NMDS2))
plot_nmds = plot_nmds + geom_point(aes(fill = Photic_zone, pch = Class_host), size = 7, alpha = 0.8)
plot_nmds = plot_nmds + scale_shape_manual(values = c(22,23))
plot_nmds = plot_nmds + scale_fill_manual(values = color_vector_zone, labels = c("Mesophotic", "Upper rariphotic", "Lower rariphotic"))
plot_nmds = plot_nmds + theme(plot.title = element_text(size=22, face="bold"), axis.text=element_text(size=12),axis.title=element_text(size=14))
plot_nmds = plot_nmds + guides(fill = guide_legend(override.aes = list(shape = 22, size = 7)))
plot_nmds = plot_nmds + annotate("text", label="2D Stress = 0.16", x = -0.2, y = 0.15, size = 5)
plot_nmds = plot_nmds + stat_contour(data = ordi.mite.na, aes(x = x, y = y, z = z, colour = after_stat(level)))
plot_nmds =  plot_nmds + scale_colour_gradient2(low = "gold", mid = "aquamarine", high = "darkblue", midpoint = 200, guide=guide_colourbar(reverse = TRUE)) 
plot_nmds = plot_nmds + labs(shape = "Sponge class", fill = "Photic zone", colour = "Depth")
plot_nmds = plot_nmds + theme_bw(base_size = 25)
plot_nmds
plot_nmds_classzone_ordi = plot_nmds


ggsave(filename = "Plot_NMDS_class_zone_ordi.pdf", 
       plot = plot_nmds_classzone_ordi, 
       device = "pdf" , 
       width = 30 , height = 25, units = "cm", 
       path = "./3_Beta_div_results")


plot_nmds = ggplot(data_nmds_samples, aes(x = NMDS1, y = NMDS2))
plot_nmds = plot_nmds + geom_point(aes(fill = Genus__celso, pch = Class_host), size = 7, alpha = 0.8)
plot_nmds = plot_nmds + scale_shape_manual(values = c(22,23))
plot_nmds = plot_nmds + scale_fill_manual(values = color_vector_genus)
plot_nmds = plot_nmds + theme(plot.title = element_text(size=22, face="bold"), axis.text=element_text(size=12),axis.title=element_text(size=14))
plot_nmds = plot_nmds + guides(fill = guide_legend(override.aes = list(shape = 22, size = 7)))
plot_nmds = plot_nmds + annotate("text", label="2D Stress = 0.16", x = -0.2, y = 0.15, size = 5)
plot_nmds = plot_nmds + stat_contour(data = ordi.mite.na, aes(x = x, y = y, z = z, colour = after_stat(level)))
plot_nmds =  plot_nmds + scale_colour_gradient2(low = "gold", mid = "aquamarine", high = "darkblue", midpoint = 200, guide=guide_colourbar(reverse = TRUE)) 
plot_nmds = plot_nmds + labs(shape = "Sponge class", fill = "Sponge genus", colour = "Depth")
plot_nmds = plot_nmds + guides(shape = guide_none(), colour = guide_none())
plot_nmds = plot_nmds + theme_bw(base_size = 25)
plot_nmds
plot_nmds_classgenus_ordi = plot_nmds


ggsave(filename = "Plot_NMDS_classgenus_ordi.pdf", 
       plot = plot_nmds_classgenus_ordi, 
       device = "pdf" , 
       width = 30 , height = 25, units = "cm", 
       path = "./3_Beta_div_results")



### 5.3.4  NMDS analysis with species instead of samples ####

data_nmds_species <- plot_ordination(physeq_compo, nmds, type = "taxa", justDF = TRUE)
data_nmds_species


sum_taxa = taxa_sums(physeq_compo)
sum_taxa
select.taxa = plot_ordination(physeq_compo, nmds, type="taxa", justDF = TRUE)
select.taxa
select.taxa.2 = cbind(select.taxa, sum_taxa)
select.taxa.2

select.taxa.4 = subset(select.taxa.2, sum_taxa > 0.1 )
select.taxa.4
datatable(select.taxa.4)

check_asv_names <- select.taxa.4[order(select.taxa.4$sum_taxa,decreasing=TRUE),]
head(check_asv_names)

head(select.taxa.4)
select.taxa.4$Phylum
levels(factor(select.taxa.4$Phylum))


colorpal_phylum = c("Archaea | Crenarchaeota" =    "darkorchid" ,
                    "Archaea | Nanoarchaeota" = "hotpink3",
                    "Bacteria | Acidobacteriota"   ="darkgoldenrod1",
                    "Bacteria | Actinobacteriota"  ="darkslateblue",
                    "Bacteria | Bacteroidota" = "tomato",
                    "Bacteria | Chloroflexi"     =  "olivedrab1",
                    "Bacteria | Dadabacteria"     = "plum1",
                    "Bacteria | Entotheonellaeota" ="orangered1",
                    "Bacteria | Myxococcota"       ="slategray2",
                    "Bacteria | NA" = "grey",
                    "Bacteria | Nitrospinota"     = "darkblue",
                    "Bacteria | Nitrospirota"     ="aquamarine4",
                    "Bacteria | PAUC34f"          ="yellow1",
                    "Bacteria | Proteobacteria"    ="dodgerblue",
                    "Bacteria | SAR324 clade(Marine group B)" = "brown",
                    "Bacteria | Spirochaetota"     ="peachpuff")



plot_nmds = ggplot(data = data_nmds_samples, aes(NMDS1, NMDS2))
plot_nmds = plot_nmds + geom_point(aes(pch = Class_host) , color = "black" , fill = "grey" , size = 7, alpha = 0.3)
plot_nmds = plot_nmds + scale_shape_manual(values = c(22,23))
plot_nmds = plot_nmds + annotate("text", label="2D Stress = 0.16", x = -0.2, y = 0.15, size = 5)
plot_nmds = plot_nmds +  geom_point(data = select.taxa.4, aes(NMDS1, NMDS2, fill = Phylum, size = sum_taxa), stroke = 1 , color = "black", alpha = 0.8, pch = 21) #+ facet_wrap(vars(Family))
plot_nmds = plot_nmds + scale_size_continuous(range = c(1, 15)) 
plot_nmds = plot_nmds  + scale_fill_manual(values = colorpal_phylum)
plot_nmds = plot_nmds + theme_bw(base_size = 25)  + theme(legend.position="right") 
plot_nmds = plot_nmds + guides(fill = guide_legend(override.aes = list(size = 5)))
plot_nmds = plot_nmds + labs(shape = "Sponge class", fill = "Sponge genus", size = "ASV relative abundance (%)", colour = "Depth")
plot_nmds = plot_nmds + guides(shape = guide_none())
plot_nmds_species = plot_nmds
plot_nmds_species


#for ASV id
plot_nmds = plot_nmds + geom_text_repel(data = select.taxa.4, aes(NMDS1, NMDS2, label = row.names(select.taxa.4)), size = 3, box.padding = 0.8)
plot_nmds
ggsave(filename = "Plot_NMDS_species.pdf", 
       plot = plot_nmds_species, 
       device = "pdf" , 
       width = 35 , height = 30, units = "cm", 
       path = "./3_Beta_div_results")


plot_nmds = ggplot(data = data_nmds_samples, aes(NMDS1, NMDS2))
plot_nmds = plot_nmds + geom_point(aes(pch = Class_host) , color = "black" , fill = "grey" , size = 7, alpha = 0.3)
plot_nmds = plot_nmds + scale_shape_manual(values = c(22,23))
plot_nmds = plot_nmds + annotate("text", label="2D Stress = 0.16", x = -0.2, y = 0.15, size = 5)
plot_nmds = plot_nmds +  geom_point(data = select.taxa.4, aes(NMDS1, NMDS2, fill = Phylum, size = sum_taxa), stroke = 1 , color = "black", alpha = 0.8, pch = 21) #+ facet_wrap(vars(Family))
plot_nmds = plot_nmds + scale_size_continuous(range = c(1, 15)) 
plot_nmds = plot_nmds  + scale_fill_manual(values = colorpal_phylum)
plot_nmds = plot_nmds + theme_bw(base_size = 25)  + theme(legend.position="right") 
plot_nmds = plot_nmds + guides(fill = guide_legend(override.aes = list(size = 5)))
plot_nmds = plot_nmds + stat_contour(data = ordi.mite.na, aes(x = x, y = y, z = z, colour = after_stat(level)))
plot_nmds = plot_nmds + scale_colour_gradient2(low = "gold", mid = "aquamarine", high = "darkblue", midpoint = 200, guide=guide_colourbar(reverse = TRUE)) 
plot_nmds = plot_nmds + labs(shape = "Sponge class", fill = "Prokaryotic phylum", size = "ASV relative abundance (%)", colour = "Depth")
plot_nmds = plot_nmds + guides(shape = guide_none(), color = guide_none())
plot_nmds
plot_nmds_species_ordi = plot_nmds


ggsave(filename = "Plot_NMDS_species_ordi.pdf", 
       plot = plot_nmds_species_ordi, 
       device = "pdf" , 
       width = 35 , height = 30, units = "cm", 
       path = "./3_Beta_div_results")




### 5.3.5 Triplot NMDS with ordisurf ####

legend_nmds_species = as_ggplot(get_legend(plot_nmds_species_ordi+ theme_bw(base_size = 12)))
legend_nmds_species

legend_nmds_genus = as_ggplot(get_legend(plot_nmds_classgenus_ordi + theme_bw(base_size = 17)))
legend_nmds_genus

legend_nmds_ordi_zone = as_ggplot(get_legend(plot_nmds_classzone_ordi+ theme_bw(base_size = 17)))
legend_nmds_ordi_zone


plot_nmds_species_ordi_f = plot_nmds_species_ordi + theme(legend.position="none")
plot_nmds_species_ordi_f

plot_nmds_classgenus_ordi_f = plot_nmds_classgenus_ordi + theme(legend.position="none") 
plot_nmds_classgenus_ordi_f

plot_nmds_classzone_ordi_f = plot_nmds_classzone_ordi + theme(legend.position="none") 
plot_nmds_classzone_ordi_f


legend_all = plot_grid(legend_nmds_ordi_zone, legend_nmds_genus ,legend_nmds_species , ncol = 3 )
legend_all

triplot_nmds_ordi = plot_grid( plot_nmds_classzone_ordi_f , plot_nmds_classgenus_ordi_f , plot_nmds_species_ordi_f ,  legend_all , labels=c("A", "B" ,"C", ""),rel_widths = c(1, 1 , 1, 2 )  ,  ncol = 2, nrow = 2 ,label_size = 20)
triplot_nmds_ordi

ggsave(filename = "Triplot_NMDS_ordi.pdf", 
       plot = triplot_nmds_ordi, 
       device = "pdf" , 
       width = 40 , height = 35, units = "cm", 
       path = "./3_Beta_div_results")



## 5.5 Multivariate test to test the effect of each factor on the overall variance of the beta-diversity ####


data_distbeta = as.matrix(distance(physeq_compo, method="bray"))
data_distbeta

metadata_subsampled <- as(sample_data(physeq_compo), "data.frame")
metadata_subsampled


### 5.5.1 Two-way permanova class*zone ####
data_permnova = adonis2(data_distbeta ~ Class_host*Photic_zone, data = metadata_subsampled)
data_permnova
datatable(data_permnova)

write.csv(as.data.frame(data_permnova), 
          file.path("./3_Beta_div_results" , "Data_twoway_permnova_class_zone.csv"))


### 5.5.2 Nested permanova class/genus ####
data_permnova = adonis2(data_distbeta ~ Class_host/Genus__celso, data = metadata_subsampled)
data_permnova

datatable(data_permnova)

write.csv(as.data.frame(data_permnova), 
          file.path("./3_Beta_div_results" , "Data_twoway_permnova_class_genus.csv"))

data_permnova = adonis2(data_distbeta ~ Class_host/Genus__celso/Species_celso, data = metadata_subsampled)
data_permnova

datatable(data_permnova)

write.csv(as.data.frame(data_permnova), 
          file.path("./3_Beta_div_results" , "Data_twoway_permnova_class_genus_species.csv"))


### 5.5.3 Pairwise adonis ####

metadata_subsampled["Class_zone"] <- paste(metadata_subsampled$Class_host, metadata_subsampled$Photic_zone)
datatable(metadata_subsampled)

data_pairwiseadonis_classzone= pairwise.adonis(data_distbeta, metadata_subsampled$Class_zone)
data_pairwiseadonis_classzone
datatable(data_pairwiseadonis_classzone)

write.csv(data_pairwiseadonis_classzone, 
          file.path("./3_Beta_div_results" , "Pairwise_adonis_classzone.csv"))

data_pairwiseadonis_genus = pairwise.adonis(data_distbeta, metadata_subsampled$Genus_host)
data_pairwiseadonis_genus
datatable(data_pairwiseadonis_genus)

write.csv(data_pairwiseadonis_genus, 
          file.path("./3_Beta_div_results" , "Pairwise_adonis_genus.csv"))


# 6. Compositional analysis ####

## 6.1. Barplot analysis at the family level  ####

physeq_aggreg = aggregate_rare(physeq_compo, level = "Family", detection = 5/100, prevalence = 0/100)
physeq_aggreg

write.csv(as.data.frame(otu_table(physeq_aggreg)), "Family_table.csv")

tax_table_family = as.data.frame(tax_table(physeq_aggreg))
tax_table_family$Family

color_families = c("Archaea | Aenigmarchaeota | Aenigmarchaeia | Aenigmarchaeales | NA"                  = "purple3",                   
                     "Archaea | Crenarchaeota | Nitrososphaeria | Nitrosopumilales | Nitrosopumilaceae"  = "lightslateblue"         ,         
                     "Archaea | Nanoarchaeota | Nanoarchaeia | Woesearchaeales | NA"                     = "mediumpurple3"    ,               
                     "Bacteria | Acidobacteriota | Subgroup 11 | NA | NA"                                = "paleturquoise3"           ,        
                     "Bacteria | Acidobacteriota | Vicinamibacteria | Vicinamibacterales | NA"           = "powderblue"       ,           
                     "Bacteria | Actinobacteriota | Acidimicrobiia | Microtrichales | Microtrichaceae"   = "wheat4"   ,               
                     "Bacteria | Bacteroidota | Bacteroidia | Cytophagales | Cyclobacteriaceae"          = "gold"              ,   
                     "Bacteria | Bacteroidota | Bacteroidia | Flavobacteriales | Flavobacteriaceae"      = "orange"           ,       
                     "Bacteria | Bacteroidota | Kapabacteria | Kapabacteriales | NA"                     = "sienna1"      ,             
                     "Bacteria | Chloroflexi | Anaerolineae | Caldilineales | Caldilineaceae"            = "chartreuse2"      ,             
                     "Bacteria | Chloroflexi | Anaerolineae | SBR1031 | A4b"                             = "lightgreen"        ,           
                     "Bacteria | Chloroflexi | Dehalococcoidia | SAR202 clade | NA"                      = "limegreen"        ,           
                     "Bacteria | Chloroflexi | TK10 | NA | NA"                                           = "greenyellow"       ,            
                     "Bacteria | Myxococcota | bacteriap25 | NA | NA"                                    = "peachpuff2"       ,            
                     "Bacteria | NA | NA | NA | NA"                                                      = "grey70"           ,       
                     "Bacteria | Nitrospinota | Nitrospinia | Nitrospinales | Nitrospinaceae"            = "pink3"            ,       
                     "Bacteria | Nitrospirota | Nitrospiria | Nitrospirales | Nitrospiraceae"            = "hotpink"          ,         
                     "Bacteria | PAUC34f | NA | NA | NA"                                                 = "yellow"           ,        
                     "Bacteria | Planctomycetota | Phycisphaerae | Phycisphaerales | Phycisphaeraceae"   = "firebrick"        ,           
                     "Bacteria | Planctomycetota | Planctomycetes | Pirellulales | Pirellulaceae"        = "orangered"        ,           
                     "Bacteria | Proteobacteria | Alphaproteobacteria | NA | NA"                         = "turquoise4"       ,            
                     "Bacteria | Proteobacteria | Alphaproteobacteria | Parvibaculales | PS1 clade"      = "mediumaquamarine"  ,                 
                     "Bacteria | Proteobacteria | Gammaproteobacteria | Burkholderiales | EC94"          = "turquoise1"         ,          
                     "Bacteria | Proteobacteria | Gammaproteobacteria | Burkholderiales | NA"            = "deepskyblue4"        ,           
                     "Bacteria | Proteobacteria | Gammaproteobacteria | Enterobacterales | Alteromonadaceae"                = "turquoise3",
                     "Bacteria | Proteobacteria | Gammaproteobacteria | Enterobacterales | Colwelliaceae"                   = "cornflowerblue",
                     "Bacteria | Proteobacteria | Gammaproteobacteria | Enterobacterales | Pseudoalteromonadaceae"          = "midnightblue",
                     "Bacteria | Proteobacteria | Gammaproteobacteria | Enterobacterales | Vibrionaceae"                    = "deepskyblue",
                     "Bacteria | Proteobacteria | Gammaproteobacteria | Gammaproteobacteria Incertae Sedis | Unknown Family"= "slategrey",
                     "Bacteria | Proteobacteria | Gammaproteobacteria | JTB23 | NA"                                         = "turquoise4",
                     "Bacteria | Proteobacteria | Gammaproteobacteria | NA | NA"                                            = "royalblue4",
                     "Bacteria | Proteobacteria | Gammaproteobacteria | Nitrosococcales | Nitrosococcaceae"                 = "cadetblue1",
                     "Bacteria | Proteobacteria | Gammaproteobacteria | Pseudomonadales | Thioglobaceae"                    = "royalblue2",
                     "Bacteria | Proteobacteria | Gammaproteobacteria | UBA10353 marine group | NA"                         = "dodgerblue2",
                     "Bacteria | SAR324 clade(Marine group B) | NA | NA | NA"                                               = "skyblue",
                     "Bacteria | Spirochaetota | Spirochaetia | Spirochaetales | Spirochaetaceae"                           = "chocolate4",
                     "Other" = "grey")
                     

# for a barplot without border the phyloseq function plot_bar need to be changed. 
# Use the plot_bar_2 function as follows

plot_bar_2 <-  function (physeq, x = "Sample", y = "Abundance", fill = NULL, title = NULL, facet_grid = NULL, border_color = NA) 
{
  mdf = psmelt(physeq)
  p = ggplot(mdf, aes_string(x = x, y = y, fill = fill))
  p = p + geom_bar(stat = "identity", position = "stack",  color = border_color)
  p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0))
  if (!is.null(facet_grid)) {
    p <- p + facet_grid(facet_grid)
  }
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}



barplot = plot_bar_2(physeq_aggreg ,"Sample_names_celso" ,fill = "Family") + theme_bw()
barplot = barplot + geom_bar(stat = "identity"  , position="stack")
barplot = barplot + scale_fill_manual(values = color_families)
barplot = barplot + theme(legend.text = element_text(size=9),
                          axis.ticks.y = element_blank(),
                          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                          legend.position="bottom",
                          legend.title = element_blank(),
                          axis.ticks.x=element_blank())
barplot = barplot + facet_nested(~Class_host+Photic_zone, scales = "free", space = "free")
barplot = barplot + guides(fill = guide_legend(ncol = 2)) 
barplot = barplot +  scale_y_percent()
barplot = barplot + theme(axis.title = element_blank())
barplot


ggsave(filename = "Barplot_family.pdf", 
       plot = barplot, 
       device = "pdf" , 
       width = 40 , height = 30, units = "cm", 
       path = "./4_Compositional_results")

## 6.2. Barplot analysis for Archaea only ####


physeq_subset_archaea = subset_taxa(physeq_subsampled, (Kingdom == "Archaea" ))
physeq_subset_archaea

datatable(tax_table(physeq_subset_archaea))

physeq_archaea_compo = transform(physeq_subset_archaea, "compositional")
physeq_archaea_compo

physeq_archaea_aggreg = aggregate_rare(physeq_archaea_compo, level = "Genus", detection = 1/100, prevalence = 0/100)
physeq_archaea_aggreg

tax_table_genus= as.data.frame(tax_table(physeq_archaea_aggreg))
tax_table_genus$Genus


color_vector_genus_archaea = c( "Archaea | Aenigmarchaeota | Aenigmarchaeia | Aenigmarchaeales | NA | NA" = "deeppink",                                      
                               "Archaea | Aenigmarchaeota | Deep Sea Euryarchaeotic Group(DSEG) | NA | NA | NA"  = "hotpink",                              
                               "Archaea | Crenarchaeota | Nitrososphaeria | Caldiarchaeales | Geothermarchaeaceae | NA" = "lightslateblue",                       
                               "Archaea | Crenarchaeota | Nitrososphaeria | Group 1.1c | NA | NA" = "purple4",                                  
                               "Archaea | Crenarchaeota | Nitrososphaeria | Marine Benthic Group A | NA | NA"  = "mediumpurple3",                             
                               "Archaea | Crenarchaeota | Nitrososphaeria | NA | NA | NA"   = "purple",                              
                               "Archaea | Crenarchaeota | Nitrososphaeria | Nitrosopumilales | Nitrosopumilaceae | Candidatus Nitrosopelagicus" = "dodgerblue2",
                               "Archaea | Crenarchaeota | Nitrososphaeria | Nitrosopumilales | Nitrosopumilaceae | Candidatus Nitrosopumilus"  = "cadetblue1", 
                               "Archaea | Crenarchaeota | Nitrososphaeria | Nitrosopumilales | Nitrosopumilaceae | Cenarchaeum"   = "royalblue3",  
                               "Archaea | Crenarchaeota | Nitrososphaeria | Nitrosopumilales | Nitrosopumilaceae | NA"  = "lightseagreen",     
                               "Archaea | Nanoarchaeota | Nanoarchaeia | Woesearchaeales | NA | NA"     = "orange2",        
                               "Archaea | Nanoarchaeota | Nanoarchaeia | Woesearchaeales | SCGC AAA011-D5 | NA" =  "gold",          
                               "Archaea | Thermoplasmatota | Thermoplasmata | Marine Group II | NA | NA"      = "limegreen",              
                               "Other" = "grey")



barplot = plot_bar_2(physeq_archaea_aggreg ,"Sample_names_celso" ,fill = "Genus") + theme_bw()
barplot = barplot + geom_bar(stat = "identity"  , position="stack")
barplot = barplot + scale_fill_manual(values = color_vector_genus_archaea)
barplot = barplot + theme(legend.text = element_text(size=8),
                          axis.ticks.y = element_blank(),
                          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                          legend.position="bottom",
                          legend.title = element_blank(),
                          axis.ticks.x=element_blank())
barplot = barplot + facet_nested(~Class_host+Photic_zone, scales = "free", space = "free")
barplot = barplot + guides(fill = guide_legend(ncol = 2)) 
barplot = barplot +  scale_y_percent()
barplot = barplot + theme(axis.title = element_blank())
barplot
barplot_archaea = barplot

ggsave(filename = "Barplot_archaea_genus.pdf", 
       plot = barplot_archaea, 
       device = "pdf" , 
       width = 35 , height = 20, units = "cm", 
       path = "./4_Compositional_results")

# 7. Discriminant analyses - Metacoder heattrees ####

### 7.1. Preparation of the datatset ####

ASV_table = read.csv(file = "ASV_table.csv" , header = TRUE , sep = ";" , row.names = 1)
dim(ASV_table)
#datatable(ASV_table)

TAX_table = read.csv(file = "Taxonomy_table.csv" , sep = ";" , header = T , row.names = 1)
dim(TAX_table)
#datatable(TAX_table)

META_table = read.csv(file = "Metadata_table3b.csv" , header = TRUE , sep = ";" , row.names = 1)
dim(META_table)

ASV_phylo = otu_table(ASV_table, taxa_are_rows = TRUE)
dim(ASV_phylo)

TAX_table = as.matrix(TAX_table)

TAX_phylo = tax_table(TAX_table)
dim(TAX_phylo)

META_phylo = sample_data(META_table)
dim(META_phylo)

physeq_raw = phyloseq(ASV_phylo, TAX_phylo, META_phylo)
physeq_raw

sample_data(physeq_raw)$is.neg <- sample_data(physeq_raw)$Sample_or_Control == "Control"

contamdf.prev05 <- isContaminant(physeq_raw, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant)

contamdf.prev05 = cbind(as.data.frame(tax_table(physeq_raw)) , contamdf.prev05)
contamdf.prev05

physeq_decontam = prune_taxa(!contamdf.prev05$contaminant, physeq_raw)
physeq_decontam

physeq_filtered = subset_taxa(physeq_decontam, 
                              (Kingdom != "Eukaryota" | is.na(Kingdom)) &
                                (Order!= "Chloroplast"| is.na(Order)) &
                                (Family != "Mitochondria"| is.na(Family)))
physeq_filtered

physeq_presubsampled = prune_samples(sample_sums(physeq_filtered)>=min_reads_threshold , physeq_filtered)
physeq_presubsampled

physeq_subsampled = subset_samples(physeq_presubsampled, ID_Samples != "CAO_420")
physeq_subsampled

physeq_compo = transform(physeq_subsampled, "compositional")
physeq_compo


## 7.2.1 Metacoder analyses class differences ####

physeq_compo_subset = prune_taxa(taxa_sums(physeq_compo) > 0.04, physeq_compo) 
physeq_compo_subset

metaco_compo <- parse_phyloseq(physeq_compo_subset)
metaco_compo

metadata_metaco_compo <- metaco_compo$data$sample_data
metadata_metaco_compo

metaco_compo$data$rel_abd <- calc_obs_props(metaco_compo, "otu_table", other_cols = T)
metaco_compo$data$rel_abd 

## Calculate per-taxon abundance
metaco_compo$data$tax_rel_abd <- calc_taxon_abund(metaco_compo, "rel_abd")
print(metaco_compo)
metaco_compo$data$tax_rel_abd 

metadata_metaco_compo$sample_id
metadata_metaco_compo$Class_host

#calculate differential data between the two groups to compare
metaco_compo$data$diff_group <- compare_groups(metaco_compo, data = "tax_rel_abd",
                                                  cols = metadata_metaco_compo$sample_id,
                                                  groups = metadata_metaco_compo$Class_host)

metaco_compo$data$diff_group

metaco_compo$data$diff_group$wilcox_p_value

metaco_compo$data$diff_group$wilcox_p_value <- p.adjust(metaco_compo$data$diff_group$wilcox_p_value, method = "fdr")
metaco_compo$data$diff_group$log2_median_ratio[metaco_compo$data$diff_group$wilcox_p_value > 0.05] <- 0
metaco_compo$data$diff_group$log2_median_ratio


Heattree_class =heat_tree(metaco_compo,
                         node_size = n_obs,
                         node_color = log2_median_ratio,
                         node_label = taxon_names,
                         tree_label = taxon_names, 
                         node_size_axis_label = "OTU count",
                         node_color_interval = c(-9, 9),
                         node_color_range = c("chartreuse","gray", "purple3"),
                         node_color_axis_label = "Log 2 ratio of median proportions",
                         layout = "davidson-harel", # The primary layout algorithm
                         initial_layout = "reingold-tilford") # The layout algorithm that initializes node locations

Heattree_class  

ggsave(filename = "Heattree_class.pdf", 
       plot = Heattree_class, 
       device = "pdf" , 
       width = 20 , height = 20, units = "cm", 
       path = "./5_Differential_results")

## 7.2.2 Metacoder analyses zone differences Demospongiae ####

physeq_subsampled_demospongiae = subset_samples(physeq_subsampled,  Class_host == "Demospongiae" )
physeq_subsampled_demospongiae

physeq_compo_demospongiae = transform(physeq_subsampled_demospongiae, "compositional")
physeq_compo_demospongiae
  
physeq_compo_subset = prune_taxa(taxa_sums(physeq_compo_demospongiae) > 0.04, physeq_compo_demospongiae) 
physeq_compo_subset

metaco_compo <- parse_phyloseq(physeq_compo_subset)
metaco_compo

metadata_metaco_compo <- metaco_compo$data$sample_data
metadata_metaco_compo

metaco_compo$data$rel_abd <- calc_obs_props(metaco_compo, "otu_table", other_cols = T)
metaco_compo$data$rel_abd 

## Calculate per-taxon abundance
metaco_compo$data$tax_rel_abd <- calc_taxon_abund(metaco_compo, "rel_abd")
print(metaco_compo)
metaco_compo$data$tax_rel_abd 

metadata_metaco_compo$sample_id
metadata_metaco_compo$Photic_zone

#calculate differential data between the two groups to compare
metaco_compo$data$diff_group <- compare_groups(metaco_compo, data = "tax_rel_abd",
                                               cols = metadata_metaco_compo$sample_id,
                                               groups = metadata_metaco_compo$Photic_zone)

metaco_compo$data$diff_group
datatable(metaco_compo$data$diff_group)
metaco_compo$data$diff_group$wilcox_p_value

#metaco_compo$data$diff_group$wilcox_p_value <- p.adjust(metaco_compo$data$diff_group$wilcox_p_value, method = "fdr")
#metaco_compo$data$diff_group$log2_median_ratio[metaco_compo$data$diff_group$wilcox_p_value > 0.05] <- 0
#metaco_compo$data$diff_group$log2_median_ratio

diverging_palette()
divergent_palette_res = c("#a6611a", "#DDDDDD", "#018571")

color_vector_zone_metacoder_ok12 = c("aquamarine", "#DDDDDD", "yellow1" )
color_vector_zone_metacoder_ok13  = c("royalblue", "#DDDDDD", "yellow1" )
color_vector_zone_metacoder_ok23 = c("aquamarine", "#DDDDDD", "royalblue" )

Heattree_zone_demospongiae = heat_tree_matrix(metaco_compo,
                                               data = "diff_group",
                                               node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
                                               node_label = taxon_names,
                                               node_color = log2_median_ratio, # A column from `obj$data$diff_table`
                                               node_color_range = color_vector_zone_metacoder_ok23, # The built-in palette for diverging data
                                               node_color_trans = "linear", # The default is scaled by circle area
                                               node_color_interval = c(-2, 2), # The range of `log2_median_ratio` to display
                                               edge_color_interval = c(-2, 2), # The range of `log2_median_ratio` to display
                                               node_size_axis_label = "Number of OTUs",
                                               node_color_axis_label = "Log2 ratio median proportions",
                                               layout = "davidson-harel", # The primary layout algorithm
                                               initial_layout = "reingold-tilford", # The layout algorithm that initializes node locations
                                               output_file = "differential_heat_tree.pdf") # Saves the plot as a pdf file

Heattree_zone_demospongiae

ggsave(filename = "Heattree_zone_demospongiae_ok23.pdf", 
       plot = Heattree_zone_demospongiae, 
       device = "pdf" , 
       width = 20 , height = 20, units = "cm", 
       path = "./5_Differential_results")

## 7.2.3 Metacoder analyses zone differences Hexactinellida ####

physeq_subsampled_hexactinellida = subset_samples(physeq_subsampled,  Class_host == "Hexactinellida" )
physeq_subsampled_hexactinellida

physeq_compo_hexactinellida = transform(physeq_subsampled_hexactinellida, "compositional")
physeq_compo_hexactinellida

physeq_compo_subset = prune_taxa(taxa_sums(physeq_compo_hexactinellida) > 0.02, physeq_compo_hexactinellida) 
physeq_compo_subset

metaco_compo <- parse_phyloseq(physeq_compo_subset)
metaco_compo

metadata_metaco_compo <- metaco_compo$data$sample_data
metadata_metaco_compo

metaco_compo$data$rel_abd <- calc_obs_props(metaco_compo, "otu_table", other_cols = T)
metaco_compo$data$rel_abd 

## Calculate per-taxon abundance
metaco_compo$data$tax_rel_abd <- calc_taxon_abund(metaco_compo, "rel_abd")
print(metaco_compo)
metaco_compo$data$tax_rel_abd 

metadata_metaco_compo$sample_id
metadata_metaco_compo$Photic_zone

#calculate differential data between the two groups to compare
metaco_compo$data$diff_group <- compare_groups(metaco_compo, data = "tax_rel_abd",
                                               cols = metadata_metaco_compo$sample_id,
                                               groups = metadata_metaco_compo$Photic_zone)

metaco_compo$data$diff_group
datatable(metaco_compo$data$diff_group)
metaco_compo$data$diff_group$wilcox_p_value

#metaco_compo$data$diff_group$wilcox_p_value <- p.adjust(metaco_compo$data$diff_group$wilcox_p_value, method = "fdr")
#metaco_compo$data$diff_group$log2_median_ratio[metaco_compo$data$diff_group$wilcox_p_value > 0.05] <- 0
#metaco_compo$data$diff_group$log2_median_ratio

diverging_palette()
divergent_palette_res = c("#a6611a", "#DDDDDD", "#018571")


color_vector_zone_metacoder_hexa = c("aquamarine", "#DDDDDD", "royalblue" )



Heattree_hexactinellida =heat_tree(metaco_compo,
                          node_size = n_obs,
                          node_color = log2_median_ratio,
                          node_label = taxon_names,
                          tree_label = taxon_names, 
                          node_size_axis_label = "OTU count",
                          node_color_interval = c(-9, 9),
                          node_color_range = color_vector_zone_metacoder_hexa,
                          node_color_axis_label = "Log 2 ratio of median proportions",
                          layout = "davidson-harel", # The primary layout algorithm
                          initial_layout = "reingold-tilford") # The layout algorithm that initializes node locations

Heattree_hexactinellida  

ggsave(filename = "Heattree_zone_hexactinellida.pdf", 
       plot = Heattree_hexactinellida, 
       device = "pdf" , 
       width = 20 , height = 20, units = "cm", 
       path = "./5_Differential_results")
