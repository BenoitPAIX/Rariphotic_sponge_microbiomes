### ALPHA-CHEMODIVERSITY ANALYSIS ########


setwd("your path here")

library(vegan)
library(ggplot2)
library(DT)

Metabo_table <- read.csv(file = "data_norm_chemodiv_filtered.csv", header = TRUE, sep = ";" , row.names = 1)
Metabo_table
dim(Metabo_table)
datatable(Metabo_table)

Metabo_table
#Metabo_table[Metabo_table < 0.5] <- 0
Metabo_table

Metadata <- read.csv(file ="Metadata_table.csv", header = TRUE, sep = ";", row.names = 1)
Metadata
dim(Metadata)
datatable(Metadata)



Shannon_chemodiv <- diversity(t(Metabo_table), "shannon")
Shannon_chemodiv
Shannon_chemodiv = as.data.frame(Shannon_chemodiv)
Shannon_chemodiv

Shannon_chemodiv_meta = cbind(Shannon_chemodiv, Metadata)
Shannon_chemodiv_meta

dim(Shannon_chemodiv_meta)
datatable(Shannon_chemodiv_meta)
Shannon_chemodiv_meta$Genus_host

#write.csv(Shannon_chemodiv_meta, "Alpha_table_LCMS_filtered.csv")


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



Shannon_chemodiv_meta$Shannon_chemodiv
Shannon_chemodiv_meta$shannon



plot_alpha = ggplot(Shannon_chemodiv_meta, aes(Photic_zone ,shannon  , fill = Photic_zone, pch = Class_host)) + guides(fill=guide_legend(title="Photic zone"))
plot_alpha = plot_alpha + geom_boxplot(alpha = 0.8, size = 1) + facet_grid( . ~ Class_host , scales="free", space = "free_x") 
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
plot_alpha_shannon_metabo = plot_alpha



ggsave(filename = "Plot_alpha_class_photic.pdf", 
       plot = plot_alpha, 
       device = "pdf" , 
       width = 20 , height = 15, units = "cm")



shapiro_data_shannon = shapiro.test(Shannon_chemodiv_meta$shannon)
shapiro_data_shannon

# p > 0.05 so normal distrib. so ANOVA test

aov_Shannon_chemodiv_lcms = aov(shannon ~ Class_host*Photic_zone, Shannon_chemodiv_meta)
aov_Shannon_chemodiv_lcms 
summary(aov_Shannon_chemodiv_lcms)


Shannon_chemodiv_meta["Class_zone"] <- paste(Shannon_chemodiv_meta$Class_host, Shannon_chemodiv_meta$Photic_zone)

Shannon_chemodiv_meta$Photic_zone


data_hsd_shannon_chemodiv = HSD.test(aov(shannon ~ Class_zone, Shannon_chemodiv_meta), "Class_zone", group=T)
data_hsd_shannon_chemodiv

data_hsd_shannon_chemodiv = HSD.test(aov(shannon ~ Class_host, Shannon_chemodiv_meta), "Class_host", group=T)
data_hsd_shannon_chemodiv

data_hsd_shannon_chemodiv = HSD.test(aov(shannon ~ Photic_zone, Shannon_chemodiv_meta), "Photic_zone", group=T)
data_hsd_shannon_chemodiv



aov_Shannon_chemodiv_lcms = aov(shannon ~ Class_host/Genus_host, Shannon_chemodiv_meta)
aov_Shannon_chemodiv_lcms 
summary(aov_Shannon_chemodiv_lcms)