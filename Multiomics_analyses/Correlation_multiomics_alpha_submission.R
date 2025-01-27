#### Correlation between alpha-diversity metrics between omcis datasets ####

library(ggplot2)
library(reshape2)
library(ggpubr)
devtools::install_github("thomasp85/patchwork")
library(patchwork)

setwd("C:/Users/bpaix/Mon Drive/2021_2024 Postdoc Naturalis/2022 Curacao subdive project/4. Lab notebook/8. Multiomics analyses/alpha_multiomics")


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



Alpha_table_all <- read.csv(file ="Alpha_table_multiomics.csv", header = TRUE, sep = ";", row.names = 1)
Alpha_table_all
Alpha_table_all$shannon_ASV

ggscatter(Alpha_table_all, x = "Values_Shannon_LCMS", y = "shannon_ASV", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Values_Shannon_LCMS", ylab = "shannon_ASV")

plot_alpha_multiomic = ggplot(Alpha_table_all, aes(Values_Shannon_LCMS ,shannon_ASV))
plot_alpha_multiomic = plot_alpha_multiomic + geom_point(aes(fill = Photic_zone, pch = Class_host), size = 4, alpha = 0.7)
plot_alpha_multiomic = plot_alpha_multiomic + scale_shape_manual(values = c(22,23))
#plot_alpha_multiomic = plot_alpha_multiomic + geom_text(aes(label = Sample_name))
plot_alpha_multiomic = plot_alpha_multiomic + scale_fill_manual(values = color_vector_zone, labels = c("Mesophotic", "Upper rariphotic", "Lower rariphotic"))
plot_alpha_multiomic = plot_alpha_multiomic + theme(plot.title = element_text(size=22, face="bold"), axis.text=element_text(size=12),axis.title=element_text(size=14))
plot_alpha_multiomic = plot_alpha_multiomic + guides(fill = guide_legend(override.aes = list(shape = 22, size = 7)))
plot_alpha_multiomic = plot_alpha_multiomic + labs(shape = "Sponge class", fill = "Photic zone")
plot_alpha_multiomic = plot_alpha_multiomic + theme_bw(base_size = 15) + theme(legend.position="bottom", legend.box="vertical", legend.margin=margin())
plot_alpha_multiomic = plot_alpha_multiomic + geom_smooth(method = "lm", level = 0.8, alpha = 0.2, color = "gray")
plot_alpha_multiomic = plot_alpha_multiomic + labs(x = "Shannon metabolomics dataset", y = "Shannon metabarcoding dataset")
plot_alpha_multiomic = plot_alpha_multiomic + theme(legend.position="none")
plot_alpha_multiomic

density_x = ggplot(Alpha_table_all, aes(x = Values_Shannon_LCMS, fill = Class_host))  
density_x = density_x + geom_density(alpha = 0.4) + theme_void()  
density_x = density_x + scale_fill_manual(values = color_vector_class)
density_x = density_x + theme(legend.position = "none")
density_x

density_y = ggplot(Alpha_table_all, aes(x = shannon_ASV, fill = Class_host)) 
density_y = density_y + geom_density(alpha = 0.4) 
density_y = density_y + theme_void() 
density_y = density_y + theme(legend.position = "none")  
density_y = density_y + scale_fill_manual(values = color_vector_class)
density_y = density_y + coord_flip()
density_y

plot_alpha_multiomic_all = density_x + plot_spacer() + plot_alpha_multiomic + density_y 
plot_alpha_multiomic_all = plot_alpha_multiomic_all +  plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))
plot_alpha_multiomic_all


ggsave(filename = "Plot_alpha_multiomics_all.pdf", 
       plot = plot_alpha_multiomic_all, 
       device = "pdf" , 
       width = 22 , height = 25, units = "cm")


