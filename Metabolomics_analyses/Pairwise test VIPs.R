# VIP univarites stats ######


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


setwd("your path here")



#### STATS VIP demospongiae / zone #####

Data_table_normalized = read.csv(file = "./data_normalized_demospongiae.csv" , header = TRUE , sep = ";" , row.names = 1)

datatable(Data_table_normalized)
dim(Data_table_normalized)

VIP_table_normalized = Data_table_normalized[, c("Sample_names","Label", 
                                                 "X2127_510.3875_1069.92", #Lyso-PC(O-C18:0)
                                                 "X1717_176.12565_89.24", #Rhodosamine
                                                 "X1638_155.07952_69.96", #C7H11N2O2  Cyclo-(Pro-Gly)?
                                                 "X406_174.12135_91.21", #C7H15N3O2 Amidinonorleucine ?
                                                 "X1186_522.38681_1217.11", #PC(O-C18:0/O-1:0)
                                                 "X1694_166.08386_71.91", #C9H12NO2 Phenylalanine ?
                                                 "X423_367.14706_131.15" , #C19H18N4O4 alkaloid
                                                 "X1724_180.09894_89.64",  #Salsolinol
                                                 "X1645_146.11527_69.93", #C7H16NO2 Methyl-isoleucine
                                                 "X1750_137.04436_94.87", #C5H5N4O Hypoxanthine
                                                 "X454_144.07941_287.22", #C5H10N3O2 Thymine ?
                                                 "X402_204.13031_89.27", #Gly-Lys
                                                 "X394_226.11945_71.84",  #Gln-Pro
                                                 "X1626_187.10489_69.52", #Ala-Pro
                                                 "X946_614.36268_994.69")] #C37H48N3O5_NI

datatable(VIP_table_normalized)
dim(VIP_table_normalized)

colnames(VIP_table_normalized) <- c("Sample_names","Photic_zone", 
                                    "Lyso_PC_O_C18_0", 
                                    "Rhodosamine", 
                                    "C7H11N2O2_Cyclo_Pro_Gly",
                                    "C7H15N3O2_Amidinonorleucine",
                                    "PC_O_C18_0_O_1_0", 
                                    "C9H12NO2_Phenylalanine",
                                    "C19H18N4O4_alkaloid",
                                    "Salsolinol",
                                    "C7H16NO2_Methyl_isoleucine",
                                    "C5H5N4O_Hypoxanthine",
                                    "C5H10N3O2", 
                                    "Gly_Lys",
                                    "Gln_Pro",
                                    "Ala_Pro",
                                    "C37H48N3O5_ni")

aov_Shannon_chemodiv_lcms = aov(Lyso_PC_O_C18_0 ~ Photic_zone, VIP_table_normalized)
aov_Shannon_chemodiv_lcms 
summary(aov_Shannon_chemodiv_lcms)

hsd_vip_demos_zone = HSD.test(aov(Lyso_PC_O_C18_0 ~ Photic_zone, VIP_table_normalized), "Photic_zone", group=T)
hsd_vip_demos_zone


aov_Shannon_chemodiv_lcms = aov(Rhodosamine ~ Photic_zone, VIP_table_normalized)
aov_Shannon_chemodiv_lcms 
summary(aov_Shannon_chemodiv_lcms)

hsd_vip_demos_zone = HSD.test(aov(Rhodosamine ~ Photic_zone, VIP_table_normalized), "Photic_zone", group=T)
hsd_vip_demos_zone

aov_Shannon_chemodiv_lcms = aov(C7H11N2O2_Cyclo_Pro_Gly ~ Photic_zone, VIP_table_normalized)
aov_Shannon_chemodiv_lcms 
summary(aov_Shannon_chemodiv_lcms)

hsd_vip_demos_zone = HSD.test(aov(C7H11N2O2_Cyclo_Pro_Gly ~ Photic_zone, VIP_table_normalized), "Photic_zone", group=T)
hsd_vip_demos_zone

aov_Shannon_chemodiv_lcms = aov(C7H15N3O2_Amidinonorleucine ~ Photic_zone, VIP_table_normalized)
aov_Shannon_chemodiv_lcms 
summary(aov_Shannon_chemodiv_lcms)

hsd_vip_demos_zone = HSD.test(aov(C7H15N3O2_Amidinonorleucine ~ Photic_zone, VIP_table_normalized), "Photic_zone", group=T)
hsd_vip_demos_zone

aov_Shannon_chemodiv_lcms = aov(PC_O_C18_0_O_1_0 ~ Photic_zone, VIP_table_normalized)
aov_Shannon_chemodiv_lcms 
summary(aov_Shannon_chemodiv_lcms)

hsd_vip_demos_zone = HSD.test(aov(PC_O_C18_0_O_1_0 ~ Photic_zone, VIP_table_normalized), "Photic_zone", group=T)
hsd_vip_demos_zone

aov_Shannon_chemodiv_lcms = aov(C9H12NO2_Phenylalanine ~ Photic_zone, VIP_table_normalized)
aov_Shannon_chemodiv_lcms 
summary(aov_Shannon_chemodiv_lcms)

aov_Shannon_chemodiv_lcms = aov(C19H18N4O4_alkaloid ~ Photic_zone, VIP_table_normalized)
aov_Shannon_chemodiv_lcms 
summary(aov_Shannon_chemodiv_lcms)

aov_Shannon_chemodiv_lcms = aov(Salsolinol ~ Photic_zone, VIP_table_normalized)
aov_Shannon_chemodiv_lcms 
summary(aov_Shannon_chemodiv_lcms)

hsd_vip_demos_zone = HSD.test(aov(Salsolinol ~ Photic_zone, VIP_table_normalized), "Photic_zone", group=T)
hsd_vip_demos_zone

aov_Shannon_chemodiv_lcms = aov(C7H16NO2_Methyl_isoleucine ~ Photic_zone, VIP_table_normalized)
aov_Shannon_chemodiv_lcms 
summary(aov_Shannon_chemodiv_lcms)

aov_Shannon_chemodiv_lcms = aov(C5H5N4O_Hypoxanthine ~ Photic_zone, VIP_table_normalized)
aov_Shannon_chemodiv_lcms 
summary(aov_Shannon_chemodiv_lcms)

aov_Shannon_chemodiv_lcms = aov(C5H10N3O2 ~ Photic_zone, VIP_table_normalized)
aov_Shannon_chemodiv_lcms 
summary(aov_Shannon_chemodiv_lcms)

aov_Shannon_chemodiv_lcms = aov(Gly_Lys ~ Photic_zone, VIP_table_normalized)
aov_Shannon_chemodiv_lcms 
summary(aov_Shannon_chemodiv_lcms)

hsd_vip_demos_zone = HSD.test(aov(Gly_Lys ~ Photic_zone, VIP_table_normalized), "Photic_zone", group=T)
hsd_vip_demos_zone

aov_Shannon_chemodiv_lcms = aov(Gln_Pro ~ Photic_zone, VIP_table_normalized)
aov_Shannon_chemodiv_lcms 
summary(aov_Shannon_chemodiv_lcms)

hsd_vip_demos_zone = HSD.test(aov(Gln_Pro ~ Photic_zone, VIP_table_normalized), "Photic_zone", group=T)
hsd_vip_demos_zone

Gln_Proaov_Shannon_chemodiv_lcms = aov(Ala_Pro ~ Photic_zone, VIP_table_normalized)
aov_Shannon_chemodiv_lcms 
summary(aov_Shannon_chemodiv_lcms)

aov_Shannon_chemodiv_lcms = aov(C37H48N3O5_ni ~ Photic_zone, VIP_table_normalized)
aov_Shannon_chemodiv_lcms 
summary(aov_Shannon_chemodiv_lcms)



library(dplyr)


mean_by_zone <- VIP_table_normalized %>% group_by(Photic_zone) %>% summarise_at(vars(-Sample_names), funs(mean(., na.rm=TRUE)))
mean_by_zone = as.matrix.data.frame(mean_by_zone)
mean_by_zone = as.data.frame(mean_by_zone)
mean_by_zone

write.csv(mean_by_zone, "Mean_table_VIP_demospongiae.csv")




#### STATS VIP Hexactinellida / zone #####



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
                                                 "X1033_794.59263_1085.76", #PC(O-C38:5) ?
                                                 "X1849_188.06838_303.97", #Imidazole derivative
                                                 "X2060_530.31959_948.86", #Lyso-PE(C22:4)
                                                 "X1615_235.16304_67.96", #C10H22N2O4
                                                 "X402_204.13031_89.27", #Glycyllysine
                                                 "X2031_570.35324_908.5", #Lyso-PC(C22:5)
                                                 "X1679_196.07246_71.6", #Cysteinyl-Glycine ??
                                                 "X2163_271.18483_1112.69", #Unidentified C17H22N2O ?
                                                 "X2339_685.43179_1257.55", #PG(31:4) ?
                                                 "X2007_502.28847_877.11", #Lyso-PE(C20:4)
                                                 "X2222_810.59565_1160.99", #PC(C38:4)
                                                 "X1725_196.07238_88.74")] #Cysteinyl-Glycine ??


datatable(VIP_table_normalized)
dim(VIP_table_normalized)

colnames(VIP_table_normalized) <- c("Sample_names","Photic_zone", 
                                    "LysoPS_C15_1",
                                    "Xanthurenic_acid",
                                    "Lyso_PC_C22:5",
                                    "Leu_Val",
                                    "Leu_Leu",
                                    "Lyso_PE_C22:2",
                                    "PC_O_C38:5_",
                                    "Imidazole_derivative",
                                    "Lyso_PE_C22_4",
                                    "C10H22N2O4",
                                    "Gly_Lys",
                                    "Lyso_PC_C22:5_2",
                                    "Cysteinyl_Glycine",
                                    "C17H22N2O",
                                    "PG_31:4",
                                    "Lyso_PE_C20:4",
                                    "PC_C38:4", 
                                    "Cysteinyl_Glycine_2")

datatable(VIP_table_normalized)
VIP_table_normalized

library(readr)
library(dplyr)

#Initialiser une liste pour stocker les résultats des tests t
resultats_ttest <- list()

# Boucle sur chaque métabolite (toutes les colonnes sauf les deux premières)
for (metabolite in colnames(VIP_table_normalized)[-c(1:2)]) {  
  # Sous-ensembles des données pour les deux groupes
  groupe1 <- VIP_table_normalized[-c(1)] %>%
    dplyr::filter(Photic_zone == "2.Upper rariphotic") %>%
    pull(!!sym(metabolite))
  
  groupe2 <- VIP_table_normalized[-c(1)] %>%
    dplyr::filter(Photic_zone == "3.Lower rariphotic") %>%
    pull(!!sym(metabolite))
  
  # Effectuer le test t
  ttest <- t.test(groupe1, groupe2, var.equal = TRUE)  # ou var.equal = FALSE
  
  # Sauvegarder le résultat
  resultats_ttest[[metabolite]] <- ttest
  
  # Afficher un résumé
  cat("\nTest t pour", metabolite, ":\n")
  print(ttest)
}


mean_by_zone <- VIP_table_normalized %>% group_by(Photic_zone) %>% summarise_at(vars(-Sample_names), funs(mean(., na.rm=TRUE)))
mean_by_zone = as.matrix.data.frame(mean_by_zone)
mean_by_zone = as.data.frame(mean_by_zone)
mean_by_zone

write.csv(mean_by_zone, "Mean_table_VIP_hexactinellida.csv")



### STATS VIP Class ####


Data_table_normalized = read.csv(file = "./data_normalized_class.csv" , header = TRUE , sep = ";" , row.names = 1)

datatable(Data_table_normalized)
dim(Data_table_normalized)



VIP_table_normalized = Data_table_normalized[, c("Sample_names","Label", 
                                                 "X2064_510.35064_957.99"  ,
                                                 "X1838_261.14089_274.29"  ,
                                                 "X2022_544.33574_897.56"  ,
                                                 "X2126_524.36706_1070.53"  ,
                                                 "X421_264.15603_127.02"  ,
                                                 "X1625_190.04718_68.03"  ,
                                                 "X1619_168.06684_69.51"  ,
                                                 "X2106_496.37124_1030.26"  ,
                                                 "X1991_542.31997_846.88"  ,
                                                 "X1752_241.15223_94.47"  ,
                                                 "X1890_541.03354_523.3"  ,
                                                 "X2030_508.3351_903.47"  ,
                                                 "X1730_231.08297_90.25"  ,
                                                 "X2127_510.3875_1069.92"  ,
                                                 "X1117_538.38277_1168.36"  ,
                                                 "X1624_198.08495_69.22"  ,
                                                 "X2192_597.33806_1141.14"  ,
                                                 "X2148_550.38314_1095.4"  ,
                                                 "X1626_187.10489_69.52"  ,
                                                 "X2054_496.33453_941.01"  ,
                                                 "X677_271.02398_521.86"  ,
                                                 "X1008_583.32127_1075.61"  ,
                                                 "X2007_502.28847_877.11"  ,
                                                 "X1884_569.04746_504.2"  ,
                                                 "X1657_189.11944_70.58"  ,
                                                 "X1839_206.04277_275.93"  ,
                                                 "X465_246.16695_320.81"  )]
ncol(VIP_table_normalized)

colnames(VIP_table_normalized) <- c("Sample_names","Sponge_class", 
                                    "Lyso_PC_C17_0",
                                    "Dipeptide_or_derivative_Glutamylisoleucine",
                                    "Lyso_PC_C20_4",
                                    "Lyso_PC_C18_0",
                                    "Dipetide_or_derivative_Glutamylvaline_Low_intensity",
                                    "Hydroyphenylglycine",
                                    "Hydroyphenylglycine_2",
                                    "Lyso_PC_O_C17_0_NO6P",
                                    "Lyso_PC_C20_5",
                                    "Unidentified_C12H21N2O3",
                                    "Aphrocallistin",
                                    "Lyso_PC_C17_1",
                                    "Methionylalanine_derivative",
                                    "Lyso_PC_O_C18_0_NO6P",
                                    "Lyso_PC_19_0",
                                    "Methylated_aminoacid_derivative",
                                    "Lyso_PI_19_0_phosphoinositol",
                                    "Lyso_PC_C20_1",
                                    "Alanylproline",
                                    "Lyso_PC_C16_0",
                                    "Aphro_frag",
                                    "Lyso_PI_18_0_phosphoinositol",
                                    "Lyso_PE_20_4",
                                    "Aphrocallistin_derivative__methoxy",
                                    "NI",
                                    "Xanthurenic_acid",
                                    "Valerylcarnitine") 
VIP_table_normalized
#Initialiser une liste pour stocker les résultats des tests t
resultats_ttest <- list()

# Boucle sur chaque métabolite (toutes les colonnes sauf les deux premières)
for (metabolite in colnames(VIP_table_normalized)[-c(1:2)]) {  
  # Sous-ensembles des données pour les deux groupes
  groupe1 <- VIP_table_normalized[-c(1)] %>%
    dplyr::filter(Sponge_class == "Demospongiae") %>%
    pull(!!sym(metabolite))
  
  groupe2 <- VIP_table_normalized[-c(1)] %>%
    dplyr::filter(Sponge_class == "Hexactinellida") %>%
    pull(!!sym(metabolite))
  
  # Effectuer le test t
  ttest <- t.test(groupe1, groupe2, var.equal = TRUE)  # ou var.equal = FALSE
  
  # Sauvegarder le résultat
  resultats_ttest[[metabolite]] <- ttest
  
  # Afficher un résumé
  cat("\nTest t pour", metabolite, ":\n")
  print(ttest)
}



mean_by_class<- VIP_table_normalized %>% group_by(Sponge_class) %>% summarise_at(vars(-Sample_names), funs(mean(., na.rm=TRUE)))
mean_by_class = as.matrix.data.frame(mean_by_class)
mean_by_class = as.data.frame(mean_by_class)
mean_by_class

write.csv(mean_by_class, "Mean_table_VIP_class.csv")


