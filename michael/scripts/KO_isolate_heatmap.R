# Load packages:
library(MetaboSignal)
library(RColorBrewer)
library(plyr)
library(reshape2)
library(ComplexHeatmap)
library(tidyverse)
library(scales)
library(igraph)

####   Link KOs to isolate (abundance)    ####
# Filter abundance table to only include KOs of interest
KO_iso_abund.1 <- as.data.frame(apply(KO2_iso_gn_ct.cast_names[,3:5533], 2, as.numeric))
colnames(KO_iso_abund.1)[] <- colnames(KO2_iso_gn_ct.cast_names)[3:5533]
rownames(KO_iso_abund.1)[] <- gsub(".eggnog.emapper.annotations.csv", "", KO2_iso_gn_ct.cast_names$Isolate)
KO_iso_abund.2 <- KO_iso_abund.1[,which(names(KO_iso_abund.1) %in% paste0('ko.', SPEgmNet2.df$target))]
# change Isolate names to new
KO_iso_abund.2.newnames <- join(IsolateNames_oldandnew, data.frame('old_name' = rownames(KO_iso_abund.2), KO_iso_abund.2))
rownames(KO_iso_abund.2.newnames) <- KO_iso_abund.2.newnames$new_name
KO_iso_abund.2.newnames <- KO_iso_abund.2.newnames[,-c(1:2)]
# Which are the most abundant KOs?
KO_iso_MostAbund.2 = data.frame('SumAbund' = colSums(KO_iso_abund.2))
KO_iso_MostAbund.2.desc <- KO_iso_MostAbund.2[order(KO_iso_MostAbund.2$SumAbund, decreasing = TRUE), , drop = FALSE] # subset of KOs which link to ALI cult metabolites
# Using a manual split:
KO_iso_MostAbund.2.desc.1 <- data.frame("KO" = rownames(KO_iso_MostAbund.2.desc), "Abundance" = KO_iso_MostAbund.2.desc)
KO_iso_MostAbund.2_0.75 <- KO_iso_MostAbund.2.desc.1[which(KO_iso_MostAbund.2.desc.1$SumAbund > 126*0.75),]
KO_iso_MostAbund.2_0.5 <- KO_iso_MostAbund.2.desc.1[which(KO_iso_MostAbund.2.desc.1$SumAbund < 126*0.75 & KO_iso_MostAbund.2.desc.1$SumAbund > 126*0.25),]
KO_iso_MostAbund.2_0.25 <- KO_iso_MostAbund.2.desc.1[which(KO_iso_MostAbund.2.desc.1$SumAbund < 126*0.25),]
KO_iso_MostAbund.2_KOgroups <- rbind(data.frame(KO_iso_MostAbund.2_0.75, "KOgroup" = rep("Frequent", nrow(KO_iso_MostAbund.2_0.75))),
                                     data.frame(KO_iso_MostAbund.2_0.5, "KOgroup" = rep("Intermediate", nrow(KO_iso_MostAbund.2_0.5))),
                                     data.frame(KO_iso_MostAbund.2_0.25, "KOgroup" = rep("Rare", nrow(KO_iso_MostAbund.2_0.25))))
KO_iso_MostAbund.2_KOgroups.1 <- KO_iso_MostAbund.2_KOgroups[order(KO_iso_MostAbund.2_KOgroups$KO),]

# Visualise abundance tables:
col_fun2 <- colorRamp2(c(0,1), c('white', 'red'))
col_fun2(seq(0,2))
row_ha <- rowAnnotation(Isolate_Cluster = as.character(IsolateClusterLabs2$cluster_longname), 
                        col = list(Isolate_Cluster = c("Strep.I" = alpha('#83E952', 0.7), "Strep.II" = alpha('#73E4D7', 0.7), 
                                                       "Veillonella" = alpha('#E6705D', 0.7), "Rothia" = alpha('#7990D7', 0.7), 
                                                       "Gemella" = alpha('#9843EB', 0.7), "Prevotella" = alpha('#7B59C2', 0.7),
                                                       "Micrococcus" = alpha('#89E598', 0.7), "Neisseria" = alpha('#88C0DC', 0.7),
                                                       "Pauljensenia" = alpha('#CCCC65', 0.7), "Staphy.+Nialla" = alpha('#E556D5', 0.7), 
                                                       "Haemophilus" = alpha('#D2E6D6', 0.7), "Granulicatella" = alpha('#79986B', 0.7), 
                                                       "Fusobacterium+Leptotrichia"= alpha('#E2CA95', 0.7),  "Cutibacterium" = alpha('#D897E2', 0.7), 
                                                       "Actinomyces" = alpha('#DBB4C8', 0.7), "Strep.III"= alpha('#E1A249', 0.7))), border = TRUE, show_legend = TRUE,
                        annotation_legend_param = list(Isolate_Cluster = list(ncol= 1,
                                                                        title = 'Isolate Cluster')))

col_ha <- HeatmapAnnotation(Frequency = KO_iso_MostAbund.2_KOgroups.1$KOgroup, 
                            col = list(Frequency = c("Frequent" = '#3288BD', "Intermediate" = '#88419D', "Rare" = '#D53E4F')), border = TRUE,
                            show_legend = TRUE, annotation_legend_param = list(Frequency = list(ncol= 1,
                                                                                           title = 'Frequency', title_position = 'topleft',
                                                                                           legend_height = unit(3, 'cm'), at = c("Frequent", "Intermediate", "Rare"))))

KO_iso_abund_KOnames <- MS_changeNames(gsub('ko.', '', colnames(KO_iso_abund.2.newnames)), 'hsa')
KO_iso_abund_isonames <- stringr::str_remove_all(rownames(KO_iso_abund.2.newnames), "\\d+|_")

Heatmap(KO_iso_abund.2.newnames, name = "KO abundance", width = unit(20, "cm"), height = unit(26, "cm"), col = col_fun2,
         column_names_gp = gpar(fontsize = 8), column_names_rot = 90, show_column_dend = FALSE,
         row_names_gp = gpar(fontsize = 8), cluster_rows = FALSE, row_names_side = 'left', show_row_dend = TRUE, row_dend_side = 'right',
         cluster_columns = TRUE, show_row_names = TRUE, show_column_names = TRUE, rect_gp = gpar(col = "black", lwd = 1),
         heatmap_legend_param = list(title = 'Abundance', at = at,
                                     legend_height = unit(3, 'cm'), grid_width = unit(1, 'cm'),
                                     title_position = 'topleft', border = 'black', labels_gp = gpar(fontsize = 10, fill = at),
                                     title_gp = gpar(fontsize = 10, font = 2)),
         row_labels = KO_iso_abund_isonames, column_labels = KO_iso_abund_KOnames,
         left_annotation = row_ha, top_annotation = col_ha)
