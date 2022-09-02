###
library(GGally)
library(ggplot2)
library(rtk)
library(edgeR)
library(vegan)
library(tidyverse)
library(dplyr)
library(phyloseq)
library(metadeconfoundR)

load(file = "~/Dropbox/lung/ThedaTill/raw_data/metadata_red.r")
load(file= "~/Dropbox/lung/ThedaTill/intermediate_files/RTK_div.r")
load(file = "~/Dropbox/lung/ThedaTill/intermediate_files/Genus_nice_meta_name.r")

colnames(RTK_div)[which(names(RTK_div) == "invsimpson")] <- "Invsimpson index"
colnames(RTK_div)[which(names(RTK_div) == "simpson")] <- "Simpson index"
colnames(RTK_div)[which(names(RTK_div) == "eveness")] <- "Evenness"
colnames(RTK_div)[which(names(RTK_div) == "shannon.x")] <- "Shannon index"
colnames(RTK_div)[which(names(RTK_div) == "richness")] <- "Richness"

## study ID and site as randome effects 
## alwasy check with order(row.names((metadata_red)))

rownames(RTK_div) <- sub("^X", "", rownames(RTK_div))
metadata_alpha <- metadata_red
metadata_alpha[,c("simpson_evenness", "Shannon_diversity","Simpson_diversity")] <-NULL
# metadata_alpha <- metadata_red[,-c(13:15)] ACHTUNG dont use it more than ones otherwise everything is messed 
metadeconf_Celtic_RTK_div <- MetaDeconfound(featureMat = RTK_div,
                                              metaMat = metadata_alpha,
                                              nnodes = 3,
                                              randomVar = list("+ (1|Study_ID) + (1|Site)", c("Study_ID", "Site")), 
                                              logfile = "~/Dropbox/lung/ThedaTill/intermediate_files/20210527_metadeconf_div.log")


plot <- BuildHeatmap(metadeconf_Celtic_RTK_div, metaVariableNames = nice.meta, d_range = "full")

heatmap <-
 plot + theme(legend.position = "none",
                 axis.text.y = element_text(size = 14),
                 axis.text.x = element_text(size = 14))
heatmap


ggsave(filename = "~/Dropbox/lung/ThedaTill/output/Updated_Figure/20210527_metadeconf_alpha_div.svg", 
       plot = plot, device = "svg", width = 4, height = 3)




