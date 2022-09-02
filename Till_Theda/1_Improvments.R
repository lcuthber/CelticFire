## Till and Theda 
## Start 11.05.2021
library(asbio)
library(picante)
library(DirichletMultinomial)
library(reshape2)
library(magrittr)
library(dplyr)
library(cowplot)
library(ggrepel)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(phyloseq)
library(devtools)
library(ggpubr)
library(rstatix)
library(patchwork)
library(ggrepel)

### Make new file for Perl script run to detect make congruence 
metadata <- read.table("~/Dropbox/lung/ThedaTill/raw_data/cf_formatted_mapfile.csv",
                       sep = ",",
                       header = TRUE, stringsAsFactors = F)
# there is an issues with 0 in the sampleID names
## separate metadata by samplesID with and w/O r in the name

metadata <- metadata[- (grep(pattern = "*r", x = metadata$SampleID)), ]
metadata$SampleID <- sub(pattern = "^0+", replace ="", x = metadata$SampleID) # remove all zeros
metadata <- rbind(metadata, metadata_endr) # rbind the two sub metadatafiles 
row.names(metadata) <- metadata$SampleID
metadata$SampleID <- NULL

metadata <- metadata[order(rownames(metadata)), ]
# check for values 
summary(metadata)
# clean up the metadata remove genus, change sex == 2 to 0 
metadata$Sex[metadata$Sex == 2 ] <- 0
metadata <- metadata[, -c(18:42)]
metadata$Ethnicity <- paste0("Ethnicity_", metadata$Ethnicity)

### create file for the perl script 
meta_make_congruence <- metadata
Pulmotypes <- read.table("~/Dropbox/lung/ThedaTill/20201117_Enterotyping/20201217_pulmotype_result_assignment_RTK.txt",
                         sep = "\t", header = T)
row.names(Pulmotypes) <- str_replace(row.names(Pulmotypes), "XX", "X" )
row.names(Pulmotypes) <- str_replace(row.names(Pulmotypes), "X", "" )

Pulmotypes <- Pulmotypes[row.names(Pulmotypes) %in% row.names(meta_make_congruence), ]

meta_make_congruence <- merge(meta_make_congruence, Pulmotypes, by = 0)
meta_make_congruence <- meta_make_congruence[, -c(34:38)]
colnames(meta_make_congruence)[ncol(meta_make_congruence)] <- "Pulmotype"
row.names(meta_make_congruence) <- meta_make_congruence$Row.names
meta_make_congruence$Row.names <- NULL
colnames(meta_make_congruence)[1] <- "Donor"
row.names(meta_make_congruence) <- paste0( "X", row.names(meta_make_congruence))
meta_make_congruence$ID <- row.names(meta_make_congruence) #put ID on first place 
meta_make_congruence <- meta_make_congruence[, c(ncol(meta_make_congruence),1: ncol(meta_make_congruence)-1)] #put ID on first place 



write.table (meta_make_congruence, 
             file = "~/Dropbox/lung/ThedaTill/make_Congruence_perl/CF_only/meta_make_congruence.tsv", 
             quote = F, sep = "\t", row.names = F) # remove rownames

con_improved <- read.table(file = "~/Dropbox/lung/ThedaTill/make_Congruence_perl/CF_only/All_congruence_20210511.tsv",
                           header = T, sep = "\t", check.names = F)


con_3 <- con_improved
con_3$SamePulmotype[con_3$SamePulmotype== "NO"] <- "Different\nCommunity"
con_3$SamePulmotype[con_3$SamePulmotype== "YES"] <- "Same\nCommunity"

con_3$SameDonor[con_3$SameDonor== "YES"] <- "Same\nDonor"
con_3$SameDonor[con_3$SameDonor== "NO"] <- "Diffrent\nDonor"

con_3$SameSampleSite[con_3$SameSampleSite== "YES"] <- "Same\nSwab Site"
con_3$SameSampleSite[con_3$SameSampleSite== "NO"] <- "Diffrent\nSwab Site"

con_3$SameAsthma[con_3$SameAsthma== "YES"] <- "Same\nAsthma Status"
con_3$SameAsthma[con_3$SameAsthma== "NO"] <- "Diffrent\nAsthma Status"

save(con_3, file = "~/Dropbox/lung/ThedaTill/intermediate_files/congurance_metadata_chiq.r")
aFrame_s <- melt(t(apply(table (con_3$SameSampleSite, con_3$SamePulmotype),
                         1, function(x) x/sum(x))))
#aFrame$value <- c(paste0(round((aFrame$value*100), digits = 1)," %"))
#print(paste0(round((aFrame$value*100), digits = 1)," %"))
heatmapconS <- ggplot (aFrame_s, aes (x = Var1, y = Var2)) +
  theme_classic () +
  theme (axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
         legend.position = "none", 
         axis.title = element_blank(), 
         axis.ticks.y = element_blank(),
         axis.text.y = element_blank(),
         axis.line.y = element_blank(), 
         plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  geom_tile (aes (fill = value),colour = "black") +
  geom_text (aes (label = scales::percent(value, accuracy = 1))) +
  scale_fill_gradient(high = "#d7301f", low = "#fff7ec") 

aFrame_s <- melt(t(apply(table (con_3$SameSampleSite, con_3$SamePulmotype),
                         1, function(x) x/sum(x))))
#aFrame$value <- c(paste0(round((aFrame$value*100), digits = 1)," %"))
#print(paste0(round((aFrame$value*100), digits = 1)," %"))
aFrame_P <- melt(t(apply(table (con_3$SameDonor, con_3$SamePulmotype),
                         1, function(x) x/sum(x))))

heatmapconP <- ggplot (aFrame_P, aes (x = Var1, y = Var2)) +
  theme_classic () +
  theme (axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
         legend.position = "none", 
         axis.title = element_blank()) +
  geom_tile (aes (fill = value),colour = "black") +
  geom_text (aes (label = scales::percent(value, accuracy = 1))) +
  scale_fill_gradient(high = "#d7301f", low = "#fff7ec") 

plot <- plot_grid(heatmapconP, heatmapconS, align = "h", rel_widths = (c(1,0.7)))
plot
ggsave(filename = "~/Dropbox/lung/ThedaTill/output/Updated_Figure/20210527_heatmapCombined_Donor_SwabSite_Pulmo_new.svg", 
       plot = plot, device = "svg", width = 3, height = 2)

#### Bulid the heatmaps
aFrame_Asthma <- melt(t(apply(table (con_3$SamePulmotype, con_3$SameAsthma),
                         1, function(x) x/sum(x))))

### Build heatmap from the metadatafile 
meta <- meta_make_congruence
meta$Pulmotype[meta$Pulmotype == 1] <- "Community\nType 1"
meta$Pulmotype[meta$Pulmotype == 2] <- "Community\nType 2"

meta$Asthma[meta$Asthma == 0] <- "Control"
meta$Asthma[meta$Asthma == 1] <- "Asthma"

save(meta, file = "~/Dropbox/lung/ThedaTill/intermediate_files/Chiqs_CF_meta.r")


aFrame_1 <- melt(t(apply(table (meta$Asthma, meta$Pulmotype),
                       1, function(x) x/sum(x))))
heatmapA <- ggplot (aFrame_1, aes (x = Var1, y = Var2)) +
 theme_classic () +
 theme (axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "none", 
        axis.title = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(), 
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
 geom_tile (aes (fill = value), colour = "black") + # add boarders around heatmap
 geom_text (aes (label = scales::percent(value, accuracy = 1))) +
 scale_fill_gradient(high = "#d7301f", low = "#fff7ec") 

print(heatmapA)
print(chisq.test (table (meta$Asthma, meta$Pulmotype)))
# data:  table(meta$Asthma, meta$Pulmotype)
# X-squared = 3.8975, df = 1, p-value = 0.04836
### if we want do fdr copy all p-values from all chiq test together and run p-adjust
ggsave(filename = "~/Dropbox/lung/ThedaTill/output/Updated_Figure/Heatmap_AirwaysType_ASTHMA.svg",
       plot = heatmapA, device = "svg", width = 2.5, height = 1.5)

### load old heatmaps from the cleaned_pulmotype_script
load(file = "~/Dropbox/lung/ThedaTill/intermediate_files/heatmapconS.r")
load(file = "~/Dropbox/lung/ThedaTill/intermediate_files/heatmapconP.r")
load(file = "~/Dropbox/lung/ThedaTill/intermediate_files/heatmapStudy.r")
load(file = "~/Dropbox/lung/ThedaTill/intermediate_files/heatmapSwabsite.r")

load(file = "~/Dropbox/lung/ThedaTill/intermediate_files/aFrame_study.r")
aFrame_study$Var2 <- as.character(aFrame_study$Var2)
aFrame_study$Var2[aFrame_study$Var2 == "SVH"] <- "Community\nType 1"
aFrame_study$Var2[aFrame_study$Var2 == "SVP"] <- "Community\nType 2"
aFrame_study$Var1 <- as.character(aFrame_study$Var1)
aFrame_study$Var1[aFrame_study$Var1 == "Celticfire"] <- "Throat swabs\nCelticfire"
  heatmapStudy <- ggplot (aFrame_study, aes (x = Var1, y = Var2)) +
  theme_classic () +
  theme (axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
         legend.position = "none",
         axis.title = element_blank(), 
         plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  geom_tile (aes (fill = value), colour = "black") + # add boarders around heatmap
  geom_text (aes (label = scales::percent(value, accuracy = 1))) +
  scale_fill_gradient(high = "#d7301f", low = "#fff7ec") 

load(file = "~/Dropbox/lung/ThedaTill/intermediate_files/aFrame.r")
aFrame$Var1 <- sub( "Left Lower Lobe","Left\nLower Lobe", aFrame$Var1)
aFrame$Var1 <- sub( "Left Upper Lobe","Left\nUpper Lobe", aFrame$Var1)

heatmapSwabsite <- ggplot (aFrame, aes (x = Var1, y = Var2)) +
  theme_classic () +
  theme (axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
         legend.position = "none",
         axis.title = element_blank(), 
         axis.ticks.y = element_blank(),
         axis.text.y = element_blank(),
         axis.line.y = element_blank(),
         plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  geom_tile (aes (fill = value), colour = "black") + # add boarders around heatmap
  geom_text (aes (label = scales::percent(value, accuracy = 1))) +
  scale_fill_gradient(high = "#d7301f", low = "#fff7ec") 

### merg all three 
plot <-
  plot_grid(
    heatmapStudy,
    heatmapSwabsite,
    heatmapA,
    align = "h",
    nrow = 1,
    rel_widths = (c(1.25, 1.2, 0.85))
  )
plot

ggsave(
  filename = "~/Dropbox/lung/ThedaTill/output/Updated_Figure/20210527_heatmaps_Community_new.svg",
  plot = plot,
  device = "svg",
  width = 4.8,
  height = 1.9
)
#### Chisq.test this is only in the CF dataset
#### install package 
#### install.packages("devtools")
# library(devtools)
# install_github("dustinfife/fifer")
meta_chisq <- meta
meta_chisq <- fastDummies::dummy_cols(meta_chisq, select_columns = "SampleSite")
colnames(meta_chisq)[which(names(meta_chisq) == "SampleSite_1")] <- "Throat"
colnames(meta_chisq)[which(names(meta_chisq) == "SampleSite_2")] <- "Left Lower Lobe"
colnames(meta_chisq)[which(names(meta_chisq) == "SampleSite_3")] <- "Left Upper Lobe"
colnames(meta_chisq)[which(names(meta_chisq) == "SampleSite_4")] <- "Lingula"

p <- chisq.test (table (meta_chisq$Asthma, meta_chisq$Pulmotype))
## X-squared = 3.8975, df = 1, p-value = 0.04836
print(chisq.test(table(meta_chisq$Throat,meta_chisq$Pulmotype)))
## X-squared = 29.533, df = 1, p-value = 5.499e-08
print(chisq.test(table(meta_chisq$`Left Lower Lobe`,meta_chisq$Pulmotype)))
## X-squared = 4.7494, df = 1, p-value = 0.02931
print(chisq.test(table(meta_chisq$`Left Upper Lobe`,meta_chisq$Pulmotype)))
## X-squared = 15.165, df = 1, p-value = 9.849e-05
print(chisq.test(table(meta_chisq$Lingula,meta_chisq$Pulmotype)))
### X-squared = 6.4064e-32, df = 1, p-value = 1



### alle alpha div plotten. Nach Asthma, Swabsite, Pulmotype 

### get different div. measurements from RTK run 
load(file = "~/Dropbox/lung/ThedaTill/intermediate_files/20210324_rtk_OTUs.r")
### get the other div.
shannon <- get.mean.diversity(rtk_OTU, div = "shannon") 
simpson   <- get.mean.diversity(rtk_OTU, div = "simpson")
invsimpson   <- get.mean.diversity(rtk_OTU, div = "invsimpson")
richness   <- get.mean.diversity(rtk_OTU, div = "richness")
eveness   <- get.mean.diversity(rtk_OTU, div = "eveness")
### transfer into dataframe
simpson <- as.data.frame(simpson)
invsimpson <- as.data.frame(invsimpson)
richness <- as.data.frame(richness)
eveness <- as.data.frame(eveness)
shannon <- as.data.frame(shannon)
#### add sample ID to Shannon use the RTK output
load(file = "~/Dropbox/lung/ThedaTill/intermediate_files/20210324_rtk_normalized_OTUs.r")
row.names(simpson) <- colnames(rtk_use) 
row.names(invsimpson) <- colnames(rtk_use) 
row.names(richness) <- colnames(rtk_use) 
row.names(eveness) <- colnames(rtk_use) 
row.names(shannon) <- colnames(rtk_use)

#### merge all the diveristy measurements in one dataframe 
merge.all <- function(shannon,simpson,invsimpson,richness,eveness, by = "row.names") {
 L <- list(shannon,simpson,invsimpson,richness,eveness)
 for (i in seq_along(L)) {
  shannon <- merge(shannon, L[[i]], by = by)
  rownames(shannon) <- shannon$Row.names
  shannon$Row.names <- NULL
 }
 return(shannon)
}
RTK_div <- merge.all(shannon,simpson,invsimpson,richness,eveness)
RTK_div$shannon.y <- NULL
meta_RTK <- merge(RTK_div, meta, by = 0)
row.names(meta_RTK) <- meta_RTK$Row.names
meta_RTK$Row.names <- NULL

save(RTK_div, file= "~/Dropbox/lung/ThedaTill/intermediate_files/RTK_div.r")
### Pimp the names for diversities
colnames(meta_RTK)[which(names(meta_RTK) == "shannon.x")] <- "shannon"
colnames(meta_RTK)[which(names(meta_RTK) == "Chao1_diversity")] <- "chao1"
colnames(meta_RTK)[which(names(meta_RTK) == "pielou_evenness")] <- "pielou"

######## first run with Pulmotypes #########
meta_RTK_Pulmotype <- meta_RTK[, c(1:5, 31, 39)]
meta_RTK_Pulmotype$Pulmotype[meta_RTK_Pulmotype$Pulmot == "Airway Community Type 1"] <- "Type 1"
meta_RTK_Pulmotype$Pulmotype[meta_RTK_Pulmotype$Pulmot == "Airway Community Type 2"] <- "Type 2"

alphaPlots <- list()
for (i in colnames(meta_RTK_SA)) {
 
 if (i == "SampleSite") {
  next
 }
 print(i)
 multitest <- add_xy_position(pairwise_wilcox_test(data = meta_RTK_SA, 
                                                   formula = formula(paste0(i, " ~ SampleSite")),  # names the varaiabele as character thingy
                                                   p.adjust.method = "fdr"))
 multitest$p.adj.signif[multitest$p.adj < 0.001] <- "***"
 multitest$p.adj.signif[multitest$p.adj < 0.01] <- "**"
 multitest$p.adj.signif[multitest$p.adj < 0.1] <- "*"
 
 alphaPlots[[i]] <- ggviolin(meta_RTK_Swabsite, 
                             x = "SampleSite", 
                             y = i, 
                             fill = "SampleSite", 
                             draw_quantiles = c(0.25, 0.5, 0.75)) +
  theme (axis.text.x = element_text (angle = 90, vjust = 0.5, hjust = 1, size = 7),
         axis.title.x = element_blank()) +
  
  
  stat_pvalue_manual(multitest, hide.ns = T)
 
}

title <- ggdraw() + draw_label(paste0("alpha diversity, richness, evenness: "))
subtitle <- ggdraw() + draw_label("pairwise wilcox.test, FDR:  = 0.1, * = 0.001, *** = 0.0001 ", size = 10)
title <- plot_grid(title, subtitle, rel_heights = c(1, 0.7), nrow = 2)
plotplot <- plot_grid(alphaPlots[[1]] + theme(legend.position = "none"), 
                      alphaPlots[[2]] + theme(legend.position = "none"), 
                      alphaPlots[[3]] + theme(legend.position = "none"), 
                      alphaPlots[[4]] + theme(legend.position = "none"), 
                      alphaPlots[[5]] + theme(legend.position = "none"),
                      alphaPlots[[6]] + theme(legend.position = "none"))
wholeDiv_Swabsite <- plot_grid(title, 
                                    plotplot, 
                                    rel_heights = c(0.1, 1), 
                                    nrow = 2)
print(wholeDiv_Swabsite)

ggsave(filename = "~/Dropbox/lung/ThedaTill/output/Updated_Figure/01_Alpha_SAMPELSITE.svg",
       plot = wholeDiv_Swabsite, device = "svg", width = 9, height = 8)

######## first run with Swabsite #########
meta_RTK_Swabsite <- meta_RTK[, c(1:5, 31, 8)]
## add real names!!!!

######## first run with Asthma #########
meta_RTK_Asthma <- meta_RTK[, c(1:5, 31, 15)]
### only include one sample per patient other we increase bias
### create loop to subset for each samplesite
meta_RTK$SampleSite[meta_RTK$SampleSite == 1] <- "Throat"
meta_RTK$SampleSite[meta_RTK$SampleSite == 2] <-  "Left Lower Lobe"
meta_RTK$SampleSite[meta_RTK$SampleSite == 3] <- "Left Upper Lobe"
meta_RTK$SampleSite[meta_RTK$SampleSite == 4] <- "Lingula"

### use meta_RTK_SA to just have Asthma and Swabsite
meta_RTK_SA$SampleSite[meta_RTK_SA$SampleSite == 1] <- "Throat"
meta_RTK_SA$SampleSite[meta_RTK_SA$SampleSite == 2] <-  "Left Lower Lobe"
meta_RTK_SA$SampleSite[meta_RTK_SA$SampleSite == 3] <- "Left Upper Lobe"
meta_RTK_SA$SampleSite[meta_RTK_SA$SampleSite == 4] <- "Lingula"


meta_RTK_SA$Pulmotype <- meta_RTK$Pulmotype


for (k in unique(meta_RTK_SA$SampleSite)) {
  # if (k == "Lingula"){next}
  meta_subset <- subset(meta_RTK_SA, SampleSite == k)# subset meta_RTK for sampleSite using k as value for samplesite
  print(k)
  for (j in c("Asthma", "Pulmotype")) {
    alphaPlots <- list()
    print(j)
    for (i in colnames(meta_subset)) {
      
      if (i %in% c("SampleSite","Asthma", "Pulmotype")) {
        next
      }
      print(i)
      multitest <- add_xy_position(pairwise_wilcox_test(data = meta_subset, 
                                                        formula = formula(paste0(i, " ~ ", j)),  # names the varaiabele as character thingy
                                                        p.adjust.method = "fdr"))
      multitest$p.adj.signif[multitest$p.adj < 0.001] <- "***"
      multitest$p.adj.signif[multitest$p.adj < 0.01] <- "**"
      multitest$p.adj.signif[multitest$p.adj < 0.1] <- "*"
      
      alphaPlots[[i]] <- ggviolin(meta_subset, 
                                  x = j, 
                                  y = i, 
                                  fill = j, 
                                  draw_quantiles = c(0.25, 0.5, 0.75)) +
        theme (axis.text.x = element_text (angle = 90, vjust = 0.5, hjust = 1, size = 7),
               axis.title.x = element_blank()) +
        
        
        stat_pvalue_manual(multitest, hide.ns = T)
      
    } 
    title <- ggdraw() + draw_label(paste0("alpha diversity, richness, evenness: ", k))
    subtitle <- ggdraw() + draw_label("pairwise wilcox.test, FDR:  = 0.1, * = 0.001, *** = 0.0001 ", size = 10)
    title <- plot_grid(title, subtitle, rel_heights = c(1, 0.7), nrow = 2)
    plotplot <- plot_grid(alphaPlots[[1]] + theme(legend.position = "none"), 
                          alphaPlots[[2]] + theme(legend.position = "none"), 
                          alphaPlots[[3]] + theme(legend.position = "none"), 
                          alphaPlots[[4]] + theme(legend.position = "none"), 
                          alphaPlots[[5]] + theme(legend.position = "none"),
                          alphaPlots[[6]] + theme(legend.position = "none"))
    wholeDiv_AS <- plot_grid(title, 
                             plotplot, 
                             rel_heights = c(0.1, 1), 
                             nrow = 2)
    print(wholeDiv_AS)
    
    ggsave(filename = paste0("~/Dropbox/lung/ThedaTill/output/Updated_Figure/01_Alpha_", j, "_", k, ".svg"),
           plot = wholeDiv_AS, device = "svg", width = 9, height = 8)
  } 
}









