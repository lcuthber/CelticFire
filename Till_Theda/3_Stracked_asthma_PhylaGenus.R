library (lme4)
library (gtools)
library (ggplot2)
library (reshape2)
library (vegan)
library (plyr)
library (lmtest)
library (orddom)
library (ggrepel)
library (stringr)
library(dplyr)
library(ggpubr)
# Use phylum and Genus 
biom_file <- import_biom("~/Dropbox/lung/ThedaTill/raw_data/otu_table_json_copy.biom", sep = "")
load(file = "~/Dropbox/lung/ThedaTill/intermediate_files/20210324_joined_genus.R")

phylum <- tax_glom(physeq = biom_file, taxrank = "Rank2")
save(phylum, file = "~/Dropbox/lung/ThedaTill/intermediate_files/20210324_biom_file_taxGlom_phylum.R")
data_biom_all <- as.data.frame(tax_table(biom_file)) # still includes NA's
data_biom_tax <- as.data.frame(tax_table(phylum)) #biome dataframe mit taxa -> damit die zuordnung klappt beim mergen
data_biom <-as.data.frame(otu_table(phylum)) #biome datafram mit otu
# umbennen der Ranks in phylum etc.
colnames(data_biom_tax)<-c("kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# genus_otu_short <- as.data.frame(otu_table(data_biom_tax))
data_biom_tax$hash <- rownames(data_biom_tax)
data_biom$hash <- rownames(data_biom)
joined_phylum <- left_join(data_biom_tax, data_biom)
colnames(joined_phylum) <- sub(pattern = "_.+", replacement = "", x = colnames(joined_phylum), )
joined_phylum$hash <- NULL
joined_phylum[, c(1,3:7) ]<- NULL
rownames(joined_phylum) <- make.names(joined_phylum$Phylum, unique = T)
rownames(joined_phylum) <- sub("p__", "", rownames(joined_phylum))
joined_phylum$Phylum <- NULL
joined_phylum <- joined_phylum[order(rownames(joined_phylum)), ]
joined_phylum <- t(joined_phylum)

save(joined_phylum, file = "~/Dropbox/lung/ThedaTill/intermediate_files/20210527_joined_phylum.r")

#Metadata
load("~/Dropbox/lung/ThedaTill/raw_data/metadata_red.r")

#Phylum tranformieren
Phylum <- joined_phylum
Phylum <- Phylum[rownames(Phylum) %in% rownames(metadata_red), ]

#Phylum <- t(Phylum)
ad.index.keep <- which(colSums(Phylum)*100/(sum(Phylum)) > 0.01) # reduce dataset 
Phylum <- Phylum[, ad.index.keep]
Phylum <- Phylum[order(row.names(Phylum)), ]
metadata_red <- metadata_red[order(row.names(metadata_red)), ]

# Phylum <- Phylum[order(row.names(Phylum)), ]
# points to - hack
row.names(Phylum)
class(Phylum)
Phylum <- as.data.frame(Phylum)
Phylum$Smpl <- row.names(Phylum)
md <- metadata_red
md$Smpl <- row.names(md)
md <- md[c("Asthma", "Smpl")] # Nur status und smp,l id rausnhemen 
Abundance <- left_join(Phylum, md)
Abundance$Smpl <- NULL

# Phylum$Smpl <-NULL
# Abundance gucken 
# DualStatus in 4 Gruppen einteil
Abundance_Asthma <- colSums(Abundance[Abundance$Asthma == 1, - ncol(Abundance) ])
Abundance_Control <- colSums(Abundance[Abundance$Asthma == 0 , - ncol(Abundance) ])
#merging
Abundance_merged <- rbind(Abundance_Asthma, Abundance_Control)
# Abundance_merged <-t(Abundance_merged)

rel_Phylum <- apply(Abundance_merged, 1, function(x) x/sum(x)) # rows 1 als prozente
# er tranformiert automatisch 
rel_Phylum <- rel_Phylum[names(sort(rowSums(rel_Phylum), decreasing = TRUE)[1:6]), ] #nur 1-6 raussuchen wir wollen die größten 6 drinnen behalten 

# delete Bacteria in front of Phylum names

# rownames(rel_Phylum) <- sub("Bacteria;",  "", rownames(rel_Phylum))
rel_Phylum <- as.data.frame(rel_Phylum)

#id Phylum erstellen row.names must be phylum colnames COPD and shit
# rel_Phylum <-t(rel_Phylum)
rel_Phylum$Phylum <- row.names(rel_Phylum)

for (i in 1:2) {
 rel_Phylum[, i] <- as.numeric(rel_Phylum[, i])
} # macht aus den abundance für die bacteria numerical zahlen

colnames(rel_Phylum) <- c("Asthma", "Control", "Phylum")
molten_rel_Phylum <- melt(data = rel_Phylum, id.vars = "Phylum")
molten_rel_Phylum$labels <- scales::percent(molten_rel_Phylum$value)
molten_rel_Phylum$variable <- str_replace((molten_rel_Phylum$variable), "Abundance_", "")
# appended_molten <- read.csv(file = "molten_Phylum_intermediate.cvs", sep = "\t")

pControls <- ggplot (data = molten_rel_Phylum, aes(x = variable, y = value, fill = Phylum)) +
 geom_bar(stat = "identity", color = "black") +
 # geom_text(aes(label = labels), vjust=1.6) +
 theme_classic() +
 theme(axis.text.x = element_text(angle = 90),
       axis.title.x = element_blank()) +
 ylab("relative abundance")
print(pControls)
# ggsave(filename = "Genus_controls_theoretical.svg", plot = pControls, device = "svg", width = 5, height = 5)
#Plot Abundance x -aches und DualStatus y- aches 

# Genus shit 

#Phylum tranformieren
Genus <- joined_genus
Genus <- Genus[rownames(Genus) %in% rownames(metadata_red), ]
#Phylum <- t(Phylum)
ad.index.keep <- which(colSums(Genus)*100/(sum(Genus)) > 0.01) # reduce dataset 
Genus <- Genus[, ad.index.keep]
Genus <- Genus[order(row.names(Genus)), ]

metadata_red <- metadata_red[order(row.names(metadata_red)),]

# Phylum <- Phylum[order(row.names(Phylum)), ]
# points to - hack
row.names(Genus)
class(Genus)
Genus <- as.data.frame(Genus)
Genus$Smpl <- row.names(Genus)
md <- metadata_red
md$Smpl <- row.names(md)
md <- md[c("Asthma", "Smpl")] # Nur status und smp,l id rausnhemen 
Abundance <- left_join(Genus, md)
Abundance$Smpl <- NULL

# Phylum$Smpl <-NULL
# Abundance gucken 
# DualStatus in 4 Gruppen einteil
Abundance_Asthma <- colSums(Abundance[Abundance$Asthma == 1, - ncol(Abundance) ])
Abundance_Control <- colSums(Abundance[Abundance$Asthma == 0 , - ncol(Abundance) ])
#merging
Abundance_merged <- rbind(Abundance_Asthma, Abundance_Control)
# Abundance_merged <-t(Abundance_merged)

rel_Genus <- apply(Abundance_merged, 1, function(x) x/sum(x)) # rows 1 als prozente
# er tranformiert automatisch 
rel_Genus <- rel_Genus[names(sort(rowSums(rel_Genus), decreasing = TRUE)[1:15]), ] #nur 1-6 raussuchen wir wollen die größten 6 drinnen behalten 

# delete Bacteria in front of Phylum names

# rownames(rel_Phylum) <- sub("Bacteria;",  "", rownames(rel_Phylum))
rel_Genus <- as.data.frame(rel_Genus)

#id Phylum erstellen row.names must be phylum colnames COPD and shit
# rel_Phylum <-t(rel_Phylum)
rel_Genus$Genus <- row.names(rel_Genus)

for (i in 1:2) {
 rel_Genus[, i] <- as.numeric(rel_Genus[, i])
} # macht aus den abundance für die bacteria numerical zahlen

colnames(rel_Genus) <- c("Asthma", "Control", "Genus")
molten_rel_Genus <- melt(data = rel_Genus, id.vars = "Genus")
molten_rel_Genus$labels <- scales::percent(molten_rel_Genus$value)
molten_rel_Genus$variable <- str_replace((molten_rel_Genus$variable), "Abundance_", "")
# appended_molten <- read.csv(file = "molten_Phylum_intermediate.cvs", sep = "\t")

pControl <- ggplot (data = molten_rel_Genus, aes(x = variable, y = value, fill = Genus)) +
 geom_bar(stat = "identity", color = "black") +
 # geom_text(aes(label = labels), vjust=1.6) +
 theme_classic() +
 theme(axis.text.x = element_text(angle = 90),
       axis.title.x = element_blank()) +
 ylab("relative abundance")
print(pControl)

plot <- plot_grid(pControls, pControl, nrow = 1, align = "hv")
ggsave(filename = "~/Dropbox/lung/ThedaTill/output/Updated_Figure/abundance_PhylaGenus_Asthma.svg",
       plot = plot, device = "svg", width = 6, height = 4)

### maybe also do that in the pulmotype ??
###### Nur plot auf Genus lvl mit den Farben von Phylum namen mit ggrepel 
data_biom_red <- data_biom_all
data_biom_red$Rank6 <- gsub( "g__", "", data_biom_red$Rank6)        
data_biom_red <- data_biom_red[data_biom_red$Rank6 %in% molten_rel_Genus$Genus, ]
which(data_biom_red$Rank6 == molten_rel_Genus)
data_biom_red$Rank6[molten_rel_Genus$Genus, "Rank6"]
x <- match(molten_rel_Genus$Genus, data_biom_red$Rank6)
y <- data_biom_red[x, ]
y$Rank2 <- gsub( "p__", "", y$Rank2)        
molten_rel_Genus$Phylum <- y$Rank2


write.table(molten_rel_Genus, file = "~/Dropbox/lung/ThedaTill/intermediate_files/molten_rel_Genus.tsv", sep = "\t", quote = F)
molten_rel_Genus <- read.table("~/Dropbox/lung/ThedaTill/intermediate_files/molten_rel_Genus.tsv", sep = "\t", header = T)
pControlx <- ggplot (data = molten_rel_Genus, aes(x = variable, y = value, fill = Phylum, width = 0.7)) +
                     ylim(-0.1,1.0) +
        geom_bar(stat = "identity", color = "black") +
        # geom_text(aes(label = labels), vjust=1.6) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, size = 14),
              axis.title.x = element_blank(),
              axis.text.y = element_text( size = 14),
              axis.title.y = element_text(size = 14),
              legend.position = "bottom", plot.margin = unit(c(1,1,1,1), "cm"), legend.text = element_text(size = 14)) +
        scale_x_discrete(limits = c("Control", "Asthma", "Genus")) +
        geom_text_repel(
                na.rm = TRUE,
                data = subset(molten_rel_Genus, variable == "Asthma"),
                label = subset(molten_rel_Genus, variable == "Asthma")$Genus,
                force = 80,
                fontface = "italic",
                size = 5,
                max.overlaps = 170, position = position_stack(vjust = 0.5),
                ylim = c(-0.1, 1.5), xlim = c(2.6,2.8), box.padding = 0.075
        ) +
        ylab("relative abundance")
print(pControlx)

# labels can easily be stacked by using position = position_stack(vjust = 0.5) in geom_text.

ggsave(filename = "~/Dropbox/lung/ThedaTill/output/Updated_Figure/abundance_Phyla-Genus_Asthma.svg",
       plot = pControlx, device = "svg", width = 3, height = 6)


