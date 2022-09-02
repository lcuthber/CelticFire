# load packages
library(DirichletMultinomial)
library(GGally)
library(ggplot2)
library(cowplot)
library(rtk)
library(edgeR)
library(vegan)
library(tidyverse)
library(dplyr)
library(phyloseq)
library(metadeconfoundR)

if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

BiocManager::install("KEGGREST")

library(KEGGREST)

metadata <- read.table("~/Dropbox/lung/ThedaTill/raw_data/cf_formatted_mapfile.csv",
                       sep = ",",
                       header = TRUE, stringsAsFactors = F)
# there is an issues with 0 in the sampleID names
## separate metadata by samplesID with and w/O r in the name
metadata_endr <- metadata[grep(pattern = "*r", x = metadata$SampleID), ]
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
metadata$Ethnicity <- paste0("Ethnicity_", metadata$Ethnicity) # makes Enthnicity 1,2 and 3
# replace Ethnicity form contiouse to pseudo binary 
for (i in unique(metadata$Ethnicity)) {
 binaryDummy <- rep(0, length(metadata$Ethnicity))
 binaryDummy[metadata$Ethnicity == i] <- 1
 metadata[[i]] <- binaryDummy
 colnames(metadata)[ncol(metadata)] <- i
}
metadata$Ethnicity <- NULL
metadata$Study <- NULL # since its only Celticfire
metadata$Study_ID <- paste0("Patient_", metadata$Study_ID) #Donor specify 

# change Samplesite names 
metadata$SampleSite <- paste0("SampleSite_", metadata$SampleSite)
for (i in unique(metadata$SampleSite)) {
 binaryDummy <- rep(0, length(metadata$SampleSite))
 binaryDummy[metadata$SampleSite == i] <- 1
 metadata[[i]] <- binaryDummy
 colnames(metadata)[ncol(metadata)] <- i
}
metadata$SampleSite <- NULL 

# Prepare Genus/ OTUs etc for analysis 
biom_file <- import_biom("~/Dropbox/lung/ThedaTill/raw_data/otu_table_json_copy.biom", sep = "")
OTUs <- as.data.frame(otu_table(biom_file))
# colnames(OTUs) <- paste0("0", colnames(OTUs))
OTUs_cleaned <- OTUs[, colnames(OTUs) %in% rownames(metadata)]
samplesize <- min(colSums(OTUs_cleaned))
rtk_OTU <-
        rtk(
                OTUs_cleaned,
                repeats = 10,
                depth = samplesize,
                ReturnMatrix = 1,
                margin = 2,
                verbose = FALSE,
                threads = 1,
                tmpdir = NULL,
                seed = 0
        )
save(rtk_OTU, file = "~/Dropbox/lung/ThedaTill/intermediate_files/20210324_rtk_OTUs.r")

# get OTU, select raremat in rtk file
rtk_use <-
        as.data.frame(rtk_OTU$raremat) # <- Otus used for the mapping of Isolates later
save(rtk_use, file = "~/Dropbox/lung/ThedaTill/intermediate_files/20210324_rtk_normalized_OTUs.r")

otu_table(biom_file) <-
        otu_table(rtk_use, taxa_are_rows = T) #normalized OTU file in biomefile zurück
otu_table(biom_file)[1:5, 1:5]
#biome datafram mit otus
#data_biom <-as.data.frame(otu_table(biom_file)) #biome datafram mit otus

# if biom file on genus level was already generated, load it. 
#Otherwise compute it and save it for later use
if (file.exists("~/Dropbox/lung/ThedaTill/intermediate_files/20210324_biom_file_taxGlom_phylum.R")) {
 load("~/Dropbox/lung/ThedaTill/intermediate_files/20210324_biom_file_taxGlom_phylum.R")
} else {
 genus <- tax_glom(physeq = biom_file, taxrank = "Rank2")
 save(genus, file = "~/Dropbox/lung/ThedaTill/intermediate_files/20210324_biom_file_taxGlom_phylum.R")
}
data_biom_all <- as.data.frame(tax_table(biom_file)) # still includes NA's
data_biom_tax <- as.data.frame(tax_table(genus)) #biome dataframe mit taxa -> damit die zuordnung klappt beim mergen
data_biom <-as.data.frame(otu_table(genus)) #biome datafram mit otu
# umbennen der Ranks in phylum etc.
colnames(data_biom_tax)<-c("kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# genus_otu_short <- as.data.frame(otu_table(data_biom_tax))
data_biom_tax$hash <- rownames(data_biom_tax)
data_biom$hash <- rownames(data_biom)
joined_genus <- left_join(data_biom_tax, data_biom)
colnames(joined_genus) <- sub(pattern = "_.+", replacement = "", x = colnames(joined_genus), )
joined_genus$hash <- NULL
joined_genus[, c(1,3:7) ]<- NULL
rownames(joined_genus) <- make.names(joined_genus$Phylum, unique = T)
rownames(joined_genus) <- sub("p__", "", rownames(joined_genus))
joined_genus$Phylum <- NULL
joined_genus <- joined_genus[order(rownames(joined_genus)), ]
joined_genus <- t(joined_genus)
rownames(joined_genus) <- sub("^X", "", rownames(joined_genus))
load(file = "~/Dropbox/lung/ThedaTill/intermediate_files/20210324_joined_phylum.R")
# load(file = "~/Dropbox/lung/ThedaTill/intermediate_files/20210324_joined_genus.R")

# Prepare the data for metdadeconfounder, remnove low abundance features
features <- joined_genus
dim(features)
ad.index.keep <- which(colSums(features)*100/(sum(features)) > 0.01)
features <- features[, ad.index.keep]
dim(features)
order(row.names(features)) #check for same order of row.nmaes in both df
#use Study as randome effect -> currentsmorker as first column
metadata <- metadata[, c(3,1,2,4:ncol(metadata))] #new first column -> splate 3 1 und 2 dann der rest
# cut out some variables in the metadata
metadata_red <- metadata[, c(1:5, 7:12, 14, 18, 22, 26:27, 30:36)] #ACHTUNG
metadata_red$Chao1_diversity <- NULL
metadata_red$Ethnicity_1 <- NULL
metadata_red$Ethnicity_2 <- NULL
metadata_red$Ethnicity_3 <- NULL 

save(metadata_red, file = "~/Dropbox/lung/ThedaTill/raw_data/metadata_red.r")
load(file = "~/Dropbox/lung/ThedaTill/raw_data/metadata_red.r")

features <- features[order(row.names(features)), ]
metadata_red <- metadata_red[order(row.names(metadata_red)), ]

metadeconf_Celtic_fullGenus <- MetaDeconfound(featureMat = features,
                                               metaMat = metadata_red,
                                               nnodes = 7,
                                               randomVar = list("+ (1|Study_ID) + (1|Site)", c("Study_ID", "Site")), 
                                               logfile = "~/Dropbox/lung/ThedaTill/intermediate_files/new_metadeconf_fullGenus.log")

# you need to change to rank_2 in the biom_file filtering step above to have 
        #phylum level information in "features"!!
metadeconf_Celtic_fullPhylum <- MetaDeconfound(featureMat = features,
                                              metaMat = metadata_red,
                                              nnodes = 3,
                                              randomVar = list("+ (1|Study_ID) + (1|Site)", c("Study_ID", "Site")), 
                                              logfile = "~/Dropbox/lung/ThedaTill/intermediate_files/20210421_metadeconf_fullPhylum.log")


nice_meta_names <- as.data.frame(nice_meta_names)
nice_meta_names$nice_meta_names[1] <- "Sex"
nice_meta_names$nice_meta_names[2] <- "Study ID"
nice_meta_names$nice_meta_names[3] <- "Site"


plot <- BuildHeatmap(metadeconf_Celtic_fullPhylum, metaVariableNames = nice_meta_names)

save(metadeconf_Celtic_fullPhylum, file = "~/Dropbox/lung/ThedaTill/output/20210421_Metadeconf_Celtic_fullPhylum_final.r")

ggsave(filename = "~/Dropbox/lung/ThedaTill/output/Updated_Figure/20210421_Metadeconf_Celtic_fullPhylum_final.svg", 
       plot = plot, device = "svg", width = 4.75, height = 3.5)

load(file = "~/Dropbox/lung/ThedaTill/intermediate_files/20210527_joined_phylum.r")
features <- joined_phylum
dim(features)
ad.index.keep <- which(colSums(features)*100/(sum(features)) > 0.01)
features <- features[, ad.index.keep]
dim(features)
dim(metadata_red)
features <- features[rownames(features) %in% rownames(metadata_red), ]
# Pimp the names with the new feature 
feature_names <- features
feature_names <- t(feature_names) # run this row.name 5 times to remove all the _
row.names(feature_names) <- sub("_", " ", row.names(feature_names))
row.names(feature_names) <- sub("\\.", " ", row.names(feature_names))
row.names(feature_names) <- sub("\\.", "", row.names(feature_names))
row.names(feature_names) <- sub("\\.", " ", row.names(feature_names))
row.names(feature_names) <- sub("_", " ", row.names(feature_names))
feature_names <- t(feature_names) # human readle names 
colnames(feature_names)[which(names(feature_names) == "Allorhizobium NeorhizobiumPararhizobium Rhizobium")] <- "Allorhizobium Neorhizobium Pararhizobium Rhizobium"

nice_feature_names <- colnames(features) # old names
nice_feature_names <- cbind(nice_feature_names, colnames(feature_names)) # cbind do make df with new and old names side by side 
save(nice_feature_names, file = "~/Dropbox/lung/ThedaTill/intermediate_files/Genus_nice_feature_name.r")

# Changes the metadata
metaVar_Names <- metadata_red
#change colnames 
colnames(metaVar_Names)[which(names(metaVar_Names) == "SampleSite_1")] <- "Throat"
colnames(metaVar_Names)[which(names(metaVar_Names) == "SampleSite_2")] <- "Left Lower Lobe"
colnames(metaVar_Names)[which(names(metaVar_Names) == "SampleSite_3")] <- "Left Upper Lobe"
colnames(metaVar_Names)[which(names(metaVar_Names) == "SampleSite_4")] <- "Lingula"
#now changes the colnames 
colnames(metaVar_Names) <- sub("_", " ",  colnames(metaVar_Names))
colnames(metaVar_Names) <- sub("_", " ",  colnames(metaVar_Names))
colnames(metaVar_Names)[which(names(metaVar_Names) == "Final 16S copies_per_ul")] <- "Final 16S copy/µL"

nice_meta_names <- colnames(metadata_red) # old names
nice_meta_names <- cbind(nice_meta_names, colnames(metaVar_Names)) # cbind do make df with new and old names side by side 
save(nice_meta_names, file = "~/Dropbox/lung/ThedaTill/intermediate_files/Genus_nice_meta_name.r")
#### here is a problem with the names, be aware need to fix
load("~/Dropbox/lung/ThedaTill/output/20210421_Metadeconf_Celtic_fullGenus_final.r")
load("~/Dropbox/lung/ThedaTill/intermediate_files/Genus_nice_feature_name.r")
load("~/Dropbox/lung/ThedaTill/intermediate_files/Genus_nice_meta_name.r")


plot <- BuildHeatmap(metadeconf_Celtic_fullGenus, featureNames = nice_feature_names)

save(plot, file = "~/Dropbox/lung/ThedaTill/intermediate_files/20210422_Metadeconfounder_fullGenus_NameNOOOOOOO_Final.r")


ggsave(filename = "~/Dropbox/lung/ThedaTill/output/20210421_metadeconf_Celtic_fullGenus_nice_Name.svg", 
       plot = plot, device = "svg", width = 8, height = 10)

# Colapse 16S data onto Isolates for more condensed metadeconfoundR input 
# Sum up abundances for all OTUS mapping to an isolate
# use the rarefied data to bin them 20210315 

Isolate <- read.table (file = "~/Dropbox/lung/ThedaTill/Isolates/BusCF_OTUs_allisolates16S.e5.pi78.blastn.sorted.tab", sep = "\t")
dim(Isolate)
length(unique(Isolate$V2))
isolate97 <- Isolate[Isolate$V10 >= 97,]

#OTUs <- read.table (file = "~/Dropbox/lung/ThedaTill/", sep = "\t", header = TRUE)
# biom_file <- import_biom("~/Dropbox/lung/ThedaTill/otu_table_json.biom", sep = "")
# OTUs <- as.data.frame(otu_table(biom_file))
# colnames(OTUs) <- paste0("X", colnames(OTUs))
# OTUs_cleaned <- OTUs[, colnames(OTUs) %in% rownames(metadata_reduced)]

load(file ="~/Dropbox/lung/ThedaTill/intermediate_files/20210324_rtk_normalized_OTUs.r")
OTUs_cleaned <- rtk_use
# for each unique name in V2 of isolate97
mapped_isolates <- as.data.frame(matrix(ncol = ncol(OTUs_cleaned)))
for (i in unique(isolate97$V2)) {
 #find all OTU names matching to isolate name i
 otuNames <- (isolate97[isolate97$V2 == i, 1])
 # extract all matching OTUs from abundance table
 OTUsSub <- OTUs[rownames(OTUs_cleaned) %in% otuNames, ]
 # sum up all rows of OTU subset to get isolate i abundance and add to mapped_isolates data.frame
 mapped_isolates[i, ] <- (colSums(OTUsSub))
}
mapped_isolates <- mapped_isolates[-1, ] # remove first row
colnames(mapped_isolates) <- colnames(OTUs_cleaned)
dim(mapped_isolates)
# all lines with the same species identifying 
# using gsub 


# rest <- OTUs_cleaned[!(row.names(OTUs_cleaned)%in%isolate97$V1), ] #collect the OTUs not matching the Isolates
# colSums(rest)
# mapped_isolates = rbind(mapped_isolates, colSums(rest))
# row.names(mapped_isolates)[118] <- "unmapped" #merge unmapped on mapped_isolates 
# 
# mapped_isolates_rtk <- rtk(mapped_isolates, repeats = 10, depth = samplesize, ReturnMatrix = 1, margin = 2,
#                            verbose = FALSE, threads = 1, tmpdir = NULL, seed = 0)
# 
# # get OTU, select raremat in rtk file 
# mapped_isolates_rtk <- as.data.frame(mapped_isolates_rtk$raremat) # <- Otus used for the mapping of Isolates later 
# save(mapped_isolates_rtk, file ="~/Dropbox/lung/ThedaTill/20201117_Enterotyping/mapped_isolates_rtk.r")

save(mapped_isolates, file ="~/Dropbox/lung/ThedaTill/intermediate_files/20210325_mapped_isolates.r")


sum(!(row.names(OTUs)%in%isolate97$V1)) # means only 700 OTUs from the 16S reads can be mapped on Isolates 
dim(isolate97)

# binn auf species name so its pretty

SpeciesNames <- gsub('^.+\\|([a-z|A-Z]+_[a-z|A-Z]+).+', '\\1', row.names(mapped_isolates)) # filter names 

names_number <- as.data.frame(summary(as.factor(SpeciesNames))) #Summery befehl als df to count the multiople species
binnedSpecies <- as.data.frame(matrix(nrow = ncol(mapped_isolates), ncol = 1)) # leerer df
rownames(binnedSpecies) <- colnames(mapped_isolates)

for (i in unique(SpeciesNames)) {
 binnedSpecies[[i]] <- colSums(mapped_isolates[which(SpeciesNames == i), ])
 
}
binnedSpecies$V1 <- NULL

binnedSpecies <- binnedSpecies[, order(colnames(binnedSpecies))]
colnames(binnedSpecies) <- paste0(rownames(names_number), "_", names_number$`summary(as.factor(SpeciesNames))`)
# für jeden uniquen species namen alle rows summiert und in den dataframe gespeichert, colnames angepasst
save(binnedSpecies, file ="~/Dropbox/lung/ThedaTill/intermediate_files/20210325_binnedSpecies_Isolates.r")
write.table(binnedSpecies, file = "~/Dropbox/lung/ThedaTill/intermediate_files/20210325_binnedSpecies_Isolates.csv")
write.table(mapped_isolates, file = "~/Dropbox/lung/ThedaTill/intermediate_files/20210325_mapped_isolates.csv", sep = "\t")


# remove low abundance isolates from analysis
features <- binnedSpecies # for the metadeconfounder
rownames(features) <- sub("^X", "", rownames(features))
# now remove all the Busselton since they dont match here 
metadata_isolates <- metadata_red
# metadata_isolates <- as.data.frame(t(metadata_isolates))
# features <- as.data.frame(t(features)) # stand 13.7 nur ein schnelles helperchen
features <- features[row.names(features)%in%row.names(metadata_isolates), ]
dim(features)
ad.index.keep <- which(colSums(features)*100/(sum(features)) > 0.01)
features <- features[, ad.index.keep]
dim(features)
features <- features[order(row.names(features)), ]
metadata_isolates <- metadata_isolates[order(row.names(metadata_isolates)), ]
#normalizing 

metadeconf_Isolates_binned <- MetaDeconfound(featureMat = features,
                                             metaMat = metadata_isolates,
                                             nnodes = 3,
                                             randomVar = list("+ (1|Study_ID) + (1|Site)", c("Study_ID", "Site")),
                                             logfile = "~/Dropbox/lung/ThedaTill/intermediate_files/20210325_metadeconf_Isolate_binned.log")

plot <- BuildHeatmap(metadeconf_Isolates_binned)
ggsave(filename = "~/Dropbox/lung/ThedaTill/output/metadeconf_Isolates_binned_First.svg",
       plot = plot, device = "svg", width = 8, height = 10)

# read cpd:KO pathways 
Iso_Func <- read.table( "~/Dropbox/lung/ThedaTill/intermediate_files/out_mapped_KO2Metabolites_20210317.r",
             sep = "\t", header = TRUE, row.names = 1)


features <- t(Iso_Func)
rownames(features) <- sub("^X", "", rownames(features))
features <- features[row.names(features)%in%row.names(metadata_isolates), ]
metadata_isolates <- metadata_isolates[row.names(metadata_isolates)%in%row.names(features), ] # 285 sampls in metadata only 280 in iso_KO


dim(features)
ad.index.keep <- which(colSums(features)*100/(sum(features)) > 0.01)
features <- features[, ad.index.keep]
dim(features)
features <- features[order(row.names(features)), ]
metadata_isolates <- metadata_isolates[order(row.names(metadata_isolates)), ]
#normalizing 
metadeconf_Isolates_Metabo <- MetaDeconfound(featureMat = features,
                                             metaMat = metadata_isolates,
                                             nnodes = 3,
                                             randomVar = list("+ (1|Study_ID) + (1|Site)", c("Study_ID", "Site")),
                                             logfile = "~/Dropbox/lung/ThedaTill/intermediate_files/20210421_metadeconf_Isolate_Metabo.log")

load( file = "~/Dropbox/lung/ThedaTill/intermediate_files/20210421_metadeconf_Isolates_Metabo.r")

# changes metbolite names 
KeggComp <- read.table("~/Dropbox/lung/ThedaTill/raw_data/KEGGCompoundNames.tsv", 
                       sep = "\t", quote = "") # So there was a problem with the 5' in the names 
# use quote = "" to get ride of this problem 

existing_metabolites <- KeggComp[KeggComp[, 1] %in% rownames(metadeconf_Isolates_Metabo$status), ]
rownames(metadeconf_Isolates_Metabo$status)[!(rownames(metadeconf_Isolates_Metabo$status)%in% KeggComp[, 1])]
x

longNames <- merge(x = metadeconf_Isolates_Metabo$status, y = existing_metabolites, by.x = 0, by.y = 1) # names are now in longNames$V2
write.table(longNames, file = "~/Dropbox/lung/ThedaTill/intermediate_files/Metabolites_long.tsv", sep = "\t", quote = F)


# longNames$long <- make.names(longNames$long, unique = TRUE)
### assigning weird long rownmaes leads to errors in buildHeatmap()
longNames$V2 <- gsub(';.*', '', longNames$V2) 
rownames(metadeconf_Isolates_Metabo$Ps) <- longNames$V2
rownames(metadeconf_Isolates_Metabo$Qs) <- longNames$V2
rownames(metadeconf_Isolates_Metabo$Ds) <- longNames$V2
rownames(metadeconf_Isolates_Metabo$status) <- longNames$V2

save(metadeconf_Isolates_Metabo, file ="~/Dropbox/lung/ThedaTill/Picrust/20210421_metadeconf_Isolates_metabolites_Final.r" )

# Chnages Names to human readable 
# Changes the metadata
metaVar_Names <- metadata_isolates
#change colnames 
colnames(metaVar_Names)[which(names(metaVar_Names) == "SampleSite_1")] <- "Throat"
colnames(metaVar_Names)[which(names(metaVar_Names) == "SampleSite_2")] <- "Left Lower Lobe"
colnames(metaVar_Names)[which(names(metaVar_Names) == "SampleSite_3")] <- "Left Upper Lobe"
colnames(metaVar_Names)[which(names(metaVar_Names) == "SampleSite_4")] <- "Lingula"
colnames(metaVar_Names)[which(names(metaVar_Names) == "simpson evenness")] <- "simpson evenness"
#now changes the colnames 
colnames(metaVar_Names) <- sub("_", " ",  colnames(metaVar_Names))
colnames(metaVar_Names) <- sub("_", " ",  colnames(metaVar_Names))
colnames(metaVar_Names)[which(names(metaVar_Names) == "Final 16S copies_per_ul")] <- "Final 16S copy/µL"

nice_meta_names <- colnames(metadata_isolates) # old names
nice_meta_names <- cbind(nice_meta_names, colnames(metaVar_Names)) # cbind do make df with new and old names side by side 
save(nice_meta_names, file = "~/Dropbox/lung/ThedaTill/intermediate_files/Metabolites_nice_meta_name.r")

load("~/Dropbox/lung/ThedaTill/Picrust/20210421_metadeconf_Isolates_metabolites_Final.r")
#write.table(nice_meta_names, "~/Dropbox/lung/ThedaTill/intermediate_files/Metabolites_nice_meta_name.tsv", sep = "\t", quote = F)
nice.meta <- read.table("~/Dropbox/lung/ThedaTill/intermediate_files/Metabolites_nice_meta_name.tsv", sep = "\t",header = T )


nice.feature <- read.table("~/Dropbox/lung/ThedaTill/intermediate_files/Metabolites_long.csv", sep = "," ,header = T )


y <- BuildHeatmap(metadeconf_Isolates_Metabo, metaVariableNames = nice.meta, d_range = "full")
y + theme(legend.position = "none",
             axis.text.y = element_text(size = 14),
             axis.text.x = element_text(size = 14))


save(plot,file ="~/Dropbox/lung/ThedaTill/intermediate_files/20210422_metadeconf_Isolates_Metabo_Names_Final.r")



ggsave(filename = "~/Dropbox/lung/ThedaTill/output/20210421_metadeconf_Isolates_metbolites_nice_Names.svg",
       plot = plot, device = "svg", width = 5, height = 7)


#read Model pathways KO
Iso_Modul <- read.table( "~/Dropbox/lung/ThedaTill/intermediate_files/out_mapped_Modules_20210316.r",
                        sep = "\t", header = TRUE, row.names = 1)

features <- t(Iso_Modul)
rownames(features) <- sub("^X", "", rownames(features))
features <- features[row.names(features)%in%row.names(metadata_isolates), ]
metadata_isolates <- metadata_isolates[row.names(metadata_isolates)%in%row.names(features), ] # 285 sampls in metadata only 280 in iso_KO
dim(features)
ad.index.keep <- which(colSums(features)*100/(sum(features)) > 0.01)
features <- features[, ad.index.keep]
dim(features)
features <- features[order(row.names(features)), ]
metadata_isolates <- metadata_isolates[order(row.names(metadata_isolates)), ]
#normalizing 
metadeconf_Isolates_KO <- MetaDeconfound(featureMat = features,
                                         metaMat = metadata_isolates,
                                         nnodes = 3,
                                         randomVar = list("+ (1|Study_ID) + (1|Site)", c("Study_ID", "Site")),
                                         logfile = "~/Dropbox/lung/ThedaTill/intermediate_files/202104_metadeconf_Isolate_KO.log")


save(metadeconf_Isolates_KO,
     file = "~/Dropbox/lung/ThedaTill/intermediate_files/metadeconf_Isolates_KO.r")

load(file = "~/Dropbox/lung/ThedaTill/intermediate_files/metadeconf_Isolates_KO.r")

# changes module names 
mod_names1 <- read.table("~/Dropbox/Uganda/Uganda_Lotus/PiCrust/binning/script/module.defs", fill = TRUE, sep = "\t", quote = "", check.names = F)
mod_names2 <- read.table("~/Dropbox/Uganda/Uganda_Lotus/PiCrust/binning/script/module.gmm.defs", fill = TRUE, sep = "\t", quote = "", check.names = F)


module_names <- rbind(mod_names1[, c(1,2)], mod_names2[, c(1,2)])
colnames(module_names) <- c("short", "long")
missing_KEGG <- read.table("~/Dropbox/lung/ThedaTill/intermediate_files/KEGG_missing_names.tsv",
                           sep = "\t", check.names = F)
#[1] "M00840" "M00841" "M00842" "M00843" "M00844" "M00845" diese 5 fehlen
colnames(missing_KEGG) <- c("short", "long")
Full_modules <- rbind(module_names, missing_KEGG)
existing_modules <- Full_modules[Full_modules[, 1] %in% rownames(metadeconf_Isolates_KO$status), ]
rownames(metadeconf_Isolates_KO$status)[!(rownames(metadeconf_Isolates_KO$status)%in% module_names[, 1])]

# Isolat_KO and module_names dont have the same length ??
longNames <- merge(x = metadeconf_Isolates_KO$Ps, y = existing_modules, by.x = 0, by.y = 1) # names are now in longNames$V2
# longNames$long <- make.names(longNames$long, unique = TRUE)
### assigning weird long rownmaes leads to errors in buildHeatmap()
rownames(metadeconf_Isolates_KO$Ps) <- longNames$long
rownames(metadeconf_Isolates_KO$Qs) <- longNames$long
rownames(metadeconf_Isolates_KO$Ds) <- longNames$long
rownames(metadeconf_Isolates_KO$status) <- longNames$long
# save(metadeconf_Isolates_KO, file ="~/Dropbox/lung/ThedaTill/Picrust/20210323_metadeconf_Isolates_KO.r" )

plot <- BuildHeatmap(metadeconf_Isolates_KO, intermedData = T)


ggsave(filename = "~/Dropbox/lung/ThedaTill/output/20210406_metadeconf_Isolates_KO.svg",
       plot = plot, device = "svg", width = 18, height = 35)


ggsave(filename = "~/Dropbox/lung/ThedaTill/output/metadeconf_Isolates_KO_First.svg",
       plot = plot, device = "svg", width = 8, height = 10)

# Binn modules names use file from https://www.genome.jp/kegg-bin/show_brite?ko00002.keg 

# KEGG_hirachie <- fromJSON(file = "~/Dropbox/lung/ThedaTill/raw_data/KEGG_Module_hirarchie.json")
# print(KEGG_hirachie)
# KEGG_hirachie <- as.data.frame(KEGG_hirachie)
#####binn the modules use chia-yu's table 

higher_Modules <- read.table("~/Dropbox/lung/ThedaTill/raw_data/Kegg_modules_pathway_myedit.txt", header = T,
                             sep = "\t")
#we want to bin on second_lvl
features <- t(features)
features_red <- features[row.names(features) %in% higher_Modules$Module, ]
features_red <- as.data.frame(features_red)
names_number <- as.data.frame(summary(as.factor(higher_Modules_red$Second_level)))#Summery befehl als df to count the multiple species
features_red <- features_red[order(row.names(features_red)), ]
row.names(higher_Modules_red) <- higher_Modules_red$Module
higher_Modules_red <- higher_Modules_red[order(row.names(higher_Modules_red)), ]
binnedModules <- as.data.frame(matrix(nrow = ncol(features), ncol = 1)) # leerer df
rownames(binnedModules) <- colnames(features) # this must be a datafram empty with the samples names
for (i in unique(higher_Modules_red$Second_level)) {
 print(i)
 # if(i == "Plant_terpenoid_biosynthesis"){
 #   next
 # }
 if (names_number[i, 1] <2) { 
  print(i)
  binnedModules[[i]] <- colSums(features_red[which(higher_Modules_red$Second_level == i), ])
  next
 }
 print(i)
 binnedModules[[i]] <- colSums(features_red[which(higher_Modules_red$Second_level == i), ])
 
}


binnedModules$V1 <- NULL
# Number are different 
binnedModules <- binnedModules[, order(colnames(binnedModules))]
colnames(binnedModules) <- paste0(rownames(names_number), "_", names_number$`summary(as.factor(higher_Modules_red$Second_level))`)
# puts counts on the modules -> suggestins change colnames in names_numbers for the future 
write.table(binnedModules, file = "~/Dropbox/lung/ThedaTill/output/22.7.13new.binned.modules.tsv", sep = "\t", quote= F)

row.names(binnedModules) <- sub("^X", "", rownames(binnedModules))
binnedModules_red <- binnedModules[row.names(metadata_isolates) %in% row.names(metadata_isolates), ]
metadata_isolates_red <- metadata_isolates[row.names(metadata_isolates)%in%row.names(metadata_isolates), ]
ad.index.keep <- which(colSums(binnedModules_red)*100/(sum(binnedModules_red)) > 0.01)
binnedModules_red <- binnedModules_red[, ad.index.keep]
dim(binnedModules_red)
binnedModules_red <- binnedModules_red[order(row.names(binnedModules_red)), ]
metadata_isolates_red <- metadata_isolates_red[order(row.names(metadata_isolates_red)), ]


metadeconf_Isolates_binned_Modules <- metadeconfoundR::MetaDeconfound(featureMat = binnedModules_red,
                                         metaMat = metadata_isolates_red,
                                         nnodes = 3,
                                         randomVar = list("+ (1|Study_ID) + (1|Site)", c("Study_ID", "Site")),
                                         logfile = "~/Dropbox/lung/ThedaTill/intermediate_files/20220713_metadeconf_Isolates_binned_Modules.log")

save(metadeconf_Isolates_binned_Modules, file = "~/Dropbox/lung/ThedaTill/intermediate_files/22.7.13.metadeconf_Isolates_binned_Modules_Final.r")


feature_names <- binnedModules_red
feature_names <- t(feature_names) # run this row.name 5 times to remove all the _
row.names(feature_names) <- sub("_", " ", row.names(feature_names))
row.names(feature_names) <- sub("_", " ", row.names(feature_names))
row.names(feature_names) <- sub("_", " ", row.names(feature_names))
row.names(feature_names) <- sub("_", " ", row.names(feature_names))
row.names(feature_names) <- sub("_", " ", row.names(feature_names))

feature_names <- t(feature_names) # human readle names 
nice_feature_names <- colnames(binnedModules_red) # old names
nice_feature_names <- cbind(nice_feature_names, colnames(feature_names)) # cbind do make df with new and old names side by side 
save(nice_feature_names, file = "~/Dropbox/lung/ThedaTill/intermediate_files/binnes_Isolates_nice_feature_name.r")

metaVar_Names <- metadata_isolates_red
#change colnames 
colnames(metaVar_Names)[which(names(metaVar_Names) == "SampleSite_1")] <- "Throat"
colnames(metaVar_Names)[which(names(metaVar_Names) == "SampleSite_2")] <- "Left Lower Lobe"
colnames(metaVar_Names)[which(names(metaVar_Names) == "SampleSite_3")] <- "Left Upper Lobe"
colnames(metaVar_Names)[which(names(metaVar_Names) == "SampleSite_4")] <- "Lingula"
#now changes the colnames 
colnames(metaVar_Names) <- sub("_", " ",  colnames(metaVar_Names))
colnames(metaVar_Names) <- sub("_", " ",  colnames(metaVar_Names))
colnames(metaVar_Names)[which(names(metaVar_Names) == "Final 16S copies_per_ul")] <- "Final 16S copy/µL"

nice_meta_names <- colnames(metadata_isolates_red) # old names
nice_meta_names <- cbind(nice_meta_names, colnames(metaVar_Names)) # cbind do make df with new and old names side by side 
save(nice_meta_names, file = "~/Dropbox/lung/ThedaTill/intermediate_files/binned_Isolates_nice_meta_name.r")

load("~/Dropbox/lung/ThedaTill/intermediate_files/binnes_Isolates_nice_feature_name.r")
load("~/Dropbox/lung/ThedaTill/intermediate_files/binned_Isolates_nice_meta_name.r")
load("~/Dropbox/lung/ThedaTill/intermediate_files/metadeconf_Isolates_binned_Modules_Final.r")



plot <- metadeconfoundR::BuildHeatmap(metadeconf_Isolates_binned_Modules, featureNames = nice_feature_names, 
                     metaVariableNames  = nice_meta_names, d_col = c("red", "white", "blue"))

save(plot, file = "~/Dropbox/lung/ThedaTill/intermediate_files/20210422_Metadeconfunder_KO_binned_Name_final.r")

     ggsave(filename = "~/Dropbox/lung/ThedaTill/output/20210421_Metadeconf_Isolates_binned_Modules_nice_Names.svg",
       plot = plot, device = "svg", width = 6, height = 8)
# dataframe mit zwei columns -> alte namen daneben neuen namen für 
# beides features and metavariablen 


#Convert metadeconfounder output into long_data_frame
long_out <- reshape2::melt(metadeconf_Isolates_KO$Ps, varnames = c("feature", "metadata"), value.name = "Ps")
long_out$Qs <- reshape2::melt(metadeconf_Isolates_KO$Qs)[, 3]
long_out$Ds <- reshape2::melt(metadeconf_Isolates_KO$Ds)[, 3]
long_out$status <- reshape2::melt(metadeconf_Isolates_KO$status)[, 3]
# long_out_signif <- subset(x = long_out, (status != "NS") & !is.na(status))


# Till 20210422
# plot heatmaps nicely
library(cowplot)
load("~/Dropbox/lung/ThedaTill/intermediate_files/20210422_metadeconf_Isolates_Metabo_Names_Final.r")
metadeconf_Isolates_Metabo <- plot
load("~/Dropbox/lung/ThedaTill/intermediate_files/20210422_Metadeconfounder_fullGenus_Name_Final.r")
Metadeconfounder_fullGenus <- plot
load("~/Dropbox/lung/ThedaTill/intermediate_files/20210422_Metadeconfunder_KO_binned_Name_final.r")
Metadeconfunder_KO_binned <- plot
rm(plot)

metadeconf_Isolates_Metabo <- metadeconf_Isolates_Metabo + theme(plot.title = element_blank(),
                                                                 plot.subtitle = element_blank(),
                                                                 axis.title.x = element_blank(),
                                                                 legend.position = "none") +
        ylab("metabololites")

Metadeconfunder_KO_binned <- Metadeconfunder_KO_binned + theme(plot.title = element_blank(),
                                                               plot.subtitle = element_blank(),
                                                               axis.title.x = element_blank(),
                                                               legend.position = "none")+
        ylab("binned KO modules")
Metadeconfounder_fullGenus <- Metadeconfounder_fullGenus + theme(plot.title = element_blank(),
                                                                 plot.subtitle = element_blank(),
                                                                 legend.position = "right",
                                                                 axis.title.x = element_blank())+
        ylab("bacterial genera")

rightPlot <- plot_grid(metadeconf_Isolates_Metabo,
                       Metadeconfunder_KO_binned, 
                       ncol = 1,
                       rel_heights = c(1.1,1),
                       align = "v")

allHeatmaps <- plot_grid(Metadeconfounder_fullGenus,
                         rightPlot,
                         ncol = 2,
                         rel_widths = c(1.7, 1), labels = c("16S gut microbiome", "isolates"))

ggsave(allHeatmaps, 
       filename = "~/Dropbox/lung/ThedaTill/output/20210422_heatmapsCombined.svg",
       width = 10, 
       height = 10, 
       device = "svg")



# compute prevalence and relative abundance per genus
genus_data <- features
testRelAbu <- colSums(genus_data)
prevalence_genus <- genus_data > 0
prevalence_genus_sum <- colSums(prevalence_genus)
marginalData <- as.data.frame(log(testRelAbu)/max(log(testRelAbu)))
colnames(marginalData) <- "abundance"
marginalData$prevalence <- prevalence_genus_sum/nrow(genus_data)

# create metadeconf heatmap and move its legend to the left
heatmap <-
        BuildHeatmap(
                metadeconf_Celtic_fullGenus,
                metaVariableNames = nice.meta, d_range = "full"
        )
heatmap <-
        heatmap + theme(legend.position = "none",
                        axis.text.y = element_text(face = "italic", size = 14),
                        axis.text.x = element_text(size = 14))
heatmap

# remove features from marginalData, that were removed by BuildHeatmap() internal filtering
dim(marginalData)
length(unique(heatmap$data$feature))
marginalData <- marginalData[rownames(marginalData) %in% heatmap$data$feature, ]

# assign order created by clustering inside BuildHeatmap() to the marginalData features
marginalData$feature <- rownames(marginalData)
marginalData$feature <- factor(as.factor(marginalData$feature),
                               levels = levels(heatmap$data$feature))
# plot a vertical bar plot for abundance
barplot_abun <- ggplot(data = marginalData, aes(x = feature, y = abundance)) +
        geom_bar(stat = "identity", width = 0.8, color ="black", fill = "white") + 
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust = 1, size = 14),
              axis.text.y = element_blank(),
              axis.title.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.title.x = element_text(angle = 90)) +
        scale_y_continuous(limits = c(0.0, 1.0), breaks =c(0.0, 0.5, 1.0),
                           labels = c(0.0, 0.5, 1.0)) +
        coord_flip()
# plot a vertical bar plot for prevalence
barplot_preva <- ggplot(data = marginalData, aes(x = feature, y = prevalence)) +
        geom_bar(stat = "identity", width = 0.8, color ="black", fill = "white") + 
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust = 1, size = 14),
              axis.text.y = element_blank(),
              axis.title.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.title.x = element_text(angle = 90)) +
        scale_y_continuous(limits = c(0.0, 1.0), breaks =c(0.0, 0.5, 1.0),
                           labels = c(0.0, 0.5, 1.0)) +
        coord_flip()
# combine all three plots into one
merged_heatmap <- plot_grid(heatmap, barplot_abun, barplot_preva, nrow = 1, align = "h", rel_widths = c(1, 0.1, 0.1))
merged_heatmap
ggsave(merged_heatmap, 
       filename = "~/Dropbox/lung/ThedaTill/output/Updated_Figure/20210629_heatmapsGenusMarginal.svg",
       width = 8, 
       height = 10, 
       device = "svg")
### Do the same with phylum
load(file = "~/Dropbox/lung/ThedaTill/output/20210421_nf_Celtic_fullPhylum_final.r")
# compute prevalence and relative abundance per genus
genus_data <- features
testRelAbu <- colSums(genus_data)
prevalence_genus <- genus_data > 0
prevalence_genus_sum <- colSums(prevalence_genus)
marginalData <- as.data.frame(
        log(testRelAbu)/max(log(testRelAbu)))
colnames(marginalData) <- "abundance"
marginalData$prevalence <- prevalence_genus_sum/nrow(genus_data)

# create metadeconf heatmap and move its legend to the left
heatmap <- BuildHeatmap(metadeconf_Celtic_fullPhylum, metaVariableNames = nice_meta_names, d_range = "full")
heatmap <- heatmap + theme(legend.position = "left", axis.text.y = element_text(face = "italic"))
heatmap

# remove features from marginalData, that were removed by BuildHeatmap() internal filtering
dim(marginalData)
length(unique(heatmap$data$feature))
marginalData <- marginalData[rownames(marginalData) %in% heatmap$data$feature, ]

# assign order created by clustering inside BuildHeatmap() to the marginalData features
marginalData$feature <- rownames(marginalData)
marginalData$feature <- factor(as.factor(marginalData$feature),
                               levels = levels(heatmap$data$feature))
# plot a vertical bar plot for abundance
barplot_abun <- ggplot(data = marginalData, aes(x = feature, y = abundance)) +
        geom_bar(stat = "identity", width = 0.8, color ="black", fill = "white") + 
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust = 1),
              axis.text.y = element_blank(),
              axis.title.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.title.x = element_text(angle = 90)) +
        scale_y_continuous(limits = c(0.0, 1.0), breaks =c(0.0, 0.5, 1.0),
                           labels = c(0.0, 0.5, 1.0)) +
        coord_flip()
# plot a vertical bar plot for prevalence
barplot_preva <- ggplot(data = marginalData, aes(x = feature, y = prevalence)) +
        geom_bar(stat = "identity", width = 0.8, color ="black", fill = "white") + 
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust = 1),
              axis.text.y = element_blank(),
              axis.title.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.title.x = element_text(angle = 90)) +
        scale_y_continuous(limits = c(0.0, 1.0), breaks =c(0.0, 0.5, 1.0),
                           labels = c(0.0, 0.5, 1.0)) +
        coord_flip()
# combine all three plots into one
merged_heatmap <- plot_grid(heatmap, barplot_abun, barplot_preva, nrow = 1, align = "h", rel_widths = c(1, 0.1, 0.1))
merged_heatmap
ggsave(merged_heatmap, 
       filename = "~/Dropbox/lung/ThedaTill/output/Updated_Figure/20210722_heatmapsPhylumMarginal.svg",
       width = 6, 
       height = 4.5, 
       device = "svg")
