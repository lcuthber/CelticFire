##Load packages 
library(phyloseq)
library(microbiome)
library(tidyverse)
library(data.table)
require(ggpubr)
require(RColorBrewer)
require(rstatix)
library(cowplot)
library(gridExtra)
library(grid)
library(doParallel)  
library(ggsci)

#now we are looking at disease status
metadata <- read.table(
 "~/Dropbox/lung/clean_metadata_20201021.csv",
 sep = "\t",
 header = TRUE,
 row.names = 1
)
metadata <- metadata[order(rownames(metadata)),]
## get Pulmotype results
all_pulmotype <-
 read.table(
  "~/Dropbox/lung/ThedaTill/20201117_Enterotyping/20201217_pulmotype_result_assignment_RTK.txt",
  sep = "\t",
  header = T
 )
row.names(all_pulmotype) <-
 str_replace(row.names(all_pulmotype), "XX", "X")
metadata_chiq <-
 metadata[rownames(metadata) %in% rownames(all_pulmotype),]

### use only the meta
load(file = "~/Dropbox/lung/ThedaTill/intermediate_files/Chiqs_CF_meta.r")

meta_chisq <- meta
meta_chisq <-
 fastDummies::dummy_cols(meta_chisq, select_columns = "SampleSite")
colnames(meta_chisq)[which(names(meta_chisq) == "SampleSite_1")] <-
 "Throat"
colnames(meta_chisq)[which(names(meta_chisq) == "SampleSite_2")] <-
 "Left Lower Lobe"
colnames(meta_chisq)[which(names(meta_chisq) == "SampleSite_3")] <-
 "Left Upper Lobe"
colnames(meta_chisq)[which(names(meta_chisq) == "SampleSite_4")] <-
 "Lingula"

p <- chisq.test (table (meta_chisq$Asthma, meta_chisq$Pulmotype))
## X-squared = 3.8975, df = 1, p-value = 0.04836
print(chisq.test(table(meta_chisq$Throat, meta_chisq$Pulmotype))) # this is müll
## X-squared = 29.533, df = 1, p-value = 5.499e-08
x <- print(chisq.test(table(
 meta_chisq$`Left Lower Lobe`, meta_chisq$Pulmotype
)))
## X-squared = 4.7494, df = 1, p-value = 0.02931
y <- print(chisq.test(table(
 meta_chisq$`Left Upper Lobe`, meta_chisq$Pulmotype
)))
## X-squared = 15.165, df = 1, p-value = 9.849e-05
z <- print(chisq.test(table(
 meta_chisq$Lingula, meta_chisq$Pulmotype
)))
#### FDR corrected p values :) 
p.adjust(c(x$p.value, y$p.value, z$p.value), method = "fdr")
# [1] 0.0439636673 0.0002954779 1.0000000000



### X-squared = 6.4064e-32, df = 1, p-value = 1
# ### This is not so nice. Instead compare all together
# meta_chisq_u <-
#   meta_chisq %>% unite("Swabsite", Throat:Lingula, remove = F)
# 
# meta_chisq_u$Swabsite[meta_chisq_u$Swabsite == "0_0_0_1"] <- "Lingula"
# meta_chisq_u$Swabsite[meta_chisq_u$Swabsite == "0_1_0_0"] <- "Left Lower Lobe"
# meta_chisq_u$Swabsite[meta_chisq_u$Swabsite == "0_0_1_0"] <- "Left Upper Lobe"
# meta_chisq_u$Swabsite[meta_chisq_u$Swabsite == "1_0_0_0"] <- "OTS"
# 
# meta_chisq_u <- subset(meta_chisq_u, Swabsite != "OTS")
# 
# print(chisq.test(table(meta_chisq_u$Swabsite, meta_chisq_u$Pulmotype))) # this is müll



### now look at the busselton vs. celtic fire only ots
metadata_chiq_x  <- metadata_chiq %>%
 filter(throatswab == 1) ## this is tidyverse!! select only Throatswabs = 1
all_pulmotype <-
 all_pulmotype[rownames(all_pulmotype) %in% rownames(metadata_chiq_x),]
order(rownames(metadata_chiq_x))
order(rownames(all_pulmotype))
metadata_chiq_x$Pulmotype <- all_pulmotype$DMM_k.2  

print(chisq.test(table(metadata_chiq_x$Study, metadata_chiq_x$Pulmotype)))
## X-squared = 2.6415e-29, df = 1, p-value = 1
load(file = "~/Dropbox/lung/ThedaTill/intermediate_files/congurance_metadata_chiq.r")

con_chiq <- con_3
con_chiq <-
 fastDummies::dummy_cols(con_chiq, select_columns = "SameDonor")
con_chiq <-
 fastDummies::dummy_cols(con_chiq, select_columns = "SamePulmotype")

print(chisq.test(table(con_chiq$SameDonor, con_chiq$SamePulmotype)))

# X-squared = 362.04, df = 1, p-value < 2.2e-16
 print(chisq.test(table(con_chiq$SameSampleSite, con_chiq$SamePulmotype)))
# X-squared = 1299.7, df = 1, p-value <
 # 2.2e-16

## now compare the OTS only in the two cohorts 
load(file = "~/Dropbox/lung/ThedaTill/raw_data/metadata_complete_Pulmo.r")
metadata_OTS <- subset(metadata_reduced, SwabSite == "OTS")
print(chisq.test(table(metadata_OTS$Study, metadata_OTS$Pulmotype))) # this is müll
# X-squared = 1.0037e-29, df = 1, p-value
# = 1
