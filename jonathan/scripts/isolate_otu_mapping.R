library(WGCNA)
library(phyloseq)
library(dplyr)
library(data.table)
library(tidyr)
library(tibble)
library(stringr)

library(ggplot2)

source("../scripts/kegg_utils.R")

#
# Load OTU tables
load("../data/Bus_CF_working_noduplicates.Rdata")

X <- list()
X[["cf_lll"]] <- t(as.matrix(otu_table(subset_samples(Bus_CF1, SwabSite=="Brushing LLL 16S"))))
X[["cf_lul"]] <- t(as.matrix(otu_table(subset_samples(Bus_CF1, SwabSite=="Brushing LUL 16S"))))
X[["cf_ots"]] <- t(as.matrix(otu_table(subset_samples(Bus_CF1, SwabSite=="OTS"))))
X[["bus_ots"]] <- t(as.matrix(otu_table(subset_samples(Bus_CF1, SwabSite=="Throat"))))
lapply(X, dim)

# load taxonomy of OTUs
taxonomy_df <- phyloseq::tax_table(Bus_CF1)@.Data %>%
  as.data.frame() %>%
  rename(OTU=OTUID)

#
# Load OTU -> isolate mapping
isolate_otu_map <- fread("../data/isolate_names/filtered_OTU_isolate.txt") %>%
  rename(label=lable)

# change to nsp naming
isolate_name_label_map <- fread("../data/isolate_names/new_label_isolates_nsp.csv") %>% # fread("../Data/new_lable_isolates.tsv") %>%
  rename(old_name=Isolate, new_name=label) %>%
  select(-Genus)

# this files contains a Genus NA for sequence HK556267.1.1522, OTU Haemophilus_15385
isolate_otu_map <- isolate_otu_map %>%
  select(OTU, isolate, percentage_identity) %>%
  rename(old_name=isolate) %>%
  left_join(isolate_name_label_map, by="old_name") %>%
  rename(isolate=new_name, Seq_name=OTU) %>%
  mutate(isolate=coalesce(isolate, old_name)) %>%
  select(-old_name)

cat(sprintf("Loaded %d OTUs and %d isolates in mapping file\n",
            length(unique(isolate_otu_map$Seq_name)),
            length(unique(isolate_otu_map$isolate))))

isolate_otu_map_w_pc_id <- taxonomy_df %>% 
  select(Seq_name, OTU) %>%
  left_join(
    isolate_otu_map, by="Seq_name"
  ) %>%
  select(-Seq_name) %>%
  replace_na(list(isolate="no_isolate")) %>%
  as_tibble()

cat(sprintf("OTU-isolate mapping contains %d OTUs and %d isolates\n",
            length(unique(isolate_otu_map_w_pc_id$OTU)),
            length(unique(isolate_otu_map_w_pc_id$isolate))))

isolate_otu_map_w_pc_id %>%
  fwrite("../data/isolate_names/otu_isolate_map.csv")
