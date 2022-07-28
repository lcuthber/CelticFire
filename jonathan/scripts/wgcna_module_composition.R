# tables showing the composition of the WGCNA modules at each site

library(dplyr)
library(data.table)
library(tidyr)
library(tibble)
library(forcats)
library(stringr)

library(ggplot2)
theme_set(theme_minimal(base_size=16))

# Load required WGCNA results
load("../data/Bus_CF_working_noduplicates.Rdata")
module_assignments <- readRDS("../results/wgcna/module_assignments.rds")
hub_otu_labels <- readRDS("../results/wgcna/hub_otu_labels.rds")
module_eigenvectors <- readRDS("../results/wgcna/module_eigenvectors.rds")
module_membership <- readRDS("../results/wgcna/module_membership.rds")
module_name_tbl <- fread("../results/wgcna/wgcna_module_names.csv") %>%
  mutate(
    module_longname=str_replace(module_longname, "Family_XIII", "Clostridiales spp.") %>%
      str_replace("Prevotella_\\d+", "Prevotella")
  )
taxonomy_df <- phyloseq::tax_table(Bus_CF1)@.Data %>%
  as_tibble() %>%
  dplyr::rename(OTU=OTUID)
X <- readRDS("../data/wgcna/otu_tables_wgcna.rds")


calc_abund_prev <- function(x) {
  outdf <- tibble(
    OTU=colnames(x),
    total_reads=colSums(x),
    perc_site_reads=100.0*colSums(x)/sum(x),
    prevalence=100.0*apply(x, 2, function(xx) sum(xx>0)/length(xx))
  )
  return(outdf)
}

otu_summaries <- lapply(X, calc_abund_prev) %>% 
  bind_rows(.id="site")

#
# create tables describing the OTUs in each module
all_module_assignments <- module_assignments %>%
  bind_rows(.id="site") %>%
  as_tibble()

all_module_membership <- module_membership %>% 
  bind_rows(.id="site") %>%
  mutate(module=str_remove(module, "ME")) %>%
  select(-pvalue)

table_data <- all_module_assignments %>%
  left_join(otu_summaries, by=c("site", "OTU")) %>%
  left_join(taxonomy_df, by="OTU") %>%
  left_join(all_module_membership, by=c("site", "OTU", "module")) %>%
  group_by(site, module) %>%
  filter(perc_site_reads > 0.001) %>%
  slice_max(correlation, n=20) %>%
  ungroup() %>%
  mutate(
    total_reads=scales::comma(total_reads, accuracy=1),
    label=sprintf("%s (%s %%)", total_reads, format(round(perc_site_reads, digits=3), nsmall=3)),
    label=str_replace(label, "\\( ", "\\("),
    prevalence=sprintf("%s %%", round(prevalence, digits=1)),
    correlation=round(correlation, digits=2)
  ) %>% 
  select(-c(Seq_name, Species, Kingdom)) %>%
  relocate(site, module, OTU, Phylum, Class, Order, Family, Genus, correlation) %>%
  select(-c(Class, Order, Family)) %>%
  rename(Reads=label, Prevalence=prevalence, Correlation=correlation)
table_data

table_data %>%
  select(-c(Reads)) %>%
  fwrite("../tmp/wgcna_module_tables_all_sites_modules.tsv", sep="\t")

output_filename <- "../tmp/wgcna_module_tables.txt"
new_file <- TRUE

for(site_name in unique(table_data$site)) {
  write(paste0(site_name, "\n"), file=output_filename, append=!new_file)
  new_file <- FALSE

  tmp_tbl <- table_data %>%
    filter(site==site_name) %>%
    select(-site)
  
  for(module_name in unique(tmp_tbl$module)) {
    write(paste0(module_name, "\n"), file=output_filename, append=TRUE)
    tmp_tbl %>%
      filter(module==module_name) %>%
      select(-module) %>%
      fwrite(output_filename, col.names=TRUE, append=TRUE, sep="\t")
    write("\n", file=output_filename, append=!new_file)
  }
  
  write("\n", file=output_filename, append=!new_file)
}
