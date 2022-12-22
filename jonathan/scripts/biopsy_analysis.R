library(dplyr)
library(tidyr)
library(readxl)
library(tibble)
library(stringr)
library(data.table)
library(phyloseq)
library(rstatix)
library(broom)

library(ggplot2)
theme_set(theme_classic(base_size=14))


kruskal_with_effsize <- function(...) {
  tryCatch({
    kruskal_pvals <- kruskal_test(...) %>% select(-method)
    kruskal_effects <- kruskal_effsize(...) %>% select(-method)
    return(kruskal_pvals %>% inner_join(kruskal_effects, by=c(".y.", "n")))
  },
  error=function(e) {
    return(tibble())
  })
  
}

# biopsy scores
biopsy_scores <- read_xlsx("../data/biopsies/H & E Celtic Fire Scoring Data v1 02_12_2022.xlsx")
biopsy_scores_valid <- biopsy_scores %>% 
  filter(`Problem with Scan/Image`==0, `Sample Adequacy`==1)

biopsy_columns <- c("Patient ID", "Bronchial Epithelium", "Denudation", "Epithelial hyperplasia",
                    "Goblet cell hyperplasia", "Epithelial metaplasia", "Basement membrane",
                    "Inflammation", "Macrophages", #  "Type of abnormality",
                    "Giant cells/Granuloma",	"Lymphocytes",	"Plasma cells",	"Eosinophils", "Neutrophils",	"Bacteria"
)

biopsy_scores_clean <- biopsy_scores_valid %>%
  select(
    `Disease (Case) or Control`,
    all_of(
      setNames(biopsy_columns, str_replace_all(biopsy_columns, " |/", "_"))
    )
  ) %>%
  mutate(
    across(
      all_of(c("Inflammation", "Macrophages", "Giant_cells_Granuloma", "Lymphocytes", "Plasma_cells", 
               "Eosinophils", "Neutrophils", "Bacteria")),
      function(x) recode(as.character(x), "0"="absent", "1"="sparse", "2"="abundant")
    ),
    across(
      all_of(c("Denudation", "Epithelial_hyperplasia",
               "Goblet_cell_hyperplasia", "Epithelial_metaplasia", "Basement_membrane")),
      function(x) recode(as.character(x), "0"="No", "1"="Yes")
    ),
    Bronchial_Epithelium=recode(as.character(Bronchial_Epithelium), "0"="Normal", "1"="Abnormal")
  ) %>%
  rename(Disease=`Disease (Case) or Control`)

#
# compare biopsy labels to disease status
biopsy_scores_clean %>% 
  pivot_longer(-c(Patient_ID, Disease)) %>%
  group_by(Disease, name, value) %>%
  tally() %>%
  ggplot(aes(x=value, y=n, fill=Disease)) +
  geom_col(position=position_dodge(0.9, "single")) +
  facet_wrap(~name, scales="free_x") +
  ylab("Number of individuals") +
  xlab("Biopsy label")

safe_fisher_test <- function(...) {
  tryCatch({
    fisher.test(...) %>% tidy()
  },
  error=function(e) {
    return(tibble())
  })
}

biopsy_scores_clean %>% 
  pivot_longer(-c(Patient_ID, Disease))  %>% 
  group_by(name) %>%
  do(safe_fisher_test(x=as.factor(.$value), y=as.factor(.$Disease))) %>%
  ungroup() %>%
  filter(!is.na(estimate)) %>%
  mutate(p_adj=p.adjust(p.value, method="fdr")) %>%
  arrange(p_adj)

#
# correlation network anaylsis for histology lbaelss
safe_cramersv <- function(...) {
  tryCatch({
    return(cramer_v(...))
  },
  error=function(e) {
    return(NA)
  }
  )
}

hist_assocs <- biopsy_scores_clean %>%
  pivot_longer(-c(Patient_ID, Disease)) %>%
  group_by(name) %>%
  mutate(n_levels=length(unique(value))) %>%
  filter(n_levels>1) %>%
  tidybayes::gather_pairs(name, value, triangle="both") %>%
  summarise(v=safe_cramersv(x=.x, y=.y)) %>%
  ungroup() 

hist_assocs %>%
  ggplot(aes(x=.col, y=.row)) +
  geom_tile(aes(fill=v)) +
  scale_fill_gradient(low="white", high="red") +
  theme(axis.text.x=element_text(angle=60, hjust=1, vjust=1))

ComplexHeatmap::Heatmap(
  hist_assocs %>%
    pivot_wider(names_from=.col, values_from=v) %>%
    column_to_rownames(".row") %>%
    as.matrix(),
  col=circlize::colorRamp2(c(0,1), c("white", "red")),
  name="Cramers_v",
  row_dend_width=unit(30, "mm"),
  column_dend_height=unit(30, "mm"),
  heatmap_width=unit(200, "mm"),
  heatmap_height=unit(200, "mm")
)
  

#
# compare biopsy labels to WGCNA eigenvectors

# WGCNA module eigenvectors
module_eigenvectors <- readRDS("../results/wgcna/module_eigenvectors.rds")

# need to map samples to individuals
map <- fread("../data/metadata/formatted_metadata.csv") %>%
  as_tibble()

sample_id_map <- map %>%
  filter(Study!="Busselton") %>%
  select(SampleID, Study_ID) %>%
  distinct() %>%
  mutate(Study_ID=as.character(Study_ID))

samples_with_indiv_IDs <- lapply(module_eigenvectors, function(x) tibble(SampleID=rownames(x))) %>%
  bind_rows(.id="site") %>%
  filter(grepl("cf_", site)) %>%
  left_join(
    sample_id_map, by="SampleID"
  )

matched_IDs <- biopsy_scores_clean %>% 
  select(Patient_ID) %>%
  mutate(Study_ID=str_remove(Patient_ID, "^0+")) %>%
  inner_join(samples_with_indiv_IDs, by="Study_ID")

tmpdata <- lapply(
  module_eigenvectors[ grepl("cf_", names(module_eigenvectors)) ],
 function(x) x %>%
   rownames_to_column("SampleID") %>%
   pivot_longer(-SampleID, names_to="module", values_to="mod_eigenvec_val")
) %>%
  bind_rows(.id="site") %>%
  left_join(matched_IDs) %>%
  left_join(biopsy_scores_clean %>% 
              pivot_longer(-c(Patient_ID, Disease), names_to="histology_label", values_to="histology_val")) 

kruskal_test_res <- tmpdata %>%
  group_by(site, module, histology_label) %>%
  do(kruskal_with_effsize(data=., mod_eigenvec_val~histology_val)) %>%
  ungroup()

kruskal_test_res %>%
  filter(magnitude!="small") %>%
  left_join(tmpdata) %>%
  unite(module, c(module, site)) %>%
  ggplot(aes(x=histology_val, y=mod_eigenvec_val)) +
  geom_boxplot() +
  geom_jitter() +
  facet_wrap(~histology_label+module, scales="free") +
  geom_hline(yintercept=0, colour="red", linetype="dotted", size=1)

#
# differential abundance analysis



