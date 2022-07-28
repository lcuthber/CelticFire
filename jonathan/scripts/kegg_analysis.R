# This script loads the N x K matrix of KO occurrances (N isolates, K KOs)
# and runs hierarchichal clustering on using the Manhattan distance between
# the isolates (i.e. phylogenetic distance). 15 clusters are identified
# by the dynamic tree cut algorithm and then KOs that are specific/close to specfic 
# for a given cluster are identified using

library(data.table)
library(dplyr)
library(tidyr)
library(tibble)
library(dendextend)
library(matrixStats)
library(dynamicTreeCut)
library(stringr)

library(ggplot2)
theme_set(theme_classic(base_size=20))
library(ggpubr)
library(ggforce)
library(drlib)
library(lemon)

source("../scripts/kegg_utils.R")
source("../scripts/ko_lookup_utils.R")

options(stringsAsFactors=FALSE)

##############################################################################
# Load isolate x KO matrix
##############################################################################

binarise_column <- function(x) { as.integer(x>0) }

# If TRUE, load emapper file without duplicated KOs from data/kegg/formatted/duplicated_KOs.csv
# otherwise check which KOs are duplicated (takes a few minutes)
LOAD_EXISTING_EMAPPER <- TRUE

if(!LOAD_EXISTING_EMAPPER) {
  
  # load KEGG matrix then encode as binary
  emapper <- fread("../data/kegg/formatted/isolate_ko_matrix.csv") %>%
    mutate_if(is.integer, binarise_column)
  
  cat("KEGG emapper contains", emapper %>% rownames() %>% unique() %>% length(), "unique isolates\n")
  cat("KEGG emapper contains", emapper %>% ncol() -1, "unique KOs\n")
  
  # change to new isolate names (use nsp)
  isolate_name_map_nsp <- fread("../data/isolate_names/new_label_isolates_nsp.csv") %>%
    rename(old_name=Isolate, new_name=label) %>%
    select(-Genus)
  
  emapper <- emapper %>%
    left_join(isolate_name_map_nsp, by=c(isolate="old_name")) %>%
    mutate(new_name=coalesce(new_name, isolate)) %>%
    select(-isolate) %>%
    rename(isolate=new_name)
  
  emapper <- emapper %>%
  #  filter(!grepl("Haemophilus_influenzae", isolate)) %>% # H. influenzae analysed separately
    column_to_rownames("isolate")
  
  ##############################################################################
  ##############################################################################
  
  ##############################################################################
  # Filter zero-variance KOs
  ##############################################################################
  
  zero_variance_KOs <- colVars(as.matrix(emapper))>0
  cat("Retaining", sum(zero_variance_KOs), "zero-variance KOS\n")
  
  emapper <- emapper[,zero_variance_KOs]
  emapper %>% dim()
  cat(ncol(emapper), "KOs remaning\n")
  
  ##############################################################################
  # remove duplicated KOs
  ##############################################################################
  
  duplicate_ko_res <- group_duplicate_columns(emapper, n_cores=22) # takes a few minutes

  # keep unique KOs for subsequent analysis
  duplicate_ko_res$new_df %>%
    rownames_to_column("isolate") %>%
    fwrite("../data/kegg/formatted/binary_emapper_matrix_no_duplicates.csv")
  
  duplicate_ko_res$tbl %>% 
    group_by(included_column) %>%
    summarise(identical_to=paste0(excluded_column, collapse="|")) %>%
    rename(representative_ko=included_column) %>%
    fwrite("../data/kegg/formatted/duplicated_KOs.csv")
  
  emapper <- duplicate_ko_res$new_df
} else {
  
  # load pre-saved results, 126 x 2,964
  emapper <- fread("../data/kegg/formatted/binary_emapper_matrix_no_duplicates.csv") %>%
    column_to_rownames("isolate")
  dim(emapper)
  
  duplicated_KOs <- fread("../data/kegg/formatted/duplicated_KOs.csv") %>%
    mutate(identical_to=str_split(identical_to, "\\|")) %>%
    unnest(identical_to)
  
}


##############################################################################
##############################################################################

###############################################################################
# Phylogenetic clustering of isolates
# (cluster using KO distances then use Dynamic Tree Cut)
##############################################################################

# complete linkage and Manhattan distance
isolate_hclust <- hclust(dist(emapper, method="manhattan"), method="complete")
dend <- isolate_hclust %>%
  as.dendrogram()
par(mar=c(0,0,0,20)+0, cex=0.65)
dend %>% plot(horiz=TRUE)

cluster_labels <- cutreeDynamic(
  isolate_hclust,
  minClusterSize=1,
  distM=as.matrix(dist(emapper, method="manhattan"))
)
names(cluster_labels) <- rownames(emapper)

# give Strep III its own cluster
strep3_isolates <- isolate_hclust %>% labels() %>% head(9)
cluster_labels[ names(cluster_labels) %in% strep3_isolates ] <- 16 

# save isolate phylogenetic clustering and cluster labels
saveRDS(isolate_hclust, "../results/kegg/isolate_phylogenetic_hclust.rds")
saveRDS(cluster_labels, "../results/kegg/isolate_cluster_labels.rds")

cluster_names <- tibble(
  cluster_longname=c(
    "Strep.I", "Strep.II", "Veillonella", "Rothia", 
    "Gemella", "Prevotella", "Micrococcus", "Neisseria", "Pauljensenia", "Staphy.+Nialla",
    "Haemophilus", "Granulicatella", "Fusobacterium+Leptotrichia",  "Cutibacterium", "Actinomyces",
    "Strep.III"
  )
) %>%
  rownames_to_column("cluster") %>%
  mutate(cluster=as.integer(cluster))

cluster_labels %>%
  as.data.frame() %>%
  rownames_to_column("isolate") %>%
  rename(cluster=".") %>%
  left_join(cluster_names) %>%
  fwrite("../results/kegg/isolate_cluster_labels.csv")

# tibble with list of isolates in each cluster
tibble(
  isolate=names(cluster_labels),
  cluster=cluster_labels
) %>%
  group_by(cluster) %>%
  summarise(
    n_isolates=length(unique(isolate)),
    all_isolates=paste0(unique(isolate), collapse=",")
  )

#
# Plot the isolate dendrogram - not paper version
cmap <- gg_color_hue(length(unique(cluster_labels))) # randomcoloR::distinctColorPalette(k=length(unique(cluster_labels2)))
cmap <- c("grey", cmap)
label_colours  <- cmap[ cluster_labels+1 ]

label_colours <- tibble(
  label=names(cluster_labels),
  cluster=cluster_labels,
  colour=label_colours
) %>%
  arrange(match(row_number(), order.dendrogram(dend)))

legend_keys <- label_colours %>%
  select(-label) %>%
  distinct() %>%
  left_join(cluster_names, by="cluster") %>%
  replace_na(list(cluster_longname="<NA>"))

label_renamer <- setNames(
  paste0(label_colours$label, "__", label_colours$cluster),
  label_colours$label
)

pdf("../plots/kegg/isolate_clusters_dynamic.pdf", width=14, height=9)
  par(mar=c(0,0,0,20)+0, cex=0.5)
  dend %>%
    set("labels", label_renamer) %>%
    color_labels(col=label_colours$colour) %>%
    color_branches(clusters=label_colours$cluster, groupLabels=paste0("cluster_", legend_keys$cluster)) %>%
    plot(horiz=TRUE)
  legend(
    "bottomleft",
    legend=legend_keys$cluster_longname,
    col=legend_keys$colour,
    pch=19,
    cex=2.3
  )
dev.off()

##############################################################################
##############################################################################

##############################################################################
# Score KOs for each cluster using Odds Ratio as summary statistic
##############################################################################

# cluster labels
y <- setNames(cluster_labels, names(cluster_labels))
cluster_tbl <- tibble(
  isolate=names(y),
  cluster=y
) %>%
  left_join(cluster_names,
            by="cluster")

cluster_sizes <- cluster_tbl %>%
  group_by(cluster, cluster_longname) %>%
  tally(name="cluster_size")

# number of isolates each KO is mapped to
n_isolates_per_KO <- tibble(
  KO=colnames(emapper),
  n_isolates=colSums(emapper)
)

y <- y[ rownames(emapper) ]
stopifnot(identical(names(y), rownames(emapper)))

# calculate Odds Ratios for each cluster
# one-vs-all approach
or_result <- list()
for(this_cluster_name in unique(y)) {
  cat(this_cluster_name, "\n")
  
  yy <- as.integer(y==this_cluster_name)
  
  # fit one vs all logistic regression model
  odds_ratios <- pbmcapply::pbmclapply(
    setNames(colnames(emapper), colnames(emapper)),
    function(ko_name) calculate_or(emapper[[ ko_name ]], yy),
    mc.cores=22
  )
  
  or_result[[ paste0("cluster_", this_cluster_name) ]] <- odds_ratios %>%
    rbindlist(idcol="KO") 
}

saveRDS(or_result, "../results/kegg/cluster_ko_odds_ratios.rds")

# combine KO scores from all clusters
odds_ratios_all_clusters <- or_result %>%
  bind_rows(.id="cluster") %>%
  mutate(cluster=as.integer(str_remove(cluster, "cluster_"))) %>%
  left_join(cluster_names, by="cluster")
stopifnot(
  length(unique(odds_ratios_all_clusters$KO))==ncol(emapper)
)

odds_ratios_all_clusters %>%
  ggplot(aes(x=Odds_ratio)) +
  geom_histogram() +
  scale_x_log10() +
  lemon::facet_rep_wrap(~cluster_longname) +
  geom_vline(xintercept=1, colour="red", size=1) +
  geom_vline(xintercept=10, colour="darkgreen", size=1) +
  ylab("Number of KOs") +
  xlab("Odds ratio")
ggsave("../plots/kegg/all_ko_ORs.pdf")

#
# sanity check plots
top_hits <- odds_ratios_all_clusters %>%
  slice_max(abs(log10(Odds_ratio)), n=20)
emapper[, top_hits$KO ] %>%
  rownames_to_column("isolate") %>%
  pivot_longer(-isolate, names_to="KO") %>%
  left_join(top_hits %>% select(KO, Odds_ratio, cluster),
            by="KO") %>%
  mutate(KO_label=sprintf(
    "%s, %s, %s",
    KO,
    ifelse(Odds_ratio>1, "positive", "negative"),
    cluster
  )) %>%
  select(-c(cluster)) %>%
  left_join(cluster_tbl %>%
              mutate(cluster=factor(
                cluster,
                levels=unique(cluster),
                labels=str_replace(unique(cluster), "cluster_", ""))),
            by="isolate") %>%
  ggplot(aes(y=isolate, x=KO_label, fill=as.factor(value))) +
  geom_tile() +
  facet_grid(rows=vars(cluster), scales="free_y", space="free",
             switch="y") +
  theme(axis.text.x=element_text(size=12, angle=90, hjust=1, vjust=0.5),
        axis.text.y=element_text(size=4),
        legend.position="none",
        panel.spacing.y=unit(0, "lines"),
        strip.placement="outside",
        strip.background.y=element_blank(),
        panel.background=element_rect(fill = NA, color = "black")) +
  scale_fill_manual(values=c("white", "red"))

# KO information
# this file is produced by scripts/parse_kegg_info.R
# ko_lookup_slim <- fread("../data/kegg/formatted/formatted_ko_lookup.tsv", sep="\t")
ko_lookup_raw <- readRDS("../data/kegg/ko_lookup_raw.rds")
ko_lookup_slim <- lapply(ko_lookup_raw, format_query) %>%
  bind_rows()

# add duplicated KOs
cluster_ko_scores <- odds_ratios_all_clusters %>%
  mutate(score=log10(Odds_ratio)) %>%
  select(cluster, KO, score)

column_order <- cluster_names %>%
  arrange(cluster_longname) %>%
  mutate(cluster=sprintf("%s__(%s)", cluster_longname, cluster)) %>%
  pull(cluster)

cluster_scores_final_table_data <- cluster_ko_scores %>%
  left_join(cluster_names, by="cluster") %>%
  mutate(cluster=sprintf("%s__(%s)", cluster_longname, cluster)) %>%
  select(-cluster_longname) %>%
  left_join(n_isolates_per_KO, by="KO") %>%
  pivot_wider(names_from=cluster, values_from=score) %>%
  left_join(duplicated_KOs %>%
              rename(KO=representative_ko, excluded_KO=identical_to) %>%
              bind_rows(
                tibble(KO=unique(.$KO), excluded_KO=unique(.$KO))
              ),
            by="KO") %>%
  rename(base_KO=KO, KO=excluded_KO) %>%
  relocate(KO, base_KO) %>%
  mutate(KO=coalesce(KO, base_KO),
         perc_isolates=100*n_isolates/126) %>%
  # left_join(ko_lookup_slim, by="KO") %>%
  relocate(KO, base_KO, n_isolates, perc_isolates) %>%
  arrange(base_KO, KO) %>%
  mutate_if(is.numeric, function(x) round(x, digits=3)) %>%
  left_join(ko_lookup_slim %>% rename(KO=entry),
            by="KO") %>%
  relocate(KO, base_KO, name, symbol, n_isolates, perc_isolates, module, pathways) %>%
  relocate(all_of(column_order), .after=last_col())

stopifnot(length(unique(cluster_scores_final_table_data$KO))==5277)
stopifnot(length(unique(cluster_scores_final_table_data$base_KO))==2964)
stopifnot(sum(is.na(cluster_scores_final_table_data %>% select(-c(symbol, module, pathways, name))))==0)

cluster_scores_final_table_data %>%
  fwrite("../results/kegg/final_KO_cluster_scores.tsv", sep="\t")

# Supplemantary Table 1 - add ubiquitous KOs
emapper_complete <- fread("../data/kegg/formatted/isolate_ko_matrix.csv")
all_KOs <- emapper_complete %>% colnames()
all_KOs <- all_KOs[ all_KOs!="isolate" ]
stopifnot(length(all_KOs)==5531)

emapper_complete %>% 
  column_to_rownames("isolate") %>% 
  t() %>%
  as.data.frame() %>%
  rownames_to_column("KO")

supp_tab_1_data <- ko_lookup_slim %>% 
  as_tibble() %>%
  rename(KO=entry) %>%
  bind_rows(
    tibble(
      KO=setdiff(all_KOs, .$KO),
      symbol="not_in_database",
      name="not_in_database",
      module="not_in_database",
      pathways="not_in_database"
    )
  ) %>%
  left_join(
    emapper_complete[,2:ncol(emapper_complete)] %>%
      do(tibble(
        KO=colnames(.),
        n_isolates=apply(., 2, function(x) sum(x>0))
        )) %>%
      mutate(perc_isolates=100*n_isolates/nrow(.)),
    by="KO"
    ) %>%
  relocate(KO, symbol, n_isolates, perc_isolates) %>%
  mutate(perc_isolates=format(perc_isolates, digits=1)) %>%
  left_join(emapper_complete %>% 
              column_to_rownames("isolate") %>% 
              t() %>%
              as.data.frame() %>%
              rownames_to_column("KO") %>%
              mutate_if(is.integer, binarise_column),
            b="KO")

stopifnot(length(unique(supp_tab_1_data$KO))==5531)
supp_tab_1_data %>%
  write_xlsx("../tables/kegg/supplementary_table_1.xlsx")

# write as excel file - Supplementary Table 2
library(writexl)

excel_sheet_spec <- list(
  c("2a All", NA, NA),		
  c("2b Uncharacterised",	"uncharacterized protein", "name"),
  c("2c Biofilm",	"Biofilm",	"pathways"),
  c("2d Antibiotic",	"antibiotic|mycin",	"pathways"),
  c("2e Toxins",	"toxin",	"name"),
  c("2f NO",	"nitric",	"name"),
  c("2g Iron", 	"iron",	"name"),
  c("2h Heme",	"heme",	"name"),
  c("2i Sphingolipid and ceramide",	"Sphingolipid|sphingolipid", "pathways"),
  c("2j Immune",	"immun",	"name"),
  c("2k CRISPR",	"CRISPR",	"name"),
  c("2l Sporulation",	"sporu",	"name")
) %>%
  do.call(rbind, .) %>%
  as_tibble() %>%
  magrittr::set_colnames(c("Subtable", "Search", "Field"))

excel_sheet_spec %>%
  fwrite("../tmp/excel_search_terms.csv")

xlsx_sheets <- list()

for(i in 1:nrow(excel_sheet_spec)) {
  cat("Sheet", i, "of", nrow(excel_sheet_spec), "\n")
  
  if(is.na(excel_sheet_spec$Field[[i]])) {
    xlsx_sheets[[ excel_sheet_spec$Subtable[[ i ]] ]] <- cluster_scores_final_table_data
  } else {
    xlsx_sheets[[ excel_sheet_spec$Subtable[[ i ]] ]] <- cluster_scores_final_table_data %>%
      filter(grepl(excel_sheet_spec$Search[[i]], !!sym(excel_sheet_spec$Field[[i]]), ignore.case=TRUE))
  }
}

write_xlsx(xlsx_sheets,
           path="../tables/kegg/supplementary_table_2.xlsx")

##############################################################################
##############################################################################

##############################################################################
# Plot the hierarchy of the KOs most assigned with each isolate cluster
##############################################################################

# largest positive KO associations for each cluster
top_hits <- odds_ratios_all_clusters %>%
  # group_by(cluster, cluster_longname) %>%
  # slice_max((Odds_ratio), n=500) 
  filter(Odds_ratio>10)
cat(sprintf("Top KO hits include %d unique KOs\n", length(unique(top_hits$KO))))

top_hit_ko_lookup_raw <- lookup_kos(unique(top_hits$KO), 20)
parsed_top_hit_info <- pbmclapply(top_hit_ko_lookup_raw, parse_brite, mc.cores=20) %>%
  bind_rows()

plot_data <- top_hits %>%
  left_join(parsed_top_hit_info, by="KO") %>%
  group_by(cluster, cluster_longname, term1, term2, term3) %>%
  tally() %>%
  ungroup() %>%
  arrange(term1, term2) %>%
  mutate(
    term1=str_remove(term1, "\\d+") %>%
      str_remove("\\d+$") %>%
      trimws(),
    term2=factor(term2, levels=unique(term2)),
    term3=factor(term3, levels=unique(term3))
  )  

# plot_data %>%
#   group_by(cluster, cluster_longname) %>%
#   do(
#     plot=ggplot(data=., aes(y=term2, x=n, fill=term1)) +
#       geom_col() +
#       facet_grid(cols=vars(cluster_longname), rows=vars(term1))
#   )
  
plot_data %>%
  filter(
    !is.na(term1),
    !term1 %in% c("Not Included in Pathway or Brite", "Brite Hierarchies")
  ) %>%
  right_join(
    cluster_sizes %>% filter(cluster_size>=5),
    by=c("cluster", "cluster_longname")
  ) %>%
  mutate(
    term1=str_replace(term1, " ", "\n")
  ) %>%
  ggplot(aes(y=term2, x=n, fill=term1)) +
  geom_col() +
  facet_grid(
    cols=vars(cluster_longname), rows=vars(term1),
    scales="free_y", space="free"
  ) +
  theme(
    axis.text.x=element_text(size=8),
    axis.text.y=element_text(size=10),
    strip.text.x=element_text(size=8),
    strip.text.y=element_text(size=8),
    legend.position="none"
  ) +
  xlab("Number of KOs") +
  ylab("KO pathway annotation (second level)") +
  scale_y_discrete(labels=function(x) str_remove(x, "\\d+ "))

ggsave("../plots/kegg/isolate_clusters_ko_pathways.png",
       height=10, width=14)

##############################################################################
##############################################################################

odds_ratios_all_clusters
  
cat(sprintf("Top KO hits include %d unique KOs\n", length(unique(odds_ratios_all_clusters$KO))))

ko_info_all <- readRDS("../data/kegg/ko_lookup_raw.rds")
parsed_ko_pathways_all <- pbmclapply(ko_info_all, parse_brite, mc.cores=20) %>%
  bind_rows()

odds_ratios_all_clusters %>%
  left_join(parsed_ko_pathways_all, by="KO") %>%
  filter(
    !is.na(term1),
    !term1 %in% c("Not Included in Pathway or Brite", "Brite Hierarchies")
  ) %>%
  mutate(term1=str_remove(term1, "\\d+$")) %>%
  right_join(
    cluster_sizes %>% filter(cluster_size>=6),
    by=c("cluster", "cluster_longname")
  ) %>%
  ggplot(aes(x=term3, y=Odds_ratio)) +
  geom_jitter(size=0.5) +
  scale_y_log10() +
  facet_grid(rows=vars(cluster_longname), cols=vars(term1), 
             scales="free_x", space="free") +
  theme(
    axis.text.x=element_text(size=8, angle=90, hjust=1, vjust=0.5)
  )
