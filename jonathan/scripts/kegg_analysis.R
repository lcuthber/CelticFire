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

options(stringsAsFactors=FALSE)

##############################################################################
# Load data
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

# save isolate phylogenetic clustering and cluster labels
saveRDS(isolate_hclust, "../results/kegg/isolate_phylogenetic_hclust.rds")
saveRDS(cluster_labels, "../results/kegg/isolate_cluster_labels.rds")
cluster_labels %>%
  as.data.frame() %>%
  rownames_to_column("isolate") %>%
  rename(cluster=".") %>%
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
cluster_names <- tibble(
  cluster_longname=c(
    "Strep.I", "Strep.II",  "Veillonella", "Rothia", 
    "Gemella", "Prevotella", "Micrococcus", "Cupriayidus", "Pauljensenia", "Staphy.+Nialla",
    "Haemophilus", "Granulicatella", "Fusobacterium+Leptotrichia",  "Cutibacterium", "Actinomyces"
  )
) %>%
  rownames_to_column("cluster") %>%
  mutate(cluster=as.integer(cluster))

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
    # color_branches(clusters=label_colours$cluster, groupLabels=paste0("cluster_", legend_keys$cluster)) %>%
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
or_result <- list()
for(this_cluster_name in unique(y)) {
  cat(this_cluster_name, "\n")
  
  yy <- as.integer(y==this_cluster_name)
  
  # fit one vs all logistic regression model
  odds_ratios <- pbmcapply::pbmclapply(
    setNames(colnames(emapper), colnames(emapper)),
    function(ko_name) calculate_or(emapper[[ ko_name ]], yy),
    mc.cores=10
  )
  
  or_result[[ paste0("cluster_", this_cluster_name) ]] <- odds_ratios %>%
    rbindlist(idcol="KO") 
}

saveRDS(or_result, "../results/kegg/cluster_ko_odds_ratios.rds")

# combine KO scores from all clusters
odds_ratios_all_clusters <- or_result %>%
  rbindlist(idcol="cluster") %>%
  mutate(cluster=as.integer(str_remove(cluster, "cluster_"))) %>%
  left_join(cluster_names, by="cluster")
stopifnot(length(unique(odds_ratios_all_clusters$KO))==length(duplicate_ko_res$non_duplicated_columns))

odds_ratios_all_clusters %>%
  ggplot(aes(x=Odds_ratio)) +
  geom_histogram() +
  scale_x_log10() +
  lemon::facet_rep_wrap(~cluster_longname) +
  geom_vline(xintercept=1, colour="red")

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
ko_lookup_slim <- fread("../data/kegg/formatted/formatted_ko_lookup.tsv", sep="\t")

# add duplicated KOs
cluster_ko_scores <- odds_ratios_all_clusters %>%
  mutate(score=log10(Odds_ratio)) %>%
  select(cluster, KO, score)

cluster_ko_scores_15clusters <- cluster_ko_scores %>%
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
  mutate(KO=coalesce(KO, base_KO)) %>%
  left_join(ko_lookup_slim, by="KO") %>%
  relocate(KO, base_KO, symbol, name, module, pathways, n_isolates) %>%
  arrange(base_KO, KO) %>%
  mutate_if(is.numeric, function(x) round(x, digits=3))

cluster_ko_scores_15clusters %>%
  fwrite("../results/kegg/cluster_ko_scores_15clusters.tsv", sep="\t")

############################################################################################################################
############################################################################################################################





