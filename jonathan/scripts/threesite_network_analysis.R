library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(phyloseq)
library(ggplot2)
theme_set(theme_bw(base_size=14))
library(ggpubr)
library(microbiome)
library(ALDEx2)
library(tidybayes)
library(rstatix)
library(dendextend)

source("../scripts/plot_utils.R")

#####################################################################################
# Load data
#####################################################################################

# load phyloseq 
load("../data/Bus_CF_working_noduplicates.Rdata")

# OTU tables
X <- list()
X[["cf_lll"]] <- t(as.matrix(otu_table(subset_samples(Bus_CF1, SwabSite=="Brushing LLL 16S"))))
X[["cf_lul"]] <- t(as.matrix(otu_table(subset_samples(Bus_CF1, SwabSite=="Brushing LUL 16S"))))
X[["cf_ots"]] <- t(as.matrix(otu_table(subset_samples(Bus_CF1, SwabSite=="OTS"))))

X$cf_lll <- X$cf_lll[ rownames(X$cf_lll)!="11595", ] # remove repeated sample

# remove individual with likely undiagnosed CF
X$cf_lul <- X$cf_lul[ rownames(X$cf_lul)!="12405", ] 
X$cf_lll <- X$cf_lll[ rownames(X$cf_lll)!="12403", ] 
X$cf_ots <- X$cf_ots[ rownames(X$cf_ots)!="12402", ]


lapply(X, dim)

all_sample_ids <- lapply(X, function(x) tibble(sample_id=rownames(x))) %>%
  bind_rows(.id="site")

sample_metadata <- sample_data(Bus_CF1) %>% 
  as_tibble() %>%
  filter(.sample %in% all_sample_ids$sample_id) %>%
  rename(sample_id=.sample)

sample_metadata %>%
  select(sample_id, Disease) %>%
  left_join(all_sample_ids, by="sample_id") %>%
  group_by(site, Disease) %>% tally()

# map samples to individuals
study_ids <- readxl::read_xlsx(
  "../data/metadata/Combined Celtic Fire and Busselton16S samples metadata for microbiome analyses correct.xlsx"
) %>%
  select("#SampleID", "Study_ID") %>%
  rename(sample_id="#SampleID") %>%
  filter(sample_id %in% all_sample_ids$sample_id) %>%
  right_join(all_sample_ids, by="sample_id")

stopifnot(all(all_sample_ids$sample_id %in% study_ids$sample_id))

# check that all study IDs only appear once
study_ids %>%
  group_by(site, Study_ID) %>% 
  tally() %>%
  pull(n) %>%
  unique()

study_ids %>%
  left_join(sample_metadata %>% select(sample_id, Disease), by="sample_id") %>%
  group_by(site, Disease) %>% 
  tally()

# taxonomic classifications
taxonomy_tbl <- speedyseq::tax_table(Bus_CF1) %>%
  as_tibble() %>%
  rename(otu=.otu)

#####################################################################################
#####################################################################################

#####################################################################################
# Bootstrap cluster co-occurrence
#####################################################################################

library(cluster)
library(ClusBoot)
library(corrr)
library(igraph)
library(ComplexHeatmap)

compute_summary_stats <- function(x) {
  tibble(
    otu=colnames(x),
    prevalence=apply(x, 2, function(xx) 100*sum(xx>0)/length(xx)),
    total_reads=colSums(x),
    total_reads_per_sample=colSums(x)/nrow(x),
    relative_total_reads=100*colSums(x)/sum(x)
  )
}

cor_dist <- function(cor_mat) as.dist(1-cor_mat)

do_cor_clust <- function(x, k) {
  out <- cor(t(x), method="spearman") %>%
    cor_dist() %>%
    hclust() %>%
    cutree(k=k)
  return(out %>% unname())
}

otu_summary <- lapply(X, compute_summary_stats) %>%
  bind_rows(.id="site") %>%
  left_join(taxonomy_tbl, by="otu")

plotted_otus <- otu_summary %>%
  group_by(site) %>%
  slice_max(total_reads, n=100) %>%
  ungroup() %>%
  pull(otu) %>%
  unique()

X_subset <- lapply(X, function(x)x %>%
                     compositions::clr() %>%
                     .[,plotted_otus])



clustboot_propr <- list()
silhouette_vals <- list()
for(k in 2:20) {
  cat(k, "\n")
  
  clustboot_res_this_it <- pbmcapply::pbmclapply(
    X_subset,
    function(x) clusboot(t(x), B=1000, clustering.func=purrr::partial(do_cor_clust, k=k)),
    mc.cores=3
  )
  
  clustboot_propr[[ as.character(k) ]] <- lapply(clustboot_res_this_it, function(x) x$proportions)
  silhouette_vals[[ as.character(k) ]] <- lapply(clustboot_res_this_it,
                                                 function(x) tibble(silhouette=boot.silhouette(x)) %>%
    mutate(cluster=1:n())
  )
}

# plot clustboot co-cluster proportion
clust_propr_mats <- pbmcapply::pbmclapply(
  unlist(clustboot_propr, recursive=FALSE),
  function(x) x %>%
     as_cordf(diagonal=diag(.)) %>%
     stretch(),
  mc.cores=10
  ) %>%
  bind_rows(.id="k.site") %>%
  separate(k.site, c("k", "site"), sep="\\.") %>%
  mutate(k=as.integer(k))
  
# clust_propr_mats %>%
#   ggplot(aes(x=x, y=y, fill=r)) +
#   geom_tile() +
#   facet_grid(rows=vars(site), cols=vars(k)) +
#   scale_fill_gradient(low="white", high="red")
#   
# booted silouette values
all_silhouette_vals <- silhouette_vals %>%
  unlist(recursive=FALSE) %>%
  bind_rows(.id="k.site") %>%
  separate(k.site, c("k", "site"), sep="\\.") %>%
  mutate(k=as.integer(k)) 

sil_plot <- all_silhouette_vals %>%
  group_by(site, k) %>%
  summarise(mean_sil=mean(silhouette, na.rm=TRUE)) %>%
  ggplot(aes(x=k, y=mean_sil, colour=site)) +
  geom_line() +
  geom_point()

# gap statistic
wrapped_do_cor_clust <- function(...) list(cluster=do_cor_clust(...))
booted_gapstats <- lapply(
  X_subset,
  function(x) clusGap(x, K.max=20, FUNcluster=wrapped_do_cor_clust, B=100))
all_booted_gapstats <- lapply(
  booted_gapstats,
  function(x) as_tibble(x$Tab) %>% mutate(k=1:n())) %>%
  bind_rows(.id="site")
gap_plot <- all_booted_gapstats  %>%
  ggplot(aes(x=k, y=gap, colour=site)) +
  geom_line() +
  geom_point()

cowplot::plot_grid(sil_plot, gap_plot, nrow=1, commmon.legend=TRUE)

all_booted_gapstats %>%
  group_by(site, k) %>%
  summarise(mean_gap=mean(gap, na.rm=TRUE)) %>%
  group_by(site) %>%
  summarise(k=k[ which.max(mean_gap) ]) %>%
  ungroup()

# best_k <- all_silhouette_vals %>%
#   group_by(site, k) %>%
#   summarise(mean_sil=mean(silhouette, na.rm=TRUE)) %>%
#   group_by(site) %>%
#   summarise(k=k[ which.max(mean_sil) ]) %>%
#   ungroup()

best_k <- tribble(
  ~site, ~k,
  "cf_lll", 9,
  "cf_lul", 9,
  "cf_ots", 7
)

clust_propr_mats %>%
  right_join(best_k, by=c("k", "site")) %>%
  ggplot(aes(x=x, y=y, fill=r)) +
  geom_tile() +
  facet_wrap(~site) +
  scale_fill_gradient(low="white", high="red")


propr_mats_best_k <- clust_propr_mats %>%
  right_join(best_k, by=c("k", "site")) %>%
  group_by(site) %>%
  do(mat=
       select(., x, y, r) %>%
       pivot_wider(names_from=y, values_from=r) %>%
       column_to_rownames("x") %>%
       as.matrix()
  ) %>%
  rowwise() %>%
  mutate(
    mat=list(mat[plotted_otus,plotted_otus]),
    final_clustering=as.dist(1-mat) %>% hclust() %>% cutree(h=0.5) %>% list()
  )

obs_cor_mats <- lapply(
  X_subset,
  function(x) x %>% cor(method="spearman")
)
obs_cor_hclust <- lapply(
  obs_cor_mats,
  function(x) x %>% cor_dist() %>% hclust()
)

# OTU clusters at each site
stopifnot(identical(best_k$site, names(obs_cor_hclust)))
otu_cluster_labels <- mapply(
  function(hc, k) hc %>%
    cutree(k=k) %>%
    as.data.frame() %>%
    rownames_to_column("otu") %>%
    rename(cluster=".") %>%
    as_tibble(),
  obs_cor_hclust,
  best_k$k,
  SIMPLIFY=FALSE
) %>%
  bind_rows(.id="site") %>%
  pivot_wider(names_from=site, values_from=cluster)

otu_cluster_labels %>%
  select(-otu) %>%
  colpair_map(function(x,y) WGCNA::randIndex(table(x, y)))

as.dendlist(lapply(obs_cor_hclust, as.dendrogram)) %>%
  cor.dendlist()

# taxonomic annotatoin for proportionality heatmap
heatmap_tax_anno_data <- taxonomy_tbl %>%
  select(otu, Phylum) %>%
  filter(otu %in% plotted_otus) %>%
  arrange(match(otu, plotted_otus)) %>%
  left_join(otu_cluster_labels, by="otu")
unique_phyla <- unique(heatmap_tax_anno_data$Phylum)
phylum_colours <- setNames(scales::brewer_pal("Set2", type="qual")(length(unique_phyla)), unique_phyla)

anno_col_palettes <- c(Phylum="Set2", cf_lll="Set1", cf_lul="Paired", cf_ots="Set1")

anno_col_list <- lapply(
  setNames(names(anno_col_palettes), names(anno_col_palettes)),
  function(x) setNames(
    scales::brewer_pal(anno_col_palettes[[x]], type="qual")(length(unique(heatmap_tax_anno_data[[x]]))), unique(heatmap_tax_anno_data[[x]])
  )
)

get_corr_pvals <- function(x) {
  out <- psych::corr.test(x, method="spearman", adjust="fdr")
  p_adj_mat <- out$p
  p_adj_mat[lower.tri(p_adj_mat)] <- 0
  p_adj_mat <- p_adj_mat + t(p_adj_mat)
  diag(p_adj_mat) <- NA
  return(p_adj_mat)
}

corr_test_pvals <- lapply(
  X_subset, get_corr_pvals
)

#
# co-occurence in clusteres for each site
make_heatmap <- function(im, var_hclust,
                         site_name, draw_pvals,
                         show_dend=TRUE, symm_colbar=TRUE, ...) {

  if(symm_colbar) {
    cpal <- circlize::colorRamp2(c(-1, 0, 1), c("Darkblue", "white", "red"))
  } else {
    cpal <- circlize::colorRamp2(c(0,1), c("white", "darkgreen"))
  }
  
  Heatmap(
    im,
    col=cpal,
    cluster_rows=var_hclust,
    cluster_columns=var_hclust,
    row_names_gp=gpar(fontsize=9),
    column_names_gp=gpar(fontsize=9),aa
    row_dend_width=unit(30, "mm"),
    column_dend_height=unit(30, "mm"),
    show_column_dend=show_dend,
    show_row_dend=show_dend,
    top_annotation=columnAnnotation(
      phylum=setNames(heatmap_tax_anno_data$Phylum, heatmap_tax_anno_data$otu),
      cluster=setNames(heatmap_tax_anno_data[[ site_name ]], heatmap_tax_anno_data$otu),
      col=list(phylum=anno_col_list$Phylum, cluster=anno_col_list[[site_name]])
    ),
    left_annotation=rowAnnotation(
      phylum=setNames(heatmap_tax_anno_data$Phylum, heatmap_tax_anno_data$otu),
      cluster=setNames(heatmap_tax_anno_data[[ site_name ]], heatmap_tax_anno_data$otu),
      show_legend=FALSE,
      col=list(phylum=anno_col_list$Phylum, cluster=anno_col_list[[site_name]])
    ),
    width=ncol(im)*unit(2.5, "mm"), 
    height=nrow(im)*unit(2.5, "mm"),
    cell_fun = function(j, i, x, y, w, h, fill) {
      if(draw_pvals) {
        if(!is.na(corr_test_pvals[[ site_name ]][i, j])) {
          if(corr_test_pvals[[ site_name ]][i, j] < 0.01) {
            grid.text('*', x, y, vjust=0.75)
          }
        }
      }
    },
    ...
  )
}

# heatmaps comparing correlation and stability
for(site_name in names(obs_cor_mats)) {
  cat(site_name, "\n")
  
  h1 <- make_heatmap(obs_cor_mats[[ site_name ]],
                     obs_cor_hclust[[ site_name ]],
                     site_name=site_name, draw_pvals=TRUE,
                     symm_colbar=TRUE, name="corr")
  h2 <- make_heatmap(propr_mats_best_k$mat[ propr_mats_best_k$site==site_name ][[1]],
                     obs_cor_hclust[[ site_name ]],
                     site_name=site_name, draw_pvals=FALSE,
                     show_dend=FALSE, symm_colbar=FALSE, name="cluster\nstability")
  
  plot_filename <- sprintf("../plots/cf_lll_vs_lul/corr_heatmap__%s.png", site_name)

  png(plot_filename, height=1050, width=1800, res=100)
  draw(h1+h2,
       merge_legend=TRUE,
       column_title=site_name,
       column_title_gp=gpar(fontsize=24))
  dev.off()
  # knitr::plot_crop(plot_filename)
}

propr_mats_best_k$propr_hclust <- lapply(propr_mats_best_k$mat, function(x) as.dist(1-x) %>% hclust())
propr_mats_best_k$propr_clustering <- lapply(propr_mats_best_k$propr_hclust, function(x) x %>% cutree(h=0.5))

# define the clusters on the stability rather than the correlation values
for(i in 1:nrow(propr_mats_best_k)) {
  site_name <- propr_mats_best_k$site[[ i ]]
  cat(site_name, "\n")
  make_heatmap(
    propr_mats_best_k$mat[[i]], propr_mats_best_k$propr_hclust[[i]],
    site_name=propr_mats_best_k$site[[i]], draw_pvals=FALSE,
    symm_colbar=FALSE
  ) -> h3
  
  pdf(sprintf("../plots/cf_lll_vs_lul/corr_heatmap__stability__%s.pdf", site_name),
      height=12, width=13)
  draw(h3, merge_legend=TRUE)
  dev.off()
}



#
# differential correlation analysis
library(boot)
library(ggraph)
library(tidygraph)

cor_mat_stat <- function(data, indices, cor.type="spearman") {
  out <- cor(data[indices,], method=cor.type) %>%
    as_cordf() %>%
    # shave() %>%
    stretch() %>%
    pull(r)
    # filter(x!=y) %>%
    # pull(r)
  # return()
}

# This returns the column and row-names in the transposed order
# only works because it's a symmetric matrix
get_obs_cor <- function(data, cor.type="spearman") {
  data %>%
    cor(method=cor.type) %>%
    as_cordf() %>%
    stretch()
}

# heatmaps comparing correlation and stability
sample_groups <- sample_metadata %>% select(sample_id, Disease)
booted_cor_adjacencies <- list()
for(site_name in names(X_subset)) {
  for(condition_name in c("all", "Asthma", "Control")) {
    
    cat(site_name, "\n")
    
    xx <- X_subset[[ site_name ]]
    
    if(condition_name!="all") {
      .ids <- sample_groups %>% filter(sample_id %in% rownames(xx), Disease==condition_name) %>% pull(sample_id)
      xx <- xx[ .ids, ]
    }
    print(dim(xx))
    
    boot_strap <- boot(data=xx,
                       statistic=cor_mat_stat, R=100, parallel="multicore", ncpus=20)
    booted_cor_adjacencies[[ paste0(site_name, "----", condition_name) ]] <- boot_strap$t %>%
      matrixStats::colQuantiles(probs=c(0.025, 0.5, 0.975), na.rm=TRUE) %>%
      as_tibble() %>%
      bind_cols(xx %>% get_obs_cor()) %>%
      relocate(x,y,r)
  }
}

make_unique_key <- function(x, y) {
  stopifnot(length(x)==length(y))
  unique_xy <- unique(c(x, y)) %>% sort()
  x_fac <- factor(x, levels=unique_xy)
  y_fac <- factor(y, levels=unique_xy)
  
  out <- rep(NA, length(x))
  out[ as.integer(x_fac)>as.integer(y_fac) ] <- sprintf("%s--%s", y[as.integer(x_fac)>as.integer(y_fac)], x[as.integer(x_fac)>as.integer(y_fac)])
  out[ as.integer(x_fac)<=as.integer(y_fac) ] <- sprintf("%s--%s", x[as.integer(x_fac)<=as.integer(y_fac)], y[as.integer(x_fac)<=as.integer(y_fac)])
  stopifnot(sum(is.na(out))==0)

  return(out)
}

#
# make adjacency matrices
# zero weight if CI includes zero, otherwise the bootstrapped value (50th percentile)
adj_df_all <- booted_cor_adjacencies %>%
  bind_rows(.id="label") %>%
  separate(label, c("site", "group"), sep="----") %>%
  filter(group!="all") %>%
  mutate(key=make_unique_key(x, y)) %>%
  select(-c(x,y)) %>%
  distinct() %>%
  mutate(weight=if_else(
    data.table::between(0, `2.5%`, `97.5%`),
    0,
    r
  )) %>%
  filter(!is.na(r))

adj_df_all %>%
  group_by(site) %>%
  summarise(sum(weight!=0))


all_adj_lists <- adj_df_all %>%
  group_by(site, group) %>%
  do(
    adj_mat=select(., key, weight) %>%
      separate(key, c("V1", "V2"), sep="--") %>%
      mutate(V1=factor(V1, levels=plotted_otus),
             V2=factor(V2, levels=plotted_otus)) %>%
      complete(V1, V2, fill=list(weight=0)) %>%
      pivot_wider(names_from=V2, values_from=weight) %>%
      column_to_rownames("V1") %>%
      as.matrix()
  )
all_adj_mats$adj_mat <- lapply(all_adj_mats$adj_mat, function(x) x[ plotted_otus,plotted_otus ])

make_symmetric <- function(x) {
  x2 <- x + t(x)
  diag(x2) <- 1
  return(x2)
}

all_adj_mats$g <- lapply(
  all_adj_mats$adj_mat,
  function(x) x %>% graph_from_adjacency_matrix(mode="max", weighted="weight")
)

g <- all_adj_mats$g[[2]]
lou <- g %>% cluster_louvain()
LO <- layout_with_fr(g, weights=E(g)$weight )
plot(lou, g, vertex.label = NA, vertex.size=5,  edge.width=E(g)$width,
     edge.arrow.size = .2, layout=LO)



make_heatmap_3 <- function(x, title) {
  im <- make_symmetric(x)
 #  im_hclust <- as.dist(1-im) %>% hclust()
  Heatmap(
    im,
    # cluster_rows=im_hclust, cluster_columns=im_hclust,
    col=circlize::colorRamp2(c(-1, 0, 1), c("Darkblue", "white", "red")),
    row_dend_width=unit(35, "mm"), column_dend_height=unit(35, "mm"),
    height=unit(7, "inch"), width=unit(7, "inch"),
    row_names_gp=gpar(fontsize=7),
    column_names_gp=gpar(fontsize=7),
    column_title=title
  )
}

all_adj_mats$hm <- mapply(make_heatmap_3, all_adj_mats$adj_mat, all_adj_mats$group)

for(row_idx in 1:nrow(all_adj_mats)) {
  site_name <- all_adj_mats$site[[ row_idx ]]
  group_name <- all_adj_mats$group[[ row_idx ]]
  cat(site_name, group_name, "\n")
  
  save_path <- sprintf("../plots/cf_lll_vs_lul/heatmap_%s_%s.pdf", site_name, group_name)
  pdf(save_path, width=12, height=12)
  draw(all_adj_mats$hm[[ row_idx ]])
  dev.off()
  knitr::plot_crop(save_path)
}

#
# bb
library(dendextend)
dend1 <- all_adj_mats$adj_mat[[1]] %>%
  make_symmetric() %>%
  dist("euclidean") %>%
  hclust() %>%
  as.dendrogram()
dend_1_clusters <- dend1 %>% dendextend::cutree(k=3)

dend2 <- all_adj_mats$adj_mat[[2]] %>%
  make_symmetric() %>%
  dist("euclidean") %>%
  hclust() %>%
  as.dendrogram()
dend2 %>% plot()

dend1_clusters <- dynamicTreeCut::cutreeDynamic(
  dend1 %>% as.hclust(),
  distM=all_adj_mats$adj_mat[[1]] %>%
    make_symmetric() %>%
    dist() %>%
    as.matrix(),
  minClusterSize=2
) %>%
  unname() %>%
  setNames(., colnames(all_adj_mats$adj_mat[[1]]))


dend2_clusters <- dynamicTreeCut::cutreeDynamic(
  dend2 %>% as.hclust(),
  distM=all_adj_mats$adj_mat[[2]] %>%
    make_symmetric() %>%
    dist() %>%
    as.matrix(),
  minClusterSize=2
) %>%
  unname() %>%
  setNames(., colnames(all_adj_mats$adj_mat[[2]]))

cbind(dend1_clusters, dend2_clusters)

WGCNA::plotDendroAndColors(
  dend2, WGCNA::labels2colors(dend2_clusters)
)

par(mfrow=c(1,1))
dend2 %>%
  colour_branches(col=WGCNA::labels2colors(dend2_clusters)[ order.dendrogram(.) ]) %>%
  plot()


library(ggraph)
library(tidygraph)

remove_unconnected_nodes <- function(g) {
  Isolated <- which(degree(g)==0)
  g = delete.vertices(g, Isolated)
  return(g)
}

intersection(all_adj_mats$g[[1]], all_adj_mats$g[[2]], keep.all.vertices=FALSE) %>% 
  remove_unconnected_nodes() %>%
  as_tbl_graph() %>%
  left_join(taxonomy_tbl %>%
              select(otu, Genus) %>%
              mutate(Genus=sub("_[^_]+$", "", Genus)),
            by=c(name="otu")) %>%
  ggraph() +
  geom_edge_link(alpha=0.5) +
  geom_node_point(aes(colour=Genus), size=4)


difference(all_adj_mats$g[[1]], all_adj_mats$g[[2]], keep.all.vertices=FALSE) %>%
  remove_unconnected_nodes() %>%
  as_tbl_graph() %>%
  left_join(taxonomy_tbl %>%
              select(otu, Genus) %>%
              mutate(Genus=sub("_[^_]+$", "", Genus)),
            by=c(name="otu")) %>%
  ggraph() +
  geom_edge_link(alpha=0.5) +
  geom_node_point(aes(colour=Genus), size=4)

difference(all_adj_mats$g[[2]], all_adj_mats$g[[1]], keep.all.vertices=FALSE) %>%
  remove_unconnected_nodes() %>%
  as_tbl_graph() %>%
  left_join(taxonomy_tbl %>%
              select(otu, Genus) %>%
              mutate(Genus=sub("_[^_]+$", "", Genus)),
            by=c(name="otu")) %>%
  ggraph() +
  geom_edge_link(alpha=0.5) +
  geom_node_point(aes(colour=Genus), size=4)





all_adj_mats$g[[1]] - all_adj_mats$g[[2]]




corr_sim_mat <- corr_sim_mat[ plotted_otus,plotted_otus ]
stopifnot(identical(colnames(corr_sim_mat), rownames(corr_sim_mat)))
corr_sim_mat <- corr_sim_mat + t(corr_sim_mat)

# library(apcluster)
# tmp <- apcluster(s=corr_sim_mat, details=TRUE)
# plot(tmp, corr_sim_mat)
# 
# ap_clusters <- lapply(tmp@clusters, function(x) tibble(otu=names(x))) %>%
#   bind_rows(.id="cluster")
# ap_clusters %>%
#   arrange(cluster)

# ap_clusters


#
# differential correlations in asthmatics vs controls - two sample
get_cor_diff <- function(data, indices, cor.type="pearson") {
  
  cor_g0 <- data[indices,] %>%
    filter(group==0) %>%
    select(-group) %>%
    correlate(method=cor.type, quiet=TRUE) %>%
    stretch()
  cor_g1 <- data[indices,] %>%
    filter(group==1) %>%
    select(-group) %>%
    correlate(method=cor.type, quiet=TRUE) %>%
    stretch()
  
  cor_g0 %>%
    inner_join(cor_g1, by=c("x", "y")) %>%
    mutate(r_diff=r.x-r.y) %>%
    pull(r_diff)
}

make_heatmap_4 <- function(x, title, name="", ...) {
  im <- x
  im_hclust <- as.dist(1-abs(im)) %>% hclust()
  Heatmap(
    im,
    name=name,
    # cluster_rows=im_hclust, cluster_columns=im_hclust,
    col=circlize::colorRamp2(c(-1, 0, 1), c("Darkblue", "white", "red")),
    row_dend_width=unit(50, "mm"), column_dend_height=unit(50, "mm"),
    height=unit(7, "inch"), width=unit(7, "inch"),
    row_names_gp=gpar(fontsize=6),
    column_names_gp=gpar(fontsize=6),
    column_title=title,
    ...
  )
}

# heatmaps comparing correlation and stability
sample_groups <- sample_metadata %>% select(sample_id, Disease)

#
# site-wise
booted_diffcorr <- list()
for(site_name in names(X_subset)) {
  cat(site_name, "\n")
  
  x <- X_subset[[ site_name ]] %>%
    as.data.frame()
  
  x_groups <- sample_groups %>%
    filter(sample_id %in% rownames(x)) %>%
    mutate(Disease_enc=recode(Disease, Control=0, Asthma=1))
  x_group <- setNames(x_groups$Disease_enc, x_groups$sample_id)[ rownames(x) ]
  
  x <- cbind(x, x_group) %>%
    as.data.frame() %>%
    rename(group=x_group)
  
  boot_strap <- boot(data=x,
                     strata=x$group,
                     statistic=get_cor_diff, R=1000, parallel="multicore", ncpus=10)
  n_nas <- sum(is.na(boot_strap$t))
  cat(sprintf("There are %d NAs in bootstrapped estimates\n", n_nas))
  
  booted_diffcorr[[ site_name ]] <- boot_strap$t %>%
    matrixStats::colQuantiles(probs=c(0.05, 0.5, 0.95), na.rm=TRUE) %>%
    as_tibble() %>%
    bind_cols(x %>% select(-group) %>% correlate(quiet=TRUE) %>% stretch() %>% select(x,y)) %>%
    relocate(x,y)
}

pdf("../plots/cf_lll_vs_lul/diffcor_heatmap_cf_lul.pdf", height=11, width=12)
booted_diffcorr[[2]] %>%
  filter(!data.table::between(0, `5%`, `95%`)) %>%
  select(x, y, `50%`) %>%
  rename(weight=`50%`) %>%
  pivot_wider(names_from=y, values_from=weight, values_fill=0) %>%
  column_to_rownames("x") %>%
  as.matrix() %>%
  make_heatmap_4("", "Control - Asthmatics", row_split=4, column_split=4)
dev.off()

#
# check which 

#
# between pairs of sites
site_pairs <- tibble(
  site=names(X_subset), value=1
) %>%
  gather_pairs(key=site, value=value) %>%
  select(.row, .col) %>%
  ungroup() %>%
  mutate_all(as.character)


booted_diff_cors <- list()

for(site_pair_idx in 1:nrow(site_pairs)) {
  
  cat(site_pair_idx, "\n")
  
  site_name_1 <- site_pairs$.row[[ site_pair_idx ]]
  site_name_2 <- site_pairs$.col[[ site_pair_idx ]]
  
  xx <- rbind(X_subset[[ site_name_1 ]], X_subset[[ site_name_2 ]]) %>%
    as.data.frame()
  xx$group <- c(rep(0, nrow(X_subset[[ site_name_1]])), rep(1, nrow(X_subset[[ site_name_2]]) ))
  
  boot_strap <- boot(data=xx,
                     statistic=get_cor_diff, R=1000,
                     strata=xx$group,
                     parallel="multicore", ncpus=20)
  
  booted_diff_cors[[ site_pair_idx ]] <- boot_strap$t %>%
    matrixStats::colQuantiles(probs=c(0.025, 0.5, 0.975), na.rm=TRUE) %>%
    as_tibble() %>%
    bind_cols(xx %>% select(-group) %>% correlate(quiet=TRUE) %>% stretch() %>% select(x,y)) %>%
    relocate(x,y)
}



booted_diff_cors[[2]] %>%
  filter(!data.table::between(0, `2.5%`, `97.5%`)) %>%
  select(x, y, `50%`) %>%
  rename(weight=`50%`) %>%
  pivot_wider(names_from=y, values_from=weight, values_fill=0) %>%
  column_to_rownames("x") %>%
  as.matrix() %>%
  make_heatmap_4("", "")

#
# differential correlations in asthmatics vs controls
diffcor_res <- booted_cor_adjacencies %>%
  bind_rows(.id="label") %>%
  separate(label, c("site", "group"), sep="----") %>%
  mutate(key=make_unique_key(x, y)) %>%
  select(-c(x,y)) %>%
  distinct()

get_cor_type <- function(cor_lower, cor_upper, eps) {
  stopifnot(length(cor_lower)==length(cor_upper))
  stopifnot(all(cor_lower<cor_upper, na.rm=TRUE))
  out <- rep(NA, length(cor_lower))
  
  if(eps==0) {
    out[ data.table::between(0, cor_lower, cor_upper) ] <- "zero"
    out[ cor_lower>0 ] <- "positive"
    out[ cor_upper<0 ] <- "negative"
  } else {
    stopifnot(eps>0)
    out[ cor_lower<eps | cor_lower>eps ] <- "zero"
    out[ cor_lower>eps ] <- "positive"
    out[ cor_upper<-eps ] <- "negative"
  }
  
  # bind_cols(cor_lower, cor_upper, out)
  return(out)
}

diff_network_tbl <- diffcor_res %>%
  filter(group!="all") %>%
  select(site, group, key, `5%`, `95%`) %>%
  mutate(cor_type=get_cor_type(`5%`, `95%`, 0.5)) %>%
  select(site, group, key, cor_type) %>%
  pivot_wider(names_from=group, values_from=cor_type) %>%
  mutate(diff_cor_type=sprintf("%s--%s", Asthma, Control)) %>%
  select(-c(Asthma, Control)) %>%
  separate(key, c("node1", "node2"), sep="--")


for(site_name in c("cf_lll", "cf_lul", "cf_ots")) {
  cat(site_name, "\n")

  g <- diff_network_tbl %>%
    filter(site==site_name,
           !diff_cor_type %in% c("zero--zero", "positive--positive", "negative--negative")) %>%
    # mutate(weight=if_else(diff_cor_type!="zero--zero", diff_cor_type, NA_character_)) %>%
    select(-site) %>%
    as_tbl_graph(directed=FALSE)
  
  g2 <- delete_edges(g, which(E(g)$diff_cor_type == "zero--zero"))
  g2 <- delete_edges(g2, which(grepl("NA", E(g2)$diff_cor_type)))
  g2 <- g2 %>%
    as_tbl_graph() %>%
    mutate(node_label=sub("_[^_]+$", "", name)) %>%
    left_join(taxonomy_tbl %>%
                select(otu, Phylum) %>%
                mutate(Phylum=sub("_[^_]+$", "", Phylum)),
              by=c(name="otu")) %>%
    left_join(otu_cluster_labels %>%
                select(otu, !!sym(site_name)) %>%
                rename(name=otu, cluster=!!sym(site_name))) %>%
    arrange(cluster)
  
  plt <- g2 %>%
    arrange(Phylum, name) %>%
    ggraph(layout="linear", circular=TRUE) +
    # geom_node_label(aes(label=node_label), alpha=0.5) +
    geom_node_point(aes(colour=Phylum), size=3) +
    geom_edge_link(aes(colour=diff_cor_type)) +
    # facet_wrap(~diff_cor_type)+
    theme(legend.position="top")# +
    ggforce::geom_mark_hull(aes(fill=Genus, x=x, y=y))

    
  arc_tbl <- plt$data %>%
    group_by(Phylum) %>%
    slice(1, n()) %>%
    ungroup() %>%
    mutate(
      r=sqrt(x**2 + y**2), theta=atan2(x, y),
      tmp=rep(c("start", "end"), 5)
    ) %>%
    select(Phylum, theta, tmp) %>%
    pivot_wider(names_from=tmp, values_from=theta) %>%
    mutate(r=seq(1, 2, length.out=5)) %>%
    mutate(start_deg=start*360/(2*pi), end_deg=end*360/(2*pi),
           reverse_startend=start-end>pi)
  plt <- plt + 
    ggforce::geom_arc(data=arc_tbl %>% filter(!reverse_startend),
                      aes(start=end, end=start, x0=0, y0=0, r=1.05, colour=Phylum),
                      size=2) +
    ggforce::geom_arc(data=arc_tbl %>% filter(reverse_startend),
                      aes(start=start, end=(2*pi)+end, x0=0, y0=0, r=1.05, colour=Phylum),
                      size=2) +
    ggrepel::geom_label_repel(aes(label=name, x=x, y=y), alpha=0.5) +
    ggtitle(site_name)
  
  ggsave(sprintf("../plots/cf_lll_vs_lul/diffcor_network_%s.png", site_name),
         dpi=450, width=9, height=9,
         plot=plt)
}



#
# correlation heatmaps for each group
group_corr_heatmap_data <- booted_cor_adjacencies %>%
  bind_rows(.id="label") %>%
  separate(label, c("site", "group"), sep="----") %>%
  mutate(r=if_else(x==y, 1, r),
         cor_type=get_cor_type(`5%`, `95%`, 0.0)
  )

heatmap_data <- group_corr_heatmap_data %>%
  mutate(r_plot=if_else(cor_type!="zero", r, 0)) %>%
  filter(group!="all") %>%
  select(x, y, group, site, r_plot) %>%
  group_by(site, group) %>%
  do(mat=pivot_wider(., names_from=y, values_from=r_plot) %>%
                     column_to_rownames("x") %>%
                     select(-c(site, group)) %>%
                     as.matrix()
  )

make_heatmap2 <- function(im,
                         show_dend=TRUE, symm_colbar=TRUE, ...) {
  
  if(symm_colbar) {
    cpal <- circlize::colorRamp2(c(-1, 0, 1), c("Darkblue", "white", "red"))
  } else {
    cpal <- circlize::colorRamp2(c(0,1), c("white", "darkgreen"))
  }
  
  Heatmap(
    im,
    col=cpal,
    # cluster_rows=var_hclust,
    # cluster_columns=var_hclust,
    row_names_gp=gpar(fontsize=9),
    column_names_gp=gpar(fontsize=9),
    row_dend_width=unit(30, "mm"),
    column_dend_height=unit(30, "mm"),
    show_column_dend=show_dend,
    show_row_dend=show_dend,
    # top_annotation=columnAnnotation(
    #   phylum=setNames(heatmap_tax_anno_data$Phylum, heatmap_tax_anno_data$otu),
    #   cluster=setNames(heatmap_tax_anno_data[[ site_name ]], heatmap_tax_anno_data$otu),
    #   col=list(phylum=anno_col_list$Phylum, cluster=anno_col_list[[site_name]])
    # ),
    # left_annotation=rowAnnotation(
    #   phylum=setNames(heatmap_tax_anno_data$Phylum, heatmap_tax_anno_data$otu),
    #   cluster=setNames(heatmap_tax_anno_data[[ site_name ]], heatmap_tax_anno_data$otu),
    #   show_legend=FALSE,
    #   col=list(phylum=anno_col_list$Phylum, cluster=anno_col_list[[site_name]])
    # ),
    width=ncol(im)*unit(2.5, "mm"), 
    height=nrow(im)*unit(2.5, "mm"),
    ...
  )
}

library(ggplotify)

heatmap_data$hclust <- lapply(heatmap_data$mat, function(x) as.dist(1-x) %>% hclust())
pw_cophcor <- as.dendlist(lapply(heatmap_data$hclust, as.dendrogram)) %>%
  cor.dendlist()
rownames(pw_cophcor) <- paste0(heatmap_data$site, "--", heatmap_data$group)
colnames(pw_cophcor) <- paste0(heatmap_data$site, "--", heatmap_data$group)

pw_cophcor %>%
  # as_cordf(diagonal=1) %>%
  corrplot::corrplot(method="number")
  stretch() %>%
  ggplot(aes(x=x, y=y, fill=r)) +
  geom_tile() +
  scale_fill_gradient2()
  
tidy_mantel <- function(d1, d2, ...) {
  out <- vegan::mantel(d1, d2, ...)
  return(tibble(estimate=out$statistic, pvalue=out$signif,
                method=paste0("Mantel:", out$method)))
}


mantel_test_res <- heatmap_data %>%
  select(site, group, mat) %>%
  unite("label", site, group, sep="--") %>%
  gather_pairs(label, mat) %>%
  rowwise() %>%
  mutate(tidy_mantel(1-.x, 1-.y)) %>%
  select(.row, .col, estimate, pvalue) %>%
  mutate(p_adj=p.adjust(pvalue, method="fdr"))

mantel_test_res %>%
  separate(.row, c(".row_site", ".row_group"), sep="--") %>%
  separate(.col, c(".col_site", ".col_group"), sep="--") %>%
  filter(.row_site==.col_site) %>%
  ggplot()

mantel_test_res %>%
  ggplot(aes(x=.row, y=.col)) +
  geom_tile(aes(fill=estimate)) +
  geom_label(aes(label=format(estimate, digits=2))) +
  scale_fill_gradient(low="white", high="red", limits=c(0,1))

library(ggplotify)
h1 <- heatmap_data$mat[[3]] %>% make_heatmap2() %>% as.ggplot()
h2 <- heatmap_data$mat[[4]] %>% make_heatmap2() %>% as.ggplot()

ggarrange(h1, h2)

#
# for each site, bootstrap the difference in (i) cophenectic correlation and (ii) 
# the mantel distance between correlation matrices
sample_groups <- sample_metadata %>% select(sample_id, Disease)

get_mantel_dist <- function(data, indices, cor.type="spearman") {
  data_0 <- data[indices,] %>% 
    filter(x_group==0) %>%
    select(-x_group)
  data_1 <- data[indices,] %>% 
    filter(x_group==1) %>%
    select(-x_group)
  
  cor_data_0 <- cor(data_0, method=cor.type)
  cor_data_1 <- cor(data_1, method=cor.type)
  
  estimate <- tryCatch({
    mantel_out <- vegan::mantel(1-cor_data_0, 1-cor_data_1)
    return(mantel_out$statistic)
    },
    error=function(e) {
      return(NA)
    })
  
  return(estimate)
}

booted_mantel_dist <- list()
for(site_name in names(X_subset)) {
  cat(site_name, "\n")
  
  x <- X_subset[[ site_name ]]
  x_groups <- sample_groups %>%
    filter(sample_id %in% rownames(x)) %>%
    mutate(Disease_enc=recode(Disease, Control=0, Asthma=1))
  x_group <- setNames(x_groups$Disease_enc, x_groups$sample_id)[ rownames(x) ]
  
  x <- cbind(x, x_group) %>%
    as.data.frame()
  
  boot_strap <- boot(data=x,
                     statistic=get_mantel_dist, R=1000, parallel="multicore", ncpus=10)
  n_nas <- sum(is.na(boot_strap$t))
  cat(sprintf("There are %d NAs in bootstrapped estimates\n", n_nas))
  
  booted_mantel_dist[[ site_name ]] <- matrixStats::colQuantiles(
   boot_strap$t, probs=c(0.05, 0.5, 0.95), na.rm=TRUE
  )
}

#
# same thing for cophenectic correlation
get_cop_cor <- function(data, indices, cor.type="spearman") {
  
  estimate <- tryCatch({
    data_0 <- data[indices,] %>% 
      filter(x_group==0) %>%
      select(-x_group)
    data_1 <- data[indices,] %>% 
      filter(x_group==1) %>%
      select(-x_group)
    
    cor_data_0 <- cor(data_0, method=cor.type)
    cor_data_1 <- cor(data_1, method=cor.type)
    
    dend_data_0 <- as.dist(1-cor_data_0) %>% hclust() %>%  as.dendrogram()
    dend_data_1 <- as.dist(1-cor_data_1) %>% hclust() %>%  as.dendrogram()
    
    return(dendextend::cor_cophenetic(dend_data_0, dend_data_1))
  },
  error=function(e) {
    return(NA)
  })
  
  return(estimate)
}

booted_coph_cor <- list()
for(site_name in names(X_subset)) {
  cat(site_name, "\n")
  
  x <- X_subset[[ site_name ]]
  x_groups <- sample_groups %>%
    filter(sample_id %in% rownames(x)) %>%
    mutate(Disease_enc=recode(Disease, Control=0, Asthma=1))
  x_group <- setNames(x_groups$Disease_enc, x_groups$sample_id)[ rownames(x) ]
  
  x <- cbind(x, x_group) %>%
    as.data.frame()
  
  boot_strap <- boot(data=x,
                     statistic=get_cop_cor, R=1000, parallel="multicore", ncpus=10)
  n_nas <- sum(is.na(boot_strap$t))
  cat(sprintf("There are %d NAs in bootstrapped estimates\n", n_nas))
  
  booted_coph_cor[[ site_name ]] <- matrixStats::colQuantiles(
    boot_strap$t, probs=c(0.05, 0.5, 0.95), na.rm=TRUE
  )
}

# plot the results
mantel_plot <- booted_mantel_dist %>%
  bind_rows() %>%
  mutate(site=names(booted_mantel_dist)) %>%
  ggplot(aes(y=`50%`, x=site)) +
  geom_col() +
  geom_errorbar(aes(ymin=`5%`, ymax=`95%`)) +
  ylim(c(0,1)) +
  ylab("Network similarity\nbetween asthmatics and controls (90% CI)") +
  xlab("Site")  +
  get_plot_theme(14)
mantel_plot
cropped_ggsave("../plots/cf_lll_vs_lul/booted_mantel_dist.pdf", height=5, width=8)

cophcor_plot <- booted_coph_cor %>%
  bind_rows() %>%
  mutate(site=names(booted_coph_cor)) %>%
  ggplot(aes(y=`50%`, x=site)) +
  geom_col() +
  geom_errorbar(aes(ymin=`5%`, ymax=`95%`)) +
  ylim(c(0,1)) +
  ylab("Cophenectic correlation\nbetween groups (90% CI)") +
  xlab("Site")


ggpubr::ggarrange(mantel_plot, cophcor_plot)  +
  get_plot_theme(14)



#
# do the same for the pairwise distances between the sites
ids_by_site_pair <- study_ids %>%
  group_by(site) %>%
  summarise(ids=list(unique(sample_id))) %>%
  gather_pairs(site, ids) %>%
  ungroup()

ids_by_site_pair$mantel_dist_res <- rep(list(NA), nrow(ids_by_site_pair))

booted_mantel_dist <- list()
for(row_idx in 1:nrow(ids_by_site_pair)) {
  cat(row_idx, "\n")
  
  site0 <- ids_by_site_pair$.row[[ row_idx ]]
  site1 <- ids_by_site_pair$.col[[ row_idx ]] 
  
  x0 <- X_subset[[ site0 ]] %>% as.data.frame() %>% mutate(x_group=0)
  x1 <- X_subset[[ site1 ]] %>% as.data.frame() %>% mutate(x_group=1)
  
  x <- bind_rows(x0, x1)
  
  boot_strap <- boot(data=x,
                     statistic=get_mantel_dist, R=1000, parallel="multicore", ncpus=10)
  n_nas <- sum(is.na(boot_strap$t))
  cat(sprintf("There are %d NAs in bootstrapped estimates\n", n_nas))
  
  ids_by_site_pair$mantel_dist_res[[ row_idx ]] <- matrixStats::colQuantiles(
    boot_strap$t, probs=c(0.05, 0.5, 0.95), na.rm=TRUE
  )
}


ids_by_site_pair %>%
  unnest_wider(mantel_dist_res) %>%
  mutate(y_label=sprintf("%s--%s", .row, .col)) %>%
  ggplot(aes(y=y_label, x=`50%`)) +
  geom_point() +
  geom_errorbarh(aes(xmin=`5%`, xmax=`95%`)) +
  xlab("Mantel similarity (90% CI)") +
  ylab("Site pair") +
  theme(text=element_text(size=20))
















  pivot_wider(names_from=group, values_from=r) %>%
  mutate(r_diff=Asthma-Control) %>%
  group_by(site) %>%
  slice_max(abs(r_diff), n=50) %>%
  ungroup() %>%
  select(site, key, r_diff) %>%
  left_join(diffcor_res %>% filter(group!="all")) %>%
  # mutate(key=str_replace(key, "--", "\n")) %>%
  ggplot(aes(y=key, x=r, colour=group)) +
  geom_point() +
  geom_errorbarh(aes(xmin=`5%`, xmax=`95%`)) +
  facet_wrap(~site, scales="free_y", ncol=1) +
  geom_vline(xintercept=0, colour="red", linetype="dashed")
ggsave("../plots/cf_lll_vs_lul/diffcor_test_plot.pdf", plot=plt)



g_tbl <- booted_cor_adjacencies$cf_lll %>%
  filter(!data.table::between(0, `5%`, `95%`)) %>%
  rename(from=x, to=y, weight=r) %>%
  filter(weight!=0) %>%
  # mutate(
  #   edge_colour=if_else(weight<0, "negative", "positive"),
  #   weight=if_else(abs(weight)>0.4, abs(weight), 0)
  # ) %>%
  as_tbl_graph(directed=FALSE) %>%
  simplify() %>%
  as_tbl_graph() %>%
  left_join(node_info, by="name") 
  
  

# graph layout
node_info <- otu_cluster_labels %>%
  filter(otu %in% names(V(g_tbl))) %>%
  rename(name=otu) %>%
  select(name, cf_lll) %>%
  mutate(node_label=sub("_[^_]+$", "", name))
cluster_sizes <- node_info %>%
  group_by(cf_lll) %>%
  tally()
sum(cluster_sizes$n)

bb <- graphlayouts::layout_as_backbone(g_tbl, keep=0.4)
E(g_tbl)$col <- F
E(g_tbl)$col[bb$backbone] <- T

set_graph_style()
# ggraph(g_tbl, layout="kk") +
ggraph(g_tbl,
       layout = "manual",
       x = bb$xy[, 1],
       y = bb$xy[, 2]) +
  geom_node_point(size=5, aes(colour=as.factor(cf_lll))) +
  geom_edge_link(aes(edge_alpha=weight, colour=edge_colour)) +
  scale_edge_colour_manual(values=c(positive="darkgreen", negative="red")) +
  geom_node_label(aes(label=node_label), alpha=0.6)

# eigengene analysis
library(WGCNA)


mod_eigengenes_comp <- lapply(
  setNames(names(X_subset), names(X_subset)),
  function(site_name) moduleEigengenes(
    X_subset[[ site_name ]],
    colors=otu_cluster_labels[[ site_name ]]
  )
)

proj_eingenes <- lapply(
  mod_eigengenes_comp,
  function(x) x$eigengenes
)

#
# check for association with asthma
library(rstatix)

kruskal_with_effsize <- function(...) {
  kruskal_pvals <- kruskal_test(...) %>% select(-method)
  kruskal_effects <- kruskal_effsize(...) %>% select(-method)
  return(kruskal_pvals %>% inner_join(kruskal_effects, by=c(".y.", "n")))
}

proj_eingenes_long <- lapply(
  proj_eingenes,
    function(x) x %>%
    rownames_to_column("sample_id") %>%
    pivot_longer(-sample_id, names_to="module")) %>%
  bind_rows(.id="site") 

module_disease_assocs <- proj_eingenes_long %>%
  left_join(sample_metadata %>% select(sample_id, Disease)) %>%
  group_by(site, module) %>%
  do(kruskal_with_effsize(data=., formula=value~Disease)) %>%
  ungroup() %>%
  mutate(p_adj=p.adjust(p, method="fdr")) %>%
  arrange(p_adj)

#
# check overlap between sites
tidy_overlap_test <- function(x, y) {
  out <- WGCNA::overlapTable(x, y)
  out$countTable %>% 
    as.data.frame() %>%
    rownames_to_column("varX") %>%
    pivot_longer(-varX, names_to="varY", values_to="count") %>%
    inner_join(
      out$pTable %>% 
        as.data.frame() %>%
        rownames_to_column("varX") %>%
        pivot_longer(-varX, names_to="varY", values_to="pvalue"),
      by=c("varX", "varY")
    ) %>% list()
}

tidy_overlap_test(otu_cluster_labels$cf_lll,
                  otu_cluster_labels$cf_lul)


module_overlap_tbl <- otu_cluster_labels %>%
  select(-otu) %>%
  colpair_map(tidy_overlap_test) %>%
  shave() %>%
  pivot_longer(-term) %>%
  unnest(value) %>%
  drop_na()  %>%
  mutate(p_adj=p.adjust(pvalue, method="fdr"))
module_overlap_tbl %>%
  mutate(facet_title=paste0(term, "--", name),
         label=if_else(
           p_adj<0.1,
           sprintf("%d\n%.2f", count, p_adj),
           NA_character_))  %>%
  ggplot(aes(x=varX, y=varY)) +
  geom_tile(aes(fill=-log10(p_adj))) +
  geom_label(aes(label=label)) +
  facet_wrap(~facet_title, scales="free") +
  scale_fill_gradient(low="white", high="red")
module_overlap_tbl %>%
  filter(p_adj<0.1)
  
#
# check for association between sites
# proj_eingenes_long <- lapply(proj_eingenes,
#        function(x) x %>%
#          rownames_to_column("sample_id") %>%
#          pivot_longer(-sample_id, names_to="module")) %>%
#   bind_rows(.id="site") 
# 
# paired_proj_eigengenes <- proj_eingenes_long %>%
#   left_join(study_ids) %>%
#   select(-sample_id) %>%
#   pivot_wider(names_from=c(site, module), names_sep="____") %>%
#   drop_na() %>%
#   pivot_longer(-Study_ID) %>%
#   separate(name, c("site", "module"), sep="____")
# 
# paired_proj_eigengenes %>%
#   group_by(site) %>%
#   do(mat=pivot_wider(., names_from=module) %>%
#        column_to_rownames("Study_ID") %>%
#        select(-site) %>%
#        as.matrix()) %>%
#   gather_pairs(key=site, value=mat) %>%
#   rowwise() %>%
#   mutate(
#     eigengene_cor=cor(.y, .x) %>% list(),
#     eigengene_cor_pvals=
#   ) -> tmp

