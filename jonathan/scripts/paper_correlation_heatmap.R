library(phyloseq)
library(data.table)
library(tidyr)
library(dplyr)
library(stringr)
library(tibble)
library(dendextend)

library(ggplot2)
library(ggdendro)
library(ggrepel)
library(ComplexHeatmap)
ht_opt$message=FALSE
library(randomcoloR)
library(RColorBrewer)

########################################################################
# utility functions
########################################################################

source("../scripts/kegg_utils.R")
source("../scripts/plot_utils.R")

get_most_abundant <- function(x, n) {
  sum_rel_abund <- colSums(x) %>% sort(decreasing=TRUE)
  return(x[,names(sum_rel_abund)[1:n]])
}

########################################################################
########################################################################

########################################################################
# load required objects
########################################################################

# OTU counts
load("../data/Bus_CF_working_noduplicates.Rdata")

X <- list()
X[["cf_lll"]] <- t(as.matrix(otu_table(subset_samples(Bus_CF1, SwabSite=="Brushing LLL 16S"))))
X[["cf_lul"]] <- t(as.matrix(otu_table(subset_samples(Bus_CF1, SwabSite=="Brushing LUL 16S"))))
X[["cf_ots"]] <- t(as.matrix(otu_table(subset_samples(Bus_CF1, SwabSite=="OTS"))))
X[["bus_ots"]] <- t(as.matrix(otu_table(subset_samples(Bus_CF1, SwabSite=="Throat"))))

module_name_tbl <- fread("../results/wgcna/wgcna_module_names.csv") %>%
  mutate(
    module_longname=str_replace(module_longname, "Family_XIII", "Clostridiales spp.") %>%
      str_replace("Prevotella_\\d+", "Prevotella")
  )
module_name_tbl$module_longname[ module_name_tbl$site=="cf_ots" & module_name_tbl$module=="turquoise" ] <- "turquoise (Heterogeneous\nmodule)"

taxonomy_df <- phyloseq::tax_table(Bus_CF1)@.Data %>%
  as.data.frame() %>%
  rename(OTU=OTUID)

wgcna_modules <- readRDS("../results/wgcna/module_assignments.rds")
wgcna_module_membership <- readRDS("../results/wgcna/module_membership.rds")

otu_trees <- readRDS("../results/wgcna/otu_trees.rds")

########################################################################
########################################################################

########################################################################
# make heatmap
########################################################################

for(SITE_NAME in names(X)[ 3 ]) {
  cat(SITE_NAME, "\n")
  
  x <- X[[ SITE_NAME ]] # plot Spearman correlation of this matrix
  
  # order OTUs by WGCNA modules
  otu_dend <- otu_trees[[  SITE_NAME  ]] %>%
    as.dendrogram()
  ordered_otus <- otu_dend %>% 
    labels()
  
  # or by WGCNA modules
  ordered_otus <- wgcna_modules[[ SITE_NAME ]] %>%
    left_join(
      taxonomy_df %>%
        select(OTU, Phylum),
      by="OTU"
    ) %>%
    mutate(
      module=forcats::fct_relevel(module, "grey", after=Inf)
    ) %>%
    arrange(module, Phylum) %>%
    pull(OTU) 
    
  
  # spearman correlation between the absolute counts
  corr_im <- x[,ordered_otus] %>%
    # get_most_abundant(n_plotted_OTUs) %>%
    cor(method="spearman")
  
  otu_summary <- tibble(
    OTU=colnames(x),
    pc_reads=100*colSums(x)/sum(x),
    prevalence=100*apply(x, 2, function(xx) sum(xx>0)/length(xx))
  ) %>%
    filter(OTU %in% rownames(corr_im)) %>%
    arrange(match(OTU, rownames(corr_im))) %>%
    left_join(wgcna_modules[[ SITE_NAME ]], by="OTU")
  
  # plot annotations: Phylum and WGCNA module
  annotation_df <- taxonomy_df %>%
    select(OTU, Phylum) %>%
    filter(OTU %in% rownames(corr_im)) %>%
    arrange(match(OTU, rownames(corr_im))) %>%
    as_tibble() %>%
    mutate(Phylum2=recode_legend_keys(OTU, Phylum, 10))
  unique_legend_keys <- unique(annotation_df$Phylum2)
  
  # top annotation (taxonomy) colour maps
  annotation_cmap <- setNames(
    scales::hue_pal()(length(unique_legend_keys)-1),#randomcoloR::distinctColorPalette(length(unique_legend_keys)),
    unique_legend_keys[ unique_legend_keys!="Other" ]
  )
  annotation_cmap[[ "Other" ]] <- "grey50"
  module_cmap <- setNames(unique(otu_summary$module), unique(otu_summary$module))
  
  top_annotation <- HeatmapAnnotation(
    Phylum=annotation_df$Phylum2 %>% as.factor() %>% forcats::fct_relevel("Other", after=Inf),
    which="column",
    col=list(Phylum=annotation_cmap),
    show_legend=TRUE,
    annotation_label="Phylum"
  )
  
  # bottom annotatoin - WGCNA modules
  otu_labels <- tibble(
    OTU=ordered_otus) %>%
    left_join(wgcna_modules[[ SITE_NAME ]], by="OTU")
  
  rle_encoding <- rle(otu_labels$module)
  rle_encoding2 <- rle_encoding
  rle_encoding2$values <- make.unique(rle_encoding$values)
  inverse.rle(rle_encoding2)
  
  tmp <- otu_labels %>%
    mutate(new_label=inverse.rle(rle_encoding2)) %>%
    group_by(new_label) %>%
    add_tally() %>%
    mutate(group_idx=row_number()) %>%
    filter(n>5) %>%
    ungroup() %>%
    mutate(x_loc=match(OTU, colnames(corr_im))) %>%
    group_by(new_label) %>%
    summarise(median_x=median(x_loc)) %>%
    mutate(new_label=str_remove(new_label, "\\.\\d+")) %>%
    left_join(module_name_tbl %>% filter(site==SITE_NAME), by=c(new_label="module")) %>%
    mutate(new_label2=str_extract(module_longname, "\\([^()]+\\)") %>%
             str_remove_all("\\(|\\)")) %>%
    filter(!is.na(new_label2))
    
  bottom_annotation <- HeatmapAnnotation(
    "WGCNA module"=otu_summary$module,
    label=anno_mark(at=tmp$median_x %>% as.integer(), labels=paste0("", tmp$new_label2),
                    which="column", side="bottom",
                    labels_gp=gpar(fontsize=10), labels_rot=30, extend=unit(100, "mm")),
    which="column",
    col=list("WGCNA module"=module_cmap [ !is.na(module_cmap) ]),
    show_legend=FALSE
  )
  
  
  # abundance and prevalence sidebars
  barplot_axis_params <- default_axis_param("row")
  barplot_axis_params$direction <- "reverse"
  
  barplot_axis_params_pc_reads <- barplot_axis_params
  barplot_axis_params_pc_reads$at <- seq(1, -5, by=-2)
  barplot_axis_params_pc_reads$labels <- c(
    # expression(10^0), expression(10^-2), expression(10^-4), expression(10^-6)
    expression(10^1), expression(10^-1), expression(10^-3), expression(10^-5)
    )# 10**seq(0, -8, by=-2)
  
  barplot_axis_params_prev <- barplot_axis_params
  barplot_axis_params_prev$at <- c(100, 50, 0)
  barplot_axis_params_prev$labels <- paste0(c(100, 50, 0), "%")
  
  BARPLOT_WIDTH <- unit(2, "cm")
  
  left_annotation <- HeatmapAnnotation(
    log10_pc_reads=anno_barplot(
      log10(otu_summary$pc_reads), axis_param=barplot_axis_params_pc_reads,
      width=BARPLOT_WIDTH
      ),
    prevalence=anno_barplot(otu_summary$prevalence, axis_param=barplot_axis_params_prev,
                            width=BARPLOT_WIDTH),
    gap=unit(1.5, "mm"),
    which="row",
    show_legend=TRUE,
    annotation_name_side="top", annotation_name_rot=30,
    annotation_name_align=TRUE
  )
  
  # make heatamp
  dendrogram_size <- 10
  
  ht_opt$ROW_ANNO_PADDING <- unit(1.5, "mm")
  
  file_name <- sprintf(
    "../plots/wgcna/annotated_corr_heatmap__%s.pdf", SITE_NAME
  )
    
  
  pdf(file_name, width=11, height=8)
  hm <- Heatmap(
    corr_im,
    col=circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
    name="Spearman\ncorrelation",
    show_row_dend=FALSE,
    show_column_dend=FALSE,
    row_dend_width=unit(dendrogram_size, "mm"),
    column_dend_height=unit(dendrogram_size, "mm"),
    cluster_rows=FALSE,#otu_trees$cf_ots,
    cluster_columns=FALSE,#otu_trees$cf_ots,
    top_annotation=top_annotation,
    left_annotation=left_annotation,
    bottom_annotation=bottom_annotation,
    show_row_names=FALSE,
    show_column_names=FALSE,
    heatmap_legend_param=list(title="Spearman\ncorrelation", grid_border="black"),
    heatmap_width=unit(18, "cm"), heatmap_height=unit(18, "cm")
  )
  draw(hm, merge_legend=TRUE)
  dev.off()
  
  knitr::plot_crop(file_name)
}


########################################################################
########################################################################

########################################################################
# annotated isolate dendrogram
########################################################################

# OTU -> isolate map
otu_isolate_map <- fread("../data/isolate_names/otu_isolate_map.csv")
sapply(X, function(x) all(otu_isolate_map$OTU %in% colnames(x)))

isolate_hclust <- readRDS("../results/kegg/isolate_phylogenetic_hclust.rds")
cluster_tbl <- fread("../results/kegg/isolate_cluster_labels.csv")
cluster_tbl$cluster[cluster_tbl$isolate=="Cupriavidus gilardii 27098_8_120"] <- NA
cluster_tbl$cluster_longname[cluster_tbl$isolate=="Cupriavidus gilardii 27098_8_120"] <- NA

#
# isolate dendrogram (clustered on KO distances)
tree_data <- isolate_hclust %>%
  as.dendrogram() %>%
  dendro_data(type="rectangle")

tree_data$labels <- tree_data$labels %>%
  left_join(cluster_tbl, by=c(label="isolate")) %>%
  # mutate(disp_label=str_trunc(label, 30))
  mutate(disp_label=gsub('[[:digit:]]+', '', label) %>% str_remove_all("_") %>% trimws("right"))
isolate_leaf_order <- tree_data$labels$label

colourbar_data <- expand_grid(
  y=c(-10,-20), x=min(tree_data$labels$x):max(tree_data$labels$x)
) %>%
  left_join(
    tree_data$labels %>%
      select(x, label),
    by="x") %>%
  left_join(cluster_tbl, by=c(label="isolate")) %>%
  mutate(
    cluster_integer=factor(cluster_longname,
                           levels=unique(cluster_longname)) %>%
      as.integer(),
    fill_colour=as.factor(cluster_integer %% 2)
  )

cluster_text_labels <- colourbar_data %>%
  select(x, label, cluster_longname) %>%
  group_by(cluster_longname) %>%
  summarise(x=mean(x)) # %>%
# mutate(cluster_longname=str_replace(cluster_longname, "\\+", "\n"))

tree_panel <- ggplot() +
  geom_segment(data=tree_data$segments, aes(x=x, y=y, xend=xend, yend=yend))+
  # geom_text(data=tree_data$labels,
  #           aes(x=x, y=y-20, label=disp_label, hjust=0),
  #           colour="black",
  #           key_glyph=draw_key_label,
  #           size=2, angle=90, hjust=1) +
  theme_dendro() +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank()) +
  labs(colour="Cluster") +
  coord_cartesian(ylim=c(-350,NA)) # + # xlim=c(6.5, 121.5), 
tree_panel

# colour palette for isolate clusters
set.seed(123456)
n <- 17
palette <- distinctColorPalette(n)
pie(rep(1, n), col=palette)

cluster_colours <- palette[-11]

rect_tops <- c(
  "Strep.I"=360, "Strep.II"=320, "Strep.III"=200, "Granulicatella"=150, "Gemella"=240,
  "Staphy.+Nialla"=250, "Micrococcus"=100, "Cutibacterium"=50,
  "Rothia"=180, "Pauljensenia"=150, "Actinomyces"=150,
  "Haemophilus"=370, "Prevotella"=370,
  "Veillonella"=200, "Fusobacterium+Leptotrichia"=630,
  "Neisseria"=150
)

label_y_pos <- c(
  "Strep.I"=450, "Strep.II"=450, "Strep.III"=450, "Granulicatella"=450, "Gemella"=450,
  "Staphy.+Nialla"=450, "Micrococcus"=450, "Cutibacterium"=450,
  "Rothia"=450, "Pauljensenia"=450, "Actinomyces"=450,
  "Haemophilus"=450, "Prevotella"=450,
  "Veillonella"=300, "Fusobacterium+Leptotrichia"=500,
  "Neisseria"=450
)

rect_locations <- colourbar_data %>%
  select(x, label, cluster_longname) %>%
  group_by(cluster_longname) %>%
  summarise(x_min=min(x)-0.05, x_max=max(x)+0.05) %>%
  mutate(y_min=0) %>%
  left_join(
    tibble(cluster_longname=names(rect_tops), y_max=rect_tops),
    by="cluster_longname"
  ) %>%
  left_join(
    tibble(cluster_longname=names(label_y_pos), label_y=label_y_pos),
    by="cluster_longname"
  ) %>%
  mutate(label_x=(x_max+x_min)/2,
         label_text=str_replace(cluster_longname, "\\+", "\n"))

tree_panel_with_rects <- tree_panel +
  geom_rect(data=rect_locations,
            aes(xmin=x_min, xmax=x_max, ymin=y_min, ymax=y_max, fill=cluster_longname),
            alpha=0.4) +
  geom_label_repel(data=rect_locations,
                   aes(x=label_x, y=label_y, label=label_text, colour=cluster_longname),
                   size=5, fill = alpha(c("white"),0.5)) +
  scale_colour_manual(values=cluster_colours, guide='none') +
  scale_fill_manual(values=cluster_colours)
tree_panel_with_rects +
  theme(legend.position="none")
cropped_ggsave("../plots/070722_meeting_plots/isolate_dendrogram.pdf",
               height=6, width=13)

#
# OTU-isolate percentage identity heatmap
get_fill_label <- function(x) {
  out <- rep(NA, length(x))
  out[ x>= 97 | is.na(x) ] <- ">=97%"
  out[ x>=99 ] <- ">=99%" 
  out[ x==100 ] <- "=100%"
  out[ x<97  | is.na(x) ] <- "<97%"
  return(factor(out, levels=c("=100%", ">=99%", ">=97%", "<97%")))
}

# only plot 100 most abundant OTUs from CF OTS samples
otu_abundances <- tibble(
  OTU=colnames(X$cf_ots),
  total_reads=colSums(X$cf_ots),
  total_relative_abundance=colSums(X$cf_ots/rowSums(X$cf_ots))
)
plotted_otus <- otu_abundances %>%
  slice_max(total_relative_abundance, n=50) %>%
  pull(OTU)

plotted_otus <- otu_isolate_map %>% 
  left_join(otu_abundances) %>%
  filter(isolate!="no_isolate") %>%
  select(OTU, total_reads) %>%
  distinct() %>%
  arrange(desc(total_reads)) %>%
  head(50) %>%
  pull(OTU)
  

pc_identity_heatmap_data <- otu_isolate_map %>%
  filter(isolate!="no_isolate") %>%
  mutate(isolate=factor(isolate, levels=isolate_leaf_order)) %>%
  complete(isolate) %>% 
  filter(!is.na(isolate)) %>%
  complete(isolate, OTU) %>%
  mutate(fill_label=get_fill_label(percentage_identity)) %>%
  filter(OTU %in% plotted_otus)

stopifnot(length(unique(pc_identity_heatmap_data$isolate))==126)

colours <- colorRampPalette(brewer.pal(9,"Blues"))(15) %>% tail(10) %>% rev()
colours <- colours[c(1, length(colours)/2, length(colours))]

cluster_vline_tbl <- cluster_tbl %>%
  mutate(isolate=factor(isolate, levels=isolate_leaf_order)) %>%
  arrange(isolate) %>%
  mutate(x=row_number())
tmpcol <- cluster_vline_tbl$cluster_longname %>% tail(-1) != cluster_vline_tbl$cluster_longname %>% head(-1)
cluster_vline_tbl <- cluster_vline_tbl %>% mutate(
  show_vline=c(FALSE, tmpcol)
)


isolate_otu_identity_heatmap <- pc_identity_heatmap_data %>%
  replace_na(list(fill_value="<97%")) %>%
  filter(!is.na(isolate) | is.na(OTU)) %>%
  ggplot(aes(y=OTU, x=isolate, fill=fill_label)) +
  geom_tile() +
  theme(axis.text.x=element_text(size=8, angle=45, hjust=1),#element_blank(),
        axis.text.y=element_text(size=6.5),
        axis.title=element_text(size=16),
        axis.title.x=element_text(vjust=7),
        panel.grid=element_blank(),
        legend.text=element_text(size=14),
        legend.box.background = element_rect(colour="black"),
        strip.placement="outside",
        strip.background=element_blank(),
        strip.text.x=element_text(angle=90, hjust=1)) +
  scale_fill_manual(values=c(colours, "white"), na.value="white") +
  labs(fill="16S rRNA sequence identity", x="Isolate") +
  scale_x_discrete(labels=function(x) str_remove_all(x, "\\d+|_")) +
  geom_vline(data=cluster_vline_tbl %>% filter(show_vline),
             aes(xintercept=x-0.5),
             linetype="dotted")

isolate_otu_identity_heatmap 




hm_grob <- grid::grid.grabExpr(draw(hm, merge_legend=TRUE))

#
# abundance and prevalence comparison plots
otu_summary_all_sites <- lapply(
  X, 
  function(x) tibble(
    OTU=colnames(x),
    total_reads=colSums(x),
    pc_reads=100*colSums(x)/sum(x),
    prevalence=100*apply(x, 2, function(xx) sum(xx>0)/length(xx))
  )
) %>%
  bind_rows(.id="site") %>%
  left_join(taxonomy_df, by="OTU") 

annotation_cmap

plotted_phyla <- names(annotation_cmap) %>% head(-1)

otu_summary_all_sites  <- otu_summary_all_sites %>%
  mutate(fill_colour=if_else(Phylum %in% plotted_phyla, Phylum, "Other"))

POINT_SIZE <- 1

otu_summary_all_sites

plt1 <- otu_summary_all_sites %>%
  filter(grepl("ots", site)) %>%
  select(site, OTU, total_reads, fill_colour) %>%
  pivot_wider(names_from=site, values_from=total_reads) %>%
  filter(cf_ots>0, bus_ots>0) %>%
  ggplot(aes(x=cf_ots, y=bus_ots)) +
  geom_point(aes(colour=fill_colour), size=POINT_SIZE) +
  scale_x_continuous(trans="log10", limits=c(c(1e-0, 1e7))) +
  scale_y_continuous(trans="log10", limits=c(c(1e-0, 1e7))) +
  xlab("Reads in CELF\nptOP samples") +
  ylab("Reads in BUS\nptOP samples") +
  theme_linedraw() +
  theme(aspect.ratio=1,
        axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  geom_abline() +
  scale_colour_manual(values=annotation_cmap) +
  labs(colour="Phylum") +
  guides(color=guide_legend(override.aes=list(size=5)))
  

plt2 <- otu_summary_all_sites %>%
  filter(grepl("cf_ots|cf_lll", site)) %>%
  select(site, OTU, total_reads, fill_colour) %>%
  pivot_wider(names_from=site, values_from=total_reads) %>%
  filter(cf_ots>0, cf_lll>0) %>%
  ggplot(aes(x=cf_ots, y=cf_lll)) +
  geom_point(aes(colour=fill_colour), size=POINT_SIZE) +
  scale_x_continuous(trans="log10", limits=c(c(1e-0, 1e7))) +
  scale_y_continuous(trans="log10", limits=c(c(1e-0, 1e7))) +
  xlab("Reads in CELF\nptOP samples") +
  ylab("Reads in CELF\nLLL samples") +
  theme_linedraw() +
  theme(aspect.ratio=1,
        axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  geom_abline()  +
  scale_colour_manual(values=annotation_cmap) +
  labs(colour="Phylum") +
  guides(color=guide_legend(override.aes=list(size=5)))


plt3 <- otu_summary_all_sites %>%
  filter(grepl("ots", site)) %>%
  select(site, OTU, prevalence, fill_colour) %>%
  pivot_wider(names_from=site, values_from=prevalence) %>%
  filter(cf_ots>0, bus_ots>0) %>%
  ggplot(aes(x=cf_ots, y=bus_ots)) +
  geom_point(aes(colour=fill_colour), size=POINT_SIZE) +
  scale_x_continuous(labels=function(x) scales::percent(x, scale=1), limits=c(c(0, 100))) +
  scale_y_continuous(labels=function(x) scales::percent(x, scale=1), limits=c(c(0, 100))) +
  xlab("Prevalence in\nCELF ptOP samples") +
  ylab("Prevalence in\nBUS ptOP samples") +
  theme_linedraw() +
  theme(aspect.ratio=1,
        axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  geom_abline() +
  scale_colour_manual(values=annotation_cmap) +
  labs(colour="Phylum") +
  guides(color=guide_legend(override.aes=list(size=5)))

plt4 <- otu_summary_all_sites %>%
  filter(grepl("cf_ots|cf_lll", site)) %>%
  select(site, OTU, prevalence, fill_colour) %>%
  pivot_wider(names_from=site, values_from=prevalence) %>%
  filter(cf_ots>0, cf_lll>0) %>%
  ggplot(aes(x=cf_ots, y=cf_lll)) +
  geom_point(aes(colour=fill_colour), size=POINT_SIZE) +
  scale_x_continuous(labels=function(x) scales::percent(x, scale=1), limits=c(c(0, 100))) +
  scale_y_continuous(labels=function(x) scales::percent(x, scale=1), limits=c(c(0, 100))) +
  xlab("Prevalence in\nCELF ptOP samples") +
  ylab("Prevalence in\nCELF LLL samples") +
  theme_linedraw() +
  theme(aspect.ratio=1,
        axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  geom_abline() +
  scale_colour_manual(values=annotation_cmap) +
  labs(colour="Phylum") +
  guides(color=guide_legend(override.aes=list(size=5)))
plt4


#
# comparing Strep KO scores
strep_mm_plotdata <- wgcna_module_membership$cf_ots %>%
  right_join(otu_isolate_map %>%
               filter(percentage_identity==100),
             by="OTU") %>%
  filter(grepl("blue|greenyellow|red", module), grepl("Streptococcus", OTU)) %>%
  left_join(cluster_tbl, by="isolate") %>%
  mutate(module=str_remove(module, "ME"))

strep_mm_plotdata %>%
  group_by(module, cluster_longname) %>%
  summarise(correlation=mean(correlation)) %>%
  ggplot(aes(x=cluster_longname, y=correlation)) +
  geom_col(position="dodge") +
  ylab("Mean module\nmembership") +
  xlab("WGCNA module") +
  labs(fill="Isolate cluster") +
  theme_bw() +
  theme(
    axis.text=element_text(size=12),
    axis.title=element_text(size=16),
    legend.title=element_text(size=16),
    legend.text=element_text(size=12)
  ) +
  lemon::facet_rep_wrap(~module)

strep_mm_plotdata %>%
  ggplot(aes(x=cluster_longname, y=correlation)) +
  geom_boxplot() +
  ylab("Mean module\nmembership") +
  xlab("WGCNA module") +
  labs(fill="Isolate cluster") +
  theme_bw() +
  theme(
    axis.text=element_text(size=12),
    axis.title=element_text(size=16),
    legend.title=element_text(size=16),
    legend.text=element_text(size=12)
  ) +
  lemon::facet_rep_wrap(~module)

LABEL_SIZE <- 28

bc_legend <- cowplot::get_legend(plt1)

#
# create composite figure and save
fig2 <- egg::ggarrange(
  tree_panel_with_rects +
    coord_cartesian(xlim=c(6.25, 121)) +
    theme(legend.position="none"),
  ggplot() + theme_void(),
  isolate_otu_identity_heatmap +
    theme(plot.margin=margin(-2, 0, 0, 0, "pt"),
          legend.position=c(0.15, 0.3)) +
    guides(fill=guide_legend(ncol=2)),
  ncol=1,
  heights=c(1.2, -0.22, 2.5)
)

fig3 <- cowplot::plot_grid(
  fig2,
  cowplot::plot_grid(
    cowplot::plot_grid(
      plotlist=lapply(
        list(plt1, plt2, plt3, plt4),
        function(x) x + theme(legend.position="none", panel.grid=element_line(size=0.01, colour="grey50"))
        ),
      ncol=2, nrow=2,
      align="vh",
      axis = "b",
      labels=c("b", "", "c", ""),
      label_size=LABEL_SIZE
    ),
    # cowplot::plot_grid(
    #   bc_legend, bc_legend,
    #   ncol=1
    # ),
    hm_grob,
    nrow=1,
    rel_widths=c(1, 1.5),
    labels=c("",  "d"),
    label_size=LABEL_SIZE
  ),
  ncol=1,
  labels=c("a"),
  label_size=LABEL_SIZE
)

ggsave(
  plot=fig3, filename="../plots/composite_figures/wgcna_and_isolate_dend.pdf",
  height=15, width=16
)
