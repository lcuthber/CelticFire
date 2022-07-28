library(dplyr)
library(tibble)
library(tidyr)
library(data.table)
library(corrr)
library(stringr)
library(forcats)

library(ggplot2)
library(ggforce)
library(ggpubr)

source("../scripts/plot_utils.R")

###############################################################################
# Module overlap plot
###############################################################################

wgcna_modules <- readRDS("../results/wgcna/module_assignments.rds")

module_name_tbl <- fread("../results/wgcna/wgcna_module_names.csv") %>%
  mutate(
    module_longname=str_replace(module_longname, "Family_XIII", "Clostridiales spp.") %>%
      str_replace("Prevotella_\\d+", "Prevotella")
  )
module_name_tbl$module_longname[ module_name_tbl$site=="cf_ots" & module_name_tbl$module=="turquoise" ] <- "turquoise (Heterogeneous\nmodule)"


module_assignments_wide <- wgcna_modules %>%
  rbindlist(idcol="label") %>%
  pivot_wider(names_from=label, values_from=module) %>%
  column_to_rownames("OTU") 

colname_pairs <- apply(
  combn(colnames(module_assignments_wide), 2), 2, function(x) x, simplify=FALSE
)
colname_pairs <- combn(colnames(module_assignments_wide), 2)

module_overlap <- purrr::map2(colname_pairs[1,],
               colname_pairs[2,],
               function(x,y) WGCNA::overlapTable(module_assignments_wide[[ x ]], module_assignments_wide[[ y ]])
)

module_overlap_pvalues <- lapply(
  module_overlap,
  function(x) x$pTable %>%
    as.data.frame() %>%
    rownames_to_column("site1_module") %>%
    pivot_longer(-site1_module, names_to="site2_module", values_to="pvalue")
) %>%
  as_tibble_col() %>%
  bind_cols(t(colname_pairs)) %>%
  magrittr::set_colnames(c("value", "site1", "site2")) %>%
  unnest(value) %>%
  mutate(pvalue=p.adjust(pvalue, method="fdr"))

overlap_axis_titles <- c(
  bus_ots="BUS ptOP",
  cf_ots="CELF ptOP",
  cf_lll="CELF LLL"
)

module_overlap_pvalues %>%
  filter(site1 %in% c("cf_ots", "bus_ots", "cf_lll"),
         site2 %in% c("cf_ots", "bus_ots", "cf_lll")) %>%
  left_join(
    module_name_tbl, by=c(site1_module="module", site1="site")
  ) %>%
  left_join(
    module_name_tbl, by=c(site2_module="module", site2="site")
  ) %>%
  mutate(pval_label=if_else(pvalue<0.05, "*", NA_character_)) %>%
  group_by(site1, site2) %>%
  do(
    plot=ggplot(data=mutate(., 
                            module_longname.x=fct_relevel(module_longname.x,
                                                          "grey", 
                                                          after=Inf),
                            module_longname.y=fct_relevel(module_longname.y,
                                                          "grey", 
                                                          after=0)),
                aes(x=module_longname.x,
                    y=module_longname.y,
                    fill=-log10(pvalue))) +
      geom_tile() +
      geom_text(aes(label=pval_label), size=5) +
      scale_fill_gradient(low="white", high="red", limits=c(0, 55),
                          guide=guide_colorbar(
                            barwidth=10, ticks.colour="black",
                            frame.colour="black"
                          )) +
      theme_minimal() +
      theme(
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11, angle=45, hjust=1, vjust=1),
        axis.title=element_text(size=16),
        strip.text.x=element_text(size=20),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14),
        panel.border=element_rect(fill=NA),
        axis.ticks=element_line()
      ) +
      coord_fixed() +
      xlab(overlap_axis_titles[[ .$site1[[1]] ]]) +
      ylab(overlap_axis_titles[[ .$site2[[1]] ]]) +
      scale_y_discrete(labels=function(x) str_replace(x, " ", "\n")) +
      scale_x_discrete(labels=function(x) str_replace_all(x, "\n", " "))
  ) %>%
  pull(plot) %>%
  rev() %>%
  ggarrange(plotlist=., nrow=1, common.legend=TRUE,
            widths=c(0.9, 1, 2.5) %>% rev(),
            labels="auto")
cropped_ggsave("../plots/supplementary_figures/wgcna_module_overlap.png",
               width=12, height=8,
               dpi=450)

  
  
# ggplot(aes(x=module_longname.x, y=module_longname.y, fill=-log10(pvalue))) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="red", limits=c(0, 55)) +
#   theme_minimal() +
#   theme(
#     axis.text.y=element_text(size=14),
#     axis.text.x=element_text(size=14, angle=45, hjust=1, vjust=1),
#     axis.title=element_text(size=16),
#     strip.text.x=element_text(size=20),
#     legend.text=element_text(size=14),
#     legend.title=element_text(size=14),
#     panel.border=element_rect(fill=NA),
#     axis.ticks=element_line()
#   ) +
#   ggforce::facet_row(~site2) +
#   coord_fixed() +
#   xlab("CELF ptOP module") +
#   ylab("BUS ptOP module") +
#   scale_x_discrete(labels=function(x) str_replace(x, "\n", " "))
# cropped_ggsave("../plots/070722_meeting_plots/ptop_overlap_pvalues.pdf")
# 
# module_overlap_pvalues %>%
#   filter(site2=="cf_ots", site1!="bus_ots") %>%
#   # mutate(
#   #   site1_module=recode(site1_module, greenyellow="green\nyellow")
#   #   ) %>%
#   left_join(
#     module_name_tbl, by=c(site1_module="module", site1="site")
#   ) %>%
#   left_join(
#     module_name_tbl, by=c(site2_module="module", site2="site")
#   ) %>%
#   mutate(site1=recode(site1, cf_lll="CELF LLL\nmodule", cf_lul="CELF LUL\nmodule")) %>%
#   filter(!grepl("LLL", site1)) %>%
#   ggplot(aes(x=module_longname.x, y=module_longname.y, fill=-log10(pvalue))) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="red", limits=c(0, 55)) +
#   theme_minimal() +
#   theme(
#     panel.spacing.x=unit(20, "pt"),
#     axis.text.y=element_text(size=14),
#     axis.text.x=element_text(size=14, angle=45, hjust=1, vjust=1),
#     axis.title=element_text(size=16),
#     strip.text.x=element_text(size=16),
#     legend.text=element_text(size=14),
#     legend.title=element_text(size=14),
#     panel.border=element_rect(fill=NA),
#     axis.ticks=element_line(),
#     strip.placement="outside"
#   ) +
#   ylab("CELF ptOP module") +
#   xlab("") +
#   ggforce::facet_row(~site1, scales="free_x", space="free", strip.position="bottom") +
#   scale_y_discrete(labels=function(x) str_replace(x, "\n", " "))
# cropped_ggsave("../plots/070722_meeting_plots/celf_overlap_pvalues.pdf")

###########################################################################
# Strep isolate/WGCNA module comparison
###########################################################################

get_fill_label <- function(x) {
  out <- rep(NA, length(x))
  out[ x>= 97 | is.na(x) ] <- ">=97%"
  out[ x>=99 ] <- ">=99%" 
  out[ x==100 ] <- "=100%"
  out[ x<97  | is.na(x) ] <- "<97%"
  return(factor(out, levels=c("=100%", ">=99%", ">=97%", "<97%")))
}
otu_isolate_map <- fread("../data/isolate_names/otu_isolate_map.csv")
cluster_tbl <- fread("../results/kegg/isolate_cluster_labels.csv")
cluster_tbl$cluster[cluster_tbl$isolate=="Cupriavidus gilardii 27098_8_120"] <- NA
cluster_tbl$cluster_longname[cluster_tbl$isolate=="Cupriavidus gilardii 27098_8_120"] <- NA

#
# Strep module membership - now in Supplementals
wgcna_modules %>%
  bind_rows(.id="site") %>%
  left_join(
    module_name_tbl, by=c("site", "module")
  ) %>%
  as_tibble() %>%
  right_join(
    otu_isolate_map %>%
      filter(percentage_identity>=99),
    by="OTU") %>%
  left_join(cluster_tbl, by="isolate") %>%
  filter(grepl("Strep", cluster_longname), grepl("ots", site)) %>%
  group_by(site) %>%
  do(tbl=table(.$module_longname, .$cluster_longname) %>% as.data.frame()) %>%
  unnest(tbl) %>%
  mutate(site=recode(site, bus_ots="BUS ptOP", cf_ots="CELF ptOP"),
         Var1=fct_relevel(
           Var1,
           c("blue (Streptococcus)", "red (Streptococcus)",
             "greenyellow (Streptococcus)", "yellow (Streptococcus)"), after=0) %>%
           fct_relevel()) %>%
  ggplot(aes(x=Var1, y=Var2, fill=Freq)) +
  geom_tile() +
  geom_label(aes(label=Freq)) +
  scale_fill_gradient(low="white", high="red",
                      guide=guide_colorbar(
                        barheight=10, ticks.colour="black",
                        frame.colour="black"
                      )) +
  xlab("WGCNA module") +
  ylab("Isolate cluster") +
  theme_minimal() +
  theme(
    axis.text=element_text(size=12),
    axis.text.x=element_text(size=12, angle=30, hjust=1, vjust=1),
    axis.title=element_text(size=16),
    strip.text.x=element_text(size=20),
    legend.text=element_text(size=14),
    legend.title=element_text(size=14),
    panel.border=element_rect(fill=NA)
  ) +
  ggforce::facet_row(~site, scales="free_x", space="free") +
  labs(fill="# OTUs") +
  scale_x_discrete(labels=function(x) str_replace(x, " ", "\n"))
cropped_ggsave("../plots/composite_figures/strep_wgcna_isolate_comp.pdf")

###########################################################################
###########################################################################

###########################################################################
# triple-decked heatmap: failed attempt
###########################################################################

X <- readRDS("../data/wgcna/otu_tables_wgcna.rds")

strep_counts <- lapply(
  X,
  function(x) x[,grepl("Streptococcus", colnames(x))] %>%
    as.data.frame() %>%
    rownames_to_column("sample_id") %>%
    pivot_longer(-sample_id, names_to="OTU", values_to="reads")) %>%
  bind_rows(.id="site") %>%
  filter(grepl("ots", site))

otus_w_isolate_clusters <- otu_isolate_map %>%
  filter(OTU %in% strep_counts$OTU, grepl("Streptococcus", isolate)) %>%
  left_join(cluster_tbl) %>%
  filter(percentage_identity>99.03) 

# %>%
#   group_by(OTU) %>%
#   summarise(n=length(unique(cluster_longname))) %>%
#   do(table(.$n) %>% as.data.frame())

bind_rows(
  strep_counts %>%
    rename(facet_row_title=site, y=sample_id, fill_value=reads, x=OTU),
  pc_identity_heatmap_data %>%
    mutate(facet_row_title="isolates") %>%
    rename(x=isolate, y=OTU, fill_value=percentage_identity)
) 

# order samples by MM in strep WGCNA modules
module_eigenvecs <- readRDS("../results/wgcna/module_eigenvectors.rds") 
sample_strep_modules <- lapply(module_eigenvecs,
                               function(x) x %>%
                                 rownames_to_column("sample_id") %>%
                                 pivot_longer(-sample_id, names_to="module", values_to="mm")) %>%
  bind_rows(.id="site") %>%
  filter(grepl("ots", site)) %>%
  right_join(
    tribble(
      ~site, ~module,
      "bus_ots", "MEblue",
      "bus_ots", "MEyellow",
      "cf_ots", "MEblue",
      "cf_ots", "MEgreenyellow",
      "cf_ots", "MEred"
    ),
    by=c("site", "module")
  ) %>%
  group_by(site, sample_id) %>%
  summarise(max_mm_module=module[ which.max(mm)],
            max_mm_val=mm[ which.max(mm)]) %>%
  arrange(desc(max_mm_val))
  
strep_counts %>%
  left_join(otus_w_isolate_clusters, by="OTU") %>%
  left_join(sample_strep_modules) %>%
  mutate(facet_row_title=sprintf("%s__%s", site, max_mm_module)) %>%
  ggplot(aes(x=OTU, y=sample_id, fill=log(reads+1))) +
  geom_tile() +
  # facet_col(~facet_row_title, space="free", scales="free_y",
  #           strip.position="left") +
  facet_grid(
    rows=vars(facet_row_title), cols=vars(cluster_longname),
    scales="free", space="free"
  ) +
  scale_fill_gradient(low="white", high="red") +
  theme_linedraw() +
  theme(
    axis.text.y=element_text(size=4),
    axis.text.x=element_text(size=4, angle=30, hjust=1, vjust=1),
    strip.background=element_blank(),
    strip.text=element_text(colour="black"),
    strip.placement="outside"
  )

mapping_duplicate_res <- otu_isolate_map %>%
  filter(isolate!="no_isolate") %>%
  pivot_wider(names_from=isolate, values_from=percentage_identity) %>%
  column_to_rownames("OTU") %>%
  group_duplicate_columns()
unique_otu_isolate_map <- mapping_duplicate_res$new_df %>%
  rownames_to_column("OTU") %>%
  pivot_longer(-OTU, names_to="isolate", values_to="percentage_identity") %>%
  filter(!is.na(percentage_identity))

wgcna_mm <- readRDS("../results/wgcna/module_membership.rds")

wgcna_mm %>%
  bind_rows(.id="site") %>%
  filter(grepl("ots", site)) %>%
  right_join(
    unique_otu_isolate_map %>%
      filter(grepl("Streptococcus", isolate)) %>%
      left_join(cluster_tbl, by="isolate"),
    by="OTU"
  ) %>%
  right_join(
    tribble(
      ~site, ~module,
      "bus_ots", "MEblue",
      "bus_ots", "MEyellow",
      "cf_ots", "MEblue",
      "cf_ots", "MEgreenyellow",
      "cf_ots", "MEred"
    ),
    by=c("site", "module")
  ) %>%
  mutate(binned_pc_ident=factor(
    get_fill_label(percentage_identity),
    levels=c("<97%", ">=97%", ">=99%", "=100%")
    )) %>%
  filter(site=="cf_ots") %>%
  ggplot(aes(y=correlation, x=binned_pc_ident, colour=module)) +
  geom_point(position=position_dodge(0.9)) +
  facet_wrap(~cluster_longname)
###########################################################################
###########################################################################

###########################################################################
# emapper visualistion (non-dupcalites only)
###########################################################################

emapper <- fread("../data/kegg/formatted/binary_emapper_matrix_no_duplicates.csv") %>%
  column_to_rownames("isolate")
duplicated_KOs <- fread("../data/kegg/formatted/duplicated_KOs.csv") %>%
  mutate(identical_to=str_split(identical_to, "\\|")) %>%
  unnest(identical_to)

isolate_hclust <- readRDS("../results/kegg/isolate_phylogenetic_hclust.rds")

isolate_hclust %>%
  as.dendrogram() %>%
  plot()

isolate_leaf_order <- isolate_hclust %>% labels()
ko_order <- colSums(emapper) %>% sort(decreasing=TRUE) %>% names()

emapper %>% 
  rownames_to_column("isolate") %>%
  pivot_longer(-isolate) %>%
  mutate(isolate=factor(isolate, levels=isolate_leaf_order),
         name=factor(name, levels=ko_order)) %>%
  ggplot(aes(x=name, y=isolate, fill=as.factor(value))) +
  geom_tile() +
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.text.y=element_text(size=8),
    axis.title=element_text(size=20),
    legend.position="none",
    panel.border=element_rect(fill=NA)
  ) +
  scale_fill_manual(values=c("white", "red")) +
  xlab("KO") + ylab("Isolate") +
  scale_y_discrete(labels=function(x) str_remove_all(x, "\\d+|_"))
cropped_ggsave("../plots/070722_meeting_plots/emapper_vis.pdf", width=14, height=10)

###########################################################################
###########################################################################

###########################################################################
# KO plots
###########################################################################

cluster_ko_scores <- readRDS("../results/kegg/cluster_ko_odds_ratios.rds") %>%
  bind_rows(.id="cluster") %>%
  mutate(score=log10(Odds_ratio)) %>%
  left_join(cluster_names %>% mutate(cluster=paste0("cluster_", cluster)),
            by="cluster")

cluster_tbl <- fread("../results/kegg/isolate_cluster_tbl.csv") 

cluster_ko_scores <- cluster_ko_scores %>%
  left_join(duplicated_KOs %>%
              rename(KO=representative_ko, excluded_KO=identical_to) %>%
              bind_rows(
                tibble(KO=unique(.$KO), excluded_KO=unique(.$KO))
              ),
            by="KO") %>%
  rename(base_KO=KO, KO=excluded_KO) %>%
  relocate(KO, base_KO) %>%
  mutate(KO=coalesce(KO, base_KO)) 

low_score_KOs <- cluster_ko_scores %>%
  filter(cluster_longname=="Strep.III") %>%
  slice_min(order_by=score, n=10) %>%
  mutate(label="low")

high_score_KOs <- cluster_ko_scores %>%
  filter(cluster_longname=="Strep.III") %>%
  slice_max(order_by=score, n=10) %>%
  mutate(label="high")

emapper %>% 
  rownames_to_column("isolate") %>%
  pivot_longer(-isolate, names_to="KO") %>%
  inner_join(bind_rows(high_score_KOs, low_score_KOs) %>%
               select(-c(cluster, cluster_longname)),
             by="KO") %>%
  left_join(cluster_tbl, by="isolate") %>%
  mutate(
    cluster_longname=str_replace(cluster_longname, "\\+", ",\n"),
    label=recode(label, high="Over\nrepresented", low="Under\nrepresented")
    ) %>%
  ggplot(aes(x=isolate, y=KO, fill=as.factor(value))) +
  geom_tile() +
  scale_fill_manual(values=c("white", "red")) +
  ylab("KO") + xlab("Isolate") +
  theme_minimal() +
  facet_grid(rows=vars(label), cols=vars(cluster_longname),
             space="free", scales="free") +
  theme(
    axis.text.x=element_text(size=8, angle=90, vjust=1, hjust=1),
    axis.text.y=element_text(size=8),
    strip.text.x=element_text(size=10, angle=90, vjust=0.5, hjust=0),
    strip.text.y=element_text(size=10),
    panel.spacing.x=unit(10, "pt"),
    panel.spacing.y=unit(20, "pt"),
    legend.position="none",
    panel.border=element_rect(fill=NA)
  ) +
  scale_x_discrete(labels=function(x) str_remove_all(x, "\\d+|_"))
cropped_ggsave("../plots/070722_meeting_plots/strep_iii_KOs.pdf", width=14, height=7)

strep_ko_comp_plotdata <- cluster_ko_scores %>%
  filter(grepl("Strep", cluster_longname)) %>%
  tidybayes::gather_pairs(key=cluster_longname, value=score) %>%
  mutate(score_diff=.x-.y,
         abs_score_diff=abs(score_diff),
         colour_label=if_else(abs_score_diff>1, "different", "not_different"))

strep_ko_comp_plotdata %>%
  do(plot=ggplot(data=., aes(x=.x, y=.y)) +
       geom_point(aes(colour=colour_label), size=1) +
       xlab(.$.col) + ylab(.$.row) +
       geom_abline(colour="red", size=0.5) +
       scale_colour_manual(values=c(different="red", not_different="grey50")) +
       theme_minimal() +
       theme(legend.position="none",
             panel.border=element_rect(fill=NA),
             text=element_text(size=16)) +
       tune::coord_obs_pred()
     ) %>%
  pull(plot) %>%
  ggarrange(plotlist=., nrow=1)
cropped_ggsave("../plots/070722_meeting_plots/strep_iii_KOs_scatter.pdf",
               width=14, height=5)

  
KOs_for_lookup <- strep_ko_comp_plotdata %>%
  filter(colour_label=="different") %>%
  pull(KO) %>% 
  unique()

# top_hit_ko_lookup_raw <- lookup_kos(KOs_for_lookup, 1)
# top_hit_ko_lookup_raw <- readRDS("../tmp/strep_ko_lookup.rds")
# get gene names
ko_info_all <- readRDS("../data/kegg/ko_lookup_raw.rds")
names(ko_info_all) <- sapply(ko_info_all, function(x) x$ENTRY)

parsed_top_hit_info <- pbmclapply(ko_info_all[ KOs_for_lookup ], parse_brite, mc.cores=20) %>%
  bind_rows() %>%
  mutate(term1=str_remove(term1, "\\d+$")) %>%
  distinct()

formatted_strep_lookup <- lapply(ko_info_all[ parsed_top_hit_info$KO ], format_query) %>% bind_rows() %>%
  as_tibble() %>%
  rename(KO=entry) %>%
  distinct()

library(ggrepel)
library(drlib)

strep_ko_comp_plotdata2 <- strep_ko_comp_plotdata %>%
  left_join(parsed_top_hit_info, by="KO") %>%
  filter(!is.na(term1), !grepl("Brite Hierarchies|Not Included in Pathway or Brite", term1)) %>%
  mutate(term1=str_remove(term1, "\\d+ ")) %>%
  left_join(
    formatted_strep_lookup %>% select(KO, symbol),
    by="KO"
  ) %>%
  mutate(
    KO_scorediff_rank=dense_rank(desc(score_diff)),
    label=if_else(KO_scorediff_rank<5, symbol, NA_character_)
  ) %>%
  select(-c(term3, term2)) %>%
  distinct() %>%
  slice_max(order_by=abs_score_diff, n=10)

unique_term1_vals <- strep_ko_comp_plotdata2$term1 %>%
  unique()
term1_cmap <- setNames(
  scales::hue_pal()(length(unique_term1_vals)),
  unique_term1_vals
) 

strep_ko_comp_plotdata2 %>% 
  do(plot=ggplot(data=.,
                 aes(x=score_diff, y=reorder(symbol, abs_score_diff))) +
       geom_col(aes(fill=term1), position="dodge") +
       theme_minimal() +
       theme(panel.border=element_rect(fill=NA),
             text=element_text(size=20),
             legend.title=element_blank()) +
       ggtitle(sprintf("%s vs %s", .$.row, .$.col)) +
       ylab("KO symbol") +
       xlab("Score difference") +
       scale_fill_manual(values=term1_cmap)
  ) %>%
  pull(plot) %>%
  ggarrange(plotlist=., nrow=1, common.legend=TRUE, legend="bottom")
cropped_ggsave(
  "../plots/070722_meeting_plots/strep_diffscore_KOs.pdf",
  width=16, height=5
)


emapper %>% 
  rownames_to_column("isolate") %>%
  pivot_longer(-isolate, names_to="KO") %>%
  left_join(cluster_tbl, by="isolate") %>%
  filter(grepl("Strep", cluster_longname)) %>%
  right_join(
    strep_ko_comp_plotdata2 %>%
      ungroup() %>%
      select(KO, symbol) %>%
      distinct(),
    by="KO"
  ) %>%
  mutate(
    cluster_longname=str_replace(cluster_longname, "\\+", ",\n"),
  ) %>%
  ggplot(aes(x=isolate, y=symbol, fill=as.factor(value))) +
  geom_tile() +
  scale_fill_manual(values=c("white", "red")) +
  ylab("KO") + xlab("Isolate") +
  theme_minimal() +
  facet_grid(cols=vars(cluster_longname),
             space="free", scales="free") +
  theme(
    axis.text.x=element_text(size=12, angle=90, vjust=0.5, hjust=1),
    axis.text.y=element_text(size=12),
    strip.text.x=element_text(size=10, angle=0, vjust=0.5, hjust=0.5),
    strip.text.y=element_text(size=10),
    panel.spacing.x=unit(10, "pt"),
    panel.spacing.y=unit(20, "pt"),
    legend.position="none",
    panel.border=element_rect(fill=NA),
    axis.ticks=element_line()
  ) +
  scale_x_discrete(labels=function(x) str_remove_all(x, "\\d+|_"))

cropped_ggsave(
  "../plots/070722_meeting_plots/strep_diffscore_KOs_emapper.pdf",
  width=16, height=7
)

kegg_graphs <- parsed_top_hit_info %>% 
  right_join(
    formatted_strep_lookup %>%
      select(KO, symbol, name) 
  ) %>%
  filter(symbol %in% c("ATP1A", "E4.3.1.4", "creC")) %>%
  select(-c(KO)) %>%
  relocate(symbol) %>%
  mutate(term1=str_remove(term1, "\\d+ ") %>% str_replace_all(" ", "\n"),
         term2=str_remove(term2, "\\d+ ") %>% str_replace_all(" ", "\n")
         ) %>%
  group_by(symbol) %>%
  do(tmp=igraph::graph_from_data_frame(.[,2:3]))

for(row_idx in 1:nrow(kegg_graphs)) {
  symbol_name <- kegg_graphs$symbol[[ row_idx ]]
  g <- kegg_graphs$tmp[[ row_idx ]]
  filename <- sprintf("../plots/070722_meeting_plots/kegg_hierarchy_%s.pdf", symbol_name)
  pdf(filename,
      width=7, height=7)
  
  V(g)$label.cex <- 1.2
  V(g)$label.color <- "black"
  g %>% 
    plot(layout=igraph::layout_as_tree)
  
  dev.off()
  
  # knitr::plot_crop(filename)
  
}

parsed_top_hit_info %>% 
  right_join(
    formatted_strep_lookup %>%
      select(KO, symbol, name) 
  ) %>%
  filter(symbol %in% c("ATP1A", "E4.3.1.4", "creC")) %>%
  select(-c(KO)) %>%
  relocate(symbol, name)%>%
  mutate(term1=str_remove(term1, "\\d+ "),
         term2=str_remove(term2, "\\d+ "),
         term3=str_remove(term3, "\\d+ ")
  ) %>%
  select(symbol, name, term1) %>% 
  distinct() %>%
  miscFuncs::latextable()
###########################################################################
###########################################################################