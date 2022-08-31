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
           fct_relevel("grey", after=Inf)) %>%
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
cropped_ggsave(
  "../plots/supplementary_figures/strep_wgcna_isolate_comp.png",
  dpi=450
)

###########################################################################
###########################################################################