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

subset_samples(Bus_CF1, ID %in% c("11593", "11595")) %>% speedyseq::tax_glom("Genus") %>% plot_bar(fill="Phylum")

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
# compare prevalence and abundance between the two sites
#####################################################################################

compute_summary_stats <- function(x) {
    tibble(
      otu=colnames(x),
      prevalence=apply(x, 2, function(xx) 100*sum(xx>0)/length(xx)),
      total_reads=colSums(x),
      total_reads_per_sample=colSums(x)/nrow(x),
      relative_total_reads=100*colSums(x)/sum(x)
  )
}

otu_summary <- lapply(X, compute_summary_stats) %>%
  bind_rows(.id="site") %>%
  left_join(taxonomy_tbl, by="otu")

for(quantity_name in c("prevalence", "total_reads_per_sample")) {
  cat(quantity_name, "\n")
  
   plot_tbl <- otu_summary %>%
    select(site, otu, !!sym(quantity_name), Phylum) %>%
    gather_pairs(key=site, value=!!sym(quantity_name)) %>%
    do(
      plot=ggplot(
        data=., aes(x=.x, y=.y, colour=Phylum)
      ) +
        geom_point(size=0.5) +
        xlab(.$.col[[1]]) + ylab(.$.row[[1]]) +
        guides(colour=guide_legend(ncol=2)) +
        get_plot_theme(10)
    )
    
  if(quantity_name=="total_reads_per_sample") {
    plot_tbl$plot <- lapply(
      plot_tbl$plot, function(x) x +
        scale_x_log10() +
        scale_y_log10() +
        geom_abline()
    )
  } else if(quantity_name=="prevalence") {
    plot_tbl$plot <- lapply(
      plot_tbl$plot,
      function(x) x + 
        scale_x_continuous(labels=scales::percent_format(scale=1)) +
        scale_y_continuous(labels=scales::percent_format(scale=1)) +
        coord_obs_pred() +
        geom_abline()
    )
  }
   
   fig <- ggarrange(
     plotlist=plot_tbl$plot, common.legend=TRUE, nrow=1,
     legend=if_else(quantity_name=="prevalence", "right", "none")
   )
  
   cropped_ggsave(sprintf("../plots/cf_lll_vs_lul/three_site_%s.png", quantity_name),
          plot=fig,
          width=14, height=6, dpi=450)
}



otu_summary %>%
  left_join(taxonomy_tbl) %>%
  ggplot(aes(x=prevalence, y=relative_total_reads, colour=Phylum)) +
  geom_point() + 
  facet_wrap(~site) +
  get_plot_theme(14) +
  xlab("Prevalence") +
  ylab("Percentage of reads") +
  scale_x_continuous(labels=scales::percent_format(scale=1)) +
  scale_y_continuous(labels=percent_format(scale=1),
                     trans="log10",
                     breaks=10**seq(-5, 1)) +
  theme(legend.position="top")
ggsave("../plots/cf_lll_vs_lul/rel_abund_and_prev_by_site.png",
       height=6, width=12, dpi=450)

# # now split by asthma status
# otu_summary_by_disease <- list()
# for(disease_label in c("Asthma", "Control")) {
#   otu_summary_by_disease[[ disease_label ]] <- lapply(
#     X,
#     function(x) x[rownames(x) %in% sample_metadata$.sample[sample_metadata$Disease==disease_label],] %>%
#       compute_summary_stats()
#   ) %>%
#     bind_rows(.id="site")
# }
# 
# otu_summary_by_disease <- otu_summary_by_disease %>% 
#   bind_rows(.id="Disease") %>%
#   left_join(taxonomy_tbl, by="otu")
# 
# plot_panels <- list()
# for(quantity_name in c("prevalence", "total_reads_per_sample")) {
#   plot_panels[[ quantity_name ]] <- otu_summary_by_disease %>%
#     select(site, otu, !!sym(quantity_name), Phylum, Disease) %>%
#     pivot_wider(names_from=Disease, values_from=!!sym(quantity_name)) %>%
#     ggplot(aes(x=Control, y=Asthma, colour=Phylum)) +
#     geom_point() +
#     xlab("Control") + ylab("Asthma") +
#     ggtitle(quantity_name) +
#     geom_abline() +
#     facet_wrap(~site)
#   
# 
#   if(quantity_name=="total_reads_per_sample") {
#     plot_panels[[ quantity_name ]] <- plot_panels[[ quantity_name ]] +
#       scale_x_log10() + scale_y_log10()
#   } else if(quantity_name=="prevalence") {
#     plot_panels[[ quantity_name ]] <- plot_panels[[ quantity_name ]] +
#       scale_x_continuous(labels=scales::percent_format(scale=1)) +
#       scale_y_continuous(labels=scales::percent_format(scale=1))
#   }
# }
# ggarrange(plotlist=plot_panels, common.legend=TRUE,
#           nrow=2)
# ggsave("../plots/cf_lll_vs_lul/prevalence_and_normed_reads_by_disease.png",
#        width=8, height=10, dpi=450)

#
# ratio of LUL/LLL abundnce at the phylum level
phylum_counts <- subset_samples(Bus_CF1, ID %in% study_ids$sample_id) %>%
  speedyseq::tax_glom("Phylum") %>%
  otu_table() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  left_join(study_ids, by="sample_id") %>%
  as_tibble() %>%
  pivot_longer(-c(sample_id, Study_ID, site), names_to="phylum")

phylum_otu_names <- taxonomy_tbl %>% 
  filter(otu %in% phylum_counts$phylum) %>%
  select(otu, Phylum)

phylum_counts <- phylum_counts %>%
  mutate(
    phylum=setNames(phylum_otu_names$Phylum, phylum_otu_names$otu)[ phylum ]
  ) %>%
  left_join(
    study_ids %>% left_join(sample_metadata %>% select(sample_id, Disease),
                            by="sample_id")
  ) %>%
  select(-sample_id)

phylum_counts %>%
  bind_rows(
    phylum_counts %>% mutate(Disease="all_samples")
  ) %>%
  group_by(phylum, Disease, site) %>%
  summarise(value=sum(value)) %>%
  ungroup() %>%
  gather_pairs(key=site, value=value) %>%
  group_by(.row, .col, Disease) %>%
  do(
    plot=ggplot(data=., aes(x=.x, y=.y, colour=phylum)) +
      geom_point(size=3, stroke=1) +
      scale_y_log10(limits=c(1, 1e7)) + 
      scale_x_log10(limits=c(1, 1e7)) +
      geom_abline() +
      xlab(.$.col[[1]]) + ylab(.$.row[[1]]) +
      ggtitle(.$Disease) +
      theme(aspect.ratio=1)
  ) %>%
  pull(plot) %>%
  ggarrange(plotlist=., nrow=3, ncol=3, common.legend=TRUE)

phylum_counts %>%
  filter(site %in% c("cf_lll", "cf_lul")) %>%
  pivot_wider(names_from=site, values_from=value) %>%
  drop_na() %>%
  ggplot(aes(x=cf_lll, y=cf_lul, colour=Disease)) +
  geom_point() +
  facet_wrap(~phylum) +
  scale_y_log10(limits=c(1, 1e7)) + 
  scale_x_log10(limits=c(1, 1e7)) +
  geom_abline()
  
# ggarrange(
#   qplot(wide_phylum_counts$cf_lll, wide_phylum_counts$cf_lul, colour=wide_phylum_counts$Disease) +
#     scale_y_log10() +  scale_x_log10() + geom_abline() + theme(legend.title=element_blank()),
#   qplot(wide_phylum_counts$cf_lll, wide_phylum_counts$cf_ots, colour=wide_phylum_counts$Disease) +
#     scale_y_log10() +  geom_abline() + scale_x_log10(),
#   qplot(wide_phylum_counts$cf_lul, wide_phylum_counts$cf_ots, colour=wide_phylum_counts$Disease)  +
#     scale_y_log10() +  geom_abline() + scale_x_log10(),
#   nrow=1,
#   common.legend=TRUE,
#   legend="top"
# )

#####################################################################################
#####################################################################################

#####################################################################################
# stacked barplots
#####################################################################################

stacked_barplot_data <- lapply(X, function(x) 100*x/rowSums(x)) %>%
  lapply(function(x) x %>% as.data.frame() %>%
           rownames_to_column("sample_id") %>%
           pivot_longer(-sample_id, names_to="otu", values_to="rel_abund")) %>%
  bind_rows(.id="site") %>%
  left_join(taxonomy_tbl %>% select(otu, Phylum), by="otu") %>%
  left_join(sample_metadata %>% select(sample_id, Disease), by="sample_id")

fig <- stacked_barplot_data %>%
  right_join(study_ids, by=c("site", "sample_id")) %>%
  ggplot(aes(x=Study_ID, y=rel_abund, fill=Phylum)) +
  geom_col(colour="black", size=0.5) +
  ggh4x::facet_grid2(rows=vars(site), cols=vars(Disease), scales="free_x", space="free_x") +
  scale_y_continuous(expand=c(0,0), labels=scales::percent_format(scale=1)) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))

# sample_order <- study_ids %>%
#   pivot_wider(names_from=site, values_from=sample_id, values_fill=NA) %>%
#   mutate(
#     across(contains("cf_"), function(x) if_else(is.na(x), ))
#     # cf_lul=if_else(!is.na(cf_lul), cf_lul, paste0("missing_", Study_ID))
#   ) 
#   left_join(sample_metadata %>% select(sample_id, Disease), by=c("cf_lll"="sample_id")) %>%
#   arrange(Disease) %>%
#   select(-Disease)

cropped_ggsave("../plots/cf_lll_vs_lul/stacked_barplots_aligned.png", plot=fig, height=8, width=18, dpi=450)

#####################################################################################
#####################################################################################


#####################################################################################
# diversity between asthmatics and controls at each site
#####################################################################################

diversity_tbl <- subset_samples(
  Bus_CF1, ID %in% all_sample_ids$sample_id
) %>%
  estimate_richness() %>%
  rownames_to_column("sample_id") %>% 
  as_tibble() %>%
  mutate(sample_id=str_remove(sample_id, "X"))

alpha_metrics <- subset_samples(
  Bus_CF1, ID %in% all_sample_ids$sample_id
) %>%
  microbiome::alpha() %>%
  rownames_to_column("sample_id") %>% 
  as_tibble()

diversity_plotdata <- diversity_tbl %>%
  select(sample_id, Observed, InvSimpson, Shannon) %>%
  left_join(
    alpha_metrics %>%
      select(sample_id, evenness_pielou, dominance_relative),
    by="sample_id"
  ) %>%
  left_join(
    study_ids,
    by="sample_id"
  ) %>%
  left_join(sample_metadata %>% select(sample_id, Disease),
            by="sample_id") %>%
  pivot_longer(-c(sample_id, Study_ID, site, Disease))


  
diversity_wilcox_test <- diversity_plotdata %>%
  group_by(site, name) %>%
  do(rstatix::wilcox_test(data=., formula=value~Disease)) %>%
  ungroup() %>%
  mutate(p_adj=p.adjust(p, method="fdr")) %>%
  left_join(
    diversity_plotdata %>%
      group_by(site, name) %>%
      do(rstatix::wilcox_effsize(data=., formula=value~Disease)) %>%
      ungroup()
  ) %>%
  mutate(
    label=sprintf("p=%.2f\nr=%.2f (%s)", p_adj, effsize, magnitude)
  )


diversity_plotdata %>%
  ggplot(aes(x=site, y=value, fill=Disease)) +
  geom_boxplot() +
  lemon::facet_rep_wrap(~name, scales="free_y") +
  ylab("Diversity value") +
  geom_label(data=diversity_wilcox_test,
            aes(x=site, label=label), y=Inf, vjust=1.5,
            alpha=0.5,
            inherit.aes=FALSE) +
  scale_y_continuous(expand=expansion(mult=c(0.1, 0.5))) +
  theme(legend.position="top") +
  get_plot_theme()
ggsave(
  "../plots/cf_lll_vs_lul/diversity_differences.png",
  height=7, width=14, dpi=450
)


#####################################################################################
#####################################################################################


#####################################################################################
# paired diversity differences between the sites
#####################################################################################

paired_wilcox_with_effsize <- function(x, y, ...) {
  p_vals <- wilcox.test(x=x, y=y, paired=TRUE, ...) %>% broom::tidy()
  eff_sizes <- effectsize::rank_biserial(x=x, y=y, paired=TRUE)
  return(p_vals %>%
           bind_cols(eff_sizes %>%
                       effectsize::interpret_r(rules="cohen1988"))
         )
}

paired_samples <- study_ids %>%
  pivot_wider(names_from=site, values_from=sample_id) %>%
  drop_na() %>%
  pivot_longer(-Study_ID, names_to="site", values_to="sample_id")

paired_samples %>%
  group_by(Study_ID) %>%
  tally() %>%
  pull(n) %>%
  table()

paired_diversity_vals <- diversity_plotdata %>%
  right_join(paired_samples, by=c("Study_ID", "sample_id", "site"))

#
# plot the diversity at each site
paired_diversity_vals %>%
  gather_pairs(key=site, value=value) %>%
  group_by(name, .row, .col) %>%
  do(
    plot=ggplot(data=., aes(x=.x, y=.y)) +
      geom_point(aes(colour=Disease)) +
      xlab(.$.col[[1]]) +
      ylab(.$.row[[1]]) +
      ggtitle(.$name) +
      geom_abline() +
      tune::coord_obs_pred()
  ) %>%
  pull(plot) %>%
  ggarrange(plotlist=., ncol=3, nrow=5, common.legend=TRUE)

# paired_diversity_vals %>%
#   filter(name!="dominance_absolute") %>%
#   ggplot(aes(x=site, y=value)) +
#   geom_point(size=0.75) +
#   geom_line(aes(group=Study_ID), alpha=0.5) +
#   geom_boxplot(alpha=0.8) + 
#   facet_wrap(Disease~name, scales="free_y", nrow=2)
# ggsave("../plots/cf_lll_vs_lul/paired_site_div_boxplots.png",
#        height=6, width=12, dpi=450)
# 
# wilcox_test_result <- paired_diversity_vals %>%
#   filter(name!="dominance_absolute") %>%
#   group_by(Disease, name) %>%
#   gather_pairs(site, value) %>%
#   do(paired_wilcox_with_effsize(.$.x, .$.y))

# kruskal-wallis test as there are three groups
kruskal_res_tbl <- paired_diversity_vals %>%
  filter(name!="dominance_absolute") %>%
  group_by(Disease, name) %>%
  do(kruskal_effsize(data=., formula=value~site)) %>%
  select(-method) %>%
  inner_join(
    paired_diversity_vals %>%
      filter(name!="dominance_absolute") %>%
      group_by(Disease, name) %>%
      do(kruskal_test(data=., formula=value~site)) %>%
      select(-method)
  ) %>%
  ungroup() %>%
  mutate(p_adj=p.adjust(p, method="fdr"))

kruskal_res_tbl %>%
  mutate(
    facet_title=sprintf(
      "%s (%s)\n\np=%.2f, eta_sq=%.2f (%s)",
      name, Disease, p_adj, effsize, magnitude
    )#paste0(name, " (", Disease, ")\np=", p_adj, "eta-sq=", effsize)
  ) %>%
  arrange(Disease, name) %>%
  mutate(facet_title=factor(facet_title, levels=unique(facet_title))) %>%
  left_join(paired_diversity_vals) %>%
  ggplot(aes(x=site, y=value)) +
  geom_point(size=0.75) +
  geom_line(aes(group=Study_ID), alpha=0.5) +
  geom_boxplot(alpha=0.8) + 
  facet_wrap(~facet_title, scales="free_y", nrow=2)

ggsave("../plots/cf_lll_vs_lul/paired_site_div_boxplots.png",
       height=6, width=13, dpi=450)


wilcox_test_result %>%
  filter(name=="InvSimpson")
  ggplot(aes(x=.row, y=.col))

 plot_panel_data <- list()
paired_wilcox_res <- list()
for(disease_label in c("Asthma", "Control")) {
  plotdata <- paired_diversity_vals
  if(disease_label!="All") {
    plotdata <- plotdata %>% filter(Disease==disease_label)
  }
  # plot_panels[[ disease_label ]] <- ggpaired(
  #   plotdata,
  #   x="site", y="value", id="Study_ID",
  #   facet.by="name", scales="free"
  # )
  # wilcox_pvalues <- plotdata %>% 
  #   group_by(name) %>%
  #   do(rstatix::wilcox_test(data=., formula=value~site, paired=TRUE)) %>%
  #   ungroup()
  # wilcox_effsize <- plotdata %>% 
  #   group_by(name) %>%
  #   do(rstatix::wilcox_effsize(data=., formula=value~site, paired=TRUE)) %>%
  #   ungroup()
  plot_panel_data[[ disease_label ]] <- plotdata
  paired_wilcox_res[[ disease_label ]] <- plotdata %>% 
    select(-sample_id) %>%
    pivot_wider(names_from=site) %>%
    group_by(name) %>%
    do(paired_wilcox_with_effsize(.$cf_lul, .$cf_lll))  %>%
    ungroup()
}

paired_wilcox_res <- paired_wilcox_res %>%
  bind_rows(.id="group") %>%
  mutate(p_adj=p.adjust(p.value, method="fdr"))

plot_panel_data %>%
  bind_rows(.id="group") %>%
  ggplot(aes(x=Disease, y=value, fill=site)) +
  geom_boxplot(position=position_dodge(width=0.9)) +
  lemon::facet_rep_wrap(~name, scales="free_y") +
  geom_point(position=position_dodge(width=0.9)) +
  geom_line(aes(group=interaction(Study_ID, name)), position=position_dodge(width=0.9)) +
  geom_label(
    data=paired_wilcox_res %>%
      mutate(label=sprintf("p=%.2f\nr=%.2f (%s)", p_adj, r_rank_biserial, Interpretation)),
    aes(x=group, label=label), y=Inf, vjust=1.5, alpha=0.5,
    inherit.aes=FALSE
  ) +
  scale_y_continuous(expand=expansion(mult=c(0.1, 0.5))) +
  theme(legend.position="top") +
  get_plot_theme(14)

# plot_panel_data %>%
#   bind_rows(.id="group") %>%
#   select(-sample_id) %>%
#   pivot_wider(names_from="site") %>%
#   ggpaired(
#     "cf_lul", "cf_lll",
#     x="Disease", y="value", id="Study_ID",
#     facet.by="name", scales="free_y"
#   )

ggsave("../plots/cf_lll_vs_lul/paired_site_div_boxplots.png",
       height=7, width=11, dpi=450)

paired_wilcox_res %>%
  bind_rows(.id="disease") %>%
  mutate(p=p.adjust(p, method="fdr")) %>%
  select(disease, name, n1, p, effsize, magnitude) %>%
  rename(n=n1) %>%
  mutate(p=format(p, digits=3), effsize=format(effsize, digits=3)) %>%
  data.table::fwrite("../tmp/paired_test_res.csv")


paired_diversity_vals %>% 
  ggplot(aes(x=Disease, y=value, fill=site)) +
  geom_boxplot() +
  lemon::facet_rep_wrap(~name, scales="free_y") +
  ylab("Diversity value")
  # geom_text(data=paired_wilcox_res,
  #           aes(label=p_adj_label), x=Disease, y=Inf, vjust=1.5, size=5,
  #           inherit.aes=FALSE)
ggsave(
  "../plots/cf_lll_vs_lul/site_diversity_differences.png",
  height=4, width=9, dpi=450
)

#####################################################################################
#####################################################################################

#####################################################################################
# PCoA plots based on UniFrac distances
#####################################################################################

get_biplot_data <- function(x, plot.axes=c(1,2)) {
  pr.coo <- x$vectors
  diag.dir <- diag(c(1,1))
  pr.coo[,plot.axes] <- pr.coo[,plot.axes] %*% diag.dir
  print(100*x$values$Relative_eig[ plot.axes ])
  return(pr.coo[,plot.axes] %>%
           as.data.frame() %>%
           rownames_to_column("sample_id") %>%
           as_tibble())
}

uf_dists <- list()
uf_dists[[ "cf_lll" ]] <- subset_samples(
    Bus_CF1, ID %in% all_sample_ids$sample_id[ all_sample_ids$site=="cf_lll" ]
) %>% UniFrac()
uf_dists[[ "cf_lul" ]] <- subset_samples(
  Bus_CF1, ID %in% all_sample_ids$sample_id[ all_sample_ids$site=="cf_lul" ]
) %>% UniFrac()

dist_mats <- lapply(X, function(x) compositions::clr(x) %>% vegan::vegdist(method="euclidean"))

lapply(dist_mats, function(x) x %>% ape::pcoa() %>% get_biplot_data(plot.axes=c(1,3))) %>%
  bind_rows(.id="site") %>%
  left_join(sample_metadata %>% select(sample_id, Disease)) %>%
  as_tibble() %>%
  ggplot(aes(x=Axis.1,
             y=Axis.3,
             colour=Disease)) +
  geom_point() +
  scale_colour_discrete() +
  facet_wrap(~site) +
  stat_ellipse()

#####################################################################################
#####################################################################################

#####################################################################################
# differential abundance analysis - asthmtics vs controls
#####################################################################################

# included_otus <- otu_summary %>%
#   select(site, otu, relative_total_reads) %>%
#   filter(relative_total_reads>0.5) %>%
#   select(otu) %>%
#   distinct()
included_otus <- otu_summary %>%
  group_by(site) %>%
  slice_max(order_by=total_reads, n=20) %>%
  ungroup() %>% select(otu) %>%
  distinct()
cat(sprintf("%d OTUs included in ALDEx2 analysis\n", length(included_otus$otu)))


run_MC_aldex2 <- TRUE
aldex2_test_results <- list()
for(site_name in names(X)) {
  cat(site_name, "\n")
  
  xx <- X[[ site_name ]]
  conditions <- sample_metadata %>%
    filter(sample_id %in% rownames(xx)) %>%
    arrange(match(sample_id, rownames(xx)))
  stopifnot(identical(conditions$sample_id, rownames(xx)))
  
  aldex_obj <- aldex.clr(t(xx), conds=conditions$Disease, useMC=run_MC_aldex2)
  
  ttest_res <- aldex_obj %>% aldex.ttest(paired.test=FALSE)
  effect_res <- aldex_obj %>% aldex.effect(CI=TRUE, verbose=FALSE, useMC=run_MC_aldex2)
  
  aldex2_test_results[[ site_name ]] <- data.frame(ttest_res, effect_res) %>%
    rownames_to_column("taxa") %>%
    as_tibble() %>%
    arrange(wi.eBH)
  
}

aldex2_test_results %>%
  bind_rows(.id="site") %>%
  rename(otu=taxa) %>%
  right_join(included_otus, by="otu") %>%
  group_by(site) %>%
  mutate(wi.eBH=p.adjust(wi.ep, method="fdr"),
         we.eBH=p.adjust(we.ep, method="fdr")
         # effsize_rank
  ) %>%
  mutate(effsize_rank=dense_rank(desc(abs(effect)))) %>%
  ungroup() %>%
  # arrange(we.eBH) %>%
  mutate(label=if_else(effsize_rank<=10, otu, NA_character_)) %>%
  ggplot(aes(x=effect, y=-log10(we.eBH))) +
  geom_point(size=2) +
  # geom_errorbar(aes(xmin=effect.low, xmax=effect.high)) +
  facet_wrap(~site) +
  xlim(c(-0.45, 0.45)) +
  xlab("Effect size") +
  ylab(expression(paste(-log[10], " p-value"))) +
  geom_vline(xintercept=0, linetype="dashed", colour="grey50") +
  # geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="grey50") +
  geom_hline(yintercept=-log10(0.1), linetype="dashed", colour="grey50") +
  get_plot_theme(13) +
  ggrepel::geom_label_repel(
    aes(label=label), max.overlaps=100, force=100
  )
cropped_ggsave("../plots/cf_lll_vs_lul/threesite_volcano_asthma_vs_control.png",
               width=12, height=4, dpi=450)

# compare effect sizes between sites
aldex2_test_results %>%
  bind_rows(.id="site") %>%
  rename(otu=taxa) %>%
  select(site, otu, effect) %>%
  pivot_wider(names_from=site, values_from=effect) %>%
  drop_na() %>%
  pivot_longer(-otu, names_to="site", values_to="effect") %>%
  gather_pairs(key=site, value=effect) %>%
  do(
    plot=ggplot(
      data=., aes(x=.x, y=.y)
    ) +
      geom_point(size=1) +
      xlab(.$.col[[1]]) + ylab(.$.row[[1]]) +
      guides(colour=guide_legend(ncol=2)) +
      get_plot_theme(12) +
      stat_cor() +
      tune::coord_obs_pred() +
      geom_abline()
  ) %>%
  pull(plot) %>%
  ggarrange(plotlist=., nrow=1)

cropped_ggsave("../plots/cf_lll_vs_lul/threesite_aldex_effectsize_agreement_asthma_vs_control.png",
               width=12, height=4, dpi=450)
#####################################################################################
#####################################################################################

#####################################################################################
# random forest two-sample test - distinguishing sites
#####################################################################################

library(tidymodels)
library(ranger)

n_cores <- 16
n_outer_folds <- 5
n_inner_folds <- 10

# returns an RF model
make_rf_mod <- partial(
  rand_forest,
  trees=100,
  engine="ranger",
  mode="classification"
)

# hparam grid
grid_rf2 <- grid_max_entropy(        
  mtry(range=c(1, 40)), 
  min_n(range=c(2, 10)),
  size=50) 


taxa_table_full <- subset_samples(
  Bus_CF1, ID %in% study_ids$sample_id
) %>%
  speedyseq::tax_glom("Genus") %>%
  otu_table() %>%
  t() %>%
  as.data.frame() 

taxa_prevalance <- apply(taxa_table_full, 2, function(x) sum(x>0))
taxa_prevalance %>% qplot()
taxa_table_full <- taxa_table_full / rowSums(taxa_table_full)
taxa_table_full <- taxa_table_full[ , taxa_prevalance>3 ]

rf_df <- taxa_table_full %>%
  rownames_to_column("sample_id") %>%
  left_join(study_ids %>% select(-Study_ID)) %>%
  column_to_rownames("sample_id")  %>%
  mutate(site=as.factor(site))

# nested CV folds
cv_resamples <- nested_cv(
  rf_df,
  outside=vfold_cv(v=n_outer_folds, repeats=1, strata=site),
  inside=vfold_cv(v=n_inner_folds, repeats=1, strata=site)
)

# run nested CV
library(doParallel)
cl <- makeForkCluster(n_cores)
registerDoParallel(cl)

outer_fold_preds <- list()
for(outer_fold_idx in 1:nrow(cv_resamples)) {
  cat(outer_fold_idx, "\n")
  
  # hparam selection on training set
  rf_fit <- tune_grid(
    workflow() %>% 
      add_model(make_rf_mod(mtry=tune(), min_n=tune())) %>% 
      add_recipe(recipe(site ~ ., data=rf_df)),
    resamples=cv_resamples$inner_resamples[[ outer_fold_idx ]],
    grid=grid_rf2,
    # metrics=metric_set(roc_aunu(truth=site, .pred_cf_lll, .pred_cf_lul, .pred_cf_ots)),
    control=control_grid(verbose=TRUE, allow_par=TRUE, parallel_over="everything")
  )
  
  # best hparam is the one that maximises the AU-ROC
  best_hparams <- rf_fit$.metrics %>%
    bind_rows(.id="inner_fold") %>%
    filter(.metric=="roc_auc") %>%
    group_by(mtry, min_n, .metric) %>%
    summarise(mean_estimate=mean(.estimate, na.rm=TRUE)) %>%
    ungroup() %>%
    group_by(.metric) %>%
    slice_max(mean_estimate, n=1, with_ties=FALSE)
  stopifnot(nrow(best_hparams)==1)
  
  # evaluate on validation set, saving predictions and labels
  outer_fold_preds[[ outer_fold_idx ]] <- make_rf_mod(
    mtry=best_hparams$mtry, min_n=best_hparams$min_n
  ) %>%
    parsnip::fit(site ~ ., data=analysis(cv_resamples$splits[[ outer_fold_idx ]])) %>%
    predict(assessment(cv_resamples$splits[[ outer_fold_idx ]]), type="prob") %>%
    bind_cols(assessment(cv_resamples$splits[[ outer_fold_idx ]]) %>% select(site)) 
}

stopCluster(cl)

all_validation_preds <- outer_fold_preds %>%
  bind_rows(.id="outer_fold")
hard_preds <- all_validation_preds %>%
  mutate(
    .pred=probably::make_class_pred(
      .pred_cf_lll, .pred_cf_lul, .pred_cf_ots,
      levels=c("cf_lll", "cf_lul", "cf_ots"),
    )
  )

# multiclass AUC (one-vs-all)
multiclass_roc_aucs <- list()
for(class_name in levels(rf_df$site)) {
  multiclass_roc_aucs[[ class_name ]] <- all_validation_preds %>%
    mutate(
      site=if_else(as.character(site)==class_name, as.character(site), paste0("not_", class_name)),
      site=factor(site, levels=c(class_name, paste0("not_", class_name)))
    ) %>% 
    group_by(outer_fold) %>%
    roc_auc(estimate=!!sym(paste0(".pred_", class_name)), truth=site) %>%
    summarise(.estimate=mean(.estimate))
    
}
multiclass_roc_aucs



roc_curve_multiclass <- all_validation_preds %>%
  roc_curve(.pred_cf_lll, .pred_cf_lul, .pred_cf_ots, truth=site) %>%
  autoplot() +
  theme(legend.position=c(0.8, 0.2)) +
  get_plot_theme(10) +
  theme(panel.spacing.x=unit(5, "mm")) +
  tune::coord_obs_pred()

n_by_class <- rf_df$site %>%
  table()

pr_baseline_prec <- as.data.frame(n_by_class /  sum(n_by_class)) %>%
  rename(site=".")

pr_curve_multiclass <- all_validation_preds %>%
  # bind_rows(
  #   all_validation_preds %>% mutate(outer_fold="pooled")
  # ) %>%
  pr_curve(.pred_cf_lll, .pred_cf_lul, .pred_cf_ots, truth=site) %>%
  autoplot() +
  theme(legend.position=c(0.8, 0.2)) +
  get_plot_theme(10) +
  theme(panel.spacing.x=unit(5, "mm")) +
  tune::coord_obs_pred() +
  geom_hline(data=pr_baseline_prec %>% rename(.level=site),
             aes(yintercept=Freq),
             linetype="dotted")

ggarrange(
  roc_curve_multiclass, pr_curve_multiclass,
  nrow=2
)
cropped_ggsave("../plots/cf_lll_vs_lul/multiclass_roc_pr_curves_site.png",
       height=5, width=10, dpi=450)

hard_preds %>%
  conf_mat(truth=site, estimate=.pred) %>%
  autoplot(type="heatmap") 

hard_preds %>%
  mcc(truth=site, estimate=.pred)

all_validation_preds %>%
  bind_rows(
    all_validation_preds %>% mutate(outer_fold="pooled")
  ) %>%
  group_by(site, outer_fold) %>%
  roc_auc(Disease, .pred_Control)

tmp <- all_validation_preds %>%
  group_by(site) %>%
  summarise(roc_out=pROC::roc(response=Disease, predictor=.pred_Asthma, levels=c("Control", "Asthma")) %>% list())
tmp$roc_out[[1]]
lapply(tmp$roc_out, pROC::auc)
lapply(tmp$roc_out, pROC::ci.auc)
pROC::roc.test(tmp$roc_out[[1]], tmp$roc_out[[2]],
               alternative="greater",
               paired=FALSE) %>% broom::tidy()
pROC::roc.test(tmp$roc_out[[1]], tmp$roc_out[[3]],
               alternative="greater",
               paired=FALSE) %>% broom::tidy()

sample_metadata %>%
  select(sample_id, Disease) %>%
  filter(sample_id %in% samples_for_rf$cf_lll) %>%
  pull(Disease) %>%
  table()

  

ggsave("../plots/cf_lll_vs_lul/rf_perf_curves.png",
       width=9, height=4, dpi=450)
# lift and gain curves
all_validation_preds %>%
  group_by(site) %>%
  gain_curve(Disease, .pred_Control) %>%
  autoplot()

all_validation_preds %>%
  group_by(site) %>%
  lift_curve(Disease, .pred_Control) %>%
  autoplot()

# hard label metrics for different thresholds
lapply(
  seq(0, 1, length.out=100),
  function(eps) all_validation_preds %>%
  group_by(site) %>%
  mutate(
    .pred=probably::make_two_class_pred(
      .pred_Control, levels=c("Control", "Asthma"), threshold=eps
    )
  ) %>%
  metric_set(accuracy, bal_accuracy, kap, mcc, f_meas, precision, recall, sens, spec, j_index)(Disease, estimate=.pred) %>%
  mutate(.threshold=eps)
) %>%
  bind_rows() %>%
  ggplot(aes(x=.threshold, y=.estimate)) +
  geom_line(aes(group=site, colour=site)) +
  facet_wrap(~.metric, scales="free_y", nrow=2) +
  theme(legend.position="top")
ggsave("../plots/cf_lll_vs_lul/hard_class_rf_metrics.png",
       width=10, height=4, dpi=450)


# all_validation_preds %>%
#   group_by(site) %>%
#   probably::threshold_perf(Disease, .pred_Control, thresholds=seq(0, 1, length.out=100))%>%
#   ggplot(aes(x=.threshold, y=.estimate)) +
#   geom_line(aes(group=site, colour=site)) +
#   facet_wrap(~.metric, scales="free_y")

tidy_cvauc <- function(...) {
  tmpout <- cvAUC::ci.cvAUC(...)
  return(tibble(
    cvAUC=tmpout$cvAUC, se=tmpout$se,
    ci_lower=tmpout$ci[[1]],  ci_upper=tmpout$ci[[2]],
    confidence=tmpout$confidence
  ))
}

all_validation_preds %>%
  group_by(site, outer_fold) %>%
  metric_set(roc_auc, pr_auc)(Disease, .pred_Control) %>%
  group_by(site, .metric) %>%
  summarise(mean_estimate=mean(.estimate)) %>%
  pivot_wider(names_from=.metric, values_from=mean_estimate)

all_validation_preds %>%
  group_by(site) %>%
  summarise(
    cvauc=list(tidy_cvauc(predictions=.pred_Asthma,
                         labels=as.integer(Disease)-1,
                         folds=outer_fold))
  ) %>%
  unnest(cvauc)





#####################################################################################
#####################################################################################

#####################################################################################
# random forest two-sample test - asthmatics vs controls
#####################################################################################

library(tidymodels)
library(ranger)

n_cores <- 16
n_outer_folds <- 10
n_inner_folds <- 10

# returns an RF model
make_rf_mod <- partial(
  rand_forest,
  trees=1000,
  engine="ranger",
  mode="classification"
)

# hparam grid
grid_rf2 <- grid_max_entropy(        
  mtry(range=c(1, 40)), 
  min_n(range=c(2, 10)),
  size=50) 

my_metrics <- metric_set(roc_auc, accuracy, pr_auc)

samples_for_rf <- study_ids %>%
  pivot_wider(names_from=site, values_from=sample_id) %>%
  drop_na()

library(doParallel)
cl <- makeForkCluster(n_cores)
registerDoParallel(cl)

outer_fold_preds <- list()
for(site_name in names(X)) {
  cat(site_name, "\n")
  outer_fold_preds[[ site_name ]] <- list()
  
  # get genus counts, removing very rare taxa
  taxa_table_full <- subset_samples(
    Bus_CF1, ID %in% samples_for_rf[[ site_name ]]
  ) %>%
    speedyseq::tax_glom("Genus") %>%
    otu_table() %>%
    t() %>%
    as.data.frame()
  taxa_prevalance <- apply(taxa_table_full, 2, function(x) sum(x>0))
  taxa_table_full <- taxa_table_full / rowSums(taxa_table_full)
  taxa_table_full <- taxa_table_full[ , taxa_prevalance>3 ]
  
  rf_df <- taxa_table_full %>%
    rownames_to_column("sample_id") %>%
    left_join(sample_metadata %>% select(sample_id, Disease),
              by="sample_id") %>%
    mutate(Disease=factor(Disease, levels=c("Control", "Asthma"))) %>%
    select(-c(sample_id))
  
  # nested CV folds
  cv_resamples <- nested_cv(
    rf_df,
    outside=vfold_cv(v=n_outer_folds, repeats=1, strata=Disease),
    inside=vfold_cv(v=n_inner_folds, repeats=1, strata=Disease)
  )
  
  for(outer_fold_idx in 1:nrow(cv_resamples)) {
    cat(outer_fold_idx, "\n")
    
    # hparam selection on training set
    rf_fit <- tune_grid(
      workflow() %>% 
        add_model(make_rf_mod(mtry=tune(), min_n=tune())) %>% 
        add_recipe(recipe(Disease ~ ., data=rf_df)),
      resamples=cv_resamples$inner_resamples[[ outer_fold_idx ]],
      grid=grid_rf2,
      metrics=my_metrics,
      control=control_grid(verbose=FALSE, allow_par=TRUE, parallel_over="everything")
    )
    
    # best hparam is the one that maximises the AU-ROC
    best_hparams <- rf_fit$.metrics %>%
      bind_rows(.id="inner_fold") %>%
      filter(.metric=="roc_auc") %>%
      group_by(mtry, min_n, .metric) %>%
      summarise(mean_estimate=mean(.estimate, na.rm=TRUE)) %>%
      ungroup() %>%
      group_by(.metric) %>%
      slice_max(mean_estimate, n=1, with_ties=FALSE)
    stopifnot(nrow(best_hparams)==1)
    
    # evaluate on validation set, saving predictions and labels
    outer_fold_preds[[site_name]][[ outer_fold_idx ]] <- make_rf_mod(
      mtry=best_hparams$mtry, min_n=best_hparams$min_n
    ) %>%
      parsnip::fit(Disease ~ ., data=analysis(cv_resamples$splits[[ outer_fold_idx ]])) %>%
      predict(assessment(cv_resamples$splits[[ outer_fold_idx ]]), type="prob") %>%
      bind_cols(assessment(cv_resamples$splits[[ outer_fold_idx ]]) %>% select(Disease)) 
  }
}
stopCluster(cl)

all_validation_preds <- outer_fold_preds %>%
  lapply(function(x) bind_rows(x, .id="outer_fold")) %>%
  bind_rows(.id="site")

roc_plot <- all_validation_preds %>%
  # bind_rows(
  #   all_validation_preds %>% mutate(outer_fold="pooled")
  # ) %>%
  group_by(site) %>%
  roc_curve(Disease, .pred_Control) %>%
  autoplot() +
  theme(legend.position=c(0.8, 0.2)) +
  get_plot_theme() +
  tune::coord_obs_pred()

all_validation_preds %>%
  bind_rows(
    all_validation_preds %>% mutate(outer_fold="pooled")
  ) %>%
  group_by(site, outer_fold) %>%
  roc_auc(Disease, .pred_Control) %>%
  group_by(site)

tmp <- all_validation_preds %>%
  group_by(site) %>%
  summarise(roc_out=pROC::roc(response=Disease, predictor=.pred_Asthma, levels=c("Control", "Asthma")) %>% list())
tmp$roc_out[[1]]
lapply(tmp$roc_out, pROC::auc)
lapply(tmp$roc_out, pROC::ci.auc)
pROC::roc.test(tmp$roc_out[[1]], tmp$roc_out[[2]],
               alternative="greater",
               paired=FALSE) %>% broom::tidy()
pROC::roc.test(tmp$roc_out[[1]], tmp$roc_out[[3]],
               alternative="greater",
               paired=FALSE) %>% broom::tidy()

sample_metadata %>%
  select(sample_id, Disease) %>%
  filter(sample_id %in% samples_for_rf$cf_lll) %>%
  pull(Disease) %>%
  table()
pr_plot <- all_validation_preds %>%
  group_by(site) %>%
  pr_curve(Disease, .pred_Control) %>%
  autoplot() +
  theme(legend.position=c(0.8, 0.2)) +
  get_plot_theme() +
  tune::coord_obs_pred() +
  geom_hline(yintercept=16/(39+16), linetype="dotted")


# lift and gain curves
all_validation_preds %>%
  group_by(site) %>%
  gain_curve(Disease, .pred_Control) %>%
  autoplot()

all_validation_preds %>%
  group_by(site) %>%
  lift_curve(Disease, .pred_Control) %>%
  autoplot()

# hard label metrics for different thresholds
lapply(
  seq(0, 1, length.out=100),
  function(eps) all_validation_preds %>%
    group_by(site) %>%
    mutate(
      .pred=probably::make_two_class_pred(
        .pred_Control, levels=c("Control", "Asthma"), threshold=eps
      )
    ) %>%
    metric_set(accuracy, bal_accuracy, kap, mcc, f_meas, precision, recall, sens, spec, j_index)(Disease, estimate=.pred) %>%
    mutate(.threshold=eps)
) %>%
  bind_rows() %>%
  ggplot(aes(x=.threshold, y=.estimate)) +
  geom_line(aes(group=site, colour=site)) +
  facet_wrap(~.metric, scales="free_y", nrow=2) +
  theme(legend.position="top")
ggsave("../plots/cf_lll_vs_lul/hard_class_rf_metrics.png",
       width=10, height=4, dpi=450)


# all_validation_preds %>%
#   group_by(site) %>%
#   probably::threshold_perf(Disease, .pred_Control, thresholds=seq(0, 1, length.out=100))%>%
#   ggplot(aes(x=.threshold, y=.estimate)) +
#   geom_line(aes(group=site, colour=site)) +
#   facet_wrap(~.metric, scales="free_y")

tidy_cvauc <- function(...) {
  tmpout <- cvAUC::ci.cvAUC(...)
  return(tibble(
    cvAUC=tmpout$cvAUC, se=tmpout$se,
    ci_lower=tmpout$ci[[1]],  ci_upper=tmpout$ci[[2]],
    confidence=tmpout$confidence
  ))
}

all_validation_preds %>%
  group_by(site, outer_fold) %>%
  metric_set(roc_auc, pr_auc)(Disease, .pred_Control) %>%
  group_by(site, .metric) %>%
  summarise(mean_estimate=mean(.estimate)) %>%
  pivot_wider(names_from=.metric, values_from=mean_estimate)

all_validation_preds %>%
  group_by(site) %>%
  summarise(
    cvauc=list(tidy_cvauc(predictions=.pred_Asthma,
                          labels=as.integer(Disease)-1,
                          folds=outer_fold))
  ) %>%
  unnest(cvauc)


ggarrange(
  roc_plot, pr_plot,
  common.legend=TRUE
)
ggsave("../plots/cf_lll_vs_lul/rf_perf_curves.png",
       width=9, height=4, dpi=450)


#####################################################################################
#####################################################################################

#####################################################################################
# MOFA
#####################################################################################

library(MOFA2)
reticulate::use_condaenv("R4")

paired_samples_wide <- paired_samples %>%
  left_join(sample_metadata %>%
              select(sample_id, Disease, CurrentSmoker),
            by="sample_id") %>%
  pivot_wider(names_from=site, values_from=sample_id)

# prepare data - agglomerate to Genus then CLR
genus_phyloseq <- Bus_CF1 %>%
  speedyseq::tax_glom("Genus") 

X_mofa <- list()
for(site_name in names(X)) {
  X_mofa[[ site_name ]] <- phyloseq::subset_samples(genus_phyloseq, ID %in% paired_samples_wide[[ site_name ]]) %>%
    otu_table() %>%
    as.data.frame() %>%
    t() %>%
    # compositions::clr() %>%
    # as.matrix() %>%
    magrittr::set_rownames(paired_samples_wide$Study_ID) %>%
    as.data.frame() 
}
lapply(X_mofa, dim)
lapply(X_mofa, class)
stopifnot(all(sapply(X_mofa, nrow)==55))

genus_prevalences <- lapply(
  X_mofa,
  function(x) tibble(genus=colnames(x), num_nonzero=apply(x, 2, function(xx) sum(xx>0)))
) %>%
  bind_rows(.id="site")
genus_prevalences %>%
  ggplot(aes(x=num_nonzero, fill=site)) +
  geom_histogram() +
  xlab("Num samples") + ylab("Num genera")
mofa_genera <- genus_prevalences %>%
  filter(num_nonzero>=3) %>%
  pull(genus) %>%
  unique()
cat(sprintf("Including %d genera in MOFA analysis\n", length(mofa_genera)))

X_mofa <- lapply(X_mofa, function(x) compositions::clr(x)[,mofa_genera] %>% t())
lapply(X_mofa, dim)

MOFAobject <- create_mofa(X_mofa)
mofa_opts <- get_default_data_options(MOFAobject)
mofa_opts$scale_views <- FALSE
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 5
MOFAobject <- prepare_mofa(MOFAobject, data_options=mofa_opts, model_options=model_opts)
MOFAobject <- run_mofa(MOFAobject)

plot_variance_explained(MOFAobject) +
  geom_label(aes(label=format(value, digits=2, scientific=FALSE)))

plot_factors(MOFAobject, factors=1:3, color_by=paired_samples_wide$Disease)


#####################################################################################
#####################################################################################


#####################################################################################
# network analysis
#####################################################################################

library(ComplexHeatmap)
library(dendextend)

# TODO: add p-values

# visualise correlation at each site
corr_mats <- lapply(X, function(x) x %>% compositions::clr() %>% cor(method="spearman"))

# only plot the most abundant OTUs
plotted_otus <- otu_summary %>%
  group_by(site) %>%
  slice_max(total_reads, n=50) %>%
  ungroup() %>%
  pull(otu) %>%
  unique()

corr_mats <- lapply(corr_mats, function(x) x[plotted_otus,plotted_otus])
otu_hclusts <- lapply(
  corr_mats,
  function(x) as.dist(1-x) %>% hclust())
cor_coph <- as.dendlist(lapply(otu_hclusts, as.dendrogram)) %>%
  cor.dendlist()
rownames(cor_coph) <- names(otu_hclusts)
colnames(cor_coph) <- names(otu_hclusts)

annotation_tbl <- taxonomy_tbl %>%
  filter(otu %in% plotted_otus) %>%
  arrange(match(otu, plotted_otus))
unique_phyla <- unique(annotation_tbl$Phylum)
phylum_colours <- setNames(brewer_pal("Set1", type="qual")(length(unique_phyla)), unique_phyla)

corr_hms <- list()
for(site_name in names(corr_mats)) {
  cat(site_name, "\n")
  corr_hms[[ site_name ]] <- Heatmap(
    corr_mats[[ site_name ]],
    col=circlize::colorRamp2(c(-1, 0, 1), c("Darkblue", "white", "red")),
    name=site_name,
    heatmap_width=unit(8, "inch"), heatmap_height=unit(8, "inch"),
    cluster_rows=otu_hclusts[[ site_name ]],
    cluster_columns=otu_hclusts[[ site_name ]],
    row_dend_width=unit(25, "mm"), column_dend_height=unit(25, "mm"),
    row_names_gp=gpar(fontsize=7), column_names_gp=gpar(fontsize=7),
    top_annotation=columnAnnotation(
      phylum=setNames(annotation_tbl$Phylum, annotation_tbl$otu),
      col=list(phylum=phylum_colours)
    ),
    left_annotation=rowAnnotation(
      phylum=setNames(annotation_tbl$Phylum, annotation_tbl$otu),
      show_legend=FALSE,
      col=list(phylum=phylum_colours)
    )
  )
  
  plot_filename <- sprintf("../plots/cf_lll_vs_lul/corr_heatmap__%s.png", site_name)
  
  png(plot_filename, height=1000, width=1000, res=100)
  draw(corr_hms[[ site_name ]], merge_legend=TRUE)
  dev.off()
  knitr::plot_crop(plot_filename)
}



matrix_trace <- function(x) sum(diag(x))
cor_mat_dist <- function(c1, c2) {
  d <- 1 - matrix_trace(c1 %*% c2) / (norm(c1, type="F") * norm(c2, type="F"))
  return(d)
}
cor_mat_dist(corr_mats$cf_lll, corr_mats$cf_lul)
cor_mat_dist(corr_mats$cf_lll, corr_mats$cf_ots)
cor_mat_dist(corr_mats$cf_lul, corr_mats$cf_ots)


vegan::mantel(1-corr_mats$cf_lll, 1-corr_mats$cf_lul)
vegan::mantel(1-corr_mats$cf_lll, 1-corr_mats$cf_ots)
vegan::mantel(1-corr_mats$cf_lul, 1-corr_mats$cf_ots)

library(bootnet)

paired_samples_wide

glommed_phyloseq <- subset_samples(Bus_CF1, ID %in% paired_samples_wide$cf_lul) %>%
  speedyseq::tax_glom("Genus") 

xx <- glommed_phyloseq %>% otu_table() %>% t() %>% as.data.frame()%>%
  compositions::clr() %>%
  as.matrix()


library(propr)

prop_mats <- lapply(X, function(x) phis(x %>% as.data.frame() %>% as.matrix()))

prop_hclusts <- lapply(prop_mats, function(x) x@matrix[plotted_otus,plotted_otus] %>% as.dist() %>% hclust())

for(site_name in names(prop_mats)) {
  cat(site_name, "\n")
  hm <- Heatmap(
    prop_mats[[ site_name ]]@matrix[plotted_otus,plotted_otus],
    col=circlize::colorRamp2(c(0, 4), c("white", "red")),
    name=site_name,
    heatmap_width=unit(8, "inch"), heatmap_height=unit(8, "inch"),
    cluster_rows=prop_hclusts[[ site_name ]],
    cluster_columns=prop_hclusts[[ site_name ]],
    row_dend_width=unit(25, "mm"), column_dend_height=unit(25, "mm"),
    row_names_gp=gpar(fontsize=7), column_names_gp=gpar(fontsize=7),
    top_annotation=columnAnnotation(
      phylum=setNames(annotation_tbl$Phylum, annotation_tbl$otu),
      col=list(phylum=phylum_colours)
    ),
    left_annotation=rowAnnotation(
      phylum=setNames(annotation_tbl$Phylum, annotation_tbl$otu),
      show_legend=FALSE,
      col=list(phylum=phylum_colours)
    )
  )
  
  plot_filename <- sprintf("../plots/cf_lll_vs_lul/propr_heatmap__%s.png", site_name)
  
  png(plot_filename, height=1000, width=1000, res=100)
  draw(hm, merge_legend=TRUE)
  dev.off()
  knitr::plot_crop(plot_filename)
}


#
# bootsrtap hclust
n_hclust_boots <- 1000
cor_dist <- function(cor_mat) as.dist(1-cor_mat)

plotted_otus <- otu_summary %>%
  group_by(site) %>%
  slice_max(total_reads, n=50) %>%
  ungroup() %>%
  pull(otu) %>%
  unique()

X_subset <- X$cf_ots %>%
  compositions::clr() %>%
  .[,plotted_otus]

obs_hclust <-  X_subset %>%
  cor(method="spearman") %>%
  cor_dist() %>%
  hclust()

booted_hclusts <- pbmcapply::pbmclapply(
  1:n_hclust_boots,
  function(...) X_subset[sample(nrow(X_subset), replace=TRUE, size=nrow(X_subset)),] %>%
    cor(method="spearman") %>%
    cor_dist() %>%
    hclust(),
  mc.cores=10
)

hclust_list <- list(list(obs_hclust), booted_hclusts) %>% unlist(recursive=FALSE)
tmp <- shipunov::Bclust(hclist=hclust_list, relative=TRUE, mc.cores=10)
plot(tmp)


#
# 
library(cluster)

do_cor_clust <- function(x) {
  out <- cor(t(x), method="spearman") %>%
    cor_dist() %>%
    hclust() %>%
    cutree(k=k)
  # return(list(cluster=out %>% unlist()))
  return(out %>% unname())
}

X_subset <- lapply(X, function(x)x %>%
                     compositions::clr() %>%
                     .[,plotted_otus])

booted_gapstats <- lapply(X_subset, function(x) clusGap(x, 20, FUNcluster=do_cor_clust, B=100))

lapply(booted_gapstats, function(x) as_tibble(x$Tab) %>% mutate(k=1:n())) %>%
  bind_rows(.id="site") %>%
  ggplot(aes(x=k, y=gap, colour=site)) +
  geom_line() +
  geom_point()

library(ClusBoot)
clustboot_res <- clusboot(t(X_subset$cf_lll), B=100, clustering.func=do_cor_clust)
rownames(clustboot_res$proportions) <- colnames(X_subset$cf_lll)
colnames(clustboot_res$proportions) <- colnames(X_subset$cf_lll)

otu_plotorder <- clustboot_res$clustering

otu_clust_dist <- as.dist(1-clustboot_res$proportions) %>% hclust()

Heatmap(
  clustboot_res$proportions,
  cluster_rows=otu_clust_dist,
  cluster_columns=otu_clust_dist
)

#
# network stability
column_vars <- as.matrix(xx) %>% matrixStats::colVars()
column_vars %>% qplot()

xx <- xx[ ,column_vars>1.5]
dim(xx)

xx %>% Rfast::standardise() %>% as.data.frame() %>% rownames_to_column("sample_id") %>% pivot_longer(-sample_id) %>% ggplot(aes(x=name, y=sample_id, fill=value)) +
  geom_tile() +
  scale_fill_gradient2()

bootnet::estimateNetwork(xx, default="cor",  corArgs = list(method = "spearman")) -> tmpres


fig <- plot(tmpres)

bootnet::bootnet(xx %>% Rfast::standardise(), nBoots=100, default="cor",  corArgs = list(method = "spearman")) -> tmpres
plot(tmpres, plot="interval") -> tmpplot

minArea <- "q2.5"
maxArea <- "q97.5"
meanVar <- "mean"
legendNcol <- 1

tmpres %>%
  summary() %>%
  ggplot(aes_string(x = 'id', y = meanVar, group = 'id', colour = 'id',ymin = minArea, ymax = maxArea, fill = "id")) + 
  geom_errorbar(position =  position_dodge(width = 0.4)) +
  geom_point(position =  position_dodge(width = 0.4)) +
  geom_line(position =  position_dodge(width = 0.4)) +
  theme_bw() + 
  theme(legend.position="none")
  xlab("Sampled nodes") + ylab("") + 
  guides(fill=guide_legend(ncol=legendNcol),colour=guide_legend(ncol=legendNcol)) + 
  scale_x_reverse(breaks = seq(0.9,0.1,by=-0.1) * ncol(x$sample$graph), labels=c(paste0(seq(90,10,by=-10),"%")),
                  limits = c( ncol(x$sample$graph)-1, 1))

#####################################################################################
#####################################################################################



#####################################################################################
# Dirichlet Multinomial Mixtures
#####################################################################################

library(DirichletMultinomial)

# prepare data - agglomerate to Genus
glommed_phyloseq <- Bus_CF1 %>%
  speedyseq::tax_glom("Class") 

paired_samples_wide <- paired_samples %>%
  left_join(sample_metadata %>%
              select(sample_id, Disease, CurrentSmoker),
            by="sample_id") %>%
  pivot_wider(names_from=site, values_from=sample_id)

X_dmn <- list()
for(site_name in names(X)) {
  X_dmn[[ site_name ]] <- phyloseq::subset_samples(glommed_phyloseq, ID %in% paired_samples_wide[[ site_name ]]) %>%
    otu_table() %>%
    as.data.frame() %>%
    t() %>%
    # compositions::clr() %>%
    # as.matrix() %>%
    magrittr::set_rownames(paired_samples_wide$Study_ID) %>%
    as.data.frame() 
}
lapply(X_dmn, dim)
lapply(X_dmn, class)
stopifnot(all(sapply(X_dmn, nrow)==55))

taxa_prevalences <- lapply(
  X_dmn,
  function(x) tibble(taxa=colnames(x), num_nonzero=apply(x, 2, function(xx) sum(xx>0)))
) %>%
  bind_rows(.id="site")
taxa_prevalences %>%
  ggplot(aes(x=num_nonzero, fill=site)) +
  geom_histogram() +
  xlab("Num samples") + ylab("Num taxa")
dmn_taxa <- taxa_prevalences %>%
  filter(num_nonzero>=3) %>%
  pull(taxa) %>%
  unique()
cat(sprintf("Including %d taxa in DMN analysis\n", length(dmn_taxa)))
X_dmn <- lapply(X_dmn, function(x) x[,dmn_taxa] %>% as.matrix())
lapply(X_dmn, dim)

#
# run DMM analysis
get_gof <- function(fit) {
  return(c("laplace"=laplace(fit), "AIC"=AIC(fit), "BIC"=BIC(fit)))
}

max_k <- 10
SEED <- 123412345

dmn_fits <- lapply(
  X_dmn,
  function(x) pbmcapply::pbmclapply(
    1:max_k,
    function(k) dmn(x, k=k, seed=SEED),
    mc.cores=10
  )
)

dmn_fit_vals <- lapply(dmn_fits,
                       function(x) lapply(x,
                                          function(xx) goodnessOfFit(xx))
                       %>% bind_rows() %>%
                         rownames_to_column("k")
) %>%
  bind_rows(.id="site") %>%
  pivot_longer(-c(site, k)) %>%
  filter(name %in% c("AIC", "BIC", "Laplace")) %>%
  mutate(k=as.integer(k))

dmn_fit_vals %>%
  ggplot(aes(x=k, y=value, colour=site,
             group=interaction(site, name))) +
  geom_point() +
  geom_line() +
  lemon::facet_rep_wrap(~name,
                        scales="free_y") +
  theme_classic(base_size=16) +
  xlab("Number of clusters") + 
  ylab("Goodness of fit metric value")

best_k <- c(cf_lll=2, cf_lul=3, cf_ots=2)

dmn_fits$cf_lll[[2]]

best_fits <- list(
  cf_lll=dmn_fits$cf_lll[[2]], cf_lul=dmn_fits$cf_lul[[3]], cf_ots=dmn_fits$cf_ots[[2]]
)

dmm_clusters <- lapply(
  best_fits,
  function(x) mixture(x, assign=TRUE) %>% as.data.frame() %>%
    rownames_to_column("Study_ID") %>%
    rename("."="cluster")
) %>%
  bind_rows(.id="site") %>%
  pivot_wider(names_from=site, values_from=cluster)

WGCNA::randIndex(table(dmm_clusters$cf_lll, dmm_clusters$cf_lul))
WGCNA::randIndex(table(dmm_clusters$cf_lll, dmm_clusters$cf_ots))
WGCNA::randIndex(table(dmm_clusters$cf_lul, dmm_clusters$cf_ots))

WGCNA::overlapTable(dmm_clusters$cf_lll, dmm_clusters$cf_lul)
WGCNA::overlapTable(dmm_clusters$cf_lll, dmm_clusters$cf_ots)
WGCNA::overlapTable(dmm_clusters$cf_lul, dmm_clusters$cf_ots)

stacked_barplot_data <- lapply(X, function(x) 100*x/rowSums(x)) %>%
  lapply(function(x) x %>% as.data.frame() %>%
           rownames_to_column("sample_id") %>%
           pivot_longer(-sample_id, names_to="otu", values_to="rel_abund")) %>%
  bind_rows(.id="site") %>%
  left_join(taxonomy_tbl %>% select(otu), by="otu") %>%
  left_join(sample_metadata %>% select(sample_id, Disease), by="sample_id")

stacked_barplot_data %>%
  group_by(site, sample_id, Class) %>%
  summarise(rel_abund=sum(rel_abund)) %>%
  ungroup() %>%
  left_join(study_ids, by=c("sample_id", "site")) %>%
  right_join(
    dmm_clusters %>% pivot_longer(-Study_ID, names_to="site", values_to="dmm_cluster"),
    by=c("Study_ID", "site")
  ) %>%
  group_by(site) %>%
  do(
    plot=ggplot(data=., aes(x=Study_ID, y=rel_abund, fill=Class)) +
      geom_col(colour="black", size=0.5) +
      ggforce::facet_row(~dmm_cluster, scales="free_x", space="free") +
      scale_y_continuous(expand=c(0,0), labels=scales::percent_format(scale=1)) +
      theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
  ) %>%
  pull(plot) %>%
  ggarrange(plotlist=., ncol=1, legend="right")
  

# try clustering by UniFrac distance
library(dendextend)
library(ComplexHeatmap)

phylo_tree <- Bus_CF1 %>% phy_tree()
phylo_dend <- phylo_tree  %>% spiralize::phylo_to_dendrogram() 
phylo_tree$tip.label %>% length()

is_tip <- phylo_tree$edge[,2] <= length(phylo_tree$tip.label)
ordered_tips <- phylo_tree$edge[is_tip, 2]
leaf_order <- phylo_tree$tip.label[ ordered_tips ]

stacked_barplot_data <- lapply(X, function(x) 100*x/rowSums(x)) %>%
  lapply(function(x) x %>% as.data.frame() %>%
           rownames_to_column("sample_id") %>%
           pivot_longer(-sample_id, names_to="otu", values_to="rel_abund")) %>%
  bind_rows(.id="site") %>%
  left_join(taxonomy_tbl %>% select(otu), by="otu") %>%
  left_join(sample_metadata %>% select(sample_id, Disease), by="sample_id")

colourbar_data <- sample_metadata %>%
  select(sample_id, Disease, CurrentSmoker) %>%
  mutate(
    Asthma=if_else(Disease=="Asthma", "Yes", "No")
  ) %>%
  select(-Disease)

sample_2_study_id <- setNames(study_ids$Study_ID, study_ids$sample_id)

unifrac_dist_mats <- list()

for(site_name in names(X)) {
  cat(site_name, "\n")
  
  dist_mat <- subset_samples(Bus_CF1, ID %in% study_ids$sample_id[ study_ids$site==site_name ]) %>%
    phyloseq::distance("wunifrac") 
  sample_dendro <- dist_mat %>%
    hclust() %>%
    as.dendrogram()

  sample_ordering <- sample_dendro %>% labels()
  
  # fig <- cowplot::plot_grid(
  #   ggdendro::ggdendrogram(data=sample_dendro) +
  #     scale_x_continuous(expand=expansion(add=c(0.5,0.5))) +
  #     theme_void() +
  #     theme(panel.background=element_rect(colour="white")) +
  #     ggtitle(site_name),
  #   colourbar_data %>%
  #     filter(sample_id %in% sample_ordering) %>%
  #     pivot_longer(-sample_id) %>%
  #     mutate(sample_id=factor(sample_id, levels=sample_ordering)) %>%
  #     ggplot(aes(x=sample_id, y=name, fill=value)) +
  #     geom_tile() +
  #     theme(axis.text.x=element_blank(),
  #           axis.title=element_blank(),
  #           legend.title=element_blank()) +
  #     scale_x_discrete(labels=function(x) sample_2_study_id[ x ]),
  #   stacked_barplot_data %>%
  #     filter(site==site_name) %>%
  #     left_join(taxonomy_tbl, by="otu") %>%
  #     group_by(site, sample_id, Phylum, Genus) %>%
  #     summarise(rel_abund=sum(rel_abund)) %>%
  #     ungroup() %>%
  #     mutate(sample_id=factor(sample_id, levels=sample_ordering)) %>%
  #     ggplot(aes(x=sample_id, y=rel_abund, fill=Phylum)) +
  #     geom_col(colour="black", size=0.5) +
  #     scale_y_continuous(expand=c(0,0), labels=scales::percent_format(scale=1)) +
  #     scale_x_discrete(labels=function(x) sample_2_study_id[ x ]) +
  #     theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  #     ylab("Relative abundance"),
  #   ncol=1,
  #   align="v",
  #   axis="rl",
  #   rel_heights=c(1, 0.15, 1)
  # )
  # 
  # ggsave(sprintf("../plots/cf_lll_vs_lul/dendro_and_relabund__%s.png", site_name),
  #        plot=fig,
  #        width=16, height=10, dpi=450)
  
  
}

tidy_mantel <- function(d1, d2, ...) {
  out <- vegan::mantel(d1, d2, ...)
  return(tibble(estmate=out$statistic, pvalue=out$signif,
                method=paste0("Mantel:", out$method)))
}

# just use the paired samples for the dendrograms
unifrac_site_dendros <- list()
unifrac_dist_mats <- list()

for(site_name in names(X)) {
  cat(site_name, "\n")
  
  dist_mat <- subset_samples(Bus_CF1, ID %in% paired_samples_wide[[ site_name ]]) %>%
    phyloseq::distance("wunifrac") %>%
    as.matrix()
  sample_renamer <- setNames(paired_samples_wide$Study_ID, paired_samples_wide[[ site_name ]])
  colnames(dist_mat) <- sample_renamer[ colnames(dist_mat) ] %>% unname()
  rownames(dist_mat) <- sample_renamer[ rownames(dist_mat) ] %>% unname()
  
  sample_dendro <- dist_mat %>%
    as.dist() %>%
    hclust() %>%
    as.dendrogram()
  unifrac_site_dendros[[ site_name ]] <- sample_dendro
  unifrac_dist_mats[[ site_name ]] <- dist_mat
}
as.dendlist(unifrac_site_dendros) %>%
  cor.dendlist()
as.dendlist(unifrac_site_dendros) %>%
  dist.dendlist()

unifrac_dist_mats <- lapply(
  unifrac_dist_mats,
  function(x) x[ paired_samples_wide$Study_ID, paired_samples_wide$Study_ID ]
)

tidy_mantel(unifrac_dist_mats$cf_lll, unifrac_dist_mats$cf_lul)
tidy_mantel(unifrac_dist_mats$cf_lll, unifrac_dist_mats$cf_ots)
tidy_mantel(unifrac_dist_mats$cf_lul, unifrac_dist_mats$cf_ots)

patient_clusters <- mapply(
  function(dend, dmat) dynamicTreeCut::cutreeDynamic(
    dend %>% as.hclust(),
    distM=dmat %>% as.matrix(),
    minClusterSize=7
  ),
  unifrac_site_dendros,
  unifrac_dist_mats,
  SIMPLIFY=FALSE
) # %>%
  # bind_rows(.id="site") %>%
  # pivot_wider(names_from=site, values_from=cluster)


ordered_clusters <- list()
for(site_name in names(unifrac_site_dendros)) {
  unifrac_site_dendros[[ site_name ]] %>% plot(main=site_name)
  colored_bars(patient_clusters[[ site_name ]][ order.dendrogram(unifrac_site_dendros[[ site_name ]]) ])
  ordered_clusters[[ site_name ]] <- patient_clusters[[ site_name ]][ order.dendrogram(unifrac_site_dendros[[ site_name ]]) ]
}
ordered_clusters %>%
  bind_cols() %>%
  colpair_map(function(x,y) WGCNA::randIndex(table(x, y)))

ordered_clusters %>%
  bind_cols() %>%
  colpair_map(function(x,y) WGCNA::overlapTable(x, y) %>% list())
#####################################################################################
#####################################################################################

#####################################################################################
# compare within-patient similarity to between-patient similarity
#####################################################################################


library(corrr)

long_unifrac_dists <- subset_samples(Bus_CF1, ID %in% study_ids$sample_id) %>%
  phyloseq::distance("wunifrac") %>%
  as.matrix() %>% 
  as_cordf() %>%
  shave() %>%
  stretch() %>%
  left_join(
    study_ids, by=c(x="sample_id")
  ) %>%
  left_join(
    study_ids, by=c(y="sample_id")
  )

# distances within the same patients
long_unifrac_dists %>%
  filter(Study_ID.x==Study_ID.y, x!=y) %>%
  mutate(x_label=if_else(
    site.x<site.y,
    paste0(site.x, "\n", site.y),
    paste0(site.y, "\n", site.x)
    )) %>%
  ggplot(aes(x=x_label, y=r)) +
  geom_boxplot()

# distances between patients
long_unifrac_dists %>%
  filter(Study_ID.x!=Study_ID.y, x!=y) %>%
  mutate(x_label=if_else(
    site.x<site.y,
    paste0(site.x, "\n", site.y),
    paste0(site.y, "\n", site.x)
  ))%>%
  ggplot(aes(x=x_label, y=r)) +
  geom_boxplot()

unifrac_distmat_all_samples <- subset_samples(Bus_CF1, ID %in% study_ids$sample_id) %>%
  phyloseq::distance("wunifrac")
all_sample_dend <- unifrac_distmat_all_samples %>%
  hclust() %>% 
  as.dendrogram()

all_sample_dend %>%
  plot()
leaf_order <- all_sample_dend %>% labels()
ordered_sample_info <- study_ids %>%
  arrange(match(sample_id, leaf_order))

site_cmap <- c(cf_lll="darkgreen", cf_lul="darkred", cf_ots="blue")
colored_bars(site_cmap[ ordered_sample_info$site ], dend=all_sample_dend)

#####################################################################################
#####################################################################################