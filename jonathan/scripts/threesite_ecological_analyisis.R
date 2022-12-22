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
library(boot)

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
# site vs site ecological analysis
#####################################################################################
library(boot)
library(broom)

alpha_metrics <- subset_samples(
  Bus_CF1, ID %in% all_sample_ids$sample_id
) %>%
  microbiome::alpha() %>%
  rownames_to_column("sample_id") %>% 
  as_tibble() %>%
  left_join(all_sample_ids, by="sample_id") %>%
  left_join(sample_metadata %>%
              select(sample_id, Disease),
            by="sample_id")

alpha_metrics_slim <- alpha_metrics %>%
  select(sample_id, site, Disease,
         observed, diversity_inverse_simpson, evenness_pielou, dominance_relative) %>%
  pivot_longer(-c(sample_id, Disease, site))
booted_ecological_metrics <- alpha_metrics_slim %>%
  group_by(site, Disease, name) %>%
  do(
    boot(.$value, statistic=function(data, indices) mean(data[indices]), R=1000) %>%
      tidy(conf.int=TRUE, conf.level=0.95)
  )
booted_ecological_metrics %>%
  ggplot(aes(y=site, x=statistic, colour=Disease)) +
  geom_point(size=2) +
  geom_errorbarh(aes(xmin=conf.low, xmax=conf.high)) +
  ggforce::facet_row(~name, scales="free_x") +
  get_plot_theme(16)

# mean_sitediff_boot <- function(data, indices) {
#   tmpdf <- data[indices,] %>%
#     group_by(site) %>%
#     summarise(value=mean(value))
#   site_values <- setNames(tmpdf$value, tmpdf$site)
#   return(
#     c(site_values[["cf_lll"]]-site_values[["cf_lul"]],
#       site_values[["cf_lll"]]-site_values[["cf_ots"]],
#       site_values[["cf_lul"]]-site_values[["cf_ots"]])
#   )
# }
# 
# booted_ecological_metrics2 <- alpha_metrics_slim %>%
#   group_by(Disease, name) %>%
#   do(
#     boot(., statistic=mean_sitediff_boot, R=1000) %>%
#       tidy(conf.int=TRUE, conf.level=0.90)
#   ) %>%
#   mutate(
#     site_pair=c("cf_lll-cf_lul", "cf_lll-cf_ots", "cf_lul-cf_ots")
#   )
# pd <- position_dodge(0.75)
# booted_ecological_metrics2 %>%
#   ungroup() %>%
#   mutate(
#     site_pair=str_remove_all(
#       site_pair, "cf_"
#     ) %>%
#       toupper() %>%
#       str_replace("-", " - "),
#     name=factor(
#       name,
#       levels=c("observed", "evenness_pielou", "diversity_inverse_simpson", "dominance_relative")
#     )
#   ) %>%
#   ggplot(aes(y=site_pair, x=statistic, colour=Disease)) +
#   geom_point(size=2, position=pd) +
#   geom_errorbarh(aes(xmin=conf.low, xmax=conf.high),  position=pd) +
#   facet_wrap(~name, scales="free_x") +
#   get_plot_theme(16) +
#   theme(panel.spacing.x=unit(2, "lines"),
#         legend.position="top") +
#   geom_vline(xintercept=0, colour="red", linetype="dotted", size=1) +
#   ylab("Site pair") +
#   xlab("Diversity measure value with 90% CI")
# cropped_ggsave(
#   "../plots/cf_lll_vs_lul/pairwise_diversity_differences.pdf",
#   height=9, width=9
# )

# just paired samlpes
mean_sitediff_boot_paired <- function(data, indices) {
  tmpdf <- data[indices,] %>%
    select(-c(Disease, name)) %>%
    pivot_longer(-Study_ID, names_to="site") %>%
    group_by(site) %>%
    summarise(value=mean(value))
  site_values <- setNames(tmpdf$value, tmpdf$site)
  return(
    c(site_values[["cf_lll"]]-site_values[["cf_lul"]],
      site_values[["cf_lll"]]-site_values[["cf_ots"]],
      site_values[["cf_lul"]]-site_values[["cf_ots"]])
  )
}

paired_alpha_metrics <- alpha_metrics_slim %>%
  left_join(study_ids, by=c("sample_id", "site")) %>%
  select(-sample_id) %>%
  pivot_wider(names_from=site) %>%
  drop_na()

booted_ecological_metrics3 <- paired_alpha_metrics %>%
  group_by(Disease, name) %>%
  do(
    boot(., statistic=mean_sitediff_boot_paired, R=100) %>%
      tidy(conf.int=TRUE, conf.level=0.90)
  ) %>%
  mutate(
    site_pair=c("cf_lll-cf_lul", "cf_lll-cf_ots", "cf_lul-cf_ots")
  )
pd <- position_dodge(0.75)
booted_ecological_metrics3 %>%
  ungroup() %>%
  mutate(
    site_pair=str_remove_all(
      site_pair, "cf_"
    ) %>%
      toupper() %>%
      str_replace("-", " - "),
    name=factor(
      name,
      levels=c("observed", "evenness_pielou", "diversity_inverse_simpson", "dominance_relative")
    )
  ) %>%
  ggplot(aes(y=site_pair, x=statistic, colour=Disease)) +
  geom_point(size=2, position=pd) +
  geom_errorbarh(aes(xmin=conf.low, xmax=conf.high),  position=pd) +
  facet_wrap(~name, scales="free_x") +
  get_plot_theme(16) +
  theme(panel.spacing.x=unit(2, "lines"),
        legend.position="top") +
  geom_vline(xintercept=0, colour="red", linetype="dotted", size=1) +
  ylab("Site pair") +
  xlab("Diversity measure value with 90% CI")
cropped_ggsave(
  "../plots/cf_lll_vs_lul/pairwise_diversity_differences.pdf",
  height=9, width=9
)

#
# difference between asthmatics and controls
get_div_diff <- function(data, indices, cor.type="pearson") {
  
  tmp <- data[indices,] %>%
    group_by(Disease) %>%
    summarise(value=mean(value)) %>%
    arrange(Disease)
  return(tmp$value[[1]] - tmp$value[[2]])
  
}

botoed_alpha_diffs <- alpha_metrics_slim %>%
  group_by(site, name) %>%
  do(
    boot(data=., statistic=get_div_diff, R=1000) %>%
      tidy(conf.int=TRUE, conf.level=0.95)
  )

botoed_alpha_diffs %>%
  mutate(name=factor(
    name,
    levels=c("observed", "evenness_pielou", "diversity_inverse_simpson", "dominance_relative")
  )) %>%
  ggplot(aes(y=site, x=statistic)) +
  geom_point(size=2) +
  geom_errorbarh(aes(xmin=conf.low, xmax=conf.high)) +
  ggforce::facet_row(~name, scales="free_x") +
  get_plot_theme(16) +
  geom_vline(xintercept=0, colour="red", linetype="dotted", size=1) +
  xlab("Asthmatics - Controls") +
  ylab("Site") +
  theme(panel.spacing.x=unit(2, "lines"),
        legend.position="top") 
cropped_ggsave("../plots/cf_lll_vs_lul/asthma_vs_control_div_boot.pdf",
               height=5, width=13)
# heatmaps comparing correlation and stability
sample_groups <- sample_metadata %>% select(sample_id, Disease)

#
# site-wise
booted_eco_diff <- list()
for(site_name in names(X_subset)) {
  cat(site_name, "\n")
  
  x <- alpha_metrics_slim %>%
    filter(site==site_name) %>%
    select(-site)
  
  boot_strap <- boot(data=x,
                     strata=x$Disease,
                     statistic=get_div_diff, R=10, parallel="multicore", ncpus=10)
  n_nas <- sum(is.na(boot_strap$t))
  cat(sprintf("There are %d NAs in bootstrapped estimates\n", n_nas))
  
  booted_diffcorr[[ site_name ]] <- boot_strap$t %>%
    matrixStats::colQuantiles(probs=c(0.05, 0.5, 0.95), na.rm=TRUE) %>%
    as_tibble() %>%
    bind_cols(x %>% select(-group) %>% correlate(quiet=TRUE) %>% stretch() %>% select(x,y)) %>%
    relocate(x,y)
}
