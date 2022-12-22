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




#
# prevalence at each site
plot_tbl <- otu_summary %>%
  select(site, otu, prevalence, Phylum) %>%
  gather_pairs(key=site, value=prevalence) %>%
  do(
    plot=ggplot(
      data=., aes(x=.x, y=.y, colour=Phylum)
    ) +
      geom_point(size=2) +
      xlab(.$.col[[1]]) + ylab(.$.row[[1]]) +
      guides(colour=guide_legend(ncol=2)) +
      get_plot_theme(18) +
      scale_x_continuous(labels=scales::percent_format(scale=1)) +
      scale_y_continuous(labels=scales::percent_format(scale=1)) +
      tune::coord_obs_pred() +
      geom_abline() +
      guides(colour=guide_legend(ncol=4))
  )
ggarrange(
  plotlist=plot_tbl$plot, common.legend=TRUE, nrow=1,
  legend="top"
)
cropped_ggsave("../plots/cf_lll_vs_lul/three_site_prevalence.pdf",
               height=8, width=18)

#
# reads per sample at each site
plot_tbl <- otu_summary %>%
  select(site, otu, total_reads_per_sample, Phylum) %>%
  gather_pairs(key=site, value=total_reads_per_sample) %>%
  do(
    plot=ggplot(
      data=., aes(x=.x, y=.y, colour=Phylum)
    ) +
      geom_point(size=2) +
      xlab(.$.col[[1]]) + ylab(.$.row[[1]]) +
      guides(colour=guide_legend(ncol=2)) +
      get_plot_theme(16) +
      scale_x_log10(limits=c(1e-3,1e5)) +
      scale_y_log10(limits=c(1e-3,1e5)) +
      geom_abline() +
      guides(colour=guide_legend(ncol=4))
  )
plot_tbl$plot[[3]]
ggarrange(
  plotlist=plot_tbl$plot, common.legend=TRUE, nrow=1,
  legend="top"
)
cropped_ggsave("../plots/cf_lll_vs_lul/three_site_read_per_sample.pdf",
               height=6, width=16)
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
        get_plot_theme(10) +
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