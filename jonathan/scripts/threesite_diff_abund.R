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
# differential abundance analysis between asthmatics and controls
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

included_otus <- otu_summary %>%
  group_by(site) %>%
  slice_max(total_reads, n=20) %>%
  ungroup() %>%
  pull(otu) %>%
  unique()
cat(sprintf("Including %d OTUs in anaylsis\n", length(included_otus)))


wilcox_with_effsize <- function(tbl, paired=FALSE, effsize_ci=FALSE, ...) {
  
  if(!paired) {
    tbl_ <- tbl %>%
      group_by(group) %>%
      summarise(xx=list(x)) %>%
      arrange(group)
    x <- tbl_$xx[[1]]
    y <- tbl_$xx[[2]]
  } else {
    tbl_ <- tbl %>%
      pivot_wider(names_from=group, values_from=x) %>%
      relocate(id)
    stopifnot(ncol(tbl_)==3)
    x <- tbl_$site1
    y <- tbl_$site2
  }

  
  p_vals <- wilcox.test(x=x, y=y, paired=paired, ...) %>% broom::tidy()
  eff_sizes <- effectsize::rank_biserial(x=x, y=y, paired=paired, ci=effsize_ci)
  return(p_vals %>%
           bind_cols(eff_sizes))#  %>%
                      # effectsize::interpret_r(rules="cohen1988"))
  
}

get_univar_wilcox_stats <- function(data, indices) {
  otu_names <- colnames(data)[ !colnames(data) %in% c("sample_id", "group", "Study_ID")]

  data_ <- data[ indices, ]

  test_res <- lapply(
    otu_names,
    function(nm_) data_ %>%
      select(sample_id, Study_ID, group, nm_) %>%
      rename(x=nm_) %>%
      mutate(tmpcol=nm_) %>%
      wilcox_with_effsize(paired=FALSE)
  ) %>% bind_rows()
  all_pvals <- setNames(p.adjust(test_res$p.value, method="fdr"), paste0(otu_names, "__pvaladj"))
  all_effsizes <- setNames(test_res$r_rank_biserial, paste0(otu_names, "__effsize"))
  return(c(all_pvals, all_effsizes))
}



X_clr_w_cond <- lapply(
  X,
  function(x) x %>% compositions::clr() %>% as.data.frame() %>%
    rownames_to_column("sample_id") %>%
    # pivot_longer(-sample_id, names_to="otu", values_to="reads_clr") %>%
    left_join(
      sample_metadata %>%
        select(sample_id, Disease),
      by="sample_id"
    ) %>%
    select(sample_id, Disease, all_of(included_otus)) %>%
    relocate(sample_id, Disease) %>%
    as_tibble())

booted_clr_wilcox_res <- list()
for(site_name in names(X_clr_w_cond)) {
    
  cat(site_name, "\n")
  
  xx <- X_clr_w_cond[[ site_name ]]

  boot_strap <- boot(data=xx,
                     statistic=get_univar_wilcox_stats,
                     R=1000,
                     strata=xx$group %>% as.factor(),
                     parallel="multicore", ncpus=20)
  
  otu_names <- colnames(xx)[ !colnames(xx) %in% c("sample_id", "Disease")]
  colnames(boot_strap$t) <- c(paste0(otu_names, "__pvaladj"), paste0(otu_names, "__effsize"))
  
  booted_clr_wilcox_res[[ site_name ]] <- boot_strap$t %>%
    matrixStats::colQuantiles(probs=c(0.05, 0.5, 0.95), na.rm=TRUE) %>%
    as.data.frame() %>%
    rownames_to_column("labelcol") %>%
    separate(labelcol, c("otu", "name"), sep="__") %>%
    as_tibble()
}

volcano_plotdata <- booted_clr_wilcox_res %>%
  bind_rows(.id="site") %>%
  pivot_wider(names_from=name, values_from=contains("%")) %>%
  group_by(site) %>%
  mutate(effsize_rank=dense_rank(desc(abs(`50%_effsize`)))) %>%
  ungroup() %>%
  mutate(
    across(contains("pvaladj"), function(x) -log10(x)),
    label=if_else(
      `95%_effsize` < 0.00 | `5%_effsize` > 0,
      otu, NA_character_
    )
  ) %>%
  left_join(taxonomy_tbl, by="otu") 
  
volcano_plotdata %>%
  ggplot(aes(x=`50%_effsize`, y=`50%_pvaladj`, colour=Phylum)) +
  geom_point() +
  # geom_errorbar(aes(ymin=`5%_pvaladj`, ymax=`95%_pvaladj`)) +
  geom_errorbarh(data=filter(volcano_plotdata, !is.na(label)),
                 aes(xmin=`5%_effsize`, xmax=`95%_effsize`)) +
  facet_wrap(~site) +
  # scale_y_continuous(
  #   trans=scales::trans_new("-log10", function(x) -log10(x), function(x) 10**(-x))
  # ) +
  ggrepel::geom_text_repel(aes(label=label), show.legend=FALSE) +
  geom_vline(data=data.frame(x=c(0.0, 0.3, -0.3, 0.5, -0.5)), aes(xintercept=x),
             colour="black", linetype="dotted") +
  geom_hline(data=data.frame(y=c(0.1, 0.05)), aes(yintercept=-log10(y)),
             colour="black", linetype="dotted") +
  ylab("-log10(p-value)") +
  xlab("Effect size with 90% CI") +
  get_plot_theme(14) +
  theme(legend.position="top")
  
cropped_ggsave("../plots/cf_lll_vs_lul/threesite_volvanoplot_bootstrap.png",
               height=6, width=12, dpi=450)

# Glutamibacter and Haemophilius OTUs
volcano_plotdata %>%
  filter(grepl("Glut|Hae", otu)) %>%
  ggplot(aes(x=`50%_effsize`, y=`50%_pvaladj`, colour=Phylum)) +
  geom_point() +
  # geom_errorbar(aes(ymin=`5%_pvaladj`, ymax=`95%_pvaladj`)) +
  geom_errorbarh(aes(xmin=`5%_effsize`, xmax=`95%_effsize`)) +
  facet_wrap(~site) +
  # scale_y_continuous(
  #   trans=scales::trans_new("-log10", function(x) -log10(x), function(x) 10**(-x))
  # ) +
  ggrepel::geom_text_repel(aes(label=otu), show.legend=FALSE) +
  geom_vline(data=data.frame(x=c(0.0, 0.3, -0.3, 0.5, -0.5)), aes(xintercept=x),
             colour="black", linetype="dotted") +
  geom_hline(data=data.frame(y=c(0.1, 0.05)), aes(yintercept=-log10(y)),
             colour="black", linetype="dotted") +
  ylab("-log10(p-value)") +
  xlab("Effect size with 90% CI") +
  get_plot_theme(14) +
  theme(legend.position="top")

cropped_ggsave("../plots/cf_lll_vs_lul/threesite_volvanoplot_bootstrap_selected_otus.png",
               height=6, width=12, dpi=450)


#####################################################################################
#####################################################################################

#####################################################################################
# paired comparison between sites - paired analysis (n=55)
#####################################################################################

get_paired_wilcoxon_estimate <- function(data, indices) {
  data_ <- data[indices,] %>%
    mutate(id=make.unique(id)) %>%
    pivot_longer(-id, names_to="site") %>%
    unnest(value)
  
  otu_names <- colnames(data_)[ !colnames(data_) %in% c("id", "site")]
  
  test_res <- lapply(
    otu_names,
    function(nm_) data_ %>%
      select(id, site, nm_) %>%
      rename(x=nm_, group=site) %>%
      wilcox_with_effsize(paired=TRUE)
  ) %>% bind_rows()
  all_pvals <- setNames(p.adjust(test_res$p.value, method="fdr"), paste0(otu_names, "__pvaladj"))
  all_effsizes <- setNames(test_res$r_rank_biserial, paste0(otu_names, "__effsize"))
  return(c(all_pvals, all_effsizes))
}

paired_samples <- study_ids %>%
  pivot_wider(names_from=site, values_from=sample_id) %>%
  drop_na() %>%
  pivot_longer(-Study_ID, names_to="site", values_to="sample_id")

paired_samples %>%
  left_join(sample_metadata %>% select(sample_id, Disease),
            by="sample_id") %>%
  group_by(site, Disease) %>%
  tally()

X_clr_paired <- X_clr_w_cond %>%
  bind_rows(.id="site") %>%
  inner_join(paired_samples, by=c("site", "sample_id")) %>%
  relocate(sample_id, Study_ID, site) %>%
  select(-Disease)

site_pairs <- tribble(
  ~site1, ~site2,
  "cf_lll", "cf_lul",
  "cf_lll", "cf_ots",
  "cf_lul", "cf_ots"
)

booted_site_wilcox_res <- list()
for(site_pair_idx in 1:nrow(site_pairs)) {
  
  cat(site_pair_idx, "\n")
  
  xx <- X_clr_paired %>%
    filter(site %in% site_pairs[site_pair_idx,] %>% c()) %>%
    rename(group=site) %>%
    group_by(Study_ID, group) %>%
    do(otu_mat=.[,4:ncol(.)]) %>%
    pivot_wider(names_from=group, values_from=otu_mat)
  colnames(xx) <- c("id", "site1", "site2")
  
  boot_strap <- boot(data=xx,
                     statistic=get_paired_wilcoxon_estimate,
                     R=100,
                     parallel="multicore", ncpus=20)
  
  otu_names <- colnames(xx$site1[[1]])
  colnames(boot_strap$t) <- c(paste0(otu_names, "__pvaladj"), paste0(otu_names, "__effsize"))
  
  booted_site_wilcox_res[[ site_pair_idx ]] <-  boot_strap$t %>%
    matrixStats::colQuantiles(probs=c(0.05, 0.5, 0.95), na.rm=TRUE) %>%
    as.data.frame() %>%
    rownames_to_column("labelcol") %>%
    separate(labelcol, c("otu", "name"), sep="__") %>%
    as_tibble() %>%
    mutate(site1=site_pairs[site_pair_idx,1]$site1, site2=site_pairs[site_pair_idx,2]$site2)
}

site_volcano_plotdata <-  booted_site_wilcox_res %>%
  bind_rows() %>%
  pivot_wider(names_from=name, values_from=contains("%")) %>%
  mutate(facet_title=sprintf("%s--%s", site1, site2)) %>%
  group_by(facet_title) %>%
  mutate(effsize_rank=dense_rank(desc(abs(`50%_effsize`)))) %>%
  ungroup() %>%
  mutate(
    across(contains("pvaladj"), function(x) -log10(x)),
    label=if_else(
      `95%_effsize` < -0.3 | `5%_effsize` > 0.3,
      otu, NA_character_
    )
  ) %>%
  left_join(taxonomy_tbl, by="otu")

site_volcano_plotdata %>%
  group_by(facet_title) %>% 
  do(
    plot=ggplot(data=.,
                aes(x=`50%_effsize`, y=`50%_pvaladj`, colour=Phylum)) +
    geom_point() +
    # geom_errorbar(aes(ymin=`5%_pvaladj`, ymax=`95%_pvaladj`)) +
    geom_errorbarh(data=filter(., !is.na(label)),
                               aes(xmin=`5%_effsize`, xmax=`95%_effsize`)) +
    # ggrepel::geom_text_repel(aes(label=label),
    #                          show.legend=FALSE, max.overlaps=100) +
    geom_vline(data=data.frame(x=c(0.0, 0.3, -0.3, 0.5, -0.5)), aes(xintercept=x),
               colour="black", linetype="dotted") +
    geom_hline(data=data.frame(y=c(0.1, 0.05)), aes(yintercept=-log10(y)),
               colour="black", linetype="dotted") +
    ylab("-log10(p-value)") +
    xlab("Effect size with 90% CI") +
    get_plot_theme(14) +
    theme(legend.position="top") +
    ggtitle(.$facet_title) +
    coord_cartesian(xlim=c(-1, 1), ylim=c(0, 10))
  ) %>%
  pull(plot) %>%
  ggarrange(plotlist=., common.legend=TRUE, nrow=3, ncol=1)

cropped_ggsave("../plots/cf_lll_vs_lul/threesite_volvanoplot_bootstrap_sitewise_paired.png",
               height=12, width=12, dpi=450)


#####################################################################################
#####################################################################################


#####################################################################################
# check biomass
#####################################################################################


pcr_log_copyno <- sample_metadata %>%
  select(sample_id, LogqPCRCopy, Disease) %>%
  mutate(LogqPCRCopy=as.numeric(LogqPCRCopy)) 


# paired analysis
biomass_pair_data <- pcr_log_copyno %>%
  right_join(paired_samples) %>%
  select(-(sample_id)) %>%
  pivot_wider(names_from=site, values_from=LogqPCRCopy) 

p1 <- biomass_pair_data %>%
  ggpaired(
    "cf_lul", "cf_ots", id="Study_ID"
  ) 
p2 <- biomass_pair_data %>%
  ggpaired(
    "cf_lul", "cf_lll", id="Study_ID",
  ) 
p3 <- biomass_pair_data %>%
  ggpaired(
    "cf_lll", "cf_ots", id="Study_ID"
  ) 
ggarrange(plotlist=lapply(list(p1, p2, p3),
                          function(x) x +
                            ylab("log PCR copyno") +
                            xlab("Sample site") +
                            ylim(c(3, 10))
                          ),
          nrow=1)
cropped_ggsave("../plots/cf_lll_vs_lul/biomass_paired_plots.png",
               width=12, height=4, dpi=450)

#####################################################################################
#####################################################################################

#####################################################################################
# Univariate linear models for difference between (i) asthmatics and controls
# and (ii) each site (with random intercept)
#####################################################################################

library(broom)
library(multidplyr)

between_eps <- function(lower, upper, eps) {
  stopifnot(length(lower)==length(upper))
  out <- rep(FALSE, length(lower))
  out[ lower > eps | upper < -eps ] <- TRUE
  return(out)
}


# collect bootstrap estimates
collect_bootreps <- function(tbl_list, col_names, probs) {
  booted_estimates <- lapply(tbl_list, function(x) x$estimate) %>%
    do.call(rbind, .) %>%
    as.matrix() 
  
  for(col_name_ in col_names) {
    sapply(
      tbl_list,
      function(x) identical(tbl_list[[1]][[col_name_]], x[[col_name_]])
    ) %>% all() %>% stopifnot()
  }
  
  return(
    tbl_list[[1]] %>%
      select(all_of(col_names)) %>%
      bind_cols(
        matrixStats::colQuantiles(booted_estimates, probs=probs, na.rm=TRUE)
      )
  )
}

included_covariates <- c("age", "Disease", "CurrentSmoker", "qPCRCN")

otu_summary <- lapply(X, compute_summary_stats) %>%
  bind_rows(.id="site") %>%
  left_join(taxonomy_tbl, by="otu")

included_otus <-  otu_summary %>%
  group_by(site) %>%
  slice_max(total_reads, n=100) %>%
  ungroup() %>%
  pull(otu) %>%
  unique()
cat(sprintf("Including %d OTUs in anaylsis\n", length(included_otus)))

# CLR normlisation
X_ols_clr <- lapply(
  X,
  function(x) x %>%
    compositions::clr() %>%
    as.data.frame() %>%
    rownames_to_column("sample_id") %>%
    as_tibble()
)

# DeSeq2 normalisation
deseq_sample_list <- lapply(
  X,
  function(x) tibble(sample_id=rownames(x)) %>%
    left_join(sample_metadata) %>%
    select(sample_id, Disease) %>%
    group_by(Disease) %>%
    summarise(x=list(unique(sample_id)))
) %>%
  lapply(function(x) setNames(x[[2]], x[[1]]))

deseq_X <- lapply(X, function(x) x %>% as.data.frame() %>% t())

X_deseq_norm <- mapply(
  function(x, s_list) metaseqR2::normalizeDeseq2(x[,do.call(c, s_list)], s_list),
  deseq_X,
  deseq_sample_list
) %>%
  lapply(function(x) x %>%
           t() %>%
           as.data.frame() %>%
           rownames_to_column("sample_id") %>%
           as_tibble())

list(clr=X_ols_clr, deseq=X_deseq_norm) %>%
  unlist(recursive=FALSE)


ols_model_data <- list(clr=X_ols_clr, deseq=X_deseq_norm) %>%
  unlist(recursive=FALSE) %>%
  bind_rows(.id="label") %>%
  separate(label, c("normalisation", "site"), sep="\\.") %>%
  left_join(sample_metadata %>%
              select(sample_id, all_of(included_covariates)),
            by="sample_id"
  ) %>%
  pivot_longer(-c(sample_id, normalisation, site, all_of(included_covariates)), names_to="otu", values_to="y") %>%
  filter(otu %in% included_otus) %>%
  mutate(qPCRCN=log(qPCRCN))

ols_model_data %>%
  filter(normalisation=="deseq") %>%
  ggplot(aes(x=y, fill=Disease)) +
  geom_histogram() +
  facet_wrap(~otu)

univar_lm_fits <- ols_model_data %>%
  filter(otu %in% included_otus) %>%
  group_by(normalisation, site, otu) %>%
  do(
    lm(
      as.formula(sprintf("y~%s", paste0(included_covariates, collapse="+"))),
      data=.
    ) %>% tidy()
  )

# univar_lm_fits %>%
#   select(normalisation, site, otu, term, estimate) %>%
#   pivot_wider(names_from=normalisation, values_from=estimate) %>%
#   filter(term=="DiseaseControl") %>%
#   ggplot(aes(x=clr, y=deseq)) +
#   geom_point() +
#   facet_wrap(~site) +
#   stat_cor(method="pearson") +
#   coord_cartesian(ylim=c(-500,500))

booted_univar_lm_fits <- pbmcapply::pbmclapply(
  1:100,
  function(x) ols_model_data %>%
    group_by(normalisation, site, otu, Disease, CurrentSmoker) %>%
    slice_sample(prop=1, replace=TRUE) %>% 
    group_by(normalisation, site, otu) %>%
    do(
      lm(
        as.formula(sprintf("y~%s", paste0(included_covariates, collapse="+"))),
        data=.
      ) %>% tidy()
    ),
  mc.cores=20
)

sapply(booted_univar_lm_fits, function(x) !is.null(x)) %>%
  all() %>%
  stopifnot()

ols_model_results <- collect_bootreps(
  booted_univar_lm_fits, c("normalisation", "site", "term", "otu"), probs=c(0.05, 0.5, 0.95))

ols_model_results %>%
  filter(term=="qPCRCN") %>%
  arrange(desc(abs(`50%`)))

ols_model_results %>%
  filter(
    term=="DiseaseControl",
    # !data.table::between(0, `5%`, `95%`)
    # between_eps(`5%`, `95%`, 0.05)
  ) %>%
  left_join(taxonomy_tbl %>%
              select(otu, Phylum)) %>%
  ggplot(aes(x=reorder(otu, `50%`), y=`50%`, colour=Phylum)) +
  geom_point() +
  geom_errorbar(aes(ymin=`5%`, ymax=`95%`)) +
  theme(axis.text.x=element_text(size=8, angle=90, hjust=1, vjust=0.5),
        legend.position="top") +
  xlab("OTU") +
  ylab("Asthma effect size\n(negative direction -> higher in asthmatics)") +
  geom_hline(yintercept=0, colour="red", linetype="dotted") +
  facet_grid(cols=vars(site), rows=vars(normalisation), scales="free_y")
  
ols_model_results %>%
  filter(
    term=="CurrentSmokerYes",
    !data.table::between(0, `5%`, `95%`)
    # between_eps(`5%`, `95%`, 0.05)
  ) %>%
  left_join(taxonomy_tbl %>%
              select(otu, Phylum)) %>%
  ggplot(aes(x=reorder(otu, `50%`), y=`50%`, colour=Phylum)) +
  geom_point() +
  geom_errorbar(aes(ymin=`5%`, ymax=`95%`)) +
  theme(axis.text.x=element_text(size=8, angle=90, hjust=1, vjust=0.5),
        legend.position="top") +
  xlab("OTU") +
  ylab("CurrentSmoker effect size") +
  geom_hline(yintercept=0, colour="red", linetype="dotted") +
  ggforce::facet_row(~site, scales="free_x", space="free")

# agreement between deseq and clr normalisation
plt <- ols_model_results %>%
  filter(
    term=="DiseaseControl"
  ) %>%
  pivot_wider(names_from=normalisation, values_from=contains("%")) %>%
  left_join(taxonomy_tbl %>%
              select(otu, Phylum)) %>%
  ggplot(aes(x=`50%_clr`, y=`50%_deseq`, colour=Phylum)) +
  geom_point(size=3) +
  # geom_errorbar(aes(ymin=`5%_deseq`, ymax=`95%_deseq`)) +
  # geom_errorbarh(aes(xmin=`5%_clr`, xmax=`95%_clr`)) +
  theme(axis.text=element_text(size=12),
        legend.position="top") +
  xlab("CLR normalisation") +
  ylab("DESeq2 normalisation") +
  geom_hline(yintercept=0, colour="red", linetype="dotted", size=1) +
  geom_vline(xintercept=0, colour="red", linetype="dotted", size=1) +
  facet_wrap(~site) +
  stat_cor(method="spearman", aes(x=`50%_clr`, y=`50%_deseq`), size=5, inherit.aes=FALSE) #  +
  # stat_smooth(method="lm", aes(x=`50%_clr`, y=`50%_deseq`), inherit.aes=FALSE)

plt + 
  coord_cartesian(ylim=c(-250, 250))

#
# site by site analysis - mixed model as repeated samples from individuals
library(lme4)
library(broom.mixed)
library(lmeresampler)

do_lmerboot <- function(mod, B, level, ...) {
  boot_out <- lmeresampler::bootstrap(model=mod, B=B, type="case", .f=fixef, resample=c(TRUE, FALSE))
  # return(boot_out$stats)
  return(confint(boot_out, level=level) %>% filter(type=="perc"))
}

lmm_data <- ols_model_data %>%
  left_join(study_ids, by=c("site", "sample_id")) %>%
  rename(id=Study_ID) %>%
  select(-sample_id) %>%
  mutate(age=scale(age)[,1], qPCRCN=scale(qPCRCN)[,1],
         id=as.factor(id))

all_lmm_mods <- lmm_data %>%
  group_by(normalisation, otu) %>%
  summarise(
    mod=lmer(
      as.formula(sprintf("y~site+%s + (1 | id)", paste0(included_covariates, collapse="+"))),
      control=lmerControl(optimizer="bobyqa")
    ) %>% list()
  ) %>%
  ungroup() %>%
  mutate(is_singular=sapply(mod, isSingular))

lmm_booted_fixefs <- pbmcapply::pbmclapply(
  setNames(all_lmm_mods$mod[ !all_lmm_mods$is_singular ], all_lmm_mods$otu[ !all_lmm_mods$is_singular ]),
  function(m) do_lmerboot(m, B=100, level=0.9),
  mc.cores=20
) %>%
  as_tibble_col() %>%
  bind_cols(all_lmm_mods %>% filter(!is_singular))

lmm_booted_fixefs %>%
  select(-mod) %>%
  unnest(value) %>%
  filter(grepl("sitecf_", term), normalisation=="clr") %>%
  mutate(
    estimate=exp(estimate),
    lower=exp(lower),
    upper=exp(upper)
  ) %>%
  left_join(
    taxonomy_tbl %>%
      select(otu, Phylum),
    by="otu"
  ) %>%
  mutate(Phylum=as.factor(Phylum)) %>%
  arrange(Phylum) %>%
  mutate(otu=factor(otu, levels=unique(otu))) %>%
  ggplot(aes(x=otu, y=estimate, colour=Phylum)) +
  geom_point() + 
  geom_errorbar(aes(ymin=lower, ymax=upper)) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  geom_hline(yintercept=1, colour="red", linetype="dotted") +
  ylab("Odds ratio (relative to LLL, 90% CI)") +
  facet_grid(rows=vars(term))


























# stratified sample of individuals
booted_lmm_params <- lmm_data %>%
  select(Disease, CurrentSmoker, id) %>%
  distinct() %>%
  group_by(Disease, CurrentSmoker) %>%
  slice_sample(prop=1, replace=TRUE) %>%
  mutate(unique_id=make.unique(id)) %>%
  left_join(lmm_data) %>%
  group_by(otu) %>%
  do(
    lmer(
      as.formula(sprintf("y~site+%s + (1 | unique_id) +  (Disease | unique_id)", paste0(included_covariates, collapse="+"))),
      data=.
    ) %>% tidy()
  )

booted_univar_lm_fits_site <- pbmcapply::pbmclapply(
  1:10,
  function(...) lmm_data %>%
    select(Disease, CurrentSmoker, id) %>%
    distinct() %>%
    group_by(Disease, CurrentSmoker) %>%
    slice_sample(prop=1, replace=TRUE) %>%
    mutate(unique_id=make.unique(id)) %>%
    left_join(lmm_data) %>%
    group_by(otu) %>%
    do(
      lmer(
        as.formula(sprintf("y~site+%s + (1 | unique_id) +  (Disease | unique_id)", paste0(included_covariates, collapse="+"))),
        data=.
      ) %>% tidy()
    ),
  mc.cores=10
)

ols_sitemodel_results <- collect_bootreps(
  booted_univar_lm_fits_site,
  c("term", "otu"), probs=c(0.05, 0.5, 0.95))

pd <- position_dodge(0.8)
ols_sitemodel_results %>%
  filter(
    grepl("site", term),
    !data.table::between(0, `5%`, `95%`)
  ) %>%
  mutate(
    term=str_remove(term, "site") %>%
      sprintf("%s / cf_lll", .)
  ) %>%
  filter(term=="cf_ots / cf_lll") %>%
  left_join(taxonomy_tbl %>%
              select(otu, Phylum)) %>%
  ggplot(aes(x=reorder(otu, `50%`), y=`50%`, colour=term)) +
  geom_point(position=pd) +
  geom_errorbar(aes(ymin=`5%`, ymax=`95%`), position=pd) +
  theme(axis.text.x=element_text(size=8, angle=90, hjust=1, vjust=0.5)) +
  xlab("OTU") +
  ylab("Site effect size") +
  geom_hline(yintercept=0, colour="red", linetype="dotted") +
  ggforce::facet_row(~Phylum, scales="free_x", space="free") +
  theme(legend.position="top")
