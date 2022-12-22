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
library(data.table)

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

#####################################################################################
# various distances between sites for asthmatics and controls
#####################################################################################

library(corrr)
library(car)
library(broom)

paired_samples <- study_ids %>%
  pivot_wider(names_from=site, values_from=sample_id) %>%
  drop_na() %>%
  pivot_longer(-Study_ID, names_to="site", values_to="sample_id") %>%
  left_join(sample_metadata %>% select(sample_id, Disease),
            by="sample_id")
paired_samples %>%
  group_by(site, Disease) %>%
  tally()

# compute all required beta diversity distance matrices
paired_phyloseq <- subset_samples(Bus_CF1, ID %in% paired_samples$sample_id)
dist_mats <- list()
dist_mats[["w_unifrac"]] <-paired_phyloseq %>%
  phyloseq::distance("wunifrac") %>%
  as.matrix()
dist_mats[["uw_unifrac"]] <- paired_phyloseq %>%
  phyloseq::distance("unifrac") %>%
  as.matrix()
dist_mats[["jaccard"]] <- vegan::vegdist(paired_phyloseq %>% otu_table() %>% t(), method="jaccard") %>%
  as.matrix()
dist_mats[["bray"]] <- vegan::vegdist(paired_phyloseq %>% otu_table() %>% t(), method="bray") %>%
  as.matrix()

long_dist_mats <- lapply(
  dist_mats,
  function(x) as_cordf(x) %>%
    shave() %>%
    stretch(na.rm=TRUE) %>%
    left_join(paired_samples, by=c(x="sample_id")) %>%
    left_join(paired_samples, by=c(y="sample_id")))

paired_distance_data <- long_dist_mats %>%
  bind_rows(.id="distance") %>%
  mutate(site_pair=make_unique_key(site.x, site.y)) %>%
  filter(Study_ID.x==Study_ID.y)
stopifnot(all(paired_distance_data$Disease.x==paired_distance_data$Disease.y))

paired_distance_data %>%
  ggplot(aes(x=site_pair, colour=Disease.x, y=r)) +
  geom_boxplot(notch=FALSE, outlier.shape=NA) +
  geom_jitter(position=position_jitterdodge()) +
  facet_wrap(~distance)

samples_for_boot <- paired_distance_data %>%
  select(Study_ID.x, Disease.x) %>%
  distinct()

# bootstrap levine effect size
levene_w_effsize <- function(y, group, center=median, effsize=TRUE, ...) {
  valid <- complete.cases(y, group)
  meds <- tapply(y[valid], group[valid], center, ...)
  resp <- abs(y - meds[group])
  anova_res <- anova(lm(resp ~ group))
  table <- anova_res[, c(1, 4, 5)]
  results <- broom::tidy(table) %>%
    filter(term=="group") 
  if(effsize) {
    effsize_tbl <- effectsize::eta_squared(anova_res, partial=FALSE, ci=FALSE) %>%
      as_tibble() %>%
      rename(term=Parameter)
    results <- results %>% 
      left_join(effsize_tbl, by="term")
  }
  return(results)
}

ci_levels <- c(0.05, 0.5, 0.95)
booted_levine_test_results <- pbmcapply::pbmclapply(
  1:1000,
  function(...) samples_for_boot %>%
    group_by(Disease.x) %>%
    slice_sample(prop=1, replace=TRUE) %>%
    ungroup() %>%
    left_join(
      paired_distance_data,
      by=c("Study_ID.x", "Disease.x")
    ) %>%
    mutate(Disease.x=as.factor(Disease.x)) %>%
    group_by(distance, site_pair) %>%
    do(
      levene_w_effsize(.$r, .$Disease.x, center=median)
    ),
    mc.cores=10
) %>% rbindlist(idcol="resample_idx") %>%
  group_by(distance, site_pair) %>%
  summarise(
    statistic_median=quantile(statistic, ci_levels[[2]]),
    eta2_median=quantile(Eta2, ci_levels[[2]]),
    statistic_lower=quantile(statistic, ci_levels[[1]]),
    eta2_lower=quantile(Eta2, ci_levels[[1]]),
    statistic_upper=quantile(statistic, ci_levels[[3]]),
    eta2_upper=quantile(Eta2, ci_levels[[3]])
  ) %>%
  ungroup()

# permutation test for p-values
obs_levine_test_res <- paired_distance_data %>%
  mutate(Disease.x=as.factor(Disease.x)) %>%
  group_by(distance, site_pair) %>%
  do(
    levene_w_effsize(.$r, .$Disease.x, center=median)
  )

permuted_levine_test_results <- pbmcapply::pbmclapply(
  1:1000,
  function(...) samples_for_boot %>%
    mutate(Disease.x=sample(Disease.x),
           Disease.y=Disease.x) %>%
    ungroup() %>%
    left_join(
      paired_distance_data %>%
        select(-c(Disease.x, Disease.y)),
      by=c("Study_ID.x")
    ) %>%
    mutate(Disease.x=as.factor(Disease.x)) %>%
    group_by(distance, site_pair) %>%
    do(
      levene_test(data=., formula=r~Disease.x, center=median)
    ),
  mc.cores=20
) %>%
  rbindlist(idcol="resample_idx") %>%
  select(resample_idx, distance, site_pair, statistic) %>%
  left_join(obs_levine_test_res %>%
              rename(obs_statistic=statistic) %>%
              select(distance, site_pair, obs_statistic)) %>%
  group_by(distance, site_pair) %>%
  summarise(
    p_value=(1+sum(statistic>obs_statistic))/(1+n())
  ) %>%
  ungroup()

permuted_levine_test_results %>%
  filter(distance=="w_unifrac") %>%
  mutate(p_adj=p.adjust(p_value, method="fdr")) %>%
  select(-p_value) %>%
  pivot_wider(names_from=distance, values_from=p_adj)

booted_levine_test_results %>%
  filter(distance=="w_unifrac") %>%
  ggplot(aes(y=site_pair, x=eta2_median)) +
  geom_point() +
  geom_errorbarh(aes(xmin=eta2_lower, xmax=eta2_upper)) + 
  facet_wrap(~distance) +
  xlab("eta2 (90% CI)") +
  geom_vline(xintercept=0.01, colour="red") +
  geom_vline(xintercept=0.06, colour="orange") +
  geom_vline(xintercept=0.14, colour="darkgreen")

#####################################################################################
#####################################################################################

#####################################################################################
# PCoA on weighted UniFrac distance
#####################################################################################

get_biplot_data <- function(x, plot.axes=c(1,2)) {
  pr.coo <- x$vectors
  diag.dir <- diag(c(1,1))
  pr.coo[,plot.axes] <- pr.coo[,plot.axes] %*% diag.dir
  return(pr.coo[,plot.axes] %>%
           as.data.frame() %>%
           rownames_to_column("sample_id") %>%
           as_tibble())
}

pcoa_res <- ape::pcoa(dist_mats$w_unifrac)

pcoa_res %>%
  get_biplot_data() %>%
  left_join(paired_samples, by="sample_id") %>%
  mutate(facet_title=sprintf("%s (%s)", Study_ID, Disease)) %>%
  ggplot(aes(x=Axis.1,
             y=Axis.2,
             colour=site)) +
  geom_text(aes(label=Study_ID)) +
  # geom_path(aes(group=Study_ID), colour="black") +
  # stat_ellipse() +
  # ggforce::geom_mark_hull(aes(group=Study_ID)) +
  facet_wrap(~facet_title)

#####################################################################################
#####################################################################################

#####################################################################################
# test on distances using linear model
#####################################################################################


# compute all required beta diversity distance matrices
slim_phyloseq <- subset_samples(Bus_CF1, ID %in% study_ids$sample_id)

dist_mats <- list()
dist_mats[["w_unifrac"]] <- slim_phyloseq %>%
  phyloseq::distance("wunifrac") %>%
  as.matrix()
dist_mats[["uw_unifrac"]] <- slim_phyloseq %>%
  phyloseq::distance("unifrac") %>%
  as.matrix()
dist_mats[["jaccard"]] <- vegan::vegdist(slim_phyloseq %>% otu_table() %>% t(), method="jaccard") %>%
  as.matrix()
dist_mats[["bray"]] <- vegan::vegdist(slim_phyloseq %>% otu_table() %>% t(), method="bray") %>%
  as.matrix()


long_dist_mats <- lapply(
  dist_mats,
  function(x) as_cordf(x) %>%
    shave() %>%
    stretch(na.rm=TRUE) %>%
    left_join(study_ids, by=c(x="sample_id")) %>%
    left_join(study_ids, by=c(y="sample_id")))

within_patient_dists <- long_dist_mats %>%
  rbindlist(idcol="distance") %>%
  as_tibble() %>%
  filter(Study_ID.x==Study_ID.y) %>%
  mutate(site_pair=make_unique_key(site.x, site.y)) %>%
  rename(Study_ID=Study_ID.x) %>%
  select(-c(site.x, site.y, Study_ID.y, x, y))

patient_metadata <- sample_metadata %>%
  select(sample_id, age, sex, Disease, BMI, CurrentSmoker) %>%
  left_join(study_ids) %>%
  select(-c(sample_id, site)) %>%
  distinct()

# fit models
library(lme4)
library(lmeresampler)
library(broom)
library(broom.mixed)
library(purrr)

# response: LLL - OTS distance
fit_model <- partial(lm, formula=y ~ age + sex + BMI + Disease + CurrentSmoker + 1)
mod1_data <-  within_patient_dists %>%
  filter(site_pair=="cf_lll--cf_ots") %>%
  left_join(patient_metadata) %>%
  rename(y=r)
mod1_data %>%
  ggplot(aes(x=y)) +
  geom_histogram()
lm_res_obs <- mod1_data %>%
  group_by(distance) %>%
  do(
    lm_out=fit_model(data=.) %>% lm.beta::lm.beta() %>% tidy()
  ) %>%
  ungroup() %>%
  unnest(lm_out)

## bootstrap
bootstrap_mod1_params <- pbmcapply::pbmclapply(
  1:1000,
  function(...) mod1_data %>%
  group_by(distance, sex, Disease) %>%
  slice_sample(prop=1, replace=TRUE) %>%
  group_by(distance) %>%
  do(
    lm_out=fit_model(data=.) %>% lm.beta::lm.beta() %>% tidy()
  ) %>%
  unnest(lm_out),
  mc.cores=20
) %>%
  rbindlist(idcol="resample_idx") %>%
  as_tibble() %>%
  group_by(distance, term) %>%
  summarise(
    estimate_median=quantile(estimate, 0.5),
    estimate_lower=quantile(estimate, 0.05),
    estimate_upper=quantile(estimate, 0.95)
  ) 

pd <- position_dodge(width=0.7)
bootstrap_mod1_params %>%
  filter(term!="(Intercept)") %>%
  ggplot(aes(x=estimate_median, y=term, colour=distance)) +
  geom_point(position=pd) +
  geom_errorbarh(aes(xmin=estimate_lower, xmax=estimate_upper),
                 position=pd) +
  geom_vline(xintercept=0, linetype="dotted", colour="red") +
  xlab("Standardised effect size (90% CI)") +
  ggtitle("response:  d(LLL,OTS)")


#
# ratio of LLL-OTS to LLL-LUL
mod2_data <- within_patient_dists %>%
  left_join(patient_metadata) %>%
  pivot_wider(names_from=site_pair, values_from=r) %>%
  drop_na() %>%
  mutate(y=log(`cf_lll--cf_ots`/`cf_lll--cf_lul`)) %>%
  select(-contains("--"))

mod2_data %>%
  ggplot(aes(x=y)) +
  geom_histogram()

lm2_res_obs <- mod2_data %>%
  group_by(distance) %>%
  do(
    lm_out=fit_model(data=.) %>% lm.beta::lm.beta() %>% tidy()
  ) %>%
  ungroup() %>%
  unnest(lm_out)


bootstrap_mod2_params <- pbmcapply::pbmclapply(
  1:1000,
  function(...) mod2_data %>%
    group_by(distance, sex, Disease) %>%
    slice_sample(prop=1, replace=TRUE) %>%
    group_by(distance) %>%
    do(
      lm_out=fit_model(data=.) %>% lm.beta::lm.beta() %>% tidy()
    ) %>%
    unnest(lm_out),
  mc.cores=10
) 

bootstrap_mod2_params[ sapply(bootstrap_mod2_params, is_tibble) ]%>%
  rbindlist(idcol="resample_idx") %>%
  as_tibble() %>%
  group_by(distance, term) %>%
  summarise(
    estimate_median=quantile(estimate, 0.5),
    estimate_lower=quantile(estimate, 0.05),
    estimate_upper=quantile(estimate, 0.95)
  ) %>%
  filter(term!="(Intercept)") %>%
  ggplot(aes(x=estimate_median, y=term, colour=distance)) +
  geom_point(position=pd) +
  geom_errorbarh(aes(xmin=estimate_lower, xmax=estimate_upper),
                 position=pd) +
  geom_vline(xintercept=0, linetype="dotted", colour="red") +
  xlab("Standardised effect size (90% CI)") +
  ggtitle(expression(paste("response:  log ", frac("d(LLL,OTS)", "d(LLL,LUL)"))))
