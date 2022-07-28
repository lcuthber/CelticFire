library(dplyr)
library(data.table)
library(tidyr)
library(tibble)
library(forcats)
library(broom)
library(stringr)

library(ggplot2)
theme_set(theme_minimal(base_size=16))
library(ggforce)

source("../scripts/kegg_utils.R")
source("../scripts/ko_lookup_utils.R")

##############################################################################
# Compare WGNCNA modules to isolate clusters using beta regression
##############################################################################

# load results from KEGG analysis
emapper <- fread("../data/kegg/formatted/binary_emapper_matrix_no_duplicates.csv") %>%
  column_to_rownames("isolate")
isolate_clusters <- fread("../results/kegg/isolate_cluster_labels.csv")
cluster_names <- fread("../results/kegg/isolate_cluster_names.csv")
duplicated_KOs <- fread("../data/kegg/formatted/duplicated_KOs.csv") %>%
  mutate(identical_to=str_split(identical_to, "\\|")) %>%
  unnest(identical_to)

# Load required WGCNA results
module_assignments <- readRDS("../results/wgcna/module_assignments.rds")
hub_otu_labels <- readRDS("../results/wgcna/hub_otu_labels.rds")
module_eigenvectors <- readRDS("../results/wgcna/module_eigenvectors.rds")
module_membership <- readRDS("../results/wgcna/module_membership.rds")
module_name_tbl <- fread("../results/wgcna/wgcna_module_names.csv") %>%
  mutate(
    module_longname=str_replace(module_longname, "Family_XIII", "Clostridiales spp.") %>%
      str_replace("Prevotella_\\d+", "Prevotella")
  )

# load taxonomy of OTUs
load("../data/Bus_CF_working_noduplicates.Rdata")
taxonomy_df <- phyloseq::tax_table(Bus_CF1)@.Data %>%
  as.data.frame() %>%
  dplyr::rename(OTU=OTUID)
rm(Bus_CF1)

# load isolate -> OTU mapping
# removing anything with <95% identity (marked as NA)
isolate_otu_map <- fread("../data/isolate_names/otu_isolate_map.csv") %>%
  filter(!is.na(percentage_identity))

emapper_isolates_no_mapping <- setdiff(isolate_clusters$isolate, unique(isolate_otu_map$isolate))
mapping_isolates_no_emapper <- setdiff(unique(isolate_otu_map$isolate), isolate_clusters$isolate)
cat(sprintf("%d emapper isolates are missing from the mapping:\n%s\n",
            length(emapper_isolates_no_mapping), paste0(emapper_isolates_no_mapping, collapse="\n")))
cat(sprintf("%d mapping isolates are missing from emapper:\n%s\n",
            length(mapping_isolates_no_emapper), paste0(mapping_isolates_no_emapper, collapse="\n")))

isolate_otu_map_slim <- isolate_otu_map %>%
  group_by(isolate) %>%
  add_tally(name="n_OTUs") %>%
  ungroup() 

isolate_map_collapsed <- isolate_otu_map_slim %>%
  group_by(isolate) %>%
  summarise(n_OTUs=length(unique(OTU)), OTU_list=paste0(unique(sort(OTU)), collapse=","))

# OTU-isolate adjacency matrix
isolate_otu_adj_mat <- isolate_otu_map_slim %>%
  select(-c(n_OTUs, percentage_identity)) %>%
  mutate(value=1) %>%
  pivot_wider(names_from=isolate, values_from=value, values_fill=0) %>%
  column_to_rownames("OTU")

# extract non-redundant set of isolates
duplicate_isolate_info <- group_duplicate_columns(isolate_otu_adj_mat, 1)

# duplicate_isolate_info$tbl %>%
#   left_join(isolate_map_collapsed, by=c(included_column="isolate")) %>%
#   fwrite("../csv/duplicate_isolate_mappings_collapsed.tsv", sep="\t")

# check that if two isolates are mapped to the same OTUs they are
# also in the same cluster
same_cluster_check <- duplicate_isolate_info$tbl %>%
  left_join(isolate_clusters, by=c(included_column="isolate")) %>%
  left_join(isolate_clusters, by=c(excluded_column="isolate")) %>%
  filter(cluster.x!=cluster.y)
stopifnot(nrow(same_cluster_check)==0)

# isolate -> OTU map having removed any isolates that map to identical OTUs
# as another isolate
unique_isolate_otu_map <- isolate_otu_map_slim %>%
  filter(isolate %in% duplicate_isolate_info$non_duplicated_columns)
stopifnot(
  length(unique(unique_isolate_otu_map$OTU))==length(unique(isolate_otu_map_slim$OTU))
) # should not have removed any OTUs in this step

# score = fraction of OTUs for a given isolate in that module
isolate_module_scores <- lapply(
  module_assignments,
  function(x) x %>%
    left_join(unique_isolate_otu_map, by="OTU") %>%
    group_by(isolate, n_OTUs, module) %>%
    tally(name="n_OTUs_in_module") %>%
    mutate(value=n_OTUs_in_module/n_OTUs)
)
saveRDS(isolate_module_scores, file="../results/beta_regression/isolate_module_scores.rds")

#
# heatmap of isolate-module scores
isolate_plot_renamer <- duplicate_isolate_info$tbl %>%
  group_by(included_column) %>%
  summarise(label_suffix=sprintf("and %d others", length(unique(excluded_column)))) %>%
  mutate(label=sprintf("%s %s", included_column, label_suffix)) %>%
  rename(isolate=included_column, isolate_longname=label) %>%
  select(-label_suffix)

isolate_module_score_plotdata <- isolate_module_scores %>%
  rbindlist(idcol="site") %>%
  # mutate(isolate=fct_relevel("no_isolate", after=Inf)) %>%
  left_join(
    isolate_clusters, by="isolate"
  ) %>%
  left_join(module_name_tbl, by=c("site", "module")) %>%
  left_join(isolate_plot_renamer, by="isolate") %>%
  mutate(isolate_longname=coalesce(isolate_longname, isolate))

isolate_module_score_plotdata %>%
  filter(grepl("ots", site)) %>%
  ggplot(aes(y=isolate_longname, x=module_longname, fill=value)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="red", limits=c(0,1),
                      guide=guide_colorbar(frame.colour="black", ticks.colour="black")) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  xlab("WGCNA module") +
  ylab("Isolate") +
  facet_grid(rows=vars(cluster), cols=vars(site),
             scales="free", space="free", switch="y") +
  theme(strip.placement="outside",
        strip.background.y=element_blank(),
        strip.text.y=element_text(size=10),
        axis.text.y=element_text(size=8)) +
  labs(fill="Isolate-module\nscore")

############################################################################################################################
############################################################################################################################

############################################################################################################################
# relate KOs to WGCNA modules via isolates - beta regression
############################################################################################################################

library(betareg)

transform_score <- function(scores, s=0.5) {
  # https://www.researchgate.net/publication/7184584_A_better_lemon_squeezer_Maximum-likelihood_regression_with_beta-distributed_dependent_variables
  n <- length(scores)
  return((scores*(n-1) + s)/n)
}

# add duplicated isolates back in
all_isolate_merger <- duplicate_isolate_info$tbl %>%
  bind_rows(tibble(
    included_column=duplicate_isolate_info$non_duplicated_columns,
    excluded_column=duplicate_isolate_info$non_duplicated_columns)
  ) %>%
  filter(included_column!="no_isolate") %>%
  rename(isolate=excluded_column, join_key=included_column)

#
# fits univariate beta regression models
fit_beta_model <- function(formula, data, ...) {
  out <- tryCatch({
    return(betareg(formula=formula, data=data, ...))
  }, error=function(e) {
    return(NULL)
  })
  return(out)
}

logit <- function(x) { log(x/(1-x)) }
logistic <- function(x) { 1 / (1 + exp(-x)) }


# covariates are fixed for each site
beta_reg_covariates_df <- emapper %>%
  rownames_to_column("isolate")

# fit univariate models for each site
all_beta_model_coefs <- list()
all_model_fit_info <- list()
n_cores <- 20

# this loop takes about 5 minutes to run on 20 cores
for(site_name in c("cf_ots", "bus_ots")) {
  cat(site_name, "\n")
  
  # format the isolate-module scores for this site
  # TODO: add any unassigned isolates back in with score=0
  beta_reg_label_df <- isolate_module_scores[[ site_name ]] %>%
    ungroup() %>%
    select(isolate, module, value)  %>%
    filter(isolate!="no_isolate") %>%
    rename(join_key=isolate) %>%
    left_join(all_isolate_merger, by="join_key") %>%
    select(-join_key) %>%
    bind_rows(
      tibble(
        isolate=beta_reg_covariates_df$isolate[ !beta_reg_covariates_df$isolate %in% .$isolate ]
      )
    ) %>% 
    pivot_wider(names_from=module, values_from=value, values_fill=0) %>%
    select(-"NA")
  
  # fit univariate beta regressions for each module
  all_module_names <- colnames(beta_reg_label_df)[ colnames(beta_reg_label_df)!="isolate" ]
  for(module_name in all_module_names) {
    cat(module_name, "\n")
    
    # contains this label (module) and all covariates (KOs)
    df <- beta_reg_label_df %>%
      select(isolate, !!sym(module_name)) %>%
      mutate(target=transform_score(!!sym(module_name))) %>%
      select(-sym(module_name)) %>%
      left_join(beta_reg_covariates_df, by="isolate") %>%
      column_to_rownames("isolate") %>%
      drop_na()
    
    cat(sprintf("Fitting model on %d isolates\n", nrow(df)))
    ko_names <- colnames(df)[ colnames(df)!="target" ]
    cat(sprintf("Fitting model for %d KOs\n", length(ko_names)))
    
    # fit the models
    beta_reg_fits <- pbmcapply::pbmclapply(
      setNames(ko_names, ko_names),
      function(ko) fit_beta_model(
        formula=as.formula(sprintf("target~%s", ko)),
        data=df %>% select(ko, target)
      ),
      mc.cores=n_cores
    )
    
    # store formatted results
    idx <- paste0(site_name, "____", module_name)
    all_beta_model_coefs[[ idx ]] <- pbmcapply::pbmclapply(
      beta_reg_fits, function(x) tidy(x, conf.int=TRUE, conf.level=0.95), mc.cores=n_cores) %>%
      rbindlist(idcol="KO")
    all_model_fit_info[[ idx ]] <- pbmcapply::pbmclapply(
      beta_reg_fits, glance, mc.cores=n_cores) %>%
      rbindlist(idcol="KO")
  }
}

# format and save results
beta_coefs <- all_beta_model_coefs %>%
  rbindlist(idcol="label") %>%
  separate(label, c("site", "module"), sep="____") %>%
  left_join(duplicated_KOs %>%
              rename(KO=representative_ko, excluded_KO=identical_to) %>%
              bind_rows(
                tibble(KO=unique(.$KO), excluded_KO=unique(.$KO))
              ),
            by="KO") %>%
  rename(base_KO=KO, KO=excluded_KO) %>%
  relocate(KO, base_KO) %>%
  mutate(KO=coalesce(KO, base_KO)) %>%
  arrange(base_KO, KO)
beta_coefs %>%
  fwrite("../results/beta_regression/beta_reg_coefs.csv")

beta_reg_GoFs <- all_model_fit_info %>%
  rbindlist(idcol="label") %>%
  separate(label, c("site", "module"), sep="____") %>%
  left_join(duplicated_KOs %>%
              rename(KO=representative_ko, excluded_KO=identical_to) %>%
              bind_rows(
                tibble(KO=unique(.$KO), excluded_KO=unique(.$KO))
              ),
            by="KO") %>%
  rename(base_KO=KO, KO=excluded_KO) %>%
  relocate(KO, base_KO) %>%
  mutate(KO=coalesce(KO, base_KO)) %>%
  arrange(base_KO, KO)
beta_reg_GoFs %>%
  fwrite("../results/beta_regression/beta_reg_GoFs.csv")

##############################################################################
##############################################################################


##############################################################################
# evaluation of model fits
##############################################################################

# model goodness of fit
beta_reg_GoFs %>%
  arrange(desc(pseudo.r.squared))

# model coefficients
all_beta_model_coefs_slim <- beta_coefs %>%
  mutate(q.value=p.adjust(p.value, method="fdr")) %>%
  filter(!term %in% c("(Intercept)", "(phi)")) %>%
  arrange(q.value, estimate) %>%
  mutate(beta_mean=logistic(estimate)) %>%
  select(site, KO, base_KO, term, module, estimate, conf.low, conf.high, std.error, q.value)
  # mutate(est_upper_ci=estimate+1.96*std.error,
  #        est_lower_ci=estimate-1.96*std.error,
  #        odds_ratio=exp(estimate),
  #        odds_ratio_upper_ci=exp(est_upper_ci),
  #        odds_ratio_lower_ci=exp(est_lower_ci)
  # )

# stronger KO-module associations when the modules are more
# phylogenetically homogeneous
all_beta_model_coefs_slim %>%
  filter(q.value<1e-20) %>%
  ggplot(aes(x=term, y=module, fill=estimate)) +
  geom_tile() +
  facet_col(~site, scales="free_y", space="free") +
  theme(axis.text.x=element_text(size=6, angle=90, hjust=1)) +
  scale_fill_gradient2(
    low="blue", mid="white", high="red", na.value="white",
    limits=c(-5,5),
    guide=guide_colorbar(frame.colour="black", ticks.colour="black")
  )

rank_by_beta_mean <- function(beta_means) {
  return(dense_rank(desc((beta_means-0.5)**2.0)))
}

# isolate-module scores with all the duplicates
isolate_module_scores_all_isolates <- lapply(
  isolate_module_scores,
  function(x) x %>%
    ungroup() %>%
    select(isolate, module, value)  %>%
    filter(isolate!="no_isolate") %>%
    rename(join_key=isolate) %>%
    left_join(all_isolate_merger, by="join_key") %>%
    select(-join_key) %>%
    bind_rows(
      tibble(
        isolate=beta_reg_covariates_df$isolate[ !beta_reg_covariates_df$isolate %in% .$isolate ]
      )
    ) %>% 
    pivot_wider(names_from=module, values_from=value, values_fill=0) %>%
    select(-"NA")
)

# sanity checks for top hits
top_hits <- all_beta_model_coefs_slim %>%
  filter(site=="bus_ots") %>%
  # mutate(ko_rank=rank_by_beta_mean(beta_mean)) %>%
  arrange(q.value, desc(abs(estimate))) %>%
  head(20) %>%
  mutate(label=sprintf("%s (%.2e,%.2e,%s)", term, estimate, q.value, module))

emapper[,top_hits$term] %>%
  rownames_to_column("isolate") %>%
  pivot_longer(-isolate, names_to="KO_name", values_to="KO_value") %>%
  left_join(
    isolate_module_scores_all_isolates$bus_ots %>%
      pivot_longer(-isolate, names_to="module")
  ) %>%
  rename(isolate_module_score=value) %>%
  mutate(KO_value=as.numeric(KO_value)) %>%
  pivot_wider(names_from=module, values_from=isolate_module_score) %>%
  pivot_wider(names_from=KO_name, values_from=KO_value) %>%
  pivot_longer(-c(isolate)) %>%
  mutate(facet_title=ifelse(grepl("K\\d+", name), "KO", "im-score")) %>%
  left_join(
    top_hits %>% select(term, label) %>% rename(name=term)
  ) %>%
  mutate(label=coalesce(label, name)) %>%
  left_join(isolate_clusters) %>%
  arrange(cluster) %>%
  mutate(isolate=factor(isolate, levels=unique(isolate))) %>%
  ggplot(aes(x=isolate, y=label, fill=value)) +
  geom_tile() +
  facet_grid(rows=vars(facet_title), cols=vars(cluster), 
             scales="free", space="free") +
  theme(axis.text.x=element_text(size=8, angle=90, hjust=1)) +
  scale_fill_gradient(low="white", high="red") 

##############################################################################
##############################################################################

##############################################################################
# search for KOs that relate to modules
##############################################################################

# attach KO information to beta coefficients
top_hits <- all_beta_model_coefs_slim %>%
  # group_by(site, module) %>%
  # slice_max(estimate, n=20) %>%
  # ungroup() %>%
  filter(conf.low>1.0) %>%
  rename(wgcna_module=module)
cat(sprintf("Top KO hits include %d unique KOs\n", length(unique(top_hits$KO))))

top_hit_ko_lookup_raw <- lookup_kos(unique(top_hits$KO), 20)
parsed_top_hit_info <- pbmclapply(top_hit_ko_lookup_raw, parse_brite, mc.cores=20)

plot_data <- parsed_top_hit_info %>% 
  bind_rows() %>%
  distinct() %>%
  inner_join(top_hits, by="KO") %>%
  group_by(site, wgcna_module, term1, term2, term3) %>%
  tally() %>%
  ungroup() %>%
  arrange(term1, term2) %>%
  mutate(
    term1=str_remove(term1, "\\d+") %>% trimws(),
    term2=factor(term2, levels=unique(term2)),
    term3=factor(term3, levels=unique(term3))
  ) %>%
  left_join(
    module_name_tbl %>%
      rename(wgcna_module=module, wgcna_module_displayname=module_longname),
    by=c("site", "wgcna_module")
  )

plot_data2 <- plot_data %>%
  filter(wgcna_module!="grey", site=="bus_ots") %>%
  mutate(
    wgcna_module_displayname=str_replace(
      wgcna_module_displayname, " \\(", "\n"
    ) %>%
      str_remove("\\)")
  ) %>%
  mutate(term1=str_remove(term1, "\\d+$"))

all_term1_vals <- unique(plot_data2$term1)

plt <- ggplot() +
  geom_col(
    data=plot_data2 %>% filter(term1==all_term1_vals[[1]]),
    aes(y=term3, x=n, fill=term2))  +
  labs(fill=all_term1_vals[[1]])

for(term1_val in tail(unique(plot_data2$term1), -1)) {
  cat(term1_val, "\n")
  plt <- plt +
    ggnewscale::new_scale_fill() +
    geom_col(
      data=plot_data2 %>% filter(term1==term1_val),
      aes(y=term3, x=n, fill=term2)) +
    labs(fill=term1_val)
    
}
plt <- plt + 
  facet_wrap(~wgcna_module_displayname) +
  theme(
    legend.title=element_text(size=10),
    legend.text=element_text(size=8)
  )

plot_data2$term1 %>% table()

plot_data2 %>%
  # filter(!term1 %in% c("Not Included in Pathway or Brite", "Brite Hierarchies")) %>%
  filter(n>1) %>%
  ggplot(aes(y=term2, x=n, fill=term1)) +
  geom_col() +
  theme(axis.text.y=element_text(size=10),
        legend.text=element_text(size=8),
        legend.position="none",
        legend.title=element_blank(),
        strip.placement="outside",
        strip.text.y.left=element_text(size=10, angle=0, hjust=1),
        strip.text.x=element_text(size=10)) +
  xlab("Number of KOs in pathway") +
  ylab("KO pathway") +
  facet_grid(rows=vars(term1), cols=vars(wgcna_module_displayname),
             scales="free_y", space="free",
             switch="y")

ggsave(
  "../plots/beta_reg/busselton_ots_wgcna_KOs.pdf",
  height=12, width=14
)

##########################################################################
# comparing the two streptococci clusters
##########################################################################

get_point_colour <- function(q1, q2, alpha) {
  stopifnot(length(q1)==length(q2))
  out <- rep("neither", length(q1))
  out[ q1<=alpha & q2<=alpha ] <- "both"
  out[ q1>=alpha & q2<=alpha ] <- "just_2"
  out[ q1<=alpha & q2>=alpha ] <- "just_1"
  return(out)
}

strep_scatter_data <- all_beta_model_coefs_slim %>%
  filter(site=="bus_ots", module %in% c("blue", "yellow")) %>%
  select(-c(std.error)) %>%
  pivot_wider(
    names_from=module,
    values_from=c(estimate, conf.low, conf.high, q.value)
    ) %>%
  mutate(
    point_colour=get_point_colour(q.value_blue, q.value_yellow, 0.1),
    point_alpha=if_else(point_colour %in% c("neither", "both"), "low", "high")
  )

strep_scatter_data_plot <- strep_scatter_data %>%
  mutate(sign_yellow=estimate_yellow>0, sign_blue=estimate_blue>0) %>%
  filter(
    !data.table::between(0, conf.low_blue, conf.high_blue) |
    !data.table::between(0, conf.low_yellow, conf.high_yellow),
  ) %>%
  filter(sign_blue!=sign_yellow)

# lookup KO pathways
strep_ko_lookup_raw <- lookup_kos(unique(strep_scatter_data_plot$KO), 20)
parsed_strep_ko_info <- pbmclapply(strep_ko_lookup_raw, parse_brite, mc.cores=20)

strep_scatter_data_plot %>%
  left_join(
    parsed_strep_ko_info %>% bind_rows(),
    by="KO"
  ) %>%
  mutate(term1=str_remove(term1, "\\d+$")) %>%
  ggplot(aes(x=estimate_blue, y=estimate_yellow, colour=term1)) +
  geom_point(size=1.5) +
  # geom_errorbarh(aes(xmin=conf.low_blue, xmax=conf.high_blue)) +
  # geom_errorbar(aes(ymin=conf.low_yellow, ymax=conf.high_yellow)) +
  geom_hline(yintercept=0, colour="red") +
  geom_vline(xintercept=0, colour="red") +
  geom_abline(colour="red")

strep_scatter_data %>%
  ggplot(aes(x=estimate_blue, y=estimate_yellow)) +
  geom_point(aes(alpha=point_alpha), size=0.5) +
  geom_hline(yintercept=0, colour="red") +
  geom_vline(xintercept=0, colour="red") +
  geom_abline() +
  scale_colour_manual(
    values=c(both="green", just_1="blue", just_2="orange", neither="grey50")
  ) +
  scale_alpha_manual(values=c(1.0, 0.2)) +
  guides(alpha="none")

##########################################################################
##########################################################################

# table for the networks to name them
# # 
