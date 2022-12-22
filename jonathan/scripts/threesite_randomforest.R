library(tidymodels)
library(ranger)

#####################################################################################
# load data
#####################################################################################



#####################################################################################
#####################################################################################

#####################################################################################
# site vs site ecological analysis
#####################################################################################

n_cores <- 16
n_outer_folds <- 10
n_inner_folds <- 10

# returns an RF model
make_rf_mod <- partial(
  rand_forest,
  trees=500,
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

# check evaluations
all_validation_preds <- outer_fold_preds %>%
  lapply(function(x) bind_rows(x, .id="outer_fold")) %>%
  bind_rows(.id="site")

#
# mean AUCs (ROC and PRC)
all_validation_preds %>%
  group_by(site, outer_fold) %>%
  metric_set(roc_auc, pr_auc)(Disease, .pred_Control) %>%
  group_by(site, .metric) %>%
  summarise(value=mean(.estimate)) %>%
  pivot_wider(names_from=.metric)

all_validation_preds %>%
  group_by(site, outer_fold) %>%
  metric_set(roc_auc, pr_auc)(Disease, .pred_Control) %>%
  ggpaired(x="site", y=".estimate", facet.by=".metric") +
  get_plot_theme(14) +
  geom_hline(
    data=data.frame(.metric=c("pr_auc", "roc_auc"), y=c(0.29, 0.5)),
    aes(yintercept=y),
    colour="red", linetype="dashed", size=2
  )

#
# hard class label metrics
hard_metrics <- lapply(
  seq(0, 1, length.out=20),
  function(eps) all_validation_preds %>%
    group_by(site, outer_fold) %>%
    mutate(
      .pred=probably::make_two_class_pred(
        .pred_Control, levels=c("Control", "Asthma"), threshold=eps
      )
    ) %>%
    metric_set(accuracy, bal_accuracy, kap, mcc, f_meas, precision, recall, sens, spec, j_index)(Disease, estimate=.pred) %>%
    mutate(.threshold=eps) %>%
    group_by(site, .metric, .threshold) %>%
    summarise(.estimate=mean(.estimate, na.rm=TRUE)) %>%
    ungroup()
) %>%
  bind_rows() 
hard_metrics


#
# random baseline
random_hard_metrics <- pbmcapply::pbmclapply(
  1:100,
  function(...) all_validation_preds %>%
  group_by(site, outer_fold) %>%
  mutate(.pred=sample(Disease)) %>%
  metric_set(
    accuracy, bal_accuracy, kap, mcc, f_meas, precision, recall, sens, spec, j_index
    )(Disease, estimate=.pred),
  mc.cores=10
)

random_baseline_plotdata <- random_hard_metrics %>%
  bind_rows(.id="resample") %>%
  group_by(site, outer_fold, .metric) %>%
  summarise(.mean_estimate=mean(.estimate),
            .estimate_lower=quantile(.estimate, p=0.05),
            .estimate_upper=quantile(.estimate, p=0.95)) %>%
  group_by(site, .metric) %>%
  summarise(.estimate=mean(.mean_estimate)) %>%
  ungroup()
  

hard_metrics %>%
  ggplot(aes(x=.threshold, y=.estimate)) +
  geom_line(aes(group=site, colour=site)) +
  facet_wrap(~.metric, scales="free_y", nrow=2) +
  theme(legend.position="top") +
  get_plot_theme(14) +
  geom_hline(
    data=random_baseline_plotdata %>% filter(site=="cf_ots"),
    aes(yintercept=.estimate)
  )


#####################################################################################
#####################################################################################

#####################################################################################
# site vs site ecological analysis
#####################################################################################

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
cl <- makeForkCluster(n_cores)
registerDoParallel(cl)

outer_fold_preds_sitemodel <- list()
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
  outer_fold_preds_sitemodel[[ outer_fold_idx ]] <- make_rf_mod(
    mtry=best_hparams$mtry, min_n=best_hparams$min_n
  ) %>%
    parsnip::fit(site ~ ., data=analysis(cv_resamples$splits[[ outer_fold_idx ]])) %>%
    predict(assessment(cv_resamples$splits[[ outer_fold_idx ]]), type="prob") %>%
    bind_cols(assessment(cv_resamples$splits[[ outer_fold_idx ]]) %>% select(site)) 
}

stopCluster(cl)

all_validation_preds <- outer_fold_preds_sitemodel %>%
  bind_rows(.id="outer_fold")
hard_preds <- all_validation_preds %>%
  mutate(
    .pred=probably::make_class_pred(
      .pred_cf_lll, .pred_cf_lul, .pred_cf_ots,
      levels=c("cf_lll", "cf_lul", "cf_ots"),
    )
  )

# multiclass AUC (one-vs-all)
multiclass_aucs <- list()
for(class_name in levels(rf_df$site)) {
  multiclass_aucs[[ class_name ]] <- all_validation_preds %>%
    mutate(
      site=if_else(as.character(site)==class_name, as.character(site), paste0("not_", class_name)),
      site=factor(site, levels=c(class_name, paste0("not_", class_name)))
    ) %>% 
    group_by(outer_fold) %>%
    roc_auc(estimate=!!sym(paste0(".pred_", class_name)), truth=site) %>%
    group_by(.metric) %>%
    summarise(.estimate=mean(.estimate))
  
}
for(class_name in levels(rf_df$site)) {
  multiclass_aucs[[ class_name ]] <- multiclass_aucs[[ class_name ]] %>%
    bind_rows(
      all_validation_preds %>%
        mutate(
        site=if_else(as.character(site)==class_name, as.character(site), paste0("not_", class_name)),
        site=factor(site, levels=c(class_name, paste0("not_", class_name)))
      ) %>% 
      group_by(outer_fold) %>%
      pr_auc(estimate=!!sym(paste0(".pred_", class_name)), truth=site) %>%
      group_by(.metric) %>%
      summarise(.estimate=mean(.estimate))
    )
}

multiclass_aucs %>%
  bind_rows(.id="site") %>%
  pivot_wider(names_from=.metric, values_from=.estimate)


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

# confusion matrix
hard_preds %>%
  conf_mat(truth=site, estimate=.pred) %>%
  .$table %>%
  as.data.frame() %>%
  mutate(cell_type=if_else(Prediction==Truth, "Correct", "Incorrect")) %>%
  ggplot(aes(x=Truth, y=Prediction, fill=Freq)) +
  geom_tile() +
  geom_label(aes(label=Freq, colour=cell_type), size=10, show.legend=FALSE) +
  get_plot_theme(18) +
  scale_fill_gradient(low="white", high="grey") +
  scale_colour_manual(values=c(Correct="darkgreen", Incorrect="red")) +
  theme(aspect.ratio=1)
cropped_ggsave("../plots/cf_lll_vs_lul/confusion_matrix_threesite.pdf",
               height=7, width=7)

#####################################################################################
#####################################################################################
