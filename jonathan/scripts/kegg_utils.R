get_duplicates <- function(x, data_matrix) {
  # returns column names in data_matrix that are identical to x
  identical_to_this <- sapply(
    colnames(data_matrix),
    function(other_ko_name) identical(
      x, data_matrix[ ,other_ko_name ]
    )
  )
  identical_to_this <- colnames(data_matrix)[ identical_to_this ]
  return(identical_to_this)
}

group_duplicate_columns <- function(df, n_cores=1) {
  cat(sprintf("Finding duplicates for dataframe with %d columns using %d cores\n", ncol(df), n_cores))
  
  # remove duplicated KOs
  duplicated_columns <- colnames(df)[ t(df) %>% duplicated() ]
  cat(sprintf("%d of %d columns are duplicates\n", length(duplicated_columns), ncol(df)))
  
  non_duplicated_columns <- colnames(df)[ !t(df) %>% duplicated() ]
  cat(sprintf("%d of %d KOs are NOT duplicates\n", length(non_duplicated_columns), ncol(df)))
  
  if(length(duplicated_columns)==0) {
    return()
  }
  
  if(length(duplicated_columns)+length(non_duplicated_columns)!=ncol(df)) {
    warning(
      sprintf(
        "# duplicates (%d) + # non-duplicates (%d) not equal to total number of columns (%d) ",
        length(duplicated_columns),
        length(non_duplicated_columns),
        ncol(df)
      ),
      call.=FALSE
    )
  }
  
  # search non-duplicated columns for the corresponding duplicate
  duplicate_column_log <- pbmcapply::pbmclapply(
    setNames(duplicated_columns, duplicated_columns),
    function(x) get_duplicates(df[[ x ]], df[,non_duplicated_columns]),
    mc.cores=n_cores
  )
  
  tbl <- tibble(
    included_column=simplify2array(duplicate_column_log),
    excluded_column=names(duplicate_column_log)
  )
  
  return(
    list(
      new_df=df[,non_duplicated_columns],
      non_duplicated_columns=non_duplicated_columns,
      duplicated_columns=duplicated_columns,
      tbl=tbl
    )
  )
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

calculate_or <- function(x, group_vec) {
  # calculate odds ratios from 2x2 contingency table
  # adding 0.5 if there are zero cells
  cont_tab <- table(x, group_vec) %>% as.matrix()
  
  zero_cells <- any(cont_tab==0)
  
  if(zero_cells) {
    # Haldane-Anscombe correction
    cont_tab <- cont_tab + 0.5
  }
  
  or <- cont_tab %>% effectsize::oddsratio() %>% as.data.frame()
  return(data.frame(or, zero_cells=zero_cells))
}

collapse_to_isolates <- function(otu_reads, isolate_otu_map_df) {
  
  # from wide OTU reads to wide isolate reads using the mapping isolate_otu_map_df
  
  otu_reads[ ,intersect(colnames(otu_reads), unique(isolate_otu_map_df$OTU)) ] %>%
    as.data.frame() %>%
    rownames_to_column("sample_id") %>%
    pivot_longer(-sample_id, names_to="OTU", values_to="reads")  %>%
    left_join(isolate_otu_map_df) %>% 
    group_by(label, sample_id) %>%
    summarise(total_reads=sum(reads)) %>%
    pivot_wider(names_from=label, values_from=total_reads) %>%
    column_to_rownames("sample_id")
}

# tidy_mantel <- function(dx, dy, ...) {
#   out  <- vegan::mantel(dx, dy, ...)
#   return(tibble(estimate=out$statistic, pvalue=out$signif))
# }
# 
# test_ko_abundance_sim <- function(x_counts, x_ko, n_boots, n_cores, ...) {
#   # mantel test between correlation distances between isolate counts
#   # and KO distances between isolates
#   
#   common_isolates <- intersect(colnames(x_counts), rownames(x_ko))
#   cat(sprintf("%d common isolates\n", length(common_isolates)))
#   
#   ko_dist <- x_ko[ common_isolates, ] %>% dist("manhattan")
#   corr_dist <- x_counts[ ,common_isolates] %>% get_cor_dist()
#   
#   obs_mantel <- tidy_mantel(ko_dist, corr_dist, ...)
#   
#   # bootstrap p-value and mantel estimate
#   resampled_rows <- lapply(
#     1:n_boots,
#     function(...) sample(1:nrow(x_counts), size=nrow(x_counts), replace=TRUE)
#   )
#   
#   bootstrap_replicates <- pbmcapply::pbmclapply(
#     resampled_rows,
#     function(idx) x_counts[ idx, ] %>%
#       get_cor_dist() %>%
#       tidy_mantel(ko_dist),
#     mc.cores=n_cores
#   )
#   
#   # format results
#   out <- bootstrap_replicates %>%
#     rbindlist(idcol="resample")
#   out2 <- tibble(ci_upper=quantile(out$estimate, 0.975),
#                  ci_lower=quantile(out$estimate, 0.025),
#                  name="estimate")
#   out3 <- tibble(ci_upper=quantile(out$pvalue, 0.975),
#                  ci_lower=quantile(out$pvalue, 0.025),
#                  name="pvalue")
#   obs_mantel %>%
#     pivot_longer(everything(), values_to="observed") %>%
#     left_join(
#       bind_rows(out2, out3), by="name"
#     )
# }

recode_legend_keys <- function(x_vals, facet_vals, eps) {
  
  tbl <- tibble(x=x_vals, facet=facet_vals) %>%
    distinct() %>%
    replace_na(list(facet="<NA>")) %>%
    group_by(facet) %>%
    add_tally(name="facet_count")
  
  tbl <- tbl %>%
    mutate(facet_label=ifelse(facet_count>eps, facet, "Other"))
  
  mapping <- setNames(tbl$facet_label, tbl$x)
  return(mapping[ x_vals ])
}

# recode_legend_keys <- function(vec, eps) {
#   vec_counts <- table(vec, useNA="always") %>%
#     as.data.frame(stringsAsFactors=FALSE) %>%
#     mutate(label=ifelse(Freq > eps, vec, "Other"))
#   label_counts <- vec_counts %>%
#     group_by(label) %>%
#     summarise(n=sum(Freq))
#   label_counts$label_new <- paste0(label_counts$label, " (", label_counts$n, ")")
#   vec_counts <- vec_counts %>%
#     left_join(label_counts)
#   mapping <- vec_counts$label_new
#   names(mapping) <- vec_counts$vec
#   vec_counts <- vec_counts %>% arrange(desc(n))
#   mapping <- factor(mapping, levels=unique(vec_counts$label_new))
#   return(sapply(vec, function(x) mapping[[x]]))
# }