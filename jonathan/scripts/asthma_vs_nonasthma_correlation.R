library(phyloseq)
library(data.table)
library(tidyr)
library(dplyr)
library(stringr)
library(tibble)
library(corrr)

library(ggplot2)
theme_set(theme_linedraw(base_size=14))

source("../scripts/kegg_utils.R")
source("../scripts/plot_utils.R")

########################################################################
# load required objects
########################################################################

# OTU counts
load("../data/Bus_CF_working_noduplicates.Rdata")

X <- list()
X[["cf_lll"]] <- t(as.matrix(otu_table(subset_samples(Bus_CF1, SwabSite=="Brushing LLL 16S"))))
X[["cf_lul"]] <- t(as.matrix(otu_table(subset_samples(Bus_CF1, SwabSite=="Brushing LUL 16S"))))
X[["cf_ots"]] <- t(as.matrix(otu_table(subset_samples(Bus_CF1, SwabSite=="OTS"))))
X[["bus_ots"]] <- t(as.matrix(otu_table(subset_samples(Bus_CF1, SwabSite=="Throat"))))

taxonomy_df <- phyloseq::tax_table(Bus_CF1)@.Data %>%
  as.data.frame() %>%
  rename(OTU=OTUID)

sample_metadata <- fread("../data/metadata/formatted_metadata.csv") %>%
  as_tibble()

# to subsample the OTUs to those used in WGCNA
otu_trees <- readRDS("../results/wgcna/otu_trees.rds")
wgcna_modules <- readRDS("../results/wgcna/module_assignments.rds")

########################################################################
########################################################################

#
# correlation structure in asthmatics vs non-asthmatics
for(SITE_NAME in names(X)) {
  cat(SITE_NAME, "\n")

  included_OTUs <- otu_trees[[ SITE_NAME ]] %>%
    as.dendrogram() %>%
    labels()
  x <- X[[ SITE_NAME ]][ , included_OTUs]
  
  grouped_sample_IDs <- sample_metadata %>%
    filter(SampleID %in% rownames(x)) %>%
    group_by(Asthma) %>%
    summarise(sample_ids=list(unique(SampleID)))
  
  grouped_sample_IDs$corrmat <- lapply(
    grouped_sample_IDs$sample_ids,
    function(ids) cor(x[ ids, ], method="spearman") %>%
      as_cordf(diagonal=1) %>%
      stretch()
  )
  
  otu_order <- wgcna_modules[[ SITE_NAME ]] %>%
    arrange(module)
  facet_order <- unique(otu_order$module)
  
  otu_summary <- tibble(
    OTU=colnames(x),
    pc_reads=100*colSums(x)/sum(x),
    prevalence=100*apply(x, 2, function(xx) sum(xx>0)/length(xx))
  ) %>%
    left_join(wgcna_modules[[ SITE_NAME ]], by="OTU")
  
  plotted_OTUs <- otu_summary %>% slice_max(pc_reads, n=250) %>% pull(OTU)
  
  fig <- grouped_sample_IDs %>%
    filter(Asthma %in% c(0,1)) %>%
    mutate(n_samples=sapply(sample_ids, length)) %>%
    select(-sample_ids) %>%
    mutate(Asthma=if_else(Asthma==1, "Asthma", "Control"),
           Asthma=sprintf("%s (n=%d)", Asthma, n_samples)) %>%
    unnest(corrmat) %>%
    filter(x %in% plotted_OTUs, y %in% plotted_OTUs) %>%
    mutate(x=factor(x, levels=otu_order$OTU),
           y=factor(y, levels=otu_order$OTU)) %>%
    left_join(
      otu_order, by=c(x="OTU")
    ) %>%
    left_join(
      otu_order, by=c(y="OTU")
    ) %>%
    mutate(
      module.x=factor(module.x, levels=facet_order),
      module.y=factor(module.y, levels=rev(facet_order)),
    ) %>%
    group_by(Asthma) %>%
    do(
      plot=ggplot(data=., aes(x=x, y=y, fill=r)) +
        geom_tile() +
        facet_grid(cols=vars(module.x), rows=vars(module.y), scales="free", space="free") +
        theme(axis.text=element_text(size=2),
              axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
              strip.background=element_blank(),
              strip.text=element_text(colour="black"),
              strip.text.x=element_text(angle=90),
              strip.text.y=element_text(angle=0),
              panel.spacing=unit(1, "pt")) +
        scale_fill_gradient2(low="blue", mid="white", high="red", limits=c(-1, 1)) +
        ggtitle(.$Asthma)
    ) %>%
    pull(plot) %>%
    ggpubr::ggarrange(plotlist=., nrow=1, common.legend=TRUE, legend="top")
    
  ggsave(sprintf("../plots/wgcna/spearcor_asthma_control__%s.pdf", SITE_NAME),
         plot=fig,
         width=16, height=9)
}

  


# #
# # differential correlation analysis
# library(DGCA)
# design_mat <- sample_metadata %>%
#   filter(SampleID %in% rownames(x)) %>%
#   arrange(match(SampleID, rownames(x))) %>%
#   select(SampleID, Asthma) %>% 
#   mutate(Asthma=if_else(Asthma==1, "Asthma", "Control")) %>%
#   mutate(value=1) %>%
#   pivot_wider(values_from=value, names_from=Asthma, values_fill=0) %>%
#   column_to_rownames("SampleID") %>%
#   as.matrix()
# 
# ddcor_res <- ddcorAll(
#   t(x),
#   design=design_mat,
#   compare=c("Asthma", "Control"),
#   corrType="pearson",
#   adjust="fdr",
#   nPerms=0
# )
# 
# plotted_OTUs <- otu_summary %>% slice_max(pc_reads, n=50) %>% pull(OTU)
# 
# top_hits <- ddcor_res %>%
#   as_tibble() %>%
#   filter(Gene1 %in% plotted_OTUs | Gene2 %in% plotted_OTUs) %>%
#   arrange(pValDiff_adj) %>%
#   head(10)
# 
# make_scatterplot <- function(x, xname, yname) {
#   ggplot(x, aes_string(x=name, y=yname))
# }
# 
# plotCors(ddcor_res, )
# 
# plotCors(inputMat=t(x) %>% log1p(), design=design_mat,
#          compare=c("Asthma", "Control"), geneA="Streptococcus_30181", geneB="Prevotella_7_9923")
#   scale_y_continuous(trans="log1p")
#   
# 
