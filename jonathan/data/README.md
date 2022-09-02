```├── data
│   ├── Bus_CF_working_noduplicates.Rdata
│   ├── isolate_names
│   │   ├── filtered_OTU_isolate.txt
│   │   ├── new_label_isolates_nsp.csv
│   │   └── otu_isolate_map.csv
│   ├── kegg
│   │   ├── formatted
│   │   │   ├── binary_emapper_matrix_no_duplicates.csv
│   │   │   ├── duplicated_KOs.csv
│   │   │   ├── isolate_ko_matrix.csv
│   │   │   └── parsed_pathway_info.csv
│   │   ├── formatted_ko_lookup.tsv
│   │   ├── kegg_pathway_categories.csv
│   │   ├── ko_annotations
│   │   │   ├── emapper_outputs
│   │   │   │   ├── new_haemophilus_isolates
│   │   │   └── reference_files
│   │   │       ├── emapper_all_binary.tsv
│   │   │       └── emapper_all.tsv
│   │   └── ko_lookup_raw.rds
│   ├── metadata
│   │   ├── Combined Celtic Fire and Busselton16S samples metadata for microbiome analyses correct.xlsx
│   │   └── formatted_metadata.csv
│   └── wgcna
│       └── otu_tables_wgcna.rds
├── markdown
│   ├── metadata_formatting.Rmd
│   └── wgcna_analysis.Rmd
├── results	
│   ├── kegg
│   │   ├── cluster_ko_odds_ratios.rds
│   │   ├── final_KO_cluster_scores.tsv
│   │   ├── isolate_cluster_labels.csv
│   │   ├── isolate_cluster_labels.rds
│   │   ├── isolate_phylogenetic_hclust.rds
│   │   └── supplementary_table_2.xlsx
│   └── wgcna
│       ├── hub_otu_labels.rds
│       ├── module_assignments.rds
│       ├── module_eigenvectors.rds
│       ├── module_membership.rds
│       ├── module_trait_assocs.rds
│       ├── otu_trees.rds
│       ├── tom_matrices.rds
│       └── wgcna_module_names.csv
└── scripts
              ├── asthma_vs_nonasthma_correlation.R
              ├── eggnog_parsing.R
              ├── isolate_otu_mapping.R
              ├── kegg_analysis.R
              ├── kegg_utils.R
              ├── ko_lookup_utils.R
              ├── paper_correlation_heatmap.R
              ├── paper_isolate_otu_mapping.R
              ├── parse_kegg_info.R
              ├── plot_utils.R
              ├── supplementary_figures.R
              ├── wgcna_module_composition.R
              └── wgcna_module_isolate_cluster_betareg.R
```
