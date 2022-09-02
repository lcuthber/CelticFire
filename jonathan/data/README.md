```
├── data
│   ├── Bus_CF_working_noduplicates.Rdata                       # decontaminated and filtered and phyloseq
│   ├── isolate_names
│   │   ├── filtered_OTU_isolate.txt                            # OTU-isolate alignment
│   │   ├── new_label_isolates_nsp.csv                          # new isolate names
│   │   └── otu_isolate_map.csv                                 # OTU-isolate sequence identity
│   ├── kegg
│   │   ├── formatted
│   │   │   ├── binary_emapper_matrix_no_duplicates.csv         # KO presence/absence in each isolate (126 x 5,531 binary matrix)
│   │   │   ├── duplicated_KOs.csv                              # list of redundant KOs
│   │   │   ├── isolate_ko_matrix.csv                           # KO prevalence in each isolate (126 x 5,531 count matrix)
│   │   │   └── parsed_pathway_info.csv                         # KO pathway categories
│   │   ├── formatted_ko_lookup.tsv                             # KO descriptions
│   │   ├── kegg_pathway_categories.csv                         # KO pathways listed according to category
│   │   ├── ko_annotations
│   │   │   ├── emapper_outputs                                 # eggnog-mapper output for 126 isolates
│   │   │   │   ├── new_haemophilus_isolates                    # eggnog-mapper output for 8 additional Haemophilus isolates (to remove)
│   │   │   └── reference_files                                 # previous versions of files for unit testing
│   │   │       ├── emapper_all_binary.tsv
│   │   │       └── emapper_all.tsv
│   │   └── ko_lookup_raw.rds                                   # KO database query for all 5,531 KOs
│   ├── metadata
│   │   ├── Combined Celtic Fire and[...].xlsx                  # Metadata prepared by Miriam
│   │   └── formatted_metadata.csv                              # metadata formatted for WGNCA analysis
│   └── wgcna
│       └── otu_tables_wgcna.rds                                # 4 OTU tables (one per site) used in WGCNA analysis
├── markdown
│   ├── metadata_formatting.Rmd                                 # formats Miriam's metadata for WGCNA
│   └── wgcna_analysis.Rmd                                      # runs WGCNA analysis
├── results	
│   ├── kegg
│   │   ├── cluster_ko_odds_ratios.rds                          # odds ratio scores for each KO-isolate cluster pair (raw)
│   │   ├── final_KO_cluster_scores.tsv                         # odds ratio scores for each KO-isolate cluster pair (formatted)
│   │   ├── isolate_cluster_labels.csv                          # cluster assigments for isolates
│   │   ├── isolate_cluster_labels.rds                          # cluster assigments for isolates
│   │   ├── isolate_phylogenetic_hclust.rds                     # dendrogram for isolate clustering using KO distance matrix
│   │   └── supplementary_table_2.xlsx                          # odds ratio scores for each KO-isolate cluster pair (manuscript version)
│   └── wgcna
│       ├── hub_otu_labels.rds                                  # hub OTUs for each WGCNA module
│       ├── module_assignments.rds                              # WGCNA module assignments
│       ├── module_eigenvectors.rds                             # WGCNA module eigenvectors
│       ├── module_membership.rds                               # WGCNA module membership
│       ├── module_trait_assocs.rds                             # associations between host traits and WGCNA modules
│       ├── otu_trees.rds                                       # WGNCA dendrograms
│       ├── tom_matrices.rds                                    # WGCNA TOM matrices
│       └── wgcna_module_names.csv                              # WGCNA module names
└── scripts
              ├── eggnog_parsing.R                              # parses lookups from KO database
              ├── isolate_otu_mapping.R                         # prepare the isolate-OTU map
              ├── kegg_analysis.R                               # run KO-based isolate clustering
              ├── kegg_utils.R                                  # utility functions for isolate clustering
              ├── ko_lookup_utils.R                             # utility functions for KO database lookups
              ├── paper_correlation_heatmap.R                   # Figure 2
              ├── parse_kegg_info.R                             # parse raw KO lookups
              ├── plot_utils.R                                  # plot saving function
              ├── supplementary_figures.R                       # Figure S3a, S3b
              ├── wgcna_module_composition.R                    # makes tables of WGNCA module composition - Tables not currently in manuscript
              └── wgcna_module_isolate_cluster_betareg.R        # not currently in manuscript - to remove
```
