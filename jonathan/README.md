The `data` subdirectory contains raw and processed files required to run the code in the `scripts` directory. These scripts assume that the current working directory is `run`.

The workflow is divided into three data-preprcessing steps (you do not have to run these as their results are already saved in the repo):

1. Parse the raw eggnog-mapper outputs into a matrix (run `scripts/eggnog_parsing.R` to produce `data/kegg/formatted/isolate_ko_matrix.csv`).
2. Lookup the KO information (run `scripts/parse_kegg_info.R` to produce `data/kegg/formatted/parsed_pathway_info.csv`).
3. Format the sample metadata (run `markdown/metadata_formatting.Rmd` to produce `results/metadata/formatted_metadata.csv`).
4. Create a formatted OTU -> isolate mapping file (run `scripts/isolate_otu_mapping.R`).

This produces the files needed to run the following analyses:

1. WGCNA analysis of the OTU abundances at the four sites (`markdown/wgcna_analysis.Rmd`, results in `results/wgcna`). 
2. Hierarchichal clustering of the isolates based on their KOs and scoring each cluster-KO pair using odds ratios (`scripts/kegg_analysis.R`, results in `results/kegg`).
3. Associate the WGCNA network modules (defined on the OTUs) with the isolate clusters using beta regression (`scripts/wgcna_module_isolate_cluster_betareg.R`).