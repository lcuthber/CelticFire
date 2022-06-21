#
# lookup KEGG information (e.g. pathways, definition, name, ...) from server
# Also parses the hierarchichal structure for the pathways into a dataframe

library(data.table)
library(dplyr)
library(tidyr)

source("../scripts/ko_lookup_utils.R")

####################################################################################
# Lookup the KO info
####################################################################################

# lookup pathways, modules, names, etc for all KOs (takes about 10 mins on 20 cores)
all_kos <- fread("../data/kegg/formatted/isolate_ko_matrix.csv") %>%
  colnames()
all_kos <- all_kos[ all_kos!="isolate" ]
ko_lookup_raw <- lookup_kos(all_kos, n_cores=20)

# format into single dataframe
cat("Loaded KO information for", length(ko_lookup_raw), "KOs\n")
ko_lookup_slim <- lapply(ko_lookup_raw, format_query) %>%
  rbindlist() %>%
  rename(KO=entry) %>%
  get_slim_query()

# save formatted KO lookup. Use tsv because some columns contain comma-separeted
# lists
ko_lookup_slim %>%
  fwrite("../data/kegg/formatted/formatted_ko_lookup.tsv", sep="\t")

####################################################################################
####################################################################################

####################################################################################
# parse the strucuture of the KEGG pathways
####################################################################################

# read csv containing pathway categories
parse_subsection_pathways <- function(x) {
  x <- x[ 2:(length(x)-1) ]
  x <- matrix(x, length(x)/2, 2, byrow=TRUE)
  ko_labels <- paste0("K0", str_extract(x[,1], "\\d+"))
  pathway_desc <- trimws(x[,2], "left")
  return(paste(ko_labels, pathway_desc))
}

raw_pathway_csv <- readLines("../data/kegg/kegg_pathway_categories.csv")

# cut the section headings
section_start_idxs <- c(which(str_detect(raw_pathway_csv, "^\\d+\\. ")), length(raw_pathway_csv)+1)
pathway_sections <- list()
for(i in 1:(length(section_start_idxs)-1)) {
  pathway_sections[[ raw_pathway_csv[[ section_start_idxs[[i]] ]] ]] <- raw_pathway_csv[ (section_start_idxs[[i]]+1):(section_start_idxs[[i+1]]-1) ]
}

# cut the subsection headings
for(j in 1:length(pathway_sections)) {
  section <- pathway_sections[[ j ]]
  subsection_start_idxs <- c(which(str_detect(section, "^\\d+\\.\\d+")), length(section)+1)
  tmp <- list()
  for(i in 1:(length(subsection_start_idxs)-1)) {
    tmp[[ section[[ subsection_start_idxs[[i]] ]] ]] <- section[ (subsection_start_idxs[[i]]+1):(subsection_start_idxs[[i+1]]-1) ]
  }
  pathway_sections[[ names(pathway_sections)[[j]] ]] <- lapply(tmp, parse_subsection_pathways)
}

names(pathway_sections) <- paste0(names(pathway_sections), "___")

pathway_sections <- unlist(pathway_sections, recursive=FALSE)
pathway_sections <- lapply(pathway_sections, as.data.frame) %>%
  rbindlist(idcol="label") %>%
  separate(label, c("section", "subsection"), "___\\.") %>%
  rename(pathways="X[[i]]")

fwrite(pathway_sections, "../data/kegg/formatted/parsed_pathway_info.csv")

####################################################################################
####################################################################################