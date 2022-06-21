library(pbmcapply)
library(KEGGREST)

chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels=FALSE)) 

lookup_kos <- function(ko_list, n_cores=1) {
  n_chunks <- (length(ko_list) %/% 10) + 1
  cat("Looking up", length(ko_list), "KOs in", n_chunks, "chunks \n")
  query <- pbmclapply(
    chunk2(ko_list, n_chunks),
    keggGet,
    mc.cores=n_cores)
  return(unlist(query, recursive=FALSE))
}

format_query <- function(q) {
  formatted_pathways <- paste0(names(q$PATHWAY), " ", q$PATHWAY, collapse=",")
  return(
    data.frame(entry=q$ENTRY,
               name=q$NAME, # definition=q$DEFINITION,
               symbol=q$SYMBOL,
               module=ifelse("MODULE" %in% names(q), q$MODULE, NA),
               pathways=formatted_pathways)
  )
}

get_slim_query <- function(df) {
  # format a raw query for single KO (df)
  df %>%
    select(KO, symbol, name, module, pathways) %>%
    # mutate(pathways=str_split(pathways, "map|map,")) %>%
    unnest(pathways) %>%
    mutate(pathways=trimws(pathways, "right", whitespace=",")) %>%
    mutate(pathways=trimws(pathways, "left", whitespace=","))
  # %>%
  #   mutate(#pathways=str_replace(pathways, "map", ""),
  #          pathways=trimws(pathways, "right", whitespace=","),
  #          pathways=ifelse(pathways!="" & pathways!=" ",
  #                          paste0("K0", pathways),
  #                          NA)) %>%
  #   mutate(pathway_KO=str_extract(pathways, "^K\\d+")) 
}