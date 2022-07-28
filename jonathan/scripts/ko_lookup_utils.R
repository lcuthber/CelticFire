library(pbmcapply)
library(KEGGREST)
library(yaml)

chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels=FALSE)) 

lookup_kos <- function(ko_list, n_cores=1) {
  # maximum number of lookups at a time in KEGGREST is 10
  if(n_cores>10) {
    n_chunks <- (length(ko_list) %/% 10) + 1
  } else {
    n_chunks <- n_cores
  }
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

parse_brite <- function(query) {
  
  brite_terms <- query$BRITE
  k_number <- query$ENTRY
  
  if(brite_terms[[1]]=="KEGG Orthology (KO) [BR:ko00001]") {
    ko_end_idx <- which(str_count(c(brite_terms, ""), "\\G ")==0)[[2]] - 1
    final_term <- brite_terms[[ ko_end_idx ]]
    if(!grepl(k_number, brite_terms[[ ko_end_idx ]])) {
      warning(sprintf("%s does not contain %s", brite_terms[[ ko_end_idx ]], k_number), .call=FALSE)
    } else{
      relevant_terms <- brite_terms[ 2:(ko_end_idx-1) ]
      relevant_terms <- relevant_terms[ !grepl(k_number, relevant_terms) ]
      hierarchy_level <- str_count(relevant_terms, "\\G ")
      
      inserted_ws <- 2 * (hierarchy_level - 1)
      inserted_ws <- lapply(inserted_ws, function(x) rep(" ", x)) %>%
        lapply(function(x) paste0(x, collapse=""))
      
      tmp_str <- paste0(
        inserted_ws, 
        relevant_terms %>% trimws() %>% str_remove_all(":"),
        ifelse(hierarchy_level==max(hierarchy_level), ": 0", ": ")) %>%
        paste0(collapse="\n")
      # writeLines(tmp_str)
      
      tmp_yaml <- paste0("\n", tmp_str) %>%
        yaml.load()
      
      parsed_yaml <- lapply(tmp_yaml, function(x) unlist(x, recursive=TRUE) %>% names())
      if(class(parsed_yaml)=="list") {
        parsed_yaml <- parsed_yaml %>% unlist()
      }
      parsed_yaml <- setNames(str_replace_all(parsed_yaml, "\\.", "__"), names(parsed_yaml))
      parsed_yaml <- paste0(names(parsed_yaml), "--", parsed_yaml, collapse="____") 
      
      tbl <- tibble(kegg_label=str_split(parsed_yaml, "____")) %>%
        unnest(kegg_label) %>%
        mutate(split_label=str_split(kegg_label, "--|__")) 
      
      max_terms <- sapply(tbl$split_label, length) %>% max()
      new_names <- paste0("term", 1:max_terms)
      
      tbl <- tbl %>%
        select(-kegg_label) %>%
        unnest_wider(split_label, names_repair=~new_names, names_sep="") %>%
        mutate(KO=k_number)
      
      return(tbl)
    }
  } else {
    cat("Bad start\n")
    tbl <- tibble()
  }
  
  return(tbl)
}