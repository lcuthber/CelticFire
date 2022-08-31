#### Metabolite-Gene Linking via Metabolic Networks ####
# Load packages:
library(MetaboSignal)
library(igraph)

# Get human metabolic network
MS_keggFinder(KEGG_database = "organism", organism_code = "hsa")[[1]] # identify organism code and search KEGG for it
hsa_paths_info <- MS_getPathIds(organism_code = "hsa")
paths <- as.matrix(read.table("hsa_pathways_edited.txt", header = TRUE, sep = "\t"))
metabo_paths <- paths[paths[, "Path_type"] == "metabolic", "Path_id"]
signaling_paths <- paths[paths[, "Path_type"] == "signaling", "Path_id"]
MS_network <- MS_keggNetwork(metabo_paths, signaling_paths, expand_genes = FALSE)

# Map genes (KOs) and metabolites (ALI media) to metabolic networks
genes_K_all <- KO2.2$K[-1]
genes_mapped_K_all <- data.frame("Mapped_gene" = unlist(MS_findMappedNodes(nodes = c(genes_K_all), MS_network)[1]))
cpd_mapped99 <- data.frame("Mapped_cpd" = unlist(MS_findMappedNodes(nodes = c(paste0("cpd:",Am_Kegg_Hmdb$Variables_99_info.KEGG)), MS_network)[1]))

# Calculate shortest paths between metabolites and genes
SPEgmNet <- MS_shortestPathsNetwork(MS_network,
                                    organism_code = 'hsa',
                                    source_nodes = as.character(cpd_mapped99[,1]),
                                    target_nodes = genes_mapped_K_all[,1],
                                    mode='out',
                                    type="all",
                                    distance_th = 'Inf', names = FALSE, export_cytoscape = FALSE)
View(SPEgmNet)
# Filter by input data
SPE_allmapped_nodes <- c(cpd_mapped99[,1], genes_mapped_K_all[,1])
SPEgmNet.1 <- SPEgmNet[which(SPEgmNet[,1] %in% SPE_allmapped_nodes & SPEgmNet[,2] %in% SPE_allmapped_nodes ),]
SPEgmNet1.5 <- SPEgmNet.1[grep('k_compound', SPEgmNet.1[,3]),] # KO -> KO removed
SPEgmNet1.5.named <- MS_changeNames(SPEgmNet1.5, 'hsa')
SPEgmNet2 <- SPEgmNet1.5[grep('cpd', SPEgmNet1.5[,1]),] # only the cpd -> KOs
SPEgmNet2.df <- as.data.frame(SPEgmNet2[,1:2])
SPEgmNet2.named <- MS_changeNames(SPEgmNet2, 'hsa') # 54 cpds to 82 genes
SPEgmNet2.named.df <- as.data.frame(SPEgmNet2.named[,1:2])

## Create Bipartite figure:
g1.arrows <- as.factor(unlist(SPEgmNet2.named[,3]))
levels(g1.arrows) <- c('->', '<->')
g1.arrows[grep('k_compound:irreversible', g1.arrows)] <- '->' 
g1.arrows[grep('k_compound:reversible', g1.arrows)] <- '<->'
SPEgmNet2.named.2 <- SPEgmNet2.named
SPEgmNet2.named.2[,3] <- g1.arrows
SPEgmNet.graph <- as.data.frame(SPEgmNet2.named.2[,1:3])
g1.2 <- graph.data.frame(SPEgmNet.graph, directed = T)
V(g1.2)$type <- V(g1.2)$name %in% SPEgmNet.graph$source
E(g1.2)$arrow #<- g1.arrows
shapee <- c('circle', 'square')
par(mar = c(1,1,1,1))
plot(g1.2, layout = layout_as_bipartite(g1.2)[,2:1],
     vertex.color = alpha(colorme2[3:4], 0.7)[as.numeric(V(g1.2)$type)+1],
     vertex.shape = shapee[as.numeric(V(g1.2)$type)+1],
     vertex.size = 8, 
     vertex.label.dist = 5,
     vertex.label.cex = 0.8,
     #vertex.label.family = "TT Arial",
     vertex.label.degree = pi*V(g1.2)$type,
     vertex.label.color = colorme2[3:4][as.numeric(V(g1.2)$type)+1],
     edge.arrow.size = 0.2,
     edge.color = 'black',
     edge.width = 1.5,
     edge.arrow.mode = E(g1.2)$type,
     asp = 5,
     margin = -0.1)
