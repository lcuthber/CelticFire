#### Import Data ####
Egfilenames <- list.files(pattern = "*.csv", 
                          path = "../emapper_outputs")
temp <- list.files(pattern = "*.csv")
EggNOG_all <- lapply(temp, function(x) as.data.frame(read.delim(x,	header	=	FALSE,	check.names	=	FALSE)))
paste(EggNOG_all[[1]][9][1:20,])
# extract all 'ko:....' from string and arrange into separate object
which(paste(EggNOG_all[[1]][9][1:20,]) %in% paste(EggNOG_all[[2]][9][1:20,]))
EgNg <- lapply(EggNOG_all, function(x){colnames(x) = x[4,]
x = x[-c(1:4),]})

# Extract KO data ####
KO <- lapply(EgNg, function(x){
  KOs.temp <- (x[["KEGG_ko"]])
  KOs <- unlist(str_split(KOs.temp, ","))
  KOs.1 <- count(data.frame('KO' = KOs), data.frame('KO' = KOs))
}) 

KO2 <- lapply(KO, `[`, c('KO'))
KO2.1 <- count(data.frame('KO' = unlist(KO2)), data.frame('KO' = unlist(KO2))) #Collect all KOs into one data.frame and remove duplicates by counting
KO2.2 <- data.frame(KO2.1, 'K' = gsub('ko:', '', KO2.1$KO))
KO2_126 <- KO2.2[which(KO2.2$n == "126"),] #Filter for the genes in all isolates
KO2_iso_gn <- KO2
KO2_iso_gn2 <- KO2_iso_gn %>% 
  set_names(seq_along(.)) %>%
  enframe %>%
  unnest(cols = c(value))
KO2_iso_gn2.2 <- as.data.frame(KO2_iso_gn2)
KO2_iso_gn2.3 <- transform(KO2_iso_gn2.2, name = as.numeric(name))
KO2_iso_gn_ct <- as.data.frame(table(KO2_iso_gn2.3[,1], KO2_iso_gn2.3[,2]))
KO2_iso_gn_ct.1 <- KO2_iso_gn_ct[which(KO2_iso_gn_ct$Var2 != ''),] # remove empty
KO2_iso_gn_ct.cast <- dcast(KO2_iso_gn_ct.1, Var1~Var2)
# Attempt to add isolate names in for numbering:
KO2_iso_gn_ct.cast_names <- data.frame('Isolate' = Egfilenames, KO2_iso_gn_ct.cast)

# KO with KO names (MS_changeNames)
KO2.2.names <- MS_changeNames(KO2.2$K, 'hsa')
KO2.3 <- data.frame('Ko' = KO2.2, 'Gene' = KO2.2.names)
KO2.3.1 <- KO2.3[-1,]; KO2.3.1 <- KO2.3.1[,-c(2)]
KO2GeneNames <- KO2.3.1
names(KO2GeneNames) <- c("Ko.KO", "KO", "Gene")

# Import Isolate data: ####
IsolateNames_oldandnew <- read.csv(file = '../old_and_new_isolate_names.csv')
IsolateClusterLabs2 <- read.csv("../isolate_cluster_labels.csv")
