library (ggplot2)
library (reshape2)
library(ggforce)
data <- read.table (file = "CelticFire/MIGA/novelty.r", header = T, sep = "\t")
ddata <- melt (data, value.name = "Significance", id.vars = "ID", variable.name = "Level")
ddata$Level <- factor (ddata$Level, levels = c ("Root", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Subspecies", "Dataset"))
ddata$Hypothesis <- "Novel"

adata <- read.table (file = "tax_novelty_miga.txt", header = T, sep = "\t")
addata <- melt (adata, value.name = "Significance", id.vars = "ID", variable.name = "Level")
addata$Level <- factor (addata$Level, levels = c ("Root", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Subspecies", "Dataset"))
addata$Hypothesis <- "Known"

laddata <- rbind (ddata, addata)
laddata <- subset (laddata, Level != "Dataset")

xx <- 12

ggplot (laddata, aes (x = Level, y = Significance)) + geom_jitter (height = 0, size = 3, alpha = 0.5, aes (color = Hypothesis)) + theme (axis.text.x = element_text (angle = -90, hjust = 0, size = xx), axis.text.y = element_text (size = xx), axis.title.y = element_text (size = xx), axis.title.x = element_text (size = xx), legend.position = c (0.2, 0.5)) + xlab ("") + ylab ("P-value") # + scale_y_continuous (trans = "log1p")

ggsave("novelty.pdf",plot=last_plot(),dpi=600)
