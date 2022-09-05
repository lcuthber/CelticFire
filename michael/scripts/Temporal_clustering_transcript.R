# Load packages
library(e1071)
library(Mfuzz)

# Data import
load("./de.R") # Formal class expression set object

# cluster
centers = 8
cl_de_8 <- mfuzz(de, centers = centers, m=2) # cmeans fuzzy clustering with 8 centers
