#####################################
#### Step 1. Load data.
#####################################
library(minfi)
library(ChAMP)
library(ENmix)
library(DMRcate)

#####load data
dataDir<-IDAT_DIR
targets <- read.metharray.sheet(dataDir)
rgSet <- read.metharray.exp(targets = targets,extended = TRUE,force = FALSE)






