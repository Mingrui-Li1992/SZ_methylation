#####################################
##### Step 4. Background, dye-bias correction and functional quantile normalization
#####################################

## Beta before Functional Normalization

message("\n Step 3: Start Functional Normalization!\n")

GmSet_FunNorm<-preprocessFunnorm(rgSet, nPCs = 2, bgCorr = TRUE, dyeCorr = TRUE, keepCN = TRUE, ratioConvert = FALSE, verbose = TRUE)


