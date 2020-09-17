#####################################
#### Step 3. Remove low quality samples according to detection P value.
#####################################

mset<-preprocessRaw(rgSet)
gmSet<-mapToGenome(mset)

tmp = getBeta(gmSet, "Illumina",offset = 100)

detP <- detectionP(rgSet)
detP<-detP[rownames(tmp),]
detPcut = DETECTION_PVALUE_CUTOFF

tmp[detP >= detPcut] <- NA 
numfail <- matrix(colMeans(is.na(tmp)))
rownames(numfail) <- colnames(detP)
colnames(numfail) <- "Failed CpG Fraction"
outfile = paste(RESULT_QC_DIR, "/Failed_CpG_fraction.txt", sep = "")
write.table(numfail,outfile,sep = "\t")

SampleCutoff= SAMPLE_CUT_OFF

RemainSample <- which(numfail < SAMPLE_CUT_OFF) 
rgSet <- rgSet[,RemainSample]
gmSet <- gmSet[,RemainSample]
detP <- detP[,RemainSample]
tmp <- tmp[,RemainSample]



