#####################################
#### Step 13. DMR selection!
#####################################

###convert beta matrix into GenomicMethylSet
grSet<-makeGenomicRatioSetFromMatrix(tmp_lm,array = "IlluminaHumanMethylationEPIC",annotation = "ilm10b4.hg19", what = "Beta")

#### differentially methylated regions (DMRs)
pheno <- pData(GmSet_XY)$Sample_Group
designMatrix <- model.matrix(~ pheno)

#########################################
myannotation <- cpg.annotate("array", grSet, arraytype = "EPIC",analysis.type="differential", design=designMatrix, coef=2)
dmrcoutput <- dmrcate(myannotation, lambda=1000, C=2)
results.ranges <- extractRanges(dmrcoutput, genome = "hg19")

