###########################
#### Filter_bad_snp_mulitHit_ChrXY_probes.
#####################################
#### Step 5. filter probes with >0.01 detP in more than 5% samples.
#####################################
ProbeCutoff = PROBE_CUT_OFF 
GmSet<-GmSet_FunNorm[rownames(gmSet),]
GmSet_detP = GmSet[rowSums(is.na(tmp)) <= ProbeCutoff*ncol(detP),]
tmp_detP = getBeta(GmSet_detP, "Illumina",offset = 100)

#####################################
###### Step 6. filter out probes with less than 3 beads in at least 5% of samples per probe.
#####################################

### mybeadcount function from ChAMP R package https://bioconductor.org/packages/release/bioc/vignettes/ChAMP/inst/doc/ChAMP.html

mybeadcount <- function(x)
{
  #select out bead count data frame
  getNBeads(x) -> nb
  locusNames <- getManifestInfo(x, "locusNames")
  bc_temp <- matrix(NA_real_, ncol = ncol(x), nrow = length(locusNames),
                    dimnames = list(locusNames, sampleNames(x)))
  
  TypeII.Name <- getProbeInfo(x, type = "II")$Name
  bc_temp[TypeII.Name, ] <- nb[getProbeInfo(x, type = "II")$AddressA,]
  
  TypeI <- getProbeInfo(x, type = "I")
  
  bc_temp->bcB
  bc_temp->bcA    
  
  bcB[TypeI$Name, ] <- nb[TypeI$AddressB,]
  bcA[TypeI$Name, ] <- nb[TypeI$AddressA,]
  
  which(bcB< BEAD_COUNTS) -> bcB3
  which(bcA< BEAD_COUNTS) -> bcA3
  bcA->bcA2
  bcB->bcB2
  bcA2[bcA3]<-NA
  bcA2[bcB3]<-NA
  
  data.frame(bcA2)->bc
  bc
}

bc = mybeadcount(rgSet)
beadCutoff = BEAD_CUTOFF 
bc2 = bc[rowSums(is.na(bc)) < beadCutoff*(ncol(bc)),]

remain_cpg<-row.names(bc2)
GmSet_BC = GmSet_detP[featureNames(GmSet_detP) %in% remain_cpg,]
tmp_BC = getBeta(GmSet_BC, "Illumina",offset = 100)

#####################################
#### Step 7. remove SNP related CpGs
#####################################

data(EPIC.manifest.hg19)
maskname <- rownames(EPIC.manifest.hg19)[which(EPIC.manifest.hg19$MASK_general == TRUE)]
GmSet_SNP = GmSet_BC[!featureNames(GmSet_BC) %in% maskname,]
tmp_SNP = getBeta(GmSet_SNP, "Illumina",offset = 100)

#####################################
##### Step 8. filter multi-hit probes
#####################################
data(multi.hit)
GmSet_multi <-GmSet_SNP[!featureNames(GmSet_SNP) %in% multi.hit$TargetID,]
tmp_multi = getBeta(GmSet_multi, "Illumina",offset = 100)

#####################################
#### Step 9. filter out all probes located in chromosome X and Y.
#####################################
annotation <- getAnnotation(GmSet_multi)
sex_probe <- rownames(annotation)[annotation$chr %in% c("chrX", "chrY")]
GmSet_XY <- GmSet_multi[!featureNames(GmSet_multi) %in% sex_probe, ]
tmp_XY<-getBeta(GmSet_XY, "Illumina",offset = 100)


