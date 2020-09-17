#####################################
#### Step 12. DMP identification
#####################################

pheno <- pData(GmSet_XY)$Sample_Group
dmp <- dmpFinder(tmp_lm, pheno = pheno, type = "categorical",qCutoff = 1, shrinkVar = FALSE)


