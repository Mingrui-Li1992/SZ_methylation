#####################################
#### Step 10. BMIQ normalization
#####################################

tmp_BMIQ<-champ.norm(beta = tmp_XY,
                     mset = GmSet_XY,
                     resultsDir = RESULT_BMIQ_DIR, 
                     method = "BMIQ",
                     plotBMIQ = TRUE,
                     arraytype = "EPIC",
                     cores = CPU_CORES_USED)


