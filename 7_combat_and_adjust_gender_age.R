############
### Step 11. ComBat after BMIQ
############

###correct for batch effect

batch1 = factor(pData(GmSet_XY)$Well)
tmp1<-ComBat.mc(tmp_BMIQ,batch = batch1,nCores = CPU_CORES_USED)
message("[ = = = = = = = = = = = = = = ]")
message("<<---Correct for Well done!--->>")
message("[ = = = = = = = = = = = = = = ]")


batch2 = factor(pData(GmSet_XY)$Array)
tmp2<-ComBat.mc(tmp1,batch = batch2,nCores = CPU_CORES_USED)
message("[ = = = = = = = = = = = = = = ]")
message("<<---Correct for Array done!--->>")
message("[ = = = = = = = = = = = = = = ]")

batch3 = factor(pData(GmSet_XY)$Slide)
tmp3<-ComBat.mc(tmp2,batch = batch3,nCores = CPU_CORES_USED)
message("[ = = = = = = = = = = = = = = ]")
message("<<---Correct for Slide done!--->>")
message("[ = = = = = = = = = = = = = = ]")


batch4 = factor(pData(GmSet_XY)$Array_Plate)
tmp4<-ComBat.mc(tmp3,batch = batch4,nCores = CPU_CORES_USED)
message("[ = = = = = = = = = = = = = = ]")
message("<<---Correct for Array_Plate done!--->>")
message("[ = = = = = = = = = = = = = = ]")



##### correct covariates: sex and age using linear model.

lm_model <- lm(t(tmp4) ~ as.factor(pData(GmSet_XY)$Sex)+pData(GmSet_XY)$Age)
tmp_lm <- t(lm_model$res)+rowMeans(tmp4)
