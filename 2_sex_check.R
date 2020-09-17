############################
########Step 2: sex check
###########################

sampleNames(rgSet)=rgSet[[1]]
mset<-preprocessRaw(rgSet)
GMSet<-mapToGenome(mset)
sex<-getSex(GMSet, cutoff = -2)
sex_out<-paste(RESULT_QC_DIR,"/sex_check_results.txt",sep="")
write.table(sex,sex_out,sep="\t")
