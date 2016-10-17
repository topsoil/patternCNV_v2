patCNV.adjustGuessedMale<-function(session.info, covg.info,guessedMales)  {

#  X.FemalevsMale.ratio<-session.info$X.FemalevsMale.ratio

  diploid.subset<-which((session.info$exon_info$Chr!="chrX" & session.info$exon_info$Chr!="chrY"  ) | session.info$exon_info$PAR==1)
  haploid.subset<-which((session.info$exon_info$Chr=="chrX" | session.info$exon_info$Chr=="chrY"  ) & session.info$exon_info$PAR==0)

#
# By default, cog.info corrects counts for MALES, when given in the sample.info file.
#  if the sex is not given and the SexQC part of (patCNV.learnExon) thinks the sample looks male,
#  then the code will set guessedMales to TRUE. The program will never override given gender (e.g. if the sample is given as FEMALE, but looks male).
#

   if(length(guessedMales)>0) {
      print(paste("adjusting Males  for samples:",paste(colnames(covg.info$exon_RPKM_mtx)[guessedMales],collapse=","),sep=""))
      covg.info$exon_RPKM_mtx[haploid.subset,guessedMales]<-  2*covg.info$exon_RPKM_mtx[haploid.subset,guessedMales]
      covg.info$exon_count_mtx[haploid.subset,guessedMales]<- 2*covg.info$exon_count_mtx[haploid.subset,guessedMales]
  }
# note covg.info passed by value
  return(covg.info)
}