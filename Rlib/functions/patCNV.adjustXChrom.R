patCNV.adjustXChrom <- function(session.info, covg.info) 
{
# The first time the coverage was adjusted, the X.FrmalevsMale.ratio was 2.0 , 
#  now the ratio was learned from the germline samples.
  X.FemalevsMale.ratio<-session.info$X.FemalevsMale.ratio




  haploidXExons<-which((session.info$exon_info$Chr=="chrX" ) & session.info$exon_info$PAR==0)

  male.idx<-c()
  if("sex" %in% colnames(session.info$file_info)) {
      trueSex<-as.character(session.info$file_info$sex)
      male.idx<-which(trueSex=="MALE")
  }
  if(length(male.idx)==0 || is.null(male.idx)) {
    return(covg.info)
  }
  total_count_vec<-covg.info$total_count_vec
  exon_count_mtx<-covg.info$exon_count_mtx
  exon_RPKM_mtx<-covg.info$exon_RPKM_mtx

  exon_count_mtx[haploidXExons,male.idx]  <-X.FemalevsMale.ratio*0.5*exon_count_mtx[haploidXExons,male.idx]
  exon_RPKM_mtx[haploidXExons,male.idx]  <-X.FemalevsMale.ratio*0.5*exon_RPKM_mtx[haploidXExons,male.idx]



  return(list(total_count_vec=total_count_vec,
              exon_count_mtx=exon_count_mtx,exon_RPKM_mtx=exon_RPKM_mtx))

  
}

