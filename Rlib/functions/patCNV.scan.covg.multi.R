patCNV.scan.covg.multi <- function(session_info, sample.type=NULL, bin_size=10, is.verbose=TRUE, is.plot=FALSE)

{
  #=========== loading configuration information
  plot_output_DIR <- session_info$DIR_info$plot_output_DIR
  txt_output_DIR <- session_info$DIR_info$txt_output_DIR

  nonXY<-which((session_info$exon_info$Chr!="chrX" & session_info$exon_info$Chr!="chrY") |  session_info$exon_info$PAR==1)


  male.sample.idx<-c()
  trueSex<-NULL
  if("sex" %in% colnames(session_info$file_info)) {
      trueSex<-as.character(session_info$file_info$sex)
      male.sample.idx<-which(trueSex=="MALE")
  }

  haploidXExons<-which((session_info$exon_info$Chr=="chrX" ) & session_info$exon_info$PAR==0)
  haploidYExons<-which((session_info$exon_info$Chr=="chrY") & session_info$exon_info$PAR==0)
  print(paste("Have ",length(haploidXExons) +length(haploidYExons)," haploid exons\n",sep=""))


  exon_bin_vec <- session_info$exon_info$exon_bin_vec
  is_capture_vec <- session_info$exon_info$is_capture_vec

  if(is.null(sample.type)) # using all the samples
  {
    wig_filename_vec <- session_info$file_info$file.name
    sample_ID_vec <- session_info$file_info$ID
    sample_name_vec <- session_info$file_info$ID
  } else {                 
    # using selected samples with given sample type
    sel_sample_idx <- which(session_info$file_info$sample.type==sample.type)
    wig_filename_vec <- session_info$file_info$file.name[sel_sample_idx]
    sample_ID_vec <- session_info$file_info$ID[sel_sample_idx]
    sample_name_vec <- session_info$file_info$ID[sel_sample_idx]
    if(!is.null(trueSex)) {
        male.sample.idx<-which(trueSex[sel_sample_idx]=="MALE")
    }
  }

  print(paste("For sample type==",sample.type,", Using this Male Sample=",wig_filename_vec[male.sample.idx],sep="\n"))
 

  N_exons <- length(exon_bin_vec)
  N_samples <- length(wig_filename_vec)
  #=========== end of loading configuration information
  
  total_count_vec <- mat.or.vec(N_samples,1)
  names(total_count_vec) <- sample_ID_vec
  exon_count_mtx <- mat.or.vec(N_exons,N_samples)
  colnames(exon_count_mtx) <- sample_ID_vec
  exon_RPKM_mtx <- mat.or.vec(N_exons,N_samples)
  colnames(exon_RPKM_mtx) <- sample_ID_vec

  
  for(k in 1:N_samples) 
  {
    tmp_wig_filename <- wig_filename_vec[k]
    if (is.verbose)
    {
      print(paste('scanning ',k,'-th sample:',sample_name_vec[k]))
    }
    count_vec <- 
      patCNV.scan.covg.single(wig_filename=tmp_wig_filename,sample_ID=sample_name_vec[k],
                              exon_bin_vec=exon_bin_vec,is_capture_vec=is_capture_vec,
                              bin_size=bin_size,is.plot=is.plot,plot_output_DIR=plot_output_DIR)
                                            
    if(k %in% male.sample.idx) {
      if(length(haploidYExons)>0) {
          count_vec[haploidYExons]<- 2.0*count_vec[haploidYExons] 
      }
      if(length(haploidXExons)>0) {
#           print(paste("correcting male sample ", sample_name_vec[k]," in scan.covg.multi\n"))
          count_vec[haploidXExons]<- 2.0*count_vec[haploidXExons]
      }
    }
    exon_count_mtx[,k] <- count_vec

#    notCNV<-which(

    total_count_vec[k] <- sum(count_vec[nonXY],na.rm=TRUE) # total bp counts

    print(paste("Total_Coverage: ",sample_name_vec[k]," : ",total_count_vec[k],sep=""))

    exon_RPKM_mtx[,k] <- (count_vec*1e9)/(bin_size*exon_bin_vec*total_count_vec[k])
      # RPKM = count_in_exon/ ( (kb) * (total_counts/1e6) )
      #      = count_in_exon/ ( (N_bins * bin_size/1e3) * total_counts/1e6)
      #      = count_in_exon / (N_bins * bin_size * total_counts/1e9)
      #      = count_in_exon*1e9 / (N_bins * bin_size * total_counts)
  }
  
  #============= return exon counts
  return(list(total_count_vec=total_count_vec,
              exon_count_mtx=exon_count_mtx,exon_RPKM_mtx=exon_RPKM_mtx))
  
}
