patCNV.scan.covg.multi.bin <- function(session_info, sample.type=NULL, 
				bin_size=10, is.verbose=TRUE)
  
{
  #=========== loading configuration information
  plot_output_DIR <- session_info$DIR_info$plot_output_DIR
  txt_output_DIR <- session_info$DIR_info$txt_output_DIR
    
  exon_bin_vec <- session_info$exon_info$exon_bin_vec
 
  
  if(is.null(sample.type)) # using all the samples
  {
    wig_filename_vec <- session_info$file_info$file.name
    sample_ID_vec <- session_info$file_info$ID  
  } else {                 
    # using selected samples with given sample type
    sel_sample_idx <- which(session_info$file_info$sample.type==sample.type)
    wig_filename_vec <- session_info$file_info$file.name[sel_sample_idx]
    sample_ID_vec <- session_info$file_info$ID[sel_sample_idx]
  }
  
  
  N_exons <- length(exon_bin_vec)
  N_samples <- length(wig_filename_vec)
  #=========== end of loading configuration information
  
  total_count_vec <- mat.or.vec(N_samples,1)
  names(total_count_vec) <- sample_ID_vec

  exon_count_list <- list()
  exon_RPKM_list <- list()
  
  for(p in 1:N_exons)
   {
     tmp.mtx <- mat.or.vec(exon_bin_vec[p],N_samples) 
     colnames(tmp.mtx) <- sample_ID_vec	
	exon_count_list[[p]]  <-  exon_RPKM_list[[p]] <- tmp.mtx
   } 	 



  
  for(k in 1:N_samples) 
  {
    tmp_wig_filename <- wig_filename_vec[k]
    tmp_sample_ID <- sample_ID_vec[k]
    if (is.verbose)
    {
      print(paste('scanning ',k,'-th sample:',tmp_sample_ID))  
    }
    single.exon.list <- 
      patCNV.scan.covg.single.bin(wig_filename=tmp_wig_filename,sample_ID=tmp_sample_ID,
                              exon_bin_vec=exon_bin_vec)
     for(p in 1:N_exons)
	   {	                                            
	    exon_count_list[[p]][,k] <- single.exon.list$binlist[[p]] 
	    total_count_vec[k] <- single.exon.list$total.covg # total bp counts

	    exon_RPKM_list[[p]][,k] <- single.exon.list$binlist[[p]] *1e9/
			(bin_size*1*total_count_vec[k])

    #exon_RPKM_mtx[,k] <- (count_vec*1e9)/(bin_size*exon_bin_vec*total_count_vec[k])
      # RPKM = count_in_exon/ ( (kb) * (total_counts/1e6) )
      #      = count_in_exon/ ( (N_bins * bin_size/1e3) * total_counts/1e6)
      #      = count_in_exon / (N_bins * bin_size * total_counts/1e9)
      #      = count_in_exon*1e9 / (N_bins * bin_size * total_counts)
	   }	
  }
  
  #============= return exon counts
  return(list(total_count_vec=total_count_vec,
              exon_count_list=exon_count_list,exon_RPKM_list=exon_RPKM_list))
  
}

