patCNV.scan.covg.multi <- function(session_info, sample.type=NULL, bin_size=10, is.verbose=TRUE, is.plot=FALSE)
  
{
  #=========== loading configuration information
  plot_output_DIR <- session_info$DIR_info$plot_output_DIR
  txt_output_DIR <- session_info$DIR_info$txt_output_DIR
    
  exon_bin_vec <- session_info$exon_info$exon_bin_vec
  is_capture_vec <- session_info$exon_info$is_capture_vec
  
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
  exon_count_mtx <- mat.or.vec(N_exons,N_samples)
  colnames(exon_count_mtx) <- sample_ID_vec
  exon_RPKM_mtx <- mat.or.vec(N_exons,N_samples)
  colnames(exon_RPKM_mtx) <- sample_ID_vec

  
  for(k in 1:N_samples) 
  {
    tmp_wig_filename <- wig_filename_vec[k]
    tmp_sample_ID <- sample_ID_vec[k]
    if (is.verbose)
    {
      print(paste('scanning ',k,'-th sample:',tmp_sample_ID))  
    }
    count_vec <- 
      patCNV.scan.covg.single(wig_filename=tmp_wig_filename,sample_ID=tmp_sample_ID,
                              exon_bin_vec=exon_bin_vec,is_capture_vec=is_capture_vec,
                              bin_size=bin_size,is.plot=is.plot,plot_output_DIR=plot_output_DIR)
                                            
    exon_count_mtx[,k] <- count_vec
    total_count_vec[k] <- sum(count_vec) # total counts
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

