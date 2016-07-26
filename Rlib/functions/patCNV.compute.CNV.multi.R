patCNV.compute.CNV.multi <- function( session_info, sample.type=NULL,
				   ref_type='average.pattern',	
                                   episl=1, small_delta=1e-5, bin_size=10,
				     zero.median.adjust=TRUE,is.verbose=FALSE)

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
    sample_name_vec <- session_info$file_info$sample.name	  
  } else {                 
    # using selected samples with given sample type
    sel_sample_idx <- which(session_info$file_info$sample.type==sample.type)
    wig_filename_vec <- session_info$file_info$file.name[sel_sample_idx]
    sample_ID_vec <- session_info$file_info$ID[sel_sample_idx]
    sample_name_vec <- session_info$file_info$sample.name[sel_sample_idx]	  
  }
    
  N_exon <- length(exon_bin_vec)
  N_sample <- length(wig_filename_vec)
  
  
  ref_mean_wigfile <- session_info$pattern$avg_wig
  ref_SD_wigfile <- session_info$pattern$var_wig
  
  #=========== end of loading configuration information
  
  
  
  
  CNV.mtx <- mat.or.vec(N_exon,N_sample)
  colnames(CNV.mtx) <- sample_ID_vec
  
  pb <- txtProgressBar(style=3,max=N_sample)
  for (k in 1:N_sample)
  {
    
    wig_filename <- wig_filename_vec[k]
    sample_ID <- sample_ID_vec[k]
    sel.sample.name <- sample_name_vec[k]
 
    #cat('\n',k,wig_filename,sample_ID,'\n',sep=' | ')
    if(is.verbose) {  cat(k,'-th of',N_sample,'\n') }
    setTxtProgressBar(pb, k)    
	   
    cnv_res <- patCNV.compute.CNV.single ( session_info, sel.sample.name, 
					ref_type=ref_type,episl=episl,small_delta=small_delta,bin_size=bin_size,
					is.verbose=is.verbose)
    CNV.mtx[,k] <- cnv_res$CNV
  }  
  
  if(zero.median.adjust) {

           for (k in 1:N_sample)
           {
             CNV.mtx[,k] <- CNV.mtx[,k] - median(CNV.mtx[,k],na.rm=TRUE)
            }
	}

   cat('\n')	

  return(
	list(CNV=CNV.mtx,sample.ID=sample_ID_vec,sample.name=sample_name_vec)
	)
  
}

