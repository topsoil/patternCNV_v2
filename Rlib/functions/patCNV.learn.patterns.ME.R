patCNV.learn.patterns.ME <- function( session_info, covg_info,
                                  sample.type=NULL, episl=1, bin_size=10,
                                  exclude_sample_name = NULL )
                              
{
  
  #============= summarizing total counts and mean/median RPKM
  total_count_vec <- covg_info$total_count_vec
  mean_RPKM_file <- session_info$Misc$mean_RPKM_file
  median_RPKM_file <- session_info$Misc$median_RPKM_file
  
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
    sel_sample_idx <- which(session_info$file_info$sample.type==sample.type &
                     !is.element(session_info$file_info$sample.name, exclude_sample_name) )
    wig_filename_vec <- session_info$file_info$file.name[sel_sample_idx]
    sample_ID_vec <- session_info$file_info$ID[sel_sample_idx]
  }
  
  RKM_vec <- bin_size*total_count_vec[sample_ID_vec]/1e9
  
  # computing RPKM based on total counts of selected samples
  # RPKM = read/ ( (kb) * (total_counts/1e6) )
  #      = read/ (bin_size/1e3 * total_counts/1e6)
  #      = read/ (bin_size*total_counts/1e9)
  
  N_exons <- length(exon_bin_vec)
  N_samples <- length(wig_filename_vec)
  cat('\n', N_samples,' samples are selected for learning patterns.\n',sep='')
  cat('\n Total counts information for selected samples: \n')
  print(total_count_vec[sample_ID_vec])
  
  
  wig.mean_RPKM_vec <- apply(covg_info$exon_RPKM_mtx[,sample_ID_vec],1,mean)
  wig.median_RPKM_vec <- apply(covg_info$exon_RPKM_mtx[,sample_ID_vec],1,median)
  write.table(x=wig.mean_RPKM_vec,file=mean_RPKM_file,quote=FALSE,row.names=FALSE,col.names=FALSE)
  write.table(x=wig.median_RPKM_vec,file=median_RPKM_file,quote=FALSE,row.names=FALSE,col.names=FALSE)
  
  #=========== end of loading configuration information
  
  
  #================ initialize mean and SD wig files
  mean_vec <- mat.or.vec(N_exons,1)
  SD_vec <- mat.or.vec(N_exons,1)
  
  mean_wig_file <- session_info$pattern$avg_wig
  SD_wig_file <- session_info$pattern$var_wig
  
  mean_wcon <- file(mean_wig_file,'w')
  SD_wcon <- file(SD_wig_file,'w')
  
  con_list <- list()
  for(z in 1:N_samples) 
  {
     con_list[[z]] <- file(wig_filename_vec[z],'r')
  }
  
  
  #========= read 1st line
  for(z in 1:N_samples) 
  { v <- readLines(con_list[[z]],1)  }
  
  writeLines(v,mean_wcon)
  writeLines(v,SD_wcon)
  
  
  #======== read each exons
  exon_x_vec <- mat.or.vec(N_samples,1)
  
  txt.pb <- txtProgressBar(min=1,max=N_exons,style=3)
  
  for (k in 1:N_exons)
  {
    setTxtProgressBar(txt.pb, k)
    
    N_bins <- exon_bin_vec[k]
    for ( j in 1:N_bins)
    {
      for(z in 1:N_samples)
      {
        exon_x_vec[z] <- log2( (as.numeric(readLines(con_list[[z]],1)) + episl)/
                              RKM_vec[z] )
      }
      est_mean <- median(exon_x_vec)    # robust estimate
      est_SD <- mad(exon_x_vec)*1.4826  # robust estimate
      writeLines(format(est_mean,digits=3),mean_wcon)
      writeLines(format(est_SD,digits=3),SD_wcon)
    }
    #========= skip exon first line
      for(z in 1:N_samples) 
        { v <- readLines(con_list[[z]],1)  }
    
        writeLines(v,mean_wcon)
        writeLines(v,SD_wcon)
  }
  
  
  
  #========= closing 
  for(z in 1:N_samples) 
  {
    close(con_list[[z]])
  }
  close(mean_wcon)
  close(SD_wcon)
  
  
}

