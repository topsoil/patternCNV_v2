patCNV.learn.patterns <- function( session_info, covg_info,
                                  sample.type=NULL, episl=1, bin_size=10,
                                  exclude_sample_name = NULL )
                              
{
  require('matrixStats')

  
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
  
  # computing RPKM based on total bp counts of selected samples
  # RPKM = read/ ( (kb) * (total_counts/1e6) )
  #      = read/ (bin_size/1e3 * total_counts/1e6)
  #      = read/ (bin_size*total_counts/1e9)
  
  N_exons <- length(exon_bin_vec)
  N_samples <- length(wig_filename_vec)
  cat('\n', N_samples,' samples are selected for learning patterns.\n',sep='')
  cat('\n Total bp counts information for selected samples: \n')
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
  
  
  #=== 
  z <- 1
  r_str <- readLines(wig_filename_vec[z])
  #r_int <- as.integer(readLines(wig_filename_vec[z]))
  r_int <- suppressWarnings(as.integer(readLines(wig_filename_vec[z])))
  # suppress the warning of "NAs introduced by coercion", induced by converting text header of wig to integer

  N_wig_bins <- length(r_int)
  N_exons <- length(exon_bin_vec)
  exon_header_vec <- which(is.na(r_int))
  exon_start_vec <- exon_header_vec + 1
  exon_end_vec <- c(exon_header_vec[2:N_exons]-1,length(r_int))
  
  wig_bin_mtx <- mat.or.vec(N_wig_bins,N_samples)
  wig_bin_mtx[,z] <- log2( (r_int + episl)/ RKM_vec[z])
  for(z in 2:N_samples) 
  {
    #r_int <- as.integer(readLines(wig_filename_vec[z]))
    r_int <- suppressWarnings(as.integer(readLines(wig_filename_vec[z])))
    # suppress the warning of "NAs introduced by coercion", induced by converting text header of wig to integer
    wig_bin_mtx[,z] <- log2( (r_int + episl)/ RKM_vec[z])
    cat(z,'-th sample: ', names(RKM_vec)[z],'\n')
    #print(proc.time()-ptm)
  }
  
  
  
  #avg_pattern_vec <-  apply(wig_bin_mtx,1,median)      # robust estimate

  median_vec <- rowMedians(wig_bin_mtx)
  avg_pattern_vec <- format(median_vec,digits=2,trim=TRUE)

  #var_pattern_vec <-  apply(wig_bin_mtx,1,mad)*1.4826  # robust estimate of SD
 

  var_pattern_vec <- format(rowMads(wig_bin_mtx)*1.4826,digits=2,trim=TRUE)

  
  
  avg_pattern_vec[exon_header_vec] <- r_str[exon_header_vec]
  var_pattern_vec[exon_header_vec] <- r_str[exon_header_vec]
  #log2( (as.numeric(readLines(con_list[[z]],1)) + episl)/
  #        RKM_vec[z] )
  
  
  #========= writing
  
  write.table(x=avg_pattern_vec,file=mean_wig_file,quote=FALSE,row.names=FALSE,col.names=FALSE)
  cat('outputing average pattern file... \n')

  write.table(x=var_pattern_vec,file=SD_wig_file,quote=FALSE,row.names=FALSE,col.names=FALSE)
  cat('outputing variability pattern file... \n')

  
  
}

