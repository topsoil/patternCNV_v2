patCNV.scan.covg.single <- function(wig_filename,sample_ID,exon_bin_vec,is_capture_vec,
                               is.plot=FALSE,log.for.plot=TRUE,plot_output_DIR,
                               bin_size=10, ylab='# of exons',xlab='log10(coverage per base)')
{
  #========= It doesn't expect users to directly use this function
  
  #r_int <- as.integer(readLines(wig_filename))
  r_int <- suppressWarnings(as.integer(readLines(wig_filename)))
  # suppressWarnings(.) suppress the warning of "NAs introduced by coercion", which is casued integer conversion of wig text header	  

  N_exons <- length(exon_bin_vec)
  exon_header_vec <- which(is.na(r_int))
  
  if(length(exon_header_vec)!=N_exons)
  {
    cat('Number of exons differs between exon key file and wig file',wig_filename,'\n')
    cat('Expected # of exons from key file:',N_exons,'\n')
    cat('Observed # of exons in wig file:',length(exon_header_vec),'\n')
    cat('Please check and re-run BAM2WIG conversion if needed. \n')
    stop('Number of exons mismatch \n')
  }
  
  exon_start_vec <- exon_header_vec + 1
  exon_end_vec <- c(exon_header_vec[2:N_exons]-1,length(r_int)) # orginal
  
    
  exon_vec <- mat.or.vec(N_exons,1)
  for (k in 1:N_exons)
  {
     exon_vec[k] <- sum(r_int[exon_start_vec[k]:exon_end_vec[k]])
  }
  
    
  #============ plot exon coverage distributions
  if(is.plot)
  {
   
    patCNV.create.DIR(plot_output_DIR)
    captured_idx <- which(is_capture_vec==1)
    non_captured_idx <- which(is_capture_vec==0)
    
    
    pdf_filename <- paste(plot_output_DIR,sample_ID,'_exon_covg.pdf',sep='')
    pdf(pdf_filename)
    
    if(log.for.plot)
    {
      hist_res <- hist(log10(exon_vec/(bin_size*exon_bin_vec)), 50,
                       main= paste('All exons ',sample_ID,sep='@'), ylab=ylab,xlab=xlab )
      hist_res <- hist(log10(exon_vec[captured_idx]/(bin_size*exon_bin_vec[captured_idx])), 50,
                       main= paste('Captured exons ',sample_ID,sep='@'), ylab=ylab, xlab=xlab )
      if(length(non_captured_idx)>0)
      {
        hist_res <- hist(log10(exon_vec[non_captured_idx]/(bin_size*exon_bin_vec[non_captured_idx])), 50,
                         main= paste('Non-captured exons ',sample_ID,sep='@'), ylab=ylab, xlab=xlab )  
      }
      
    } else {
      
      hist_res <- hist((exon_vec/(bin_size*exon_bin_vec)), 50,
                       main= paste('All exons ',sample_ID,sep='@'), ylab=ylab,xlab=xlab )
      hist_res <- hist((exon_vec[captured_idx]/(bin_size*exon_bin_vec[captured_idx])), 50,
                       main= paste('Captured exons ',sample_ID,sep='@'), ylab=ylab, xlab=xlab )
      if(length(non_captured_idx)>0)
      {
        hist_res <- hist((exon_vec[non_captured_idx]/(bin_size*exon_bin_vec[non_captured_idx])), 50,
                       main= paste('Non-captured exons ',sample_ID,sep='@'), ylab=ylab, xlab=xlab )
      }
    }
    dev.off()
  } # end of plot
  
  #============= return exon counts
  return(exon_vec)
  
}

