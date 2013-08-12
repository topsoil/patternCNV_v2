
patCNV.compute.CNV.single <- function( session_info, sample_name,
			   ref_type='average.pattern',
                           episl=1, small_delta=1e-5, bin_size=10,
			   is.verbose=TRUE)

# ref_type: either "basic.paired", "Germline" or "average.pattern"

# ref_avg_wigfile:    reference or average pattern wig file
# var_wigfile:	      varability pattern wig file 	
{
   
  sample.idx <- which(session_info$file_info$sample.name==sample_name)
  wig_filename <- session_info$file_info$file.name[sample.idx]
  sample_ID <- session_info$file_info$ID[sample.idx]
  
  exon_bin_vec <- session_info$exon_info$exon_bin_vec
  is_capture_vec <- session_info$exon_info$is_capture_vec

  if(ref_type=='average.pattern')
	{ref_avg_wigfile <- session_info$pattern$avg_wig}
  var_wigfile <- session_info$pattern$var_wig
  
  if(ref_type=='Germline'|ref_type=='basic.paired')	{
		sel.sample.subjectID <- session_info$file_info$subject.ID[sample.idx]
		sel.sample.germline.idx <- 
			which(session_info$file_info$subject.ID==sel.sample.subjectID &
				session_info$file_info$sample.type=='Germline')
		if (!length(sel.sample.germline.idx))	{
			stop(paste('Germline sample of subject',sel.sample.subjectID,'cannot be located'))		}
		ref_avg_wigfile <- session_info$file_info$file.name[sel.sample.germline.idx]
		germline_sample_ID <- session_info$file_info$ID[sel.sample.germline.idx]
		germline.exon_vec <- 
              patCNV.scan.covg.single(ref_avg_wigfile,germline_sample_ID,
			exon_bin_vec,is_capture_vec,is.plot=FALSE, bin_size=bin_size)

		  germline.total_count <- sum(germline.exon_vec,na.rm=TRUE)
	  	  germline.RPKM_deno <- bin_size*germline.total_count/1e9
	}

  if(is.verbose) {  cat('\n Processing',sample_ID,':\n',wig_filename,'\n',sep=' ')  }

  
  exon_vec <- 
      patCNV.scan.covg.single(wig_filename,sample_ID,exon_bin_vec,is_capture_vec,
                             is.plot=FALSE, bin_size=bin_size)
      
  total_count <- sum(exon_vec,na.rm=TRUE)  
  

   N_exons <- length(exon_bin_vec)
    
  
  RPKM_deno <- bin_size*total_count/1e9


  
  # RPKM = read/ ( (kb) * (total_counts/1e6) )
  #      = read/ (bin_size/1e3 * total_counts/1e6)
  #      = read/ (bin_size*total_counts/1e9)
  
  #================ initialize exon-summary vectors
  CNV_vec <- mat.or.vec(N_exons,1)
  #CNV_SD_vec <- mat.or.vec(N_exons,1)
   
  case_wcon <- file(wig_filename,'r')
    
  mean_wcon <- file(ref_avg_wigfile,'r')

  if(ref_type!='basic.paired')  { SD_wcon <- file(var_wigfile,'r') }
  
    
  #========= read 1st line
  v <- readLines(case_wcon,1)
  v <- readLines(mean_wcon,1)
  
  if(ref_type!='basic.paired')  { v <- readLines(SD_wcon,1) }
   
  #======== read each exons (average.pattern)
   
  if(ref_type=='average.pattern')
 {

  for (k in 1:N_exons)
  {
    N_bins <- exon_bin_vec[k]
    
    y_vec <- mat.or.vec(N_bins,1) # from case
    mean_vec <- mat.or.vec(N_bins,1)  # from ref_avg_wig
    SD_vec <- mat.or.vec(N_bins,1)    # from var_wig
    
    weight_vec <- mat.or.vec(N_bins,1) # weighted mean
    
    for ( j in 1:N_bins)
    {
      y_vec[j] <-  log2(
        (as.numeric(readLines(case_wcon,1)) + episl)/ RPKM_deno )
      
      mean_vec[j] <- as.numeric(readLines(mean_wcon,1))
      SD_vec[j] <- as.numeric(readLines(SD_wcon,1))
      
     }
  
    #========= computing CNV 
    weight_vec <- 1/(SD_vec+small_delta)^2
    CNV_vec[k] <- sum( (y_vec-mean_vec) * weight_vec )/sum(weight_vec,na.rm=TRUE)
 #   CNV_SD_vec[k] <- sqrt(1/sum(weight_vec))
     
    
    #========= skip exon first line
    v <- readLines(case_wcon,1)
    v <- readLines(mean_wcon,1)
    v <- readLines(SD_wcon,1)
  }
  
  } # end of if(ref_type=='average.pattern')


  #======== read each exons (Germline)
   
  if(ref_type=='Germline')
 {

  for (k in 1:N_exons)
  {
    N_bins <- exon_bin_vec[k]
    
    y_vec <- mat.or.vec(N_bins,1) # from case
    mean_vec <- mat.or.vec(N_bins,1)  # from ref_avg_wig
    SD_vec <- mat.or.vec(N_bins,1)    # from var_wig
    
    weight_vec <- mat.or.vec(N_bins,1) # weighted mean
    
    for ( j in 1:N_bins)
    {
      y_vec[j] <-  log2(
        (as.numeric(readLines(case_wcon,1)) + episl)/ RPKM_deno )
      
      #mean_vec[j] <- as.numeric(readLines(mean_wcon,1))
      mean_vec[j] <-  log2(
        (as.numeric(readLines(mean_wcon,1)) + episl)/ germline.RPKM_deno )
      SD_vec[j] <- as.numeric(readLines(SD_wcon,1))
      
     }
  
    #========= computing CNV 
    weight_vec <- 1/(SD_vec+small_delta)^2
    CNV_vec[k] <- sum( (y_vec-mean_vec) * weight_vec )/sum(weight_vec,na.rm=TRUE)
 #   CNV_SD_vec[k] <- sqrt(1/sum(weight_vec))
     
    
    #========= skip exon first line
    v <- readLines(case_wcon,1)
    v <- readLines(mean_wcon,1)
    v <- readLines(SD_wcon,1)
  }
  
  } # end of if(ref_type=='Germline')



  #======== read each exons (basic.paired: somatic vs. germline, no pattern information is used at all)
   
  if(ref_type=='basic.paired')
 {

  for (k in 1:N_exons)
  {
    N_bins <- exon_bin_vec[k]
    
    y_vec <- mat.or.vec(N_bins,1) # from case
    mean_vec <- mat.or.vec(N_bins,1)  # from ref_avg_wig
    
    #weight_vec <- mat.or.vec(N_bins,1) # weighted mean    
    for ( j in 1:N_bins)
    {
      y_vec[j] <-  log2(
        (as.numeric(readLines(case_wcon,1)) + episl)/ RPKM_deno )
      
      #mean_vec[j] <- as.numeric(readLines(mean_wcon,1))
      mean_vec[j] <-  log2(
        (as.numeric(readLines(mean_wcon,1)) + episl)/ germline.RPKM_deno )
      
     }
  
    #========= computing CNV 
    CNV_vec[k] <- sum( (y_vec-mean_vec) )/N_bins

     
    
    #========= skip exon first line
    v <- readLines(case_wcon,1)
    v <- readLines(mean_wcon,1)
   }
  
  } # end of if(ref_type=='basic.paired')
  

  #========= closing 
  
  close(case_wcon)
  close(mean_wcon)

  if(ref_type!='basic.paired')  {   close(SD_wcon)  }
  
  return(
	list(CNV=CNV_vec,sample.ID=sample_ID,sample.name=sample_name)
		)
  
}

