patCNV.scan.covg.single.bin <- function(wig_filename,sample_ID,exon_bin_vec)

{
  
  #r_int <- as.integer(readLines(wig_filename))
   
  exon.list <- list()   

  r_int <- suppressWarnings(as.integer(readLines(wig_filename)))
  # suppressWarnings(.) suppress the warning of "NAs introduced by coercion", which is casued integer conversion of wig text header	  

  N_exons <- length(exon_bin_vec)
  exon_header_vec <- which(is.na(r_int))
  exon_start_vec <- exon_header_vec + 1
  exon_end_vec <- c(exon_header_vec[2:N_exons]-1,length(r_int)) # orginal

  #exon.list$total.covg <- sum(r_int,na.rm=TRUE)
  exon.list$total.covg <- median(r_int,na.rm=TRUE)*length(r_int)
  exon.list$binlist <- list()
  
  for (k in 1:N_exons)
  {
     exon.list$binlist[[k]] <- r_int[exon_start_vec[k]:exon_end_vec[k]]
  }
  
  
  #============= return exon list
  return( exon.list)
  
}

