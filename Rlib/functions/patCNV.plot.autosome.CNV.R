patCNV.plot.autosome.CNV <- function( session_info, cnv_res, sample_name, chr_ver='hg19',
                                capture.only=TRUE, min_ref_avg_RPKM=3,ref_avg_type='median',
                                ylim=c(-3,3),cex=0.6,	xlab='Mb',ylab='log2-ratio',
                                color_vec=c('red','blue','green','orange','brown','purple','black'),
				                        ... )
{
   if(chr_ver=='hg19')
   {
     patCNV.data(data_name='hg19.Chr_length',is.verbose=FALSE)
     chr_length_vec <- hg19.Chr_length_vec  
   }
   
    
   sel_sample_idx <- which(cnv_res$sample.name==sample_name)
   
   if(ref_avg_type=='median')	
     {ref_avg_RPKM <- unlist(read.delim(session_info$Misc$median_RPKM_file,header=FALSE))}

   if(ref_avg_type=='mean')	
     {ref_avg_RPKM <- unlist(read.delim(session_info$Misc$mean_RPKM_file,header=FALSE))} 
    
   if(ref_avg_type=='none')
     { ref_avg_RPKM <- mat.or.vec(nrow(session_info$exon_info),1)+123 }

  chr_color_vec <- rep(color_vec,ceiling(22/length(color_vec))) # repeat color patterns
  
  chr_start_vec <- chr_length_vec*0
  chr_start_vec[1] <- 0
  for(k in 2:22)
  {
    chr_start_vec[k] <- chr_start_vec[k-1] + chr_length_vec[k-1]
  }
  
 
    #sel_sample_idx <- j
    
    
    plot(-1e3,-1e3,xlim=c(0,(chr_start_vec[22]+chr_length_vec[22])/1e6),
         ylim=ylim,xlab='Mb',ylab='log2-ratio',
         main=paste('Whole exome CNV of',cnv_res$sample.ID[sel_sample_idx]))
    for(k in 1:22)
    {
      sel_chr <- paste('chr',k,sep='')
        if(capture.only)
        {
          sel_chr_idx <- which(session_info$exon_info$Chr==sel_chr &
                                ref_avg_RPKM>=min_ref_avg_RPKM & 
				session_info$exon_info$InCapture==1)  
        } else {
          sel_chr_idx <- which(session_info$exon_info$Chr==sel_chr & 	ref_avg_RPKM>=min_ref_avg_RPKM)
        } # exon indices in selected Chr
        
      order_chr_idx <- sel_chr_idx[order(session_info$exon_info$Start[sel_chr_idx])]
      gnm_pos <- (session_info$exon_info$Start[order_chr_idx]+chr_start_vec[k])/1e6
      lines(gnm_pos, cnv_res$CNV[order_chr_idx,sel_sample_idx], type='p',
            cex=cex,col=chr_color_vec[k] )  
    }

   
  
  
}
