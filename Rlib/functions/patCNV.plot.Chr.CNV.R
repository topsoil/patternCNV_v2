patCNV.plot.Chr.CNV <- function( cnv_res, session_info, sample_name, 
				sel_chr='chr1', chr_range=NULL,
				capture.only=TRUE, min_ref_avg_RPKM=3,ref_avg_type='median',
				cex=0.65,col='steelblue3',ylim=c(-3,3),
				xlab='Mb',ylab='log2-ratio',...)
#===== ... is for plot()
{
	
	sel.sample.idx <- which(cnv_res$sample.name==sample_name)
	if(!length(sel.sample.idx))
	 { stop(paste(sample_name,'cannot be located in input','\n')) }	

   if(ref_avg_type=='median')	
     {ref_avg_RPKM <- unlist(read.delim(session_info$Misc$median_RPKM_file,header=FALSE))}

   if(ref_avg_type=='mean')	
     {ref_avg_RPKM <- unlist(read.delim(session_info$Misc$mean_RPKM_file,header=FALSE))} 
    
   
      if(capture.only)
      {
        sel_chr_idx <- which(session_info$exon_info$Chr==sel_chr &
                    ref_avg_RPKM>=min_ref_avg_RPKM & session_info$exon_info$InCapture==1)  
      } else {
        sel_chr_idx <- which(session_info$exon_info$Chr==sel_chr & 			ref_avg_RPKM>=min_ref_avg_RPKM)
      } # exon indices in selected Chr

      if(!is.null(chr_range))
	{
	   range_idx <- which(session_info$exon_info$Start[sel_chr_idx]>=chr_range[1] &
				session_info$exon_info$Stop[sel_chr_idx]<=chr_range[2])	
	   sel_chr_idx <- sel_chr_idx[range_idx] 	
	}	      

      order_chr_idx <- sel_chr_idx[order(session_info$exon_info$Start[sel_chr_idx])]
      gnm_pos <- (session_info$exon_info$Start[order_chr_idx])/1e6
      plot(gnm_pos, cnv_res$CNV[order_chr_idx,sel.sample.idx], type='p',
		main=paste(sel_chr,' log2-ratio CNV @ ',cnv_res$sample.ID[sel.sample.idx]), 
		cex=cex,col=col,xlab=xlab,ylab=ylab,ylim=ylim,...)  
  
}
