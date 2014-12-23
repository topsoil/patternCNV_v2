patCNV.coverage.QC <- function(  session_info,covg_info,sel_genes,legend.layout="topleft",legend.cex=0.5)
# if legend.layout=='none', not plotting legend
  
{
  
   sel.exon.idx <- 
	which(is.element(session_info$exon_info$Genes,sel_genes)) 
	
   comb_exon_names <- paste(paste(session_info$exon_info$Genes[sel.exon.idx],session_info$exon_info$Chr[sel.exon.idx]),session_info$exon_info$Start[sel.exon.idx],sep=':')
	
	   N_sample <- length(covg_info$total_count_vec)
	   N_sel_exon <- length(comb_exon_names)
	   sel.rpkm.mtx <- mat.or.vec(N_sel_exon,N_sample)
   
	  for(j in 1:N_sample)
	  {
	    sel.rpkm.mtx [,j] <- covg_info$exon_RPKM_mtx[sel.exon.idx,j]
	  }  


   samplecolor_vec <- rainbow(N_sample)

	par(oma = c(10, 0, 0, 0))
 	plot(-1,-1,xlim=c(1,N_sel_exon),ylim=c(0,1.6*max(sel.rpkm.mtx)),
		xlab='',ylab='RPKM per exon',xaxt="none")
	axis(1, at=seq(1,N_sel_exon),
	labels=comb_exon_names, las = 2)
	for(j in 1:N_sample)
	{
	  lines(sel.rpkm.mtx[,j],type='b',col=samplecolor_vec[j],pch=j%%26,lwd=1.5,cex=0.8)  
	}
	grid()
  if(legend.layout!='none') 
  {
    legend(x=legend.layout,colnames(covg_info$exon_RPKM_mtx),
           cex=legend.cex,col=samplecolor_vec,pch=seq(1,N_sample)%%26,lty=1,lwd=1)  
  }
	
	
	par(oma = c(0, 0, 0, 0)) # restore
} # end of QC function
  