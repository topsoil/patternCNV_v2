patCNV.chr.coverage.QC <- function( covg_info, session_info,legend.layout="topright",plot.cex=2,legend.cex=1)
  
  
{
   N_sample <- length(covg_info$total_count_vec)
   chr.rpkm.mtx <- mat.or.vec(22,N_sample)

   for(k in 1:22)
   {
	  tmp.chr <- paste('chr',k,sep='')
	  sel.exon.idx <- which(session_info$exon_info$Chr==tmp.chr) 
	  for(j in 1:N_sample)
	  {
	    chr.rpkm.mtx[k,j] <- sum(covg_info$exon_RPKM_mtx[sel.exon.idx,j],na.rm=TRUE)
	  }  
	}  

   samplecolor_vec <- rainbow(N_sample)
	
 	plot(-1,-1,xlim=c(1,22),ylim=c(0,max(chr.rpkm.mtx,na.rm=TRUE)*1.5),xaxp=c(1,22,21),
		xlab='Chr index',ylab='RPKM sum per Chr')
	for(j in 1:N_sample)
	{
	  lines(chr.rpkm.mtx[,j],type='b',col=samplecolor_vec[j],pch=j,lwd=1.5,cex=plot.cex)  
	}
	grid()
	legend(x=legend.layout,colnames(covg_info$exon_RPKM_mtx),
	       cex=legend.cex,col=samplecolor_vec,pch=seq(1,N_sample),lty=1,lwd=1)

} # end of QC function
  