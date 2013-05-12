patCNV.estimate.FDR <- function( cnv_res, session_info, null_model, is.plot=FALSE )

{
  
  cnv.mtx <- as.matrix(cnv_res$CNV)
  N_sample <- ncol(cnv.mtx)
  N_exon <- nrow(cnv.mtx)
  
  res <- cnv_res
  res$LFDR.mtx <- mat.or.vec(N_exon,N_sample) + 1
  res$qval.mtx <- mat.or.vec(N_exon,N_sample) + 1
  res$rawp.mtx <- mat.or.vec(N_exon,N_sample) + 1
  
  colnames(res$LFDR.mtx)=colnames(res$qval.mtx)=colnames(res$rawp.mtx)=cnv_res$sample.ID
  
  if(null_model$type=='all')
  {
    for (k in 1:N_sample)
    {
        sel_vec <- cnv.mtx[,k]
        #p_tmp <- plaplace(q=sel_vec,location=Null_model$location,scale=Null_model$scale)
        #p_vec <- 2*pmin(1-p_tmp,p_tmp)
	p_vec <- patCNV.lap.pval(q=sel_vec,location=null_model$location,scale=null_model$scale)
        
        fdr_res <- fdrtool(x=p_vec,statistic='pvalue',plot=FALSE,verbose=FALSE)
        
        res$rawp.mtx[,k] <- p_vec
        res$qval.mtx[,k] <- fdr_res$qval
        res$LFDR.mtx[,k] <- fdr_res$lfdr
    }  
      
      if(is.plot)       { hist(p_vec,50) }
    return(res)
  }
    
  if(null_model$type=='autosome')
  {
    for (k in 1:N_sample)
    {
      for(j in 1:22)
      {
        sel_idx <- which(session_info$exon_info$Chr==null_model$Chr[j])
        sel_vec <- cnv.mtx[sel_idx,k]
        #p_tmp <- plaplace(q=sel_vec,location=null_model$location[j],scale=null_model$scale[j])
        #p_vec <- 2*pmin(1-p_tmp,p_tmp)
        p_vec <- patCNV.lap.pval(q=sel_vec,location=null_model$location[j],scale=null_model$scale[j])
        fdr_res <- fdrtool(x=p_vec,statistic='pvalue',plot=FALSE,verbose=FALSE)
        
        res$rawp.mtx[sel_idx,k] <- p_vec
        res$qval.mtx[sel_idx,k] <- fdr_res$qval
        res$LFDR.mtx[sel_idx,k] <- fdr_res$lfdr
      }
      
      if(is.plot)
      { hist(p_vec,50) }
      
    }
    
    return(res)
  }
  
}


