patCNV.fit.null.model <- function(ctrl.cnv_res, session_info,
                            type='autosome',is.plot=FALSE)

{

  cnv.mtx <- ctrl.cnv_res$CNV	  

  if(type=='all')
  {
    res <- list()
    res$type='all'
    v1 <- as.vector(cnv.mtx[which(session_info$exon_info$is_capture_vec==1),])
    plot(density(v1),
         xlim=c(-1,1),main='all',lwd=2,xlab='',ylab='')
    res$location <- median(v1)
    res$scale <- mad(v1)
    res$SD <- sd(v1)
    return(res)
  }
  
  if(type=='autosome')
  {
    
  
  res <- list()
  res$type='autosome'
  res$Chr <- paste('chr',seq(1,22),sep='')
  res$location <- mat.or.vec(22,1)
  res$scale <- mat.or.vec(22,1) + 1
  res$SD <- mat.or.vec(22,1) + 1
  
  #====== begin plotting
  if(is.plot)
  {
    par(mfrow = c(4, 6))
    
    for ( k in 1:22)
    {
      tmp_chr <- res$Chr[k]
      plot(density(as.vector(cnv.mtx[which(session_info$exon_info$Chr==tmp_chr),])),
           xlim=c(-1,1),ylim=c(0,5),main=tmp_chr,lwd=2,xlab='',ylab='')
    }  
  } #===== end plotting
  
  
  for ( k in 1:22)
  {
    tmp_chr <- res$Chr[k]
    v1 <- as.vector(cnv.mtx[which(session_info$exon_info$Chr==tmp_chr),])
    res$location[k] <- median(v1)
    res$scale[k] <- mad(v1)
    res$SD[k] <- sd(v1)
  }

  return(res)
  
  }
  
}


