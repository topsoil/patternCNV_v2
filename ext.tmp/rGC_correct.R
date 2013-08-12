rGC_correct <- function(chr,loc,counts,gc,sample.name='sample_1',
                        N_gc_interval=1e3,robust.q=0,
                        plot_DIR='',is.pdf=TRUE,plot.cex=0.6)
{
  
  require(limma)
  require(gplots)
   
  chr <- as.vector(chr)
  loc <- as.vector(loc)
  counts <- as.vector(counts)
  gc <- as.vector(gc)
  
  ori_median_counts <- median(counts,na.rm=TRUE)
 # robust_interval <- min(abs(ori_median_counts-
 #            quantile(counts,probs=c(robust.q,1-robust.q),na.rm=TRUE)))
  
  robust_interval <- range(quantile(counts,probs=c(robust.q,1-robust.q),na.rm=TRUE))
  sel_bin_idx <- which( counts<=(robust_interval[2]) &
                         counts>=(robust_interval[1]) &
                         !is.na(counts) )
  
  lo.gc <- gc[sel_bin_idx]
  lo.count <- counts[sel_bin_idx]
  lo.fit <- loessFit( y=lo.count, x=lo.gc)
  
    
  gcmin <- min(gc,na.rm=TRUE)
  gcmax <- max(gc,na.rm=TRUE)
  gc.interval <- c(seq(gcmin,gcmax,(gcmax-gcmin)/N_gc_interval))
  gc.shift <- mat.or.vec(length(gc.interval)-1,1)
  for(j in 1:N_gc_interval)
  {
    tmp.sel.idx <- which(lo.gc>=gc.interval[j] & lo.gc<gc.interval[j+1])
    gc.shift[j] <- 
      median(lo.fit$fitted[tmp.sel.idx])
  }
  
  gc.loess.model <- loess(gc.shift~gc.interval[1:N_gc_interval])
  corrected_vec0 <- 
    counts - predict(gc.loess.model,gc)  
  #corrected_vec <- pmax(corrected_vec,0)
  corrected_vec <- 
    (corrected_vec0 - median(corrected_vec0,na.rm=TRUE)) +  ori_median_counts
  
  mdf.lo.fit <- loessFit( y=corrected_vec[sel_bin_idx], x=lo.gc)
  mdf.gc.shift <- mat.or.vec(length(gc.interval)-1,1)
  for(j in 1:N_gc_interval)
  {
    tmp.sel.idx <- which(lo.gc>=gc.interval[j] & lo.gc<gc.interval[j+1])
    mdf.gc.shift[j] <- 
      median(mdf.lo.fit$fitted[tmp.sel.idx])
  }
  mdf.gc.loess.model <- loess(mdf.gc.shift~gc.interval[1:N_gc_interval])
  
  #=========== 
  
  max.y <- quantile(corrected_vec,0.99,na.rm=TRUE)*3.5
  
  
  pdf_filename <- paste(plot_DIR,'/',sample.name,'_GC_coverage.pdf',sep='')
  if(is.pdf){   pdf(pdf_filename)   }
  smoothScatter(lo.gc,lo.count,
                ylim=c(robust_interval[1],robust_interval[2]),
                xlab='GC content',ylab='Coverage',
                main=paste('GC-coverage (before correction) \n',sample.name,''))
  lines(gc.interval[1:N_gc_interval],gc.shift,col='red',lwd=2.5,type='l')
  
  smoothScatter(gc[sel_bin_idx],corrected_vec[sel_bin_idx],
                ylim=c(robust_interval[1],robust_interval[2]),
                xlab='GC content',ylab='Coverage',
                main=paste('GC-coverage (after correction) \n',sample.name,''))
  lines(gc.interval[1:N_gc_interval],mdf.gc.shift,col='red',lwd=2.5,type='l')
  
  smoothScatter(lo.gc,lo.count,
                ylim=c(0,max.y),
                xlab='GC content',ylab='Coverage',
                main=paste('GC-coverage (before correction) \n',sample.name,''))
  lines(gc.interval[1:N_gc_interval],gc.shift,col='red',lwd=2.5,type='l')
  
  smoothScatter(gc[sel_bin_idx],corrected_vec[sel_bin_idx],
                ylim=c(0,max.y),
                xlab='GC content',ylab='Coverage',
                main=paste('GC-coverage (after correction) \n',sample.name,''))
  lines(gc.interval[1:N_gc_interval],mdf.gc.shift,col='red',lwd=2.5,type='l')
  
  for(chr.idx in 1:22)
  {
    par(mfrow=c(2,1))
    tmp.chr <- paste('chr',chr.idx,sep='')
    sel.pos.idx <- which(chr==tmp.chr)
    plot(loc[sel.pos.idx]/1e6,counts[sel.pos.idx],cex=plot.cex,
         xlab='position (Mb)',ylab='read counts',
         main=paste(tmp.chr,'(before correction)'))
    abline(h=ori_median_counts,lwd=3,col='red',lty=3)
    plot(loc[sel.pos.idx]/1e6,corrected_vec[sel.pos.idx],cex=plot.cex,
         xlab='position (Mb)',ylab='read counts',
         main=paste(tmp.chr,'(after correction)'))
    abline(h=ori_median_counts,lwd=3,col='red',lty=3)
  }
  
  par(mfrow=c(1,1))
  
  
  if(is.pdf){   dev.off()   }
  
  sel.nan.idx <- which(!is.na(counts))
  return(list(chr=chr[sel.nan.idx],loc=loc[sel.nan.idx],
              c.count=corrected_vec[sel.nan.idx],o.count=counts[sel.nan.idx],
              gc=gc[sel.nan.idx],median=ori_median_counts,
              sample.name=sample.name))
  
}