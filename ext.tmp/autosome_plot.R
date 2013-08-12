autosome_plot <- function(res_list, plot_DIR='plot_output',
                          odd.col='orange', even.col='skyblue',plot.cex=0.6,
                          is.chr.labeled=TRUE,
                          odd.y.pos=-10,even.y.pos=-38,label.cex=0.85)
  
{
  chr <- res_list$chr
  gc <- res_list$gc
  counts <- res_list$c.count
  
  color_vec <- as.character(chr)
  odd.chr <- paste('chr',seq(1,22,2),sep='')
  even.chr <- paste('chr',seq(2,22,2),sep='')
  color_vec[which(is.element(color_vec,odd.chr))] <- odd.col
  color_vec[which(is.element(color_vec,even.chr))] <- even.col 
  
  png_filename <- paste(plot_DIR,'/',res_list$sample.name,'_autosome.png',sep='')
  png(png_filename,width=3e3,height=1e3,res=180)
  
  up.y <- quantile(gc.res$c.count,0.995,na.rm=TRUE)*2.5
  plot(counts,cex=plot.cex,
       col=color_vec, ylim = c(-50,up.y),
       xlab='Bin index',ylab='Coverage (reads)',
       main=paste('coverage of',res_list$sample.name))
  abline(h=res_list$median,lwd=3,col='red',lty=3)
  
  if(is.chr.labeled)
  {
    for(chr.idx in 1:22)
    {
      tmp.chr <- paste('chr',chr.idx,sep='')
      sel.pos.idx <- which(chr==tmp.chr)
      median_pos <- median(sel.pos.idx)
      if(is.element(tmp.chr,odd.chr))
      { tmp.col <- odd.col 
        y.tmp.pos <- odd.y.pos
      } else {
        tmp.col <- even.col
        y.tmp.pos <- even.y.pos
      }
      text(x=median_pos,y=y.tmp.pos,labels=chr.idx,col=tmp.col,cex=label.cex)
    }
  }
  
  dev.off()
  
}