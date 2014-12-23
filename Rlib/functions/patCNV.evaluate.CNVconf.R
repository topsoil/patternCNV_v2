patCNV.evaluate.CNVconf <- function ( session.info, cnv.vec, SNR.vec, 
                                      sample.ID = NULL, N.SNR.quantiles = 20)
  #SNR.vec <- pattern.res$SNR.dB
{
  SNR.GE0_Autosome.idx <- which(SNR.vec >= 0 & 
                                  !is.element(session.info$exon_info$Chr, 
                                              c("chrX", "chrY", "chrM", "X", "Y", "M")))
  

  
  
  
  
  SNR.CNV.upperMAD <- mat.or.vec(N.SNR.quantiles, 1)
  SNR.CNV.median <- mat.or.vec(N.SNR.quantiles, 1)
  SNR.median <- mat.or.vec(N.SNR.quantiles, 1)
  SNR.CNV.lowerMAD <- mat.or.vec(N.SNR.quantiles, 1)
  
  # stop here
  tmp.SNR.vec <- SNR.vec[SNR.GE0_Autosome.idx] 
  tmp.SNR.vec <- tmp.SNR.vec + (1e-6)*runif(n=length(tmp.SNR.vec))
  tmp.CNV.vec <- cnv.vec[SNR.GE0_Autosome.idx]
  tmp.CNV.vec <- tmp.CNV.vec + (1e-6)*rnorm(n=length(tmp.CNV.vec))
  
  SNR.quantile.interval <- quantile(tmp.SNR.vec, 
                                    seq(0, 1, length.out = N.SNR.quantiles + 1))
  
  feature.SNR.interval.idx <- findInterval( tmp.SNR.vec, 
                                            SNR.quantile.interval, rightmost.closed = TRUE)
  
  
  for(q in 1 : N.SNR.quantiles){
    sel.element.idx <- which(feature.SNR.interval.idx == q)
    
    SNR.median[q]<- median( tmp.SNR.vec[sel.element.idx] )
    SNR.CNV.median[q] <- median( tmp.CNV.vec[sel.element.idx] )
    
    upper.idx <- intersect( which(tmp.CNV.vec > SNR.CNV.median[q]), 
                            sel.element.idx)
    SNR.CNV.upperMAD[q] <- median(tmp.CNV.vec[upper.idx] - SNR.CNV.median[q])
    
    lower.idx <- intersect( which(tmp.CNV.vec < SNR.CNV.median[q]), 
                            sel.element.idx)
    SNR.CNV.lowerMAD[q] <- median( SNR.CNV.median[q] - tmp.CNV.vec[lower.idx] )
  }
  
  sspline.medianCNV <- 
    smooth.spline( SNR.median,  (SNR.CNV.median ))
  sspline.upperCNV.SD <- 
    smooth.spline( SNR.median,  (1.4826 *SNR.CNV.upperMAD) )
  sspline.lowerCNV.SD <- 
    smooth.spline( SNR.median,  (1.4826 *SNR.CNV.lowerMAD) )
  
  
  
  lwd.cex <- 2.5
  point.cex <- 0.6
  point.color <- "darkgrey"
  smoothScatter( tmp.SNR.vec, 
                 2*2^tmp.CNV.vec, 
                 xlim = c( 2, max(tmp.SNR.vec)),
                 ylim = c(0,5), cex = point.cex, col = point.color,
                 ylab = "copy number estimates",
                 xlab = "Signal to Noise Ratio (SNR)\n 20*log10(avg.coverage/SD.coverage)",
                 main = paste("SNR vs. copy number estimate plot \n", sample.ID ) )
  abline( h = 1 , lty = 2, lwd = lwd.cex, col = "green")
  abline( h = 2 , lty = 2, lwd = lwd.cex, col = "black")
  abline( h = 3 , lty = 2, lwd = lwd.cex, col = "orange")
  abline( h = 4 , lty = 2, lwd = lwd.cex, col = "red")
  
  tmp.minSNR <- quantile(tmp.SNR.vec, 2e-3)
  tmp.maxSNR <- quantile(tmp.SNR.vec, 1 - 2e-3)
  tmp.axis.SNR <- seq(tmp.minSNR, tmp.maxSNR, length.out = 100)
  #lines(tmp.axis.SNR, 2*2^predict(sspline.medianCNV, tmp.axis.SNR)$y, col = "black")
  
  tmp.axis.CNVupper <- predict(sspline.medianCNV, tmp.axis.SNR)$y + 
    2*predict(sspline.upperCNV.SD, tmp.axis.SNR)$y
  
  tmp.axis.CNVlower <- predict(sspline.medianCNV, tmp.axis.SNR)$y - 
    2*predict(sspline.lowerCNV.SD, tmp.axis.SNR)$y
  
  lines(tmp.axis.SNR, 2*2^(tmp.axis.CNVupper), col = "blue")
  lines(tmp.axis.SNR, 2*2^(tmp.axis.CNVlower), col = "blue")
  legend("topright", "95% confidence interval", col = "blue", lwd = 2)
  
  
  #==== compute confidence interval 
  small.mad.offset <- 0.1
  tmp.upper.zscore <- 
    ( cnv.vec - predict(sspline.medianCNV, SNR.vec)$y ) / 
    (predict(sspline.upperCNV.SD, SNR.vec)$y + small.mad.offset)
  
  tmp.lower.zscore <- 
    ( cnv.vec - predict(sspline.medianCNV, SNR.vec)$y ) / 
    (predict(sspline.lowerCNV.SD, SNR.vec)$y + small.mad.offset)
  
  #cor(tmp.upper.zscore, tmp.lower.zscore)
  
  tmp.zscore <- tmp.upper.zscore
  tmp.NegZ.idx <- which(tmp.upper.zscore < 0 )
  tmp.zscore[ tmp.NegZ.idx ] <- tmp.lower.zscore[ tmp.NegZ.idx ]
  
  tmp.1tail.pval <- 1 - pnorm(tmp.zscore)
  pval.vec <- pmin(tmp.1tail.pval, 1-tmp.1tail.pval)*2
  
  #hist(tmp.1tail.pval)
  #hist(tmp.pval)
  
  upper95.cnv.vec <- rep( NA, length(cnv.vec) )
  upper95.cnv.vec <- cnv.vec + 
    2 * predict(sspline.upperCNV.SD, SNR.vec)$y
  
  lower95.cnv.vec <- rep( NA, length(cnv.vec) )
  lower95.cnv.vec <- cnv.vec - 
    2 * predict(sspline.lowerCNV.SD, SNR.vec)$y
  
  CNV.conf.mtx <- 
    cbind( CNV.log2ratio = cnv.vec, 
           CNV.95upper = upper95.cnv.vec, CNV.95lower = lower95.cnv.vec,
           pval = pval.vec, SNR.db = SNR.vec )
  
  
  return(as.data.frame(CNV.conf.mtx))
  
}