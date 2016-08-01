patCNV.sample.QC <- function(session.info, covg.info,
                             output.prefix = c("Germline", "Somatic"),
                             outlier.min.median.corval = 0.85,
                             heatmap.text.cex = 0.35, rpkm.density.bw = 0.1,
                             txt.output.DIR = NULL, plot.output.DIR = NULL)
{
  N.sample <- length(covg.info$total_count_vec)
  sampleID.vec <- colnames(covg.info$exon_RPKM_mtx)
  
  if(is.null(txt.output.DIR)){
    txt.output.DIR <- session.info$DIR_info$txt_output_DIR
  }  
  
  if(is.null(plot.output.DIR)){
    plot.output.DIR <- session.info$DIR_info$plot_output_DIR
  }  
  
  
  sample.cormtx <- cor(covg.info$exon_RPKM_mtx)
  sample.median.corval <- apply(sample.cormtx,1,median, na.rm = TRUE)
  
  outlier.sampleSet <- names(which(sample.median.corval < outlier.min.median.corval))
  
  if(length(outlier.sampleSet)>0){
    cat("patternCNV warning: likely outliers:\n",
        paste(outlier.sampleSet,"\n"),"\n")  }
  
  sample.outlier.idx <- match(outlier.sampleSet, sampleID.vec)
  barplot.color <- rep("grey",length(sample.median.corval))
  barplot.color[sample.outlier.idx] <- "red"
  
  is.outlier <- is.element(sampleID.vec, outlier.sampleSet )
  mch.median.corval <- sample.median.corval[sampleID.vec] 
  
  sample.QC.mtx <- cbind(median.corval = mch.median.corval, is.outlier, 
                         total.million.BP.cvg = round(covg.info$total_count_vec/1e6,digits = 1))
  rownames(sample.QC.mtx) <- sampleID.vec
  
  
  approx.x <- seq(-5, 12, 0.5)
  approx.density.mtx <- mat.or.vec( length(approx.x), N.sample)
  rownames(approx.density.mtx) <- round(approx.x,digits = 2)
  colnames(approx.density.mtx) <- sampleID.vec
  
  pdf(paste(plot.output.DIR, output.prefix, "_QC_plot.pdf",sep = ""))
  
  
  hm <- heatmap.2(sample.cormtx,trace="none",col=bluered(25), 
            cexRow = heatmap.text.cex, cexCol = heatmap.text.cex,
            labRow = sampleID.vec, 
            labCol = sampleID.vec, 
            main = "pair-wise sample correlations")
  
  order.sample.median.corval.idx <- order(sample.median.corval, decreasing = TRUE)
  barplot(sample.median.corval[order.sample.median.corval.idx], 
          names.arg = sampleID.vec[order.sample.median.corval.idx],
          col = barplot.color[order.sample.median.corval.idx],
          xlab="median sample correlation", horiz = TRUE,las=1, cex.names=0.35,
          main = "median of pair-wise correlation (outlier in red)")
  
  par(mar=c(5.1, 8 ,4.1 ,2.1))
  barplot(covg.info$total_count_vec/1e6, names.arg = sampleID.vec,
          xlab="total bp covg (Million)", horiz = TRUE,las=1, cex.names=0.35,
          col=barplot.color, main = "total coverage (outlier in red)")
  par(mar=c(5.1, 4.1 ,4.1 ,2.1))
  
  
  plot( density(log2(covg.info$exon_RPKM_mtx[,1]+0.1), bw = rpkm.density.bw),
        xlim=c(-5, 12),xlab="log2(RPKM)", 
        main=paste("log2(RPKM) densities of",N.sample,"samples"), 
        col="white", ylim=c(0,0.39), lwd = 2 )
  for(j in 1:N.sample){
    tmp.density <- density(log2(covg.info$exon_RPKM_mtx[,j]+0.1), 
                           bw = rpkm.density.bw)
    
    lines( tmp.density, col=rainbow(N.sample)[j], lwd = 2)
    
    approx.density <- approxfun(tmp.density)
    approx.density.mtx[,j] <- approx.density(approx.x)
    approx.density.mtx[which(is.na(approx.density.mtx[,j])),j] <- 0
    
    lines( approx.x, approx.density(approx.x),
           col = rainbow(N.sample)[j], lwd = 1.5, pch = j, type = "p")
  }
  legend("topright", sampleID.vec, col = rainbow(N.sample), 
         lwd = 1, pch = 1 : N.sample, cex = 0.6)  
  
  heatmap.2(t(approx.density.mtx), main = "heatmap of log2(RPKM) densities",
            trace="none", Colv = FALSE, col = bluered(25), 
            cexCol = 0.7, cexRow = 0.7)
  
  dev.off()
  
  write.table(x = as.data.frame(sample.QC.mtx),
              file = paste(plot.output.DIR, output.prefix, "_QC_table.txt",sep = ""), 
              quote = FALSE, sep = "\t", row.names = TRUE)

  
  rownames(sample.cormtx) <- colnames(sample.cormtx) <- sampleID.vec
  write.table(sample.cormtx[rev(hm$rowInd), hm$colInd],
              file = paste(plot.output.DIR, output.prefix, "_heatmap_cor_matrix.txt",sep = ""),
              quote = FALSE, sep = "\t", row.names = TRUE)

  return(as.data.frame(sample.QC.mtx))
}
