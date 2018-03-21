patCNV.plotSimpleCNV <- function( chr.vec, pos.vec, cnv.vec, 
                                  text.cex = 0.8, cnv.cex = 0.6,
                                  ylim = NULL, ylab = "CNV log2-ratio",
                                  chr.color =  rep(c("magenta","dodgerblue"),12),
                                  ... )
{
  cnv.color.idx <- chr.vec
  cnv.color.idx[chr.vec == "chrX" | chr.vec == "X"] <- "chr23"
  cnv.color.idx[chr.vec == "chrY" | chr.vec == "Y"] <- "chr24"
  cnv.color.idx <- as.numeric(gsub( "chr", "", cnv.color.idx))
  
  #order.idx <- order(cnv.color.idx*1e9 + pos.vec*0.05, decreasing = FALSE) # order according to chr & pos
  order.idx <- order(cnv.color.idx, pos.vec, decreasing = FALSE) # order according to chr & pos
  
  chr.Str <- c( seq(1,22), "X", "Y")
  
  ordered.cnv.vec <- cnv.vec[order.idx]
  ordered.chr.vec <- cnv.color.idx[order.idx]
  
  if(is.null(ylim)){ ylim <- range(ordered.cnv.vec) }
  outliersLess<-which(ordered.cnv.vec< ylim[1])
  outliersMore<-which(ordered.cnv.vec> ylim[2])
  

  plot(ordered.cnv.vec, 
       col = chr.color[ordered.chr.vec], ylim = ylim, ylab = ylab,
       cex = cnv.cex, ...)
  abline(h = 0, lty = 2, lwd = 2, col = "black")
  points(seq(1,length(ordered.cnv.vec))[outliersLess],rep(ylim[1],length(outliersLess)),col = "blue", pch=6)
  points(seq(1,length(ordered.cnv.vec))[outliersMore],rep(ylim[2],length(outliersMore)),col = "red", pch=2)


  for(chr.idx in 1 : 24){
    tmp.chr.Set.idx <- which(ordered.chr.vec == chr.idx )
    tail.pos <- tail(tmp.chr.Set.idx, 1)
    if(chr.idx < 24)
    { abline(v = tail.pos, lty = 2, lwd = 2, col = "grey") } # no end line for chrY
    tmp.middle.pos <- range(tmp.chr.Set.idx)[1] + diff(range(tmp.chr.Set.idx))/2
    
    text( x = tmp.middle.pos, y = ylim[1]*(0.85 + 0.1*(chr.idx%%2)), cex = text.cex,
          labels = chr.Str[chr.idx])
    text( x = tmp.middle.pos, y = ylim[2]*(0.85 + 0.1*(chr.idx%%2)), cex = text.cex,
          labels = chr.Str[chr.idx])
    
  }
  
}
