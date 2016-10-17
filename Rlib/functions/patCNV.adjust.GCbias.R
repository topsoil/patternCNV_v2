#===  adjust GC-coverage bias from original RPKM coverage

patCNV.adjust.GCbias <- function(session.info, covg.info,
                                 output.prefix = c("Germline", "Somatic"),
                                 sspline.df = 3, GC.slope.range = c(0.3, 0.7),
                                 GC.train.min.log2RPKM = 1.5,
                                 small.delta = 1e-2,
                                 GC.hist.break = 50,
                                 GC.splinePlot.ylim = c(2, 3.8),
                                 txt.output.DIR = NULL, plot.output.DIR = NULL)
{
  
  N.sample <- length(covg.info$total_count_vec)
  sampleID.vec <- colnames(covg.info$exon_RPKM_mtx)
#  colnames.vec<-make.names(session.info$file_info$ID,unique=T)
  sampleids<-session.info$file_info$ID[match(sampleID.vec,session.info$file_info$ID)]
  
  if(is.null(txt.output.DIR)){
    txt.output.DIR <- session.info$DIR_info$txt_output_DIR
  }  
  
  if(is.null(plot.output.DIR)){
    plot.output.DIR <- session.info$DIR_info$plot_output_DIR
  }  
  
  
  GC.slope.vec <- mat.or.vec(N.sample, 1)
  GC.maxmin.Difflog2R.vec <- mat.or.vec(N.sample, 1)
  
  names(GC.slope.vec) <- names(GC.maxmin.Difflog2R.vec) <- sampleID.vec
  
  GC.vec <- (session.info$exon_info$GC_Content)
  
  LR.sdelta <- small.delta
  log2RPKM.mtx <- log2(covg.info$exon_RPKM_mtx + LR.sdelta)
  log2.median.vec <- apply(log2RPKM.mtx, 1, median, na.rm = TRUE)
  train.exon.idx <- which((log2.median.vec >= GC.train.min.log2RPKM) & ((session.info$exon_info$Chr!="chrX" & session.info$exon_info$Chr!="chrY") | session.info$exon_info$PAR==1))
  print(paste("training GC exons=",length(train.exon.idx),"\n",sep=""))
  crct.log2RPKM.mtx <- log2RPKM.mtx
  for(k in 1:N.sample){
    tmp.GC.train <- GC.vec[train.exon.idx]
    tmp.Cvg.train <- log2RPKM.mtx[train.exon.idx,k]
    tmp.median.train <- median(tmp.Cvg.train)
    
    sspline.res <- smooth.spline( tmp.GC.train, tmp.Cvg.train, df = sspline.df)
    
    GC.slope.vec[k] <- 
      diff(predict( sspline.res, GC.slope.range )$y) / diff(GC.slope.range)
    
    tmp.GC.cvg.range <- (range(predict( sspline.res, seq(GC.slope.range[1], GC.slope.range[2], length.out = 100) )$y))
    GC.maxmin.Difflog2R.vec[k] <- tmp.GC.cvg.range[2] - tmp.GC.cvg.range[1]

    tmp.GC <- GC.vec
    tmp.Cvg <- log2RPKM.mtx[,k]

    GC.predict.Cvg <- predict(sspline.res, tmp.GC)$y
    crct.log2RPKM.mtx[, k ] <- ( tmp.Cvg / GC.predict.Cvg ) * tmp.median.train
    print(paste(sampleids[k]," baseline shift= ",mean(crct.log2RPKM.mtx[train.exon.idx, k ]-tmp.Cvg.train),", slope=",GC.slope.vec[k],sep=""))
  }
  
  crct.RPKM.mtx <- 2^crct.log2RPKM.mtx
  
  
  
  
  #pdf("GC_9samples_per_page.pdf")
  pdf(paste(plot.output.DIR, output.prefix, "_GCcovg_summary_plot.pdf",sep = ""))
  
  hist(GC.vec,GC.hist.break, xlab = "GC content", main = "GC distribution")
  
  par(mfrow = c(3,3))
  #par(mar=c(1.5, 1.5 ,4.5 ,1.5))
  for(j in 1:N.sample){
    
    plot(smooth.spline( GC.vec[train.exon.idx], 
                        log2(covg.info$exon_RPKM_mtx[train.exon.idx,j]+LR.sdelta),
                        df = sspline.df),
#         col = "black", lwd = 2, ylim = GC.splinePlot.ylim, type = "l", xlab = "", ylab = "", 
         col = "black", lwd = 2, type = "l", xlab = "", ylab = "", 
         main=paste("",
                    sampleids[j],"\n",
                    "\n GC slope =", round(GC.slope.vec[j],digits = 2)
         ))
    
  }
  par(mfrow = c(1,1))
  
  for(j in 1 : N.sample){
    smoothScatter( GC.vec[train.exon.idx], 
                   log2(covg.info$exon_RPKM_mtx[train.exon.idx,j]+LR.sdelta), 
                   xlab="GC content", ylab="log2(RPKM)",
                   main=paste("",
                              sampleids[j],"\n",
                              "\n GC slope =", round(GC.slope.vec[j],digits = 2) ) )
    lines(smooth.spline( GC.vec[train.exon.idx], 
                         log2(covg.info$exon_RPKM_mtx[train.exon.idx,j]+LR.sdelta),df = sspline.df),
          col="red", lwd=2)
  }
  
  dev.off()
  
  covg.info$exon_RPKM_mtx <- crct.RPKM.mtx
  covg.info$GC.slope <- GC.slope.vec
  
  write.table( x = as.data.frame(
    list(BP.covg.In.million = covg.info$total_count_vec/1e6,
         GC.covg.bias.slope = GC.slope.vec)
  ),
  file = paste(plot.output.DIR, output.prefix, "_GC_Covg_table.txt",sep = ""), 
  quote = FALSE, sep = "\t", row.names = TRUE)
  
  return( covg.info ) # GC.adjust
}
