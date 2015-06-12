patCNV.exon.callCNV <- function(  session.info, covg.info, pattern.list,  
                                  CNV.type = c("Germline", "Somatic"),
                                  txt.output.DIR = NULL, plot.output.DIR = NULL, 
                                  feature.name.sep = " ",
                                  small.delta = 1e-2, SNR.dB.cut = 10, N.SNR.quantiles = 20)
  # output pval.mtx, CNV.mtx
{

  N.sample <- length(covg.info$total_count_vec)
  N.feature <- nrow(covg.info$exon_RPKM_mtx)

  if(is.null(txt.output.DIR)){
    txt.output.DIR <- session.info$DIR_info$txt_output_DIR
  }

  if(is.null(plot.output.DIR)){
    plot.output.DIR <- session.info$DIR_info$plot_output_DIR
  }


  sampleID.vec <- colnames(covg.info$exon_RPKM_mtx)
  featureID.vec <- paste(      session.info$exon_info$Chr,
                               session.info$exon_info$Start,
                               session.info$exon_info$Genes,
                               sep = feature.name.sep )

  CNV.mtx <- mat.or.vec(N.feature, N.sample)
  pval.mtx <- mat.or.vec(N.feature, N.sample) + 1
  colnames(CNV.mtx) <- colnames(pval.mtx) <- sampleID.vec
  rownames(CNV.mtx) <- rownames(pval.mtx) <- featureID.vec

  for( k in 1 : N.sample){
    individual.CNV.vec <-
      log2(covg.info$exon_RPKM_mtx[,k] + small.delta) -
      log2(pattern.list$median.RPKM + small.delta)

    CNV.mtx[ , k] <- individual.CNV.vec

    individual.sample.ID <- sampleID.vec[k]
    cat("processing", k, "-th sample:", individual.sample.ID, "\n")

    png(paste(plot.output.DIR, individual.sample.ID, "_CNV.png", sep = ""),
        width = 3e3, height = 1.5e3, res = 250)

    SNR.GEthreshold.idx <- which(pattern.list$SNR.dB >= SNR.dB.cut)
    patCNV.plotSimpleCNV(chr.vec = session.info$exon_info$Chr[SNR.GEthreshold.idx],
                         pos.vec = session.info$exon_info$Start[SNR.GEthreshold.idx],
                         cnv.vec = individual.CNV.vec[SNR.GEthreshold.idx],
                         main = paste("CNV plot of", individual.sample.ID),
                         ylim = c(-3, 3))
    dev.off()



    png(paste(plot.output.DIR, individual.sample.ID, "_SNR_vs_CNV.png", sep = ""),
        width = 2e3, height = 2e3, res = 300)

    CNV.conf.list <- patCNV.evaluate.CNVconf( session.info = session.info,
                                              cnv.vec = individual.CNV.vec,
                                              SNR.vec = pattern.list$SNR.dB,
                                              sample.ID = individual.sample.ID,
                                              N.SNR.quantiles = N.SNR.quantiles)

    dev.off()

    singleCNV.txt.filename <- paste(txt.output.DIR, individual.sample.ID,
                                    "_CNV.txt", sep = "")

    if( CNV.type == "Germline" ){
      pval.mtx[ , k] <- CNV.conf.list$pval
      txt.out.mtx <- cbind( chr = session.info$exon_info$Chr,
                            start.pos = session.info$exon_info$Start,
                            stop.pos = session.info$exon_info$Stop,
                            gene = session.info$exon_info$Genes,
                            CNV.conf.list )
    }  else { # somatic
      txt.out.mtx <- cbind( chr = session.info$exon_info$Chr,
                            start.pos = session.info$exon_info$Start,
                            stop.pos = session.info$exon_info$Stop,
                            gene = session.info$exon_info$Genes,
                            CNV.log2ratio = CNV.conf.list$CNV.log2ratio,
                            SNR.dB = CNV.conf.list$SNR.db )
    }
    write.table(x = txt.out.mtx, file = singleCNV.txt.filename,
                quote = FALSE, row.names = FALSE, sep = "\t")

  } #  for( k in 1 : N.sample)

  CNVmtx.txt.filename <- paste(txt.output.DIR, CNV.type,
                               "_CNV_matrix.txt", sep = "")
  Pvalmtx.txt.filename <- paste(txt.output.DIR, CNV.type,
                                "_pval_matrix.txt", sep = "")
  if( CNV.type == "Germline" ){
    write.table(x = CNV.mtx, file = CNVmtx.txt.filename,
                quote = FALSE, row.names = TRUE, sep = "\t")
    write.table(x = pval.mtx, file = Pvalmtx.txt.filename,
                quote = FALSE, row.names = TRUE, sep = "\t")

    return(list(CNV.mtx = CNV.mtx,
                pval.mtx = pval.mtx))

  } else { # somatic
    write.table(x = CNV.mtx, file = CNVmtx.txt.filename,
                quote = FALSE, row.names = TRUE, sep = "\t")
    return(list(CNV.mtx = CNV.mtx))
  }

}
