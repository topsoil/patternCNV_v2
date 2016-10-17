patCNV.exon.callCNV <- function(  session.info, covg.info, pattern.list,  
                                  CNV.type = c("Germline", "Somatic"),
                                  txt.output.DIR = NULL, plot.output.DIR = NULL, 
                                  feature.name.sep = " ",
                                  small.delta = 1e-2, SNR.dB.cut = 10, N.SNR.quantiles = 20,predicted.sex=NULL,known.sex=NULL)
  # output pval.mtx, CNV.mtx
  # output exon_count_mtx, exon_RPKM_mtx
  # return(list(total_count_vec=total_count_vec,
  #          exon_count_mtx=exon_count_mtx,exon_RPKM_mtx=exon_RPKM_mtx))
{


 sampleID.vec <- colnames(covg.info$exon_RPKM_mtx)
#  colnames.vec<-make.names(session.info$file_info$ID,unique=T)
  sampleids<-session.info$file_info$ID[match(sampleID.vec,session.info$file_info$ID)]


  N.sample <- length(covg.info$total_count_vec)
  N.feature <- nrow(covg.info$exon_RPKM_mtx)

  diploidExons<-which((session.info$exon_info$Chr!="chrX" & session.info$exon_info$Chr!="chrY"  ) | session.info$exon_info$PAR==1)
  haploidXExons<-which((session.info$exon_info$Chr=="chrX" ) & session.info$exon_info$PAR==0)

  if(is.null(txt.output.DIR)){
    txt.output.DIR <- session.info$DIR_info$txt_output_DIR
  }

  if(is.null(plot.output.DIR)){
    plot.output.DIR <- session.info$DIR_info$plot_output_DIR
  }



  featureID.vec <- paste(      session.info$exon_info$Chr,
                               session.info$exon_info$Start,
                               session.info$exon_info$Genes,
                               sep = feature.name.sep )

  CNV.mtx <- mat.or.vec(N.feature, N.sample)
  pval.mtx <- mat.or.vec(N.feature, N.sample) + 1
#  colnames(CNV.mtx) <- colnames(pval.mtx) <- sampleID.vec
  colnames(CNV.mtx) <- colnames(pval.mtx) <- sampleids
  rownames(CNV.mtx) <- rownames(pval.mtx) <- featureID.vec

  exon.raw.count.mtx <- covg.info$exon_count_mtx
  exon.nmlz.RPKM.mtx <- covg.info$exon_RPKM_mtx
  rownames(exon.raw.count.mtx) <- rownames(exon.nmlz.RPKM.mtx) <- featureID.vec
  
  for( k in 1 : N.sample){
    individual.CNV.vec <-
      log2(covg.info$exon_RPKM_mtx[,k] + small.delta) -
      log2(pattern.list$median.RPKM + small.delta)

#      print(paste("Sample ",sampleID.vec[k], ", meanX=",mean(covg.info$exon_RPKM_mtx[haploidXExons,k]),
#      ", mean X model=",mean(pattern.list$median.RPKM[haploidXExons]),
#      "mean diploid model=", mean(pattern.list$median.RPKM[diploidExons]),", mean diploid sample=",mean(covg.info$exon_RPKM_mtx[diploidExons,k]),"\n"))

    CNV.mtx[ , k] <- individual.CNV.vec

    individual.sample.ID <- sampleids[k]
    cat("processing", k, "-th sample:", individual.sample.ID, "\n")

    png(paste(plot.output.DIR, individual.sample.ID, "_CNV.png", sep = ""),
        width = 3e3, height = 1.5e3, res = 250)

    SNR.GEthreshold.idx <- which(pattern.list$SNR.dB >= SNR.dB.cut)

    sexString<-getSexTitleString(predicted.sex,known.sex,k)

    patCNV.plotSimpleCNV(chr.vec = session.info$exon_info$Chr[SNR.GEthreshold.idx],
                         pos.vec = session.info$exon_info$Start[SNR.GEthreshold.idx],
                         cnv.vec = individual.CNV.vec[SNR.GEthreshold.idx],
                         main = paste("CNV plot of", individual.sample.ID,"\n",sexString),
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
  RawCount.mtx.txt.filename <- paste(txt.output.DIR, CNV.type,
                               "_raw_count_matrix.txt", sep = "")
  RPKMmtx.txt.filename <- paste(txt.output.DIR, CNV.type,
                                "_RPKM_matrix.txt", sep = "")
  if( CNV.type == "Germline" ){
    write.table(x = CNV.mtx, file = CNVmtx.txt.filename,
                quote = FALSE, row.names = TRUE, sep = "\t")
    write.table(x = pval.mtx, file = Pvalmtx.txt.filename,
                quote = FALSE, row.names = TRUE, sep = "\t")
	
	write.table(x = exon.raw.count.mtx, file = RawCount.mtx.txt.filename,
                quote = FALSE, row.names = TRUE, sep = "\t")
				
    write.table(x = exon.nmlz.RPKM.mtx, file = RPKMmtx.txt.filename,
                quote = FALSE, row.names = TRUE, sep = "\t")
	
  
    return(list(CNV.mtx = CNV.mtx,
                pval.mtx = pval.mtx))

  } else { # somatic
    write.table(x = CNV.mtx, file = CNVmtx.txt.filename,
                quote = FALSE, row.names = TRUE, sep = "\t")
				
	write.table(x = exon.raw.count.mtx, file = RawCount.mtx.txt.filename,
                quote = FALSE, row.names = TRUE, sep = "\t")
				
    write.table(x = exon.nmlz.RPKM.mtx, file = RPKMmtx.txt.filename,
                quote = FALSE, row.names = TRUE, sep = "\t")
	
    return(list(CNV.mtx = CNV.mtx))
  }

}


