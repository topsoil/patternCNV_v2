patCNV.learnExonPattern <- function( session.info, covg.info,
                                     sample.QC.list,
                                     check.gender = TRUE, sex.minimum.RPKM = 0.1,
                                     small.delta = 1e-2)
{

  sample.outlier.idx <-
    match( rownames(sample.QC.list)[which(sample.QC.list$is.outlier == 1)],
           colnames(covg.info$exon_RPKM_mtx) )

  if(check.gender){
    sex.chrX.idx <- which(session.info$exon_info$Chr=="chrX" )
    sex.chrY.idx <- which(session.info$exon_info$Chr=="chrY" )


    chrX.mean.RPKM <- apply( covg.info$exon_RPKM_mtx[sex.chrX.idx,], 1, mean)
    chrY.mean.RPKM <- apply( covg.info$exon_RPKM_mtx[sex.chrY.idx,], 1, mean)

    reliable.sex.chrX.idx <-  sex.chrX.idx[which(chrX.mean.RPKM>=sex.minimum.RPKM)]
    reliable.sex.chrY.idx <-  sex.chrY.idx[which(chrY.mean.RPKM>=sex.minimum.RPKM)]

    sample.chrX.meanCvg <- apply(  (covg.info$exon_RPKM_mtx[reliable.sex.chrX.idx,]), 2, mean)
    sample.chrY.meanCvg <- apply(  (covg.info$exon_RPKM_mtx[reliable.sex.chrY.idx,]), 2, mean)


    YvsX.ratio.cut <- 0.45

    sample.YvsX.ratio <- (sample.chrY.meanCvg/sample.chrX.meanCvg)
    sample.IsMale <- (sample.YvsX.ratio > YvsX.ratio.cut)
    male.sample.idx <- which(sample.IsMale)
    female.sample.idx <- which(!sample.IsMale)

    chrX.RPKM.mtx <- as.matrix(covg.info$exon_RPKM_mtx[sex.chrX.idx,])
    if(!is.null(male.sample.idx)){
      chrX.RPKM.mtx[,male.sample.idx] <-  2 * chrX.RPKM.mtx[,male.sample.idx]
    }
    chrX.median.vec <- apply( (chrX.RPKM.mtx ), 1, median)
    chrX.MAD.vec <- apply( (chrX.RPKM.mtx ), 1, mad)


    if(is.null(male.sample.idx)){
      chrY.median.vec <- NULL
      chrY.MAD.vec <- NULL
    } else {
      chrY.RPKM.mtx <- 2 * as.matrix(covg.info$exon_RPKM_mtx[sex.chrY.idx, male.sample.idx])
      chrY.median.vec <- apply(  chrY.RPKM.mtx , 1, median)
      chrY.MAD.vec <- apply(  chrY.RPKM.mtx , 1, mad )
    }


  } # check gender and adjust exon patterns by genders

  if(length(sample.outlier.idx) == 0) {
    median.vec <- apply((covg.info$exon_RPKM_mtx ), 1, median)
    MAD.vec <- apply((covg.info$exon_RPKM_mtx ), 1, mad)
  } else {
    median.vec <- apply((covg.info$exon_RPKM_mtx[,-sample.outlier.idx] ), 1, median)
    MAD.vec <- apply((covg.info$exon_RPKM_mtx[,-sample.outlier.idx] ), 1, mad)
  }



  median.vec[sex.chrX.idx] <- chrX.median.vec
  median.vec[sex.chrY.idx] <- chrY.median.vec

  MAD.vec[sex.chrX.idx] <- chrX.MAD.vec
  MAD.vec[sex.chrY.idx] <- chrY.MAD.vec

  SNR.ratio.vec <- (median.vec) / (1.4826 * (MAD.vec) + small.delta )
  if(length(which(is.na(SNR.ratio.vec)))>0) {	# fill SNR == NA with min.SNR 
	  min.SNR.ratio.val <- min(SNR.ratio.vec, na.rm=TRUE)
	  SNR.ratio.vec[which(is.na(SNR.ratio.vec))] <- min.SNR.ratio.val
  }
  
  SNR.dB.vec <- 20 * log10( SNR.ratio.vec )

  mean.count.vec <- apply(covg.info$exon_count_mtx,1,mean, na.rm=TRUE)

  patternMtx <- cbind( median.RPKM = median.vec,
                       MAD.RPKM = MAD.vec,
                       mean.raw.BPcount = mean.count.vec,
                       SNR.ratio = SNR.ratio.vec,
                       SNR.dB = SNR.dB.vec)

  return( as.data.frame(patternMtx) )
}
