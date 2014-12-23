patCNV.exon.segmentCNV <- function( session.info, CNV.mtx, pattern.list,  
                                    SNR.dB.cut = 10, is.plot = TRUE,
                                    plot.ylim=c(-3,3), plot.cex=0.6,
                                    txt.output.DIR = NULL, plot.output.DIR = NULL,
                                    ...)
  # ... is for segment() in DNAcopy package  
  
{
  
  require(DNAcopy)
  
  N.sample <- ncol(CNV.mtx)
  N.feature <- nrow(CNV.mtx)
  
  if(is.null(txt.output.DIR)){
    txt.output.DIR <- session.info$DIR_info$txt_output_DIR
  }  
  
  if(is.null(plot.output.DIR)){
    plot.output.DIR <- session.info$DIR_info$plot_output_DIR
  }  
  
  sampleID.vec <- colnames(CNV.mtx)
  
  output_suffix <- "_CNV_seg_"
  sel_idx <-  which(pattern.list$SNR.dB >= SNR.dB.cut)
  
  chrom <- session.info$exon_info$Chr[sel_idx]
  maploc <- session.info$exon_info$Start[sel_idx]
  maploc.end <- session.info$exon_info$Stop[sel_idx]
  
  for(smpl_idx in 1 : N.sample)
  {
    
    sel_sample_ID <- sampleID.vec[smpl_idx]
    
    cat("processing", smpl_idx, "-th sample:", sel_sample_ID, "\n")
    
    
    tmp_output_bed_file <- 
      paste(txt.output.DIR, sel_sample_ID, output_suffix,'.bed',sep='')
    tmp_output_txt_file <- 
      paste(txt.output.DIR, sel_sample_ID, output_suffix,'.txt',sep='')
    
    cna.res_My <- CNA(CNV.mtx[sel_idx, smpl_idx], chrom, maploc,
                      data.type=c("logratio"), sampleid=sel_sample_ID)
    seg.res_My <- segment(cna.res_My, verbose=0, ...)
    
    seg_mtx <- seg.res_My$output
    bed_mtx <- cbind(as.character(seg_mtx$chrom),seg_mtx$loc.start,
                     maploc.end[seg.res_My$segRows$endRow],seg_mtx$seg.mean) # end of segment is modified according to end of exon
    write.table(x=bed_mtx,file=tmp_output_bed_file,quote=FALSE,
                row.names=FALSE,col.names=FALSE,sep='\t')
    
    txt_mtx <- cbind(as.character(seg_mtx$chrom),seg_mtx$loc.start,
                     maploc.end[seg.res_My$segRows$endRow],seg_mtx$seg.mean,seg_mtx$num.mark) # end of segment is modified according to end of exon
    write.table(x=txt_mtx,file=tmp_output_txt_file,quote=FALSE,
                row.names=FALSE,col.names=FALSE,sep='\t')
    #the only difference between .bed and .txt files is "seg_mtx$num.mark"
    
    if (is.plot)
    {
      
      CNV_seg_pdf_filename <- 
        paste(plot.output.DIR, sel_sample_ID,output_suffix,'.pdf',sep='')
      pdf(CNV_seg_pdf_filename)
      plot(seg.res_My, ylim=plot.ylim, cex=plot.cex)
      plot(seg.res_My, ylim=plot.ylim,
           cex=plot.cex, plot.type="chrombysample",xmaploc=TRUE)
      #xmaploc=TRUE x-axis should be according to chr position rather than bin number
      dev.off()  
    } # if (is.plot)
    
  } # end of for(smpl_idx
}