patCNV.segment.CNV <- function( cnv_res, session_info, 
				capture.only=TRUE,ref_avg_type='median',min_ref_avg_RPKM=3,
                                is.plot=TRUE, 
				plot.ylim=c(-3,3),plot.cex=0.6,output_suffix='CNV_seg',...)
# ... is for segment() in DNAcopy package  
 
{
  
  require(DNAcopy)
  patCNV.create.DIR(session_info$DIR$txt_output_DIR)
  
  if(ref_avg_type=='median')	
     {ref_avg_RPKM <- unlist(read.delim(session_info$Misc$median_RPKM_file,header=FALSE))}

  if(ref_avg_type=='mean')	
     {ref_avg_RPKM <- unlist(read.delim(session_info$Misc$mean_RPKM_file,header=FALSE))} 
	
  
  sample_ID_vec <- cnv_res$sample.ID
  N_sample <- length(sample_ID_vec)
  cnv.mtx <- as.matrix(cnv_res$CNV)  

  if(capture.only==TRUE){ 
                          sel_idx <- which(ref_avg_RPKM>=min_ref_avg_RPKM & 
                          session_info$exon_info$is_capture_vec==1) } else {
                          sel_idx <- which(ref_avg_RPKM>=min_ref_avg_RPKM)
                        }
  sel_idx <- setdiff(sel_idx,
  which(session_info$exon_info$Chr=='chrX' | session_info$exon_info$Chr=='chrY')) # exclude sex Chr
  
  chrom <- session_info$exon_info$Chr[sel_idx]
  maploc <- session_info$exon_info$Start[sel_idx]
  maploc.end <- session_info$exon_info$Stop[sel_idx]
  
  for(smpl_idx in 1:N_sample)
  {
    
    sel_sample_ID <- sample_ID_vec[smpl_idx]
    tmp_output_bed_file <- 
      paste(session_info$DIR$txt_output_DIR,sel_sample_ID,output_suffix,'.bed',sep='')
    tmp_output_txt_file <- 
      paste(session_info$DIR$txt_output_DIR,sel_sample_ID,output_suffix,'.txt',sep='')
    
    cna.res_My <- CNA(cnv.mtx[sel_idx,smpl_idx], chrom, maploc,
                      data.type=c("logratio"),sampleid=sel_sample_ID)
    seg.res_My <- segment(cna.res_My, verbose=0,...)
    
    seg_mtx <- seg.res_My$output
    #bed_mtx <- cbind(as.character(seg_mtx$chrom),seg_mtx$loc.start,seg_mtx$loc.end,seg_mtx$seg.mean)
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
        paste(session_info$DIR$plot_output_DIR,sel_sample_ID,output_suffix,'.pdf',sep='')
      pdf(CNV_seg_pdf_filename)
      plot(seg.res_My,ylim=plot.ylim,cex=plot.cex)
      dev.off()  
    } # if (is.plot)
    
  } # end of for(smpl_idx
} # end of CNV_segment function
  