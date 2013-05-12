patCNV.export.CNV.tables <- function(FDR_res, covg_info, session_info, 
				ref_avg_type='median', min_ref_avg_RPKM=3,
                               output_suffix='_CNV_table.txt',capture.only=TRUE,
                               FDR_type='localFDR',FDR_threshold=1e-3)
# refcovg_type = 'mean' or 'median'
# FDR_type = 'localFDR' or 'qvalue'
{
  cnv.mtx <- FDR_res$CNV
  sample_ID_vec <- FDR_res$sample.ID
  N_sample <- length(sample_ID_vec)
  
  patCNV.create.DIR(session_info$DIR$txt_output_DIR)
  
  
 if(ref_avg_type=='median')	
     {ref_avg_RPKM <- unlist(read.delim(session_info$Misc$median_RPKM_file,header=FALSE))}

   if(ref_avg_type=='mean')	
     {ref_avg_RPKM <- unlist(read.delim(session_info$Misc$mean_RPKM_file,header=FALSE))} 

  for(smpl_idx in 1:N_sample)
  {
    sel_sample_ID <- sample_ID_vec[smpl_idx]
    tmp_output_file <- paste(session_info$DIR$txt_output_DIR,sel_sample_ID,output_suffix,sep='')
    
    if(FDR_type=='localFDR')
    { sel_idx <- unique(which(FDR_res$LFDR.mtx[,smpl_idx] <= FDR_threshold)) }

    if(FDR_type=='qvalue')
    { sel_idx <- unique(which(FDR_res$qval.mtx[,smpl_idx] <= FDR_threshold)) }
    
    if(capture.only==TRUE)
    { non_cap_idx <- which(session_info$exon_info$InCapture==0)
      sel_idx <- setdiff(sel_idx,non_cap_idx) }
    
    inefficient_idx <- which(ref_avg_RPKM < min_ref_avg_RPKM)
    sel_idx <- setdiff(sel_idx,inefficient_idx)
    
    chr_vec   <- session_info$exon_info$Chr[sel_idx]
    start_vec <- session_info$exon_info$Start[sel_idx]
    stop_vec <- session_info$exon_info$Stop[sel_idx]
    gene_vec <- session_info$exon_info$Genes[sel_idx]
    cnv_vec <- cnv.mtx[sel_idx,sel_sample_ID]
    rawpval_vec <- FDR_res$rawp.mtx[sel_idx,sel_sample_ID]
    LFDR_vec <- FDR_res$LFDR.mtx[sel_idx,sel_sample_ID]
    qval_vec <- FDR_res$qval.mtx[sel_idx,sel_sample_ID]
    
    #=== for the two following matrices, columns are selected according to sample_ID,
    #    in order to prevent incorrect mapping
    count_vec <- covg_info$exon_count_mtx[sel_idx,sel_sample_ID] 
    RPKM_vec <- covg_info$exon_RPKM_mtx[sel_idx,sel_sample_ID]
    
    ref_avg_RPKM_vec <- ref_avg_RPKM[sel_idx]
    InCapture_vec <- session_info$exon_info$InCapture[sel_idx]
    
    out.mtx <- cbind(chr_vec,start_vec,stop_vec,gene_vec,
                     cnv_vec,rawpval_vec,LFDR_vec,qval_vec,count_vec,RPKM_vec,
			ref_avg_RPKM_vec,InCapture_vec)


#	print(dim(out.mtx))	# debugging
#	print(head(count_vec))
#	print(head(out.mtx))

   if(ref_avg_type=='median') {
	    colnames(out.mtx) <- c('Chr','exon.start','exon.stop','gene.symbol',
                           'CNV.log2ratio','raw.pval','localFDR','qval',
                           'observed.bp.counts','observed.RPKM',
				'average.RPKM(median)','is.captured')
		}

   if(ref_avg_type=='mean') {
	    colnames(out.mtx) <- c('Chr','exon.start','exon.stop','gene.symbol',
                           'CNV.log2ratio','raw.pval','localFDR','qval',
                           'observed.bp.counts','observed.RPKM',
				'average.RPKM(mean)','is.captured')
		}

    
    write.table(x=out.mtx,file=tmp_output_file,
                quote=FALSE,sep='\t',row.names=FALSE)
  }
  
}

