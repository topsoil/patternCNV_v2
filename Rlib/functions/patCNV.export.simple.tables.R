patCNV.export.simple.tables <- function(session_info, covg_info, cnv_res,
                                     output_suffix='_simple_table.txt',capture.only=TRUE,
                                     Log2R_threshold=0.5)
  # refcovg_type = 'mean' or 'median'
  # FDR_type = 'localFDR' or 'qvalue'
{
  cnv.mtx <- cnv_res$CNV
  sample_ID_vec <- cnv_res$sample.ID
  N_sample <- length(sample_ID_vec)
  
  patCNV.create.DIR(session_info$DIR$txt_output_DIR)
  
   
  for(smpl_idx in 1:N_sample)
  {
    sel_sample_ID <- sample_ID_vec[smpl_idx]
    tmp_output_file <- paste(session_info$DIR$txt_output_DIR,sel_sample_ID,output_suffix,sep='')
       
    sel_idx <- unique(which(abs(cnv.mtx[,smpl_idx]) >= Log2R_threshold)) 
    
    if(capture.only==TRUE)
    { non_cap_idx <- which(session_info$exon_info$InCapture==0)
      sel_idx <- setdiff(sel_idx,non_cap_idx) }
    
    
    chr_vec   <- session_info$exon_info$Chr[sel_idx]
    start_vec <- session_info$exon_info$Start[sel_idx]
    stop_vec <- session_info$exon_info$Stop[sel_idx]
    gene_vec <- session_info$exon_info$Genes[sel_idx]
    cnv_vec <- cnv.mtx[sel_idx,sel_sample_ID]
        
    #=== for the two following matrices, columns are selected according to sample_ID,
    #    in order to prevent incorrect mapping
    count_vec <- covg_info$exon_count_mtx[sel_idx,sel_sample_ID] 
    RPKM_vec <- covg_info$exon_RPKM_mtx[sel_idx,sel_sample_ID]
    
    InCapture_vec <- session_info$exon_info$InCapture[sel_idx]
    
    out.mtx <- cbind(chr_vec,start_vec,stop_vec,gene_vec,
                     cnv_vec,count_vec,RPKM_vec,
                     InCapture_vec)
    
    
    #	print(dim(out.mtx))	# debugging
    #	print(head(count_vec))
    #	print(head(out.mtx))
    
    
    colnames(out.mtx) <- c('Chr','exon.start','exon.stop','gene.symbol',
                             'CNV.log2ratio',
                             'observed.bp.counts','observed.RPKM',
                             'is.captured')
    
    
    write.table(x=out.mtx,file=tmp_output_file, quote=FALSE,sep='\t',row.names=FALSE)
  }
  
}

