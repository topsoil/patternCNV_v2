patCNV.Gene.Heatmap <- function( cnv_res, session_info, sel_genes,
			         ref_avg_type='median',min_ref_avg_RPKM=3, capture.only=TRUE,
                                 logr.cut=1.5, font.cex=0.7, heatmap.margin=10,...)
#===== ... is for plot()
{
   if(ref_avg_type=='median')	
     {ref_avg_RPKM <- unlist(read.delim(session_info$Misc$median_RPKM_file,header=FALSE))}

   if(ref_avg_type=='mean')	
     {ref_avg_RPKM <- unlist(read.delim(session_info$Misc$mean_RPKM_file,header=FALSE))} 

  require(gplots)
  if(capture.only)
  {
    sel_exon_idx <- which(ref_avg_RPKM>=min_ref_avg_RPKM & session_info$exon_info$InCapture==1)  
  } else {
    sel_exon_idx <- which(ref_avg_RPKM>=min_ref_avg_RPKM)
  } 
  
  sel_gene_idx <- sel_exon_idx[
			which( is.element(session_info$exon_info$Genes[sel_exon_idx], 
                                   sel_genes) )]
  
  exon_label <- paste(session_info$exon_info$Genes[sel_gene_idx],
        paste(session_info$exon_info$Chr[sel_gene_idx],
        session_info$exon_info$Start[sel_gene_idx],sep=':'),sep=' @ ')
  
  heatmap.2(cnv_res$CNV[sel_gene_idx,],density.info='none',
            trace='none',
            labRow=exon_label,
            col=bluered(35),symbreaks=TRUE,margins=c(heatmap.margin,heatmap.margin),
	    cexRow=font.cex,cexCol=font.cex,	
    breaks=seq(-logr.cut,logr.cut,length.out=36),...)
}
  