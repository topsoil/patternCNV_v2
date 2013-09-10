patCNV.load.config <- function( config_file )
{

 #=== checking and loading require libraries  
 patCNV.load.Rlib(lib_name='DNAcopy',lib_type='bioconductor')
 patCNV.load.Rlib(lib_name='RCurl',lib_type='CRAN')
 patCNV.load.Rlib(lib_name='fdrtool',lib_type='CRAN')
# patCNV.load.Rlib(lib_name='VGAM',lib_type='CRAN')
 patCNV.load.Rlib(lib_name='matrixStats',lib_type='CRAN')
 patCNV.load.Rlib(lib_name='gplots',lib_type='CRAN')

 #=== loading configurations  
 source(config_file)

 #=== DIR information
 DIR_info <- list()
 DIR_info$plot_output_DIR <- patCNV.DIR.str(plot_output_DIR)
 DIR_info$txt_output_DIR <- patCNV.DIR.str(txt_output_DIR)
  
 
 #=== detecting existence of important files
 if(!patCNV.file.exists(exon_key_file))
 {
   str1 <- '\n Error! Exon key file cannot be located. \n Please check \"exon_key_file\" in \"'
   str2 <- '\" and re-run the code. \n'
   stop( paste(str1,config_file,str2,sep='') )
   return()
 }
  
 if(!patCNV.file.exists(sample_info_file))
 {
   str1 <- '\n Error! sample information file cannot be located. \n Please check \"sample_info_file\" in \"'
   str2 <- '\" and re-run the code. \n'
   stop( paste(str1,config_file,str2,sep='') )
   return()
 }
 
 #=== load exon_key
 cat('\n Loading exon information from \"', exon_key_file, '\"',sep='')
 exon_info <- read.delim(exon_key_file,stringsAsFactors=FALSE)
 
 if(substr(exon_info$Chr[1],1,1)!='c'){
   exon_info$Chr <- paste('chr',exon_info$Chr,sep='')   # if not starting with 'chr'
 }
   
 
 exon_info$exon_bin_vec <- as.numeric(exon_info$Bin_Count)
 exon_info$is_capture_vec <- as.numeric(exon_info$InCapture)
 cat('\n Number of exons in capture:',length(which(exon_info$is_capture_vec==1)))
 cat('\n Number of exons not in capture:',length(which(exon_info$is_capture_vec==0)))
 cat('\n')
 
 #=== load sample_info
 cat('\n Loading sample information from \"', sample_info_file, '\"',sep='')
 file_info <- read.delim(sample_info_file,stringsAsFactors=FALSE)
 file_info$ID <- paste(file_info$sample.name,'(',file_info$sample.type,')',sep='') # unique ID = name + type
 cat('\n Number of samples:', nrow(file_info))
 cat('\n Sample type distributions:')
 print(table(file_info$sample.type))
 cat('\n')
 
 for(k in 1:nrow(file_info))
 {
   tmp_file <- file_info$file.name[k]
   if(!patCNV.file.exists(tmp_file))
   {
     cat('\n Error! Input file \"',tmp_file,'\" cannot be located',sep='')
   }
 }
 cat('\n')
 
 res_list <- list(exon_info=exon_info,file_info=file_info,DIR_info=DIR_info)


 if(exists("avg_pattern_file")) 	
	{ res_list$pattern$avg_wig <- avg_pattern_file } else 
		{ res_list$pattern$avg_wig <- 'avg_pattern.wig' }

 if(exists("var_pattern_file")) 	
	{ res_list$pattern$var_wig <- var_pattern_file } else 
		{ res_list$pattern$var_wig <- 'var_pattern.wig' }

 if(exists("median_RPKM_file")) 	
	{ res_list$Misc$median_RPKM_file <- median_RPKM_file } else 
		{ res_list$Misc$median_RPKM_file <- 'median_RPKM.txt' }

 if(exists("mean_RPKM_file")) 	
	{ res_list$Misc$mean_RPKM_file <- mean_RPKM_file } else 
		{ res_list$Misc$mean_RPKM_file <- 'mean_RPKM.txt' }

 
 return(res_list)
 
}

