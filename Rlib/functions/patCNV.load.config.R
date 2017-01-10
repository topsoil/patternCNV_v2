patCNV.load.config <- function( config_file )
{

 #=== checking and loading require libraries  
 patCNV.load.Rlib(lib_name='DNAcopy',lib_type='bioconductor')
 patCNV.load.Rlib(lib_name='RCurl',lib_type='CRAN')
 patCNV.load.Rlib(lib_name='fdrtool',lib_type='CRAN')
# patCNV.load.Rlib(lib_name='VGAM',lib_type='CRAN')
 patCNV.load.Rlib(lib_name='matrixStats',lib_type='CRAN')
 patCNV.load.Rlib(lib_name='gplots',lib_type='CRAN')

 options(stringsAsFactors = FALSE)
 
 #=== loading configurations  
 source(config_file)

 if(! exists("build")) {
     build="hg19"
     print("WARNING: Build not defined, assuming hg19")
 } else {
   if(! (build=="hg19" || build=="hg38")) {
   	stop('\n Error! build can only be hg19 or hg38')
   	return()
   }
}
if(! exists("minSexExons")) {
     minSexExons=30
}

if(!exists("learn.Xratio")) {
# Should set it to TRUE for deriving the X scaling factor from the reference set.
    learn.Xratio=FALSE
}
# If ratio is less than than, it is definitely not a male(ratio should be 1.0) (e.g. a Female or Klifelter)
if(!exists("YoverX.ratio.cut")) {
# former default was a ratio of 0.45, but that was before a PAR correction.
    YoverX.ratio.cut<-0.7
}
# If the Y over X ratio (YXX YXXX) is less than this, it may be Klinefelter
if(!exists("YoverX.Klinefelter.cut")) {
    YoverX.Klinefelter.cut<-0.7
}
if(!exists("X.FemalevsMale.ratio")) {
   X.FemalevsMale.ratio<-2.0
}

 #=== DIR information
 DIR_info <- list()
 DIR_info$plot_output_DIR <- patCNV.DIR.str(plot_output_DIR)
 DIR_info$txt_output_DIR <- patCNV.DIR.str(txt_output_DIR)
  
	patCNV.create.DIR( DIR_info$plot_output_DIR )
	patCNV.create.DIR( DIR_info$txt_output_DIR )

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
 file_info <- read.delim(sample_info_file,stringsAsFactors=FALSE,sep="\t", colClasses=c(rep("character", 6)))
 if("sex" %in% colnames(file_info)) {
    file_info$sex<-as.character(file_info$sex)
 }
 file_info$ID <- paste(file_info$sample.name,'(',file_info$sample.type,')',sep='') # unique ID = name + type
# New field.. should be file_info$sex
 cat('\n Number of samples:', nrow(file_info))
 cat('\n Sample type distributions:')
 print(table(file_info$sample.type))
 if("sex" %in% colnames(file_info)) {
     cat('\n Sample sex distributions:')
     print(table(file_info$sex))
     cat('\n')
 }

 for(k in 1:nrow(file_info))
 {
   tmp_file <- file_info$file.name[k]
   if(!patCNV.file.exists(tmp_file))
   {
     cat('\n Error! Input file \"',tmp_file,'\" cannot be located',sep='')
   }
 }
 cat('\n')


print(paste("nsamples=",length(file_info$ID),sep=""))
germlines.idx=which(file_info$sample.type=="Germline")
germlines<-file_info$ID[germlines.idx]
nGermline<-length(germlines)
print(paste("ngermlines=",nGermline,sep=""))
print(paste("nprobes=",length(exon_info$Chr),sep=""))
print(paste("nprobes.chrX=",length(which(exon_info$Chr=="chrX")),sep=""))
print(paste("nprobes.chrY=",length(which(exon_info$Chr=="chrY")),sep=""))
exon_info_mask<-matrix(data=c(1),ncol=length(file_info$ID),nrow=length(exon_info$Chr))
colnames(exon_info_mask)<-file_info$ID
#print(paste("dim(exon_info_mask)=",dim(exon_info_mask),sep=""))
if(exists("exclude_training_bed"))
{
# It's OK not to specify this file, but if specified, must be present.
	if(!patCNV.file.exists(exclude_training_bed)) {
           {
	       stop(paste("\n Error!Exclude training bed file cannot be located.n Please check for the existence (of even an empty) ", exclude_training_bed," bed file",sep=""))
	       return()  
           }
          exclude_regions<-cbind(Chr=c(),Start=c(),Stop=c(),sampleid=c())
          exclude_regions>-read.delim(exclude_training_bed,stringsAsFactors=FALSE,header=F)
	  colnames(exclude_regions)<-c("Chr","Start","Stop","sampleid")
	  exclude_regions$Chr=as.character(exclude_regions$Chr)
	  exclude_regions$Start=1+as.numeric(exclude_regions$Start)
	  exclude_regions$Stop=as.numeric(exclude_regions$Stop)
	  exclude_regions$sampleid=as.character(exclude_regions$sampleid)

	  exon_info$Chr<-as.character(exon_info$Chr)
	  exon_info$Start<-as.numeric(exon_info$Start)
	  exon_info$Stop<-as.numeric(exon_info$Stop)
	  exclude_regions$Start<-as.numeric(exclude_regions$Start)
	  exclude_regions$Stop<-as.numeric(exclude_regions$Stop)
# Just multiply coverage by this mask to omit masked out exons.

          for(i in 1:length(exclude_regions$Chr)) {
                samp.index<-which(germlines %in% exclude_regions$sample_id[i])
                exon_info_mask[unlist(which(exon_info$Chr==exclude_regions$Chr[i] & (
                (exon_info$Start>=exclude_regions$Start[i] & exon_info$Start<=exclude_regions$Stop[i])
                 | (exon_info$Stop>=exclude_regions$Start[i] & exon_info$Stop<=exclude_regions$Stop[i]
                 | (exon_info$Start<=exclude_regions$Start[i]  & exon_info$Stop>=exclude_regions$Stop[i] )
                 | (exon_info$Start>=exclude_regions$Start[i]  & exon_info$Stop<=exclude_regions$Stop[i])
		 | (exon_info$Start<=exclude_regions$Start[i]  & exon_info$Stop>=exclude_regions$Start[i])
		 | (exon_info$Start<=exclude_regions$Stop[i]  & exon_info$Stop>=exclude_regions$Stop[i])
		)))),samp.index]<-0
          }
       }		 
}


# Definitions of the PAR regions for build 37(hg19) .. These regions are always diploid
#chrY:10001-2649520 and chrY:59034050-59363566
#chrX:60001-2699520 and chrX:154931044-155260560
# for build 38 (hg20)
#chrY:10,000-2,781,479 and chrY:56,887,902-57,217,415
#chrX:10,000-2,781,479 and chrX:155,701,382-156,030,895
PARstartsX<-c(60001,154931044)
PARstartsY<-c(10001,59034050)
PARstopsX<-c(2699520,155260560)
PARstopsY<-c(2649520,59363566)

if(build == "hg38") {
   PARstartsX<-c(10000,155701382)
   PARstartsY<-c(10000,56887902)
   PARstopsX<-c(2781479,156030895)
   PARstopsY<-c(2781479,57217415)
}

PAR<-rep(0,length(exon_info$Chr))
PAR[which((exon_info$Chr=="chrX" & ((exon_info$Start>=PARstartsX[1] & exon_info$Stop<=PARstopsX[1]) 
				 | (exon_info$Start>=PARstartsX[2] & exon_info$Stop<=PARstopsX[2]) )) |
	   (exon_info$Chr=="chrY" & ((exon_info$Start>=PARstartsY[1] & exon_info$Stop<=PARstopsY[1]) 
	   | (exon_info$Start>=PARstartsY[2] & exon_info$Stop<=PARstopsY[2]) ))) ] <-1

exon_info$PAR=PAR



 res_list <- list(exon_info=exon_info,file_info=file_info,DIR_info=DIR_info,
 	     exon_info_mask=exon_info_mask,minSexExons=minSexExons,learn.Xratio=learn.Xratio,YoverX.ratio.cut=YoverX.ratio.cut,
	     YoverX.Klinefelter.cut=YoverX.Klinefelter.cut,X.FemalevsMale.ratio=X.FemalevsMale.ratio)


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

