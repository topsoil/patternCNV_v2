patCNV.Somatic.SexCheck <- function( session.info, covg.info,
                                     sample.QC.list,
                                     sex.minimum.RPKM = 0.1,
				     text.output.DIR = 'cnv-txt') {
# New Module to perform sex QC.
# The difficulty is that if the sex was given as MALE, the counts were multiplied by 2.0 (X.FemalevsMale.ratio)


#  source does not behave like "include", delayed execution makes it a problem.
#######################
#######################
#######################
################################################################################################################################
####BEGIN of DO NOT CHANGE CODE HERE, change in patCNV.SexCheck.INCLUDE.R file.. and paste back                         ############################
################################################################################################################################





 sampleID.vec <- colnames(covg.info$exon_RPKM_mtx)
#  colnames.vec<-make.names(session.info$file_info$ID,unique=T)
#  sampleids<-session.info$file_info$ID[match(sampleID.vec,colnames.vec)]


  sample.NOToutlier.idx<- which(! (sampleID.vec %in% rownames(sample.QC.list)[which(sample.QC.list$is.outlier == 1)]))
  NSamples<-length(sampleID.vec) 
  X.FemalevsMale.ratio<-as.numeric(session.info$X.FemalevsMale.ratio)

  wasAbleToSexCheck<-FALSE


  male.sample.idx<-c()
  female.sample.idx<-c()

  diploid.subset<-which((session.info$exon_info$Chr!="chrX" & session.info$exon_info$Chr!="chrY"  ) | session.info$exon_info$PAR==1)
  haploid.subset<-which((session.info$exon_info$Chr=="chrX" | session.info$exon_info$Chr=="chrY"  ) & session.info$exon_info$PAR==0)


  is.known.female.vec<-rep(FALSE,NSamples)
  is.known.male.vec<-rep(FALSE,NSamples)
  trueSex<-rep("NA",NSamples)

  if("sex" %in% colnames(session.info$file_info)){
  # The count data is in the same order as the file_info ..
      trueSex<-as.character(session.info$file_info$sex[match(sampleID.vec,session.info$file_info$ID)])
      male.sample.idx<-which(as.character(trueSex)=="MALE")
      is.known.male.vec[male.sample.idx]<-TRUE
      female.sample.idx<-which(as.character(trueSex)=="FEMALE")
      is.known.female.vec[female.sample.idx]<-TRUE
  }
  female.NOToutlier.idx<-female.sample.idx[which(female.sample.idx %in% sample.NOToutlier.idx)]
  male.NOToutlier.idx<-male.sample.idx[which(male.sample.idx %in% sample.NOToutlier.idx)]


# Set to NA known CNV on specific samples, to get 'cleaner' pattern
# This exon_info_mask should cover both Somatic and Germline samples 
  OKflag.raw<-session.info$exon_info_mask
  samples.include<-which(sampleID.vec %in% colnames(OKflag.raw))
  OKflag<-OKflag.raw[,samples.include,drop=FALSE]
  if(length(samples.include)!=length(sampleID.vec)) {
      print("BUG: exon columns are not a subset of exon_info_mask\nDIAGNOSTIC INFO FOLLOWS\n----------------------------------\n")
      print(paste("colnames for OKflag",colnames(OKflag.raw),sep="\n"))
      print(paste("colnames for exon_RPKM_mtx",sampleID.vec,sep="\n"))
      print(paste("dim(OKflag)=",dim(OKflag),"\n------------------------------------",sep=""))
      exit(-1)
  }

  if(length(which(OKflag==0))>0){
    print(paste("Masked out =",length(which(OKflag==0))," exons\n",sep=""))
    OKflag[which(OKflag==0)]<-NA
    exon_RPKM_mtx.masked<-OKflag*covg.info$exon_RPKM_mtx
    exon_count_mtx.masked<-OKflag*covg.info$exon_count_mtx
  } else {
    exon_RPKM_mtx.masked<-covg.info$exon_RPKM_mtx
    exon_count_mtx.masked<-covg.info$exon_count_mtx
  }

  median.vec <- apply((exon_RPKM_mtx.masked[,sample.NOToutlier.idx,drop=F] ), 1, median,na.rm=T)
  MAD.vec <- apply((exon_RPKM_mtx.masked[,sample.NOToutlier.idx,drop=F] ), 1, mad,na.rm=T)

# undo male correction,(doubling) so can perform sex check and better learn Sex Ratio
   if(length(male.sample.idx)>0){
      exon_RPKM_mtx.masked[haploid.subset,male.sample.idx]<- 0.5*exon_RPKM_mtx.masked[haploid.subset,male.sample.idx]
      exon_count_mtx.masked[haploid.subset,male.sample.idx]<- 0.5*exon_count_mtx.masked[haploid.subset,male.sample.idx]
#      print("Undid Haploid correction")
  }


#
# check gender and adjust exon patterns by gender  
#            Only reason not to do this is if one has an all female cohort (or all male) .. or panel with no Y coverage.
#

  sex.chrX.idx <- which(session.info$exon_info$Chr=="chrX")
  sex.chrY.idx <- which(session.info$exon_info$Chr=="chrY")



# Even if there are no probes with data.. once check.gender check is requested, must output a report.. so precompute empty QC values.

     sex.predicted<-rep("NA",NSamples)

     notPAR.X.subset<-which(session.info$exon_info$PAR[sex.chrX.idx]==0)

     cat(paste("Number of haploid X probes =",length(notPAR.X.subset),"\n"))


     n.Xprobes<-rep(length(notPAR.X.subset),NSamples)
     n.Xprobes.OKcov<-rep(0,NSamples)
     notPAR.Y.subset <- which(session.info$exon_info$PAR[sex.chrY.idx]==0 )
     n.Yprobes<-rep(length(notPAR.Y.subset),NSamples)
     n.Yprobes.OKcov<-rep(0,NSamples)

# Load preliminary coverage .. it may be null.

     sample.diploid.meanCvg <- apply(exon_RPKM_mtx.masked[diploid.subset,,drop=F], 2, mean,na.rm=T)

     if(length(notPAR.Y.subset)>0){
             sample.chrY.meanCvg <- apply(exon_RPKM_mtx.masked[sex.chrY.idx[notPAR.Y.subset],,drop=F], 2, mean,na.rm=T)
     } else {
             sample.chrY.meanCvg<-rep(NA,ncol(exon_RPKM_mtx.masked))
     }
     if(length(notPAR.X.subset)>0){
             sample.chrX.meanCvg <- apply(exon_RPKM_mtx.masked[sex.chrX.idx[notPAR.X.subset],,drop=F], 2, mean,na.rm=T)
     } else {
             sample.chrX.meanCvg<-rep(NA,ncol(exon_RPKM_mtx.masked))
     }

     sample.LooksMale<-c()



     if(length(notPAR.X.subset)==0){
     	  sample.chrX.meanCvg<-rep(NA,ncol(exon_RPKM_mtx.masked))
     } 

# Enough probes to Do Sex Check & (optionally) learn X scaling factor for males.
        if(length(notPAR.X.subset)>0){
            # vector of probe means
             sample.chrX.meanCvg <- apply(exon_RPKM_mtx.masked[sex.chrX.idx[notPAR.X.subset],,drop=F], 2, mean,na.rm=T)
             if(length(female.NOToutlier.idx)>0) {
                 OK.X.notPar.probes<-which(apply(exon_RPKM_mtx.masked[sex.chrX.idx[notPAR.X.subset],female.NOToutlier.idx,drop=F], 1, mean,na.rm=T)>=sex.minimum.RPKM)
             } else {
#         use the sex-corrected scaled data for X.
                 if(length(sample.NOToutlier.idx)>0) {
                     OK.X.notPar.probes<-which(apply(covg.info$exon_RPKM_mtx[sex.chrX.idx[notPAR.X.subset],sample.NOToutlier.idx,drop=F], 1, mean,na.rm=T)>=sex.minimum.RPKM)
                 } else {
                   OK.X.notPar.probes<-c()
                 }
             }
          }
          if(length(notPAR.Y.subset)>0){
             if(length(male.NOToutlier.idx)>0) {
                  OK.Y.notPar.probes <- which(apply( exon_RPKM_mtx.masked[sex.chrY.idx[notPAR.Y.subset],male.NOToutlier.idx,drop=F], 1, mean,na.rm=T)>=sex.minimum.RPKM)
             } else {
             # We do not know who is male, so we will underestimate Y counts.
                  if(length(sample.NOToutlier.idx)>0) {
                     OK.Y.notPar.probes<-which(apply( exon_RPKM_mtx.masked[sex.chrY.idx[notPAR.Y.subset],sample.NOToutlier.idx,drop=F], 1, mean,na.rm=T)>=sex.minimum.RPKM)
                  } else {
                     OK.Y.notPar.probes<-c()
                  }
             }
          }

       # Per Sample mean coverage, eliminating NA .. also evaluate sex for outliers.

          sample.Xoverdiploid.ratio <- NA
          sample.Yoverdiploid.ratio <- NA
          if(length(notPAR.X.subset)>0 && length(OK.X.notPar.probes)>0){
	        reliable.sex.chrX.idx <-  sex.chrX.idx[notPAR.X.subset[OK.X.notPar.probes]]
	        n.Xprobes.OKcov<-rep(length(reliable.sex.chrX.idx),NSamples)
	    # Per Sample mean coverage, eliminating NA .. also evaluate sex for outliers.
	        sample.chrX.meanCvg <- apply(exon_RPKM_mtx.masked[reliable.sex.chrX.idx,,drop=F], 2, mean,na.rm=T)
	        sample.Xoverdiploid.ratio <- (sample.chrX.meanCvg/sample.diploid.meanCvg)
          }
          if(length(notPAR.Y.subset)>0 && length(OK.Y.notPar.probes)>0){
	        reliable.sex.chrY.idx <-  sex.chrY.idx[notPAR.Y.subset[OK.Y.notPar.probes]]
	        n.Yprobes.OKcov<-rep(length(reliable.sex.chrY.idx),NSamples)
	    # Per Sample mean coverage, eliminating NA .. also evaluate sex for outliers.
	        sample.chrY.meanCvg <- apply(exon_RPKM_mtx.masked[reliable.sex.chrY.idx,,drop=F], 2, mean,na.rm=T)
	        sample.Yoverdiploid.ratio <- (sample.chrY.meanCvg/sample.diploid.meanCvg)
          }

          if(length(notPAR.X.subset)>=session.info$minSexExons && length(notPAR.Y.subset)>=session.info$minSexExons){
#   refine exons with enough coverage to perform sex check.
             if(length(OK.X.notPar.probes)>=session.info$minSexExons && length(OK.Y.notPar.probes)>=session.info$minSexExons){

# YoverX.ratio.cut <-0.45 : Male, X and Y counts should be the same. Female, X counts >> Y counts (should be almost no Y counts)
# old version of PatternCNV(12/07/2015) had ratio threshold at 0.45
# YvsX.Klinefelter.cut=0.3
	        sample.YoverX.ratio <- (sample.chrY.meanCvg/sample.chrX.meanCvg)
	        sample.Yoverdiploid.ratio <- (sample.chrY.meanCvg/sample.diploid.meanCvg)
# possibly a XYY (male with extra Y .. aka Jacobs)
	        sample.LooksXYYvector<-   sample.Yoverdiploid.ratio>0.7 & sample.Yoverdiploid.ratio<=1.3 & sample.Xoverdiploid.ratio<=0.7 & sample.Xoverdiploid.ratio>0.3 
		sample.LooksXYY<-which(sample.LooksXYYvector)
# possibly a XXYY 
	        sample.LooksXXYYvector<-sample.Yoverdiploid.ratio>0.7 & sample.Xoverdiploid.ratio>0.7
		sample.LooksXXYY<-which(sample.LooksXYYvector)
# Y counts highenough to be a male .. but not high enough to be XYY
	        sample.LooksMalevector <- sample.Yoverdiploid.ratio > 0.3 & sample.Yoverdiploid.ratio <=0.7  &  sample.Xoverdiploid.ratio<=0.7 & sample.Xoverdiploid.ratio>0.3 & !(sample.LooksXYYvector | sample.LooksXXYYvector)
	        sample.LooksMale <- which(sample.LooksMalevector)
# Y counts Too low to be a male
	        sample.LooksFeMale <- which(sample.Yoverdiploid.ratio <=0.1 & sample.YoverX.ratio<0.1  & sample.Xoverdiploid.ratio>0.7 & sample.Xoverdiploid.ratio<=1.3 )
# possibly a Klinefelter (XXY or XXXY)
	        sample.LooksKlinefelter<-which(sample.YoverX.ratio <= session.info$YoverX.Klinefelter.cut & sample.YoverX.ratio>0.1 & sample.Xoverdiploid.ratio>0.7 & sample.Yoverdiploid.ratio>0.3 & sample.Yoverdiploid.ratio<=0.7)
# possibly a Turner Syndrome (female (no Y) with a single X)
	        sample.LooksTurner<-which(sample.YoverX.ratio <=0.1 & sample.Yoverdiploid.ratio<=0.1  &  sample.Xoverdiploid.ratio<=0.7 & sample.Xoverdiploid.ratio>0.3)
# possibly a XXX (female (no Y) with 3 X)
	        sample.LooksXXX<-which(sample.YoverX.ratio <=0.1 & sample.Yoverdiploid.ratio<=0.1  &  sample.Xoverdiploid.ratio>1.3)


	        sex.predicted[sample.LooksMale]<-"MALE"
	        sex.predicted[sample.LooksFeMale]<-"FEMALE"
# female (no Y) with a single X
	        sex.predicted[sample.LooksTurner]<-"TURNER(X-female)"
# XXY XXXY
                sex.predicted[sample.LooksKlinefelter]<-"KLINEFELTER(XXY/XXXY)"
	        sex.predicted[sample.LooksXXX]<-"XXX"
	        sex.predicted[sample.LooksXYY]<-"XYY"



# XFemalevsMale.ratio<-2
# res_list <- list(exon_info=exon_info,file_info=file_info,DIR_info=DIR_info,exon_info_mask,minSexExons,learn.Xratio,YoverX.ratio.cut,YoverX.Klinefelter.cut,X.FemalevsMale.ratio)
	    
	        OKmale.sample.idx<-which(sex.predicted=="MALE" & trueSex=="MALE")
	        OKfemale.sample.idx<-which(sex.predicted=="FEMALE" & trueSex=="FEMALE")

#          (optionally) learn Male vs Female ratio for X.. Should be 2.0, but could vary 

	        if(session.info$learn.Xratio==TRUE && length(OKmale.sample.idx)>0 && length(OKfemale.sample.idx)>0){
	      	   X.FemalevsMale.ratio<-mean(sample.chrX.meanCvg[OKfemale.sample.idx])/mean(sample.chrX.meanCvg[OKmale.sample.idx])
	        }
		wasAbleToSexCheck<-TRUE
		 print("wasAbleToSexCheck with X&Y probes\n")
             } # End of Enough probes to perform Gender Check and (optionally) learn Correction
       } # End of Enough probes to perform Gender Check and (optionally) learn Correction
       

       if((!wasAbleToSexCheck) && length(notPAR.X.subset)>=session.info$minSexExons){
       		# No Y Chromosome probes, but can compare diploid coverage to X coverage to determine if Male or Female.


# not too low to be a male
	        sample.LooksMaleVec <- sample.Xoverdiploid.ratio>0.3 & sample.Xoverdiploid.ratio<=0.7
	        sample.LooksMale <- which(sample.LooksMaleVec)
# Too low to be a male 
	        sample.LooksFeMale <- which(sample.Xoverdiploid.ratio <=1.3 & sample.Xoverdiploid.ratio > 0.7)
# possibly an XXX or  XXXY_Klinefelter (Can only detect XXXY with this method but not XXY)
	        sample.LooksKlinefelter<-which(sample.Xoverdiploid.ratio>1.3)
	        sample.LooksTurner<-which(sample.LooksMaleVec & is.known.female.vec)


	        sex.predicted[sample.LooksMale]<-"MALE"
	        sex.predicted[sample.LooksFeMale]<-"FEMALE"
	        sex.predicted[sample.LooksKlinefelter]<-"KLINEFELTER_or_XXX(XXXY/XXX)"
# Note "Turner" overwrites MALE .. so code has to remain in this order.
	        sex.predicted[sample.LooksTurner]<-"TURNER(X)"

	    
	        OKmale.sample.idx<-which(sex.predicted=="MALE" & trueSex=="MALE")
	        OKfemale.sample.idx<-which(sex.predicted=="FEMALE" & trueSex=="FEMALE")

#          (optionally) learn Male vs Female ratio for X.. Should be 2.0, but could vary 

	        if(session.info$learn.Xratio==TRUE && length(OKmale.sample.idx)>0){
	      	   X.FemalevsMale.ratio<-mean(sample.diploid.meanCvg)/mean(sample.chrX.meanCvg[OKmale.sample.idx])
	        }
		 print("wasAbleToSexCheck with X probes only\n")
		wasAbleToSexCheck<-TRUE

        }


       if((!wasAbleToSexCheck) && length(notPAR.Y.subset)>=session.info$minSexExons){
       		# No X Chromosome probes, but can compare diploid coverage to Y coverage to determine if Male or Female.
# not too low to be a male
	        sample.LooksMaleVec <- sample.Yoverdiploid.ratio>0.3 & sample.Yoverdiploid.ratio<=0.7
	        sample.LooksMale <- which(sample.LooksMaleVec)
# Too low to be a male 
	        sample.LooksFeMale <- which(sample.Yoverdiploid.ratio <=0.3)
	        sample.LooksXYY<-which(sample.Yoverdiploid.ratio>0.7)
	        sex.predicted[sample.LooksMale]<-"MALE"
	        sex.predicted[sample.LooksFeMale]<-"FEMALE"
	        sex.predicted[sample.LooksXYY]<-"XYY"
	    
	        OKmale.sample.idx<-which(sex.predicted=="MALE" & trueSex=="MALE")
	        OKfemale.sample.idx<-which(sex.predicted=="FEMALE" & trueSex=="FEMALE")

#          (optionally) learn Male vs Female ratio for X.. Should be 2.0, but could vary 
		 print("wasAbleToSexCheck with Y probes only\n")
		wasAbleToSexCheck<-TRUE
        }

   if(wasAbleToSexCheck){
      guessedMales<-which((as.character(sex.predicted)=="MALE" ) &  (!(is.known.female.vec | is.known.male.vec)))
   } else {
      guessedMales<-c()
    }



#######################
#######################
#######################
################################################################################################################################
####END of DO NOT CHANGE CODE HERE, change in INCLUDE file.. and paste back                         ############################
################################################################################################################################






  sample.sex.QC=data.frame(sample=sampleID.vec,sex.supplied=trueSex,
                  sex.predicted=sex.predicted,
                  chrX.mean=sample.chrX.meanCvg,
		  chrY.mean=sample.chrY.meanCvg,
		  Diploid.mean=sample.diploid.meanCvg,
	          n.Xprobes=n.Xprobes,
		  n.Xprobes.OKcov= n.Xprobes.OKcov,
		  n.Yprobes= n.Yprobes,
		  n.Yprobes.OKcov)
 rownames(sample.sex.QC)<-NULL

# Write out sex QC table
  write.table(x = as.data.frame(sample.sex.QC),file = paste(text.output.DIR, "/Somatic_sex_QC_table.txt",sep = ""), 
    quote = FALSE, sep = "\t", row.names = FALSE)




  retObj<-list()
  attr(retObj,"guessedMales")<-guessedMales
  attr(retObj,"sex.predicted")<-sex.predicted
  attr(retObj,"sex.known")<-trueSex
  return(retObj)
}
