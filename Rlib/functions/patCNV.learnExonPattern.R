patCNV.learnExonPattern <- function( session.info, covg.info,
                                     sample.QC.list,
                                     check.gender = TRUE, sex.minimum.RPKM = 0.1,
                                     small.delta = 1e-2,
				     text.output.DIR = 'cnv-txt')

{
# Changes:
#     12/09/2015 : Hugues Sicotte
#                  -Changed the meaning of signal-to-noise (was too low by 1.48)
#                  -Added SEX QC file
#                  -X-rescaling and Sex-check now ignore PAR1 and PAR2 regions.
#                  -User can specify BED file with known CNV in training sample, to ignore during training.
# 		   -remove outliers during X chromosome rescaling and sex check.
#                  -sex now required to be specified by user.
#     7/06/2016
#                  - disabled check.gender (it's always true)
#                  - synchronized code between learExonPattern and Somatic.SexCheck
#


     print("Running Germline Sex Check\n")

################################################################################################################################
####BEGIN of DO NOT CHANGE CODE HERE, change in patCNV.SexCheck.INCLUDE.R file.. and paste back                         ############################
################################################################################################################################
# Include this file in patCNV.Somatic.SexCheck.R  and patCNV.learnExonPattern.R

 sampleID.vec <- colnames(covg.info$exon_RPKM_mtx)
#  colnames.vec<-make.names(session.info$file_info$ID,unique=T)
#  sampleids<-session.info$file_info$ID[match(sampleID.vec,session.info$file_info$ID)]


  sample.NOToutlier.idx<- which(! (sampleID.vec %in% rownames(sample.QC.list)[which(sample.QC.list$is.outlier == 1)]))
  NSamples<-length(sampleID.vec) 
  X.FemalevsMale.ratio<-as.numeric(session.info$X.FemalevsMale.ratio)

  wasAbleToSexCheck<-FALSE

  male.sample.idx<-c()
  female.sample.idx<-c()

  diploid.subset<-which((session.info$exon_info$Chr!="chrX" & session.info$exon_info$Chr!="chrY"  ) | session.info$exon_info$PAR==1)
  haploid.subset<-which((session.info$exon_info$Chr=="chrX" | session.info$exon_info$Chr=="chrY"  ) & session.info$exon_info$PAR==0)


#  print(paste("sampleID.vec=",paste(sampleID.vec,sep="",collapse=","),sep=""))

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
#  print(paste(trueSex,sep="",collapse=","))


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
#  median.vec <- apply((exon_RPKM_mtx.masked ), 1, median,na.rm=T)
#  median.vec <- apply((exon_RPKM_mtx.masked[,,drop=F] ), 1, median,na.rm=T)

  MAD.vec <- apply((exon_RPKM_mtx.masked[,sample.NOToutlier.idx,drop=F] ), 1, mad,na.rm=T)
#  MAD.vec <- apply((exon_RPKM_mtx.masked[,,drop=F] ), 1, mad,na.rm=T)
#  MAD.vec <- apply((exon_RPKM_mtx.masked ), 1, mad,na.rm=T)

#  print("MAD.vec")
#  print(MAD.vec)

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
  print(paste("notPAR.X.subset=",length(notPAR.X.subset),sep=""))
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
       

     print("About to check if can perform X-only sex check")
     if((!wasAbleToSexCheck) && length(notPAR.X.subset)>=session.info$minSexExons){
     	        print("performing X-only sex check")
       		# No Y Chromosome probes, but can compare diploid coverage to X coverage to determine if Male or Female.


# not too low to be a male
	        sample.LooksMaleVec <- sample.Xoverdiploid.ratio>0.3 & sample.Xoverdiploid.ratio<=0.7
	        sample.LooksMale <- which(sample.LooksMaleVec)
		print(paste("sample.LooksMale=",paste(sample.LooksMale,sep="",collapse=","),sep=""))
# Too low to be a male 
	        sample.LooksFeMale <- which(sample.Xoverdiploid.ratio <=1.3 & sample.Xoverdiploid.ratio > 0.7)
		print(paste("sample.LooksFeMale=",paste(sample.LooksFeMale,sep="",collapse=","),sep=""))
# possibly an XXX or  XXXY_Klinefelter (Can only detect XXXY with this method but not XXY)
	        sample.LooksKlinefelter<-which(sample.Xoverdiploid.ratio>1.3)
	        sample.LooksTurner<-which(sample.LooksMaleVec & is.known.female.vec)


	        sex.predicted[sample.LooksMale]<-"MALE"
	        sex.predicted[sample.LooksFeMale]<-"FEMALE"
	        sex.predicted[sample.LooksKlinefelter]<-"KLINEFELTER_or_XXX(XXXY/XXX)"
# Note "Turner" overwrites MALE .. so code has to remain in this order.
	        sex.predicted[sample.LooksTurner]<-"TURNER(X)"

	    	print(paste("sex.predicted=",paste(sex.predicted,sep="",collapse=","),sep=""))
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
      guessedMales<-which((as.character(sex.predicted)=="MALE") &  (!(is.known.female.vec | is.known.male.vec)))
      print(paste("is.known.female.vec=",paste(is.known.female.vec,sep="",collapse=","),sep=""))
      print(paste("is.known.male.vec=",paste(is.known.male.vec,sep="",collapse=","),sep=""))
      print(paste("!is.known.vec=",paste( (!(is.known.male.vec |is.known.female.vec)),collapse=","),sep=""))

      print(paste("guessedMales=",paste(guessedMalessep="",collapse=","),sep=""))
   } else {
      guessedMales<-c()
    }


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
  write.table(x = as.data.frame(sample.sex.QC),file = paste(text.output.DIR, "/Germline_sex_QC_table.txt",sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)





# Compute the sex-adjusted Pattern.

  male.index<-sort(c(male.sample.idx,guessedMales))
  if(wasAbleToSexCheck){
     guessedFemales<-which(as.character(sex.predicted)=="FEMALE" & (! (is.known.female.vec | is.known.male.vec)))
  } else {
     guessedFemales<-c()
  }

  female.index<-sort(c(female.sample.idx,guessedFemales))

# Rescale X chromosome (not in PAR region) counts for males.
# At least some probes that can be corrected.

  if(length(sex.chrX.idx)>0) {
     chrX.RPKM.mtx.masked <- as.matrix(exon_RPKM_mtx.masked[sex.chrX.idx,])
     if(length(male.index)>0 && !is.null(male.index) && length(notPAR.X.subset)>0){
# Now multiple males by 2
	          chrX.RPKM.mtx.masked[notPAR.X.subset,male.index] <-  2*chrX.RPKM.mtx.masked[notPAR.X.subset,male.index]
     }
# Compute median for X&Y  probes
     chrX.median.vec<-NULL
     chrX.MAD.vec<-NULL

     if(length(female.index)>0) {
        chrX.median.vec <- apply( (chrX.RPKM.mtx.masked[,female.index,drop=F] ), 1, median,na.rm=T)
        chrX.MAD.vec <- apply( (chrX.RPKM.mtx.masked[,female.index,drop=F] ), 1, mad,na.rm=T)
     } else {
        if(length(male.index)>0) {
           chrX.median.vec <- apply( (chrX.RPKM.mtx.masked[,male.index,drop=F] ), 1, median,na.rm=T)
           chrX.MAD.vec <- apply( (chrX.RPKM.mtx.masked[,male.index,drop=F] ), 1, mad,na.rm=T)
        }
     }
# Overwrite the X&Y chromosome count with the adjusted valued.
#
     if(!is.null(chrX.median.vec)) {
         median.vec[sex.chrX.idx] <- chrX.median.vec
         MAD.vec[sex.chrX.idx] <- chrX.MAD.vec
     }
   }

    if(length(sex.chrY.idx)>0) {
       chrY.RPKM.mtx.masked <- as.matrix(exon_RPKM_mtx.masked[sex.chrY.idx,])
       chrY.median.vec <- NULL
       chrY.MAD.vec <- NULL	   
       if(length(male.index)>0 && !is.null(male.index) && length(notPAR.Y.subset)>0){
           chrY.median.vec <- apply(  chrY.RPKM.mtx.masked[,male.index,drop=F] , 1, median,na.rm=T)
      	   chrY.MAD.vec <- apply(  chrY.RPKM.mtx.masked[,male.index,drop=F] , 1, mad,na.rm=T )
       }

# Overwrite the X&Y chromosome count with the adjusted valued.
#
       if(!is.null(chrY.median.vec)){
          median.vec[sex.chrY.idx] <- chrY.median.vec
          MAD.vec[sex.chrY.idx] <- chrY.MAD.vec
       }
    }



# Rescale back by gender or guessed gender.
  if(length(notPAR.X.subset)>0){
          exon_count_mtx.masked[sex.chrX.idx[notPAR.X.subset],male.index]<- X.FemalevsMale.ratio * exon_count_mtx.masked[sex.chrX.idx[notPAR.X.subset],male.index]
  }
  if(length(notPAR.Y.subset)>0){
          exon_count_mtx.masked[sex.chrY.idx[notPAR.Y.subset],male.index]<-2*exon_count_mtx.masked[sex.chrY.idx[notPAR.Y.subset],male.index]
  }




#
# HS: on 12/08/2015 removed Extra 1.4825   SNR.ratio.vec <- (median.vec) / (1.4826 * (MAD.vec) + small.delta )
#
  SNR.ratio.vec <- (median.vec) / (MAD.vec + small.delta )

#  print("raw SNR.ratio.vec .. before fixing NA or taking log10\n")
#  print(SNR.ratio.vec)

  if(length(which(is.na(SNR.ratio.vec)))>0){	# fill SNR == NA with min.SNR 
	  min.SNR.ratio.val <- min(SNR.ratio.vec, na.rm=TRUE)
	  SNR.ratio.vec[which(is.na(SNR.ratio.vec))] <- min.SNR.ratio.val
  }


  SNR.dB.vec <- 20 * log10( SNR.ratio.vec )


  mean.count.vec <- apply(exon_count_mtx.masked,1,mean, na.rm=TRUE)

  patternMtx <- cbind( median.RPKM = median.vec,
                       MAD.RPKM = MAD.vec,
                       mean.raw.BPcount = mean.count.vec,
                       SNR.ratio = SNR.ratio.vec,
                       SNR.dB = SNR.dB.vec)

  patternMtxFrame<- as.data.frame(patternMtx)
  attr(patternMtxFrame,"X.FemalevsMale.ratio")<-X.FemalevsMale.ratio
  print(paste("Guessed Males=",sampleID.vec[guessedMales],sep="",collapse="\n"))
  attr(patternMtxFrame,"guessedMales")<-guessedMales
  attr(patternMtxFrame,"sex.predicted")<-sex.predicted
  attr(patternMtxFrame,"sex.known")<-trueSex

  return(patternMtxFrame)
}
