# patternCNV code to perform somatic CNV detection
rm(list=ls())

arg <- commandArgs(TRUE) 

if(length(arg) != 3){
	stop("ERROR! Incorrect number of arguments. \nUSAGE: Rscript somatic.cnv.R [working dir] [patternCNV Rlib path] [config.ini path]")
}

setwd(arg[1])
.libPaths(arg[2])
config.ini.path <- arg[3]

print("#=== step 0. loading R funcs")
source(paste(arg[2],"/patCNV.load.package.R",sep=""), chdir=TRUE)

print("#=== step 1. initialize session")
somatic.sessionInfo <- patCNV.load.config(config.ini.path)

germline_count <- length(which(somatic.sessionInfo$file_info$sample.type=="Germline"))
germline_count

germline_samples <- subset(somatic.sessionInfo$file_info$sample.name,somatic.sessionInfo$file_info$sample.type=="Germline")
germline_samples

# check if there are enough germline samples with no duplicate samples
if(germline_count >= 3 & !anyDuplicated(germline_samples)){

	print("#=== step 2(a) germline coverage")
	germline.covg.res <- patCNV.scan.covg.multi(session_info=somatic.sessionInfo, sample.type='Germline', is.plot=FALSE)

	print("#=== step 2(b) germline coverage QC")
	germline.QC.res <- patCNV.sample.QC(session.info=somatic.sessionInfo, covg.info=germline.covg.res, output.prefix="Germline") 

	print("#=== step 2(c) germline GC correction")
	adj.germline.covg.res <- patCNV.adjust.GCbias(session.info=somatic.sessionInfo, covg.info=germline.covg.res, output.prefix="Germline")

	print("#=== step 2(d) pattern training and SNR summary")
	exon.pattern.res <- patCNV.learnExonPattern(session.info=somatic.sessionInfo, covg.info=adj.germline.covg.res, sample.QC.list=germline.QC.res,check.gender=TRUE)
	guessedMales<- attr(exon.pattern.res,"guessedMales")
	print(paste("returning from function: guessedMales=",guessedMales,sep="",collapse="\n"))
	predicted.sex<-attr(exon.pattern.res,"sex.predicted")
	known.sex<-attr(exon.pattern.res,"sex.known")

	if(somatic.sessionInfo$learn.Xratio==TRUE) {
             somatic.sessionInfo$X.FemalevsMale.ratio<-as.numeric(attr(exon.pattern.res,"X.FemalevsMale.ratio"))
	     print(paste("somatic.sessionInfo$X.FemalevsMale.ratio=",somatic.sessionInfo$X.FemalevsMale.ratio,"\n",sep=""))
	     adj.germline.covg.res<-patCNV.adjustXChrom(session.info=somatic.sessionInfo, covg.info=adj.germline.covg.res) 
        }
	if(length(guessedMales)>0) {
	   print(paste("Adjusting males N=",length(guessedMales)," guessed Males with X.FemalevsMale.ratio=",attr(exon.pattern.res,"X.FemalevsMale.ratio"),"\n",sep=""))
	   adj.germline.covg.res.Xcorrected<-patCNV.adjustGuessedMale(session.info=somatic.sessionInfo,covg.info=adj.germline.covg.res,guessedMales=attr(exon.pattern.res,"guessedMales"))
	} else {
	   adj.germline.covg.res.Xcorrected<- adj.germline.covg.res
	}

	print("#=== step 2(e) exon-level CNV calling; generate table and matrix summary")
	germline.exonCNV.res <- patCNV.exon.callCNV(session.info=somatic.sessionInfo, covg.info=adj.germline.covg.res.Xcorrected, 
			     pattern.list=exon.pattern.res, CNV.type="Germline",predicted.sex=predicted.sex,known.sex=known.sex)


	print("#=== step 2(f) exon-level CNV segmentation; generate txt and pdf results per sample")
	patCNV.exon.segmentCNV(session.info=somatic.sessionInfo, CNV.mtx=germline.exonCNV.res$CNV.mtx, pattern.list=exon.pattern.res)



	print("#=== step 3(a) somatic coverage")
	somatic.covg.res <- patCNV.scan.covg.multi(session_info=somatic.sessionInfo, sample.type='Somatic')

	print("#=== step 3(b) somatic GC-coverage bias correction")
	adj.somatic.covg.res <- patCNV.adjust.GCbias(session.info=somatic.sessionInfo, covg.info=somatic.covg.res, output.prefix="Somatic")

	if(somatic.sessionInfo$learn.Xratio==TRUE) {
# readjust X chromosome coverage for known males after learning the ratio.. Do not recommend this.
	     adj.somatic.covg.res<-patCNV.adjustXChrom(session.info=somatic.sessionInfo, covg.info=adj.somatic.covg.res) 
        }


        print("#=== step 3(c) somatic coverage QC")	 
	somatic.QC.res <- patCNV.sample.QC(session.info=somatic.sessionInfo, covg.info=adj.somatic.covg.res, output.prefix="Somatic") 


	print("#=== step 3(d) somatic sex-Check")
	sexcheck<-patCNV.Somatic.SexCheck(session.info=somatic.sessionInfo, covg.info=adj.somatic.covg.res,sample.QC.list=somatic.QC.res)
	somatic.guessedMales<- attr(sexcheck,"guessedMales")
	somatic.predicted.sex<-attr(sexcheck,"sex.predicted")
	somatic.known.sex<-attr(sexcheck,"sex.known")

	if(length(somatic.guessedMales)>0) {
	   print("Adjusting somatic X&Y counts for guessed males\n")
	   adj.somatic.covg.res.Xcorrected<-patCNV.adjustGuessedMale(session.info=somatic.sessionInfo,covg.info=adj.somatic.covg.res,guessedMales=somatic.guessedMales)
	} else {
	   adj.somatic.covg.res.Xcorrected<-adj.somatic.covg.res
        }
	attr(adj.somatic.covg.res.Xcorrected,"guessedMales")<-somatic.guessedMales


	print("#=== step 3(f) somatic exonCNV calling; text and figure results")
	somatic.CNV.res <- patCNV.exon.callCNV(session.info=somatic.sessionInfo, covg.info=adj.somatic.covg.res.Xcorrected, pattern.list=exon.pattern.res, CNV.type = "Somatic",
			           predicted.sex=somatic.predicted.sex,known.sex=somatic.known.sex)

	print("#=== step 3(g) somatic exonCNV segmentation; text and figure results")
	patCNV.exon.segmentCNV(session.info=somatic.sessionInfo, CNV.mtx=somatic.CNV.res$CNV.mtx, pattern.list=exon.pattern.res)






}else{
	print("You need at least 3 unique Germline samples to run this version of PatternCNV. Running the old version now which performs somatic calling in a pair-wise fashion...")

	.libPaths(paste(arg[2],"/old_somatic_version",sep=""))

	source(paste(arg[2],"/old_somatic_version/patCNV.load.package.R",sep=""),chdir=TRUE)
	somatic.session_info <- patCNV.load.config(config.ini.path)
	germline.covg_info <- patCNV.scan.covg.multi(session_info=somatic.session_info,sample.type='Germline',is.plot=FALSE)
	somatic.covg_info <- patCNV.scan.covg.multi(session_info=somatic.session_info,sample.type='Somatic',is.plot=FALSE)
	somatic.cnv <- patCNV.compute.CNV.multi(session_info=somatic.session_info,ref_type='basic.paired',sample.type='Somatic')
	somatic.cnv$sample.name
	# generate autosome CNV plots
	for(samp in somatic.cnv$sample.name)
	{
		png(paste(somatic.session_info$DIR$plot_output_DIR,"/",samp,".png",sep=""), height=1e3,width=3e3)
		patCNV.data('hg19.Chr_length')
		patCNV.plot.autosome.CNV(session_info=somatic.session_info,sample_name=samp,chr_length_info=hg19.Chr_length_vec,cnv_res=somatic.cnv,ref_avg_type='none')
		dev.off()
	}
	# generate CNV table
	patCNV.export.simple.tables(cnv_res=somatic.cnv,covg_info=somatic.covg_info,Log2R_threshold=0,session_info=somatic.session_info)


}

