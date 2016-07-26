# patternCNV code to perform germline CNV detection
rm(list=ls())

arg <- commandArgs(TRUE) 

if(length(arg) != 3){
	stop("ERROR! Incorrect number of arguments. \nUSAGE: Rscript germline.cnv.R [working dir] [patternCNV Rlib path] [config.ini path]")
}

setwd(arg[1])
.libPaths(arg[2])
config.ini.path <- arg[3]

print("#=== step 0. loading R funcs")
source(paste(arg[2],"/patCNV.load.package.R",sep=""), chdir=TRUE)

print("#=== step 1. initialize session")
germline.sessionInfo <- patCNV.load.config(config.ini.path)

germline_count <- length(which(germline.sessionInfo$file_info$sample.type=="Germline"))
germline_count
if(germline_count < 3){
	stop("ERROR! You need at least 3 Germline samples to run this version of PatternCNV.")
}

print("#=== step 2(a) germline coverage")
germline.covg.res <- patCNV.scan.covg.multi(session_info=germline.sessionInfo, sample.type='Germline', is.plot=FALSE)

print("#=== step 2(b) germline coverage QC")
germline.QC.res <- patCNV.sample.QC(session.info=germline.sessionInfo, covg.info=germline.covg.res, output.prefix="Germline") 

#print("#=== step 2(c) germline GC correction")
adj.germline.covg.res <- patCNV.adjust.GCbias(session.info=germline.sessionInfo, covg.info=germline.covg.res, output.prefix="Germline")

print("#=== step 2(d) pattern training and SNR summary")
exon.pattern.res <- patCNV.learnExonPattern(session.info=germline.sessionInfo, covg.info=adj.germline.covg.res, sample.QC.list=germline.QC.res,check.gender=T)

guessedMales<- attr(exon.pattern.res,"guessedMales")
predicted.sex<-attr(exon.pattern.res,"sex.predicted")
known.sex<-attr(exon.pattern.res,"sex.known")

if(germline.sessionInfo$learn.Xratio==TRUE) {
     germline.sessionInfo$X.FemalevsMale.ratio<-as.numeric(attr(exon.pattern.res,"X.FemalevsMale.ratio"))
     print(paste("germline.sessionInfo$X.FemalevsMale.ratio=",germline.sessionInfo$X.FemalevsMale.ratio,"\n"))
     adj.germline.covg.res<-patCNV.adjustXChrom(session.info=germline.sessionInfo, covg.info=adj.germline.covg.res)
}
if(length(attr(exon.pattern.res,"guessedMales"))>0) {
     adj.germline.covg.res.Xcorrected<-patCNV.adjustGuessedMale(session.info=germline.sessionInfo,covg.info=adj.germline.covg.res,guessedMales=attr(exon.pattern.res,"GuessedMales"))
} else {
     adj.germline.covg.res.Xcorrected<- adj.germline.covg.res
}








print("#=== step 2(e) exon-level CNV calling; generate table and matrix summary")
germline.exonCNV.res <- patCNV.exon.callCNV(session.info=germline.sessionInfo, covg.info=adj.germline.covg.res.Xcorrected, 
		     pattern.list=exon.pattern.res, CNV.type="Germline",predicted.sex=predicted.sex,known.sex=known.sex)

print("#=== step 2(f) exon-level CNV segmentation; generate txt and pdf results per sample")
patCNV.exon.segmentCNV(session.info=germline.sessionInfo, CNV.mtx=germline.exonCNV.res$CNV.mtx, pattern.list=exon.pattern.res)

