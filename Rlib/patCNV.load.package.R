
patCNV.install.DIR <- getwd()

#========== Nov. 18. 2014 (ver 1.5)
#=== added functions
# patCNV.sample.QC (.)
# patCNV.learnExonPattern (.)
# patCNV.adjust.GCbias (.)
# patCNV.plotSimpleCNV (.)
# patCNV.evaluate.CNVconf (.)
# patCNV.exon.callCNV (.)
# patCNV.exon.segmentCNV (.)

#========== Aug. 9. 2013 (ver 1.0)

function_vec <- c(
		  "patCNV.misc.functions.R",
		  "patCNV.load.config.R",

		  "patCNV.scan.covg.single.R",
		  "patCNV.scan.covg.multi.R", 		  
		  "patCNV.coverage.QC.R",
		  "patCNV.chr.coverage.QC.R",
 
      
	      "patCNV.learn.patterns.R",
		  "patCNV.learn.patterns.ME.R",


		  "patCNV.compute.CNV.single.R",
		  "patCNV.compute.CNV.multi.R",			
		  "patCNV.fit.null.model.R" ,
		  "patCNV.estimate.FDR.R",

      
			"patCNV.export.simple.tables.R",
      
		   "patCNV.export.CNV.tables.R",
		   "patCNV.plot.Chr.CNV.R",
		   "patCNV.plot.autosome.CNV.R",
		   "patCNV.Gene.Heatmap.R",
		   "patCNV.segment.CNV.R",
		   
		   #====================== v1.5
			"patCNV.sample.QC.R",
			"patCNV.learnExonPattern.R",
			"patCNV.adjust.GCbias.R",
			"patCNV.plotSimpleCNV.R",
			"patCNV.evaluate.CNVconf.R",
			"patCNV.exon.callCNV.R",
			"patCNV.exon.segmentCNV.R"
                                            			)		 	

for (i in 1:length(function_vec))
{
  tmp_func <- paste(patCNV.install.DIR,'/functions/',function_vec[i],sep='')
  source(tmp_func)
}

cat('\n ========================================================================= \n')
cat('\n R functions of pattern-CNV (ver 1.5) have been successfully loaded! \n')
cat('\n To cite software \"Pattern-CNV\" in publications use: \n')
cat('  Wang C, Evans JM, Bhagwate AV et al. \"PatternCNV: a versatile tool for detecting \n  copy number changes from exome sequencing data\". ')
cat('Bioinformatics. 2014 Sep 15;30(18):2678-80. \n')
cat('\n ========================================================================= \n')