
patCNV.install.DIR <- getwd()

#========== Aug. 9. 2013

function_vec <- c(
		  "patCNV.misc.functions.R",
		  "patCNV.load.config.R",

		  "patCNV.scan.covg.single.R",
		  "patCNV.scan.covg.multi.R", 		  
		  "patCNV.coverage.QC.R",
		  "patCNV.chr.coverage.QC.R",
		  
      
      
	          "patCNV.learn.patterns.R",

		  "patCNV.compute.CNV.single.R",
		  "patCNV.compute.CNV.multi.R",			
		  "patCNV.fit.null.model.R" ,
		  "patCNV.estimate.FDR.R",

      
      "patCNV.export.simple.tables.R",
      
		   "patCNV.export.CNV.tables.R",
		   "patCNV.plot.Chr.CNV.R",
		   "patCNV.plot.autosome.CNV.R",
		   "patCNV.Gene.Heatmap.R",
		   "patCNV.segment.CNV.R"
                                            			)		 	

for (i in 1:length(function_vec))
{
  tmp_func <- paste(patCNV.install.DIR,'/functions/',function_vec[i],sep='')
  source(tmp_func)
}

cat('\n ========================================================= \n')
cat('\n R functions of pattern-CNV (ver 1.0) have been successfully loaded! \n')
cat('\n ========================================================= \n')
