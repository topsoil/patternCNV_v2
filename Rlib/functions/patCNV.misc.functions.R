#========================================================================
patCNV.DIR.str <- function( DIR_name )
{
  tmp_str <- ''
  if (substr(DIR_name,nchar(DIR_name),nchar(DIR_name))=="/")  
   { tmp_str <- DIR_name}  else {
     tmp_str <- paste(substr(DIR_name,1,nchar(DIR_name)),'/',sep='')
   }
  return(tmp_str)
}


#=========================================================================
patCNV.create.DIR <- function( DIR_name )
{
 if (!file.exists(DIR_name) ){
    #  cat('\"',DIR_name,'\" does not exist, creating directory...\n',sep='')
      dir.create(file.path(DIR_name),showWarnings=FALSE)
    }
}

#=========================================================================
patCNV.load.Rlib <- function( lib_name, lib_type=c('CRAN','Bioconductor') )
{

  if(tolower(lib_type)=='cran')	{ 
	
	lib_flag <- library(package=lib_name,logical.return=TRUE,character.only=TRUE)
	if(!lib_flag)
	  {
	    	cat(paste('\n library ',lib_name,' cannot be located. installing...\n',sep=''))
		install.packages(pkgs=lib_name)	
	 	stop(paste('\n please re-run the code\n',sep=''))
	    
	  }
	} 

  if(tolower(lib_type)=='bioconductor')	{ 
	
	lib_flag <- library(package=lib_name,logical.return=TRUE,character.only=TRUE)
	if(!lib_flag)
	  {
		stop(paste('\n library ',lib_name,' cannot be located. installing...\n',sep=''))
	    source("http://bioconductor.org/biocLite.R")
	    biocLite(pkgs=lib_name)	
	    cat(paste('\n please re-run the code\n',sep=''))
	  }
	} 

}


patCNV.data <- function(data_name, DIR_name=patCNV.install.DIR)
{
	tmp.file <- paste(DIR_name,'/data/',data_name,sep='')
 	cat(tmp.file,'... loaded \n')
	load(tmp.file,parent.frame(n = 1))
        cat('\n')	
}



patCNV.lap.pval <- function (q, location = 0, scale = 1) 
{
  if (!is.numeric(scale) | scale <=0) 
    { stop(" input scale must be positive number") }
  z.score <- (q - location)/scale
  L <- max(length(q), length(location), length(scale))
  q <- rep(q, length.out = L)
  location <- rep(location, length.out = L)
  scale <- rep(scale, length.out = L)
  p_tmp <- ifelse(q < location, 0.5 * exp(z.score), 1 - 0.5 * exp(-z.score))
  2*pmin(1-p_tmp,p_tmp)
}


