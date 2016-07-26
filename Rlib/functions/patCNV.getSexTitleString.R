getSexTitleString<- function(predicted.sex,known.sex,k) {
# Supply a vector of "N" predicted and known sex. This function takes the "k-th" sample
    sexString<-""
    knownSex<-""
    predictedSex<-""
    if(length(predicted.sex)>=k) {
      predictedSex<-predicted.sex[k]
    }
    if(length(known.sex)>=k){
      knownSex<-known.sex[k]
    }
    if((nchar(knownSex)==0 || knownSex=="NA" || is.na(knownSex) || is.null(knownSex)) && (predictedSex=="NA" || nchar(predictedSex)==0 || is.na(predictedSex) || is.null(predictedSex))) {
       sexString<-"Sample of unspecified Gender\nNot able to perform Sex QC"
    } else if((!is.na(knownSex)) && knownSex=="MALE" ) {
       if((!is.na(predictedSex)) && predictedSex=="MALE") {
       	   sexString<-"Known MALE sample has rescaled X chromosome\nPasses Sex QC"
       } else if(is.na(predictedSex) || predictedSex=="NA" || nchar(predictedSex)==0 ) {
       	   sexString<-"Known MALE sample has rescaled X chromosome\nNot able to perform Sex QC"
       } else {
       	   sexString<-paste("Known MALE sample has rescaled X chromosome\nFailed Sex QC, looks like ",predictedSex,sep="")
       }
    } else if((!is.na(knownSex)) && (!is.na(predictedSex)) && knownSex!="MALE" && knownSex!="FEMALE" && predictedSex=="MALE") {
       sexString<-paste("Samples of unspecified sex(",knownSex,") has rescaled X chromosome\nLooks like MALE by Sex QC",sep="")
    } else if(is.na(predictedSex) || predictedSex=="NA" || nchar(predictedSex)==0 ) {
      sexString<-paste("Sample of Sex ",knownSex," : Not able to perform Sex QC",sep="")
    } else if(is.na(knownSex) || nchar(knownSex)==0 || knownSex=="NA" || is.null(knownSex) || (knownSex!="MALE" && knownSex!="FEMALE") ) {
      sexString<-paste("Unknown Sex sample (",knownSex,") predicted to be ",predictedSex," by Sex QC",sep="")
    } else if(predictedSex==knownSex) {
      sexString<-paste("Known ",knownSex," sample : Passes Sex QC",sep="")
    } else if(knownSex!=predictedSex) {
      sexString<-paste("Failed Sex QC: sample specified as ",knownSex," looks like it has sex ",predictedSex,sep="")
    } 
    return(sexString)
}
