partAbund <-
function(comm,countThresh=NULL){
### F. Guillaume Blanchet - Décembre 2013
###########################################
	### General check
	if(any(is.na(comm))){
		stop("There are NAs in comm")
	}
	
	if(is.data.frame(comm)){
		print("'comm' was converted to a matrix")
		comm<-as.matrix(comm)
	}
	
	if(is.null(countThresh)){
		countThresh<-max(comm)
	}
	
	### Construct result object
	commPart<-array(dim=c(nrow(comm),ncol(comm),countThresh))
	
	### Calculates all partial-abundances community matrices
	for(i in 1:countThresh){
		commtmp<-comm
		tmp<-which(comm>(i-1),arr.ind=TRUE)
		commtmp[tmp]<-i
		commPart[,,i]<-commtmp
	}
	
	class(commPart)<-"partAbund"
	
	return(commPart)
}
