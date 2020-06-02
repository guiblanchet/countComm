is.euclidean <-
function(distMat){
### F. Guillaume Blanchet - Décembre 2013
#########################################
	### General checks
	if(!inherits(distMat,"dist")){
		stop("distMat should be an object of class 'dist'")
	}
	if (any(distMat < .Machine$double.eps)){
		warning("Zero distance(s)")
	}
	
	### Calculate PCoA on distMat
	eigVal<-cmdscale(distMat,k=1,eig=TRUE)$eig
	neig<-length(eigVal)
	
	### Check of there are some "significant" negative values
	eucl <- !any(eigVal < -sqrt(.Machine$double.eps))
	
	return(eucl)
}
