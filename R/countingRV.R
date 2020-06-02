countingRV <-
function(comm,partAbund,transfo=NULL,dist=NULL,vegdist=NULL,designdist=NULL,logbase=NULL,...){
### F. Guillaume Blanchet - Décembre 2013
#########################################
	### General checks
	if(!inherits(partAbund,"partAbund")){
		stop("partAbund needs to be a 'partAbund' class object")
	}
	if(!(nrow(comm)==dim(partAbund)[1] & ncol(comm)==dim(partAbund)[2])){
		stop("The first two dimensions of partAbund needs to be the same as the dimension of comm")
	}
	if(length(which(c(!is.null(dist),!is.null(vegdist),!is.null(designdist))))>=2){
		stop("Only one of 'dist', 'vegdist' and 'designdist' argument should be defined")
	}
	if(!is.matrix(comm)){
		print("'comm' was converted to a matrix")
		comm<-as.matrix(comm)
	}
	
	### Match argument
#	transfo<-match.arg(transfo)
#	dist<-match.arg(dist)
#	vegdist<-match.arg(vegdist)
	
	### Counting procedure
	res<-numeric()
	if((is.null(transfo) & is.null(dist))){
		for(i in 1:dim(partAbund)[3]){
			res[i]<-RV(comm,partAbund[,,i])
		}
	}
	
	if(!is.null(transfo)){
		commTransfo<-decostand(comm,method=transfo)
		for(i in 1:dim(partAbund)[3]){
			res[i]<-RV(commTransfo,decostand(partAbund[,,i],method=transfo,...))	
		}
	}
	
	if(!is.null(dist)){
		Dist<-dist(comm,method=dist)
		nsp<-ncol(comm)
		if(is.euclidean(Dist)){
			commDist<-cmdscale(Dist,k=nsp-1)$points
			for(i in 1:dim(partAbund)[3]){
				res[i]<-RV(commDist,cmdscale(dist(partAbund[,,i],method=dist),k=nsp-1)$points)
			}
		}else{
			commDist<-cmdscale(Dist,k=nsp-1,add=TRUE)$points
			for(i in 1:dim(partAbund)[3]){
				res[i]<-RV(commDist,cmdscale(dist(partAbund[,,i],method=dist),k=nsp-1,add=TRUE)$points)
			}
			print("A Cailliez correction was applied")
		}
	}
	
	if(!is.null(vegdist)){
		Dist<-vegdist(comm,method=vegdist)
		nsp<-ncol(comm)
		if(is.euclidean(Dist)){
			commDist<-cmdscale(Dist,k=nsp-1)$points
			for(i in 1:dim(partAbund)[3]){
				res[i]<-RV(commDist,cmdscale(vegdist(partAbund[,,i],method=vegdist),k=nsp-1)$points)
			}
		}else{
			commDist<-cmdscale(Dist,k=nsp-1,add=TRUE)$points
			for(i in 1:dim(partAbund)[3]){
				res[i]<-RV(commDist,cmdscale(vegdist(partAbund[,,i],method=vegdist),k=nsp-1,add=TRUE)$points)
			}
			print("A Cailliez correction was applied")
		}
	}

	if(!is.null(designdist)){
		Dist<-designdist(comm,method=designdist,...)
		nsp<-ncol(comm)
		if(is.euclidean(Dist)){
			commDist<-cmdscale(Dist,k=nsp-1)$points
			for(i in 1:dim(partAbund)[3]){
				res[i]<-RV(commDist,cmdscale(designdist(partAbund[,,i],method=designdist,...),k=nsp-1)$points)
			}
		}else{
			commDist<-cmdscale(Dist,k=nsp-1,add=TRUE)$points
			for(i in 1:dim(partAbund)[3]){
				res[i]<-RV(commDist,cmdscale(designdist(partAbund[,,i],method=designdist,...),k=nsp-1,add=TRUE)$points)
			}
			print("A Cailliez correction was applied")
		}
	}
	
	if(!is.null(transfo) & !is.null(vegdist)){
		if(transfo == "log" & vegdist == "altGower"){
			Dist<-vegdist(decostand(comm,method=transfo,logbase=logbase),method=vegdist)
			nsp<-ncol(comm)
			if(is.euclidean(Dist)){
				commDist<-cmdscale(Dist,k=nsp-1)$points
				for(i in 1:dim(partAbund)[3]){
				res[i]<-RV(commDist,cmdscale(vegdist(partAbund[,,i],method=vegdist,...),k=nsp-1)$points)
				}
			}else{
				commDist<-cmdscale(Dist,k=nsp-1,add=TRUE)$points
				for(i in 1:dim(partAbund)[3]){
					res[i]<-RV(commDist,cmdscale(vegdist(partAbund[,,i],method=vegdist,...),k=nsp-1,add=TRUE)$points)
				}
				print("A Cailliez correction was applied")
			}		
		}
	}
	
	return(res)
}
