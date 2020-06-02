regquadppp <-
function(data,ndivision=2,bottom=TRUE,left=TRUE){
### 
### Description:
### 
### This function builds a series of regularly distributed quadrats separating the full sampling area.
### 
### Arguments:
### 
### data : A ppp object.
### ndivision : Vector. The number of divisions used to split the sampling area
### bottom : Logical. Whether a tree that is directly on a horizontal quadrat line goes to the top or bottom quadrat. Default is TRUE (A tree will go in the bottom quadrat).
### left : Logical. Whether a tree that is directly on a vertical quadrat line goes to the left or right quadrat. Default is TRUE (A tree will go in the left  quadrat).
### 
### Details:
###
### ndivision needs to be an integer larger than one. As many divisons are made in X and Y so that the proportion between a quadrat and the sampling area is constant.
###
### Values :
###
### A vector counting the number of individuals in each grid cell
### 
### F. Guillaume Blanchet - June 2011, modified December 2011
################################################################################
	#CC# The general checks
	if(class(data)!="ppp"){
		stop("X needs to be a 'ppp' object")
	}
	
	#CC# Select column with X and Y coordinates and species
	xcoord<-data$x
	ycoord<-data$y
	
	#CC# find the extent in X and Y of the sampled area
	xrange<-data$window$xrange
	yrange<-data$window$yrange
	
	#CC# Correct for the number of division
	ndivision<-ndivision+2
	
	#CC# Construct a vector where quadrates numbers are assigned to each individuals
	nndivision<-length(ndivision)
	quad<-matrix(NA,ncol=nndivision,nrow=length(xcoord))
	
	#CC# Loop for each ndivision
	for(i in 1:nndivision){
		xquad<-seq(xrange[1],xrange[2],length=ndivision[i])
		yquad<-seq(yrange[1],yrange[2],length=ndivision[i])
		
		nside<-length(xquad)-1
		
		if(bottom & left){
			xquad[1]<-xquad[1]-0.01
			yquad[1]<-yquad[1]-0.01
		}
		if(!bottom & left){
			xquad[1]<-xquad[1]-0.01
			yquad[length(yquad)]<-yquad[length(yquad)]+0.01
		}
		if(bottom & !left){
			xquad[length(yquad)]<-xquad[length(yquad)]+0.01
			yquad[1]<-yquad[1]-0.01
		}
		if(!bottom & !left){
			xquad[length(yquad)]<-xquad[length(yquad)]+0.01
			yquad[length(yquad)]<-yquad[length(yquad)]+0.01
		}
		
		#CC# Counter
		quadNum<-1
		#CC# Loop to find individuals X coordinates
		for(j in 1:nside){
			#CC# Loop to find individuals Y coordinates
			for(k in 1:nside){
				#CC# Associate individuals with quadrates
				if(bottom & left){
					quad[which(xcoord > xquad[j] & xcoord <= xquad[j+1] & ycoord > yquad[k] & ycoord <= yquad[k+1]),i]<-quadNum
				}
				if(!bottom & left){
					quad[which(xcoord > xquad[j] & xcoord <= xquad[j+1] & ycoord >= yquad[k] & ycoord < yquad[k+1]),i]<-quadNum
				}
				if(bottom & !left){
					quad[which(xcoord >= xquad[j] & xcoord < xquad[j+1] & ycoord > yquad[k] & ycoord <= yquad[k+1]),i]<-quadNum
				}
				if(!bottom & !left){
					quad[which(xcoord >= xquad[j] & xcoord < xquad[j+1] & ycoord >= yquad[k] & ycoord < yquad[k+1]),i]<-quadNum
				}
				quadNum<-quadNum+1
			}
		}
		#CC# Reinitialize the counter
		quadNum<-1
		
	}
	
	return(quad)
}
