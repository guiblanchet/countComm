commSimul<-function(nquad,nsp, nind,patchrad=0.2,SAD="lognormal",sdlog=5){
###
### Function to simulate a community following the procedure proposed in Blanchet et al. (submitted).
###
### This function simulates a community in a square of unit size. 
###
### Arguments :
###
### nquad : Number of quadrates in which to divide the sampling area. This should be equal to any integer squared. (See Details).
### nsp : Integer defining the number of species to be simulated.
### nind : Integer defining the maximum number of individuals for each species. Note that if SAD = "uniform", this value will be the same for each species. (See details).
### patchrad : A numeric value defining the maximum radius of a patch for a group of individuals of a single species. (See details).
### SAD : Character string defining the type species-abundance distribution for the community. Either "lognormal", "bstick", "uniform"
### sdlog : This argument is active only if SAD= "lognormal". It defines the standard deviation of the lognormal distribution.
### 
### Details
###
### The number of quadrates in the sampling area (a square of unit size) is defined by dividing the area horizontally and vertically so that each quadrate is proportionally equal to the size of the sampling area. For this reason, nquad should be an integer that can be obtain from any solution of x^2. If it is not the case, a value will be calculated that approximate the desired value by rounding the square root nquad and putting it to the power of 2.
###
### Depending on the way individuals are distributed in the sampling area (spatagg), the number of individuals per species defined by 'nind' may varies because the method use to sample individual relies on random point patterns and the value given is the average of the random point process. 
### 
### Note that when defining the patch radius through patchrad, the patch radius is sampled from a uniform distribution ranging from 0.01 to patchrad. Also, since all simulations are carried out in a square unit size, defining patchrad with a value larger than 1 is ill advised. The argument patchrad defines patch radius for individuals of each species simulated so that the patch size associated to different species within the same community can vary within the predefine range.
###
### F. Guillaume Blanchet - Juin 2014
##########################################################################################
	SAD<-match.arg(SAD)
	
	### General checks
	if(!any(SAD==c("lognormal","bstick","uniform"))){
		stop("'SAD' should be either 'lognormal', 'bstick', 'uniform'")
	}
	
	if(patchrad > 1){
		warning("'patchrad' is larger than 1")
	}
	
	ndiv<-sqrt(nquad)
	ndivRound<-round(ndiv)
	if(ndivRound^2!=ndiv^2){
		ndiv<-ndivRound
		warning("The number of quadrats define by 'nquad' will be different than the number of rows in the final community matrix")
	}

	comm<-matrix(0,nrow=ndiv^2,ncol=nsp)
	
	for(i in 1:nsp){
		### This loop is to ensure that all species have at least one individuals
		n<-0
		while(n==0){
			### Defining parameters of the distribution
			if(SAD=="lognormal"){
				probability<-dlnorm(1:nind,sdlog=sdlog)
				nindivi<-sample(nind,1,prob=probability)
			}
			if(SAD=="bstick"){
				probability<- bstick(nind)
				nindivi<-sample(nind,1,prob=probability)
			}
			if(SAD=="uniform"){
				nindivi<-nind
			}
			
			### Parameter for the Matern cluster process
			if(nindivi == 1){
				intensePoisson<-1
				meanIndCluster<-1
			}
			if(nindivi > 1){
				intensePoisson<-sample(nindivi,1)
				meanIndCluster<-round(nindivi/intensePoisson)
			}
			
			radius<-runif(1,0.01,patchrad)
			sp<-rMatClust(intensePoisson,radius,meanIndCluster)
			n<-sp$n
		}
		
		### Associate individuals to sites (quadrats)
		spquad<-regquadppp(sp,ndiv-1)
		indCount<-summary(as.factor(spquad),maxsum=sp$n+1)
		comm[as.numeric(names(indCount)),i]<-indCount
	}
	### Final results
	colnames(comm)<-paste("sp",1:nsp,sep="")
	rownames(comm)<-paste("quad",1:nquad,sep="")
	return(comm)
}