# Filename: PreHap.R
# Version : $Id: PreHap.R,v 1.6 2004/09/04 21:22:40 sblay Exp $

# HapAssoc- Inference of trait-haplotype associations in the presence of uncertain phase
# Copyright (C) 2003  K.Burkett, B.McNeney, J.Graham

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

########################################################################

pre.hapassoc <- function(dat,numSNPs,maxMissingGenos=1, pooling.tol=0.05, 
                        zero.tol=1/(2*nrow(dat)*10),method="default"){

  if(method=="PHASE")
    haplos.list<-RecodeHaplosPHASE(dat,numSNPs,maxMissingGenos)
  else
    haplos.list<-RecodeHaplos(dat,numSNPs,maxMissingGenos)

  haplotest<-FALSE; ID.check<-rep(FALSE,length(haplos.list$ID))

  # Starting matrices, some rows/columns will be deleted if there
  # are missing haplotypes

  newhaploDM <- haplos.list$haploDM
  newnonHaploDM <- haplos.list$nonHaploDM
  nonHaploDMnames<-names(haplos.list$nonHaploDM)
  newhaploMat <- haplos.list$haploMat
  newID <- haplos.list$ID

  #Run the usual EM algorithm using just the haplo information to
  #get estimates of haplotype frequencies and initial weights
  if(method=="PHASE")
  {
    emres<-EMnullPHASE() #in phase case, emres just contains gamma
    newwt<-haplos.list$wt  #In the phase case, we already have initial estimates given
  }
  else
  {
    emres<-EMnull(haplos.list)
    newwt<-emres$wts
  }

  zero.ind<-emres$gamma<zero.tol #flag haplos with zero frequency
  initGamma<-emres$gamma[!zero.ind]
  haplos.names<-names(initGamma)<-names(emres$gamma[!zero.ind])
  zeroFreqHaplos<-names(emres$gamma[zero.ind])

  if(sum(zero.ind)>0) { #then non-existent haplos need to be removed
    haplotest<-TRUE
    newhaploDM<-newhaploDM[,!zero.ind]

    #We only want rows that sum to two - others must have involved
    #haplotypes with estimated frequency of zero
    finalMatInd<- (rowSums(newhaploDM) == 2)

    newhaploDM<-newhaploDM[finalMatInd,]
    newhaploMat<-newhaploMat[finalMatInd,]
    newnonHaploDM<-data.frame(newnonHaploDM[finalMatInd,])
    names(newnonHaploDM)<-nonHaploDMnames
    newwt<-newwt[finalMatInd]
    newID<-newID[finalMatInd]
      
    # Re-calculate weights -- 2 loops!! there must be a better way
    uniqueIDs<-unique(newID)
    IDsum<-rep(0,length(uniqueIDs))
    for(i in 1:length(uniqueIDs)) {
      IDsum[i]<-sum(newwt[newID==uniqueIDs[i]])
    }
    for(i in 1:length(newwt)) {
      newwt[i] <- newwt[i]/IDsum[uniqueIDs==newID[i]]
    }
  }

  pooling.ind<-initGamma<pooling.tol #flag rare haplos 
  pooled.haplos<-"no pooled haplos"
  if(sum(pooling.ind)>1) { #then pooling to be done *in design matrix only*
    pooled.haplos<-haplos.names[pooling.ind]
    pooledDMcol<-rowSums(newhaploDM[,pooling.ind])
    newhaploDM<-data.frame(newhaploDM[,!pooling.ind],pooled=pooledDMcol)
  }
  
# The following code converts the columns of newnonHaploDM 
# to the right types to fix the effect of the previous 
# matrix conversions:                         -Sigal
j<-NULL
for(i in 1:length(dat)) j<-c(j, class(dat[[i]]))
for(i in 1:length(newnonHaploDM)) {
 if (is.numeric(dat[[i]])) 
 newnonHaploDM[[i]] <- as.numeric(as.character(newnonHaploDM[[i]]))
 else if(is.character(dat[[i]]))
 newnonHaploDM[[i]] <- as.character(newnonHaploDM[[i]])
}


  return(list(nonHaploDM=newnonHaploDM, haploDM=newhaploDM,
              haploMat=newhaploMat, wt=newwt, ID=newID, 
              haplotest=haplotest, initGamma=initGamma,
              zeroFreqHaplos=zeroFreqHaplos,pooledHaplos=pooled.haplos))
}


## Other functions called in CheckHaplos

########################################################################


isMissing <- function(test.vec){
  
  flag <- vector(length=max(test.vec))
  names(flag) <- c(1:length(flag))

  for (i in 1:length(flag)){if (sum(test.vec==i)==0){flag[i] <-TRUE} }
  return(flag)
}
  
isIn <- function(vec,vecElements) {
  #Return a vector of same length as vec. Element i of returned vector is
  #FALSE if vec[i] is not in the vector vecElements and TRUE if it is.
  flag<-rep(FALSE,length(vec))
  for(i in 1:length(vecElements)) {
    flag<-(flag | vec==vecElements[i])
  }
  return(flag)
}
  

EMnull<-function(haplos.list, gamma=FALSE, maxit=100, tol=1/(2*sum(haplos.list$wt)*100)){

 haploMat <- haplos.list$haploMat

 ID <- haplos.list$ID
 # N<-ID[length(ID)]
 N<-sum(haplos.list$wt)
 wts<-haplos.list$wt
 
 # Initial gamma values, if no gamma specified, calculate gamma values based
 # on augmented dataset.

 # We should avoid using the design matrix from an additive risk model to
 # help with haplotype frequency calculations. The following used to use
 # gamma <- (t(haplos)%*%wts)/(2*N)   where haplos is the additive model
 # design matrix as a computational trick to get haplotype frequencies.
 # Now use the tapply function and the haplotypes in haploMat to sum wts .

 if (gamma==FALSE){
    allHaps<-c(haploMat[,1],haploMat[,2])
    allWts<-c(wts,wts)
    gamma<-tapply(allWts,allHaps,sum)/(2*N)
 }

 gammadiff<-1
 it<-1
 num.prob<-vector(length=nrow(haploMat))
 
 # The EM loop

 while ( (it<maxit) && (gammadiff>tol) ){
   
        # multiplicative constant = 2 for heterozyg 1 for homozyg
        haplo.probs<-rep(1,nrow(haploMat))+isMultiHetero(haplos.list)
        haplo.probs <- haplo.probs*gamma[haploMat[,1]]*gamma[haploMat[,2]]
        num.prob<-haplo.probs

	# E step: Calculate the weights for everyone
	# Use the ID to determine the number of pseudo-individuals in the
	# denominator probability

	for (i in 1:nrow(haploMat)){               
		pseudo.index<-ID==ID[i]
                wts[i] <- num.prob[i]/sum(num.prob[pseudo.index])
        }

	# M step: Find new estimates using weighted haplotype counts
        allWts<-c(wts,wts)
	gammaNew <- tapply(allWts,allHaps,sum)/(2*N)
   	gammadiff<-max(abs(gamma-gammaNew), na.rm=TRUE)#maximum diff
   	gamma<-gammaNew

        it<-it+1
 }
 if(gammadiff>tol) 
   warning(paste("no convergence in EMnull after ",maxit,"iterations\n"))
 
 results <- list(gamma=gamma, wts=wts)

 return(results)

}

isMultiHetero<-function(haplos.list) {
  return(as.numeric(haplos.list$haploMat[,1] != haplos.list$haploMat[,2]))
}

############################################################
EMnullPHASE<-function()
{
	f=file("temp.out_freqs")
	rawdat=scan(file=f,list(x=""))

	#A little bit of subtlety.  One would expect that the length required for these
	#vectors would be 2^numSNPs.  If we did not screen out some rows in
	#handleMissingsPHASE, this would be true, however, because we do typically
	#screen out some individuals, this changes.  Why?  Because some haplotypes
	#with really low frequency (ie essentially 0) would only have appeared in
	#the data from pseudo-individuals reconstructed from rows with missing data.  As
	#we typically eliminate these rows with a lot of missing data, these pseudo-
	#individuals are not recreated, so PHASE may have less than 2^numSNPs haplotypes.
	#(divide by 4 since 4 columns and -1 since first row is header)
	gammanames=rep("h",length(rawdat$x)/4-1)
	gammavals=rep(0,length(rawdat$x)/4-1)

	i=5 #the first 4 entries are just the column headers in the PHASE *.out_freqs file
	while(!is.na(rawdat$x[i]))
	{
		idx=as.numeric(rawdat$x[i])
		i=i+1
		gammanames[idx]=paste(gammanames[idx],as.character(rawdat$x[i]),sep="")
		i=i+1
		gammavals[idx]=as.numeric(rawdat$x[i])
		i=i+1
		i=i+1 #skip over the 4th entry in the file, which is the std.err.
	}
	names(gammavals)=gammanames
	return(list(gamma=gammavals))
}
