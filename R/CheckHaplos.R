# Filename: CheckHaplos.R
# Version : $Id: CheckHaplos.R,v 1.4 2004/01/28 02:24:41 mcneney Exp $

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

CheckHaplos <- function(haplos.list, numSNPs, pooling.tol=0.05, 
                        zero.tol=1/(2*sum(haplos.list$wt)*10)){

  haplo.vec <- makeHaploLab(0:(2^numSNPs-1), numSNPs)
  names.haplos<-paste("n",haplo.vec,sep="")
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
  emres<-EMnull(haplos.list)
  newwt<-emres$wts
  
  zero.ind<-emres$gamma<zero.tol #flag haplos with zero frequency
  initGamma<-emres$gamma[!zero.ind]
  names(initGamma)<-names(emres$gamma[!zero.ind])
  zeroFreqHaplos<-names(emres$gamma[zero.ind])

  if(sum(zero.ind)>0) { #then non-existent haplos need to be removed
    haplotest<-TRUE
    newhaploDM<-newhaploDM[,!zero.ind]
    haplo.vec<-haplo.vec[!zero.ind]
  }
  pooling.ind<-initGamma<pooling.tol #flag rare haplos 
  pooled.haplos<-"no pooled haplos"
  if(sum(pooling.ind)>1) { #then there is pooling to be done
    pooled.haplos<-haplo.vec[pooling.ind]
    #Recode pooled haplos in haploMat as "pooled"
    newhaploMat[isIn(newhaploMat[,1],pooled.haplos),1]<-"pooled"
    newhaploMat[isIn(newhaploMat[,2],pooled.haplos),2]<-"pooled"
    #Sum pooled columns of haploDM
    pooledDMcol<-apply(newhaploDM[,pooling.ind],1,sum)
    newhaploDM<-data.frame(newhaploDM[,!pooling.ind],pooled=pooledDMcol)
    pooled.freq<-sum(initGamma[pooling.ind])
    temNames<-names(initGamma[!pooling.ind])
    initGamma<-c(initGamma[!pooling.ind],pooled.freq)
    names(initGamma)<-c(temNames,"pooled")
  }
    
  if(sum(zero.ind)>0 | sum(pooling.ind)>1) {
    #Some rows of haploDM may now be duplicates. Find the
    #ones that aren't as these will be in the final matrices output by
    #CheckHaplos
    uniqueIDs<-unique(newID)
    finalMatInd<-NULL
    for (i in 1:length(uniqueIDs) ){
      IDind<-(newID==uniqueIDs[i])
      duplicatesInd<-duplicated(newhaploDM[IDind,]) #built-in R func duplicated
      finalMatInd<-c(finalMatInd,!duplicatesInd)
    }
    #Also, we only want rows that sum to two - others must have involved
    #haplotypes with estimated frequency of zero
    finalMatInd<-( finalMatInd & (rowSums(newhaploDM) == 2) )

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

    ID.check <- isMissing(newID) # check if an ID has been completely deleted
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

 haplos <- haplos.list$haploDM
 haplos.names <- colnames(haplos)

 ID <- haplos.list$ID
 # N<-ID[length(ID)]
 N<-sum(haplos.list$wt)
 wts<-haplos.list$wt
 
 # Initial gamma values, if no gamma specified, calculate gamma values based
 # on augmented dataset.
 num<-nrow(haplos) 

 if (gamma!=FALSE) { names(gamma)<-haplos.names }
 if (gamma==FALSE){
	gamma <- (t(haplos)%*%wts)/(2*N)
 	names(gamma)<-haplos.names
 }

 gammadiff<-1
 it<-1
 num.prob<-vector(length=nrow(haplos))
 
 # The EM loop

 while ( (it<maxit) && (gammadiff>tol) ){
   
        # multiplicative constant = 2 for heterozyg 1 for homozyg
        haplo.probs<-rep(1,nrow(haplos))+as.numeric(isMultiHetero(haplos.list))
	for (i in 1:length(gamma)){
          haplo.probs <- haplo.probs*gamma[i]^haplos[,i]
        }
        num.prob<-haplo.probs

	# E step: Calculate the weights for everyone
	# Use the ID to determine the number of pseudo-individuals in the 
	# denominator probability

	for (i in 1:nrow(haplos)){               
		pseudo.index<-ID==ID[i]
                wts[i] <- num.prob[i]/sum(num.prob[pseudo.index])
        }

	# M step: Find new estimates using weighted haplotype counts
	gammaNew <- (t(haplos)%*%wts)/(2*N)
   	gammadiff<-max(abs(gamma-gammaNew), na.rm=TRUE)#maximum diff
   	gamma<-gammaNew

        it<-it+1
 }
 if(it==maxit) 
   warning(paste("no convergence in EMnull after ",maxit,"iterations\n"))
 gammaNames<-dimnames(gamma)[[1]]
 gamma<-as.vector(gamma)
 names(gamma)<-gammaNames
 
 results <- list(gamma=gamma, wts=wts)

 return(results)

}

isMultiHetero<-function(haplos.list) {
  return(as.numeric(haplos.list$haploMat[,1] != haplos.list$haploMat[,2]))
}
