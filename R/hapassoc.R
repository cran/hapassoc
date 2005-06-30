# Filename: hapassoc.R
# Version : $Id: hapassoc.R,v 1.21 2005/06/30 21:40:13 sblay Exp $

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

hapassoc<-function(form, haplos.list, baseline="missing", family=binomial(), 
             freq=FALSE, maxit=50, tol=0.001, ...){

 environment(form)<-environment()#set envir of formula to envir w/i hapassoc function
 column.subset <- colnames(haplos.list$haploDM)!=baseline
 hdat <- cbind(haplos.list$nonHaploDM, haplos.list$haploDM[,column.subset])
 # We used to use model.response to extract the response variable, but
 # this lead to problems in calculating the residuals in the pYgivenX function.
 # It is better to fit the model with the glm function below and then 
 # extract the response from the fitted model object - Brad Nov.2/04
 # response<-model.response(model.frame(form,data=hdat))  
 colnames(hdat)<- c(colnames(haplos.list$nonHaploDM),
                    colnames(haplos.list$haploDM[,column.subset]))
 ID <- haplos.list$ID
 # N<-ID[length(ID)]
 N<-sum(haplos.list$wt)
 wts<-haplos.list$wt
 
 # Get the haplotype columns
 haplos<-haplos.list$haploDM
 haploMat<-haplos.list$haploMat
 allHaps<-c(haploMat[,1],haploMat[,2]) #Needed later in hapassoc loop for wt calcs
 haplos.names<-names(haplos.list$initFreq)

 # Initial freq values, if no freq specified use initFreq

 if (freq!=FALSE) {
   names(freq)<-haplos.names
 } else {
   freq<-haplos.list$initFreq
 }

 # Initial beta values calculated from augmented data set
 # The ... in the following call to glm allows user to pass other args to glm

 regr<-glm(form, family=family, data=hdat, weights=wts, ...) 
 response<-regr$y #Change added Nov.2/04 to extract response from fitted model
 beta<-regr$coef
 fits<-regr$fitted.values
 betadiff<-1
 it<-1
 num.prob<-vector(length=nrow(hdat))
 
 # The hapassoc loop

 while ( (it<maxit) && (betadiff>tol) ){
   
        # Vector of P(Y|X)*P(X) probability 
	# If the person is unaffected, P(Y=0) is 1-fitted value
	# otherwise it is the fitted value 
   
        # Multiplicative const for haplo probs: 1 for homo, 2 for het
        haplo.probs<-rep(1,nrow(haplos))+isMultiHetero(haplos.list)
        haplo.probs <- haplo.probs*freq[haploMat[,1]]*freq[haploMat[,2]]

	phi<-mlPhi(regr) #Compute ML estimate of dispersion param
        if(is.null(phi)) { #no converergence in ml estimate of phi
           break() #hapassoc will throw a warning of non-convergence
        }
        num.prob<-pYgivenX(response,fits,phi,family)*haplo.probs

	# E step: Calculate the weights for everyone
	# Use the ID to determine the number of pseudo-individuals in the 
	# denominator probability

        wts<-vector("double", length(wts))
        wts<-.C("getWts", as.integer(nrow(hdat)), as.integer(ID),
                 wts = as.double(wts), as.double(as.vector(num.prob)), 
                 PACKAGE="hapassoc")
        wts <- wts$wts

	# M step: Find new estimates using GLM and weighted haplotype counts

	# Find the new betas using old betas as starting value
        regr <- glm(form, family=family, data=hdat, weights=wts,
                    control=glm.control(epsilon=1e-08),start=beta)
        betaNew<-regr$coef
   	fits<-regr$fitted.values
	betadiff<-max( max(abs(beta-betaNew), na.rm=TRUE),
	 max(abs(beta-betaNew)/(0.1+abs(betaNew)), na.rm=TRUE) )   	
	beta<-betaNew

	# Find the new freqs, weighted sum of haplotypes
        freq <- tapply(c(wts,wts),allHaps,sum)/(2*N)

        it<-it+1
 }
 
 if(betadiff>tol) { #did not converge
    warning(paste("No convergence in hapassoc in",it,"iterations\n")) 
    ans<-list(converged=FALSE)
    class(ans)<-"hapassoc"
    return(ans)
 }


 
 #freq is currently an array which causes problems in the calculations below
 freq <- as.matrix(freq) 

 EMresults <- list(beta=beta, gamma=freq, fits=fits, wts=wts, 
                   glm.final.fit=regr,dispersion=phi, response=response)
 names(haplos.list)[names(haplos.list)=="freq"] <- "gamma"
 names(EMresults)[names(EMresults)== "freq"] <- "gamma"
 var.est <- EMvar(haplos.list, EMresults, family)

 ans<-list(it=it, beta=beta, freq=freq, fits=fits, wts=wts, ID=ID,
           var=var.est, dispersionML=phi, family=family, response=response,
           converged=TRUE)

 class(ans)<-"hapassoc"
 return(ans)

}

## Other functions called in hapassoc:

########################################################################

mlPhi<-function(myfit,maxit=30,tol=1e-5) {
  if(myfit$family$family=="binomial" || myfit$family$family=="poisson") {
    return(1) #dispersion set to one
  }
  if(myfit$family$family=="gaussian") {
    return(summary(myfit)$deviance/sum(myfit$prior.weights))
  }
  if(myfit$family$family=="Gamma") { #Need to do Newton-Raphson
    return(mlPhiGamma(myfit,maxit,tol))
  }
  #Else, we don't support the family passed in
  stop(paste("ML estimation of dispersion not supported for family",
                myfit$family$family))
}

########################################################################

mlPhiGamma<-function(myfit,maxit,tol) {
  
  # Need to do Newton-Raphson, use moment estimate of phi to get
  # starting value for N-R 

  dev<-myfit$deviance
  n<-sum(myfit$prior.weights)
  phiMoment<-summary(myfit)$dispersion #moment estimate of phi
  ## Looking for root of score equation as a function of x=1/phi
  ## Function f is the score equation and fp is its derivative
  f<-function(x,n,dev) { return(dev+2*n*(digamma(x)-log(x))) }
  fp<-function(x,n) { return(2*n*(trigamma(x)-1/x)) }
  diff<-1; i<-1  # initialization
  xold<-1/phiMoment #starting value for N-R 
  while(i<=maxit && diff>tol) {
    xnew<-xold - f(xold,n,dev)/fp(xold,n)
    if(is.na(xnew)) {
       warning("No convergence for ML estimate of Gamma scale parameter")
       return(NULL)
    }
    diff<-abs(xnew-xold); xold<-xnew; i<-i+1
  }
  if(i>maxit) 
    warning("No convergence for ML estimate of Gamma scale parameter")
  return(1/xnew)
}

########################################################################

pYgivenX<-function(y,mu,phi,family){

  #Calculate P(Y|X) to be used in the weights. We use P(Y|X) in 
  #expressions like P(Y_i|X_i^j)P(X_I^j)/sum_k{ P(Y_i|X_i^k)P(X_i^k) }
  #so factors that are the same for all X_i^j can be ignored.
  #The deviance resids can be used to get P(Y|X) up to constants:
  #Assuming weights of 1 
  #dev.resid = 2*phi*(l(y,phi) - l(mu,phi)) where l(mu,phi) is the log-lik
  #using mean mu and l(y,phi) is the log-likelihood in the saturated
  #model using y as the mean.
  #So   exp(-dev.resid/(2*phi)) = P(y|x)*exp(-l(y,phi))
  #Call the dev.resid function with a weight of 1 

  return(exp(-1*family$dev.resid(y,mu,wt=1)/(2*phi)))
} 


## EMvar Functions
 
########################################################################

EMvar<-function(haplos.list, results, family)  {
 # Get the results of the last iteration
 beta<-results$beta[!is.na(results$beta)]
 betanames<-names(beta)
 num.beta<-length(beta)

 gammanames<-rownames(results$gamma)
 gamma<-as.vector(results$gamma)
 num.gamma<-length(gamma)

 weights<-results$wt
 ID <- haplos.list$ID

 final.regr <- results$glm.final.fit

 # Set up the Y vectors and P vector of fitted values
 # one for the augmented data set, one for just the missing
 
 missing<-(weights<1)

 yi.full<-results$response
 fits.full<-results$fits
 
 yi.mis<-yi.full[missing]
 fits.mis<-fits.full[missing]

 # Set up the matrices needed to find the variance:

   ## Design matrices

   Xreg <- model.matrix(final.regr)
   colnames(Xreg)<-betanames
   Xreg.mis<-Xreg[missing,]

   ## Haplotype matrix
   #Instead of assuming the glm model is additive in the haplotypes, 
   #make the Xgen matrix from scratch here. Used to be
   #   Xgen<-as.matrix(haplos.list$haploDM)
   haploMat<-haplos.list$haploMat
   #each application of outer in next line returns a matrix of T/F with
   # (i,j)th element true if haploMat[i,k]==gammanames[j]; k=1,2
   #Adding these T/F matrices coerces to 1/0 s first before adding
   Xgen<-outer(haploMat[,1],gammanames,"==")+outer(haploMat[,2],gammanames,"==")

   Xgen.mis <- as.matrix(Xgen[missing,])
 
   ## Weight matrices

   W.full<-diag(weights)
   W.mis<-diag(weights[missing])

   ## likelihood derivative matrices

   H.mis<-diag(yi.mis-fits.mis)/results$dispersion

   ## Haplotype frequency and derivative matrices

   ones.mat<-matrix(1,nrow=(num.gamma-1), ncol=(num.gamma-1))
   num.haplo<-vector(length=num.gamma)
   for (i in 1:num.gamma){
	num.haplo[i]<-sum(weights*Xgen[,i])}
   der.gamma<-1/gamma
   der2.gamma<-1/gamma^2

   G1<-diag(der.gamma)
   G<-G1[,1:(num.gamma-1)]
   G[num.gamma,]<-G1[num.gamma, num.gamma]

   ## Block diagonal matrix

   ID.mis <- ID[missing]
   i <- 1
   B<-matrix(0,nrow=sum(missing), ncol=sum(missing) )

   while (i<length(ID.mis)){
	pseudo.index<-ID.mis==ID.mis[i]
	ones.vec<-rep(1, sum(pseudo.index)*sum(pseudo.index))
	block<-matrix(ones.vec, nrow=sum(pseudo.index), ncol=sum(pseudo.index))
	block.size<-sum(pseudo.index)
	B[i:(i+block.size-1), i:(i+block.size-1)]<-block
	i<-i+sum(pseudo.index)
   }

 # Calculate the complete data expected information (Ic)

   ## Top block

   Ic.reg<-solve(summary(final.regr)$cov.unscaled)/results$dispersion

   ## Middle block if dispersion was estimated -- Ic.phi requires Sphi

   Sphi<-SPhi(final.regr,results$dispersion) #will be NULL if phi not est'd
   if(!is.null(Sphi)) { 
     Ic.phi<-IPhi(final.regr,Sphi,results$dispersion)
   } else { 
     Ic.phi<-NULL 
   }

   ## Bottom block

   Ic.cov<- diag(num.haplo[1:(num.gamma-1)]*der2.gamma[1:(num.gamma-1)])+der2.gamma[num.gamma]*num.haplo[num.gamma]*ones.mat

   ## Full block diagonal matrix (combine top, middle and bottom blocks)

   if(!is.null(Ic.phi)) {
     Ic<-matrix(0,nrow=(nrow(Ic.reg)+1+nrow(Ic.cov)),
                ncol=(nrow(Ic.reg)+1+nrow(Ic.cov))) #+1 is for phi
     Ic[1:num.beta,1:num.beta]<-Ic.reg
     Ic[(num.beta+1),(num.beta+1)]<-Ic.phi
     Ic[(num.beta+2):(num.beta+num.gamma),
        (num.beta+2):(num.beta+num.gamma)]<-Ic.cov
   } else {
     Ic<-matrix(0,nrow=(nrow(Ic.reg)+nrow(Ic.cov)),
                ncol=(nrow(Ic.reg)+nrow(Ic.cov)))
     Ic[1:num.beta,1:num.beta]<-Ic.reg
     Ic[(num.beta+1):(num.beta+num.gamma-1),
        (num.beta+1):(num.beta+num.gamma-1)]<-Ic.cov
   }
     

 # Calculate the missing data information (Imis), only over missing

   ## S_beta
 
   Sbeta<-H.mis%*%Xreg.mis

   ## S_phi

   if(!is.null(Sphi)) {
     Sphi<-Sphi[missing]
   }
 
   ## S_gamma
 
   Sgamma<-Xgen.mis%*%G

   ## Score matrix

   S<-cbind(Sbeta,Sphi,Sgamma)

   ## Calculate Imis from the other matrices

   Imis<-t(S)%*%(W.mis-W.mis%*%B%*%W.mis)%*%S  

 # Calculate the information and var/cov matrix

 Info<-Ic-Imis
 varEM<-solve(Info)

 gammanames <- paste("f", gammanames, sep="")
 if(!is.null(Sphi)) {
   colnames(varEM)<-c(betanames, "phi", gammanames[1:(num.gamma-1)])
 }else {
   colnames(varEM)<-c(betanames, gammanames[1:(num.gamma-1)])
 }
 rownames(varEM)<-colnames(varEM)


 return(varEM)

}


## Other functions called in EMvar:

########################################################################

SPhi<-function(myfit,phi) {
  if(myfit$family$family=="binomial" || myfit$family$family=="poisson") {
    return(NULL) #dispersion set to one so score not defined
  }
  if(myfit$family$family=="gaussian") {
    return(SPhiGaussian(myfit,phi))
  }
  if(myfit$family$family=="Gamma") { 
    return(SPhiGamma(myfit,phi))
  }
}

########################################################################

SPhiGaussian<-function(myfit,phi) {
  return( (myfit$y-myfit$fitted)^2/(2*phi^2) - 1/(2*phi) )
}

########################################################################

SPhiGamma<-function(myfit,phi) {
  x<-1/phi
  y<-myfit$y
  mu<-myfit$fitted.values
  return( (digamma(x)-log(x)+(y-mu)/mu-log(y/mu))*x^2 )
}

########################################################################

IPhi<-function(myfit,score,phi) {
  if(myfit$family$family=="binomial" || myfit$family$family=="poisson") {
    return(NULL) #dispersion set to one so phi not estimated
  }
  if(myfit$family$family=="gaussian") {
    return(IPhiGaussian(myfit,score,phi))
  }
  if(myfit$family$family=="Gamma") {
    return(IPhiGamma(myfit,score,phi))
  }
}

########################################################################

IPhiGaussian<-function(myfit,score,phi) {
  wts<-myfit$prior.weights
  n<-sum(wts)
  return( 2*(t(score)%*%wts)/phi + n/(2*phi^2) )
}

########################################################################

IPhiGamma<-function(myfit,score,phi) {
  x<-1/phi
  wts<-myfit$prior.weights
  n<-sum(wts)
  return( 2*x*(t(score)%*%wts) + n*x^4*(trigamma(x) - 1/x) )
}

  
