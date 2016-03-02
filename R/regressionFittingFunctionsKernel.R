#' Nearest Neighbors Regression
#'
#' This function is typically not directly used by the user; it is used inside  \code{\link{fitFlexCoDE}}
#'
#' @param kernelX matrix with covariates that will be used for training
#' @param responses matrix where each column is a response for the training data
#' @param extra list with one component named nNeighVec, which contains a vetor with different number of neighbors; the function will choose the best value among them
#'
#' @return object of the class NN containing information need to perform prediction on new points
#' @export
regressionFunction.NNKernel=function(kernelX,responses,extra=NULL)
{
  # Both x and responses are matrices
  n=dim(kernelX)[1]
  random=sample(1:n)
  nTrain=round(0.7*n)
  #kernelTrainTrain=kernelX[random[1:nTrain],random[1:nTrain]]
  responsesTrain=responses[random[1:nTrain],]
  kernelValidationTrain=kernelX[random[-c(1:nTrain)],random[1:nTrain],drop=FALSE]
  responsesValidation=responses[random[-c(1:nTrain)],,drop=FALSE]


  nNeighVec=extra$nn
  if(is.null(nNeighVec))
    nNeighVec=round(seq(1,dim(responsesTrain)[1],length.out = 100))


  error=matrix(NA,length(nNeighVec),dim(responsesTrain)[2])
  for(ii in 1:length(nNeighVec))
  {
    predictedValidation=t(apply(-kernelValidationTrain,1,function(xx) {
      nearest=sort(xx,index.return=T)$ix[1:nNeighVec[ii]]
      return(colMeans(responsesTrain[nearest,,drop=FALSE]))
    }))
    error[ii,]=colMeans((predictedValidation-responsesValidation)^2)
  }

  bestNN=apply(error,2,function(xx){
    nNeighVec[which.min(xx)]
  })
  rm(kernelValidationTrain)
  gc(verbose = FALSE)

  regressionObject=NULL
  regressionObject$bestNN=bestNN
  regressionObject$responsesTrain=responses
  class(regressionObject)="NNKernel"
  rm(responsesTrain)
  gc(verbose=FALSE)
  return(regressionObject)
}


#' Print function for object of the class NN
#'
#' This function is typically not directly used by the user; it is used inside  \code{\link{print.FlexCoDE}}, the print
#' method for the class FlexCoDE
#'
#' @param regressionObject of the class NN
#' @param bestI optimal number of expansion coefficients
#' @param nameCovariates name of the covariates
#'
#' @return prints characteristics of the regressions that were fitted
#'
print.NNKernel=function(regressionObject,bestI,nameCovariates)
{
  cat(paste("Number of neighbors chosen for each fitted regression:",paste(regressionObject$bestNN[1:bestI],collapse=", "),"\n"))
}


#' Print function for object of the class NN
#'
#' This function is typically not directly used by the user; it is used inside  \code{\link{print.FlexCoDE}}, the print
#' method for the class FlexCoDE
#'
#' @param regressionObject of the class NN
#' @param bestI optimal number of expansion coefficients
#' @param nameCovariates name of the covariates
#'
#' @return prints characteristics of the regressions that were fitted
#'
print.NNKernel=function(regressionObject,bestI,nameCovariates)
{
  cat(paste("Number of neighbors chosen for each fitted regression:",paste(regressionObject$bestNN[1:bestI],collapse=", "),"\n"))
}


#' Spectral Series Regression
#'
#' This function is typically not directly used by the user; it is used inside  \code{\link{fitFlexCoDE}}
#'
#' @param kernelX gram matrix that will be used for training
#' @param responses matrix where each column is a response for the training data
#' @param extra list with one component called nXMax and contains a single integer number that describes what is the maximum number of spectral basis functions with be used
#'
#' @return object of the class Series containing information needed to perform prediction on new points
#' @export
#'
regressionFunction.SeriesKernel=function(kernelX,responses,extra=NULL)
{
  # Both x and responses are matrices
  n=dim(kernelX)[1]
  random=sample(1:n)
  nTrain=round(0.7*n)
  kernelMatrix=kernelX[random[1:nTrain],random[1:nTrain]]
  responsesTrain=responses[random[1:nTrain],]
  kernelMatrixValidationTrain=kernelX[random[-c(1:nTrain)],random[1:nTrain]]
  responsesValidation=responses[random[-c(1:nTrain)],,drop=FALSE]



  nXBestEps=matrix(NA,ncol(responsesValidation))


  nAll = dim(kernelMatrix)[1]

  if(!is.null(extra$nXMax))
  {
    nXMax=min(nrow(kernelMatrix)-15,extra$nXMax)
  } else {
    nXMax=nrow(kernelMatrix)-15
  }
  p=10

  Omega=matrix(rnorm(nAll*(nXMax+p),0,1),nAll,nXMax+p)
  Z=kernelMatrix%*%Omega
  Y=kernelMatrix%*%Z
  Q=qr(x=Y)
  Q=qr.Q(Q)
  B=t(Q)%*%Z%*%solve(t(Q)%*%Omega)
  eigenB=eigen(B)
  lambda=eigenB$values
  U=Q%*%eigenB$vectors

  basisX=Re(sqrt(nAll)*U[,1:nXMax])
  eigenValues=Re((lambda/nAll)[1:nXMax])

  coefficientsMatrix=1/nAll*t(responsesTrain)%*%basisX

  m=nrow(kernelMatrixValidationTrain) # New

  basisXNew=kernelMatrixValidationTrain %*% basisX
  basisXNew=1/n*basisXNew*matrix(rep(1/eigenValues,m),m,ncol(basisX),byrow=T)

  errorsForEachReg=apply(as.matrix(1:nrow(coefficientsMatrix)),1,function(ii)
  {
    coeff=coefficientsMatrix[ii,]
    errors=apply(as.matrix(1:length(coeff)),1,function(kk)
    {
      betaHat=t(coeff[1:kk,drop=FALSE])
      predictedY=basisXNew[,1:kk,drop=F]%*%t(betaHat)
      return(mean((predictedY-responsesValidation[,ii])^2))
    })
    nXBest=(1:length(coeff))[which.min(errors)]
    bestError=min(errors)
    return(c(bestError,nXBest))
  })

  nXBestEps=errorsForEachReg[2,]


  rm(basisXNew)
  gc(verbose=FALSE)

  p=10

  nAll=nrow(kernelX)
  Omega=matrix(rnorm(nAll*(nXMax+p),0,1),nAll,nXMax+p)
  Z=kernelX%*%Omega
  Y=kernelX%*%Z
  Q=qr(x=Y)
  Q=qr.Q(Q)
  B=t(Q)%*%Z%*%solve(t(Q)%*%Omega)
  eigenB=eigen(B)
  lambda=eigenB$values
  U=Q%*%eigenB$vectors

  basisX=Re(sqrt(nAll)*U[,1:nXMax])
  eigenValues=Re((lambda/nAll)[1:nXMax])

  coefficientsMatrix=1/nAll*t(responses)%*%basisX


  rm(kernelX,kernelMatrixValidationTrain,kernelMatrix)
  gc(verbose = FALSE)

  regressionObject=NULL
  regressionObject$nXMax=nXMax
  regressionObject$bestNX=nXBestEps # vector with best cutoffs
  regressionObject$coefficientsMatrix=coefficientsMatrix
  regressionObject$basisX=basisX
  regressionObject$eigenValues=eigenValues
  class(regressionObject)="SeriesKernel"
  rm(responsesTrain,coefficientsMatrix)
  gc(verbose=FALSE)
  return(regressionObject)
}


#' Print function for object of the class Series
#'
#' This function is typically not directly used by the user; it is used inside  \code{\link{print.FlexCoDE}}, the print
#' method for the class FlexCoDE
#'
#' @param regressionObject of the class Series
#' @param bestI optimal number of expansion coefficients
#' @param nameCovariates name of the covariates
#'
#' @return prints characteristics of the regressions that were fitted
#'
print.SeriesKernel=function(regressionObject,bestI,nameCovariates)
{
  cat(paste("Number of expansion coefficients chosen for each fitted regression:",paste(regressionObject$bestNX[1:bestI],collapse = ", "),"\n"))
  cat("\n")
  cat(paste("Best epsilon for spectral decomposition:",regressionObject$bestEps,"\n"))

}


#' Title test
#'
#' @param kernelX
#' @param responses
#' @param extra
#'
#' @return r
#' @export
#'
#' @import kernlab
#'
regressionFunction.SDMKernel=function(kernelX,responses,extra=NULL)
{
  # Both x and responses are matrices

  C=extra$C
  if(is.null(C))
    C=exp(seq(log(0.05),log(0.5),length.out = 5))

  eps=extra$eps
  if(is.null(eps))
    eps=seq(0.1,1,length.out = 5)


  # fittedReg=list()
  # fittedReg[[1]]=1
  # errors=matrix(NA,length(C),length(eps))
  # for(i in 2:ncol(responses))
  # {
  #   for(j in 1:length(C))
  #   {
  #     for(k in 1:length(eps))
  #     {
  #       fit=try(kernlab::ksvm(x=kernelX,y=responses[,i],type="eps-svr",kernel="matrix",cross=2,epsilon=eps[k],C=C[j]),silent = TRUE)
  #       if(class(fit)=="try-error")
  #         next;
  #       errors[j,k]=fit@error
  #     }
  #   }
  #   which=which(errors==min(errors),arr.ind = TRUE)
  #   fittedReg[[i]]=kernlab::ksvm(x=kernelX,y=responses[,i],type="eps-svr",kernel="matrix",cross=2,epsilon=eps[which[2]],C=C[which[1]])
  # }
  #


  nCores=extra$nCores
  if(is.null(nCores))
    nCores=1

  cl <- parallel::makeCluster(nCores)
  doParallel::registerDoParallel(cl)


  fittedReg <- foreach(i=2:ncol(responses)) %dopar% {
    errors=matrix(NA,length(C),length(eps))
    for(j in 1:length(C))
    {
      for(k in 1:length(eps))
      {
        fit=try(kernlab::ksvm(x=kernelX,y=responses[,i],type="eps-svr",kernel="matrix",cross=2,epsilon=eps[k],C=C[j]),silent = TRUE)
        if(class(fit)=="try-error")
          next;
        errors[j,k]=fit@error
      }
    }
    which=which(errors==min(errors),arr.ind = TRUE)
    return(kernlab::ksvm(x=kernelX,y=responses[,i],type="eps-svr",kernel="matrix",cross=2,epsilon=eps[which[2]],C=C[which[1]]))
  }

  parallel::stopCluster(cl)

  fittedReg=append(1,fittedReg)

  regressionObject=NULL
  regressionObject$fittedReg=fittedReg
  class(regressionObject)="SDMKernel"
  rm(responses)
  gc(verbose=FALSE)
  return(regressionObject)
}



#' Print function for object of the class SDMKernel
#'
#' This function is typically not directly used by the user; it is used inside  \code{\link{print.FlexCoDE}}, the print
#' method for the class FlexCoDE
#'
#' @param regressionObject of the class SDMKernel
#' @param bestI optimal number of expansion coefficients
#' @param nameCovariates name of the covariates
#'
#' @return prints characteristics of the regressions that were fitted
#'
print.SDMKernel=function(regressionObject,bestI,nameCovariates)
{
  cat(paste("Number of expansion coefficients chosen for each fitted regression:",paste(regressionObject$bestNX[1:bestI],collapse = ", "),"\n"))
  cat("\n")
  cat(paste("Best epsilon for spectral decomposition:",regressionObject$bestEps,"\n"))
}
