#' Predict NN regression
#'
#' This function is typically not directly used by the user; it is used inside  \code{\link{fitFlexCoDE}}
#'
#' @param object object of the class NN
#' @param xNew matrix with covariates where prediction will be calculated
#' @param maxTerms maximum number of expansion coefficients
#'
#' @return returns matrix where element (i,j) contains the estimate of the j-th expansion coefficient for the j-th sample
#' @export
#'
predict.NNKernel=function(object,kernelNewTrain,maxTerms=NULL)
{
  if(class(object)!="NNKernel")
    stop("Object has wrong class, should be NNKernel")


  if(!is.null(maxTerms))
  {
    maxTerms=min(maxTerms,length(object$bestNN))
  } else {
    maxTerms=length(object$bestNN)
  }


  predictedValidation=t(apply(-kernelNewTrain,1,function(xx) {
    responsePredict=rep(NA,maxTerms)
    for(kk in 1:maxTerms)
    {
      nearest=sort(xx,index.return=T)$ix[1:object$bestNN[kk]]
      responsePredict[kk]=mean(object$responsesTrain[nearest,kk,drop=FALSE])
    }
    return(responsePredict)
  }))
  rm(kernelNewTrain)
  gc(verbose = FALSE)

  return(predictedValidation)

}


#' Predict Series regression
#'
#' This function is typically not directly used by the user; it is used inside  \code{\link{fitFlexCoDE}}
#'
#' @param object object of the class Series
#' @param xNew matrix with covariates where prediction will be calculated
#' @param maxTerms maximum number of expansion coefficients
#'
#' @return returns matrix where element (i,j) contains the estimate of the j-th expansion coefficient for the j-th sample
#' @export
#'
predict.SeriesKernel=function(object,kernelNewTrain,maxTerms=NULL)
{
  if(class(object)!="SeriesKernel")
    stop("Object has wrong class, should be SeriesKernel")



  if(!is.null(maxTerms))
  {
    maxTerms=min(maxTerms,length(object$bestNX))
  } else {
    maxTerms=length(object$bestNX)
  }

  m=dim(kernelNewTrain)[1] # New
  n=dim(kernelNewTrain)[2] # Old
  basisX=kernelNewTrain %*% object$basisX[,1:object$nXMax,drop=FALSE]
  basisX=1/n*basisX*matrix(rep(1/object$eigenValues[1:object$nXMax],m),m,object$nXMax,byrow=T)

  predictedValidation=matrix(NA,nrow(kernelNewTrain),maxTerms)
  for(kk in 1:maxTerms)
  {
    #predictedValidation[,kk]=basisX[,1:object$bestNX[kk],drop=F]%*%object$coefficientsMatrix[kk,1:object$bestNX[kk],drop=F]
    predictedValidation[,kk]=object$coefficientsMatrix[kk,1:object$bestNX[kk],drop=F]%*%t(basisX[,1:object$bestNX[kk],drop=F])
  }

  gc(verbose = FALSE)
  return(predictedValidation)

}



predict.SDMKernel=function(object,kernelNewTrain,maxTerms=NULL)
{
  if(class(object)!="SDMKernel")
    stop("Object has wrong class, should be SDMKernel")


  if(!is.null(maxTerms))
  {
    maxTerms=min(maxTerms,length(object$fittedReg))
  } else {
    maxTerms=length(object$fittedReg)
  }

  predictedValidation=matrix(NA,nrow(kernelNewTrain),maxTerms)
  predictedValidation[,1]=1
  for(i in 2:maxTerms)
  {
    predictedValidation[,i]=predict(object$fittedReg[[i]],kernlab::as.kernelMatrix(kernelNewTrain[,kernlab::SVindex(object$fittedReg[[i]])]))
  }

  rm(kernelNewTrain)
  gc(verbose = FALSE)

  return(predictedValidation)
}

