#' Predict NN regression
#'
#' This function is typically not directly used by the user; it is
#' used inside \code{\link{fitFlexCoDE}}
#'
#' @param object object of the class NN
#' @param xNew matrix with covariates where prediction will be calculated
#' @param maxTerms maximum number of expansion coefficients
#'
#' @return Returns a matrix where element (i, j) contains the estimate
#'   of the j-th expansion coefficient for the i-th sample
#' @export
predict.NN <- function(object, xNew, maxTerms = NULL) {
  n_predict <- nrow(xNew)
  n_train <- nrow(object$xTrain)
  maxTerms <- min(maxTerms, length(object$bestNN))

  nns <- FNN::knnx.index(object$xTrain, xNew, k = n_train)

  preds <- matrix(NA, n_predict, maxTerms)
  for (ii in seq_len(n_predict)) {
    for (jj in seq_len(maxTerms)) {
      preds[ii, jj] <- mean(object$responsesTrain[nns[ii, seq_len(object$bestNN[jj])], jj])
    }
  }

  return(preds)
}

#' Predict NW regression
#'
#' This function is typically not directly used by the user; it is used inside  \code{\link{fitFlexCoDE}}
#'
#' @param object object of the class NW
#' @param xNew matrix with covariates where prediction will be calculated
#' @param maxTerms maximum number of expansion coefficients
#'
#' @return returns matrix where element (i,j) contains the estimate of the j-th expansion coefficient for the j-th sample
#' @export
#'
predict.NW=function(object,xNew,maxTerms=NULL)
{
  if(class(object)!="NW")
    stop("Object has wrong class, should be NW")

  distanceNewTrain=fields::rdist(xNew,object$xTrain)

  if(!is.null(maxTerms))
  {
    maxTerms=min(maxTerms,length(object$bestEps))
  } else {
    maxTerms=length(object$bestEps)
  }

  predictedValidation=matrix(NA,nrow(xNew),maxTerms)
  for(kk in 1:maxTerms)
  {

    weights=exp(-(distanceNewTrain/(object$bestEps[kk]))^2)
    weights=weights/rowSums(weights)
    isThereNa=apply(weights,1,function(xx){sum(is.na(xx))>0})
    weights[isThereNa,]=1/ncol(weights)

    predictedValidation[,kk]=weights%*%object$responsesTrain[,kk,drop=F]
  }


  rm(distanceNewTrain)
  gc(verbose = FALSE)

  return(predictedValidation)
}


#' Predict SpAM regression
#'
#' This function is typically not directly used by the user; it is used inside  \code{\link{fitFlexCoDE}}
#'
#' @param object object of the class SpAM
#' @param xNew matrix with covariates where prediction will be calculated
#' @param maxTerms maximum number of expansion coefficients
#'
#' @return returns matrix where element (i,j) contains the estimate of the j-th expansion coefficient for the j-th sample
#' @export
#'
predict.SpAM=function(object,xNew,maxTerms=NULL)
{
  if(class(object)!="SpAM")
    stop("Object has wrong class, should be SAM")

  if(!is.null(maxTerms))
  {
    maxTerms=min(maxTerms,length(object$fittedReg))
  } else {
    maxTerms=length(object$fittedReg)
  }

  predictedValidation=apply(as.matrix(1:maxTerms),1,function(xx)
  {
    out.tst = try(predict(object$fittedReg[[xx]]$fit,xNew),silent = TRUE)
    if(class(out.tst)=="try-error")
    {
      return(NA)
    }
    predicted=out.tst[[1]][,object$fittedReg[[xx]]$bestCol]
    return(predicted)
  })

  whichNa=sapply(predictedValidation,function(x)any(is.na(x)))
  if(all(!whichNa))
  {
    return(predictedValidation)
  }
  whichSelect=rep(F,length(predictedValidation))
  whichSelect[1:(which.max(whichNa)-1)]=T
  predictedValidation=subset(predictedValidation,subset = whichSelect)
  predictedValidation=sapply(predictedValidation,function(x)x)
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
predict.Series=function(object,xNew,maxTerms=NULL)
{
  if(class(object)!="Series")
    stop("Object has wrong class, should be Series")

  distanceNewTrain=fields::rdist(xNew,object$xTrain)
  kernelNewOld=exp(-distanceNewTrain^2/(4*object$bestEps))

  if(!is.null(maxTerms))
  {
    maxTerms=min(maxTerms,length(object$bestNX))
  } else {
    maxTerms=length(object$bestNX)
  }


  m=dim(kernelNewOld)[1] # New
  n=dim(kernelNewOld)[2] # Old
  basisX=kernelNewOld %*% object$basisX[,1:object$nXMax,drop=FALSE]
  basisX=1/n*basisX*matrix(rep(1/object$eigenValues[1:object$nXMax],m),m,object$nXMax,byrow=T)

  predictedValidation=matrix(NA,nrow(distanceNewTrain),maxTerms)
  for(kk in 1:maxTerms)
  {
    #predictedValidation[,kk]=basisX[,1:object$bestNX[kk],drop=F]%*%object$coefficientsMatrix[kk,1:object$bestNX[kk],drop=F]
    predictedValidation[,kk]=object$coefficientsMatrix[kk,1:object$bestNX[kk],drop=F]%*%t(basisX[,1:object$bestNX[kk],drop=F])
  }

  gc(verbose = FALSE)
  return(predictedValidation)
}



#' Predict Lasso regression
#'
#' This function is typically not directly used by the user; it is used inside  \code{\link{fitFlexCoDE}}
#'
#' @param object object of the class Lasso
#' @param xNew matrix with covariates where prediction will be calculated
#' @param maxTerms maximum number of expansion coefficients
#'
#' @return returns matrix where element (i,j) contains the estimate of the j-th expansion coefficient for the j-th sample
#' @export
#'
predict.Lasso=function(object,xNew,maxTerms=NULL)
{
  if(class(object)!="Lasso")
    stop("Object has wrong class, should be Lasso")

  if(!is.null(maxTerms))
  {
    maxTerms=min(maxTerms,ncol(object$coefficients))
  } else {
    maxTerms=ncol(object$coefficients)
  }


  predictedValidation=xNew%*%object$coefficients[-1,1:maxTerms,drop=FALSE]
  predictedValidation=predictedValidation+matrix(object$coefficients[1,1:maxTerms,drop=FALSE],nrow(predictedValidation),ncol(predictedValidation),byrow = TRUE)

  return(predictedValidation)
}


#' Predict XGBoost
#'
#' This function is typically not directly used by the user; it is used inside  \code{\link{fitFlexCoDE}}
#'
#' @param object object of the class XGBoost
#' @param xNew matrix with covariates where prediction will be calculated
#' @param maxTerms maximum number of expansion coefficients
#'
#' @return returns matrix where element (i,j) contains the estimate of the j-th expansion coefficient for the j-th sample
#' @export
#'
predict.XGBoost=function(object,xNew,maxTerms=NULL)
{
  if(class(object)!="XGBoost")
    stop("Object has wrong class, should be XGBoost")

  if(!is.null(maxTerms))
  {
    if(maxTerms==1)
      return(matrix(1,nrow(xNew),1))
    maxTerms=min(maxTerms,length(object$fittedReg)+1)
  } else {
    maxTerms=length(object$fittedReg)+1
  }

  predictedValidation=apply(as.matrix(2:maxTerms),1,function(xx)
  {
    predicted=predict(object$fittedReg[[xx-1]]$fit,xNew)
    return(predicted)
  })

  return(cbind(1,predictedValidation))
}



#' Predict Forest regression
#'
#' This function is typically not directly used by the user; it is used inside  \code{\link{fitFlexCoDE}}
#'
#' @param object object of the class Forest
#' @param xNew matrix with covariates where prediction will be calculated
#' @param maxTerms maximum number of expansion coefficients
#'
#' @return returns matrix where element (i,j) contains the estimate of the j-th expansion coefficient for the j-th sample
#' @export
#'
predict.Forest=function(object,xNew,maxTerms=NULL)
{
  if(class(object)!="Forest")
    stop("Object has wrong class, should be Forest")

  if(!is.null(maxTerms))
  {
    maxTerms=min(maxTerms,length(object$fittedReg))
  } else {
    maxTerms=length(object$fittedReg)
  }

  predictedValidation=apply(as.matrix(1:maxTerms),1,function(xx)
  {
    colnames(xNew)=NULL
    predicted = predict(object$fittedReg[[xx]]$fit,newdata=xNew)
    return(predicted)
  })

  return(predictedValidation)
}


