
#' Title
#'
#' @param x
#' @param responses
#' @param extra
#'
#' @return x
#' @export
regressionFunction.NN=function(x,responses,extra=NULL)
{
  # Both x and responses are matrices
  n=dim(x)[1]
  random=sample(1:n)
  nTrain=round(0.7*n)
  xTrain=x[random[1:nTrain],]
  responsesTrain=responses[random[1:nTrain],]
  xValidation=x[random[-c(1:nTrain)],,drop=FALSE]
  responsesValidation=responses[random[-c(1:nTrain)],,drop=FALSE]

  distanceValidationTrain=fields::rdist(xValidation,xTrain)

  nNeighVec=extra$nn
  if(is.null(nNeighVec))
    nNeighVec=round(seq(1,dim(responsesTrain)[1],length.out = 100))


  error=matrix(NA,length(nNeighVec),dim(responsesTrain)[2])
  for(ii in 1:length(nNeighVec))
  {
    predictedValidation=t(apply(distanceValidationTrain,1,function(xx) {
      nearest=sort(xx,index.return=T)$ix[1:nNeighVec[ii]]
      return(colMeans(responsesTrain[nearest,,drop=FALSE]))
    }))
    error[ii,]=colMeans((predictedValidation-responsesValidation)^2)
  }

  bestNN=apply(error,2,function(xx){
    nNeighVec[which.min(xx)]
  })
  rm(distanceValidationTrain)
  gc(verbose = FALSE)

  regressionObject=NULL
  regressionObject$bestNN=bestNN
  regressionObject$xTrain=x
  regressionObject$responsesTrain=responses
  class(regressionObject)="NN"
  rm(xTrain,responsesTrain)
  gc(verbose=FALSE)
  return(regressionObject)
}

#' Title
#'
#' @param x
#' @param responses
#' @param extra
#'
#' @return x
#' @export
#'
regressionFunction.SpAM=function(x,responses,extra=NULL)
{
  # Both x and responses are matrices
  n=dim(x)[1]
  random=sample(1:n)
  nTrain=round(0.7*n)
  xTrain=x[random[1:nTrain],]
  responsesTrain=responses[random[1:nTrain],]
  xValidation=x[random[-c(1:nTrain)],,drop=FALSE]
  responsesValidation=responses[random[-c(1:nTrain)],,drop=FALSE]


  sVec=extra$sVec
  if(is.null(sVec))
    sVec=round(seq(1,14,length.out = 6))


  fittedReg=apply(as.matrix(1:ncol(responsesTrain)),1,function(ii){

    bestCol=bestError=rep(NA,length(sVec))
    for(s in 1:length(bestCol))
    {
      fit=try(SAM::samQL(xTrain,responsesTrain[,ii,drop=FALSE],p = sVec[s]),silent = TRUE)
      if(class(fit)=="try-error")
      {
        bestError[s]=Inf
        next;
      }
      out = predict(fit,xValidation)
      out$values[is.na(out$values)]=0
      errorPerLambda=colMeans((out[[1]]-matrix(responsesValidation[,ii,drop=FALSE],nrow(out[[1]]),ncol(out[[1]])))^2)
      bestCol[s]=which.min(errorPerLambda)
      bestError[s]=min(errorPerLambda)
    }

    bestS=sVec[which.min(bestError)]
    fit=SAM::samQL(x,responses[,ii,drop=FALSE],p = bestS)

    object=NULL
    object$fit=fit
    object$bestS=bestS
    object$bestCol=bestCol[which.min(bestError)]
    object$whichCovariates=fit$func_norm[,object$bestCol]>1e-23

    return(object)
  })

  regressionObject=NULL
  regressionObject$fittedReg=fittedReg
  class(regressionObject)="SpAM"
  rm(fittedReg,responsesTrain,xTrain,xValidation,responsesValidation)
  gc(verbose=FALSE)
  return(regressionObject)
}

#' Title
#'
#' @param x
#' @param responses
#' @param extra
#'
#' @return X
#' @export
#'
regressionFunction.Series=function(x,responses,extra=NULL)
{
  # Both x and responses are matrices
  n=dim(x)[1]
  random=sample(1:n)
  nTrain=round(0.7*n)
  xTrain=x[random[1:nTrain],]
  responsesTrain=responses[random[1:nTrain],]
  xValidation=x[random[-c(1:nTrain)],,drop=FALSE]
  responsesValidation=responses[random[-c(1:nTrain)],,drop=FALSE]

  distanceValidationTrain=fields::rdist(xValidation,xTrain)
  distanceTrainTrain=fields::rdist(xTrain,xTrain)

  epsGrid=extra$eps
  if(is.null(epsGrid))
    epsGrid=seq(median(distanceValidationTrain^2)/100,median(distanceValidationTrain^2),length.out = 20)


  error=matrix(NA,length(epsGrid))
  nXBestEps=matrix(NA,length(epsGrid),ncol(responsesValidation))
  for(ii in 1:length(epsGrid))
  {
    kernelMatrix=exp(-distanceTrainTrain^2/(4*epsGrid[ii]))
    kernelMatrixValidationTrain=exp(-distanceValidationTrain^2/(4*epsGrid[ii]))
    nAll = dim(kernelMatrix)[1]

    if(!is.null(extra$nXMax))
    {
      nXMax=min(nrow(xTrain)-15,extra$nXMax)
    } else {
      nXMax=nrow(xTrain)-15
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

    m=nrow(xValidation) # New

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
    error[ii]=mean(errorsForEachReg[1,])
    nXBestEps[ii,]=errorsForEachReg[2,]
    rm(basisXNew)
    gc(verbose=FALSE)
  }

  bestEps=epsGrid[which.min(error)]
  bestNX=nXBestEps[which.min(error),]

  kernelMatrix=exp(-distanceTrainTrain^2/(4*bestEps))
  nAll = dim(kernelMatrix)[1]

  nXMax=max(bestNX)
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

  basisX=Re(sqrt(nAll)*U[,1:nXMax,drop=F])
  eigenValues=Re((lambda/nAll)[1:nXMax,drop=F])

  coefficientsMatrix=1/nAll*t(responsesTrain)%*%basisX

  rm(distanceValidationTrain,distanceTrainTrain,kernelMatrixValidationTrain,kernelMatrix)
  gc(verbose = FALSE)

  regressionObject=NULL
  regressionObject$nXMax=nXMax
  regressionObject$bestNX=bestNX # vector with best cutoffs
  regressionObject$bestEps=bestEps
  regressionObject$xTrain=xTrain
  regressionObject$coefficientsMatrix=coefficientsMatrix
  regressionObject$basisX=basisX
  regressionObject$eigenValues=eigenValues
  class(regressionObject)="Series"
  rm(xTrain,responsesTrain,basisX,eigenValues,coefficientsMatrix)
  gc(verbose=FALSE)
  return(regressionObject)
}


#' Title
#'
#' @param x
#' @param responses
#' @param extra
#'
#' @return x
#' @export
#'
regressionFunction.Lasso=function(x,responses,extra=NULL)
{
  # Both x and responses are matrices

  coeffs=apply(responses,2,function(yy){
    lasso.mod = glmnet::glmnet(x,yy, alpha =1)
    cv.out = glmnet::cv.glmnet(x, yy, alpha =1)
    bestlam = cv.out$lambda.min
    coeff = coefficients(lasso.mod , s = bestlam)
    return(as.numeric(coeff))
  })

  regressionObject=NULL
  regressionObject$coefficients=coeffs
  class(regressionObject)="Lasso"
  return(regressionObject)
}


#' Title
#'
#' @param x
#' @param responses
#' @param extra
#'
#' @return x
#' @export
#'
regressionFunction.Forest=function(x,responses,extra=NULL)
{

  # Both x and responses are matrices
  n=dim(x)[1]
  random=sample(1:n)
  nTrain=round(0.7*n)
  xTrain=x[random[1:nTrain],]
  responsesTrain=responses[random[1:nTrain],]
  xValidation=x[random[-c(1:nTrain)],,drop=FALSE]
  responsesValidation=responses[random[-c(1:nTrain)],,drop=FALSE]

  p0Vec=extra$p0Vec
  if(is.null(p0Vec))
    p0Vec=round(seq(1,ncol(xTrain),length.out = 5))



  fittedReg=apply(as.matrix(1:ncol(responsesTrain)),1,function(ii){
    error=rep(NA,length(p0Vec))
    for(s in 1:length(p0Vec))
    {
      ajuste = randomForest::randomForest(x=xTrain,y=responsesTrain[,ii,drop=FALSE],mtry=p0Vec[s],importance = FALSE)
      predito = predict(ajuste, newdata = xValidation)
      error[s]=mean((predito-responsesValidation[,ii,drop=FALSE])^2)
    }
    bestP0=p0Vec[which.min(error)]
    #ajuste = randomForest(x=xTrain,y=responsesTrain[,ii,drop=FALSE],mtry=bestP0,importance = TRUE)
    ajuste = randomForest::randomForest(x=x,y=responses[,ii,drop=FALSE],mtry=bestP0,importance = TRUE)
    object=NULL
    object$fit=ajuste
    object$importance=ajuste$importance[,1]
    object$bestP0=bestP0
    object$errors=ajuste$mse
    gc(verbose=FALSE)
    return(object)
  })

  regressionObject=NULL
  regressionObject$fittedReg=fittedReg
  class(regressionObject)="Forest"
  return(regressionObject)
}

