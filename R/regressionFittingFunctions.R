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
