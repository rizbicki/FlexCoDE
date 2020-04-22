#' Nearest Neighbors Regression
#'
#' This function is typically not directly used by the user; it is used inside  \code{\link{fitFlexCoDE}}
#'
#' @param x matrix with covariates that will be used for training
#' @param responses matrix where each column is a response for the training data
#' @param extra list with one component named nNeighVec, which contains a vetor with different number of neighbors; the function will choose the best value among them
#'
#' @return object of the class NN containing information need to perform prediction on new points
#' @export
regressionFunction.NN <- function(x, responses, extra = NULL) {
  n_obs <- nrow(x)
  n_basis <- ncol(responses)
  n_train <- round(0.7 * n_obs)

  random <- sample(n_obs)
  train_ids <- random[1:n_train]
  x_train <- x[train_ids, , drop = FALSE]
  z_train <- responses[train_ids, , drop = FALSE]

  n_validation <- n_obs - n_train
  x_validation <- x[-train_ids, , drop = FALSE]
  z_validation <- responses[-train_ids, , drop = FALSE]

  nns <- FNN::knnx.index(x_train, x_validation, k = n_train)

  nNeighVec <- extra$nn
  if (is.null(nNeighVec)) {
    nNeighVec <- round(seq(1, n_train, length.out = 100))
  }
  n_k <- length(nNeighVec)

  errors <- matrix(0.0, n_basis, n_k)
  for (ii in seq_len(n_validation)) {
    for (kk in 1:n_k) {
      z_predict <- colMeans(z_train[nns[ii, seq_len(nNeighVec[kk])], , drop = FALSE])
      errors[, kk] <- errors[, kk] + (z_predict - z_validation[ii, ]) ^ 2
    }
  }
  bestNN <- nNeighVec[apply(errors, 1, which.min)]

  return(structure(list(bestNN = bestNN,
                        xTrain = x,
                        responsesTrain = responses),
                   class = "NN"))
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
print.NN=function(regressionObject,bestI,nameCovariates)
{
  cat(paste("Number of neighbors chosen for each fitted regression:",paste(regressionObject$bestNN[1:bestI],collapse=", "),"\n"))
}



#' Nadaraya Watson Regression
#'
#' This function is typically not directly used by the user; it is used inside  \code{\link{fitFlexCoDE}}
#'
#' @param x matrix with covariates that will be used for training
#' @param responses matrix where each column is a response for the training data
#' @param extra list with one component named epsGrid, which contains a vetor with different bandwidths; the function will choose the best value among them
#'
#' @return object of the class NN containing information need to perform prediction on new points
#' @export
regressionFunction.NW=function(x,responses,extra=NULL)
{
  # Both x and responses are matrices
  n=dim(x)[1]
  random=sample(1:n)
  nTrain=round(0.7*n)
  xTrain=x[random[1:nTrain],,drop=FALSE]
  responsesTrain=responses[random[1:nTrain],,drop=FALSE]
  xValidation=x[random[-c(1:nTrain)],,drop=FALSE]
  responsesValidation=responses[random[-c(1:nTrain)],,drop=FALSE]

  distanceValidationTrain=fields::rdist(xValidation,xTrain)

  epsGrid=extra$epsGrid
  if(is.null(epsGrid))
    epsGrid=seq(median(distanceValidationTrain)/100,median(distanceValidationTrain),length.out = 20)



  error=matrix(NA,length(epsGrid),dim(responsesTrain)[2])
  for(ii in 1:length(epsGrid))
  {
    weights=exp(-(distanceValidationTrain/(epsGrid[ii]))^2)


    isThereNa=apply(weights,1,function(xx){sum(is.na(xx))>0})
    weights[isThereNa,]=1/ncol(weights)
    weights=weights/ncol(weights)

    predictedValidation=weights%*%responsesTrain

    error[ii,]=colMeans((predictedValidation-responsesValidation)^2)
  }

  bestEps=apply(error,2,function(xx){
    epsGrid[which.min(xx)]
  })
  rm(distanceValidationTrain)
  gc(verbose = FALSE)

  regressionObject=NULL
  regressionObject$bestEps=bestEps
  regressionObject$xTrain=x
  regressionObject$responsesTrain=responses
  class(regressionObject)="NW"
  rm(xTrain,responsesTrain)
  gc(verbose=FALSE)
  return(regressionObject)
}

#' Print function for object of the class NW
#'
#' This function is typically not directly used by the user; it is used inside  \code{\link{print.FlexCoDE}}, the print
#' method for the class FlexCoDE
#'
#' @param regressionObject of the class NW
#' @param bestI optimal number of expansion coefficients
#' @param nameCovariates name of the covariates
#'
#' @return prints characteristics of the regressions that were fitted
#'
print.NW=function(regressionObject,bestI,nameCovariates)
{
  cat(paste("Bandiwdth chosen for each fitted regression:",paste(regressionObject$bestEps[1:bestI],collapse=", "),"\n"))
}


#' SpAM Regression (Sparse Additive Model)
#'
#' This function is typically not directly used by the user; it is used inside  \code{\link{fitFlexCoDE}}
#'
#' @param x matrix with covariates that will be used for training
#' @param responses matrix where each column is a response for the training data
#' @param extra list with one component named sVec, which contains a vetor with different number of splins; the function will choose the best value among them. The list can also contain a component named
#' nCores which contains the number of cores to be used for parallel computing. Default is one.
#'
#'
#' @return object of the class SpAM containing information needed to perform prediction on new points
#'
#' @import SAM
#' @import doParallel
#' @import foreach
#'
#' @export
regressionFunction.SpAM=function(x,responses,extra=NULL)
{
  # Both x and responses are matrices
  n=dim(x)[1]
  random=sample(1:n)
  nTrain=round(0.7*n)
  xTrain=x[random[1:nTrain],,drop=FALSE]
  responsesTrain=responses[random[1:nTrain],,drop=FALSE]
  xValidation=x[random[-c(1:nTrain)],,drop=FALSE]
  responsesValidation=responses[random[-c(1:nTrain)],,drop=FALSE]


  sVec=extra$sVec
  if(is.null(sVec))
    sVec=round(seq(1,14,length.out = 6))


  nCores=extra$nCores
  if(is.null(nCores))
    nCores=1

  cl <- parallel::makeCluster(nCores)
  doParallel::registerDoParallel(cl)


  fittedReg <- foreach(ii=1:ncol(responsesTrain)) %dopar% {

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

  }

  parallel::stopCluster(cl)



  regressionObject=NULL
  regressionObject$fittedReg=fittedReg
  class(regressionObject)="SpAM"
  rm(fittedReg,responsesTrain,xTrain,xValidation,responsesValidation)
  gc(verbose=FALSE)
  return(regressionObject)
}


#' Print function for object of the class SpAM
#'
#' This function is typically not directly used by the user; it is used inside  \code{\link{print.FlexCoDE}}, the print
#' method for the class FlexCoDE
#'
#' @param regressionObject of the class SpAM
#' @param bestI optimal number of expansion coefficients
#' @param nameCovariates name of the covariates
#'
#' @return prints characteristics of the regressions that were fitted
#'
print.SpAM=function(regressionObject,bestI,nameCovariates)
{

  bestS=sapply(regressionObject$fittedReg, function(x)x$bestS)[1:bestI]
  cat(paste("Number of splines chosen for each fitted regression:",paste(bestS,collapse = ", "),"\n"))

  cat("\n")

  if(length(regressionObject$fittedReg[[1]]$whichCovariates)==1)
  {
    #bestS=t(sapply(regressionObject$fittedReg, function(x)x$whichCovariates)[1:bestI])
    return();
  } else {
    bestS=t(sapply(regressionObject$fittedReg, function(x)x$whichCovariates)[,1:bestI])
  }
  freq=colMeans(bestS)

  table=data.frame(covariate=order(freq,decreasing = TRUE),frequency=sort(freq,decreasing =  TRUE))
  cat(paste("How many times each covariate was selected: \n"))
  print(table)

  if(is.null(nameCovariates))
    nameCovariates=1:nrow(table)


  table$covariate=factor(table$covariate,levels=table$covariate)
  ggplot2::ggplot(table, ggplot2::aes(x=as.factor(covariate),y=frequency))+
    ggplot2::geom_bar(position="dodge",stat="identity") +
    ggplot2::coord_flip()+ggplot2::ylab("Frequency")+ggplot2::xlab("Covariate")+
    ggplot2::scale_x_discrete(labels=nameCovariates[as.numeric(as.character(table$covariate))])
}


#' Spectral Series Regression
#'
#' This function is typically not directly used by the user; it is used inside  \code{\link{fitFlexCoDE}}
#'
#' @param x matrix with covariates that will be used for training
#' @param responses matrix where each column is a response for the training data
#' @param extra list with two components: the first is named epsGrid, which contains a vetor with different number of bandwidths to be used in the gaussian kernel; the function will choose the best value among them; the second is called nXMax and contains a single integer number that describes what is the maximum number of spectral basis functions with be used
#'
#' @return object of the class Series containing information needed to perform prediction on new points
#' @export
regressionFunction.Series=function(x,responses,extra=NULL)
{
  # Both x and responses are matrices
  n=dim(x)[1]
  random=sample(1:n)
  nTrain=round(0.7*n)
  xTrain=x[random[1:nTrain],,drop=FALSE]
  responsesTrain=responses[random[1:nTrain],,drop=FALSE]
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

    basisX=Re(sqrt(nAll)*U[,1:nXMax,drop=FALSE])
    eigenValues=Re((lambda/nAll)[1:nXMax])

    coefficientsMatrix=1/nAll*t(responsesTrain)%*%basisX

    m=nrow(xValidation) # New

    basisXNew=kernelMatrixValidationTrain %*% basisX
    basisXNew=1/n*basisXNew*matrix(rep(1/eigenValues,m),m,ncol(basisX),byrow=T)

    errorsForEachReg=apply(as.matrix(1:nrow(coefficientsMatrix)),1,function(ii)
    {
      coeff=coefficientsMatrix[ii,,drop=FALSE]
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
print.Series=function(regressionObject,bestI,nameCovariates)
{
  cat(paste("Number of expansion coefficients chosen for each fitted regression:",paste(regressionObject$bestNX[1:bestI],collapse = ", "),"\n"))
  cat("\n")
  cat(paste("Best epsilon for spectral decomposition:",regressionObject$bestEps,"\n"))
}


#' Lasso Regression
#'
#' This function is typically not directly used by the user; it is used inside  \code{\link{fitFlexCoDE}}
#'
#' @param x matrix with covariates that will be used for training
#' @param responses matrix where each column is a response for the training data
#' @param extra this argument is ignored in this function
#'
#' @return object of the class Lasso containing information needed to perform prediction on new points
#' @export
regressionFunction.Lasso=function(x,responses,extra=NULL)
{
  # Both x and responses are matrices

  coeffs=apply(responses[,-1],2,function(yy){
    lasso.mod = glmnet::glmnet(x,yy, alpha =1)
    cv.out = glmnet::cv.glmnet(x, yy, alpha =1)
    bestlam = cv.out$lambda.min
    coeff = coefficients(lasso.mod , s = bestlam)
    return(as.numeric(coeff))
  })
  coeffs <- cbind(c(1,ncol(x)),coeffs)

  regressionObject=NULL
  regressionObject$coefficients=coeffs
  class(regressionObject)="Lasso"
  return(regressionObject)
}

#' Print function for object of the class Lasso
#'
#' This function is typically not directly used by the user; it is used inside  \code{\link{print.FlexCoDE}}, the print
#' method for the class FlexCoDE
#'
#' @param regressionObject of the class Lasso
#' @param bestI optimal number of expansion coefficients
#' @param nameCovariates name of the covariates
#'
#' @return prints characteristics of the regressions that were fitted
#'
print.Lasso=function(regressionObject,bestI,nameCovariates)
{


  whichCoefficients=t(apply(regressionObject$coefficients[-1,1:bestI],1,function(x)(abs(x)>1e-20)))

  freq=colMeans(whichCoefficients)

  table=data.frame(covariate=order(freq,decreasing = TRUE),frequency=sort(freq,decreasing =  TRUE))
  cat(paste("How many times each covariate was selected: \n"))
  print(table)

  if(is.null(nameCovariates))
    nameCovariates=1:nrow(table)


  table$covariate=factor(table$covariate,levels=table$covariate)
  ggplot2::ggplot(table, ggplot2::aes(x=as.factor(covariate),y=frequency))+
    ggplot2::geom_bar(position="dodge",stat="identity") +
    ggplot2::coord_flip()+ggplot2::ylab("Frequency")+ggplot2::xlab("Covariate")+
    ggplot2::scale_x_discrete(labels=nameCovariates[as.numeric(as.character(table$covariate))])


}


#' Forest Regression
#'
#' This function is typically not directly used by the user; it is used inside  \code{\link{fitFlexCoDE}}
#'
#' @param x matrix with covariates that will be used for training
#' @param responses matrix where each column is a response for the training data
#' @param extra list with one components named p0Vec, which contains a vetor with different number of variables randomly sampled as candidates at each split of the forest regression (aka mtry in randomForest package); the function will choose the best value among them;  one component named ntree which contains the number of tree to be used by the forest (default is 500), and  one component named maxnodes which contains the number of f terminal nodes trees in the forest can have (if not given, trees are grown to the maximum possible). The list can also contain a component named
#' nCores which contains the number of cores to be used for parallel computing. Default is one.
#'
#' @import randomForest
#'
#' @return object of the class Forest containing information needed to perform prediction on new points
#' @export
regressionFunction.Forest=function(x,responses,extra=NULL)
{

  # Both x and responses are matrices
  n=dim(x)[1]
  random=sample(1:n)
  nTrain=round(0.7*n)
  xTrain=x[random[1:nTrain],,drop=FALSE]
  responsesTrain=responses[random[1:nTrain],,drop=FALSE]
  xValidation=x[random[-c(1:nTrain)],,drop=FALSE]
  responsesValidation=responses[random[-c(1:nTrain)],,drop=FALSE]

  p0Vec=extra$p0Vec
  if(is.null(p0Vec))
    p0Vec=round(seq(1,ncol(xTrain),length.out = 5))

  ntree=extra$ntree
  if(is.null(ntree))
    ntree=500

  maxnodes=extra$maxnodes

  nCores=extra$nCores
  if(is.null(nCores))
    nCores=1

  cl <- parallel::makeCluster(nCores)
  doParallel::registerDoParallel(cl)


  fittedReg <- foreach(ii=1:ncol(responsesTrain)) %dopar% {
    error=rep(NA,length(p0Vec))
    for(s in 1:length(p0Vec))
    {
      ajuste = randomForest::randomForest(x=xTrain,
                                          y=responsesTrain[,ii,drop=FALSE],
                                          mtry=p0Vec[s],
                                          importance = FALSE)
      colnames(xValidation)=NULL
      predito = predict(ajuste, newdata = xValidation)
      error[s]=mean((predito-responsesValidation[,ii,drop=FALSE])^2)
    }
    bestP0=p0Vec[which.min(error)]
    #ajuste = randomForest(x=xTrain,y=responsesTrain[,ii,drop=FALSE],mtry=bestP0,importance = TRUE)
    ajuste = randomForest::randomForest(x=x,y=responses[,ii,drop=FALSE],mtry=bestP0,importance = TRUE,ntree=ntree,maxnodes=maxnodes)
    object=NULL
    object$fit=ajuste
    object$importance=ajuste$importance[,1]
    object$bestP0=bestP0
    object$errors=ajuste$mse
    gc(verbose=FALSE)
    return(object)
  }

  parallel::stopCluster(cl)

  regressionObject=NULL
  regressionObject$fittedReg=fittedReg
  class(regressionObject)="Forest"
  return(regressionObject)
}

#' Print function for object of the class Forest
#'
#' This function is typically not directly used by the user; it is used inside  \code{\link{print.FlexCoDE}}, the print
#' method for the class FlexCoDE
#'
#' @param regressionObject of the class Forest
#' @param bestI optimal number of expansion coefficients
#' @param nameCovariates name of the covariates
#'
#' @return prints characteristics of the regressions that were fitted
#'
print.Forest=function(regressionObject,bestI,nameCovariates)
{

  if(length(regressionObject$fittedReg[[1]]$importance)==1)
  {
    importance=t(sapply(regressionObject$fittedReg,function(x)x$importance)[1:bestI])

  } else {
    importance=sapply(regressionObject$fittedReg,function(x)x$importance)[,1:bestI]
  }
  freq=rowMeans(importance)

  table=data.frame(covariate=order(freq,decreasing = TRUE),frequency=sort(freq,decreasing =  TRUE))
  cat(paste("Average Importance of each covariate: \n"))
  print(table)


  if(is.null(nameCovariates))
    nameCovariates=1:nrow(table)

  table$covariate=factor(table$covariate,levels=table$covariate)
  ggplot2::ggplot(table, ggplot2::aes(x=as.factor(covariate),y=frequency))+
    ggplot2::geom_bar(position="dodge",stat="identity") +
    ggplot2::coord_flip()+ggplot2::ylab("Average Importance")+ggplot2::xlab("Covariate")+
    ggplot2::scale_x_discrete(labels=nameCovariates[as.numeric(as.character(table$covariate))])



}






#' XGBoost
#'
#' This function is typically not directly used by the user; it is used inside  \code{\link{fitFlexCoDE}}
#'
#' @param x matrix with covariates that will be used for training
#' @param responses matrix where each column is a response for the training data
#' @param extra list with one components named p0Vec, which contains a vetor with different number of variables randomly sampled as candidates at each split of the forest regression (aka mtry in randomForest package); the function will choose the best value among them;  one component named ntree which contains the number of tree to be used by the forest (default is 500), and  one component named maxnodes which contains the number of f terminal nodes trees in the forest can have (if not given, trees are grown to the maximum possible). The list can also contain a component named
#' nCores which contains the number of cores to be used for parallel computing. Default is one.
#'
#' @import xgboost
#'
#' @return object of the class XGBoost containing information needed to perform prediction on new points
#' @export
#'
regressionFunction.XGBoost=function(x,responses,extra=NULL)
{

  # Both x and responses are matrices

  ninter=extra$ninter
  if(is.null(ninter))
    ninter=500

  maxnodes=extra$maxnodes

  nCores=extra$nCores
  if(is.null(nCores))
    nCores=1

  cl <- parallel::makeCluster(nCores)
  doParallel::registerDoParallel(cl)


  fittedReg <- foreach(ii=2:ncol(responses)) %dopar% {
    bst <- xgboost::xgboost(names(x),data =x,label=responses[,ii,drop=FALSE],
                              nthread = 1,
                            nround = ninter,
                            objective="reg:squarederror",
                            verbose=FALSE)

    object=NULL
    object$fit=bst
    object$importance=xgboost::xgb.importance(colnames(x), model = bst)
    gc(verbose=FALSE)
    return(object)
  }

  parallel::stopCluster(cl)

  regressionObject=NULL
  regressionObject$fittedReg=fittedReg
  class(regressionObject)="XGBoost"
  return(regressionObject)
}

#' Print function for object of the class XGBoost
#'
#' This function is typically not directly used by the user; it is used inside  \code{\link{print.FlexCoDE}}, the print
#' method for the class FlexCoDE
#'
#' @param regressionObject of the class XGBoost
#' @param bestI optimal number of expansion coefficients
#' @param nameCovariates name of the covariates
#'
#' @return prints characteristics of the regressions that were fitted
#'
print.XGBoost=function(regressionObject,bestI,nameCovariates)
{


  if(length(regressionObject$fittedReg[[1]]$importance$Gain)==1)
  {
    importance=t(sapply(regressionObject$fittedReg,function(x)x$importance$Gain)[1:(bestI-1)])

  } else {
    importance=sapply(regressionObject$fittedReg,function(x)x$importance$Gain)[,1:(bestI-1),drop=FALSE]
  }
  freq=rowMeans(importance)
  if(!is.null(regressionObject$names.covariates))
  {
    freq=freq[match(regressionObject$fittedReg[[1]]$importance$Feature,
                    nameCovariates)]
  } else {
    freq=freq[order(regressionObject$fittedReg[[1]]$importance$Feature)]
  }


  table=data.frame(covariate=order(freq,decreasing = TRUE),frequency=sort(freq,decreasing =  TRUE))
  cat(paste("Average Importance of each covariate: \n"))
  print(table)


  if(is.null(nameCovariates))
    nameCovariates=1:nrow(table)

  table$covariate=factor(table$covariate,levels=table$covariate)
  ggplot2::ggplot(table, ggplot2::aes(x=as.factor(covariate),y=frequency))+
    ggplot2::geom_bar(position="dodge",stat="identity") +
    ggplot2::coord_flip()+ggplot2::ylab("Average Importance")+ggplot2::xlab("Covariate")+
    ggplot2::scale_x_discrete(labels=nameCovariates[as.numeric(as.character(table$covariate))])



}



