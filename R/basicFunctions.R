#' FlexCoDE Fit Conditional Density Estimation via Regression
#'
#' @param xTrain Covariates x used to train the model (one observation per row)
#' @param zTrain Responses z used to train the model  (matrix with one column; one observation per row)
#' @param xValidation Covariates x used to tune the model (one observation per row; same number of columns as xTrain)
#' @param zValidation Responses z used to tune the model  (matrix with one column; one observation per row)
#' @param xTest Covariates x used to estimate risk of final model (one observation per row; same number of columns as xTrain). Default is NULL
#' @param zTest Responses z used to estimate risk of final model  (matrix with one column; one observation per row). Default is NULL
#' @param nIMax Maximum possible number of components of the series expansion (that is, the function will find the best I<nIMax). Default is 100
#' @param regressionFunction a function indicating which regression method will be used to estimate the expansion coefficients. Currently can be one of
#' @param regressionFunction.extra extra parameters to be sent to regression function; see the regression you want to use to check what are the available parameters
#' @param system Basis for z. Current options are "Fourier", "Cosine" and "discrete". Default is "Fourier"
#' @param chooseDelta Should delta, the cutoff to remove spurious bumps, be chosen?
#' @param deltaGrid Grid of threshold deltas (betwen 0 and 0.5). Default value is seq(0,0.4,0.05).
#' @param verbose Should we print what we are doing? Default is TRUE.
#'
#' @return Returns the fitted estimated conditional density, and object of the class FlexCoDE
#' @export
fitFlexCoDE=function(xTrain,zTrain,xValidation,zValidation,xTest=NULL,zTest=NULL,nIMax=length(zTrain),regressionFunction,regressionFunction.extra=NULL,system="Fourier",deltaGrid=seq(0,0.4,0.05),chooseDelta=TRUE,verbose=TRUE)
{
  objectCDE=NULL
  objectCDE$zMax=max(zTrain)
  objectCDE$zMin=min(zTrain)
  zTrain=(zTrain-objectCDE$zMin)/(objectCDE$zMax-objectCDE$zMin)

  class(objectCDE)="FlexCoDE"

  if(verbose) print("Transforming Response")
  responseFourier=calculateBasis(zTrain,nIMax,system)
  if(verbose) print("Fitting Regression Functions")
  regressionObject=regressionFunction(x=xTrain,responses=responseFourier,extra=regressionFunction.extra)
  objectCDE$nIMax=nIMax
  objectCDE$system=system
  objectCDE$zTrain=zTrain
  objectCDE$xTrain=xTrain
  objectCDE$regressionObject=regressionObject
  rm(regressionObject,xTrain,zTrain,responseFourier)
  gc(verbose = FALSE)

  zValidation=(zValidation-objectCDE$zMin)/(objectCDE$zMax-objectCDE$zMin)

  basisZValidation=calculateBasis(zValidation,objectCDE$nIMax,objectCDE$system) # returns matrix length(z)xnIMax with the basis for z

  if(verbose) print("Tuning Number of Expansion Coefficients (I)")
  coefficientsXValidation=predict(objectCDE$regressionObject,xValidation)
  term1=1/2*colMeans(coefficientsXValidation^2)
  term1=cumsum(term1)

  term2=colMeans(coefficientsXValidation*basisZValidation[,1:ncol(coefficientsXValidation),drop=F])
  term2=cumsum(term2)
  objectCDE$errors=term1-term2
  objectCDE$bestI=which.min(objectCDE$errors)
  objectCDE$bestError=min(objectCDE$errors)

  if(verbose) print("Choosing optimal cutoff Delta")
  delta=chooseDelta(objectCDE, xValidation,zValidation,deltaGrid)
  objectCDE$bestDelta=delta

  if(!is.null(xTest)&!is.null(zTest))
  {
    if(verbose) print("Estimating risk on test set")
    error=estimateErrorFlexCoDE(objectCDE,xTest,zTest,se=TRUE)
    objectCDE$estimatedRisk=error
  }

  if(objectCDE$bestI==objectCDE$nIMax)
    warning("bestI=nIMax, try increasin nIMax if you want to improve performance")

  return(objectCDE)
}

#' Title Choose threshold value
#'
#' @param objectCDE An object of the class FlexCoDE with a fitted CDE, typically fitted used fitFlexCoDE beforehand
#' @param xValidation Covariates x used to validate (tune) the model (one x observation per row).
#' @param zValidation Responses z used to validate (tune) the model  (matrix with 1 column). Each row corresponds to a row of the xValidation argument
#'
#' @return Best Delta
chooseDelta = function(objectCDE, xValidation,zValidation,deltaGrid=seq(0,0.4,0.05))
{
  if(class(objectCDE)!='FlexCoDE')
    stop("objectCDE should be of class FlexCoDE")
  error=rep(NA,length(deltaGrid))
  cat("\n Progress Bar:\n")
  for(ii in 1:length(deltaGrid))
  {
    cat(paste(c(rep("|",ii),rep(" ",length(deltaGrid)-ii),"|\n"),collapse=""))
    objectCDE$bestDelta=deltaGrid[ii]
    estimateErrors=estimateErrorFlexCoDE(objectCDE=objectCDE,xTest=xValidation,zTest=zValidation,se=FALSE)
    error[ii]=estimateErrors
  }
  #plot(error)
  whichMin=(1:length(error))[error==min(error)]
  bestDelta=deltaGrid[max(whichMin)]
  return(bestDelta)
}


#' Estimate error (risk) of FlexCoDE object via test set
#'
#' @param objectCDE is an object of the class FlexCoDE
#' @param xTest Covariates x of the sample used to test the model (one observation per row)
#' @param zTest Response z of the sample used to test the model (one observation per row)
#' @param se Should standard error be computed? Default is TRUE
#'
#' @return Estimated error (with SE if desired)
#' @export
#'
estimateErrorFlexCoDE=function(objectCDE=objectCDE,xTest,zTest,se=TRUE)
{
  zGrid=seq(objectCDE$zMin[1],objectCDE$zMax[1],length.out=500)

  predictedComplete=predictFlexCoDE(objectCDE,xNew = xTest,B=length(zGrid))
  predictedComplete=predictedComplete$CDE*(objectCDE$zMax-objectCDE$zMin)

  colmeansComplete=colMeans(predictedComplete^2)
  sSquare=mean(colmeansComplete)

  n=length(zTest)
  predictedObserved=apply(as.matrix(1:n),1,function(xx) { index=which.min(abs(zTest[xx]-zGrid))
  return(predictedComplete[xx,index])
  })
  likeli=mean(predictedObserved)

  if(!se)
    return(1/2*sSquare-likeli)

  # Bootstrap
  output=NULL
  output$mean=1/2*sSquare-likeli

  boot=1000
  meanBoot=apply(as.matrix(1:boot),1,function(xx){
    sampleBoot=sample(1:n,replace=T)

    predictedCompleteBoot=predictedComplete[sampleBoot,]
    zTestBoot=zTest[sampleBoot]

    colmeansComplete=colMeans(predictedCompleteBoot^2)
    sSquare=mean(colmeansComplete)

    predictedObserved=apply(as.matrix(1:n),1,function(xx) { index=which.min(abs(zTestBoot[xx]-zGrid))
    return(predictedCompleteBoot[xx,index])
    })
    likeli=mean(predictedObserved)
    return(1/2*sSquare-likeli)
  })
  output$seBoot=sqrt(var(meanBoot))
  return(output)


}


#' Evaluates the estimated  density of new observations (testing points)
#'
#' @param objectCDE Object of the class "FlexCoDE", typically fitted used \code{\link{fitFlexCoDE}} beforehand
#' @param xNew Matrix with nTest rows and same number of columns as xTrain, containing x's for which the estimates are desired.
#' @param B Number of point where f(z|x) will be evaluated (on the z scale). This will be equally spaced between zMin and zMax
#'
#' @return The return value is an object with the following components
#' \item{z}{Points where the density was evaluate}
#' \item{CDE }{Matrix with value of the density at points z. Each row corresponds to a different observation x (i-th row of CDE corresponds to i-th row of xTest).}
#' @export
#'
predict.FlexCoDE=function(objectCDE,xNew,B=1000)
{
  if(class(objectCDE)!="FlexCoDE")
    stop("Object should be of type FlexCoDE")
  zGrid=seq(from=0,to=1,length.out=B)

  if(is.null(objectCDE$bestI))
    objectCDE$bestI=objectCDE$nIMax

  coeff=predict(objectCDE$regressionObject,xNew,maxTerms=objectCDE$bestI)

  basisZNew=calculateBasis(zGrid,objectCDE$bestI,objectCDE$system) # returns matrix length(z)xnIMax with the basis for z

  estimates=coeff%*%t(basisZNew)

  binSize=(1)/(B+1)

  delta=ifelse(!is.null(objectCDE$bestDelta),objectCDE$bestDelta,0)


  estimates=t(apply(estimates,1,function(xx).normalizeDensity(binSize,xx,delta)))

  estimates=estimates/(objectCDE$zMax-objectCDE$zMin)
  returnValue=NULL
  returnValue$CDE=estimates
  returnValue$z=seq(from=objectCDE$zMin,to=objectCDE$zMax,length.out=B)


  return(returnValue)
}

#' Print object of classe FlexCoDE
#'
#' @param objectCDE Object of the class "FlexCoDE", typically fitted used \code{\link{fitFlexCoDE}} beforehand
#'
#' @return returns information regarding the fitted model
#' @export
#'
print.FlexCoDE=function(objectCDE)
{
  if(class(objectCDE)!="FlexCoDE")
    stop("Object should be of class FlexCoDE")
  cat("FlexCoDE - Flexible Conditional Density Estimator \n \n")
  cat(paste("Regression Method Used:",class(objectCDE$regressionObject)),"\n")

  cat(paste("Best Number of Expansion Coefficients Sected:",(objectCDE$bestI)),"\n")

  cat(paste("Basis used:",objectCDE$system,"\n"))

  if(!is.null(objectCDE$estimatedRisk)) cat(paste("Estimated risk on validation set: ",objectCDE$estimatedRisk$mean," (",objectCDE$estimatedRisk$seBoot,")","\n",sep=""))
}
