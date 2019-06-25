#' FlexCoDE Fit Conditional Density Estimation via Regression
#'
#' @param xTrain Covariates x used to train the model (one observation per row)
#' @param zTrain Responses z used to train the model  (matrix with one column; one observation per row)
#' @param xValidation Covariates x used to tune the model (one observation per row; same number of columns as xTrain)
#' @param zValidation Responses z used to tune the model  (matrix with one column; one observation per row)
#' @param xTest Covariates x used to estimate risk of final model (one observation per row; same number of columns as xTrain). Default is NULL
#' @param zTest Responses z used to estimate risk of final model  (matrix with one column; one observation per row). Default is NULL
#' @param nIMax Maximum possible number of components of the series expansion (that is, the function will find the best I<nIMax). Default is 100
#' @param n_grid Number of grid points to evaluate estimated densities. Default is 1000
#' @param regressionFunction a function indicating which regression method will be used
#' to estimate the expansion coefficients.
#' Currently can be one of regressionFunction.NN, regressionFunction.NW,
#' regressionFunction.SpAM, regressionFunction.Series,
#' regressionFunction.Lasso, regressionFunction.Forest or
#' regressionFunction.XGBoost.
#' Type ?regressionFunction.XX to find out more about method XX.
#' @param regressionFunction.extra extra parameters to be sent to
#' regression function; see the regression you want
#'  to use to check what are the available parameters.
#'  The argument nCores which contains the number of cores to be used
#'  for parallel computing. Default is one.
#' @param system Basis for z. Current options are "Fourier", "Cosine" and "discrete". Default is "Fourier"
#' @param chooseDelta Should delta, the cutoff to remove spurious bumps, be chosen? Default is TRUE
#' @param deltaGrid Grid of threshold deltas (betwen 0 and 0.5). Default value is seq(0,0.4,0.05).
#' @param chooseSharpen Should alpha, the parameter to sharpen the final estimate, be chosen? Default is FALSE
#' @param sharpenGrid Grid of sharpen parameters alpha. Default value is seq(0.01,10,length.out = 20).
#' @param zMin Minimum value z assumes. Default is min(zTrain).
#' @param zMax Maximum value z assumes. Default is max(zTrain).
#' @param verbose Should we print what we are doing? Default is FALSE.
#'
#' @return Returns the fitted estimated conditional density, and object of the class FlexCoDE. The return value is an object with the following components:
#' \item{zMin, zMax}{Minimum and maximum value of z}
#' \item{nIMax}{Maximum number of expansion coefficients (user input). Default is minimum between 25 and number of training samples.}
#' \item{system}{Basis used for expanding the response}
#' \item{zTrain}{zTrain (user input)}
#' \item{xTrain}{xTrain (user input)}
#' \item{regressionObject}{Object with fitted regressions. Class and content depend on which regression method was chosen by user}
#' \item{errors}{Estimated errors for each value of I (number of expansion coefficients) using validation set}
#' \item{bestI}{Optimal number of I according to validation set}
#' \item{bestError}{Estimated error of model with bestI expansion terms according to validation set}
#' \item{bestDelta}{Optimal value of threshold delta according to validation set}
#' \item{bestAlpha}{Optimal value of alpha according to validation set}
#' \item{estimatedRisk}{(If user provides xTest and zTest) Estimated risk (error) according to test set)}
#'
#' @example ../testPackage.R
#'
#' @export
fitFlexCoDE=function(xTrain,zTrain,xValidation,zValidation,xTest=NULL,zTest=NULL,
                     nIMax=min(25,length(zTrain)),regressionFunction,
                     regressionFunction.extra=NULL,
                     system="Fourier",
                     deltaGrid=seq(0,0.45,length.out = 15),chooseDelta=TRUE,
                     sharpenGrid=seq(0.01,10,length.out = 20),chooseSharpen=FALSE,
                     zMin=NULL,zMax=NULL,n_grid=1000,verbose=FALSE)
{
  if(!is.matrix(xTrain))
    xTrain=as.matrix(xTrain)


  if(!is.matrix(xValidation))
    xValidation=as.matrix(xValidation)


  if(!is.null(xTest))
  {
    if(!is.matrix(xTest))
      xTest=as.matrix(xTest)

    if(!is.matrix(zTest))
      zTest=as.matrix(zTest)

  }


  if(!is.matrix(zTrain))
    zTrain=as.matrix(zTrain)

  if(!is.matrix(zValidation))
    zValidation=as.matrix(zValidation)



  objectCDE=NULL

  if (is.null(zMin)) {
    zMin <- min(zTrain)
  }
  objectCDE$zMin <- zMin

  if (is.null(zMax)) {
    zMax <- max(zTrain)
  }
  objectCDE$zMax <- zMax

  zTrain <- box_transform(zTrain, zMin, zMax)

  class(objectCDE)="FlexCoDE"
  objectCDE$verbose=verbose

  if(verbose) print("Transforming Response")
  responseFourier=calculateBasis(zTrain,nIMax,system)
  if(verbose) print("Fitting Regression Functions")
  regressionObject=regressionFunction(x=xTrain,responses=responseFourier,extra=regressionFunction.extra)
  objectCDE$nIMax=nIMax
  objectCDE$system=system
  objectCDE$zTrain=zTrain
  objectCDE$xTrain=xTrain
  objectCDE$regressionObject=regressionObject

  objectCDE$covariateNames=colnames(xTrain)

  rm(regressionObject,xTrain,zTrain,responseFourier)
  gc(verbose = FALSE)

  basisZValidation=calculateBasis(box_transform(zValidation, zMin, zMax),
                                  objectCDE$nIMax,objectCDE$system) # returns matrix length(z)xnIMax with the basis for z

  if(verbose) print("Tuning Number of Expansion Coefficients (I)")
  coefficientsXValidation=predict(objectCDE$regressionObject,xValidation)

  # Find level of expansion which minimizes the loss
  term1 <- 1/2*colMeans(coefficientsXValidation^2)
  term2 <- colMeans(coefficientsXValidation * basisZValidation[, 1:ncol(coefficientsXValidation), drop = FALSE])

  levels <- attr(basisZValidation, "levels")
  uniq_levels <- unique(levels)
  objectCDE$errors <- sapply(uniq_levels, function(level) {
    return(sum(term1[levels <= level] - term2[levels <= level]))
  })
  best_level <- uniq_levels[which.min(objectCDE$errors)]
  objectCDE$bestI <- max(which(levels <= best_level))
  objectCDE$bestError <- min(objectCDE$errors)

  objectCDE$n_grid=n_grid
  if (chooseDelta) {
    if (verbose) {
      print("Choosing optimal cutoff Delta")
    }

    objectCDE$bestDelta <- chooseDelta(objectCDE, xValidation, zValidation,
                                       deltaGrid, n_grid = n_grid)
  } else {
    objectCDE$bestDelta <- 0.0
  }

  if (chooseSharpen) {
    if (verbose) {
      print("Choosing optimal sharpen parameter alpha")
    }

    objectCDE$bestAlpha <- chooseSharpen(objectCDE, xValidation, zValidation,
                                         sharpenGrid, n_grid = n_grid)
  } else {
    objectCDE$bestAlpha <- 1.0
  }


  if(!is.null(xTest)&!is.null(zTest))
  {
    if(verbose) print("Estimating risk on test set")
    error <- estimateError(objectCDE, xTest, zTest, se = TRUE,n_grid=n_grid)
    objectCDE$estimatedRisk=error
  }

  if(objectCDE$bestI==objectCDE$nIMax)
    warning("\n the optimal I found was exactly nIMax; try increasing nIMax if you want to improve performance")


  return(objectCDE)
}

#' Fits FlexZBoost
#'
#' Wrapper for fitFlexCoDE for the special case of FlexZBoost
#'
#' @param xTrain Covariates x used to train the model (one observation per row)
#' @param zTrain Responses z used to train the model  (matrix with one column; one observation per row)
#' @param xValidation Covariates x used to tune the model (one observation per row; same number of columns as xTrain)
#' @param zValidation Responses z used to tune the model  (matrix with one column; one observation per row)
#' @param xTest Covariates x used to estimate risk of final model (one observation per row; same number of columns as xTrain). Default is NULL
#' @param zTest Responses z used to estimate risk of final model  (matrix with one column; one observation per row). Default is NULL
#' @param chooseSharpen Should alpha, the parameter to sharpen the final estimate, be chosen? Default is TRUE
#' @param ... additional arguments to fitFlexCoDE
#'
#' @return Returns the fitted estimated conditional density, and object of the class FlexCoDE. The return value is an object with the following components:
#' \item{zMin, zMax}{Minimum and maximum value of z}
#' \item{nIMax}{Maximum number of expansion coefficients (user input). Default is minimum between 25 and number of training samples.}
#' \item{system}{Basis used for expanding the response}
#' \item{zTrain}{zTrain (user input)}
#' \item{xTrain}{xTrain (user input)}
#' \item{regressionObject}{Object with fitted regressions. Class and content depend on which regression method was chosen by user}
#' \item{errors}{Estimated errors for each value of I (number of expansion coefficients) using validation set}
#' \item{bestI}{Optimal number of I according to validation set}
#' \item{bestError}{Estimated error of model with bestI expansion terms according to validation set}
#' \item{bestDelta}{Optimal value of threshold delta according to validation set}
#' \item{bestAlpha}{Optimal value of alpha according to validation set}
#' \item{estimatedRisk}{(If user provides xTest and zTest) Estimated risk (error) according to test set)}
#'
#'
#' @export
FlexZBoost=function(xTrain,zTrain,xValidation,zValidation,
                    xTest=NULL,zTest=NULL,chooseSharpen=TRUE,...)
{
  return(fitFlexCoDE(xTrain=xTrain,zTrain=zTrain,xValidation=xValidation,
                     zValidation=zValidation,
                     xTest=xTest,zTest=zTest,chooseSharpen=chooseSharpen,
                     regressionFunction = regressionFunction.XGBoost,
                     ...))
}

#' Choose threshold value to remove spurius bumps
#'
#' This function is typically not directly used by the user; it is used inside  \code{\link{fitFlexCoDE}}
#'
#' @param objectCDE An object of the class FlexCoDE with a fitted CDE,
#'   typically fitted used \code{\link{fitFlexCoDE}} beforehand
#' @param X Covariates used to validate (tune) the model (one x
#'   observation per row).
#' @param Z Responses used to validate (tune) the model (matrix with 1
#'   column). Each row corresponds to a row of the xValidation
#'   argument
#'
#' @return Value of delta for bump removal which minimizes CDE loss
chooseDelta <- function(objectCDE, X, Z, delta_grid = seq(0.0, 0.4, 0.05),n_grid=n_grid) {
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }

  preds <- predict(objectCDE, X, process = FALSE,B=n_grid)

  bin_size <- diff(preds$z)[1]
  estimates <- t(apply(preds$CDE, 1, function(xx) {
    return(normalize_density(bin_size, xx))
  }))

  errors <- rep(NA, length(delta_grid))
  for (ii in seq_along(delta_grid)) {
    new_estimates <- t(apply(estimates, 1, function(xx) {
      tmp <- remove_bumps(bin_size, xx, delta = delta_grid[ii])
      return(normalize_density(bin_size, tmp))
    }))

    errors[ii] <- cde_loss(new_estimates, preds$z, Z)
  }

  return(delta_grid[max(which.min(errors))])
}


#' Choose sharpen parameter of a conditional density estimator
#'
#' Chooses optimal alpha for the new conditional density
#' estimator f(z|x) <- f^alpha(z|x).
#'
#' @param objectCDE An object of the class FlexCoDE with a fitted CDE,
#'   typically fitted used \code{\link{fitFlexCoDE}} beforehand
#' @param X Covariates used to validate (tune) the model (one x
#'   observation per row).
#' @param Z Responses used to validate (tune) the model (matrix with 1
#'   column). Each row corresponds to a row of the xValidation
#'   argument
#' @param alpha_grid Grid of values of alpha
#'
#' @return Value of alpha which minimizes CDE loss
chooseSharpen <- function(objectCDE, X, Z, alpha_grid=seq(0.01,10,length.out = 20),n_grid=n_grid) {
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }

  preds <- predict(objectCDE, X, process = TRUE,B=n_grid)
  bin_size <- diff(preds$z)[1]


  errors <- rep(NA, length(alpha_grid))
  for (ii in seq_along(alpha_grid)) {
    new_estimates=(preds$CDE)^alpha_grid[ii]
    # reescale so that it is a proper density
    new_estimates=(new_estimates/rowSums(new_estimates))/bin_size
    errors[ii] <- cde_loss(new_estimates, preds$z, Z)
  }

  return(alpha_grid[max(which.min(errors))])
}


#' Calculate CDE loss
#'
#' @param pred a matrix of conditional density estimates; rows
#'   correspond to densities on z_grid and columns correspond to
#'   observations z_test.
#' @param z_grid a vector of grid points at which pred is evaluated.
#' @param z_test a vector of observed z-values.
#'
#' @return The estimated CDE loss for the predict estimates
cde_loss <- function(pred, z_grid, z_test) {
  z_grid <- as.matrix(z_grid)

  z_min <- apply(z_grid, 2, min)
  z_max <- apply(z_grid, 2, max)
  z_delta <- prod(z_max - z_min) / nrow(z_grid)

  integrals <- z_delta * sum(pred ^ 2) / nrow(pred)

  nn_ids <- cbind(1:nrow(z_test), FNN::knnx.index(z_grid, z_test, k = 1))
  likeli <- mean(pred[nn_ids])

  return(integrals - 2 * likeli)
}

#' Estimate error (risk) of FlexCoDE object via test set
#'
#' @param obj is an object of the class FlexCoDE typically fitted used
#'   \code{\link{fitFlexCoDE}}
#' @param x_test Covariates for the sample used to test the model (one
#'   observation per row)
#' @param z_test Response for the sample used to test the model (one
#'   observation per row)
#' @param se Boolean flag determining if bootstrap standard error are
#'   computed. Default is TRUE
#' @param n_boot Number of bootstrap samples used to calculate
#'   standard errors. Default is 500
#' @param n_grid Number of grid points to evaluate conditional
#'   density. Default is 500.
#'
#' @return Estimated error (with SE if se = TRUE)
#' @export
estimateError <- function(obj, x_test, z_test, se = TRUE, n_boot = 500,
                          n_grid = 1000) {
  x_test <- as.matrix(x_test)
  z_test <- as.matrix(z_test)

  pred_obj <- predict(obj, x_test, B=n_grid)
  cde_test <- pred_obj$CDE
  z_grid <- pred_obj$z

  loss <- cde_loss(cde_test, z_grid, z_test)

  if (!se) {
    return(loss)
  }

  # Bootstrap standard errors
  boot_losses <- replicate(n_boot, {
    boot_ids <- sample(nrow(z_test), replace = TRUE)

    cde_boot <- cde_test[boot_ids, , drop = FALSE]
    z_boot <- z_test[boot_ids, , drop = FALSE]

    return(cde_loss(cde_boot, z_grid, z_boot))
  })

  return(list(mean = loss,
              se_boot = sqrt(var(boot_losses))))
}

#' Evaluates the estimated  density of new observations (testing points) of a "FlexCoDE" object
#'
#' @param objectCDE Object of the class "FlexCoDE", typically fitted used \code{\link{fitFlexCoDE}} beforehand
#' @param xNew Matrix with nTest rows and same number of columns as xTrain, containing x's for which the estimates are desired.
#' @param B Number of point where f(z|x) will be evaluated (on the z scale). This will be equally spaced between zMin and zMax
#' @param predictionBandProb Either a number indicating the probability for the highest predictive density region desired  or FALSE if bands are not desired. Default is FALSE
#'
#' @return The return value is an object with the following components
#' \item{z}{Points where the density was evaluate}
#' \item{CDE}{Matrix with value of the density at points z. Each row corresponds to a different observation x (i-th row of CDE corresponds to i-th row of xTest).}
#' \item{th}{(If predictionBandProb is not FALSE) Threshold values for each estimated density. The region where estimated densities are above these values have the approximate coverage probability desired. See  \code{\link{plot.FlexCoDE}} for ploting these regions.}
#'
#' @example ../predict.R
#'
#' @export
#'
predict.FlexCoDE <- function(obj, xNew, B = NULL, predictionBandProb = FALSE, process = TRUE) {
  if (!is.matrix(xNew)) {
    xNew <- as.matrix(xNew)
  }
  if(is.null(B))
    B=obj$n_grid
  z_grid <- seq(0.0, 1.0, length.out = B)

  if (!is.null(obj$bestI)) {
    n_basis <- obj$bestI
  } else {
    n_basis <- obj$nIMax
  }

  coeff <- predict(obj$regressionObject, xNew, maxTerms = n_basis)

  z_basis <- calculateBasis(z_grid, n_basis, obj$system)

  estimates <- tcrossprod(coeff, z_basis)

  if (!is.null(obj$bestDelta)) {
    delta <- obj$bestDelta
  } else {
    delta <- 0.0
  }

  if (!is.null(obj$bestAlpha)) {
    alpha <- obj$bestAlpha
  } else {
    alpha <- 1.0
  }


  binSize <- 1 / (B + 1)

  if (process) {
    estimates <- t(apply(estimates, 1, function(xx) {
      return(post_process(binSize, xx, delta = delta,alpha=alpha))
    }))
  }

  estimates <- estimates / (obj$zMax - obj$zMin)

  if (!predictionBandProb) {
    return(list(CDE = estimates,
                z = seq(obj$zMin, obj$zMax, length.out = B)))
  }

  th <- matrix(NA, nrow(estimates), 1)
  for (ii in 1:nrow(estimates)) {
    th[ii] = .findThresholdHPD((obj$zMax - obj$zMin) / B,
                               estimates[ii, ], predictionBandProb)
  }

  return(list(CDE = estimates,
              z = seq(obj$zMin, obj$zMax, length.out = B),
              th = th))
}

#' Print object of classe FlexCoDE
#'
#' @param objectCDE Object of the class "FlexCoDE", typically fitted used \code{\link{fitFlexCoDE}} beforehand
#'
#' @return returns information regarding the fitted model
#'
#' @export
#'
print.FlexCoDE=function(objectCDE)
{
  if(class(objectCDE)!="FlexCoDE")
    stop("Object should be of class FlexCoDE")
  cat("FlexCoDE - Flexible Conditional Density Estimator \n \n \n")
  cat("####### Caracteristic of the fitted CDE:\n\n")
  cat(paste("Regression Method Used:",class(objectCDE$regressionObject)),"\n")

  cat(paste("Best Number of Expansion Coefficients Sected:",(objectCDE$bestI)),"\n")

  cat(paste("Basis used:",objectCDE$system,"\n"))

  if(!is.null(objectCDE$estimatedRisk)) cat(paste("Estimated risk on test set: ",objectCDE$estimatedRisk$mean," (se: ",objectCDE$estimatedRisk$seBoot,")","\n",sep=""))

  cat("\n")
  cat("####### Caracteristic of the fitted regression:\n\n")
  print(objectCDE$regressionObject,bestI=objectCDE$bestI,nameCovariates=objectCDE$covariateNames)

}


#' Plots examples of estimated densities together with real response
#'
#' @param objectCDE Object of the class "FlexCoDE", typically fitted used \code{\link{fitFlexCoDE}} beforehand
#' @param xTest Covariates x of the sample used to test the model (one observation per row)
#' @param zTest Response z of the sample used to test the model (one observation per row)
#' @param nPlots Number of desired densities to be ploted (which will be picked at random). Default is minimum between 8 and number of testing points
#' @param fontSize Font size of axis labels and legend
#' @param lineWidth Line width of the curves to be ploted
#' @param predictionBandProb Either a number indicating the probability for the highest predictive density region desired  or FALSE if bands are not desired. Default is FALSE
#' @param lineWidthPred Line width of the prediction bands to be ploted
#'
#' @return Plot with estimated densities
#'
#' @examples # See \code{\link{fitFlexCoDE}}
#'
#' @export
#'
plot.FlexCoDE=function(objectCDE,xTest,zTest,nPlots=min(nrow(xTest),9),fontSize=12,lineWidth=1,predictionBandProb=FALSE,lineWidthPred=0.6)
{
  if(!is.matrix(xTest))
    xTest=as.matrix(xTest)


  if(is.null(xTest))
    stop("Please provide xTest")


  if(is.null(zTest))
    stop("Please provide zTest")

  if(class(objectCDE)!="FlexCoDE")
    stop("objectCDE needs to be of class FlexCoDE")
  if(objectCDE$verbose)  print("Calculating predicted values")
  predictedValues=predict(objectCDE,xTest,B=objectCDE$n_grid,predictionBandProb=predictionBandProb)


  randomOrder=sample(1:nrow(xTest),nPlots,replace=FALSE)
  if(objectCDE$verbose) print("Creating plots")


  data=data.frame(x=predictedValues$z,y=predictedValues$CDE[randomOrder[1],],dataPoint=rep(1,length(predictedValues$z)),vertical=zTest[randomOrder[1]])
  if(nPlots>1)
  {
    for(i in 2:nPlots)
    {
      dataB=data.frame(x=predictedValues$z,y=predictedValues$CDE[randomOrder[i],],dataPoint=rep(i,length(predictedValues$z)),vertical=zTest[randomOrder[i]])
      data=rbind(data,dataB)
    }
  }

  g=ggplot2::ggplot(data,ggplot2::aes(x=x,y=y))+ggplot2::geom_line(size=lineWidth,color=2)+ggplot2::xlab("Response")+
    ggplot2::ylab("Estimated Density")+
    ggplot2::geom_vline(ggplot2::aes(xintercept=vertical),size=lineWidth)+
    ggplot2::theme(axis.title=ggplot2::element_text(size=fontSize,face="bold"))+ ggplot2::facet_wrap(~ dataPoint)
  print(g)

  if(predictionBandProb==FALSE)
    return()

  eps=0.35
  k=nrow(xTest)
  plot(x=1:k,y=zTest,main="",ylab="Prediction Region",cex.main=1.4,cex.axis=1.4,cex.lab=1.4,cex=1.5,col=1,xaxt="n",xlim=c(0.5,k+0.5),pch=16,ylim=c(objectCDE$zMin,objectCDE$zMax),xlab="Sample",bty="l")
  for(ii in 1:k)
  {
    whichLarger=predictedValues$CDE[ii,]>predictedValues$th[ii]
    runs=rle(whichLarger>0)
    nRuns=length(runs$values)

    cumulative=cumsum(runs$lengths)
    for(jj in 1:nRuns)
    {
      if(runs$values[jj]==TRUE)
      {
        if(jj==1)
        {
          lower=objectCDE$zMin
          upper=predictedValues$z[cumulative[jj]]
          lines(c(ii,ii),c(lower,upper),col=1,lwd=lineWidthPred)
          lines(c(ii-eps,ii+eps),c(lower,lower),col=1,lwd=lineWidthPred)
          lines(c(ii-eps,ii+eps),c(upper,upper),col=1,lwd=lineWidthPred)
          next;
        }
        #points(rep(ii,sum(whichLarger)),predicted$z[whichLarger],pch=18,cex=0.9,col=2)
        lower=predictedValues$z[cumulative[jj-1]]
        upper=predictedValues$z[cumulative[jj]]
        lines(c(ii,ii),c(lower,upper),col=1,lwd=lineWidthPred)

        lines(c(ii-eps,ii+eps),c(lower,lower),col=1,lwd=lineWidthPred)
        lines(c(ii-eps,ii+eps),c(upper,upper),col=1,lwd=lineWidthPred)
      }
    }
  }

  points(x=1:k,y=zTest,main="",ylab="Estimate",cex.main=1.4,cex.axis=1.4,cex.lab=1.4,cex=1.5,col=1,xaxt="n",xlim=c(0.5,k+0.5),pch=16,
         xlab="Sample")


}




#' Plots examples of estimated densities together with real response
#'
#' @param objectCDE_binded Object of the class "FlexCoDE_binded", typically obtained using \code{\link{bindFlexCoDE}} beforehand
#' @param xTest Covariates x of the sample used to test the model (one observation per row)
#' @param zTest Response z of the sample used to test the model (one observation per row)
#' @param nPlots Number of desired densities to be ploted (which will be picked at random). Default is minimum between 8 and number of testing points
#' @param fontSize Font size of axis labels and legend
#' @param lineWidth Line width of the curves to be ploted

#' @return Plot with estimated densities
#'
#' @examples # See \code{\link{bindFlexCoDE}}
#'
#' @export
plot.FlexCoDE_binded=function(objectCDE_binded,xTest,zTest,nPlots=min(nrow(xTest),9),fontSize=12,lineWidth=1)
{

  if(is.vector(xTest))
    xTest=as.matrix(xTest)

  if(is.null(xTest))
    stop("Please provide xTest")


  if(is.null(zTest))
    stop("Please provide zTest")

  if(class(objectCDE_binded)!="FlexCoDE_binded")
    stop("objectCDE_binded needs to be of class FlexCoDE_binded")
  if(objectCDE_binded[[1]]$verbose)  print("Calculating predicted values")

  predictedValues=list()
  for(b in 1:length(objectCDE_binded))
  {
    predictedValues[[b]]=predict(objectCDE_binded[[b]],xTest,B=500)
  }


  namesEstimators=sapply(objectCDE_binded, function(x)
    class(x$regressionObject))
  randomOrder=sample(1:nrow(xTest),nPlots,replace=FALSE)

  if(objectCDE_binded[[1]]$verbose)  print("Creating plots")

  x=c(sapply(predictedValues, function(x)x$z))
  y=c(sapply(predictedValues, function(x)x$CDE[randomOrder[1],]))
  data=data.frame(x=x,y=y,Estimator=as.factor(rep(namesEstimators,each=length(predictedValues[[1]]$z))),dataPoint=1,vertical=zTest[randomOrder[1]])
  if(nPlots>1)
  {
    for(i in 2:nPlots)
    {
      y=c(sapply(predictedValues, function(x)x$CDE[randomOrder[i],]))
      dataB=data.frame(x=x,y=y,Estimator=as.factor(rep(namesEstimators,each=length(predictedValues[[1]]$z))),dataPoint=i,vertical=zTest[randomOrder[i]])
      data=rbind(data,dataB)
    }
  }

  ggplot2::ggplot(data,ggplot2::aes(x=x,y=y,color=Estimator))+ggplot2::geom_line(size=lineWidth)+ggplot2::xlab("Response")+
    ggplot2::ylab("Estimated Density")+
    ggplot2::geom_vline(ggplot2::aes(xintercept=vertical),size=lineWidth)+
    ggplot2::theme(axis.title=ggplot2::element_text(size=fontSize,face="bold"))+ ggplot2::facet_wrap(~ dataPoint)+
    ggplot2::theme(legend.direction = "horizontal",legend.position = "top",legend.title=ggplot2::element_text(size=16,face="bold"),legend.text=ggplot2::element_text(size=fontSize),axis.title=ggplot2::element_text(size=fontSize,face="bold"))

}


#' Binds together objects of the class "FlexCoDE"
#'
#' @param objectCDE1 An object of the class FlexCoDE with a fitted CDE, typically fitted used \code{\link{fitFlexCoDE}} beforehand
#' @param objectCDE2 An object of the class FlexCoDE with a fitted CDE, typically fitted used \code{\link{fitFlexCoDE}} beforehand
#' @param ... other objects of the class FlexCoDE with a fitted CDE, typically fitted used \code{\link{fitFlexCoDE}} beforehand
#'
#'
#' @return list with all objects combined. Result is of the class "FlexCoDE_binded"
#' @example ../testPackageBind.R
#'
#' @export
#'
bindFlexCoDE=function(objectCDE1,objectCDE2,...)
{
  returnValue=append(list(objectCDE1,objectCDE2),list(...))
  class(returnValue)="FlexCoDE_binded"
  return(returnValue)
}



#' Finds best linear combination of several FlexCoDE estimates
#'
#' @param objectCDE_binded An object of the class FlexCoDE_binded with a fitted CDE, typically fitted used \code{\link{bindFlexCoDE}} beforehand
#' @param xValidation Covariates x used to validate (tune) the model (one x observation per row).
#' @param zValidation Responses z used to validate (tune) the model  (matrix with 1 column). Each row corresponds to a row of the xValidation argument
#' @param xTest Covariates x used to estimate risk of final model (one observation per row; same number of columns as xTrain). Default is NULL
#' @param zTest Responses z used to estimate risk of final model  (matrix with one column; one observation per row). Default is NULL
#' @return Returns an object of the class "combinedFlexCoDE" which contains the weights best linear combination of the input models, together with all fitted models
#'
#' @example ../testPackageCombined.R
#'
#' @export
combineFlexCoDE=function(objectCDE_binded,xValidation,zValidation,xTest=NULL,zTest=NULL)
{


  if(is.vector(xValidation))
    xValidation=as.matrix(xValidation)


  if(class(objectCDE_binded)!="FlexCoDE_binded")
    stop("Class of objectCDE_binded should be FlexCoDE_binded")
  predictedValues=list()
  for(b in 1:length(objectCDE_binded))
  {
    predictedValues[[b]]=predict(objectCDE_binded[[b]],xValidation,B=500)
  }

  grid=predictedValues[[1]]$z
  estimatesValidation=lapply(predictedValues, function(x)x$CDE)


  width=grid[2]-grid[1]
  nModels=length(estimatesValidation)

  B=matrix(0,nModels,nModels)
  for(i in 1:nModels)
  {
    for(j in 1:nModels)
    {
      B[i,j]=mean(width*rowSums(estimatesValidation[[i]]*estimatesValidation[[j]]))
    }
  }

  whichZ=apply(as.matrix(zValidation),1,function(x){
    which.min(abs(x-grid))
  })
  b=rep(NA,nModels)
  for(i in 1:nModels)
  {
    m=estimatesValidation[[i]]
    b[i]=mean(m[(1:nrow(m)) + nrow(m) * (whichZ - 1)])
  }

  weights=quadprog::solve.QP(Dmat=B, dvec=b, Amat=t(rbind(1,diag(nModels))), bvec=c(1,rep(0,nModels)), meq=1, factorized=FALSE)
  weights$solution[weights$solution<0]=0
  weights$solution=weights$solution/sum(weights$solution)



  returnValue=list(objectCDEs=objectCDE_binded,weights=weights$solution)
  class(returnValue)="combinedFlexCoDE"

  rm(objectCDE_binded)
  gc(verbose = FALSE)

  if(!is.null(xTest)&!is.null(zTest))
  {
    if(returnValue$objectCDEs[[1]]$verbose) print("Estimating risk on test set")
    error=estimateError(returnValue,xTest,zTest,se=TRUE)
    returnValue$estimatedRisk=error
  }


  return(returnValue)
}

#' Evaluates the estimated  density of new observations (testing points) of a "combinedFlexCoDE" object
#'
#' @param objectCDE Object of the class "combinedFlexCoDE", typically fitted used \code{\link{combineFlexCoDE}} beforehand
#' @param xNew Matrix with nTest rows and same number of columns as xTrain, containing x's for which the estimates are desired.
#' @param B Number of point where f(z|x) will be evaluated (on the z scale). This will be equally spaced between zMin and zMax
#'
#' @return The return value is an object with the following components
#' \item{z}{Points where the density was evaluate}
#' \item{CDE }{Matrix with value of the density at points z. Each row corresponds to a different observation x (i-th row of CDE corresponds to i-th row of xTest).}
#' @param predictionBandProb Either a number indicating the probability for the highest predictive density region desired  or FALSE if bands are not desired. Default is FALSE
#' @export
#'
#' @examples # See \code{\link{combineFlexCoDE}}
predict.combinedFlexCoDE=function(objectCombined,xNew,B=1000,predictionBandProb=FALSE)
{
  if(class(objectCombined)!="combinedFlexCoDE")
    stop("objectCombined should be of class combinedFlexCoDE")


  if(is.vector(xNew))
    xNew=as.matrix(xNew)

  predictedValues=list()
  for(b in 1:length(objectCombined$objectCDEs))
  {
    predictedValues[[b]]=predict(objectCombined$objectCDEs[[b]],xNew,B=500)
  }

  grid=predictedValues[[1]]$z
  estimatesValidation=lapply(predictedValues, function(x)x$CDE)

  predictedValuesFinal=matrix(0,nrow(estimatesValidation[[1]]),ncol(estimatesValidation[[1]]))
  for(b in 1:length(estimatesValidation))
  {
    predictedValuesFinal=predictedValuesFinal+estimatesValidation[[b]]*objectCombined$weights[b]
  }

  returnValue=NULL
  returnValue$CDE=predictedValuesFinal
  returnValue$z=grid


  if(predictionBandProb==FALSE)
    return(returnValue)


  th=matrix(NA,nrow(returnValue$CDE),1)
  for(i in 1:nrow(returnValue$CDE))
  {

    th[i]=.findThresholdHPD((objectCombined$zMax-objectCombined$zMin)/B,returnValue$CDE[i,],predictionBandProb)


  }

  returnValue$th=th
  return(returnValue)

  th=matrix(NA,nrow(returnValue$CDE),2)
  for(i in 1:nrow(returnValue$CDE))
  {
    interval=.findThresholdSymmetricMode((objectCombined$zMax-objectCombined$zMin)/B,
                                         returnValue$CDE[i,],
                                         predictionBandProb)
    intervalExtended=interval[1]:interval[2]
    for(k in 1:length(intervalExtended))
    {
      if(returnValue$CDE[i,intervalExtended][k]==0)
      {
        interval[1]=interval[1]+1
      } else {
        break;
      }
    }
    for(k in length(intervalExtended):1)
    {
      if(returnValue$CDE[i,intervalExtended][k]==0)
      {
        interval[2]=interval[2]-1
      } else {
        break;
      }
    }
    th[i,1]=returnValue$z[interval[1]]
    th[i,2]=returnValue$z[interval[2]]
  }
  returnValue$th=th
  return(returnValue)

}

#' Print object of classe combinedFlexCoDE
#'
#' @param objectCDE Object of the class "combinedFlexCoDE", typically fitted used \code{\link{combineFlexCoDE}} beforehand
#'
#' @return returns information regarding the fitted model
#' @export
#'
print.combinedFlexCoDE=function(objectCombined)
{
  if(class(objectCombined)!="combinedFlexCoDE")
    stop("objectCombined should be of class combinedFlexCoDE")


  cat("Object of class combinedFlexCoDE containing",length(objectCombined$weights),"fitted FlexCoDE regression estimators with weights \n ",objectCombined$weights,"\n respectively \n")
  cat("\n Estimators use the following regression methods respectively: \n")
  for(i in 1:length(objectCombined$weights))
  {
    cat(class(objectCombined$objectCDEs[[i]]$regressionObject),"\n")
  }

  if(!is.null(objectCombined$estimatedRisk))
  {
    cat(paste("Estimated risk on test set: ",objectCombined$estimatedRisk$mean," (se: ",objectCombined$estimatedRisk$seBoot,")","\n",sep=""))
  }

  cat("\n \n ############################## \n")
  cat("############################## \n \n")
  cat("\n Regression fits are the following: \n ")
  for(i in 1:length(objectCombined$weights))
  {
    cat("\n ############################## \n")
    cat("\n Fit ",i,":\n",sep = "")
    print(objectCombined$objectCDEs[[i]])
  }
}


#' Plots examples of estimated densities together with real response
#'
#' @param objectCombined Object of the class "combinedFlexCoDE", typically fitted used \code{\link{combineFlexCoDE}} beforehand
#' @param xTest Covariates x of the sample used to test the model (one observation per row)
#' @param zTest Response z of the sample used to test the model (one observation per row)
#' @param nPlots Number of desired densities to be ploted (which will be picked at random). Default is minimum between 8 and number of testing points
#' @param fontSize Font size of axis labels and legend
#' @param lineWidth Line width of the curves to be ploted
#' @param predictionBandProb Either a number indicating the probability for the highest predictive density region desired  or FALSE if bands are not desired. Default is FALSE
#' @param lineWidthPred Line width of the prediction bands to be ploted
#'
#'
#' @return Plot with estimated densities
#' @export
#'
#' @examples # See \code{\link{combineFlexCoDE}}
plot.combinedFlexCoDE=function(objectCombined,xTest,zTest,nPlots=min(nrow(xTest),9),fontSize=12,lineWidth=1,predictionBandProb=FALSE,lineWidthPred=0.6)
{

  if(is.vector(xTest))
    xTest=as.matrix(xTest)


  if(is.null(xTest))
    stop("Please provide xTest")


  if(is.null(zTest))
    stop("Please provide zTest")

  if(class(objectCombined)!="combinedFlexCoDE")
    stop("objectCDE needs to be of class combinedFlexCoDE")
  if(objectCombined$objectCDEs[[1]]$verbose)  print("Calculating predicted values")
  predictedValues=predict(objectCombined,xTest,B=500,predictionBandProb=predictionBandProb)

  randomOrder=sample(1:nrow(xTest),nPlots,replace=FALSE)
  if(objectCombined$objectCDEs[[1]]$verbose) print("Creating plots")

  data=data.frame(x=predictedValues$z,y=predictedValues$CDE[randomOrder[1],],dataPoint=rep(1,length(predictedValues$z)),vertical=rep(zTest[randomOrder[1]],length(predictedValues$z)))
  if(nPlots>1)
  {
    for(i in 2:nPlots)
    {
      dataB=data.frame(x=predictedValues$z,y=predictedValues$CDE[randomOrder[i],],dataPoint=rep(i,length(predictedValues$z)),vertical=rep(zTest[randomOrder[i]],length(predictedValues$z)))
      data=rbind(data,dataB)
    }
  }

  g=ggplot2::ggplot(data,ggplot2::aes(x=x,y=y))+
    ggplot2::geom_line(size=lineWidth,color=2)+ggplot2::xlab("Response")+
    ggplot2::ylab("Estimated Density")+
    ggplot2::geom_vline(ggplot2::aes(xintercept=vertical),
                        size=lineWidth,color=2)+
    ggplot2::theme(axis.title=ggplot2::element_text(size=fontSize,face="bold"))+ ggplot2::facet_wrap(~ dataPoint)
  print(g)

  if(predictionBandProb==FALSE)
    return()

  eps=0.35
  k=nrow(xTest)
  plot(x=1:k,y=zTest,main="",ylab="Prediction Region",cex.main=1.4,
       cex.axis=1.4,cex.lab=1.4,cex=1.5,col=1,xaxt="n",
       xlim=c(0.5,k+0.5),pch=16,ylim=c(objectCombined$zMin,objectCombined$zMax),
       xlab="Sample",bty="l")
  for(ii in 1:k)
  {
    whichLarger=predictedValues$CDE[ii,]>predictedValues$th[ii]
    runs=rle(whichLarger>0)
    nRuns=length(runs$values)

    cumulative=cumsum(runs$lengths)
    for(jj in 1:nRuns)
    {
      if(runs$values[jj]==TRUE)
      {
        if(jj==1)
        {
          lower=objectCombined$zMin
          upper=predictedValues$z[cumulative[jj]]
          lines(c(ii,ii),c(lower,upper),col=1,lwd=lineWidthPred)
          lines(c(ii-eps,ii+eps),c(lower,lower),col=1,lwd=lineWidthPred)
          lines(c(ii-eps,ii+eps),c(upper,upper),col=1,lwd=lineWidthPred)
          next;
        }
        #points(rep(ii,sum(whichLarger)),predicted$z[whichLarger],pch=18,cex=0.9,col=2)
        lower=predictedValues$z[cumulative[jj-1]]
        upper=predictedValues$z[cumulative[jj]]
        lines(c(ii,ii),c(lower,upper),col=1,lwd=lineWidthPred)

        lines(c(ii-eps,ii+eps),c(lower,lower),col=1,lwd=lineWidthPred)
        lines(c(ii-eps,ii+eps),c(upper,upper),col=1,lwd=lineWidthPred)
      }
    }
  }

  points(x=1:k,y=zTest,main="",ylab="Estimate",cex.main=1.4,cex.axis=1.4,cex.lab=1.4,cex=1.5,col=1,xaxt="n",
         xlim=c(0.5,k+0.5),pch=16,xlab="Sample")


}
