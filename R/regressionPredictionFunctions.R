predict.NN=function(object,xNew,maxTerms=NULL)
{
  if(class(object)!="NN")
    stop("Object has wrong class, should be NN")

  distanceNewTrain=fields::rdist(xNew,object$xTrain)

  if(!is.null(maxTerms))
  {
    maxTerms=min(maxTerms,length(object$bestNN))
  } else {
    maxTerms=length(object$bestNN)
  }


  predictedValidation=t(apply(distanceNewTrain,1,function(xx) {
    responsePredict=rep(NA,maxTerms)
    for(kk in 1:maxTerms)
    {
      nearest=sort(xx,index.return=T)$ix[1:object$bestNN[kk]]
      responsePredict[kk]=mean(object$responsesTrain[nearest,kk,drop=FALSE])
    }
    return(responsePredict)
  }))
  rm(distanceNewTrain)
  gc(verbose = FALSE)

  return(predictedValidation)
}

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
