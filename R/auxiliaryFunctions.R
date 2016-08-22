
#' Calculate Basis functions for new observations
#'
#' @param z elements where basis is going to be calculated
#' @param nIMax how many functions should be calculated
#' @param system Basis system to be used. Currently, either "cosine" or "Fourier"
#'
#' @return Test
calculateBasis=function(z,nIMax,system)
{
  if(system=="cosine")
  {
    basisZ=apply(as.matrix(1:(nIMax-1)),1,function(xx)sqrt(2)*cos(xx*pi*z))
    basisZ=cbind(rep(1,length(z)),basisZ)
    return(basisZ)
  }
  if(system=="Fourier")
  {
    if(nIMax==1)
      return(matrix(1,length(z),1))

    sinBasisZ=apply(as.matrix(1:round((nIMax)/2)),1,function(xx) sqrt(2)*sin(2*xx*pi*z))
    cosBasisZ=apply(as.matrix(1:round((nIMax)/2)),1,function(xx) sqrt(2)*cos(2*xx*pi*z))
    basisZ=matrix(NA,length(z),2*round((nIMax)/2))
    basisZ[,seq(from=1,length.out=dim(sinBasisZ)[2],by=2)]=sinBasisZ
    basisZ[,seq(from=2,length.out=dim(cosBasisZ)[2],by=2)]=cosBasisZ
    basisZ=cbind(rep(1,length(z)),basisZ)
    basisZ=basisZ[,1:nIMax]
    return(basisZ)
  } else if(system=="discrete")
  {
    eps=0.0001
    #basisZ=apply(as.matrix(1:(nZ)),1,function(xx)as.numeric((xx>z-1/2)&(xx<=z+1/2)))
    basisZ=apply(as.matrix(1:(nIMax)),1,function(xx)as.numeric(((xx-1)/nIMax+eps<z) & (z<(xx)/nIMax-eps)))

    return(basisZ)
  }    else
  {
    stop("System of Basis not known")
  }
}


.normalizeDensity=function(binSize,estimates,delta=0)
{
  # internal function of renormalization of density
  estimates=matrix(estimates,1,length(estimates))
  if(all(estimates<=0)) estimates=matrix(1,1,length(estimates))
  estimatesThresh=estimates
  #th=0
  th=1e-6
  estimatesThresh[estimatesThresh<th]=0
  
  if(sum(estimatesThresh)==0)
    return(matrix(ncol(estimatesThresh),1,ncol(estimatesThresh)))

  if(sum(as.vector(binSize*estimatesThresh))>1)
  {
    maxDensity=max(estimates)
    minDensity=0
    newXi=(maxDensity+minDensity)/2
    eps=1
    ii=1
    while(ii<=500)
    {

      #estimatesNew=apply(as.matrix(estimates),2,function(xx)max(0,xx-newXi))

      estimatesNew=estimates
      estimatesNew[estimatesNew-newXi<0]=0
      estimatesNew[estimatesNew-newXi>0]=(estimatesNew-newXi)[estimatesNew-newXi>0]


      area=sum(as.vector(binSize*estimatesNew))
      eps=abs(1-area)
      if(eps<0.001) break; # level found
      if(1>area) maxDensity=newXi
      if(1<area) minDensity=newXi
      newXi=(maxDensity+minDensity)/2
      ii=ii+1
    }
    #estimatesNew=apply(as.matrix(estimates),2,function(xx)max(0,xx-newXi))

    estimatesNew=estimates
    estimatesNew[estimatesNew-newXi<0]=0
    estimatesNew[estimatesNew-newXi>0]=(estimatesNew-newXi)[estimatesNew-newXi>0]

    runs=rle(as.vector(estimatesNew)>0)
    nRuns=length(runs$values)
    jj=1
    area=lower=upper=NULL
    if(nRuns>2)
    {
      for(ii in 1:nRuns)
      {
        if(runs$values[ii]==FALSE) next;
        whichMin=1
        if(ii>1)
        {
          whichMin=sum(runs$lengths[1:(ii-1)])
        }
        whichMax=whichMin+runs$lengths[ii]
        lower[jj]=whichMin # lower interval of component
        upper[jj]=whichMax # upper interval of component
        area[jj]=sum(as.vector(binSize*estimatesNew[whichMin:whichMax])) # total area of component
        jj=jj+1
      }

      delta=min(delta,max(area))
      for(ii in 1:length(area))
      {
        if(area[ii]<delta)
          estimatesNew[lower[ii]:upper[ii]]=0
      }
      estimatesNew=estimatesNew/(binSize*sum(as.vector(estimatesNew)))
    }

    return(estimatesNew)
  }
  estimatesNew=as.vector(1/binSize*estimatesThresh/sum(as.vector(estimatesThresh)))

  runs=rle(estimatesNew>0)
  nRuns=length(runs$values)
  jj=1
  area=lower=upper=NULL
  if(nRuns>2)
  {
    for(ii in 1:nRuns)
    {
      if(runs$values[ii]==FALSE) next;
      whichMin=1
      if(ii>1)
      {
        whichMin=sum(runs$lengths[1:(ii-1)])
      }
      whichMax=whichMin+runs$lengths[ii]
      lower[jj]=whichMin # lower interval of component
      upper[jj]=whichMax # upper interval of component
      area[jj]=sum(as.vector(binSize*estimatesNew[whichMin:whichMax])) # total area of component
      jj=jj+1
    }
    delta=min(delta,max(area))
    for(ii in 1:length(area))
    {
      if(area[ii]<delta)
        estimatesNew[lower[ii]:upper[ii]]=0
    }
    estimatesNew=estimatesNew/(binSize*sum(as.vector(estimatesNew)))
  }

  return(estimatesNew)

}


.findThresholdHPD=function(binSize,estimates,confidence)
{
  estimates=as.vector(estimates)
  maxDensity=max(estimates)
  minDensity=min(estimates)
  newCut=(maxDensity+minDensity)/2
  eps=1
  ii=1
  while(ii<=1000)
  {
    prob=sum(binSize*estimates*(estimates>newCut))
    eps=abs(confidence-prob)
    if(eps<0.0000001) break; # level found
    if(confidence>prob) maxDensity=newCut
    if(confidence<prob) minDensity=newCut
    newCut=(maxDensity+minDensity)/2
    ii=ii+1
  }
  return(newCut)
}

.findThresholdSymmetricMode=function(binSize,estimates,confidence)
{
  estimates=as.vector(estimates)
  mode=which.max(estimates)
  maxInverval=length(estimates)
  minInverval=1
  newCut=round(maxInverval/2)
  eps=1
  ii=1
  while(ii<=1000)
  {
    whichRegion=round(max(c(1,mode-newCut)):min(c(length(estimates),mode+newCut)))
    prob=sum(binSize*estimates[whichRegion])
    eps=abs(confidence-prob)
    if(eps<0.0000001) break; # level found
    if(confidence>prob) minInverval=newCut
    if(confidence<prob) maxInverval=newCut
    newCut=(maxInverval+minInverval)/2
    ii=ii+1
  }
  return(c(whichRegion[1],whichRegion[length(whichRegion)]))
}
