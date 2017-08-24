# FlexCoDE

Implements *FlexCoDE*, Flexible Conditional Density Estimator, in R.

More on FlexCoDE: Izbicki, R.; Lee, A.B. [Converting High-Dimensional Regression to High-Dimensional Conditional Density Estimation](https://projecteuclid.org/euclid.ejs/1499133755). Eletronic Journal of Statistics, 2017


To install the package, run

```R
# install.packages("devtools")
devtools::install_github("rizbicki/FlexCoDE")
```

A simple example:

```R
# generate data
n=1000
d=10
data=matrix(NA,n,d+1)
data[,1:d]=matrix(rnorm(n*d),n,d)
data[,d+1]=data[,1]+rnorm(n,0,0.1)

# determine sample sizes
nTrain=round(0.7*n)
nValidation=round(0.25*n)
nTest=n-nTrain-nValidation

# split data
randomIndex=sample(1:n)
xTrain=data[randomIndex[1:nTrain],1:d]
xValidation=data[randomIndex[(nTrain+1):(nTrain+nValidation)],1:d]
xTest=data[randomIndex[(nTrain+nValidation+1):n],1:d]
zTrain=data[randomIndex[1:nTrain],d+1]
zValidation=data[randomIndex[(nTrain+1):(nTrain+nValidation)],d+1]
zTest=data[randomIndex[(nTrain+nValidation+1):n],d+1]

# Fit nearest neighbors FlexCoDE
fit=fitFlexCoDE(xTrain,zTrain,xValidation,zValidation,xTest,zTest,
            nIMax = 20,regressionFunction = regressionFunction.NN)
fit$estimatedRisk
print(fit)
plot(fit,xTest,zTest)

# Plot estimated curves evaluates on new test points, and
# compare with true conditional density
predictedValues=predict(fit,xTest,B=500)
plot(predictedValues$z,predictedValues$CDE[1,])
lines(predictedValues$z,dnorm(predictedValues$z,xTest[1,1],0.1),col=2)

# Fit sparse additive FlexCoDE
fit=fitFlexCoDE(xTrain,zTrain,xValidation,zValidation,xTest,zTest,
            nIMax = 30,regressionFunction = regressionFunction.SpAM)
fit$estimatedRisk
print(fit)
plot(fit,xTest,zTest)

# Plot estimated curves evaluates on new test points, and
# compare with true conditional density
predictedValues=predict(fit,xTest,B=500)
plot(predictedValues$z,predictedValues$CDE[1,])
lines(predictedValues$z,dnorm(predictedValues$z,xTest[1,1],0.1),col=2)


# Fit sparse additive FlexCoDE using 4 cores (i.e., using parallel computing)
fit=fitFlexCoDE(xTrain,zTrain,xValidation,zValidation,xTest,zTest,
            nIMax = 30,regressionFunction = regressionFunction.SpAM,
            regressionFunction.extra=list(nCores=4))
fit$estimatedRisk
print(fit)
plot(fit,xTest,zTest)

```
