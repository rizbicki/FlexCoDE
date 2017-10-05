# FlexCoDE

Implements *FlexCoDE*, Flexible Conditional Density Estimator, in R.

FlexCode  is a general-purpose methodology for converting any conditional mean point estimator of $z$ to a conditional {\em density} estimator $f(z \vert x)$, where $x$  represents the covariates. The key idea is to expand the unknown function $f(z \vert x)$ in an orthonormal basis $\{\phi_i(z)\}_{i}$:

$$f(z|x)=\sum_{i}\beta_{i }(x)\phi_i(z). $$

By the orthogonality property, the expansion coefficients are just conditional means 

$$\beta_{i }(x) =  \mathbb{E}\left[\phi_i(z)|x\right] \equiv \int f(z|x)   \phi_i(z) dz.$$

These coefficients can easily be estimated from data by regression. 

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

# FlexZBoost

FlexZBoost is a particular realization of FlexCode, where we use XGBoost  for the regression part, as these techniques scale well for massive data.  There is an additional tuning parameter
in FlexZBoost: an exponent $\alpha$ that we use to sharpen the computed density estimates $\widehat{f}(z|x)$, according to $\widetilde{f}(z|x) \propto \widehat{f}(z|x)^\alpha$.


A simple example of fitting flexZBoost using the data generated above:

```R
fit=flexZBoost(xTrain,zTrain,xValidation,zValidation,xTest,zTest,
            nIMax = 30)
fit$bestAlpha
fit$estimatedRisk
plot(fit,xTest,zTest)
```
