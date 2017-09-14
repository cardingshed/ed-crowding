ts.null.corr <- function(r,len, alpha.1,sigma.1=1,alpha.2,sigma.2=1) {
  # r = observed correlation coefficient between two series
  # alpha.1, sigma.1 are the AR(1) parameter and sd of series 1
  # alpha.2, sigma.2 are the AR(1) parameter and sd of series 2
  # len = length of series
  #
  # Function returns the probability of observing a correlation coefficient
  # equal to or greater than r when comparing two independent time-series
  # of lenght len with the input characteristics.
  #
  # check inputs - call to arima.sim requires stationarity and alpha non-zero
  if (alpha.1 <= 0 | alpha.1 >= 1) stop("alpha.1 must be > 0 and < 1")
  if (alpha.2 <= 0 | alpha.2 >= 1) stop("alpha.2 must be > 0 and < 1")
  # simulate 1000 AR(1) series with alpha.1 and sigma.1 parameters
  # each with len values. series.a is a lenx1000 matrix
  series.a <- replicate(1000,{
    as.vector(arima.sim(list(order=c(1,0,0),ar=alpha.1),n=len,sd=sigma.1)) })
  # similar simulation for alpha.2 and sigma.2
  series.b <- replicate(1000,{
    as.vector(arima.sim(list(order=c(1,0,0),ar=alpha.2),n=len,sd=sigma.2)) })
  # next create pairwise correlations of the series.a and series.b columns
  # e.g. cor.matrix[1,5] is the correlation between series.a[,1] and series.b[,5]
  # this value is repeated in cor.matrix[5,1] so must extract upper triangle only
  cor.matrix <- cor(series.a,series.b)
  cor.values <- data.frame(val= cor.matrix[upper.tri(cor.matrix)])
  # there are (1000^2 - 1000)/2 values in the upper triangle of cor.matrix
  probr <- sum(abs(cor.values)>=r)/499500
  # frequency a correlation >= r was obtained in the sample
  probr
}

### Example usage:
### (1) Compare two series without autocorrelation
# set.seed(1000)
# ts.null.corr(r=0.197,len=100,alpha.1=0.01,alpha.2=0.01)
# [1] 0.04912312
### (2) Interpret observed Pearson Correlation Coefficient of 0.5
### between two series with extensive autocorrelation
# ts.null.corr(r=0.5,len=100,alpha.1=0.7,alpha.2=0.8)
# [1] 0.003119119
### (3) Practical example with experimental values
# ts.null.corr(r=0.63,alpha.1=0.846,sigma.1=1,alpha.2=0.801,sigma.2=1,len=86)
# [1] 0.001813814
