
setwd('/Users/yuwang/Desktop')
x = read.csv("./fuses.csv", header=F)$V1

KK = 2 # the number of components
w     = 1/2                         #Assign equal weight to each component to start with
mu    = rnorm(1, mean(x), sd(x))   #Random cluster centers randomly spread over the support of the data
lambda = rexp(1, rate=1)
sigma = sd(x)                       #Initial standard deviation

# Plot the initial guess for the density
#xx = seq(-8,11,length=200)
#yy = w*dnorm(xx, mu[1], sigma) + (1-w)*dnorm(xx, mu[2], sigma)
#plot(xx, yy, type="l", ylim=c(0, max(yy)), xlab="x", 
#     ylab="Initial density", lwd=2)
#points(x, rep(0,n), col=cc.true)

## The actual MCMC algorithm starts here
# Priors
aa  = rep(1,KK)  # Uniform prior on w
eta = 0          # Mean 0 for the prior on mu_k
tau = 1          # Standard deviation 5 on the prior for mu_l
dd  = 2
qq  = 1
a = 1            # Gamma prior with shape parameter and rate parameter
b = 1            # both equal to 1 on lambda
           
# Number of iterations of the sampler
rrr   = 10000
burn  = 1000


# Storing the samples
cc.out    = array(0, dim=c(rrr, n))
w.out     = rep(0, rrr)
mu.out    = array(0, dim=rrr)
sigma.out = rep(0, rrr)
lambda.out = rep(0, rrr)
logpost   = rep(0, rrr)

# MCMC iterations
for(s in 1:rrr){
  # Sample the indicators
  cc = rep(0,n)
  for(i in 1:n){
    v = rep(0,KK)
    v[1] = log(w) + dexp(x[i], rate=lambda, log=TRUE)  #Compute the log of the weights
    v[2] = log(1-w) + dnorm(x[i], mu, sigma, log=TRUE)  #Compute the log of the weights
    v = exp(v - max(v))/sum(exp(v - max(v)))
    cc[i] = sample(1:KK, 1, replace=TRUE, prob=v)
  }
  
  # Sample the weights
  w = rbeta(1, aa[1] + sum(cc==1), aa[2] + sum(cc==2))
  
  # Sample the means
  n2    = sum(cc==2)
  xsum2 = sum(x[cc==2])
  tau2.hat = 1/(n2/sigma^2 + 1/tau^2)
  mu.hat  = tau2.hat*(xsum2/sigma^2 + eta/tau^2)
  mu  = rnorm(1, mu.hat, sqrt(tau2.hat))
  
  # Sample the lambda
  a.star = a + sum(cc==1)
  b.star = b + sum(x[cc==1])
  lambda = rgamma(1, shape=a.star, rate=b.star)
  
  # Sample the variances
  dd.star = dd + n2/2
  qq.star = qq + sum((x[cc==2] - mu)^2)/2
  sigma = sqrt(rinvgamma(1, dd.star, qq.star))
  
  # Store samples
  cc.out[s,]   = cc
  w.out[s]     = w
  mu.out[s]   = mu
  lambda.out[s] = lambda
  sigma.out[s] = sigma
  for(i in 1:n){
    if(cc[i]==1){
      logpost[s] = logpost[s] + log(w) + dexp(x[i], lambda, log=TRUE)
    }else{
      logpost[s] = logpost[s] + log(1-w) + dnorm(x[i], mu, sigma, log=TRUE)
    }
  }
  logpost[s] = logpost[s] + dbeta(w, aa[1], aa[2],log = T)
  logpost[s] = logpost[s] + dnorm(mu, eta, tau, log = T)
  logpost[s] = logpost[s] + dgamma(lambda, shape=1, rate=1, log=T)
  logpost[s] = logpost[s] + log(dinvgamma(sigma^2, dd, 1/qq))
  if(s/500==floor(s/500)){
    print(paste("s =",s))
  }
}