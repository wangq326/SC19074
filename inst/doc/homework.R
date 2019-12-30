## -----------------------------------------------------------------------------
# Create two groups of points from two different populations
datatest <- rbind(matrix(rnorm(100, sd = 0.3), ncol = 2),
           matrix(rnorm(100, mean = 1, sd = 0.3), ncol = 2))
colnames(datatest) <- c("x", "y")

# Visualize the data
plot(datatest)

## -----------------------------------------------------------------------------
# Use K-means clustering to divide the whole group into two clusters
kmeanscluster <- kmeans(datatest, 2)

# Visualize the clusters
plot(datatest, col = kmeanscluster$cluster)
points(kmeanscluster$centers, col = 1 : 2, pch = 8, cex = 2)

## -----------------------------------------------------------------------------
y <- cor(state.x77)
library(knitr)
kable(head(y), format = "html")

## -----------------------------------------------------------------------------
# generate n r.v. from f(x) with sigma2
rRayleigh <- function(n, sigma2){
  j<-k<-0
  y <- numeric(n)
  while (k < n) {
  u <- runif(1)
  j <- j + 1
  x <- rgamma(1,shape = 2,scale = sigma2) # random variate from g
  if (exp(-(x^2-2*x+1)/(2*sigma2)) > u) {
    #we accept x
    k <- k + 1
    y[k] <- x
    }
  }
  hist(y, prob = TRUE)  
  y <- seq(0, 20, .002)
  lines(y, y*exp(-y^2/(2*sigma2))/sigma2) 
}

## -----------------------------------------------------------------------------
rRayleigh(100000,1)

## -----------------------------------------------------------------------------
rRayleigh(100000,4)

## -----------------------------------------------------------------------------
rRayleigh(100000,0.25)

## -----------------------------------------------------------------------------
rRayleigh(100000,16)

## -----------------------------------------------------------------------------
normlocationmixture <- function(n, p){
  x1 <- rnorm(n, 0, 1)
  x2 <- rnorm(n, 3, 1)
  r <- rbinom(n, 1, p)
  z <- r*x1+(1-r)*x2
  hist(z, probability = TRUE)
  z <- seq(-5, 10, .001)
  lines(z,p*dnorm(z,0,1)+(1-p)*dnorm(z,3,1))
}

## -----------------------------------------------------------------------------
normlocationmixture(10000, 0.75)

## -----------------------------------------------------------------------------
p <- seq(0.05,1,.05)
for (i in c(1:20)){
  normlocationmixture(10000,p[i])
  print(p[i])
}

## -----------------------------------------------------------------------------
Wdsampling1 <- function(Sigma, d, n){
  library(MASS)
  L <- chol(Sigma) # Sigma=LL'
  T <- matrix(0, nrow = d, ncol = d)
  for (i in (1:d)){
    T[i,i] <- sqrt(rchisq(1,n-i+1))
    for (j in 1:i-1){
      T[i,j] <- rnorm(1)
    }
  } 
  S <- L%*%T%*%t(T)%*%t(L)  #Bartlett's decomposition
  S
}

## -----------------------------------------------------------------------------
Wdsampling1(diag(3)+1, 3, 5)

## -----------------------------------------------------------------------------
Wdsampling2 <- function(Sigma, d, n){
  library(MASS)
  Id <- diag(d)  
  Md <- numeric(d)
  x <- mvrnorm(n, Md, Id)  #x[i,]~N(0,Id)
  L <- chol(Sigma)
  for (i in (1:n)){
      A <- x[i,]%*%t(x[i,])
  } # generate a wishart(Id,d,n) sample
  S <- L%*%A%*%t(L)
  S
}

## -----------------------------------------------------------------------------
Wdsampling2(diag(3)+1, 3, 5)

## -----------------------------------------------------------------------------
m <- 1e4   # number of samplings
t <- runif(m, min=0, max=pi/3)   # generate samples from U(0,pi/3)
theta.hat <- mean(sin(t)) * pi/3  # calculate theta.hat
print(c(theta.hat,0.5))  # check it vs the true value

## -----------------------------------------------------------------------------
n <- 1e4   # number of samplings
x <- runif(n/2, min=0, max=1)   # generate n/2 samples from U(0,1)
y <- 1-x
u <- exp(-x)/(1+x^2)
v <- exp(-y)/(1+y^2)
theta1.hat <- (mean(u)+mean(v))/2  
sd1.hat <- sd((u+v)/2)
# estimate it with antithetic variables 
x1 <- runif(n, min=0, max=1)
u1 <- exp(-x1)/(1+x1^2)
theta2.hat <- mean(u1)
sd2.hat <- sd(u1)
# estimate it without variance reduction
re.sd <- (sd2.hat-sd1.hat)/sd2.hat  #the approximate reduction in varianceas
print(re.sd)

## ----echo=TRUE----------------------------------------------------------------
m <- 10000
theta.hat <- se <- numeric(5)

g <- function(x) {
  exp(-x - log(1+x^2)) * (x > 0) * (x < 1)
}

x <- runif(m) #using f0
fg <- g(x)
theta.hat[1] <- mean(fg)
se[1] <- sd(fg)

x <- rexp(m, 1) #using f1
fg <- g(x) / exp(-x)
theta.hat[2] <- mean(fg)
se[2] <- sd(fg)

x <- rcauchy(m) #using f2
i <- c(which(x > 1), which(x < 0))
x[i] <- 2 #to catch overflow errors in g(x)
fg <- g(x) / dcauchy(x)
theta.hat[3] <- mean(fg)
se[3] <- sd(fg)

u <- runif(m) #f3, inverse transform method
x <- - log(1 - u * (1 - exp(-1)))
fg <- g(x) / (exp(-x) / (1 - exp(-1)))
theta.hat[4] <- mean(fg)
se[4] <- sd(fg)

u <- runif(m) #f4, inverse transform method
x <- tan(pi * u / 4)
fg <- g(x) / (4 / ((1 + x^2) * pi))
theta.hat[5] <- mean(fg)
se[5] <- sd(fg)

res <- rbind(theta=round(theta.hat,3), se=round(se,3))
colnames(res) <- paste0('f',0:4)
knitr::kable(res, format = "html",align='c')

## -----------------------------------------------------------------------------
M <- 10000
k <- 5
N <- 50

a <- (seq(k+1)-1)/k

g<-function(x)exp(-x)/(1+x^2)*(x>0)*(x<1)

kcof <- function(i){
  ans <- (1-exp(-1))/(exp(-(i-1)/k)-exp(-i/k))
  ans
}

st.im <- function(i){
  u <- runif(M/k)   # inverse transformation method
  v <- - log(exp(-(i-1)/k)-(exp(-(i-1)/k)-exp(-i/k))*u)
  fg <- g(v)/(kcof(i)*exp(-v)/(1-exp(-1)))
  fg
}
est <- matrix(0,N,2)
for (i in 1:N){
  for (j in 1:k){
    uu <- st.im(j)
    est[i,1] <- est[i,1]+mean(uu)
    est[i,2] <- est[i,2]+sd(uu)
  }
}
ans <- rbind(apply(est,2,mean),apply(est,2,sd))
colnames(ans) <- c('mean','sd')
library(knitr)
knitr::kable(ans,format='html')

## -----------------------------------------------------------------------------
M <- 10000
k <- 5
N <- 50

a <- numeric(k+1)
for (l in 2:k)
  a[l]=-log(1-(l-1)*(1-exp(-1))/k)
a[k+1] <- 1
a[1] <- 0
# divide the real line into k intervals

g<-function(x)exp(-x)/(1+x^2)*(x>0)*(x<1)
# integrated function

st.im <- function(lower,upper){
  u <- runif(M/k)   # inverse transformation method
  v <- -log(exp(-lower)-(1-exp(-1))*u/k)
  fg <- g(v)/(k*exp(-v)/(1-exp(-1)))
  fg
}
# samples from interval [lower,upper)

est <- matrix(0,N,2)
for (i in 1:N){
  for (j in 1:k){
    uu <- st.im(a[j],a[j+1])
    est[i,1] <- est[i,1]+mean(uu)
    est[i,2] <- est[i,2]+sd(uu)
  }
}
ans <- rbind(apply(est,2,mean),apply(est,2,sd))
colnames(ans) <- c('mean','sd')
library(knitr)
knitr::kable(ans,format='html')

## -----------------------------------------------------------------------------
ecp.chi2 <- function(n, m, v, alpha){
  # coverage probability of CI,sample size n,replicate m,confidence level 1-alpha
  ecp <- numeric(m)
  for(i in 1:m){
    x <- rchisq(n, df = v) # sample from chi(v)
    lcl <- mean(x) - qt(1-alpha/2, n-1) * sd(x) / sqrt(n) # lower conﬁdence limit
    ucl <- mean(x) + qt(1-alpha/2, n-1) * sd(x) / sqrt(n) # upper conﬁdence limit
    if(lcl <= v){
        if(v <= ucl){
          ecp[i] <- 1  
        }
    } # cp for a MC experiment replicate
  }
  cp <- mean(ecp)
  cp
}
ecp <- ecp.chi2(20, 1000, 2, 0.05)
ecp

## -----------------------------------------------------------------------------
n <- 20 
alpha <- .05 
UCL <- replicate(1000, expr = { 
  x <- rchisq(n, df = 2) 
  (n-1) * var(x) / qchisq(alpha, df = n-1) 
  })  
ecp.eg4 <- mean(UCL > 4) 

## -----------------------------------------------------------------------------
sk <- function(x) {
  #computes the sample skewness coeff. 
  xbar <- mean(x) 
  m3 <- mean((x - xbar)^3) 
  m2 <- mean((x - xbar)^2) 
  return(m3 / m2^1.5) 
} 

MC_sk <- function(n, m){
  #a MC experiment of m replicate,with n sample from population
  #estimate skewness of each replicate
  MC_skewness <- numeric(m)
  for (i in 1:m){
    replicate_MC <- rnorm(n)
    MC_skewness[i] <- sk(replicate_MC)
  }
  MC_skewness
}

level <- c(0.025, 0.05, 0.95, 0.975)
n <- 20
m <- 1000
q_sk <- quantile(MC_sk(n, m), level)
cv <- qnorm(level, mean = 0, sd = sqrt(6/n))
knitr::kable(t(cbind(q_sk,cv)), format = 'html', caption = 'Quantile')

sd_hat <- sqrt(level*(1-level)/(n*dnorm(q_sk,mean=0,sd=sqrt(6*(n-2)/(n+1)/(n+3)))^2)) # Compute the standard error of the estimates
sd_hat

## -----------------------------------------------------------------------------
set.seed(1107)
sk <- function(x) {
  #computes the sample skewness 
  xbar <- mean(x) 
  m3 <- mean((x - xbar)^3) 
  m2 <- mean((x - xbar)^2) 
  return(m3 / m2^1.5) 
} 
Powersktest1 <- function(a){
  n <- 20 #sample sizes
  cv <- qnorm(.975, 0, sqrt(6*(n-2)/(n+1)/(n+3))) 
  #crit. values for each n
  p.reject <- numeric(length(a)) #to store sim. results 
  m <- 10000 #num. repl. each sim.
  for (i in 1:length(a)) { 
    sktests <- numeric(m) #test decisions 
    for (j in 1:m) { 
      x <- rbeta(n,2,a[i]) #test decision is 1 (reject) or 0
      sktests[j] <- as.integer(abs(sk(x)) >= cv) 
      } 
    p.reject[i] <- mean(sktests) #proportion rejected 
    }
return(p.reject)
}
a <- seq(2.5,10,.1)
plot(cbind('alpha 2'=a,'power'=Powersktest1(a)))

## -----------------------------------------------------------------------------
Powersktest2 <- function(a2){
  a1 <- 2
  n <- 20 #sample sizes
  cv <- qnorm(.975, 0, sqrt(6*(n-2)/(n+1)/(n+3))) #crit. values for each n
  p.reject <- numeric(length(v)) #to store sim. results 
  m <- 10000 #num. repl. each sim.
  for (i in 1:length(v)) { 
    sktests <- numeric(m) #test decisions 
    for (j in 1:m) { 
      x <- rt(n,v[i]) #test decision is 1 (reject) or 0
      sktests[j] <- as.integer(abs(sk(x)) >= cv) 
      } 
    p.reject[i] <- mean(sktests) #proportion rejected 
    }
return(p.reject)
}

v <- seq(0.5,10,.1)
plot(cbind(v,'power'=Powersktest2(v)))

## -----------------------------------------------------------------------------
alpha <- .1 
n <- 30 
m <- 2500 
epsilon <- c(seq(0, .15, .01), seq(.15, 1, .05)) 
N <- length(epsilon) 
pwr <- numeric(N) #critical value for the skewness test 
cv <- qnorm(1-alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
for (j in 1:N) { 
  #for each epsilon 
  e <- epsilon[j] 
  sktests <- numeric(m) 
  for (i in 1:m) { 
    #for each replicate 
    a <- sample(c(1, 500), replace = TRUE, size = n, prob = c(1-e, e)) 
    x <- rbeta(n, a, a) 
    sktests[i] <- as.integer(abs(sk(x)) >= cv) 
    } 
  pwr[j] <- mean(sktests) } 
#plot power vs epsilon 
plot(epsilon, pwr, type = "b", xlab = bquote(epsilon), ylim = c(0,1)) 
abline(h = .1, lty = 3) 
se <- sqrt(pwr * (1-pwr) / m) 
#add standard errors 
lines(epsilon, pwr+se, lty = 3) 
lines(epsilon, pwr-se, lty = 3)

## -----------------------------------------------------------------------------
alpha <- .1 
n <- 30 
m <- 2500 
epsilon <- c(seq(0, .15, .01), seq(.15, 1, .05)) 
N <- length(epsilon) 
pwr <- numeric(N) #critical value for the skewness test 
cv <- qnorm(1-alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
for (j in 1:N) { 
  #for each epsilon 
  e <- epsilon[j] 
  sktests <- numeric(m) 
  for (i in 1:m) { 
    #for each replicate 
    v <- sample(c(1, 20), replace = TRUE, size = n, prob = c(1-e, e)) 
    x <- rt(n, v) 
    sktests[i] <- as.integer(abs(sk(x)) >= cv) 
    } 
  pwr[j] <- mean(sktests) } 
#plot power vs epsilon 
plot(epsilon, pwr, type = "b", xlab = bquote(epsilon), ylim = c(0,1)) 
abline(h = .1, lty = 3) 
se <- sqrt(pwr * (1-pwr) / m) 
#add standard errors 
lines(epsilon, pwr+se, lty = 3) 
lines(epsilon, pwr-se, lty = 3)

## -----------------------------------------------------------------------------
t1e <- numeric(3)
n <- 200 
m <- 1000 

## -----------------------------------------------------------------------------
pvalues <- replicate(m, expr = { 
  #simulate under alternative mu1 
  x <- rchisq(n, df = 1) 
  ttest <- t.test(x, alternative = "two.sided", mu = 1) 
  ttest$p.value     
  })
t1e[1] <- mean(pvalues <= .05) 
t1e[1]

## -----------------------------------------------------------------------------
pvalues <- replicate(m, expr = { 
  #simulate under alternative mu1 
  x <- runif(n, max = 2,min = 0) 
  ttest <- t.test(x, alternative = "two.sided", mu = 1) 
  ttest$p.value     
  })
t1e[2] <- mean(pvalues <= .05) 
t1e[2]

## -----------------------------------------------------------------------------
pvalues <- replicate(m, expr = { 
  #simulate under alternative mu1 
  x <- rexp(n, rate = 1) 
  ttest <- t.test(x, alternative = "two.sided", mu = 1) 
  ttest$p.value     
  })
t1e[3] <- mean(pvalues <= .05) 
t1e[3]

## -----------------------------------------------------------------------------
distribution <- c('chisq(1)','unif(0,2)','exp(1)')
rbind(distribution,'type 1 error' = t1e)

## -----------------------------------------------------------------------------
library(bootstrap)
data(scor)
plot(scor) #the scatter plots for each pair of test scores

## -----------------------------------------------------------------------------
cor(scor) #the sample correlation matrix

## -----------------------------------------------------------------------------
set.seed(9986)
library(boot) #for boot function

cor_12 <- function(x, i){ 
  #want correlation of columns 1 and 2 
  cor(x[i,1], x[i,2])
}
obj1 <- boot(data = scor, statistic = cor_12, R = 2000)

cor_34 <- function(x, i){ 
  cor(x[i,3], x[i,4])
}
obj2 <- boot(data = scor, statistic = cor_34, R = 2000)

cor_35 <- function(x, i){ 
  cor(x[i,3], x[i,5])
}
obj3 <- boot(data = scor, statistic = cor_35, R = 2000)

cor_45 <- function(x, i){ 
  cor(x[i,4], x[i,5])
}
obj4 <- boot(data = scor, statistic = cor_45, R = 2000)

bootstrap_cor <- c(mean(obj1$t),mean(obj2$t),mean(obj3$t),mean(obj4$t))
s <- c(sd(obj1$t),sd(obj2$t),sd(obj3$t),sd(obj4$t))
s_c <- cor(scor)
sample_cor=c(cor12=s_c[1,2],cor34=s_c[3,4],cor35=s_c[3,5],cor45=s_c[4,5])

rbind(sample_cor,bootstrap_cor,se_estimation=s)

## -----------------------------------------------------------------------------
library(boot) #for boot and boot.ci 
set.seed(666)

mean.boot <- function(dat, ind) { 
  #function to compute the statistic mean
  y <- dat[ind] 
  mean(y)
}

M=200
cr <- ml <- mr <- matrix(0,nrow = M,ncol = 3)
for (i in 1:M) {
  norm_s <- rnorm(100)
  obj5 <- boot(norm_s, statistic = mean.boot, R = 2000)
  a1 <- boot.ci(obj5, type = c("norm","basic","perc"))
  
  # the empirical coverage rates for the sample mean
  cr[i,1] <- ifelse(a1[["normal"]][2]<0&a1[["normal"]][3]>0,1,0)
  cr[i,2] <- ifelse(a1[["basic"]][4]<0&a1[["basic"]][5]>0,1,0)
  cr[i,3] <- ifelse(a1[["percent"]][4]<0&a1[["percent"]][5]>0,1,0)
  
  # the proportion of times that the conﬁdence intervals miss on the left
  ml[i,1] <- (a1[["normal"]][3]<0)
  ml[i,2] <- (a1[["basic"]][5]<0)
  ml[i,3] <- (a1[["percent"]][5]<0)
  
  # the porportion of times that the conﬁdence intervals miss on the right
  mr[i,1] <- (a1[["normal"]][2]>0)
  mr[i,2] <- (a1[["basic"]][4]>0)
  mr[i,3] <- (a1[["percent"]][4]>0)
  
}

coverage_rate <- c(norm=mean(cr[,1]),basic=mean(cr[,2]),percent=mean(cr[,3]))
p_on_left <- c(norm=mean(ml[,1]),basic=mean(ml[,2]),percent=mean(ml[,3]))
p_on_right <- c(norm=mean(mr[,1]),basic=mean(mr[,2]),percent=mean(mr[,3]))

rbind(coverage_rate,p_on_left,p_on_right)

## -----------------------------------------------------------------------------
library(boot) 
set.seed(9986)

sk.boot <- function(dat,ind) {
  # function to compute the statistic skewness 
  x <- dat[ind]
  xbar <- mean(x) 
  m3 <- mean((x - xbar)^3) 
  m2 <- mean((x - xbar)^2) 
  return(m3 / m2^1.5) 
} 


M=200
crnorm <- crchisq <- matrix(0,nrow = M,ncol = 3)
for (i in 1:M) {
  norm_s <- rnorm(100)
  obj6 <- boot(norm_s, statistic = sk.boot, R = 2000) # bootstrap
  a2 <- boot.ci(obj6, type = c("norm","basic","perc")) # bootstrap ci for 3 methods
  
  # the empirical coverage rates for the sample sk
  crnorm[i,1] <- (a2[["normal"]][2]<0&a2[["normal"]][3]>0)
  crnorm[i,2] <- (a2[["basic"]][4]<0&a2[["basic"]][5]>0)
  crnorm[i,3] <- (a2[["percent"]][4]<0&a2[["percent"]][5]>0)
}

for (i in 1:M) {
  chi5_s <- rchisq(100,df=5)
  obj7 <- boot(chi5_s, statistic = sk.boot, R = 2000)
  a3 <- boot.ci(obj7, type = c("norm","basic","perc"))
  
  # the empirical coverage rates for the sample sk
  crchisq[i,1] <- (a3[["normal"]][2]<4/sqrt(10)&a3[["normal"]][3]>4/sqrt(10))
  crchisq[i,2] <- (a3[["basic"]][4]<4/sqrt(10)&a3[["basic"]][5]>4/sqrt(10))
  crchisq[i,3] <- (a3[["percent"]][4]<4/sqrt(10)&a3[["percent"]][5]>4/sqrt(10))
}

norm_coverage_rate <- c(norm=mean(crnorm[,1]),basic=mean(crnorm[,2]),percent=mean(crnorm[,3]))
chi5_coverage_rate <- c(norm=mean(crchisq[,1]),basic=mean(crchisq[,2]),percent=mean(crchisq[,3]))

rbind(norm_coverage_rate,chi5_coverage_rate)

## -----------------------------------------------------------------------------
library(bootstrap)
data(scor)

prop.var.1pc <- function(dat, ind){
  # the proportion of variance explained by the ﬁrst principal component
  sigma <- cov(dat[ind,])
  lambda <- eigen(sigma)$values
  theta <- lambda[1]/sum(lambda)
  return(theta)
}

jack.bias.and.se <- function(dat, theta){
  # function to get the jackknife estimates of bias and standard error of theta
  # dat is the sample data;theta is the interseted parameter
  n <- dim(dat)[1]
  theta.hat <- theta(dat, 1:n)
   
  theta.jack <- numeric(n) 
  for(i in 1:n){ 
    theta.jack[i] <- theta(dat, (1:n)[-i])
  } 
  bias.jack <- (n-1)*(mean(theta.jack)-theta.hat) 
  se.jack <- sqrt((n-1)*mean((theta.jack-theta.hat)^2))
  round(c(original=theta.hat,bias.jack=bias.jack,se.jack=se.jack),3)
}

jack.bias.and.se(scor,prop.var.1pc)

## -----------------------------------------------------------------------------
library(DAAG)
attach(ironslag) 
a <- seq(10, 40, .1) #sequence for plotting fits

####par(mfrow=c(2,2))

L1 <- lm(magnetic ~ chemical) 
plot(chemical, magnetic, main="Linear", pch=16) 
yhat1 <- L1$coef[1] + L1$coef[2] * a 
lines(a, yhat1, lwd=2)

L2 <- lm(magnetic ~ chemical + I(chemical^2)) 
plot(chemical, magnetic, main="Quadratic", pch=16) 
yhat2 <- L2$coef[1] + L2$coef[2] * a + L2$coef[3] * a^2 
lines(a, yhat2, lwd=2)

L3 <- lm(log(magnetic) ~ chemical) 
plot(chemical, magnetic, main="Exponential", pch=16) 
logyhat3 <- L3$coef[1] + L3$coef[2] * a 
yhat3 <- exp(logyhat3) 
lines(a, yhat3, lwd=2)

L4 <- lm(magnetic ~ chemical + I(chemical^2) + I(chemical^3)) 
plot(chemical, magnetic, main="Cubic polynomial", pch=16) 
yhat4 <- L4$coef[1] + L4$coef[2] * a  + L4$coef[3] * a^2 + L4$coef[4] * a^3 
lines(a, yhat4, lwd=2)

## -----------------------------------------------------------------------------
n <- length(magnetic) 
e1 <- e2 <- e3 <- e4 <- numeric(n)
# for n-fold cross validation 
# fit models on leave-one-out samples 
for (k in 1:n) { 
  y <- magnetic[-k] 
  x <- chemical[-k]
  J1 <- lm(y ~ x) 
  yhat1 <- J1$coef[1] + J1$coef[2] * chemical[k] 
  e1[k] <- magnetic[k] - yhat1

  J2 <- lm(y ~ x + I(x^2)) 
  yhat2 <- J2$coef[1] + J2$coef[2] * chemical[k] + J2$coef[3] * chemical[k]^2 
  e2[k] <- magnetic[k] - yhat2

  J3 <- lm(log(y) ~ x) 
  logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[k] 
  yhat3 <- exp(logyhat3) 
  e3[k] <- magnetic[k] - yhat3

  J4 <- lm(y ~ x + I(x^2) + I(x^3)) 
  yhat4 <- J4$coef[1] + J4$coef[2] * chemical[k] + J4$coef[3] * chemical[k]^2 + J4$coef[4] * chemical[k]^3
  e4[k] <- magnetic[k] - yhat4
  } 
MSE <- c(Linear=mean(e1^2), Quadratic=mean(e2^2), Exponential=mean(e3^2), Cubic_ploy=mean(e4^2)) 
TSS <- var(y)*n
SSE <- MSE*n
R2 <- 1-SSE/TSS
p <- c(1,2,1,3)
ADR2 <- 1-(1-R2)*(n-1)/(n-p-1)

MSE <- signif(MSE,3)
R2 <- signif(R2,3)
ADR2 <- signif(ADR2,3)
rbind(MSE, R2, ADR2)

## -----------------------------------------------------------------------------
# count5 statistic
count5stat <- function(x, y) { 
  X <- x - mean(x) 
  Y <- y - mean(y) 
  outx <- sum(X > max(Y)) + sum(X < min(Y)) 
  outy <- sum(Y > max(X)) + sum(Y < min(X)) # return 1 (reject) or 0 (do not reject H0) 
  return(max(c(outx, outy)))
  }

# generate samples with two different sample sizes
x <- rnorm(30,1,1)
y <- rnorm(25,1,1)


R <- 999
z <- c(x, y)
K <- 1:(length(x)+length(y))
n<-length(x)
set.seed(12345)
reps <- numeric(R)
t0 <- count5stat(x,y)
for (i in 1:R){ 
  # permutation R times
  k <- sample(K, size = n, replace = FALSE) 
  x1 <- z[k]
  y1 <- z[-k] 
  reps[i] <- count5stat(x1,y1) 
} 
p <- mean(abs(c(t0, reps)) >= abs(t0)) 
print(p)

## -----------------------------------------------------------------------------
library(Ball)
library(boot)
library(bootstrap)
library(MASS)

## -----------------------------------------------------------------------------
dCov <- function(x, y){
  x <- as.matrix(x)
  y <- as.matrix(y)
  n <- nrow(x)
  m <- nrow(y)
  if (n != m || n < 2) stop("Sample sizes must agree")
  if (! (all(is.finite(c(x, y))))) stop("Data contains missing or infinite values")
  Akl <- function(x) {
    d <- as.matrix(dist(x))
    m <- rowMeans(d)
    M <- mean(d)
    a <- sweep(d, 1, m)
    b <- sweep(a, 2, m)
    b + M
    }
  A <- Akl(x)
  B <- Akl(y)
  sqrt(mean(A * B))
}
ndCov2 <- function(z, ix, dims) {
  #dims contains dimensions of x and y
  p <- dims[1]
  q <- dims[2]
  d <- p + q
  x <- z[ , 1:p] #leave x as is
  y <- z[ix, -(1:p)] #permute rows of y
  return(nrow(z) * dCov(x, y)^2)
}

## -----------------------------------------------------------------------------
mu <- c(0,0)
I <- diag(1,2)
m <- 100
set.seed(12345)
R <- 99

pow <- function(n,model){
  p.values <- matrix(NA,m,2)
  for(i in 1:m){
    x <- mvrnorm(n,mu,I)
    e <- mvrnorm(n,mu,I)
    if(model==1) y <- x/4+e
    if(model==2) y <- x/4*e
    z<-cbind(x,y)
    boot.obj <- boot(data = z, statistic = ndCov2, R = R,sim = "permutation", dims = c(2, 2))
    tb <- c(boot.obj$t0, boot.obj$t)
    p.values[i,1] <- mean(tb>=tb[1])
    p.values[i,2] <- bcov.test(x,y,R=R,seed=i*123)$p.value
    }
  alpha <- 0.05;
  pow2 <- colMeans(p.values<alpha)
  return(pow2)
}

## -----------------------------------------------------------------------------
N <- c(10,20,30,50,75,100)
power1 <- matrix(0,6,2)
power2 <- matrix(0,6,2)
for (i in 1:6) {
  power1[i,] <- pow(N[i],1)
  power2[i,] <- pow(N[i],2)
}
plot(N,power1[,1],type = "l",col = "black",ylab = "power",ylim = c(0,1),main = "Power Comparison : Y=X/4+e")
lines(N,power1[,2],col = "red")
legend("bottomright",legend=c("Ball covariance","Distance correlation"),
       col=c("red","black"),lty=1,lwd=1)  

plot(N,power2[,1],type = "l",col = "black",ylab = "power",ylim = c(0,1),main = "Power Comparison : Y=X/4*e")
lines(N,power2[,2],col = "red")
legend("bottomright",legend=c("Ball covariance","Distance correlation"),
       col=c("red","black"),lty=1,lwd=1)  

## -----------------------------------------------------------------------------
lap.f <- function(x){
  # prop density function of Laplace distribution
  return(1/2*exp(-abs(x)))
}

## -----------------------------------------------------------------------------
random.walk.Me <- function(sigma, x0, N){ 
  # function to generate a random walk metropolis chain
  x <- numeric(N) 
  x[1] <- x0 
  u <- runif(N) 
  k <- 0 
  for (i in 2:N){ 
    y <- rnorm(1, x[i-1], sigma) 
    if (u[i] <= (lap.f(y)/lap.f(x[i-1]))){
      x[i] <- y 
    }else{ 
      x[i] <- x[i-1] 
      k <- k + 1 } 
    } 
  return(list(x=x, k=k)) 
} 

N <- 2000 
sigma <- c(.05, .5, 2, 16)
x0 <- 25 
rw1 <- random.walk.Me(sigma[1], x0, N) 
rw2 <- random.walk.Me(sigma[2], x0, N) 
rw3 <- random.walk.Me(sigma[3], x0, N) 
rw4 <- random.walk.Me(sigma[4], x0, N)

## -----------------------------------------------------------------------------
inverse.F <- function(p){
  if(p>=0.5){
     perc <- -log(1-2*abs(p-0.5))
  }else{perc <- log(1-2*abs(p-0.5))}
  return(perc) 
}
perc1 <- inverse.F(0.025)
perc2 <- inverse.F(0.975)

## -----------------------------------------------------------------------------
rw.k <- c(rw1$k, rw2$k, rw3$k, rw4$k)
rate.acceptance <- (N-rw.k)/N
rbind(sigma,rate.acceptance)

## -----------------------------------------------------------------------------
x <- runif(1000,0.1,100)
le <- log(exp(x))
el <- exp(log(x))
sum(x==le) # number of x=log(exp(x))
sum(x==el) # number of x=exp(log(x))
sum(le==el)# number of log(exp(x))=exp(log(x))

## -----------------------------------------------------------------------------
eq <- numeric(3)
for (i in 1:1000) {
  eq[1] <- eq[1] + all.equal(x[i],le[i])
  eq[2] <- eq[2] + all.equal(x[i],el[i])
  eq[3] <- eq[3] + all.equal(el[i],le[i])
}
print(eq)

## -----------------------------------------------------------------------------
cupper <- function(k,a){
  return(sqrt(a^2*k/(k+1-a^2)))
}
f1 <- function(u){
  (1+u^2/(k-1))^(-k/2)
}
f2 <- function(u){
  (1+u^2/k)^(-(k+1)/2)
}
sol1 <- function(a){
  # the toot of sol1 is A(k)
  2*gamma(k/2)/(sqrt(pi*(k-1))*gamma((k-1)/2))*integrate(f1,0,cupper(k-1,a))$value-2*gamma((k+1)/2)/(sqrt(pi*k)*gamma(k/2))*integrate(f2,0,cupper(k,a))$value
}   

kt <- c(4:25,100)
n <- length(kt)
A1 <- A2 <- numeric(n)
for (i in 1:n) {
  k <- kt[i]
  A1[i] <- uniroot(sol1,c(0.5,sqrt(k)/2+1))$root
}

## Exercise 11.4
sol2 <- function(a){
   # the toot of sol2 is A(k)
   1-pt(cupper(k-1,a),k-1)-1+pt(cupper(k,a),k)
  }
for (i in 1:n) {
   k <- kt[i]
   A2[i] <- uniroot(sol2,c(1e-5,sqrt(k)-1e-5))$root
 }

cbind(df.t=kt,root.ex11.5=A1,root.ex11.4=A2)

## -----------------------------------------------------------------------------
set.seed(999)
library(stats4)

dll <- function(x){
  # root of dll is the p and q that maximum loglikeli
  t <- c(n.ob[5],n.ob[6],n.ob[3],n.ob[1]-n.ob[5],n.ob[2]-n.ob[6],n.ob[4])
  r <- 1-sum(x)
  f1 <- 2*t[1]/x[1]-2*t[3]/r+t[6]/x[1]+t[4]/x[1]-t[4]/r-t[5]/r
  f2 <- 2*t[2]/x[2]-2*t[3]/r+t[6]/x[2]-t[4]/r+t[5]/x[2]-t[5]/r
  c(F1=f1,F2=f2)
}

loglikeli <- function(p){
  # loglikelihood function
  return(2*n.ob[5]*log(p[1])+2*n.ob[6]*log(p[2])+2*n.ob[3]*log(1-p[1]-p[2])+(n.ob[1]-n.ob[5])*log(2*p[1]*(1-p[1]-p[2]))+(n.ob[2]-n.ob[6])*log(2*p[2]*(1-p[1]-p[2]))+n.ob[4]*log(2*p[1]*p[2]))
}

N <- 10000

theta <- matrix(0,N,2)
ll <- numeric(N)
n.ob <- c(28,24,41,70,NaN,NaN)
n.ob[5] <- sample(1:n.ob[1],1)
n.ob[6] <- sample(1:n.ob[2],1)
n <- sum(n.ob[1:4])
library(rootSolve)

for (i in 1:N) {
  theta[i,] <- multiroot(dll,start=c(0.1,0.2))$root
  ll[i] <- loglikeli(theta[i,])
  w <- c(theta[i,],1-sum(theta[i,]))
  o <- numeric(6)
  o[1] <- n*w[1]^2
  o[2] <- n*w[2]^2
  o[3] <- n*w[3]^2
  o[4] <- n*2*w[1]*w[3]
  o[5] <- n*2*w[2]*w[3]
  o[6] <- n*2*w[1]*w[2]
  n.ob <- c(o[1]+o[4],o[2]+o[5],o[3],o[6],o[1],o[2])
}

print(theta[9990:10000,])

index <- 1:100
plot(index,exp(ll[index]),ylab="max-likelihood",type="l")
index <- (N-1000):N
plot(index,exp(ll[index]),ylab="max-likelihood",type="l")

## -----------------------------------------------------------------------------
data("mtcars")
formulas <- list(mpg~disp, mpg~I(1/disp), mpg~disp+wt,mpg~I(1/disp)+wt)
lm.loop <- list()
for (i in 1:4) {
  lm.loop[[i]] <- lm(formulas[[i]],data = mtcars)
}
lm.lapply <- lapply(formulas, lm,data=mtcars)
lm.loop
lm.lapply

## -----------------------------------------------------------------------------
boot.sample <- function(i){
  rows <- sample(1:nrow(mtcars), rep = TRUE) 
  mtcars[rows, ]
}

lm.loop2 <- list()
for (i in 1:10) {
  lm.loop2[[i]] <- lm(mpg~disp,data = boot.sample(i))
}

lm.lapply2 <- lapply(lapply(1:10, boot.sample),lm,formula=mpg~disp)

lm.loop2
lm.lapply2

## -----------------------------------------------------------------------------
rsq <- function(mod) summary(mod)$r.squared

r2.ex3.loop <- lapply(lm.loop, rsq)
r2.ex3.lapply <- lapply(lm.lapply,rsq)

r2.ex4.loop <- lapply(lm.loop2, rsq)
r2.ex4.lapply <- lapply(lm.lapply2, rsq)

r2.ex3 <- cbind(model=as.character(formulas),r2.ex3.loop,r2.ex3.lapply)
r2.ex4 <- cbind(r2.ex4.loop,r2.ex4.lapply)

r2.ex3
r2.ex4

## -----------------------------------------------------------------------------
set.seed(123)
trials <- replicate(100,t.test(rpois(10,10),rpois(7,10)),simplify = FALSE) 

# Use sapply() and an anonymous function to extract the p-value from every trial.
pv1 <- sapply(trials, function(x) x$p.value)
# Extra challenge: get rid of the anonymous function by using [[ directly.
pv2 <- sapply(trials,"[[",3)

cbind(pv1,pv2)

## -----------------------------------------------------------------------------
library(parallel)
mcsapply <- function(X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE, chunk.size = NULL){
  # 计算可用线程数，并设置并行使用线程数
  no_cores <- detectCores() - 1
  # 初始化
  cl <- makeCluster(no_cores)
  ans <- parSapply(cl, X, FUN, ..., simplify = TRUE,
            USE.NAMES = TRUE, chunk.size = NULL)
  stopCluster(cl)
  return(ans)
}


# example
mcsapply(trials,function(x) x$p.value)

## -----------------------------------------------------------------------------
N <- 2000 
sigma <- c(.05, .5, 2, 16)
x0 <- 25

## -----------------------------------------------------------------------------
library(Rcpp)

#dir_cpp <- '../Rcpp/' # Can create source file in Rstudio 
#sourceCpp("rwme.cpp")

sourceCpp(code = '
          #include <Rcpp.h>
          using namespace Rcpp;
          
          //[[Rcpp::export]]
          NumericMatrix rwme(double sigma,double xo,int N){
          NumericMatrix x(N,2);
          x(0,0) = xo; 
          x(1,0) = 1;
          NumericVector u = runif(N); 
          for (int i = 1; i < N ;i++){ 
            double y = as<double>(rnorm(1, x(i-1,0), sigma));
            double t = exp(-abs(y))/exp(-abs(x(i-1,0)));
            if (u[i] > t){
              x(i,0) = x(i-1,0); 
              x(i,1) = 0;
            }
            else{ 
              x(i,0) = y;
              x(i,1) = 1;} 
          };
          return x;
          }')



## -----------------------------------------------------------------------------
lap.f <- function(x){
  # prop density function of Laplace distribution
  return(1/2*exp(-abs(x)))
}
randomwalk <- function(sigma, x0, N){ 
  # function to generate a random walk metropolis chain
  x <- matrix(0,N,2) 
  x[1,1] <- x0 
  u <- runif(N) 
  for (i in 2:N){ 
    y <- rnorm(1, x[i-1], sigma) 
    if (u[i] <= (lap.f(y)/lap.f(x[i-1]))){
      x[i,1] <- y 
      x[i,2] <- 1
    }else{ 
      x[i,1] <- x[i-1] 
      x[i,2] <- 0} 
    } 
  return(x) 
} 

## -----------------------------------------------------------------------------
rw1 <- randomwalk(sigma[1], x0, N) 
rw2 <- randomwalk(sigma[2], x0, N) 
rw3 <- randomwalk(sigma[3], x0, N) 
rw4 <- randomwalk(sigma[4], x0, N)

## -----------------------------------------------------------------------------
rw5 <- rwme(sigma[1], x0, N)
rw6 <- rwme(sigma[2], x0, N) 
rw7 <- rwme(sigma[3], x0, N) 
rw8 <- rwme(sigma[4], x0, N)

## -----------------------------------------------------------------------------
acR <- c(sum(rw1[,2]==1),sum(rw2[,2]==1),sum(rw3[,2]==1),sum(rw4[,2]==1))
acC <- c(sum(rw5[,2]==1),sum(rw6[,2]==1),sum(rw7[,2]==1),sum(rw8[,2]==1))
rbind(sigma,acR,acC)

## -----------------------------------------------------------------------------
####par(mfrow=c(2,2))
qqplot(rw1[rw1[,2]==1,1],rw5[rw5[,2]==1,1],xlab = "rwR",ylab = "rwC",main = expression("Q-Q plot for" ~~ {sigma==0.05}))
qqplot(rw2[rw2[,2]==1,1],rw6[rw6[,2]==1,1],xlab = "rwR",ylab = "rwC",main = expression("Q-Q plot for" ~~ {sigma==0.5}))
qqplot(rw3[rw3[,2]==1,1],rw7[rw7[,2]==1,1],xlab = "rwR",ylab = "rwC",main = expression("Q-Q plot for" ~~ {sigma==2}))
qqplot(rw4[rw4[,2]==1,1],rw8[rw8[,2]==1,1],xlab = "rwR",ylab = "rwC",main = expression("Q-Q plot for" ~~ {sigma==16}))

## -----------------------------------------------------------------------------
library(microbenchmark) 
ts <- microbenchmark(rwR1 <- randomwalk(sigma[1], x0, N) ,rwR2 <- randomwalk(sigma[2], x0, N) ,rwR3 <- randomwalk(sigma[3], x0, N) ,rwR4 <- randomwalk(sigma[4], x0, N) ,rwC1 <- rwme(sigma[1], x0, N) ,rwC2 <- rwme(sigma[2], x0, N) ,rwC3 <- rwme(sigma[3], x0, N)  ,rwC4 <- rwme(sigma[4], x0, N))
summary(ts)[,c(1,3,5,6)] 

