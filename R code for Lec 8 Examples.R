#######################################################
    ###       R code for Lec 8 Examples                 ###
    #######################################################
    
### Example 1 (Gibbs sampler: Bivariate distribution)

    #initialize constants and parameters
    N <- 5000               #length of chain
    burn <- 1000            #burn-in length
    X <- matrix(0, N, 2)    #the chain, a bivariate sample

    rho <- -.75             #correlation
    mu1 <- 0
    mu2 <- 2
    sigma1 <- 1
    sigma2 <- .5
    s1 <- sqrt(1-rho^2)*sigma1
    s2 <- sqrt(1-rho^2)*sigma2

    ###### generate the chain #####

    X[1, ] <- c(mu1, mu2)            #initialize

    for (i in 2:N) {
        x2 <- X[i-1, 2]
        m1 <- mu1 + rho * (x2 - mu2) * sigma1/sigma2
        X[i, 1] <- rnorm(1, m1, s1)
        x1 <- X[i, 1]
        m2 <- mu2 + rho * (x1 - mu1) * sigma2/sigma1
        X[i, 2] <- rnorm(1, m2, s2)
    }

    b <- burn + 1
    x <- X[b:N, ]

    # compare sample statistics to parameters
    colMeans(x)
    cov(x)
    cor(x)

    plot(x, main="", cex=.5, xlab=bquote(X[1]),
         ylab=bquote(X[2]), ylim=range(x[,2]))
         
### Example 2  A Bayesian inference example: Body temperature data

  bodytemp<-read.table("bodytemp.txt",header=T)
   y<-bodytemp$temp
   bary<-mean(y); n<-length(y)
   Iterations<-5000
   mu0<-0; s0<-100; a0<-0.001; b0<-0.001
   theta <- matrix(nrow=Iterations, ncol=2)
   cur.mu<-0; cur.tau<-2; cur.s<-sqrt(1/cur.tau)
   for (t in 1:Iterations){
       w<- s0^2/( cur.s^2/n+ s0^2 )
       m <- w*bary + (1-w)*mu0
       s <- sqrt( w/n ) * cur.s
       cur.mu <- rnorm( 1, m, s )
       a <- a0 + 0.5*n
       b <- b0 + 0.5 * sum( (y-cur.mu)^2 )
       cur.tau <- rgamma( 1, a, b )
       cur.s <- sqrt(1/cur.tau)
       theta[t,]<-c( cur.mu, cur.s)
   }


mcmc.output<-theta
apply(mcmc.output[-(1:1000),],2,mean)
#compare to true value: 98.25, 0.542
apply(mcmc.output[-(1:1000),],2,sd)
#compare to true value: 0.06456, 0.06826   
         
#####diagnostics plots
erg.mean<-function(x){ # compute ergodic mean 
        n<-length(x)
        result<-cumsum(x)/cumsum(rep(1,n))
  }
  
par( mfrow=c(3,2), xaxs='r', yaxs='r', bty='l' , cex=0.8)
iter<-1500
burnin<-500
index<-1:iter
index2<-(burnin+1):iter

plot(index, theta[index,1], type='l', ylab='Values of mu', xlab='Iterations', main='(a) Trace Plot of mu')
plot(index, theta[index,2], type='l', ylab='Values of sigma', xlab='Iterations', main='(b) Trace Plot of sigma')


ergtheta0<-erg.mean( theta[index,1] )
ergtheta02<-erg.mean( theta[index2,1] )
ylims0<-range( c(ergtheta0,ergtheta02) )

ergtheta1<-erg.mean( theta[index,2] )
ergtheta12<-erg.mean( theta[index2,2] )
ylims1<-range( c(ergtheta1,ergtheta12) )

step<-10
index3<-seq(1,iter,step)
index4<-seq(burnin+1,iter,step)

plot(index3 , ergtheta0[index3], type='l', ylab='Values of mu', xlab='Iterations', main='(c) Ergodic Mean Plot of mu', ylim=ylims0)
lines(index4, ergtheta02[index4-burnin], col=2, lty=2)


plot(index3, ergtheta1[index3], type='l', ylab='Values of sigma', xlab='Iterations', main='(d) Ergodic Mean Plot of sigma', ylim=ylims1)
lines(index4, ergtheta12[index4-burnin], col=2, lty=2)



acf(theta[index2,1], main='Autocorrelations Plot for mu')
acf(theta[index2,2], main='Autocorrelations Plot for sigma')       
   
### Example 3 Slice Gibbs Sampler 
wais<-read.table("wais.txt",header=T)
y<-wais$senility; x<-wais$wais; n<-length(y)
positive<- y==1 
Iterations<-55000
mu.beta<-c(0,0); s.beta<-c(100,100) 
beta <- matrix(nrow=Iterations, ncol=2) 
acc.prob <- 0 
current.beta<-c(0,0); u<-numeric(n)  
for (t in 1:Iterations){  
  eta<-current.beta[1]+current.beta[2]*x
    U<-exp(y*eta)/(1+exp(eta))
    u<-runif( n, rep(0,n), U) 
    logitu<-log( u/(1-u) )
    logitu1<-  logitu[positive]
    logitu2<- -logitu[!positive]
    
    l0<- max( logitu1 - current.beta[2]*x[positive] )
    u0<- min( logitu2 - current.beta[2]*x[!positive] )
    unif.random<-runif(1,0,1)
    fa<- pnorm(l0, mu.beta[1], s.beta[1])
    fb<- pnorm(u0, mu.beta[1], s.beta[1])
    current.beta[1] <- qnorm( fa + unif.random*(fb-fa), mu.beta[1], s.beta[1])
    
    l1<- max( (logitu1 - current.beta[1])/x[positive] )
    u1<- min( (logitu2 - current.beta[1])/x[!positive] )
    unif.random<-runif(1,0,1)
    fa<- pnorm(l1, mu.beta[2], s.beta[2])
    fb<- pnorm(u1, mu.beta[2], s.beta[2])
    current.beta[2] <- qnorm( fa + unif.random*(fb-fa), mu.beta[2], s.beta[2])
    beta[t,]<-current.beta  
}

apply(beta[-(1:15000),],2,mean)
apply(beta[-(1:15000),],2,sd)

## convergence diagnositcs plots

par( mfrow=c(3,2), xaxs='r', yaxs='r', bty='l' , cex=0.8)
iter<-55000
burnin<-15000
index<-seq(1,iter,50)
index2<-(burnin+1):iter

plot(index, beta[index,1], type='l', ylab='Values of beta0', xlab='Iterations', main='(a) Trace Plot of beta0')
plot(index, beta[index,2], type='l', ylab='Values of beta1', xlab='Iterations', main='(b) Trace Plot of beta1')

iter<-55000
burnin<-15000
index<-seq(1,iter,1)
index2<-(burnin+1):iter

ergbeta0<-erg.mean( beta[index,1] )
ergbeta02<-erg.mean( beta[index2,1] )
ylims0<-range( c(ergbeta0,ergbeta02) )

ergbeta1<-erg.mean( beta[index,2] )
ergbeta12<-erg.mean( beta[index2,2] )
ylims1<-range( c(ergbeta1,ergbeta12) )

step<-50
index3<-seq(1,iter,step) 
index4<-seq(burnin+1,iter,step)

plot(index3 , ergbeta0[index3], type='l', ylab='Values of beta0', xlab='Iterations', main='(c) Ergodic Mean Plot of beta0', ylim=ylims0)
lines(index4, ergbeta02[index4-burnin], col=2, lty=2)


plot(index3, ergbeta1[index3], type='l', ylab='Values of beta1', xlab='Iterations', main='(d) Ergodic Mean Plot of beta1', ylim=ylims1)
lines(index4, ergbeta12[index4-burnin], col=2, lty=2)




lag.to.print<-900
acf1<-acf(beta[index2,1], main='Autocorrelations Plot for beta0', lag.max=lag.to.print, plot=FALSE)
acf2<-acf(beta[index2,2], main='Autocorrelations Plot for beta1', lag.max=lag.to.print, plot=FALSE)

acf.index<-seq(1,lag.to.print,20)

plot( acf1[acf.index], main='Auto-correlations for beta0' )
plot( acf2[acf.index], main='Auto-correlations for beta1' )
     
         
### Example  4 (Gelman-Rubin method of monitoring convergence)

    Gelman.Rubin <- function(psi) {
        # psi[i,j] is the statistic psi(X[i,1:j])
        # for chain in i-th row of X
        psi <- as.matrix(psi)
        n <- ncol(psi)
        k <- nrow(psi)

        psi.means <- rowMeans(psi)     #row means
        B <- n * var(psi.means)        #between variance est.
        psi.w <- apply(psi, 1, "var")  #within variances
        W <- mean(psi.w)               #within est.
        v.hat <- W*(n-1)/n + B/n+(B/(n*k))     #upper variance est.
        r.hat <- v.hat / W             #G-R statistic
        return(r.hat)
        }

    normal.chain <- function(sigma, N, X1) {
        #generates a Metropolis chain for Normal(0,1)
        #with Normal(X[t], sigma) proposal distribution
        #and starting value X1
        x <- rep(0, N)
        x[1] <- X1
        u <- runif(N)

        for (i in 2:N) {
            xt <- x[i-1]
            y <- rnorm(1, xt, sigma)     #candidate point
            r1 <- dnorm(y, 0, 1) * dnorm(xt, y, sigma)
            r2 <- dnorm(xt, 0, 1) * dnorm(y, xt, sigma)
            r <- r1 / r2
            if (u[i] <= r) x[i] <- y else
                 x[i] <- xt
            }
        return(x)
        }

    sigma <- 0.2     #parameter of proposal distribution
    k <- 4          #number of chains to generate
    n <- 15000      #length of chains
    b <- 1000       #burn-in length

    #choose overdispersed initial values
    x0 <- c(-10, -5, 5, 10)

    #generate the chains
    X <- matrix(0, nrow=k, ncol=n)
    for (i in 1:k)
        X[i, ] <- normal.chain(sigma, n, x0[i])
    #trace plots
     plot(1:n,X[1,],type="l")
     lines(1:n,X[2,],type="l",col=2)
     lines(1:n,X[3,],type="l",col=3)
     lines(1:n,X[4,],type="l",col=4)


    #compute diagnostic statistics
    psi <- t(apply(X, 1, cumsum))
    for (i in 1:nrow(psi))
        psi[i,] <- psi[i,] / (1:ncol(psi))
    print(Gelman.Rubin(psi))

    #plot psi for the four chains
    par(mfrow=c(2,2))
    for (i in 1:k)
        plot(psi[i, (b+1):n], type="l",
            xlab=i, ylab=bquote(psi))
    par(mfrow=c(1,1)) #restore default

    #plot the sequence of R-hat statistics
    rhat <- rep(0, n)
    for (j in (b+1):n)
        rhat[j] <- Gelman.Rubin(psi[,1:j])
    plot(rhat[(b+1):n], type="l", xlab="", ylab="R")
    abline(h=1.1, lty=2)


### Example  5 (Coal mining disasters)
    # Gibbs sampler for the coal mining change point

    library(boot)     #for coal data
    data(coal)
    year <- floor(coal)
    y <- table(year)
    plot(y)  #a time plot

    y <- floor(coal[[1]])
    y <- tabulate(y)
    y <- y[1851:length(y)]
    
    #y<-dput(as.numeric(y))
    n=112
   y<-c(4, 5, 4, 1, 0, 4, 3, 4, 0, 6, 3, 3, 4, 0, 2, 6, 3, 3, 5, 4, 
      5, 3, 1, 4, 4, 1, 5, 5, 3, 4, 2, 5, 2, 2, 3, 4, 2, 1, 3, 2, 2, 
      1, 1, 1, 1, 3, 0, 0, 1, 0, 1, 1, 0, 0, 3, 1, 0, 3, 2, 2, 0, 1, 
      1, 1, 0, 1, 0, 1, 0, 0, 0, 2, 1, 0, 0, 0, 1, 1, 0, 2, 3, 3, 1, 
      1, 2, 1, 1, 1, 1, 2, 3, 3, 0, 0, 0, 1, 4, 0, 0, 0, 1, 0, 0, 0, 
      0, 0, 1, 0, 0, 1, 0, 1) 
  data=list("n","y")
  parameters <- c("k","b")
  inits = function() {list(b=c(0,0),k=50)}     
  
  # model can be done in R
  # coal<-function(){for( i in 1 : n ) {
  #    y[i] ~ dpois(mu[i])
  #   log(mu[i]) <- b[1] + step(i - k) * b[2]
  #   }
  #  for (j in 1:2) {
  #   b[j] ~ dnorm( 0.0,1.0E-6)
  #    }
  #   k ~ dunif(1,n)
  #    }
  # write.model(coal,"coal.bug")
   
  
  coal.sim <- bugs(data, inits, parameters,"coal.bug", n.chains=3, n.iter=10000,n.burnin=2000,n.thin=1,codaPkg=T,DIC=F,bugs.directory="D:/WinBUGS14")
  attach.bugs(coal.sim)
  print(coal.sim)
  par(mfrow=c(2,1))
  plot(density(b[,1]),xlab="beta1")
  plot(density(b[,2]),xlab="beta2")
 
 # use CODA package to do further analysis
   
  coal.out<-read.bugs(coal.sim)
  library(coda)
  autocorr.plot(coal.out)
  gelman.plot(coal.out)

