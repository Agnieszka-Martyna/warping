library(stats4)
library(methods)
# creating the class signal.
# 
rm(list = ls())
# Create a class with three fields.
# peak: A vector that contains the position of peaks in the range of spectra.
# variance: A vector that contains the variances for different peaks.
# intensity: A vector that contains the intensity of the peaks.
# Window: A vector that contains the number of retention time to be affected
#         around the main peak.
# lowerlimit: A variable that set the minimum value for retention time.
# upperlimit: A variable that set the maximum value for retention time.

gc <- setRefClass("gc",fields = list(peak = "numeric",
                                   variance = "numeric",
                                   intensity = "numeric",
                                   window = "numeric",
                                   lowerlimit = "numeric",
                                   upperlimit = "numeric"),
                     methods = list(
                       print = function()
                       {
                         if(!any(length(peak) == length(variance) | 
                                 length(variance) == length(intensity) |
                                 length(peak == length(intensity)))) 
                           stop("vector peak, variance, and intensity must be same length")
                          
                         # number of peaks
                         np <- length(peak)
                         # create the vector for all the range of retention time
                          x <- rep(0,length(seq(lowerlimit,upperlimit,1)))
                          
                          # for each peak identify the window to be changed
                          for(j in 1:np)
                            {
                            if(peak[j]-window[j]>0)
                            {
                              ini<- peak[j]-window[j]
                              fin<- peak[j]+window[j]
                              # obtain position to be changed
                              positionToChange <- seq(ini,fin,1)
                              
                            x[ini:fin]<-x[ini:fin] + sapply(positionToChange,function(pos)
                              {
                                # generate a value from normal dist. with mean equal to
                                # intensity
                                abs(rnorm(1,intensity[j],sqrt(variance[j])))/exp(abs(peak[j]-pos)/2)
                                # exponential decay when the retention time differ from 
                                # the position of the peak
                                
                              })                  
                              
                            }
                          }
                        print = x 
                       },
                       range = function()
                       {
                         range = seq(lowerlimit,upperlimit,1)
                       }
                     ))



set.seed(1)
v1 <- c(50,110,140,175,190,230,250)
u <- round(rnorm(7,5,2),0)
u
v2 <- v1 + u
v2

gc1<-gc(peak =v1,
               variance =c(2,2,1,1,1,1,1),
               intensity = c(10,35,55,35,180,700,400),
               window = c(3,3,3,5,5,10,10),
               lowerlimit = 0,
               upperlimit = 500)

# signal shifted by adding a random normal vector to the peak vectors
gc2<-gc(peak =v2,
                variance =c(2,2,1,1,1,1,1),
                intensity = c(10,35,55,35,180,700,400),
                window = c(3,3,3,5,5,10,10),
                lowerlimit = 0,
                upperlimit = 500)

B1 <- gc1$print()
# signal shifted
T1 <- gc2$print()

plot(gc1$range(), B1, type = "l", col = "blue")
lines(gc1$range(),T1, col = "green")


Pprime <- function(P,xi,xi1)
{
  
  
}

# benefit function 
f <- function(P,T1,xi,xi1)
{
  
  # Use the warping points to find the unwarped points and obtain pj's
  # calculate the Pprime using pj's in order to compare with the target
  # Calculate correlation between Pprime and T in a segment.
  ini<- which(T$x == xi)
  fin <-which(T$x == xi1)
  
  # Transform P into Pw(Pprime)
  
  # if Pw and T1 are the same length
  a <- T1$f[ini:fin]
  b <- Pw$f[ini:fin]
  output <- cov(a,b)
  
}

align <- function(P,Target,m,t)
{
  # F1: matrix containing the cumulated benefit function
  
  # Pre-aligming length of chromatogram
  Lp<- max(P$x)-min(P$x)
  
  # Post-aligning length of chromatogram and length of 
  # target chromatogram
  Lt <- max(Target$x)-min(Target$x)
  
  # Calculate number of sections for P
  N <- round(Lp/m)
  
  # calculate difference in mean section length between P and T
  d <- round(Lt/N) - m
  
  F1 <- matrix(0.0, nrow = (N+1), ncol = (Lt+1))
  U <- matrix(0.0, nrow = (N+1), ncol = (Lt+1))
  
  for (i in 0:N)
  {
    for(x in 0:Lt)
    {
      F1[i,x] <- -Inf  
    }
  }
  
  F1[N,0] = 0
  
  for (i in (N-1): 0)
  {
    xstart <- max(i*(m+d-t),Lt-(N-i)*(m+d+t))
    xend <- min(i*(m+d + t),Lt-(N-i)*(m+d-t))
    for (x in xstart:xend)
    {
      for (u in (d-t):(d+t))
      {
        fsum <- F1[i+1,x+m+u] + f(x,x+m+u)
        if (fsum > F1[i,x]) then
        F1[i,x] <- fsum
        U[i,x] <- u
      }
    }
  }
  
  x(0) = 0
  for(i in 0:N-1)
  {
    u[i] <- U[i,x[i]]
    x[i+1] <- x[i]+m+u(i)
  }
  return(x = x, u=u,F1=F1,U=U)
}
