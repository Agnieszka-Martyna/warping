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


Compare2Gc <- function (P,Target,m,t)
{
  Lp <- max(P)-min(P) 
  Lt <- max(Target)-min(Target)
  N <- round(Lp/m)
  x <- 0:length(P)
  delta <- round(Lt/N - m)
  Interval <- seq(delta-t,delta+t,1)
  nwarping <- length(Interval)
  for (ui in 1:nwarping)
  {
    xadj <- sapply(1:(N-1),function(i){ x[i] + m + Interval[ui] })
    cov()
                                                                                                                                                                                                                                                
  }
  
  
}


