source("GC_Sim.R")
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



