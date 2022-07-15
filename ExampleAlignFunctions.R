pathGCfunctions <- "E:/University of Glasgow/AgnieskaPRoject/GC_MS"
setwd(pathGCfunctions)
source(paste0(pathGCfunctions,"/AlignFunctions.R"))


# load signals

# artificial signals 1 channel
T1 = c(0,1,1,2,5,6,5,5,4,4,3,2,1,1,0,0)
X1 = c(0,0,1,1,1,2,3,5,5,6,6,5,3,1,1,0)
X2  = c(0,0,1,1,1,2,3,5,5,6,6,5,2,1,1,0)

# artificial signals 2 channel GC-MS
MT <- cbind(T1,T1)
MX <- cbind(X1,X2)

Seg = 5
Slack = 1
Options = c(0,1,1,0,0)


# First we align X1 to T1 and X2 to T1 separately
# using the function

WX1 <- cow(T1,X1,Seg,Slack,Options)
WX1$Warping  # nodes position before and after warping
WX1$XWarped  # Warped signal for 


WX2 <- cow(T1,X2,Seg,Slack,Options)
WX2$Warping # nodes position before and after warping
WX2$XWarped # Warped signal

# GC-MS aligning
# if we have a GCMS signal with two channels 
# represented by 16 rows x 2 columns
MX
# and a Target GCMS signal
MT

# We can align MX to MT that would obtainn as a result
# the alignment of X1 and X2
cbind(as.vector(WX1$XWarped),as.vector(WX2$XWarped) )

WMX <- alignGCMS(MT,MX, Seg, Slack)
WMX
WMX$W

sum(round(WMX$W - as.matrix(cbind(as.vector(WX1$XWarped),as.vector(WX2$XWarped))),2))

# WMX$W contain the aligned matrix we are looking for