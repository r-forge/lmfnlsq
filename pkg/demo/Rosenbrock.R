#   LMFnlsqtest.m          Constrained Rosenbrock's valey
##########################################################
#   The script solves a testing problem of the Rosenbrock's function by
#   minimization of of a sum of squares of residuals and a curve fitting.
#   Example:
#   A user may run the script multiply changing only few parameters:
#       iprint      as a step in displaying intermediate results,
#       ScaleD      diagonal scale matrix, and
#       Trace       a control variable for storing intermediate data.

# Miroslav Balda
# miroslav AT balda DOT cz
#   2008-08-18  v 1.1   Modified for analytical gradient
#   2009-01-06  v 1.2   updated for modified function LMFnlsq
setwd('/home/mk/Desktop/z4/LMFnlsq article')
source('LMFnlsq.R')
cat('Rosenbrock without constrains\n')
par( mfrow = c( 1, 3 ) )
ipR<-5
sD<-c()
xy<-1
x0 <- c(-1.2, 1) #   Usual starting point for Rosenbrock valey
for (k in 1:3){             #   Cycle for analytical | finite differences gradient
    if (k==1){             #   EXAMPLE 1:  Rosenbrock without constrains
        r <- 0
        r2 <- 0
        gr <- 'AG'    #   Analytical gradient
        ros <-  function(x) matrix(c(10*(x[2]-x[1]^2), 1-x[1]),2,1,byrow=F)
        jac <-  function(x) matrix(c(-20*x[1], 10, -1, 0),2,2,byrow=T)
        cat('Analytical gradient\n')
        # With analytical Jacobian matrix
tR <- system.time( tmp <- LMFnlsq(ros,x0,list(Display = ipR, Jacobian = jac, ScaleD = sD, Trace = xy )) )

xf <- tmp$xf
ssq <- tmp$ssq
cnt <- tmp$cnt
loops <- tmp$loops
XY <- tmp$XY
  }  else if (k==2) {               #   EXAMPLE 2:  Rosenbrock with constraints
        cat('Rosenbrock with constrains\n')
        gr <- 'FDG'   #   Finite difference approx. of gradient
        r <- 0.5
        w <- 1000
        d <- function(x) t(x) %*% x-r^2 #    delta of squares of position and radius
        ros <-  function(x) matrix(c(10*(x[2]-x[1]^2), 1-x[1], (r>0)*(d(x)>0)*d(x)*w),3,1,byrow=F)

        cat('Gradient from finite differences\n')
        # With finite difference Jacobian matrix
        tR <- system.time( tmp <- LMFnlsq(ros,x0,list(Display = ipR, ScaleD = sD, Trace = xy )) )
xf <- tmp$xf
ssq <- tmp$ssq
cnt <- tmp$cnt
loops <- tmp$loops
XY <- tmp$XY
    } else if (k==3) {               #   EXAMPLE 3:  Rosenbrock with upper/lower boundaries  +-0.5
        cat('Rosenbrock with upper/lower boundaries  +-0.5\n')
        gr <- 'FDG'   #   Finite difference approx. of gradient
        r2 <- 0.5
        w2 <- 1000
        ros <-  function(x) matrix(c(10*(x[2]-x[1]^2), 1-x[1], (x[2]>r2)*(x[2]-r2)*w2, (x[2]< -r2)*(x[2]+r2)*w2, (x[1]>r2)*(x[1]-r2)*w2, (x[1]< -r2)*(x[1]+r2)*w2 ),6,1,byrow=F)
        cat('Gradient from finite differences\n')
        # With finite difference Jacobian matrix
        tR <- system.time( tmp <- LMFnlsq(ros,x0,list(Display = ipR, ScaleD = sD, Trace = xy )) )
xf <- tmp$xf
ssq <- tmp$ssq
cnt <- tmp$cnt
loops <- tmp$loops
XY <- tmp$XY
    }# if

    R <- sqrt(t(xf) %*% xf)
    cat('\n  Distance from the origin R =',R,',   R^2 = ',R^2,'\n')
cat('time: ',(proc.time()-tR)[1],' sec\n')

    if (xy){                               #   Saved sequence [x(1), x(2)] 

if (r>0) theTitle <- 'Constrained Rosenbrock valley' else theTitle <- 'Rosenbrock valley'
theTitle <- paste(theTitle, gr)

x<-seq(-2,2,by=.1)
y<-seq(-2,2,by=.1)
Y<-matrix(rep(x,length(x)),nrow=length(x),byrow=T)
X<-matrix(rep(y,length(y)),nrow=length(y),byrow=F)
Z <- (100*(Y-X^2)^2 - (1-X)^2)   #   Rosenbrock's function
image(x,y,Z,xlim=c(-2,2),ylim=c(-2,2), xlab='',ylab='',axes=FALSE)
par(new=T)
contour(x,y,Z,drawlabels=F,nlevels=30,xlim=c(-2,2),ylim=c(-2,2), xlab='X1',ylab='X2',main=theTitle)
par(new=T)
plot(x0[1],x0[2],type='p', lwd=2,xlim=c(-2,2),ylim=c(-2,2),xlab='',ylab='',axes=FALSE)#   starting point
par(new=T)
plot(xf[1],xf[2],type='p', lwd=2,xlim=c(-2,2),ylim=c(-2,2),xlab='',ylab='',axes=FALSE)#   terminal point
par(new=T)
plot(c(x0[1],XY[1,]),c(x0[2],XY[2,]),type='l',xlim=c(-2,2),ylim=c(-2,2), lwd=3,xlab='',ylab='',axes=FALSE)#   iteration path

if (k==2){ # circle = feasible domain
fi<-seq(0,2*pi,by=pi/18)
par(new=T)
plot(cos(fi)*r,sin(fi)*r,xlim=c(-2,2),ylim=c(-2,2),type='l',xlab='',ylab='',axes=FALSE,col='white')   #   circle 
}
if (k==3){ # square = feasible domain
fi<-seq(0,2*pi,by=pi/18)
par(new=T)
segments(-r2,-r2,r2,-r2,col='white')#   square 
segments(r2,r2,-r2,r2,col='white')
segments(-r2,-r2,-r2,r2,col='white')
segments(r2,r2,r2,-r2,col='white')
}
    }
    
} 





