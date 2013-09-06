setwd('/home/mk/Desktop/z4/LMFnlsq article')
source('LMFnlsq.R')

#                       EXAMPLE 4:  Curve fit of decaying exponential
print('Exponential fit y(x) = c1 + c2*exp(c3*x)')

cR <- c(1,2,-1)
x <-  matrix(seq(0,3,by=0.1),ncol=1)         #   column vector of independent variable values

sx<-dim(x)[1]

y <- cR[1] + cR[2]*exp(cR[3]*x) + 0.1*matrix(rnorm(sx),ncol=1)  # dependent variable

                        #   Initial estimates:    
c1 <- y[length(y)]            #   c1 <- y(x->inf)
c2 <- y[1]-c1           #   c2 for x=0
c3 <- Re(mldivide(x[2:(length(x)-1)], log(abs(y[2:(length(y)-1)]-c1)/c2)))  #   evaluated c3 

res <- function(cR)  Re(cR[1] + cR[2]*exp(cR[3]*x) - y)   #   anonym. funct. for residuals

tR <- system.time(tmp <- LMFnlsq(res,c(c1,c2,c3),list(Display=-1)) )# without displ. lambda
C <- tmp$xf
ssq <- tmp$SS
cnt <- tmp$cnt

plot(x,y,xlim=c(0,3),ylim=c(0,3),main=expression(paste('Regression by f(x) = c',scriptscriptstyle(1),' + c',scriptscriptstyle(2),' ', plain(e)^{c[3]*x})))
grid()
par(new=T)
plot(x,res(C)+y,type='l',xlim=c(0,3),ylim=c(0,3),main='',col='red',xlab='',ylab='',axes=FALSE)

cat('time: ',(tR[1]),' sec\n')
