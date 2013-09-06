library('Matrix')
library('pracma')

LMFnlsq<-function(FUN,xc,LMFoptions=NA)
{
# R code by Jose Gama 2013
# Based on:
# LMFnlsq - Solution of nonlinear least squares
# http://www.mathworks.com/matlabcentral/fileexchange/17534-lmfnlsq-solution-of-nonlinear-least-squares
# M. Balda,
# Institute of Thermomechanics,
# Academy of Sciences of The Czech Republic,
# balda AT cdm DOT cas DOT cz
# miroslav AT balda DOT cz
#
# Reference:
# Fletcher, R., (1971): A Modified Marquardt Subroutine for Nonlinear Least
# Squares. Rpt. AERE-R 6799, Harwell

#######################################################################
# FUN,JAC,x,r,epsx
getAv<-function(FUN,JAC,x,r,epsx)
{
#        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  Calculate A, v, r

aRg <- formals(JAC)
naRg <-length(aRg)
if (naRg==1) J <- JAC(x) 
if (naRg==4) J <- JAC(FUN,r,x,epsx) 
#if (is.function(JAC)) J = JAC(x) else J = JAC(FUN,r,x,epsx)
##J = JAC(FUN,r,x,epsx)

A <- t(J) %*% J
v <- t(J) %*% r
list(A=A,v=v)
}
# --------------------------------------------------------------------

finjac<-function(FUN,r,x,epsx)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  numerical approximation to Jacobi matrix
{
rc <- matrix(r,ncol=1)
lx <- length(x)
J  <- matrix(0,nrow=length(r),ncol=lx)
for (k in 1:lx)
{
    dx <-0.25 * epsx[k]
    if (dx==0) dx<-2^(-52)
    xd <- x
    xd[k] <- xd[k]+dx
    rd <- FUN(xd)
#   ~~~~~~~~~~~~~~~~~~~
    J[,k]<-((matrix(rd,ncol=1)-rc) / dx)
}
J
}
# --------------------------------------------------------------------

 printit<-function(ipr,cnt,res,SS,x,dx,l,lc)
{
#        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  Printing of intermediate results
#  ipr <  0  do not print lambda columns
#      =  0  do not print at all
#      >  0  print every (ipr)th iteration
#  cnt = -1  print out the header
#         0  print out second row of results
#        >0  print out first row of results
if (ipr !=0)
   if (cnt<0)                 #   table header
	{
      cat('\n')
      nch <- 50+(ipr>0)*25
      paste(rep('*',nch),sep='',collapse='')
      cat('itr	nfJ	SUM(r^2)	x	dx\n')
      if (ipr>0) cat('	l	lc\n')
      paste(rep('*',nch),sep='',collapse='')
      cat('\n')
	}
   else if ((cnt %% ipr)==0) #   iteration output
	{
          f<-'#12.4e'
          if (ipr>0)
	  cat(cnt,'\t',res,'\t',SS,'\t', x[1],'\t',dx[1],'\t',l,'\t',lc,'\n')
          else cat(cnt,'\t',res,'\t',SS,'\t', x[1],'\t',dx[1],'\n')
          for (k in 2:length(x)) cat('		',x[k],'\t',dx[k],'\n')
      }
}
#######################################################################

# Options
LMFoptionsDefault <- list(Display = 1,Jacobian = finjac,MaxIter  = 0,ScaleD = c(),FunTol = 1e-7,XTol = 1e-7,Printf = printit, Trace = 0, Lambda   = 0)
#Display  = 0         #   no print of iterations
#Jacobian = 'finjac'  #   finite difference Jacobian approximation
#MaxIter  = 0         #   maximum number of iterations allowed
#ScaleD   = c()        #   automatic scaling by D = diag(diag(J'*J))
#FunTol   = 1e-7      #   tolerace for final function value
#XTol     = 1e-7      #   tolerance on difference of x-solutions
#Printf   = 'printit' #   disply intermediate results
#Trace    = 0         #   don't save  intermediate results
#Lambda   = 0         #   start with Newton iteration
printFlag=TRUE
if (any(is.na(LMFoptions))) LMFoptions <- LMFoptionsDefault else {
if (!is.numeric(LMFoptions$Display)) LMFoptions$Display <- LMFoptionsDefault$Display
if (!is.function(LMFoptions$Jacobian)) LMFoptions$Jacobian <- LMFoptionsDefault$Jacobian
if (!is.numeric(LMFoptions$MaxIter)) LMFoptions$MaxIter <- LMFoptionsDefault$MaxIter
if (!is.numeric(LMFoptions$ScaleD)) LMFoptions$ScaleD <- LMFoptionsDefault$ScaleD
if (!is.numeric(LMFoptions$FunTol)) LMFoptions$FunTol <- LMFoptionsDefault$FunTol
if (!is.numeric(LMFoptions$XTol)) LMFoptions$XTol <- LMFoptionsDefault$XTol
if (!is.function(LMFoptions$Printf)) LMFoptions$Printf <- LMFoptionsDefault$Printf else printFlag=FALSE
if (!is.numeric(LMFoptions$Trace)) LMFoptions$Trace <- LMFoptionsDefault$Trace
if (!is.numeric(LMFoptions$Lambda)) LMFoptions$Lambda <- LMFoptionsDefault$Lambda
}

#               INITIATION OF SOLUTION
#               **********************
x <- matrix(xc,ncol=1)
n <- dim(x)[1]
epsx <- matrix(LMFoptions$XTol,ncol=1)
le <- dim(epsx)[1]
if (le==1) epsx<-c(epsx) * matrix(1,nrow=n,ncol=1) else if (le != n) stop(paste('Dimensions of vector epsx ',le,'!=',n))

epsf  <- matrix(LMFoptions$FunTol,ncol=1)
ipr   <- LMFoptions$Display
JAC   <- LMFoptions$Jacobian
maxit <- LMFoptions$MaxIter    #   maximum permitted number of iterations
if (maxit==0) maxit<-100*n
#printf<- LMFoptions$Printf

r <- FUN(x)
tmp <- getAv(FUN,JAC,x,r,epsx)
A <- tmp$A
v <- tmp$v
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SS <- t(r) %*% r
res <- 1
cnt <- 0
trcXY <- LMFoptions$Trace      #   iteration tracing
if ( trcXY){
    XY <- matrix(0,nrow=n,ncol=maxit)
    XY[,1] <- x
    } else  XY <- c()

D <- LMFoptions$ScaleD #   CONSTANT SCALE CONTROL D
if (!is.null(D)) D <- matrix(D,ncol=1)
if (length(D)==0) {
    D=diag(A)              #   automatic scaling
} else {
    ld=length(D)
    if (ld==1)
        D=abs(D) * matrix(1,nrow=n,ncol=1) #   scalar of unique scaling
    else if (ld != n) stop(paste('wrong number of scales D, lD =',ld))
}
D <- matrix(D,ncol=1)
D[which(D<=0)]<-1
T <- matrix(sqrt(D),ncol=1)

Rlo<-0.25
Rhi<-0.75
l<-LMFoptions$Lambda
lc<-1
if (printFlag) LMFoptions$Printf(ipr,-1)       #   Table header
dx <- matrix(0,nrow=n,ncol=1)

#               SOLUTION
#               ********    MAIN ITERATION CYCLE
while (TRUE) #                   ********************
{
    if (printFlag) LMFoptions$Printf(ipr,cnt,res,SS,x,dx,l,lc)
    cnt <- cnt+1
    if (trcXY) { XY[,cnt+1] <- x}
    d = matrix(diag(A),ncol=1)
    s = matrix(0,nrow=n,ncol=1)
#                           INTERNAL CYCLE
    while (TRUE) #               ~~~~~~~~~~~~~~
{
        while (TRUE)
        {
            UA2 <- upper.tri(A,1)#triu(A,1)
            diag(UA2)<-FALSE
            UA <-A
            UA[-which(UA2)]<-0
            
            Utmp <- matrix(d+c(l) * D,ncol=1)
	    ##Utmp <- mdiag(c(Utmp))
	    Utmp <- diag(c(Utmp))
	    #print(Utmp)

            A <- t(UA)+UA+ Utmp	
#check to see if A is positive definite
	    U<-NA
            try(U <- chol(A), silent =TRUE); if (any(is.na(U))) p <- 1 else p<-0

            #~~~~~~~~~~~~~~~
            if (p==0) break
            l <- 2 * l
            if (l==0) l<-1
            }
        dx <- pracma::mldivide(U,pracma::mldivide(t(U),v)) # vector of x increments
        vw <- t(dx) %*% v
        fin <- -1
        if (vw<=0) break        #   The END

        for (n2 in 1:n)
{
            z <- d[n2] * dx[n2]

            if (n2>1) z<-A[n2,1:n2-1] %*% dx[1:n2-1]+z
            if (n2<n) z<-t(A[(n2+1):n,n2]) %*% dx[(n2+1):n]+z
            s[n2] <- 2 * v[n2]-z
        }
        dq <- t(s) %*% dx
        s  <- x-dx
        rd <- FUN(s)
#            ~~~~~~~~~~~~
        res <- res+1
        SSP <- t(rd) %*% rd
        dS  <- SS-SSP
        fin <- 1
        if (all((abs(dx)-epsx)<=0) | (res>=maxit) | (abs(dS)<=epsf)) break #   The END
        fin<-0
        if (dS>=Rlo*dq) break
        A <- U
        y <- 0.5
        z <- 2*vw-dS
        if (z>0) y<-vw/z
        if (y>.5) y<-0.5
        if (y<.1) y<-0.1
        if (l==0)
{
            y <- 2*y
            for (i in 1:n) A[i,i] <- 1/A[i,i]
            for (i in 2:n)
		{
                ii <- i-1
                for (j in 1:ii) A[j,i] <- -A[j,j:ii] %*% A[j:ii,i] %*% A[i,i]
		}
            for (i in 1:n)
                for (j in i:n) A[i,j] <- abs(A[i,j:n] %*% matrix(A[j,j:n],ncol=1) )#t(A[j,j:n])
            l <- 0
            tr <- matrix(diag(A),nrow=1) %*% D
            #str(tr)
            for (i in 1:n)
		{
                z <- t(A[1:i,i]) %*% T[1:i]+z
                #str(z)
                if (i<n)
			{
                    ii <- i+1
                    z  <- A[i,ii:n] %*% T[ii:n]+z
                	}
                z <- z * T[i]
                if (z>l) l<-z
            }
            if (tr<l) l<-tr
            l <- 1/l
            lc <- l
        }
        l <- l / y
        #str(l)
        if (dS>0){ dS <- -1e300; break}
    } #  while            INTERNAL CYCLE LOOP
#                           ~~~~~~~~~~~~~~~~~~~

    if (fin) break
    if (dS>Rhi*dq)
	{
        l<-l/2
        if (l<lc) l<-0
    	}
    SS<-SSP;  x<-s;  r<-rd
tmp <- getAv(FUN,JAC,x,r,epsx)
A <- tmp$A
v <- tmp$v
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
} # while                 MAIN ITERATION CYCLE LOOP
#                           *************************

if (fin>0)
    if (dS>0)
{
        SS <- SSP
        x  <- s
}
if (ipr !=0)
{
    if (printFlag) LMFoptions$Printf(sign(ipr),cnt,res,SS,x,dx,l,lc)
}
xf <- x

if (trcXY) XY[,cnt+2] <- x
XY<-XY[,1:(cnt+2)]

if (res>=maxit) cnt=-maxit
list(xf=xf, SS=SS, cnt=cnt, res=res, XY=XY) # xf, SS, cnt, res, XY
}

