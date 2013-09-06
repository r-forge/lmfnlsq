\name{LMFnlsq}
\alias{LMFnlsq}
\title{Levenberg-Marquard-Fletcher method for curve fitting}
\description{\code{LMFnlsq} implements Fletcher's improvement of the Levenberg-Marquard method for solving non-linear least squares problems.}
\usage{
LMFnlsq(FUN,xc,LMFoptions=NA)
}
\arguments{
  \item{FUN}{\code{FUN} is a function to evaluate the residuals, from y (experimental data values): residuals = FUN(x) - y}
  \item{xc}{\code{xc} is a numeric vector with the initial guesses}
  \item{LMFoptions}{\code{LMFoptions} is a list with options}
}
\details{
LMFoptions:
Display  = 0         #   do not print iterations
Jacobian = 'finjac'  #   finite difference Jacobian approximation
MaxIter  = 0         #   maximum number of iterations allowed
ScaleD   = c()        #   automatic scaling by D = diag(diag(J'*J))
FunTol   = 1e-7      #   tolerace for final function value
XTol     = 1e-7      #   tolerance on difference of x-solutions
Printf   = 'printit' #   disply intermediate results
Trace    = 0         #   don't save  intermediate results
Lambda   = 0         #   start with Newton iteration
}
\value{
  \item{xf}{\code{xf} is a final solution approximation}
  \item{SS}{\code{SS} is the sum of squares of residuals}
  \item{Cnt}{\code{Cnt} is the number of iterations}
  \item{res}{\code{res} is the number of calls of FUN and Jacobian matrix}
  \item{XY}{\code{XY} are the points of intermediate results in iterations}  
}
\author{Jose Gama}
\examples{
\dontrun{demo(Rosenbrock)}
}
\keyword{programming}
