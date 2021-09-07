\name{clse}
\alias{clse}
\title{
Constrained Least-Squares Estimates
}
\description{
This function computes the constrained least-squares estimates when a subset of or all of the regression coefficients are constrained to be non-negative, as described in Fraser and Massam (1989).
}
\usage{
clse(dat.obj)
}

\arguments{
\item{dat.obj}{
A list with the following format, \code{list(y, mat, lam)}. Here \code{y} is the response vector, \code{mat} is the design matrix for the regression, and \code{lam} is a vector with the length that matches the number of columns in \code{mat}. The values of \code{lam} is either 0 or 1, with 0 means unconstrained and 1 means the corresponding regression coefficient is constrained to be non-negative.  
}
}

\value{
The returned value is a list with format, \code{list(dat.obj, beta.vec, yhat)}. Here \code{dat.obj} is the input of the function, \code{beta.vec} gives the estimated regression coefficient, and \code{yhat} is the vector for the fitted response values.  
}
\references{
Fraser, D. A. S. and H. Massam (1989). A mixed primal-dual bases algorithm for regression
under inequality constraints. Application to concave regression. \emph{Scandinavian Journal
of Statistics} 16, 65-74.
}


\examples{
#generate a dataset for illustration.
x=seq(1,10,,100)
y=x^2+rnorm(length(x))
#generate spline bases.
tmp=MIC.splines.basis.fast(x=x, df = 10, knots = NULL, boundary.knots=NULL,
type="Is",degree = 3,delta=0.001,eq.alloc=FALSE)
#plot the spline bases.
plot(tmp)
#generate the data object for the clse function.
dat.obj=list(y=y, mat=cbind(1, tmp$mat), lam=c(0, rep(1, ncol(tmp$mat))))
#fit clse.
fit=clse(dat.obj=dat.obj)
#visualize fitted results.
plot(x, y, pch=16)
lines(x, fit$yhat, lwd=3, col=2)

}


