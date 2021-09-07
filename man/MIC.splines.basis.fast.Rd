\name{MIC.splines.basis.fast}
\alias{MIC.splines.basis.fast}

\title{Generating MIC Spline Bases}
\description{
This function provides C implementation for the computing of monotonic spline bases, including M-splines, I-splines, and C-splines, denoted by MIC splines. The definitions of the spline bases are described in Meyer (2008).
}
\usage{
MIC.splines.basis.fast(x, df = NULL, knots = NULL, boundary.knots = NULL, 
type = "Ms", degree = 3, delta = 0.01, eq.alloc = FALSE)
}

\arguments{
  \item{x}{A numeric vector for the data to generate spline bases for.
}
  \item{df}{The degree of freedom, which equals to the number of interior knots plus the spline degree.
}
  \item{knots}{A vector for the interior knots.
  }
  \item{boundary.knots}{ The values for the left and right boundary points.
}
  \item{type}{The type of splines to be computed. \code{"Ms"} stands for M-splines, \code{"Is"} stands for I-splines, \code{"IsN"} stands for I-splines without normalization, and \code{"Cs"} stands for C-splines. 
}
  \item{degree}{The degree for the M-splines. I-splines are based on the integration of the M-splines, and C-splines are based on the integration of the I-splines. 
}
  \item{delta}{A numeric value that is used to set the bin width for numerical integration. Usually it is set to a small number.}
  \item{eq.alloc}{A logic variable, which is true if using equal spacing for the interior knots, and is false if using equal quantiles for the interior knots.  
}
}

\value{
A list with format, \code{list(mat, x, ...)}. Here \code{mat} is the matrix for the spline bases, \code{x} is the vector for the data, and the rest of the items are carrying the information from the arguments. 
}
\references{
Meyer, M. C. (2008). Inference using shape-restricted regression splines. \emph{The Annals of
Applied Statistics} 2, 1013-1033.
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
}

