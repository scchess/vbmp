\name{plotDiagnostics}
\alias{plotDiagnostics}

\title{ VBMP covariance diagnostics plot }
\description{
plot the evolution of convergence diagnostics: lower-bound, 
predictive likelihood, out-of-bound test error and theta params (when inferred)  
}
\usage{
plotDiagnostics(obj)
}

\arguments{
  \item{obj}{ an object inheriting from class \code{VBMP.obj}, usually 
    the result of a call to \code{vbmp} }
}
\seealso{ See Also as \code{\link{vbmp}}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }