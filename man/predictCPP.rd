\name{predictCPP}
\alias{predictCPP}

\title{ VBMP predict functions parameters }
\description{
 Obtains estimates of class posterior probabilities from a fitted VBMP object
}
\usage{
predictCPP(obj, X.TEST=NULL)
}

\arguments{
  \item{obj}{ an object inheriting from class \code{VBMP.obj}, usually
    the result of a call to \code{vbmp} }
  \item{X.TEST}{ optionally, matrix in which to look for variables with which to predict. If omitted, the fitted predictors are used. }
}
\seealso{ See Also as \code{\link{vbmp}}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }