\name{vbmp}
\alias{vbmp}

\title{Variational Bayesian Multinomial Probit Regression with Gaussian Process Priors.}
\description{
  Used to fit a Multinomial Probit Regression model, specified by giving the 
  matrix design \code{X}, the associated response variables \code{t.class}, kernel type and covariate 
  scaling parameters. Covariance paramters can be inferred from the data.
}
\usage{
vbmp(X, t.class, X.TEST, t.class.TEST, theta, control = list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{ Feature matrix for parameter 'estimation'}
  \item{t.class}{Target values, integer number used for class labels. }
  \item{X.TEST}{ Feature matrix to compute out-of-sample (test) prediction errors and likelihoods }
  \item{t.class.TEST}{ Target values for test data }
  \item{theta}{ The covariance function parameters (e.g. scaling coefficients for each dimension) }
  \item{control}{A list of control parameters. See Details}
}
\details{
 In this implementation a single covariance function is shared across all classes.
 Compute the predictive posteriors on the test set and the associated likelihood 
 and test errors at each iteration.

 The control argument is a list that can supply any of the following
 components:
\describe{
  \item{InfoLevel}{ 0 to suppress tracing ( > 0  to print different levels
          of monitoring information) }
  \item{sFILE.TRACE}{ File name where to redirect output (default NULL) }
  \item{bThetaEstimate}{ if covariance parameter estimation switched on. Defaults to FALSE (switched off) }
  \item{sKernelType}{ Kernel function used in training and predicting. 
   Currently implemented kernels are Gaussian (\code{"gauss"}), 
   Cauchy (\code{"cauchy"}), Laplace (\code{"laplace"}),
   Polynomial (\code{"poly"}), Homogeneous polynomial (\code{"hpoly"}),
   'Thin-plate' spline (\code{"tps"}), 'linear' spline (\code{"lsp"}) and
   Inner product(\code{"iprod"}).
   Defaults to \code{"gauss"}. }
  \item{maxIts}{ Maximum number of variational EM steps to take. 
     Defaults to 50.}
  \item{Thresh}{Convergence threshold on marginal likelihood lowerbound. 
     Defaults to 1e-4.}
  \item{method}{Integral computation method: "quadrature" (Gaussian quadrature)
     or "classic"(simple sampler).  Defaults to "quadrature".}     
  \item{nNodesQuad}{Number of nodes used for quadrature. Defaults to 49.}     
  \item{nSampsTG}{Number of samples used in obtaining mean of truncated 
     Gaussian. Defaults to 1000.}
  \item{nSampsIS}{Number of samples used in the importance sampler. 
     Defaults to 1000.}
  \item{nSmallNo}{Small number used to prevent numerical problems 
    (ill-conditioned covariance matrix). Defaults to 1e-10.}
  \item{parGammaTau,parGammaSigma}{The location and scale parameters of the 
     Gamma prior over covaraince params. Default to 1e-6.}
  \item{bMonitor}{TRUE to collect monitor convergence diagnostics at each
     iteration. Defaults to FALSE.}
  \item{bPlotFitting}{TRUE to plot test performance results at each iteration 
      during model estimation (if TRUE it forces bMonitor to TRUE). Defaults to FALSE.}
 }

}

\value{
\code{vbmp} returns an object of class "VBMP.obj". 
An object of class "VBMP.obj" is a list containing at least the following components: 
  \item{Kc}{Number of classes}
  \item{Ptest}{Matrix of multinomial class predictive posterior probabilities for the test data }
  \item{X}{Feature matrix}
  \item{invPHI}{Inverse of the Kernel matrix}
  \item{Y}{Matrix of auxiliary variables}
  \item{M}{Matrix of GP random variables}            
  \item{theta}{covariance kernel hyperparameters (estimates computed during 
       model fitting, if inferred}            
  \item{sKernelType }{Kernel function used in training and predicting}            
  \item{Test.Err}{Out-of-Sample Percent Prediction error estimates computed 
   during model fitting (0-1 error loss).}
  \item{PL}{Predictive Likelihood estimates computed during model fitting}
  \item{LOWER.BOUND}{Lower bound estimates computed during model fitting}            
  
}
%\seealso{ See Also as \code{\link{covParams,lowerBound,plotDiagnostics,
% predClass,predError, predLik}}

\references{ 
 Girolami M, Rogers S, \emph{Variational Bayesian Multinomial Probit Regression with Gaussian Process Priors}, Neural Computation 18, 1790-1817 (2006).
 Lama N, Girolami M \emph{vbmp: Variational Bayesian Multinomial Probit Regression for multi-class classification in R}, Bioinformatics 24(1):135-136 (2008). 
 \url{http://bioinformatics.oxfordjournals.org/cgi/content/short/btm535v1}
 }
\author{ N Lama \email{nicola.lama@unina2.it}, MA Girolami \email{girolami@dcs.gla.ac.uk} }
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{

## -----------------------------------------------------------------------------
## EXAMPLE 1 - Theta estimate with synthetic data
## -----------------------------------------------------------------------------
## Samples of 2-D data points drawn from three nonlinearly separable
## classes which take the form of two annular rings and one zero-centered
## Gaussian are used in this little illustrative example. 
genSample <- function(n, noiseVar=0) {
    ## class 1 and 2 (x ~ U(0,1))
    u <- 4. * matrix(runif(2*n), nrow=n, ncol=2) - 2.;
    i <- which(((u[, 1]^2 + u[, 2]^2) > .1) & ((u[, 1]^2 + u[, 2]^2) < .5) );
    j <- which(((u[, 1]^2 + u[, 2]^2) > .6) & ((u[, 1]^2 + u[, 2]^2) < 1) );
    X <- u[c(i, j),];
    t.class <- c(rep(1, length(i)),rep(2, length(j)));
    ## class 3 (x ~ N(0,1))
    x <- 0.1 * matrix(rnorm(2*length(i)), ncol=2, nrow=length(i) );
    k <- which((x[, 1]^2 + x[, 2]^2) < 0.1);
    X <- rbind(X, x[k, ]);
    t.class <- c(t.class, rep(3, length(k)));
    ## add random coloumns
    if (noiseVar>0) X <- cbind(X, matrix(rnorm(noiseVar*nrow(X)), ncol=noiseVar, nrow=nrow(X)));
    structure( list( t.class=t.class, X=X), class="MultiNoisyData");
}

set.seed(123); ## Init random number generator

## Generate training and test samples as an independent 
## test set to assess out-of-sample prediction error 
## and predictive likelihoods.
nNoisyInputs <- 0;       ## number of additional noisy input parameters
Ntest <- Ntrain <- 500;  ## sample sizes
dataXt.train <- genSample(Ntrain, nNoisyInputs);
dataXt.test  <- genSample(Ntest,  nNoisyInputs);

\dontrun{ 
 theta <- runif(ncol(dataXt.train$X));
 res <- vbmp( dataXt.train$X, dataXt.train$t.class,
        dataXt.test$X, dataXt.test$t.class, theta, 
         control=list(bThetaEstimate = T, 
         bPlotFitting=T, maxIts=50));
}

## set theta params (previously estimated) 
theta <- c(0.09488309, 0.16141604);   
## Fit the vbmp
res <- vbmp( dataXt.train$X, dataXt.train$t.class,
        dataXt.test$X, dataXt.test$t.class, theta, 
        control=list(maxIts=5));
## print out-of-sample error estimate
predError(res);

\dontrun{
## ----------------------------------------------------------
## EXAMPLE 2 - BRCA12 genomic data
## ----------------------------------------------------------
## Leave-one-out (LOO) cross-validation prediction error of the probabilistic 
## Gaussian process classifier used in Zsofia Kote-Jarai et al. 
## Clin Cancer Res 2006;12(13);3896-3901

  if(any(installed.packages()[,1]=="Biobase")) {
    library("Biobase");
    data("BRCA12");
    brca.y <- BRCA12$Target.class;
    brca.x <- t(exprs(BRCA12));
  } else {
    print("Deprecated.....");
    load(url("http://www.dcs.gla.ac.uk/people/personal/girolami/pubs_2005/VBGP/BRCA12.RData"));
    brca.y <- as.numeric(BRCA12$y);
    brca.x <- as.matrix(BRCA12[,-1]);
  }
  
  sKernelType <- "iprod";  ## Covariance function type
  Thresh <- 1e-8;  ## Iteration threshold
  InfoLevel <- 1;
  theta <- rep(1.0, ncol(brca.x));
  ITER.THETA <- 24;
  n     <- nrow(brca.x) ;
  Kfold <- n; # number of folds , if equal to n then LOO
  samps <- sample(rep(1:Kfold, length=n), n, replace=FALSE); 
  res   <- rep(NA, n);
  print(paste("LOO crossvalidation started...... (",n,"steps)"));
  for (x in 1:Kfold) {
      cat(paste(x,", ",sep="")); flush.console();
      resX <- vbmp( brca.x[samps!=x,], brca.y[samps!=x], 
                    brca.x[samps==x,], brca.y[samps==x], 
                    theta,  control=list(bThetaEstimate=F, 
                    bPlotFitting=F, maxIts=ITER.THETA, 
                    sKernelType=sKernelType, Thresh=Thresh));    
      res[samps==x] <- predClass(resX); 
  }
  print("(end)");
  print(paste("Crossvalidated error rate", round(sum(res!=brca.y)/n,2)));
}   


}
%
% res <- vbmp( dataXt.train$X, dataXt.train$t.class, dataXt.test$X, 
%      dataXt.test$t.class, theta, control=list(InfoLevel=3, maxIts=3));
% res <- vbmp( dataXt.train$X, dataXt.train$t.class, dataXt.test$X, 
%      dataXt.test$t.class, theta, control=list(bPlotFitting=T,InfoLevel=3, 
%      bThetaEstimate=T, maxIts=100));
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\seealso{ See Also as \code{\link{predictCPP}}, \code{\link{covParams}}, 
 \code{\link{lowerBound}}, \code{\link{predError}}, 
 \code{\link{predLik}}, \code{\link{predClass}}}

\keyword{ models }

