`predictCPP` <-
function(obj, X.TEST=NULL) {
   ## Computes the predictive posteriors as defined in Section 4.5 of the paper.
   ## first two equations at page 1798
   con <- list(InfoLevel=0, sFILE.TRACE=NULL, bThetaEstimate=FALSE,
        sKernelType="gauss", maxIts=50, Thresh=1e-4, tmpSave=NULL,
        nSampsTG=1000, nSampsIS=1000, nSmallNo=1e-10, parGammaSigma=1e-6,
        parGammaTau=1e-6, bPlotFitting=FALSE);
   con[names(obj$con)] <- obj$con;
   X  <- obj$X;
   if (length(X.TEST)==0) X.TEST <- X;
   Y <- obj$Y;
   M <- obj$M;
   Kc <- ncol(Y);        ## Identify the number of classes
   N  <- nrow(X); Kd <- ncol(X);   ## Get number of samples and dimension of data
   theta <- covParams(obj);
   Ntest <- nrow(X.TEST);     ## Number of test points
   invPHI <- obj$invPHI;
   PHItest <- computeKernel(X, X.TEST, con$sKernelType, theta);
   PHItestSelf <- computeKernel(X.TEST, X.TEST, con$sKernelType, theta);
   Res <- t(crossprod(Y, invPHI)%*%PHItest);
   S <- (diag(PHItestSelf) - diag(crossprod(PHItest, invPHI)%*%PHItest));
   predictive.likelihood <- 0.;
   if (Kc > 2) {
      Ptest <- matrix(1., nrow=Ntest, ncol=Kc);
      u     <- rnorm(con$nSampsTG);
      for (n in 1:Ntest) {
         for (i in 1:Kc) {
            pp <- rep(1., con$nSampsTG);
            for (j in ((1:Kc)[-i])) {
               pp <- pp * safeNormCDF(u + (Res[n, i] - Res[n, j])/(sqrt(1.+S[n])));
            }
            Ptest[n, i] <- mean(pp);
         }
      }
   } else {
      stop("Multinomial only code....")
   }
   Ptest <- t(apply(Ptest, 1, function(x) {x/sum(x)})); ## JUST IN CASE
   Ptest;
}

