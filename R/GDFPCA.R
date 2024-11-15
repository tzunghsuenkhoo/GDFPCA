#' Generalized dynamic principal component analysis 
#' @importFrom fda mean.fd center.fd is.fd
#' @importFrom gdpc auto.gdpc 
#' @importFrom stats fitted
#' @importFrom funData funData2fd simFunData
#' @param X         An fd object (see fda package) 
#' @param p         Integer. The number of GDFPC (generalized dynamic functional principal components) kept to approximate the original process.  Default is 1. 
#' @param center    Logical. If TRUE, centers the fd object. Default = TRUE
#' @param auto_comp Logical. If TRUE, computes further GDFPCs until the proportion of explained variance reaches expl_var.
#'                  Otherwise, use p. Default = FALSE.
#' @param expl_var  A number between 0 and 1. Desired proportion of explained variance (used only if auto_comp==TRUE). Default is 0.8
#' @param niter_max Integer. Maximum number of iterations. Default is 500
#' 
#' @return  \item{Xhat}{The approximation of X with GDPCs.}
#'   \item{GDFPC}{Estimated \eqn{p \times (n+k)} matrix of GDPCs with k lags}
#'   \item{GDPCA_coefs}{Estimated coefficients}
#'   \item{NMSE}{The normalized mean squared error with p number of GDPCs}
#'   \item{VAR}{The numerical explained variance with p number of GDPCs}
#'
#' 
#' @export
#' @references Khoo, T. H., Dabo, I. M., Pathmanathan, D., & Dabo-Niang, S. (2024). Generalized functional dynamic principal component analysis. arXiv preprint arXiv:2407.16024.
#' Peña, D., & Yohai, V. J. (2016). Generalized dynamic principal components. Journal of the American Statistical Association, 111(515), 1121-1131.
#'
#' @examples 
#'      # Create 50 sample paths of a Wiener process by using its truncated Karhunen-Loève expansion 
#'      
#'      argvals = seq(0, 1, length = 100)
#'      efW <- funData::simFunData(seq(0, 1, length = 100),M = 100 ,eFunType = "Wiener", eValType = "wiener", N = 50)
#'      efW.fd <- funData::funData2fd(efW$simData)
#'      
#'      # Plot the functional process
#'      plot(efW.fd)
#'      
#'      # Perform GDFPCA on the process
#'      Res <- GDFPCA(X = efW.fd, p = 1)
#'      
#'      # Plot the approximated process
#'      plot(Res$Xhat)
#'      
#'      # Plot 1st GDFPC
#'      plot(Res$GDFPC)
#'      
GDFPCA <- function(X, p = 1, center = TRUE, auto_comp = FALSE, expl_var = 0.8,  niter_max = 500){
  # Check the class of X
  if(is.fd(X) == FALSE){
    stop("X is not an fd object")
  }
  #Get the mean function (scores) for reconstruction purpose
  mean_fd = fda::mean.fd(X)
  mean_X  = rep(NA,dim(X$coefs)[1])
  for (i in 1:dim(X$coefs)[2]){
    mean_X = cbind(mean_X,mean_fd$coefs)
  }
  mean_X  = mean_X[,-1]
  mean_X  = t(mean_X)
  # Center X if center == TRUE
  if(center==TRUE){
    X = fda::center.fd(X)
  }
  # Calculate GDFPCs and its approximation
  X.gdpca        <- gdpc::auto.gdpc(t(X$coefs), num_comp = p, auto_comp = FALSE, niter_max = 500, expl_var = 0.8,
                                    tol = 1e-4, k_max = 10, normalize = 1,
                                    ncores = 1, verbose = FALSE)
  Xhat.gdpca     <- stats::fitted(X.gdpca, num_comp = p)
  Xhat.gdpca.fd  <- fda::fd(coef = t(Xhat.gdpca + mean_X), X$basis)
  
  # Normalized mean squared error
  nmse.gdpca <- sum((Xhat.gdpca - t(X$coefs))^2)/(dim(X$coefs)[1]*dim(X$coefs)[2])
  # Empirical explained variance
  var.gdpca = (1 - sum( (Xhat.gdpca - t(X$coefs))**2 ) / sum(X$coefs**2))*100
  
  lst <-  list(Xhat        = Xhat.gdpca.fd,
               GDFPC       = X.gdpca[[1]]$f,
               GDPCA_coefs = X.gdpca,
               NMSE        = nmse.gdpca,
               VAR         = var.gdpca
  )
  return(lst)  
}

