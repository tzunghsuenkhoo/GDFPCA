# Generalized dynamic principal component analysis 
#' Title
#'
#' @param X - An fd object (see fda package) 
#' @param p - the number of GDFC (generalized dynamic functional principal components) 
#' @param SCALE_ID - scaling of the coefficients. Default = FALSE
#' @param center - "center=TRUE" centers the fd object. Default = TRUE
#'
#' @return 
#' @export
#'
#' @examples
GDFPCA <- function(X, p, SCALE_ID = FALSE, center=TRUE){
  #Get the means for reconstruction purpose
  mean_fd = mean.fd(X)
  mean_X = rep(NA,dim(X$coefs)[1])
  for (i in 1:dim(X$coefs)[2]){
    mean_X = cbind(mean_X,mean_fd$coefs)
  }
  mean_X = mean_X[,-1]
  mean_X = t(mean_X)
  # X is fd object 
  if(center==TRUE){
    X = center.fd(X)
  }
  X.gdpca    <- auto.gdpc(scale(t(X$coefs), center = TRUE, scale = SCALE_ID), auto_comp = FALSE, normalize = 1, num_comp = p, niter_max = 1000)
  
  Xhat.gdpca <- scale(fitted(X.gdpca, num_comp = p), center = FALSE, scale = `if`(SCALE_ID,1/X.pca$scale,FALSE))
  Xhat.gdpca <- scale(Xhat.gdpca, center = -X.pca$center, scale = FALSE)
  
  Xhat.gdpca.fd  <- fd(coef = t(Xhat.gdpca + mean_X), X$basis)
  
  nmse.gdpca <- sum((Xhat.gdpca - t(X$coefs))^2)/(dim(X$coefs)[1]*dim(X$coefs)[2])
  
  # Empirical explained variance
  var.gdpca = (1 - sum( (Xhat.gdpca - t(X$coefs))**2 ) / sum(X$coefs**2))*100
  
  lst <-  list(RMSE = nmse.gdpca,
               VAR  = var.gdpca,
               XHAT = Xhat.gdpca.fd,
               PCAs = X_GDPCA = X.gdpca
  )
  return(lst)  
}

