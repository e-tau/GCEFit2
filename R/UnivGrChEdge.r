#################
## Function for univariate extended Gram-Charlier and Edgeworth

# Edg.univ - Edgeworth first term, any n
# hermite_polynomial - recursive formula
# GC.univ - univariate GC - provide cumulants up to any order

# ExpSquGC - Expanded Squared GC density
# psi2ck - used in EXpSquGC
# psi2_constants - provide only the constants (b_k's)
# psi2_direct - used to build the squared form (not expanded)
# LLEGC2 - LL for squared GC - not extended
# LLEGC3 - LL for squared GC - extended
# Edge.dens.squ - Edgeworth squared
# LLEdge2 - Log likelihood function for squared (extended) Edgeworth
# CF_Ext_GC - Cornish Fisher for extended GG
# EM_EGC  - EM for Gyorgy's procedure Extended GC



#' Univariate Edgeworth density
#'
#' Fits the Edgeworth approximation
#'
#'
#' @param w a vector of data,
#' @param n number of terms in the  mean
#' @param cum c(muX,varX,skewX,kurtX) the cumulants of the elements of the sum
#'
#' @return Edgeworth approximation evaluated at sample points w. w can be considered
#' a mean of \eqn{n_0} independent copies of a rv X.
#'
#' @details
#'
#' Computes the Edgeworth density at the points  \eqn{w = \bar{x}_{n_0}= (X_1+ \dots + X_{n_0}}) where
#' \eqn{X_1, \dots, X_{n_0}} are i.i.d. rv with mean \eqn{\mu_X} and
#' variance \eqn{\sigma_X^2}
#' \deqn{
#' f_{{W}^{\left( n_0\right) }}\left( {w}\right) = \left( 1+\sum_{k=1}^{2}\frac{n_0^{-k/2}}{k!}g_{{Y},k}\left(
#' {z}\right) \right) \varphi \left( {w}|\mu_X,{\sigma }^2_{%
#'{X}}\right)
#'  }
#' where
#' \deqn{
#' g_{{Y},1}\left( {z}\right) =\frac{{\kappa }_{%
#' {Y},3}{H}_{3}\left( {{z}}\right) }{6}, \qquad
#' g_{{Y},2}\left( {z}\right) =\frac{{\kappa }_{%
#'  {Y},4}{H}_{4}\left( {{z}}\right) }{12}+\frac{{\kappa }_{{Y},3}^{ 2
#'  }{H}_{6}\left( {{z}}\right) }{36}.
#' }
#' Here \eqn{Y = (X-\mu_X)/\sigma_X} and \eqn{Z=\sqrt{n_0}(w-\mu_X)/\sigma_X}
#'
#' @examples
#' ## GH distribution example
#' param <- c(mu=0, delta= 2, alpha= 1.8, beta= -0.91, lambda=2)
#' ## generate GH random data
#' set.seed(123)
#' xGH <- GeneralizedHyperbolic::rghyp(1000,param=param)
#' xGH <- sort(xGH)
#'
#' ## sample cumulants values
#' Scum <- MultiStatM::SampleMomCum(xGH,r=4,centering=FALSE,scaling=FALSE)$estCum.r
#' ## Mean, variance, sekwness and excess kurtosis
#' ScumGH <- c(Scum[1],Scum[2],Scum[3]/Scum[2]^(3/2),Scum[4]/Scum[2]^(4/2))
#' ## Fit Edgeworth density n=1
#' EdgeFit <- Edge.dens(xGH, n=1, ScumGH)
#' ## plot results
#' hist(xGH, freq=FALSE, breaks=20, main = "xGH data")
#' lines(xGH,EdgeFit, col="red",lwd=2)
#'
#'
#' @export
Edge.dens <- function(w, n=1, cum) {
  # note: w is a sum or a mean of i.i.d X
  muX <- cum[1]
  varX  <- cum[2]
  skewX <- cum[3]
  kurtX <- cum[4]


           scale <- sqrt(n) / sqrt(varX)
           z <- (w - muX) * scale

           fi.z <- 1 + n^(-1/2) * skewX * hermite_polynomial(3, z)/ 6 +
             n^(-1) * (kurtX * hermite_polynomial(4, z)/ 12 + skewX^2 * hermite_polynomial(6, z)/ 36) / 2

           Edg.dens <- scale* fi.z * stats::dnorm(z)


  return(Edg.dens)
}



#' Univariate Gram Charlier density
#'
#' Fits the Gram Charlier approximation of any order to the
#' points x using the density
#' \deqn{
#' f_{{X}}\left( {x}\right) =\left(
#' 1+\sum_{k=3}^{\infty }\frac{1}{k!}{B}_{k}\left( 0,0,%
#' {\kappa }_{{Y},3},\ldots {\kappa }_{%
#' {Y},k}\right) {H}_{k}\left( {y}|{I}%
#'  \right) \right) \varphi \left( {x}|{\mu}_X,{\sigma }_X^2\right).
#' }
#'
#'
#' @param x Numeric vector of evaluation points.
#' @param cum Numeric vector of STANDARDIZED cumulants. The first element is the mean,
#'   the second the variance, and higher-order entries correspond to
#'   higher-order cumulants. Must contain at least three elements.
#' @param beta Optional numeric vector of SHRINKAGE coefficients for the generalized
#'   Charlier expansion. If `NULL`, all coefficients are set to 1.
#'   Must have length `length(cum) - 2`.
#' @param normalize Logical; if `TRUE`, the resulting density is numerically
#'   normalized using the trapezoidal rule over `x`.
#'
#' @return A numeric vector containing the approximated density evaluated
#'   at the points in `x`.
#'
#' @details
#' The function uses `GC.Design()` to construct the generalized Charlier
#' design matrix associated with the supplied cumulants. The density is then
#' obtained by multiplying the Gaussian reference density by the polynomial
#' correction term
#'
#' \deqn{
#' q(x) = 1 + H(x)\beta.
#' }
#'
#' If `normalize = TRUE`, the density is rescaled so that its numerical
#' integral over the supplied grid `x` equals one.
#'
#' @examples
#' ## GH distribution example
#' param <- c(mu=0, delta= 2, alpha= 1.8, beta= -0.91, lambda=2)
#' ## generate GH random data
#' set.seed(123)
#' xGH <- GeneralizedHyperbolic::rghyp(1000,param=param)
#' xGH <- sort(xGH)
#'
#' ## sample cumulants values
#' Scum <- MultiStatM::SampleMomCum(xGH,r=4,centering=FALSE,scaling=FALSE)$estCum.r
#' ## Mean, variance, sekwness and excess kurtosis
#' ScumGH <- c(Scum[1],Scum[2],Scum[3]/Scum[2]^(3/2),Scum[4]/Scum[2]^(4/2))
#' ## Fit Edgeworth density n=1
#' GrChFit <- GrCh.dens(xGH, ScumGH)
#' ## plot results
#' hist(xGH, freq=FALSE, breaks=20, main = "xGH data")
#' lines(xGH,GrChFit, col="red",lwd=2)
#'
#' @seealso
#' \code{\link{GC.Design}}
#'
#' @export
GrCh.dens <- function(x, cum, beta = NULL, normalize = FALSE) {
  
  cum <- as.numeric(cum)
  k <- length(cum)
  
  if (k < 3L) {
    stop("cum must contain at least mean, variance, and kappa_3.")
  }
  
  muX <- cum[1L]
  VX  <- cum[2L]
  
  if (!is.finite(VX) || VX <= 0) {
    stop("cum[2] must be a positive variance.")
  }
  
  H <- GC.Design(x, cum)
  
  if (is.null(beta)) {
    beta <- rep(1, k - 2L)
  } else {
    beta <- as.numeric(beta)
    if (length(beta) != k - 2L) {
      stop("beta must have length length(cum) - 2.")
    }
  }
  
  q <- 1 + as.numeric(H %*% beta)
  
  dens <- q * stats::dnorm(x, mean = muX, sd = sqrt(VX))
  
  if (normalize) {
    z <- sum(diff(x) * (head(dens, -1L) + tail(dens, -1L)) / 2)
    if (is.finite(z) && z > 0) {
      dens <- dens / z
    }
  }
  
  dens
}






#' Squared Gram-Charlier density of order K
#' 
#' Fits the squared polynomial Gram-Charlier density of order K
#' 
#' \deqn{
#' f_{{X}}\left( {x}\right) =\left(
#' 1+\sum_{k=3}^{\infty }\frac{1}{k!}{B}_{k}\left( 0,0,%
#' {\kappa }_{{Y},3},\ldots {\kappa }_{%
#' {Y},k}\right) {H}_{k}\left( {y}|{I}%
#'  \right) \right)^2 \varphi \left( {x}|{\mu}_X,{\sigma }_X^2\right).
#' }
#'
#' The parameters are estimated by [GrCh.squ.MLE()] or by [GrCh.minE()].
#' 
#' @param x     vector of evaluation points
#' @param  param  c(mu, s2, alpha3, alpha4, ..., alphaK)
#'
#'
#' @return the values of squared GC evaluated at x
#'
#' @examples
#' ## GH distribution example
#' param <- c(mu=0, delta= 2, alpha= 1.8, beta= -0.91, lambda=2)
#' ## generate GH random data
#' set.seed(123)
#' xGH <- GeneralizedHyperbolic::rghyp(1000,param=param)
#' xGH <- sort(xGH)
#' ## MLE estimation
#' MLE <- GrCh.squ.MLE(data=xGH,start=c(0,1,0,3), method="BFGSR")
#' ## fit the squared Gram Charlier
#' GCSqu <- GrCh.dens.squ(xGH, MLE$estimate)
#' ## plot results
#' hist(xGH, freq=FALSE, breaks=20, main = "xGH data")
#' lines(xGH, GCSqu,  col = "#1E90FF", lwd = 2, lty = 2)
#'
#' @export
GrCh.dens.squ <- function(x, param) {
  
  K <- length(param)
  
  if (K < 2) {
    stop("param must contain at least mu and s2")
  }
  
  mu <- param[1]
  s2 <- param[2]
  
  # Variance must be positive
  if (!is.finite(s2) || s2 <= 0) {
    return(rep(NA_real_, length(x)))
  }
  
  # Standardized variable
  z <- (x - mu) / sqrt(s2)
  
  # Build alpha vector: alpha1 = alpha2 = 0
  alpha <- numeric(K)
  if (K > 2) {
    alpha[3:K] <- param[3:K]
  }
  
  # Normalization constant
  tau_val <- 1
  if (K > 2) {
    for (i in 3:K) {
      tau_val <- tau_val + alpha[i]^2 / factorial(i)
    }
  }
  tau_val <- 1 / tau_val
  
  if (!is.finite(tau_val) || tau_val <= 0) {
    return(rep(NA_real_, length(x)))
  }
  
  # Polynomial part: 1 + sum alpha_i H_i(z)/i!
  poly <- rep(1, length(x))
  if (K > 2) {
    for (i in 3:K) {
      poly <- poly + alpha[i] * hermite_polynomial(i, z) / factorial(i)
    }
  }
  
  # Density
  dens <- tau_val * (poly^2) * stats::dnorm(x, mean = mu, sd = sqrt(s2))
  
  # safeguard
  dens[!is.finite(dens)] <- NA_real_
  
  return(dens)
}






#' Squared Edgeworth density
#'
#' Fits the Edgeworth approximation with squared polynomial part.
#' Only the case \eqn{n_0=1} is considered. Parameters are estimated
#' by [Edge.squ.MLE()]
#'
#' @param w data vector
#' @param cum  c(muX,varX,skewX,kurtX)
#'
#' @return density evaluate at w
#'
#'
#' @export
Edge.dens.squ <- function(w,cum) {
  n=1
  muX <- cum[1]
  varX  <- cum[2]
  skewX <- cum[3]
  kurtX <- cum[4]

  tau <- ( 1 + skewX^2 / 6 + kurtX^2 / 6 + skewX^4 * 5 / 9 )^(-1)
           scale <- sqrt(n) / sqrt(varX)
           z <- (w - muX) * scale

           fi.z <- 1 + n^(-1/2) * skewX * hermite_polynomial(3, z)/ 6 +
             n^(-1) * (kurtX * hermite_polynomial(4, z)/ 12 + skewX^2 * hermite_polynomial(6, z)/ 36) / 2

           fi.z <- tau * fi.z^2

           Edg.dens <- scale* fi.z * stats::dnorm(z)


  return(Edg.dens)
}




#' Integrated  Gram Charlier
#'
#' Provide the integral for the Gram Charlier [GrCh.dens()] density truncated
#' at k=4 over the region (-infinity,x)
#'
#' @param x vector of upper limits of integration
#' @param cum vector (mean, variance, skewness, ex-kurtosis)
#'
#' @return the vector of estimated probabilities
#'
#' @examples
#' ## GH distribution example
#' param <- c(mu=0, delta= 2, alpha= 1.8, beta= -0.91, lambda=2)
#' ## generate GH random data
#' set.seed(123)
#' xGH <- GeneralizedHyperbolic::rghyp(1000,param=param)
#' xGH <- sort(xGH)
#'
#' ## sample cumulants values
#' Scum <- MultiStatM::SampleMomCum(xGH,r=4,centering=FALSE,scaling=FALSE)$estCum.r
#' ## Mean, variance, sekwness and excess kurtosis
#' ScumGH <- c(Scum[1],Scum[2],Scum[3]/Scum[2]^(3/2),Scum[4]/Scum[2]^(4/2))
#' ## estimate probabilities on (-infinity,q)
#' q <- c(-1, -2, 0, 1)
#' GrCh.Int(q, ScumGH)
#'
#' @export
GrCh.Int <- function(x,cum){
  sigma <- sqrt(cum[2])
  mu <- cum[1]
  z <- (x-mu)/sigma
  Fx <- stats::pnorm(z) - stats::dnorm(z) * (cum[3] * EQL:: hermite(z,2) / 6 +
    cum[4] * EQL:: hermite(z,3) / 24 )
  return(Fx)
}


#' Integrated  Edgeworth density
#'
#' Provide the integral for the Edgeworth density truncated at k=4
#'
#' Provide the integral for the Edgeworth [Edge.dens()] density truncated
#' at k=4 over the region (-infinity,x)
#'
#' @param x vector of upper limits of integration
#' @param cum vector (mean, variance, skewness, ex-kurtosis)
#'
#' @return the vector of estimated probabilities
#'
#' @examples
#' ## GH distribution example
#' param <- c(mu=0, delta= 2, alpha= 1.8, beta= -0.91, lambda=2)
#' ## generate GH random data
#' set.seed(123)
#' xGH <- GeneralizedHyperbolic::rghyp(1000,param=param)
#' xGH <- sort(xGH)
#'
#' ## sample cumulants values
#' Scum <- MultiStatM::SampleMomCum(xGH,r=4,centering=FALSE,scaling=FALSE)$estCum.r
#' ## Mean, variance, sekwness and excess kurtosis
#' ScumGH <- c(Scum[1],Scum[2],Scum[3]/Scum[2]^(3/2),Scum[4]/Scum[2]^(4/2))
#' ## estimate probabilities on (-infinity,q)
#' q <- c(-1, -2, 0, 1)
#' Edge.Int(q, ScumGH)
#'
#' @export
Edge.Int <- function(x,cum){
  sigma <- sqrt(cum[2])
  mu <- cum[1]
  z <- (x-mu)/sigma
  Fx <- stats::pnorm(z) -  stats::dnorm(z) * ( cum[3] * EQL:: hermite(z,2) / 6 +
    cum[4] * EQL:: hermite(z,3) / 24 + (cum[3])^2 * EQL:: hermite(z,5) / 72 )
  return(Fx)
}

#' MLE estimation of squared Gram Charlier density
#'
#' it uses maxLik to produce MLE estimates for [GrCh.dens.squ()]
#'
#' @param data vector of data
#' @param start c(mu, sig2, alpha3, alpha4)
#' @param method optimization methods available in maxLik
#' @param ... Additional arguments passed to maxLik()
#'
#' @return an object from maxLik
#'
#' @examples
#' ## GH distribution example
#' param <- c(mu=0, delta= 2, alpha= 1.8, beta= -0.91, lambda=2)
#' ## generate GH random data
#' set.seed(123)
#' xGH <- GeneralizedHyperbolic::rghyp(1000,param=param)
#' xGH <- sort(xGH)
#' ## MLE estimation
#' MLE <- GrCh.squ.MLE(data=xGH,start=c(0,1,0,3), method="BFGSR")
#' ## fit the squared Gram Charlier
#' GCSqu <- GrCh.dens.squ(xGH, MLE$estimate)
#' ## plot results
#' hist(xGH, freq=FALSE, breaks=20, main = "xGH data")
#' lines(xGH, GCSqu,  col = "#1E90FF", lwd = 2, lty = 2)
#'
#'
#' @export
GrCh.squ.MLE <- function(data, start, method = "BFGSR", ...) {

  if (!requireNamespace("maxLik", quietly = TRUE)) {
    stop("Package 'maxLik' is required.")
  }

  fit <- maxLik::maxLik(
    logLik = LLGC2_K,
    start = start,
    method = method,
    data = data,
    ...
  )

  return(fit)
}



#' MLE estimation of squared Edgeworth density
#'
#' it uses maxLik to produce MLE estimates for [Edge.dens.squ()]
#'
#' @param data vector of data
#' @param start c(mu, sig2, alpha3, alpha4)
#' @param method optimization methods available in maxLik
#' @param ... Additional arguments passed to maxLik()
#'
#' @return an object from maxLik
#'
#' @examples
#' ## GH distribution example
#' param <- c(mu=0, delta= 2, alpha= 1.8, beta= -0.91, lambda=2)
#' ## generate GH random data
#' set.seed(123)
#' xGH <- GeneralizedHyperbolic::rghyp(1000,param=param)
#' xGH <- sort(xGH)
#' ## MLE estimation
#' MLE <- Edge.squ.MLE(data=xGH,start=c(0,1,0,3), method="BFGSR")
#' ## fit the squared Edgeworth
#' EdgeSqu <- Edge.dens.squ(xGH, MLE$estimate)
#' ## plot results
#' hist(xGH, freq=FALSE, breaks=20, main = "xGH data")
#' lines(xGH, EdgeSqu,  col = "#1E90FF", lwd = 2, lty = 2)
#'

#'
#' @export
Edge.squ.MLE <- function(data, start, method = "BFGSR", ...) {

  if (!requireNamespace("maxLik", quietly = TRUE)) {
    stop("Package 'maxLik' is required.")
  }

  fit <- maxLik::maxLik(
    logLik = LLEdge2,
    start = start,
    method = method,
    data = data,
    ...
  )

  return(fit)
}




#' Integrated Squared Gram Charlier density 
#' 
#' Provides the integral for any order k of the suared polynomial density
#' 
#' \deqn{
#' f_{{X}}\left( {x}\right) =\left(
#' 1+\sum_{k=3}^{\infty }\frac{1}{k!}{B}_{k}\left( 0,0,%
#' {\kappa }_{{Y},3},\ldots {\kappa }_{%
#' {Y},k}\right) {H}_{k}\left( {y}|{I}%
#'  \right) \right)^2 \varphi \left( {x}|{\mu}_X,{\sigma }_X^2\right).
#' }
#'
#' The parameters are estimated by [GrCh.squ.MLE()] or [GrCh.minE()]
#' 
#' 
#' @param x vector of values
#' @param param vector c(mu, sigma2, a3, ... ak)
#' 
#' @export
GrCh.Int.squ <- function(x,param){
  k=length(param)
  mu <- param[1]
  sigma2 <- param[2]
  #alpha3 <- param[3]
  #alpha4 <- param[4]
  Ck <- HermSqu2(param)
  sigma <- sqrt(sigma2)
  z <- (x-mu)/sigma
  C28 <- 0
  for (i in 2:(2*k)){
    C28 <- C28 + Ck[i] * EQL:: hermite(z,i-1) * stats::dnorm(z)
  }
  Fx <- stats::pnorm(z) - Ck[1] * stats::dnorm(z) - C28
  return(Fx)
}

