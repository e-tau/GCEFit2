



# Univariate Hermyte polynomial
#
# @param x a vector of data
# @param n order of polynomial
#
# @return The vector of polynomial values
hermite_polynomial <- function(n, x) {
  if (n == 0) {
    return(1)
  } else if (n == 1) {
    return(x)
  } else if (n == 2) {
    return(x^2 - 1)
  } else {
    return(x * hermite_polynomial(n - 1, x) - (n - 1) * hermite_polynomial(n - 2, x))
  }
}


## Corrected psi squared from page 73 - formula 3.4
## chooose the numbero f k's
psi2ck <- function(x, mu, sig2,alpha_3, alpha_4,k=8) {
  # Calculate tau (as defined in your previous equations)
  tau <- 1 / (1 + alpha_3^2 / 6 + alpha_4^2 / 24)


  # constants
  C0 <- 1 + (alpha_3^2/6 + alpha_4^2/24)
  C1 <- (1/3) * alpha_3 * alpha_4
  C2 <- (alpha_3^2/2 + alpha_4^2/6)
  C3 <- ( alpha_3 / 3 +  alpha_3 * alpha_4/ 2)
  C4 <- (alpha_4 / 12 + alpha_3^2 / 4 + alpha_4^2 / 8)
  C5 <-  alpha_3 * alpha_4 / 6
  C6 <- ((1/36) * alpha_3^2  + alpha_4^2 / 36 )
  C7 <- (1/(72)) * alpha_3 * alpha_4
  C8 <- (1/24^2) * alpha_4^2

  Const <- c(C1,C2,C3,C4,C5,C6,C7,C8)

  #Compute the terms of the formula
  z <- (x-mu)/sqrt(sig2)

  gz <- C0

  for (i in 1:k){
    # gz <- gz + Const[i]* hermite_polynomial(i, z)/sig2^(i/2) ## origni
    gz <- gz + Const[i]* hermite_polynomial(i, z)
  }
  gz <- tau*gz
  return(gz)
}


## provide only the constants (b_k's)
psi2_constants <- function(mu, sig2,alpha_3, alpha_4,k=8) {
  # Calculate tau (as defined in your previous equations)
  tau <- 1 / (1 + alpha_3^2 / 6 + alpha_4^2 / 24)


  # constants
  C0 <- 1 + (alpha_3^2/6 + alpha_4^2/24)
  C1 <- (1/3) * alpha_3 * alpha_4
  C2 <- (alpha_3^2/2 + alpha_4^2/6)
  C3 <- ( alpha_3 / 3 +  alpha_3 * alpha_4/ 2)
  C4 <- (alpha_4 / 12 + alpha_3^2 / 4 + alpha_4^2 / 8)
  C5 <-  alpha_3 * alpha_4 / 6
  C6 <- ((1/36) * alpha_3^2  + alpha_4^2 / 36 )
  C7 <- (1/(72)) * alpha_3 * alpha_4
  C8 <- (1/24^2) * alpha_4^2

  Const <- c(C1,C2,C3,C4,C5,C6,C7,C8) *tau

  return(Const)
}



psi2_direct <- function(x,mu,sd, alpha_3, alpha_4) {
  # Calculate tau (as defined in your previous equations)
  tau <- 1 / (1 + alpha_3^2 / 6 + alpha_4^2 / 24)

  # Calculate the Hermite polynomials using the recursive function

  z <- (x-mu)/sd

  H3 <- hermite_polynomial(3, z)
  H4 <- hermite_polynomial(4, z)


  # Compute the terms of the formula
  term1 <- tau*(1 + alpha_3/6 * H3 + alpha_4 * H4 /24 )^2


  # Sum the terms and multiply by 1/tau

  return(term1)
}


# Log-likelihood for squared (direct) Gram-Charlier
#
# Builds the LL for the squared (not extended) GC. Version used
# in the literature
#
# Used to feed the maxLik function
#
# @param data vector of data
# @param param (mean, var, alpha3, alpha4)
#
#
#
# @return Scalar log-likelihood
# @export
LLEGC2 <- function(param, data) {

  mu  <- param[1]
  s2  <- param[2]
  k3  <- param[3]
  k4  <- param[4]

  # Variance must be positive
  if (!is.finite(s2) || s2 <= 0)
    return(-1e12)

  z <- (data - mu) / sqrt(s2)

  # Hermite expansion
  poly <- 1 +
    k3 * hermite_polynomial(3, z) / 6 +
    k4 * hermite_polynomial(4, z) / 24

  tau <- (1 + k3^2 / 6 + k4^2 / 24)^(-1)

  if (!is.finite(tau) || tau <= 0)
    return(-1e12)

  # Compute in log-space (numerically stable)
  log_base  <- stats::dnorm(data, mean = mu, sd = sqrt(s2), log = TRUE)
  log_poly2 <- 2 * log(abs(poly) + 1e-12)  # small safeguard
  log_tau   <- log(tau)

  ll <- sum(log_base + log_poly2 + log_tau)

  if (!is.finite(ll))
    return(-1e12)

  return(ll)
}


# Log likelihood function for squared (extended) Edgeworth
#
# Builds the LL for the squared (not extended) Edgeworth (consider only (n=1)
#
# Used to feed the maxLik function
#
# @param data vector of data
# @param param (mean, var, alpha3, alpha4) Note var not sd
#
#
# @export
LLEdge2 <- function(data,param){
  de <- Edge.dens.squ(data, param)
  LL <- sum(log(de))
  return(LL)
}



