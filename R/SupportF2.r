## builds the polynomial part of GC allowing shrinkage coefficients
GC.Design <- function(x, cum) {
  
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
  
  z <- (x - muX) / sqrt(VX)
  
  kappa <- cum
  kappa[1L] <- 0
  kappa[2L] <- 0
  
  B <- MultiStatM::Cum2Mom(kappa, Type = "Univariate")
  
  orders <- 3:k
  
  H <- sapply(orders, function(i) {
    B[i] * hermite_polynomial(i, z) / factorial(i)
  })
  
  if (length(orders) == 1L) {
    H <- matrix(H, ncol = 1L)
  }
  
  colnames(H) <- paste0("H", orders)
  
  H
}



# Log-likelihood for squared Gram-Charlier of order K
#
# param = c(mu, s2, alpha3, alpha4, ..., alphaK)
# data  = vector of observations
#
# Returns scalar log-likelihood
LLGC2_K <- function(param, data) {
  
  K <- length(param)
  
  if (K < 2) {
    stop("param must contain at least mu and s2")
  }
  
  mu <- param[1]
  s2 <- param[2]
  
  # Variance must be positive and finite
  if (!is.finite(s2) || s2 <= 0) {
    return(-1e12)
  }
  
  # Standardized data
  z <- (data - mu) / sqrt(s2)
  
  # Build alpha vector:
  # alpha[1] = alpha1 = 0
  # alpha[2] = alpha2 = 0
  # alpha[3], ..., alpha[K] taken from param[3:K]
  alpha <- numeric(K)
  if (K > 2) {
    alpha[3:K] <- param[3:K]
  }
  
  # Normalizing constant
  tau_val <- 1
  if (K > 2) {
    for (i in 3:K) {
      tau_val <- tau_val + alpha[i]^2 / factorial(i)
    }
  }
  tau_val <- 1 / tau_val
  
  if (!is.finite(tau_val) || tau_val <= 0) {
    return(-1e12)
  }
  
  # Polynomial term: 1 + sum_{i=3}^K alpha_i H_i(z)/i!
  poly <- rep(1, length(data))
  if (K > 2) {
    for (i in 3:K) {
      poly <- poly + alpha[i] * hermite_polynomial(i, z) / factorial(i)
    }
  }
  
  # Log-density in log-space
  log_base  <- stats::dnorm(data, mean = mu, sd = sqrt(s2), log = TRUE)
  log_poly2 <- 2 * log(abs(poly) + 1e-12)   # safeguard
  log_tau   <- log(tau_val)
  
  ll <- sum(log_base + log_poly2 + log_tau)
  
  if (!is.finite(ll)) {
    return(-1e12)
  }
  
  return(ll)
}


## --------------------------
HermMom2 <- function(params, cum){
  ## left side 
  ## params[1] <- MLE$estimate[1]
  ## params[2] <- log(MLE$estimate[2])
  K <- length(params)
  Bvec <- NULL
  theta <- NULL
  theta[1] <- cum[1]- params[1]
  theta[2] <- cum[2] - exp(params[2]) # + (cum[1]- params[1])^2
  
  Bvec <- theta[1]
  Bpoly <- bell_complete(1,Bvec)
  Bvec <- c(Bvec, theta[2] )
  Bpoly <- c(Bpoly, bell_complete(2,Bvec))
  for (k in 3:(K)) {
    Bvec <- c(Bvec,cum[k]) #*(1+theta[2])^((k+2)/2))
    Bpoly <- c(Bpoly, bell_complete(k,Bvec))
  }
  return(Bpoly)
}


#####################
HermSqu2 <- function(params){
  K <- length(params)
  alpha <-c(1,0,0,params[3:K]) # dim K+1
  tau_val <- 1
  for (k in 4:(K+1)) {
    tau_val <- tau_val+ alpha[k]^2/factorial(k-1)
  }
  tau_val <- tau_val^(-1)
  
  T <- hermite_tensor(K)
  beta <- NULL
  for (ell in 1:(2*K)) {
    T_ell <- factorial(ell)*T[,,ell+1]
    tilde_c_ell <- tau_val*as.numeric( t(alpha) %*% T_ell %*% alpha )
    beta <-c( beta, tilde_c_ell) 
  }
  beta <- beta / factorial(1:(2*K))
  return(beta)
}

################################# bell polynomial
bell_complete <- function(n, x) {
  # x should be a vector: x1, x2, ..., xn
  B <- numeric(n + 1)
  B[1] <- 1  # B0 = 1
  
  for (m in 1:n) {
    B[m + 1] <- 0
    for (k in 1:m) {
      B[m + 1] <- B[m + 1] +
        choose(m - 1, k - 1) * x[k] * B[m - k + 1]
    }
  }
  
  return(B[n + 1])
}
#################################   ####################################### 

hermite_tensor <- function(K) {
  # Allocate tensor of size (K+1) x (K+1) x (2K+1)
  T <- array(0, dim = c(K+1, K+1, 2*K+1))
  
  for (m in 0:K) {
    for (n in 0:K) {
      # Max r is limited by min(m,n)
      rmax <- min(m, n)
      for (r in 0:rmax) {
        ell <- m + n - 2*r
        # Only store ell in range 0..2K
        if (ell >= 0 && ell <= 2*K) {
          T[m+1, n+1, ell+1] <- T[m+1, n+1, ell+1] +
            choose(m, r) * choose(n, r) * factorial(r)/factorial(m)/factorial(n)  
        }
      }
    }
  }
  
  return(T)
}

###################
HermSqu <- function(params){
  K <- length(params)
  alpha <-c(1,0,0,params[3:K]) # dim K+1
  tau_val <- 1
  for (k in 4:(K+1)) {
    tau_val <- tau_val+ alpha[k]^2/factorial(k-1)
  }
  tau_val <- tau_val^(-1)
  
  T <- hermite_tensor(K)
  beta <- NULL
  for (ell in 1:(K)) {
    T_ell <- factorial(ell)*T[,,ell+1]
    tilde_c_ell <- tau_val*as.numeric( t(alpha) %*% T_ell %*% alpha )
    beta <-c( beta, tilde_c_ell) 
  }
  return(beta)
}


HermMom <- function(params, cum){
  ####### left side 
  Bvec <- NULL
  theta <- NULL
  theta[1] <- params[1]
  theta[2] <- params[2]
  
  Bvec <- theta[1]
  Bpoly <- bell_complete(1,Bvec)
  Bvec <- c(Bvec, theta[2] )
  Bpoly <- c(Bpoly, bell_complete(2,Bvec))
  for (k in 1:(K-2)) {
    Bvec <- c(Bvec,cum[k]*(1+theta[2])^((k+2)/2)) 
    Bpoly <- c(Bpoly, bell_complete(k+2,Bvec))
  }
  return(Bpoly)
}



