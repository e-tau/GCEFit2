#' Positive shrinkage for Gram-Charlier correction factors
#'
#' Fits shrinkage coefficients for a truncated Gram-Charlier expansion so that
#' the resulting correction factor is nonnegative on a user-supplied grid.
#'
#' The cumulant vector is assumed to have the form
#'
#' \deqn{
#'   cum = (\kappa_1, \kappa_2, \kappa_3, \ldots, \kappa_K),
#' }
#'
#' where \eqn{\kappa_1} is the mean, \eqn{\kappa_2} is the variance, and
#' \eqn{\kappa_3,\ldots,\kappa_K} are standardized higher-order cumulant
#' quantities entering the Gram-Charlier correction factor
#'
#' \deqn{
#'   q_\beta(y)
#'   =
#'   1 + \sum_{k=3}^K \beta_k
#'   \frac{\kappa_k}{k!} H_k(y),
#'   \qquad
#'   y = \frac{x - \kappa_1}{\sqrt{\kappa_2}}.
#' }
#'
#' The function solves the weighted quadratic shrinkage problem
#'
#' \deqn{
#'   \min_\beta \sum_{k=3}^K w_k (1 - \beta_k)^2
#' }
#'
#' subject to
#'
#' \deqn{
#'   q_\beta(y_i) \ge 0
#' }
#'
#' on the grid supplied through `grid_x`, together with box constraints
#' \eqn{0 \le \beta_k \le 1}. Larger weights make shrinkage of the
#' corresponding coefficient less attractive. Individual coefficients may also
#' be fixed exactly using `fixed_orders` and `fixed_beta`.
#'
#' @param grid_x Numeric vector. Grid of `x` values where nonnegativity of the
#'   Gram-Charlier correction factor is enforced.
#' @param cum Numeric vector of cumulant inputs. The required format is
#'   `cum = c(mean, variance, skewness, excess_kurtosis, ...)`.
#'   The first entry is the mean, the second entry is the variance, and entries
#'   from the third onward are standardized higher-order quantities used in the
#'   Hermite expansion.
#' @param fixed_orders Optional integer vector. Orders whose shrinkage
#'   coefficients are fixed rather than optimized. For example,
#'   `fixed_orders = c(3, 4)` fixes the skewness and kurtosis shrinkage
#'   coefficients.
#' @param fixed_beta Optional numeric vector with the same length as
#'   `fixed_orders`. Values assigned to the fixed coefficients. If `NULL`,
#'   all fixed coefficients are set to 1.
#' @param weights Optional numeric vector of length `length(cum) - 2`.
#'   The entries correspond to orders `3:length(cum)`. Larger weights penalize
#'   shrinkage more heavily. If `NULL`, all weights are set to 1.
#' @param tol Numeric tolerance used in feasibility checks.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{beta}{
#'     Named numeric vector of fitted shrinkage coefficients. The entries are
#'     named `beta_3`, `beta_4`, ..., `beta_K`. Values close to 1 indicate
#'     little or no shrinkage. Values close to 0 indicate strong shrinkage.
#'     If the optimization is infeasible, non-fixed entries are returned as
#'     `NA`.
#'   }
#'   \item{orders}{
#'     Integer vector containing the Hermite orders associated with `beta`,
#'     namely `3:length(cum)`.
#'   }
#'   \item{feasible}{
#'     Logical value. `TRUE` means that the returned coefficients satisfy
#'     nonnegativity of the correction factor on `grid_x`, up to numerical
#'     tolerance. `FALSE` indicates that no feasible solution was found or that
#'     the fixed constraints were incompatible with positivity on the grid.
#'   }
#'   \item{min_factor}{
#'     Minimum value of the fitted correction factor over `grid_x`.
#'     A value near zero means that at least one grid point is binding.
#'     Negative values indicate numerical or optimization failure.
#'   }
#'   \item{q_grid}{
#'     Numeric vector containing the fitted Gram-Charlier correction factor
#'     evaluated on `grid_x`.
#'   }
#'   \item{grid_x}{
#'     The grid supplied by the user. Returned for convenience.
#'   }
#'   \item{objective_value}{
#'     Value of the quadratic objective returned by `quadprog::solve.QP`.
#'     Smaller values indicate less total weighted shrinkage. This component is
#'     absent when all coefficients are fixed or when the solver fails.
#'   }
#'   \item{call}{
#'     The matched function call.
#'   }
#' }
#'
#' @details
#' The positivity constraint is imposed only on the supplied grid. For numerical
#' safety, it is usually advisable to fit on a dense grid and then verify the
#' fitted correction factor on a finer grid.
#'
#' The distinction between high weights and fixed coefficients is important.
#' Large weights create a soft preference for preserving selected coefficients,
#' whereas `fixed_orders` imposes exact equality constraints. For example, using
#' very large weights on orders 3 and 4 keeps skewness and kurtosis nearly
#' unchanged unless positivity requires otherwise. By contrast,
#' `fixed_orders = c(3, 4)` forces `beta_3 = beta_4 = 1`; this can make the
#' problem infeasible.
#'
#' @examples
#' cum <- c(
#'   0.0,   # mean
#'   1.0,   # variance
#'   0.8,   # skewness
#'   1.5,   # excess kurtosis
#'   0.7,   # 5th standardized cumulant-like term
#'   0.4,   # 6th standardized cumulant-like term
#'   0.25,  # 7th standardized cumulant-like term
#'   0.15   # 8th standardized cumulant-like term
#' )
#'
#' x_grid <- seq(
#'   cum[1] - 4 * sqrt(cum[2]),
#'   cum[1] + 4 * sqrt(cum[2]),
#'   length.out = 1001
#' )
#'
#' # Example 1: all coefficients are shrunk with equal weights
#' fit_all <- GrCh.shrink(
#'   grid_x = x_grid,
#'   cum = cum
#' )
#'
#' fit_all$beta
#' fit_all$feasible
#' fit_all$min_factor
#'
#' # Example 2: skewness and kurtosis are strongly protected
#' # through very large weights, but are not fixed exactly
#' fit_weighted <- GrCh.shrink(
#'   grid_x = x_grid,
#'   cum = cum,
#'   weights = c(1e4, 1e4, 1, 1, 1, 1)
#' )
#'
#' fit_weighted$beta
#' fit_weighted$feasible
#' fit_weighted$min_factor
#'
#' # Example 3: skewness and kurtosis shrinkage coefficients
#' # are fixed exactly at beta_3 = beta_4 = 1
#' fit_fixed <- GrCh.shrink(
#'   grid_x = x_grid,
#'   cum = cum,
#'   fixed_orders = c(3, 4),
#'   fixed_beta = c(1, 1)
#' )
#'
#' fit_fixed$beta
#' fit_fixed$feasible
#' fit_fixed$min_factor
#'
#' @export
GrCh.shrink <- function(grid_x,
                        cum,
                        fixed_orders = integer(0),
                        fixed_beta = NULL,
                        weights = NULL,
                        tol = 1e-10) {
  
  grid_x <- as.numeric(grid_x)
  cum <- as.numeric(cum)
  
  if (length(cum) < 3L) {
    stop("`cum` must contain at least mean, variance, and kappa_3.", call. = FALSE)
  }
  
  if (!is.finite(cum[2L]) || cum[2L] <= 0) {
    stop("`cum[2]` must be a positive variance.", call. = FALSE)
  }
  
  H <- GC.Design(grid_x, cum)
  
  p <- ncol(H)
  orders <- seq.int(3L, length(cum))
  
  if (is.null(weights)) {
    weights <- rep(1, p)
  } else {
    weights <- as.numeric(weights)
    
    if (length(weights) != p) {
      stop("`weights` must have length `length(cum) - 2`.", call. = FALSE)
    }
    
    if (any(!is.finite(weights)) || any(weights <= 0)) {
      stop("All `weights` must be finite and strictly positive.", call. = FALSE)
    }
  }
  
  if (length(fixed_orders) > 0L) {
    fixed_orders <- as.integer(fixed_orders)
    fixed_idx <- match(fixed_orders, orders)
    
    if (anyNA(fixed_idx)) {
      stop("`fixed_orders` must be among `3:length(cum)`.", call. = FALSE)
    }
    
    if (is.null(fixed_beta)) {
      fixed_beta <- rep(1, length(fixed_orders))
    } else {
      fixed_beta <- as.numeric(fixed_beta)
      
      if (length(fixed_beta) != length(fixed_orders)) {
        stop("`fixed_beta` must have the same length as `fixed_orders`.", call. = FALSE)
      }
      
      if (any(!is.finite(fixed_beta))) {
        stop("All values in `fixed_beta` must be finite.", call. = FALSE)
      }
    }
  } else {
    fixed_idx <- integer(0)
    fixed_beta <- numeric(0)
  }
  
  free_idx <- setdiff(seq_len(p), fixed_idx)
  
  offset <- rep(1, length(grid_x))
  
  if (length(fixed_idx) > 0L) {
    offset <- offset + as.numeric(H[, fixed_idx, drop = FALSE] %*% fixed_beta)
  }
  
  if (length(free_idx) == 0L) {
    beta_full <- rep(NA_real_, p)
    beta_full[fixed_idx] <- fixed_beta
    names(beta_full) <- paste0("beta_", orders)
    
    q_grid <- offset
    
    return(list(
      beta = beta_full,
      orders = orders,
      feasible = all(q_grid >= -tol),
      min_factor = min(q_grid),
      q_grid = q_grid,
      grid_x = grid_x,
      call = match.call()
    ))
  }
  
  H_free <- H[, free_idx, drop = FALSE]
  q <- length(free_idx)
  
  Amat <- cbind(
    t(H_free),
    diag(q),
    -diag(q)
  )
  
  bvec <- c(
    -offset,
    rep(0, q),
    rep(-1, q)
  )
  
  Dmat <- diag(weights[free_idx], q, q)
  dvec <- weights[free_idx]
  
  Dmat <- Dmat + diag(1e-12, q)
  
  sol <- tryCatch(
    quadprog::solve.QP(
      Dmat = Dmat,
      dvec = dvec,
      Amat = Amat,
      bvec = bvec
    ),
    error = function(e) NULL
  )
  
  if (is.null(sol)) {
    beta_full <- rep(NA_real_, p)
    
    if (length(fixed_idx) > 0L) {
      beta_full[fixed_idx] <- fixed_beta
    }
    
    names(beta_full) <- paste0("beta_", orders)
    
    return(list(
      beta = beta_full,
      orders = orders,
      feasible = FALSE,
      min_factor = NA_real_,
      q_grid = rep(NA_real_, length(grid_x)),
      grid_x = grid_x,
      call = match.call()
    ))
  }
  
  beta_free <- pmin(1, pmax(0, sol$solution))
  
  beta_full <- rep(NA_real_, p)
  beta_full[fixed_idx] <- fixed_beta
  beta_full[free_idx] <- beta_free
  names(beta_full) <- paste0("beta_", orders)
  
  q_grid <- 1 + as.numeric(H %*% beta_full)
  
  list(
    beta = beta_full,
    orders = orders,
    feasible = all(q_grid >= -tol),
    min_factor = min(q_grid),
    q_grid = q_grid,
    grid_x = grid_x,
    objective_value = sol$value,
    call = match.call()
  )
}


## Least squares estimation

objective_function_BellH_K <- function(params, cum) {
  K <- length(params)
  #################  left side 
  HMT <- HermMom(params, cum)
  #################  right side 
  HSM <-  HermSqu(params)
  return(sum((HMT-HSM[1:K])^2))
}


obj_fun_K <- function(params, cum) {
  K <- length(params)
  #################  left side 
  HMT <- HermMom2(params, cum)
  #################  right side 
  HSM <-  HermSqu(params)
  return(sum((HMT-HSM[1:K])^2))
}

Min_K <- function(start, cum) {
  
  start[2] <- log(start[2])  
  
  result <- stats::optim(
    par = start,
    fn = obj_fun_K,
    cum = cum,
    method = "BFGS"
  )
  
  converged <- result$convergence == 0
  
  if (!converged) {
    warning("Optimization failed to converge: ", result$message)
  }
  
  best_params <- result$par
  best_params[2] <- exp(result$par[2])
  param_names <- c("muG", "sigma2")
  if (length(best_params) > 2) {
    param_names <- c(
      param_names,
      paste0("alpha", 3:length(best_params))
    )
  }
  names(best_params) <- param_names
  
  return(list(
    converged = converged,
    objective_value = result$value,
    par = best_params
  ))
}

Min_K_Gy <- function(start = NULL, cum) {
  
  K <- length(cum)
  
  mu.X <- cum[1] 
  sig2.X <- cum[2]
  
  if (is.null(start)) {
    theta_1 <- mu.X / sqrt(sig2.X)
    theta_2 <- 0
    start <- c(theta_1, theta_2, cum[3:K])
  }
  
  lower_bounds <- c(-Inf, 0, rep(-Inf, K - 2))
  upper_bounds <- rep(Inf, K)
  
  res <- stats::optim(
    par = start,
    fn = objective_function_BellH_K,
    cum = cum[3:K],
    lower = lower_bounds,
    upper = upper_bounds,
    method = "L-BFGS-B"
  )
  
  converged <- res$convergence == 0
  
  if (!converged) {
    warning("Optimization failed to converge: ", res$message)
  }
  
  if (res$par[2] < 0) res$par[2] <- 0
  
  best_params <- res$par
  best_params[2] <- sig2.X / (1 + res$par[2])
  best_params[1] <- mu.X - sqrt(best_params[2]) * res$par[1]
  
  names(best_params) <- c("muG", "sigma2", paste0("alpha", 3:K))
  
  list(
    converged = converged,
    objective_value = res$value,
    par = best_params,
    type = "Simplified"
  )
}

Min_K_adaptive <- function(start, cum,
                           shrink = 0.9,
                           min_shrink = 0.05,
                           max_attempts = 50,
                           verbose = TRUE) {
  
  if (length(cum) < 3) {
    return(Min_K(start, cum))
  }
  
  original_cum <- cum
  factor <- 1
  
  for (attempt in seq_len(max_attempts)) {
    
    cum_try <- original_cum
    cum_try[-c(1, 2)] <- factor * original_cum[-c(1, 2)]
    
    fit <- tryCatch(
      Min_K(start, cum_try),
      warning = function(w) {
        invokeRestart("muffleWarning")
      },
      error = function(e) NULL
    )
    
    if (!is.null(fit) && isTRUE(fit$converged)) {
      fit$cum_used <- cum_try
      fit$shrink_factor <- factor
      fit$attempts <- attempt
      return(fit)
    }
    
    if (verbose) {
      message(
        "Attempt ", attempt,
        " failed; shrinking higher cumulants by factor ",
        signif(factor * shrink, 4)
      )
    }
    
    factor <- factor * shrink
    
    if (factor < min_shrink) break
  }
  
  warning("Adaptive optimization failed to converge.")
  
  return(list(
    converged = FALSE,
    objective_value = NA_real_,
    par = NA,
    cum_used = cum_try,
    shrink_factor = factor,
    attempts = attempt
  ))
}

### WRAPPER for MIN
#' Least square estimation for Extended Gram Charlier
#'
#' Provides a unified interface to different estimation procedures based on
#' cumulant matching. The function allows the user to choose between a full
#' optimization, an adaptive strategy for difficult cases, or a simplified
#' estimation method.
#'
#' @param start Optional numeric vector of starting values for the optimization.
#'   The vector should contain parameters in the form
#'   \code{c(muG, sigma2, alpha3, ..., alphaK)}.
#'   If \code{NULL}, a default starting value is internally constructed
#'   (only used for \code{type = "Simplified"}).
#'
#' @param cum Numeric vector of cumulants. Typically includes mean, variance,
#'   skewness, kurtosis, and possibly higher-order cumulants.
#'
#' @param type Character string specifying the estimation method. One of:
#'   \describe{
#'     \item{"Full"}{Standard minimization of the cumulant-based objective
#'       function using \code{Min_K}.}
#'     \item{"Adaptive"}{Adaptive procedure that progressively shrinks higher-order
#'       cumulants if optimization fails, based on \code{Min_K_adaptive}.}
#'     \item{"Simplified"}{Reduced estimation procedure using
#'       \code{Min_K_Gy}, with internally constructed starting values if not provided.}
#'   }
#'
#' @param ... Additional arguments passed to the underlying estimation function.
#'   For example, parameters controlling the adaptive procedure such as
#'   \code{shrink}, \code{max_attempts}, etc.
#'
#' @return A list with components:
#'   \describe{
#'     \item{converged}{Logical indicating whether the optimization converged.}
#'     \item{objective_value}{Value of the objective function at the optimum.}
#'     \item{par}{Estimated parameters.}
#'     \item{type}{Character string indicating the method used.}
#'     \item{...}{Additional elements returned by the underlying method
#'       (e.g., shrink factor for adaptive procedure).}
#'   }
#'
#' @examples
#' \dontrun{
#' # Full estimation
#' fit_full <- GrCh.minE(start = start, cum = Scum, type = "Full")
#'
#' # Adaptive estimation
#' fit_adapt <- GrCh.minE(start = start, cum = Scum,
#'                           type = "Adaptive", shrink = 0.9)
#'
#' # Simplified estimation (no start needed)
#' fit_simpl <- GrCh.minE(cum = Scum, type = "Simplified")
#' }
#'
#'
#' @export
GrCh.minE <- function(start = NULL, cum,
                 type = c("Full", "Adaptive", "Simplified"),
                 ...) {
  
  type <- match.arg(type)
  
  if (type %in% c("Full", "Adaptive") && is.null(start)) {
    stop("For type = 'Full' or 'Adaptive', start must be provided.")
  }
  
  fit <- switch(
    type,
    
    "Full" = Min_K(
      start = start,
      cum = cum
    ),
    
    "Adaptive" = Min_K_adaptive(
      start = start,
      cum = cum,
      ...
    ),
    
    "Simplified" = Min_K_Gy(
      start = start,
      cum = cum
    )
  )
  
  fit$type <- type
  
  fit
}


## Algorithm

#' Select the best Gram–Charlier density approximation via multiple estimation methods
#'
#' This function implements a selection procedure to identify the best-fitting
#' Gram–Charlier density approximation among several estimation strategies.
#' It combines different cumulant-based estimation methods and evaluates their
#' performance based on the distance between the model-implied cumulative
#' distribution function (CDF) and the empirical CDF.
#'
#' @param start Numeric vector of starting values for the optimization procedures.
#'   Typically of the form \code{c(muG, sigma2, alpha3, ..., alphaK)}.
#'
#' @param cum Numeric vector of cumulants, usually including mean, variance,
#'   skewness, kurtosis, and possibly higher-order cumulants.
#'
#' @param data Numeric vector of observed data used to evaluate the fitted models.
#'
#' @param distance Character string specifying the distance measure used to compare
#'   the fitted CDF with the empirical CDF. One of:
#'   \describe{
#'     \item{"L1"}{Mean absolute difference between model and empirical CDF.}
#'     \item{"KS"}{Kolmogorov distance (maximum absolute difference).}
#'   }
#'
#' @param verbose Logical; if \code{TRUE}, prints a summary table of results.
#'
#' @details
#' The procedure consists of three main steps:
#'
#' \strong{Step 1: Estimation}  
#' Four candidate parameter vectors are obtained using different estimation methods:
#' \itemize{
#'   \item \code{Min_K}: full cumulant matching via nonlinear optimization.
#'   \item \code{Min_K_adaptive}: adaptive version of \code{Min_K} with shrinking
#'         higher-order cumulants in case of convergence issues.
#'   \item \code{Min_K_Gy} with internally generated starting values.
#'   \item \code{Min_K_Gy} using the user-provided starting values.
#' }
#'
#' \strong{Step 2: Density evaluation}  
#' For each candidate parameter vector, the Gram–Charlier density is evaluated
#' at the observed data points using \code{GrCh.dens.squ}. A log-density score
#' and diagnostic quantities (e.g., minimum density value) are also computed.
#'
#' \strong{Step 3: CDF comparison}  
#' The model-implied CDF is computed using \code{GrCh.Int.squ} and compared to
#' the empirical CDF. The comparison is based on either the L1 distance or the
#' Kolmogorov distance.
#'
#' The candidate solution minimizing the selected distance is returned as the best fit.
#'
#' @return A list containing:
#'   \describe{
#'     \item{best_method}{Name of the estimation method yielding the best fit.}
#'     \item{best_param}{Estimated parameters corresponding to the best method.}
#'     \item{best_fit}{Full output of the selected estimation procedure.}
#'     \item{scores}{Data frame summarizing performance metrics for all candidates,
#'       including convergence status, CDF distance, and density diagnostics.}
#'     \item{all_fits}{List of all fitted models.}
#'     \item{distance}{Distance measure used for model selection.}
#'   }
#'
#' @seealso \code{\link{Min_K}}, \code{\link{Min_K_adaptive}},
#'   \code{\link{Min_K_Gy}}, \code{\link{GrCh.dens.squ}},
#'   \code{\link{GrCh.Int.squ}}
#'
#' @export
GrCh.fitting <- function(start, cum, data,
                         distance = c("L1", "KS"),
                         verbose = TRUE) {
  
  distance <- match.arg(distance)
  
  candidates <- list()
  
  # 1. Full Min_K
  candidates[["Full"]] <- tryCatch(
    Min_K(start = start, cum = cum),
    error = function(e) list(converged = FALSE, par = NA, error = e$message)
  )
  
  # 2. Adaptive Min_K
  candidates[["Adaptive"]] <- tryCatch(
    Min_K_adaptive(start = start, cum = cum),
    error = function(e) list(converged = FALSE, par = NA, error = e$message)
  )
  
  # 3. Simplified Min_K_Gy with default start
  candidates[["Simplified_default"]] <- tryCatch(
    Min_K_Gy(start = NULL, cum = cum),
    error = function(e) list(converged = FALSE, par = NA, error = e$message)
  )
  
  # 4. Simplified Min_K_Gy with provided start
  candidates[["Simplified_start"]] <- tryCatch(
    Min_K_Gy(start = start, cum = cum),
    error = function(e) list(converged = FALSE, par = NA, error = e$message)
  )
  
  x <- sort(data)
  F_emp <- stats::ecdf(data)(x)
  
  evaluate_candidate <- function(fit, name) {
    
    param <- fit$par
    
    if (all(is.na(param))) {
      return(data.frame(
        method = name,
        converged = FALSE,
        cdf_distance = NA_real_,
        log_density_score = NA_real_,
        min_density = NA_real_
      ))
    }
    
    dens <- tryCatch(
      GrCh.dens.squ(data, param),
      error = function(e) rep(NA_real_, length(data))
    )
    
    log_density_score <- if (all(is.finite(dens)) && all(dens > 0)) {
      sum(log(dens))
    } else {
      NA_real_
    }
    
    min_density <- suppressWarnings(min(dens, na.rm = TRUE))
    
    F_model <- tryCatch(
      sapply(x, function(z) GrCh.Int.squ(z, param)),
      error = function(e) rep(NA_real_, length(x))
    )
    
    diff <- F_model - F_emp
    
    cdf_distance <- if (all(is.finite(diff))) {
      
      if (distance == "L1") {
        mean(abs(diff))
      } else if (distance == "KS") {
        max(abs(diff))
      }
      
    } else {
      NA_real_
    }
    
    data.frame(
      method = name,
      converged = isTRUE(fit$converged),
      cdf_distance = cdf_distance,
      log_density_score = log_density_score,
      min_density = min_density
    )
  }
  
  scores <- do.call(
    rbind,
    Map(evaluate_candidate, candidates, names(candidates))
  )
  
  scores <- scores[order(scores$cdf_distance), ]
  
  best_method <- scores$method[which.min(scores$cdf_distance)]
  best_fit <- candidates[[best_method]]
  
  if (verbose) {
    print(scores)
  }
  
  list(
    best_method = best_method,
    best_param = best_fit$par,
    best_fit = best_fit,
    scores = scores,
    all_fits = candidates,
    distance = distance
  )
}
