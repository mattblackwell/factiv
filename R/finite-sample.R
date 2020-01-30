##' Estimate main effect IV ratios for 2^K factorial experiments
##'
##' This function estimates the ratio of the effect of treatment
##' assignment on the outcome to the effect of treatment assignment on
##' treatment uptake in 2^K factorial experiments. The approach uses
##' finite sample asymptotic inference to generate confidence
##' intervals. 
##'
##' 
##' @title Finite-Sample IV Estimation of 2^K Factorial Design
##' @param formula formula specification of the factorial design with
##'   noncompliance. The right-hand side of the formula should have
##'   two components separated by the \code{|} symbol, with the first
##'   component containing the K binary treatment variables and the
##'   second component containing the K binary instruments associated
##'   with each treatment variable. The order of the variables in the
##'   formula must match. 
##' @param data a data.frame on which to apply the \code{formula}.
##' @param ways the highest degree of interaction effect to estimate. 
##' @param level the confidence level required.
##' @param joint_compliers a logical indicating if all of the
##'   estimated effects should be calculated for compliers on all
##'   factors. By default, this is \code{FALSE} and the each factorial
##'   effect is estimated among the compliers for the factors in that
##'   effect.  
##' @return A list of class \code{iv_finite_factorial} that contains the following
##'   components: 
##' \item{tau}{a vector of estimated effect ratios for each factor.}
##' \item{tau_cis}{a matrix of confidence intervals for each effect
##'   ratio. This matrix has 4 columns because it is possible to have
##'   disjoint confidence intervals in this method.}
##' \item{tau_y}{a vector of the estimated effects on the outcome.}
##' \item{v_tau_y}{the estimated sample variances of the effects on
##'   the outcome.}
##' \item{tau_d}{a vector of the estimated effects on treatment
##'   uptake.} 
##' \item{v_tau_y}{the estimated sample variances of the effects on
##'   treatment uptake.}
##' \item{level}{the confidence level of \code{tau_cis}.}
##' @author Matt Blackwell
##' @references 
##'
##' Matthew Blackwell and Nicole Pashley (2020) "Noncompliance in
##'   Factorial Experiments." Working paper.
##'
##' @examples
##' data(newhaven)
##'
##' out <- iv_finite_factorial(turnout_98 ~ inperson + phone | inperson_rand
##'   + phone_rand, data = newhaven, ways = 2)
##'
##' out
##'
##' joint <- iv_finite_factorial(turnout_98 ~ inperson + phone |
##'   inperson_rand + phone_rand, data = newhaven, ways = 2,
##'   joint_compliers = TRUE)
##'
##' joint
##' 
##' @export
##' @importFrom stats model.matrix model.response

iv_finite_factorial <- function(formula, data, subset, ways = 1, level = 0.95,
                                joint_compliers = FALSE) {
  cl <- match.call(expand.dots = TRUE)
  mf <- match.call(expand.dots = FALSE)

  m <- match(
    x = c("formula", "data", "subset"),
    table = names(mf),
    nomatch = 0L
  )
  mf <- mf[c(1L, m)]
  mf[[1L]] <- as.name("model.frame")
  mf$drop.unused.levels <- TRUE

  ## must be valid formula
  formula <- Formula::as.Formula(formula)
  stopifnot(length(formula)[1] == 1L, length(formula)[2] %in% 1:2)
  if (inherits(try(terms(formula), silent = TRUE), "try-error")) {
    stop("cannot use dot '.' in formulas")
  }

  mt_d <- terms(formula, data = data, rhs = 1)
  attr(mt_d, "intercept") <- 0
  mt_z <- terms(formula, data = data, rhs = 2)
  attr(mt_z, "intercept") <- 0

  ## add to mf call
  mf$formula <- formula

  #  finally evaluate model.frame, create data matrix
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms") # terms object

  Y <- model.response(mf, "numeric")
  D <- model.matrix(mt_d, mf)
  Z <- model.matrix(mt_z, mf)
  out <- iv_finite_fit(Y, D, Z, level = level, ways = ways,
                       joint_compliers = joint_compliers)
  out$call <- cl
  out$terms <- mt
  class(out) <- "iv_finite_factorial"
  return(out)
}

iv_finite_fit <- function(y, d, z, level = 0.95, ways = 1, joint_compliers = TRUE) {
  K <- dim(z)[2]
  N <- length(y)

  dz_vals <- rep(list(c(1, 0)), 2 * K)
  dd_grid <- expand.grid(dz_vals)[, 1:K]
  zz_grid <- expand.grid(dz_vals)[, (K + 1):(2 * K)]
  zz_str_grid <- do.call(paste0, zz_grid)
  ps_grid <- expand.grid(rep(list(c("a", "n", "c")), K))
  ps_type <- 2 + zz_grid - dd_grid + 2 * dd_grid * zz_grid
  ps_dict <- list("a", c("n", "c"), "n", c("a", "c"))

  g <- expand.grid(rep(list(c(1, -1)), times = K))
  colnames(g) <- colnames(z)
  J <- nrow(g)
  n_eff <- sum(choose(K, 1:ways))
  V <- matrix(0, ncol = 3 ^ K, nrow = K)
  colnames(V) <-  do.call(paste0, ps_grid)
  all_comps <- which(rowSums(ps_grid == rep("c", K)) == K)
  for (k in 1:K) {
    if (!joint_compliers) {
      cvec <- which(ps_grid[, k] == "c")
      V[k, cvec] <- 1
    } else {
      V[k, all_comps] <- 1
    }
  }
  if (ways > 1) {
    for (k in 2:ways) {
      combs <- combn(1:K, k)
      contr_int <- matrix(NA, nrow = J, ncol = ncol(combs))
      colnames(contr_int) <- rep("", ncol(combs))
      Vk <- matrix(0, nrow = ncol(combs), ncol = 3 ^ K)
      for (j in 1:ncol(combs)) {
        this_comb <- g[, combs[, j], drop = FALSE]
        contr_int[, j] <- apply(this_comb, 1, prod)
        colnames(contr_int)[j] <- paste0(colnames(z)[combs[, j]],
                                         collapse = ":")
        if (!joint_compliers) {
          fact_comps <- rowSums(ps_grid[, combs[, j]] == rep("c", k)) == k
          Vk[j, fact_comps] <- 1
        } else {
          Vk[j, all_comps] <- 1
        }
        
      }
      g <- cbind(g, contr_int)
      V <- rbind(V, Vk)
    }
  }
  g <- g / 2 ^ (K - 1)

  A <- matrix(1, nrow = nrow(dd_grid), ncol = nrow(ps_grid))
  for (k in 1:K) {
    k_mat <- sapply(ps_dict[ps_type[,k]], function(x) 1 * (ps_grid[,k] %in% x))
    A <- A * t(k_mat)
  }
  colnames(A) <- do.call(paste0, ps_grid)
  Aw <- V %*% solve(crossprod(A)) %*% t(A)

  z_grid <- expand.grid(rep(list(c(1, 0)), times = K))
  z_grid_str <- do.call(paste0, z_grid)
  z_str <- apply(z, 1, paste0, collapse = "")

  d_grid <- expand.grid(rep(list(c(1, 0)), times = K))
  d_grid_str <- do.call(paste0, d_grid)
  d_str <- apply(d, 1, paste0, collapse = "")
  Dtilde <- matrix(0, nrow = N, ncol = J)
  for (j in 1:J) {
    Dtilde[d_str == d_grid_str[j], j] <- 1
  }
  
  s_z <- array(NA, dim = c(J + 1, J + 1, J))
  Ybar <- rep(NA, times = J)
  Dbar <- matrix(NA, nrow = J, ncol = J)
  Q_z <- array(0, dim = c(ncol(g) + nrow(V), J + 1, J))
  vcv <- matrix(0, ncol = ncol(g) + nrow(V), nrow = ncol(g) + nrow(V))
  for (j in 1:J) {
    jj <- which(z_str == z_grid_str[j])
    Ybar[j] <- mean(y[jj])
    Dbar[, j] <- colMeans(Dtilde[jj, ])
    s_z[,, j] <- var(cbind(y[jj], Dtilde[jj, ]))
    this_A <- Aw[, which(zz_str_grid == z_grid_str[j])]
    Q_z[1:ncol(g), 1, j] <- unlist(g[j, ])
    Q_z[1:ncol(g), -1, j] <- 0
    Q_z[(ncol(g) + 1):dim(Q_z)[1], 1, j] <- 0
    Q_z[(ncol(g) + 1):dim(Q_z)[1], -1, j] <- this_A
    vcv <- vcv + (1 / length(jj)) * (Q_z[,, j] %*% s_z[,, j] %*% t(Q_z[,, j]))
  }
  # Create science matrices -----------------------------

  tau_y <- t(g) %*% Ybar
  tau_d <- Aw %*% c(Dbar)

  tau <- tau_y / tau_d
  
  v_tau_y <- diag(vcv)[1:n_eff]
  v_tau_d <- diag(vcv)[n_eff + (1:n_eff)]
  c_tau_yd <- diag(vcv[1:n_eff, n_eff + (1:n_eff)])
  names(tau) <- names(tau_y) <- names(tau_d) <- colnames(g)
  names(v_tau_y) <- names(v_tau_d) <- names(c_tau_yd) <- colnames(g)
  
  alpha <- (1 - level) / 2
  qq <- qnorm(alpha)

  aa <- tau_d ^ 2 - qq ^ 2 * v_tau_d
  bb <- -2 * (tau_d * c(tau_y) - qq ^ 2 * c_tau_yd)
  cc <- c(tau_y ^ 2 - qq ^ 2 * v_tau_y)

  ## ignoring aa being exactly 0
  deter <- bb ^ 2 - 4 * aa * cc
  closed <- which(aa > 0 & deter > 0)
  infinite <- which(aa <= 0 & deter <= 0)
  disjoint <- which(aa <= 0 & deter > 0)

  moe <- rep(NA, K)
  moe[deter > 0] <- sqrt(deter[deter > 0]) / (2 * aa[deter > 0])
  cntr <- bb / (-2 * aa)
  tau_cis <- matrix(NA, nrow = n_eff, ncol = 4)
  tau_cis[closed, 1:2] <- cbind(cntr - moe, cntr + moe)[closed, ]
  tau_cis[infinite, 1] <- -Inf
  tau_cis[infinite, 2] <- Inf
  tau_cis[disjoint, 1] <- -Inf
  tau_cis[disjoint, 4] <- Inf
  tau_cis[disjoint, 2:3] <- cbind(cntr + moe, cntr - moe)[disjoint, ]
  rownames(tau_cis) <- names(tau)
  colnames(tau_cis) <- c("ci_1_lower", "ci_1_upper", "ci_2_lower", "ci_2_upper")
  return(list(tau = tau, tau_cis = tau_cis,
              tau_y = tau_y, v_tau_y = v_tau_y,
              tau_d = tau_d, v_tau_d = v_tau_d,
              level = level))
}

#' @export
print.iv_finite_factorial <- function(x, ...) {
  cat("\nCall:\n")
  print(x$call)

  cat("\nRatio of Y-Z:D-Z effects:\n")
  out <- cbind(x$tau, x$tau_cis)
  perc <- paste0(format(100 * x$level, trim = TRUE, scientific = FALSE,
                        digits = 3), "%")
  ci1 <- apply(format(x$tau_cis[, 1:2], digits = 3),
               1, paste, collapse = ", ", sep = "")
  ci1 <- paste0("(", ci1, ")")
  disj <- which(!is.na(x$tau_cis[, 3]))
  if (length(disj)) {
    ci2 <- apply(format(x$tau_cis[disj, 3:4, drop = FALSE], digits = 3),
                 1, paste, collapse = ", ", sep = "")
    ci2 <- paste0("(", ci2, ")")
    ci1[disj] <- paste(ci1[disj], ci2, sep = " U ")
  }
  out <- cbind(format(x$tau, digits = 3), ci1)
  colnames(out) <- c("Estimate", paste0(perc, " Confidence Interval"))
  rownames(out) <- names(x$tau)
  print(out, quote = FALSE)
  cat("\n")
  invisible(x)
}

calculate_rho_hat <- function(d, z) {
  K <- dim(d)[2]
  N <- dim(d)[1]
  
  dz_vals <- rep(list(c(1, 0)), 2 * K)
  ps_grid <- expand.grid(rep(list(c("a", "n", "c")), K))
  d_grid <- expand.grid(dz_vals)[, 1:K]
  z_grid <- expand.grid(dz_vals)[, (K + 1):(2 * K)]
  R <- nrow(z_grid)
  
  d_grid_str <- do.call(paste0, d_grid)
  z_grid_str <- do.call(paste0, z_grid)
  Dtilde <- matrix(0, nrow = N, ncol = R)
  Ztilde <- matrix(0, nrow = N, ncol = R)
  for (r in 1:R) {
    Dtilde[d_str == d_grid_str[r], r] <- 1
    Ztilde[z_str == z_grid_str[r], r] <- 1
  }

  d_grid <- expand.grid(rep(list(c(1,0)), times = k))
  d_grid_str <- do.call(paste0, d_grid)
  ps_type <- 2 + z_grid - d_grid + 2 * d_grid * z_grid
  ps_dict <- list("a", c("n", "c"), "n", c("a", "c"))
  A <- matrix(1, nrow = nrow(d_grid), ncol = nrow(ps_grid))
  for (k in 1:K) {
    k_mat <- sapply(ps_dict[ps_type[,k]], function(x) 1 * (ps_grid[,k] %in% x))
    A <- A * t(k_mat)
  }
  colnames(A) <- do.call(paste0, ps_grid)
  rownames(A) <- paste0(do.call(paste0, d_grid), "_", do.call(paste0, z_grid))
  rho <- rep(NA, times = nrow(ps_grid))
  names(rho) <- colnames(A)
  f_dz <- colSums(Ztilde * Dtilde) / colSums(Ztilde)
  
  Dnorm <- Ztilde * sweep(Dtilde, 2, f_dz)
  s_dz <- colSums(Dorm ^ 2) / colSums(Ztilde)
  return(list(est = f_dz, Ztilde = Ztilde, Dtilde = Dtilde, s_dz = s_dz))
}
