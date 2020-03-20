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
##' @param subset subset of the data to pass to estimation.
##' @param level the confidence level required.
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
##'   + phone_rand, data = newhaven)
##'
##' out
##'
##' joint <- iv_finite_factorial(turnout_98 ~ inperson + phone |
##'   inperson_rand + phone_rand, data = newhaven)
##'
##' joint
##' 
##' @export
##' @importFrom stats model.matrix model.response

iv_finite_factorial <- function(formula, data, subset, level = 0.95) {
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
  out <- iv_finite_fit(Y, D, Z, level = level)
  out$call <- cl
  out$terms <- mt
  class(out) <- "iv_finite_factorial"
  return(out)
}

iv_finite_fit <- function(y, d, z, level = 0.95) {
  K <- dim(z)[2]
  N <- length(y)
  ways <- K
  
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
  for (k in 1:K) {
    cvec <- which(ps_grid[, k] == "c")
    V[k, cvec] <- 1
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
        fact_comps <- rowSums(ps_grid[, combs[, j]] == rep("c", k)) == k
        Vk[j, fact_comps] <- 1        
      }
      g <- cbind(g, contr_int)
      V <- rbind(V, Vk)
    }
  }
  g_J <- g[,n_eff]
  G_J <- diag(g_J)
  
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
  R <- matrix(0, nrow = N, ncol = J)
  for (j in 1:J) {
    R[d_str == d_grid_str[j], j] <- 1
  }

  H <- y * R
  s_z <- array(NA, dim = c(2 * J, 2 * J, J))
  Hbar <- matrix(NA, nrow = J, ncol = J)
  Rbar <- matrix(NA, nrow = J, ncol = J)
  Q_z <- array(0, dim = c(3 * n_eff - 1, 2 * J, J))
  vcv <- matrix(0, ncol = 3 * n_eff - 1, nrow = 3 * n_eff - 1)
  theta <- rep(0, times = 3 * n_eff - 1)
  t_ind <- 1:n_eff
  tc_ind <- n_eff + 1:(n_eff - 1)
  d_ind <- 2 * n_eff - 1 + 1:n_eff

  for (j in 1:J) {
    jj <- which(z_str == z_grid_str[j])
    Hbar[, j] <- colMeans(H[jj, ])
    Rbar[, j] <- colMeans(R[jj, ])
    s_z[,, j] <- var(cbind(H[jj, ], R[jj, ]))
    this_A <- Aw[, which(zz_str_grid == z_grid_str[j])]
    this_tc <- t(g[,-n_eff]) %*% G_J * g_J[j]
    Q_z[t_ind, 1:J, j] <- unlist(g[j, ])
    Q_z[t_ind, -(1:J), j] <- 0
    Q_z[tc_ind, 1:J, j] <- this_tc
    Q_z[tc_ind, -(1:J), j] <- 0
    Q_z[d_ind, 1:J, j] <- 0
    Q_z[d_ind, -(1:J), j] <- this_A
    theta <- theta + Q_z[, , j] %*% c(Hbar[,j], Rbar[,j])
    vcv <- vcv + (1 / length(jj)) * (Q_z[,, j] %*% s_z[,, j] %*% t(Q_z[,, j]))
  }
  # Create science matrices -----------------------------

  
  tau <- theta[t_ind]
  tau_c <- theta[tc_ind]
  delta <- theta[d_ind] ##Aw %*% c(Dbar)

  phi <- tau / delta
  gamma <- tau_c / delta[n_eff]
  
  v_tau <- diag(vcv)[t_ind]
  v_delta <- diag(vcv)[d_ind]
  c_tau_delta <- diag(vcv[t_ind, d_ind])
  v_tau_c <- diag(vcv)[tc_ind]
  v_delta_c <- rep(v_delta[n_eff], times = length(gamma))
  c_tau_delta_c <- c(vcv[tc_ind, d_ind[n_eff]])
  
  names(phi) <- names(tau) <- names(delta) <- colnames(g)
  names(v_tau) <- names(v_delta) <- names(c_tau_delta) <- colnames(g)
  names(tau_c) <- names(gamma) <- colnames(g)[-n_eff]
  
  alpha <- (1 - level) / 2
  qq <- qnorm(alpha)

  num <- c(tau, tau_c)
  den <- c(delta, rep(delta[n_eff], length(tau_c)))
  v_num <- c(v_tau, v_tau_c)
  v_den <- c(v_delta, v_delta_c)
  c_num_den <- c(c_tau_delta, c_tau_delta_c)
  aa <- den ^ 2 - qq ^ 2 * v_den
  bb <- -2 * (den * num - qq ^ 2 * c_num_den)
  cc <- c(num ^ 2 - qq ^ 2 * v_num)

  ## ignoring aa being exactly 0
  deter <- bb ^ 2 - 4 * aa * cc
  closed <- which(aa > 0 & deter > 0)
  infinite <- which(aa <= 0 & deter <= 0)
  disjoint <- which(aa <= 0 & deter > 0)

  moe <- rep(NA, K)
  moe[deter > 0] <- sqrt(deter[deter > 0]) / (2 * aa[deter > 0])
  cntr <- bb / (-2 * aa)
  tau_cis <- matrix(NA, nrow = length(num), ncol = 4)
  tau_cis[closed, 1:2] <- cbind(cntr - moe, cntr + moe)[closed, ]
  tau_cis[infinite, 1] <- -Inf
  tau_cis[infinite, 2] <- Inf
  tau_cis[disjoint, 1] <- -Inf
  tau_cis[disjoint, 4] <- Inf
  tau_cis[disjoint, 2:3] <- cbind(cntr + moe, cntr - moe)[disjoint, ]
  rownames(tau_cis) <- c(names(phi), names(gamma))
  colnames(tau_cis) <- c("ci_1_lower", "ci_1_upper", "ci_2_lower", "ci_2_upper")
  mcafe_cis <- tau_cis[t_ind,]
  scafe_cis <- tau_cis[tc_ind,]
  mcafe_est <- cbind(tau, delta, phi)
  rownames(mcafe_est) <- names(phi)
  colnames(mcafe_est) <- c("itt_y", "itt_d", "mcafe")
  scafe_est <- cbind(tau_c, rep(delta[n_eff], length(tau_c)), gamma)
  rownames(scafe_est) <- names(gamma)
  colnames(scafe_est) <- c("itt_y", "pr_joint_c", "scafe")

  return(list(mcafe_est = mcafe_est, scafe_est = scafe_est,
              mcafe_cis = mcafe_cis, scafe_cis = scafe_cis,
              theta = theta, vcov = vcv, level = level))
}

#' @export
print.iv_finite_factorial <- function(x, ...) {
  cat("\nCall:\n")
  print(x$call)

  cat("\nMarginalized-complier factorial effects:\n")
  print(x$mcafe_est)

  cat("\nSupercomplier factorial effects:\n")
  print(x$scafe_est)
  invisible(x)
}


#' @export
summary.iv_finite_factorial <- function(x, ...) {
  cat("\nCall:\n")
  print(x$call)
  
  
  cis <- rbind(x$mcafe_cis, x$scafe_cis)
  mcafe_out <- cbind(x$mcafe_est[,"mcafe"], x$mcafe_cis)
  perc <- paste0(format(100 * x$level, trim = TRUE, scientific = FALSE,
                        digits = 3), "%")
  ci1 <- apply(format(cis[, 1:2], digits = 3),
               1, paste, collapse = ", ", sep = "")
  ci1 <- paste0("(", ci1, ")")
  disj <- which(!is.na(cis[, 3]))
  if (length(disj)) {
    ci2 <- apply(format(cis[disj, 3:4, drop = FALSE], digits = 3),
                 1, paste, collapse = ", ", sep = "")
    ci2 <- paste0("(", ci2, ")")
    ci1[disj] <- paste(ci1[disj], ci2, sep = " U ")
  }

  mcafe_ci <- ci1[1:nrow(x$mcafe_cis)]
  scafe_ci <- ci1[-(1:nrow(x$mcafe_cis))]
  cat("\nMarginalized-complier factorial effects:\n")
  mcafe_out <- cbind(format(x$mcafe_est[,3], digits = 3), mcafe_ci)
  
  colnames(mcafe_out) <- c("Estimate", paste0(perc, " Confidence Interval"))
  rownames(mcafe_out) <- rownames(x$mcafe_est)
  print(mcafe_out, quote = FALSE)

  cat("\nSupercomplier factorial effects:\n")
  scafe_out <- cbind(format(x$scafe_est[,3], digits = 3), scafe_ci)
  colnames(scafe_out) <- c("Estimate", paste0(perc, " Confidence Interval"))
  rownames(scafe_out) <- rownames(x$scafe_est)
  print(scafe_out, quote = FALSE)

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
