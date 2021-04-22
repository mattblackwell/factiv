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
##'   two components separated by the `|` symbol, with the first
##'   component containing the K binary treatment variables and the
##'   second component containing the K binary instruments associated
##'   with each treatment variable. The order of the variables in the
##'   formula must match.
##' @param data a data.frame on which to apply the `formula`.
##' @param subset subset of the data to pass to estimation.
##' @param level the confidence level required.
##' @return A list of class `iv_finite_factorial` that contains the
##' following components:
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
##' \item{level}{the confidence level of `tau_cis`.}
##' @author Matt Blackwell
##' @references
##'
##' Matthew Blackwell and Nicole Pashley (2021) "Noncompliance in
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
##' @importFrom stats model.matrix model.response terms sd var qnorm
##' pt 
##' @importFrom utils combn

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

  y <- model.response(mf, "numeric")
  d <- model.matrix(mt_d, mf)
  z <- model.matrix(mt_z, mf)
  out <- iv_finite_fit(y, d, z, level = level)
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
  eff_names <- rep(NA, n_eff)
  eff_names[1:K] <- colnames(d)
  for (k in 1:K) {
    cvec <- which(ps_grid[, k] == "c")
    V[k, cvec] <- 1
  }
  count <- K + 1
  if (ways > 1) {
    for (k in 2:ways) {
      combs <- combn(1:K, k)
      contr_int <- matrix(NA, nrow = J, ncol = ncol(combs))
      colnames(contr_int) <- rep("", ncol(combs))
      Vk <- matrix(0, nrow = ncol(combs), ncol = 3 ^ K)
      for (j in seq_len(ncol(combs))) {
        this_comb <- g[, combs[, j], drop = FALSE]
        contr_int[, j] <- apply(this_comb, 1, prod)
        colnames(contr_int)[j] <- paste0(colnames(z)[combs[, j]],
                                         collapse = ":")
        eff_names[count] <- paste0(colnames(d)[combs[, j]],
                                   collapse = ":")
        count <- count + 1
        fact_comps <- rowSums(ps_grid[, combs[, j]] == rep("c", k)) == k
        Vk[j, fact_comps] <- 1
      }
      g <- cbind(g, contr_int)
      V <- rbind(V, Vk)
    }
  }
  g_J <- g[, n_eff]
  G_J <- diag(g_J)

  g <- g / 2 ^ (K - 1)
  A <- matrix(1, nrow = nrow(dd_grid), ncol = nrow(ps_grid))
  for (k in 1:K) {
    k_mat <- sapply(ps_dict[ps_type[, k]],
                    function(x) 1 * (ps_grid[, k] %in% x))
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
  t_ind <- 1:(n_eff - 1)
  tc_ind <- c(n_eff:(2 * n_eff - 1))
  d_ind <- (2 * n_eff):(3 * n_eff - 1)

  for (j in 1:J) {
    jj <- which(z_str == z_grid_str[j])
    n_j <- length(jj)
    Hbar[, j] <- colMeans(H[jj, ])
    Rbar[, j] <- colMeans(R[jj, ])
    s_z[, , j] <- var(cbind(H[jj, ], R[jj, ]))
    this_A <- Aw[, which(zz_str_grid == z_grid_str[j])]
    this_tc <- t(g) %*% G_J * g_J[j]
    Q_z[t_ind, 1:J, j] <- unlist(g[j, -n_eff])
    Q_z[t_ind, -c(1:J), j] <- 0
    Q_z[tc_ind, 1:J, j] <- this_tc
    Q_z[tc_ind, -c(1:J), j] <- 0
    Q_z[d_ind, 1:J, j] <- 0
    Q_z[d_ind, -c(1:J), j] <- this_A
    theta <- theta + Q_z[, , j] %*% c(Hbar[, j], Rbar[, j])
    vcv <- vcv + (1 / n_j) * (Q_z[, , j] %*% s_z[, , j] %*% t(Q_z[, , j]))
  }
  # Create science matrices -----------------------------

  out <- list()
  out$num_ind <- c(t_ind, tc_ind)
  out$den_ind <- c(d_ind[-n_eff], rep(d_ind[n_eff], times = length(tc_ind)))
  eff_names <- c(eff_names[-n_eff], eff_names)
  tau_cis <- fieller_cis(theta, vcv, out$num_ind, out$den_ind, eff_names, level)

  taus <- theta[out$num_ind] / theta[out$den_ind]

  out$mcafe_est <- taus[t_ind]
  out$pcafe_est <- taus[tc_ind]

  names(out$mcafe_est) <- eff_names[t_ind]
  names(out$pcafe_est) <- eff_names[tc_ind]

  out$mcafe_cis <- tau_cis[t_ind, ]
  out$pcafe_cis <- tau_cis[tc_ind, ]
  out$level <- level
  out$theta <- theta
  out$vcov <- vcv
  return(out)
}

fieller_cis <- function(theta, vcv, num_inds, den_inds, eff_names, level) {
  alpha <- (1 - level) / 2
  qq <- qnorm(alpha)
  K <- length(num_inds)

  num <- theta[num_inds]
  den <- theta[den_inds]
  v_num <- diag(vcv)[num_inds]
  v_den <- diag(vcv)[den_inds]
  c_num_den <- diag(vcv[num_inds, den_inds])

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
  rownames(tau_cis) <- eff_names
  colnames(tau_cis) <- c("ci_1_lower", "ci_1_upper", "ci_2_lower", "ci_2_upper")
  return(tau_cis)
}

#' @export
print.iv_finite_factorial <- function(x, ...) {
  cat("\nCall:\n")
  print(x$call)

  cat("\nMarginalized-complier factorial effects:\n")
  print(x$mcafe_est)

  cat("\nPerfect-complier factorial effects:\n")
  print(x$pcafe_est)
  invisible(x)
}


#' @export
summary.iv_finite_factorial <- function(object, ...) {
  cat("\nCall:\n")
  print(object$call)


  cis <- rbind(object$mcafe_cis, object$pcafe_cis)
  mcafe_out <- cbind(object$mcafe_est, object$mcafe_cis)
  perc <- paste0(format(100 * object$level, trim = TRUE, scientific = FALSE,
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

  mcafes <- seq_along(object$mcafe_est)
  mcafe_ci <- ci1[mcafes]
  pcafe_ci <- ci1[-mcafes]
  cat("\nMarginalized-complier factorial effects:\n")
  mcafe_out <- cbind(format(object$mcafe_est, digits = 3), mcafe_ci)

  colnames(mcafe_out) <- c("Estimate", paste0(perc, " Confidence Interval"))
  rownames(mcafe_out) <- names(object$mcafe_est)
  print(mcafe_out, quote = FALSE)

  cat("\nPerfect-complier factorial effects:\n")
  pcafe_out <- cbind(format(object$pcafe_est, digits = 3), pcafe_ci)
  colnames(pcafe_out) <- c("Estimate", paste0(perc, " Confidence Interval"))
  rownames(pcafe_out) <- names(object$pcafe_est)
  print(pcafe_out, quote = FALSE)

  invisible(object)
}


##' Tidy summarizes information about the components of a model.
##'
##'
##' @title Tidy an iv_finite_factorial object
##' @param x An `iv_factorial` object produced by a call to
##'   [factiv::iv_finite_factorial()]
##' @param conf.level The confidence level to use for the confidence
##'   interval. Must be strictly greater than 0 and less than 1.
##'   Defaults to 0.95, which corresponds to a 95 percent confidence
##'   interval.
##' @param ... Additional arguments. Not used. Needed to match generic
##'   signature only.
##' @return A [tibble::tibble()] with columns:
##'
##'   \item{term}{The name of the effect term.}
##'
##'   \item{estimand}{Which complier effect being
##'   estimated.}
##'
##'   \item{estimate}{The estimated value of the effect.}
##' 
##'   \item{ci_1_lower}{Lower bound for the first interval of the
##'   Fieller confidence regresion for the estimate.}
##' 
##'   \item{ci_1_upper}{Upper bound for the first interval of the
##'   Fieller confidence regresion for tshe estimate.}
##' 
##'   \item{ci_2_lower}{Lower bound for the second interval of the
##'   Fieller confidence regresion for the estimate. Only non-`NA`
##'   when the confidence region is disjoint.}
##' \item{ci_2_upper}{Upper
##'   bound for the second interval of the Fieller confidence
##'   regresion for the estimate. Only non-`NA` when the confidence
##'   region is disjoint.}
##' @author Matt Blackwell
##' @export
tidy.iv_finite_factorial <- function(x, conf.level = 0.95, ...) {

  tms <- c(names(x$mcafe_est), names(x$pcafe_est))
  estimands <- c(rep("MCAFE", length(x$mcafe_est)),
                 rep("PCAFE", length(x$pcafe_est)))
  ests <- c(x$mcafe_est, x$pcafe_est)
  cis <- fieller_cis(x$theta, x$vcov, x$num_ind, x$den_ind, tms,
                     level = conf.level)
  ret <- tibble::tibble(term = tms,
                        estimand = estimands,
                        estimate = ests)
  ret <- dplyr::bind_cols(ret, as.data.frame(cis))
  return(ret)
}


##' @importFrom generics tidy
##' @export
generics::tidy
