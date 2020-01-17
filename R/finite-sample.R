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
  contr_grid <- expand.grid(rep(list(c(1, -1)), times = K)) / 2 ^ (K - 1)
  z_grid <- expand.grid(rep(list(c(1, 0)), times = K))
  J <- nrow(z_grid)
  z_grid_str <- do.call(paste0, z_grid)
  z_str <- apply(z, 1, paste0, collapse = "")

  # Create science matrices -----------------------------
  Ztilde <- matrix(0, nrow = N, ncol = J)
  Dtilde <- matrix(0, nrow = N, ncol = K * J)
  for (j in 1:J) {
    Dcols <- seq(j, K * J, by = J)
    Ztilde[z_str == z_grid_str[j], j] <- 1
    Dtilde[z_str == z_grid_str[j], Dcols] <- d[z_str == z_grid_str[j], ]
  }
  Ytilde <- Ztilde * y

  Ybar <- colSums(Ytilde) / colSums(Ztilde)
  Dbar <- matrix(colSums(Dtilde) / colSums(Ztilde), nrow = J, ncol = K)

  tau_y <- t(contr_grid) %*% Ybar
  tau_d <- diag(t(contr_grid) %*% Dbar)

  tau <- tau_y / tau_d
  r <- N / 2 ^ K
  Ynorm <- Ztilde * sweep(Ytilde, 2, Ybar)
  s_y <- colSums(Ynorm ^ 2) / (colSums(Ztilde) - 1)
  s_d <- matrix(NA, nrow = J, ncol = K)
  c_yd <-  matrix(NA, nrow = J, ncol = K)
  for (k in 1:K) {
    cols <- 1:J + J * (k - 1)
    Dnorm <- Ztilde * sweep(Dtilde[, cols], 2, Dbar[, k])
    s_d[, k] <- colSums(Dnorm ^ 2) / (colSums(Ztilde) - 1)
    c_yd[, k] <- colSums(Ynorm * Dnorm) / (colSums(Ztilde) - 1)
  }

  v_tau_y <- t(contr_grid ^ 2 / r) %*% s_y
  v_tau_d <- diag(t(contr_grid ^ 2 / r) %*% s_d)
  c_tau_yd <- diag(t(contr_grid ^ 2 / r) %*% c_yd)
  names(tau) <- names(tau_y) <- names(tau_d) <- colnames(z)
  names(v_tau_y) <- names(v_tau_d) <- names(c_tau_yd) <- colnames(z)
  
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
  tau_cis <- matrix(NA, nrow = K, ncol = 4)
  tau_cis[closed, 1:2] <- cbind(cntr - moe, cntr + moe)[closed, ]
  tau_cis[infinite, 1] <- -Inf
  tau_cis[infinite, 2] <- Inf
  tau_cis[disjoint, 1] <- -Inf
  tau_cis[disjoint, 4] <- Inf
  tau_cis[disjoint, 2:3] <- cbind(cntr + moe, cntr - moe)[disjoint, ]
  rownames(tau_cis) <- colnames(z)
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
