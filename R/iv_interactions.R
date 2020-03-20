##' Estimates principal stratum-specific effects and interactions in a
##' 2^K factorial experiment
##'
##' This function estimates treatment effects for 2^K factorial
##' experiments in the face of noncompliance on all factors. A
##' monotonicity assumption is assumed for both treatment-instrument
##' pairs, along with treatment exclusion. See Blackwell (2017) for
##' more details on those assumptions.
##'
##'
##' The procedure uses iterative generalized method of moments (GMM)
##' to  estimate both the proportions of each compliance class (also
##' known as principal strata) and the average potential outcomes
##' within those classes. It also provides estimates of several
##' one-way, joint, and interactive treatment effects within these
##' classes. 
##'
##' Under the above assumptions, the compliance classes are the
##' product of the compliance classes for each treatment-instrument
##' pair. For instance, \code{"cc"} is the class that would comply
##' with both treatments, \code{"ca"} is the class that would comply
##' with the first treatment and always take the second treatment, and
##' \code{"cn"} is the class that would comply with the first
##' treatment and never take the second treatment. Finally, note that
##' treatment effects are only well-defined for compliance classes for
##' which there is compliance on at least one treatment. 
##' 
##' @title IV Estimation of 2^K Factorial Design
##' @param formula formula specification of the factorial design with
##'   noncompliance. The right-hand side of the formula should have
##'   two components separated by the \code{|} symbol, with the first
##'   component containing the K binary treatment variables and the
##'   second component containing the K binary instruments associated
##'   with each treatment variable. The order of the variables in the
##'   formula must match. 
##' @param data A data.frame on which to apply the \code{formula}.
##' @param subset subset of the data to pass to estimation.
##' @param method character indiciating if the estimator should be
##' \code{"lm"} using the least squares approach (default) or
##' \code{"cmd"} to estimate via efficent minimum distance estimator. 
##' @param level the confidence level required.
##' @return A list of class \code{iv_factorial} that contains the following
##'   components: 
##' \item{rho}{vector of estimated compliance class
##'   probabilities.}
##' \item{psi}{vector of the estimated conditional mean of the outcome
##'   within the compliance classes.}
##' \item{vcov}{estimated asymptotic variance matrix of the combined
##'   \code{rho} and \code{psi} parameters.}
##' \item{scafe_est}{vector of estimated main effects of each factor among
##'   supercompliers.}
##' \item{scafe_se}{vector of estimated standard errors for the
##'   estimated effects in \code{tau}.}
##' \item{scafe_cis}{a matrix of confidence intervals for the SCAFE
##' estimates.}
##' \item{level}{the confidence level of \code{scafe_cis}.}
##' @author Matt Blackwell
##' @references Matthew Blackwell (2017) Instrumental Variable Methods
##'   for Conditional Effects and Causal Interaction in Voter
##'   Mobilization Experiments, Journal of the American Statistical
##'   Association, 112:518, 590-599,
##'   \url{https://doi.org/10.1080/01621459.2016.1246363}
##'
##' Matthew Blackwell and Nicole Pashley (2020) "Noncompliance in
##'   Factorial Experiments." Working paper.
##'
##' @examples
##' data(newhaven)
##'
##' out <- iv_factorial(turnout_98 ~ inperson + phone | inperson_rand
##'   + phone_rand, data = newhaven)
##'
##' summary(out)
##' 
##' @export
##' @importFrom stats model.matrix model.response

iv_factorial <- function(formula, data, subset, method = "lm", level = 0.95) {
  cl <- match.call(expand.dots = TRUE)
  mf <- match.call(expand.dots = FALSE)
  stopifnot(method %in% c("lm", "cmd"))
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
  K <- dim(Z)[2]
  if (method == "lm") {
    out <- factiv_lm_fit(Y, D, Z)
  } else {
    out <- factiv_cmd_fit(Y, D, Z)
  }
  effs <- psi_to_tau(out$psi, out$rho, K, out$vcov)
  out$scafe_est <- effs$tau
  out$scafe_se <- effs$tau_se
  out$level <- level
  alpha <- (1 - level) / 2
  qq <- qnorm(alpha)
  out$scafe_cis <- matrix(NA, nrow = length(out$scafe_est), ncol = 2)
  out$scafe_cis[, 1] <- out$scafe_est - qq * out$scafe_se
  out$scafe_cis[, 2] <- out$scafe_est + qq * out$scafe_se
  rownames(out$scafe_cis) <- colnames(D)
  colnames(out$scafe_cis) <- c("ci_lower", "ci_upper")  
  names(out$scafe_est) <- names(out$scafe_se) <- colnames(D)
  out$scafe_se[out$scafe_se == 0] <- NA
  class(out) <- "iv_factorial"
  out$call <- cl
  out$df.residual <- nrow(D) - ncol(out$vcov)
  return(out)
}




factiv_lm_fit <- function(y, d, z) {
  K <- dim(d)[2]
  N <- dim(d)[1]
  J <- 2 ^ K
  if (K != dim(z)[2]) stop("d/z dims do not match")
  
  dz_vals <- rep(list(c(1, 0)), 2 * K)
  ps_grid <- expand.grid(rep(list(c("a", "n", "c")), K))
  dz_d_grid <- expand.grid(dz_vals)[, 1:K]
  dz_z_grid <- expand.grid(dz_vals)[, (K + 1):(2 * K)]
  dz_z_grid_str <- do.call(paste0, dz_z_grid)
  dz_d_grid_str <- do.call(paste0, dz_d_grid)
  ps_type <- 2 + dz_z_grid - dz_d_grid + 2 * dz_d_grid * dz_z_grid
  ps_dict <- list("a", c("n", "c"), "n", c("a", "c"))
  
  A <- matrix(1, nrow = nrow(dz_d_grid), ncol = nrow(ps_grid))
  for (k in 1:K) {
    k_mat <- sapply(ps_dict[ps_type[,k]], function(x) 1 * (ps_grid[,k] %in% x))
    A <- A * t(k_mat)
  }
  colnames(A) <- do.call(paste0, ps_grid)
  rownames(A) <- paste0(dz_d_grid_str, "_", dz_z_grid_str)
  Aw <- solve(crossprod(A)) %*% t(A)
  
  B <- matrix(0, nrow = nrow(dz_d_grid), ncol = nrow(dz_d_grid))
  rownames(B) <- rownames(A)
  hold <- list(c("n", "c"), c("a", "c"))
  for (j in 1:nrow(unique(dz_d_grid))) {
    this_str <- unique(dz_d_grid_str)[j]
    this_strata <- expand.grid(hold[unlist(unique(dz_d_grid)[j, ]) + 1])
    s_names <- do.call(paste0, this_strata)
    grab_rows <- dz_d_grid_str == this_str
    B[grab_rows, grab_rows] <- A[grab_rows, s_names]
    colnames(B)[grab_rows] <- paste0(this_str, "_", s_names)
  }
  Bw <- solve(crossprod(B)) %*% t(B)

  z_grid <- expand.grid(rep(list(c(1, 0)), times = K))
  z_grid_str <- do.call(paste0, z_grid)
  z_str <- apply(z, 1, paste0, collapse = "")

  d_grid <- expand.grid(rep(list(c(1, 0)), times = K))
  d_grid_str <- do.call(paste0, d_grid)
  d_str <- apply(d, 1, paste0, collapse = "")

  ## calculate H and R data
  R <- matrix(0, nrow = N, ncol = J)
  for (j in 1:J) {
    R[d_str == d_grid_str[j], j] <- 1
  }
  H <- y * R
  s_z <- array(NA, dim = c(2 * J, 2 * J, J))
  Hbar <- matrix(NA, nrow = J, ncol = J)
  Rbar <- matrix(NA, nrow = J, ncol = J)
  tot_p <- ncol(A) + ncol(B)
  Q_z <- array(0, dim = c(tot_p, 2 * J, J))
  r_ind <- 1:ncol(A)
  h_ind <- (ncol(A) + 1) : tot_p
  vcv <- matrix(0, ncol = tot_p, nrow = tot_p)
  theta <- rep(0, times = tot_p)
  for (j in 1:J) {
    this_z <- z_grid_str[j]
    jj <- which(z_str == this_z)
    Hbar[, j] <- colMeans(H[jj, ])
    Rbar[, j] <- colMeans(R[jj, ])
    jjj_A <- grep(paste0("_", this_z), colnames(Aw))
    jjj_B <- grep(paste0("_", this_z), colnames(Bw))
    Q_z[r_ind, 1:J, j] <- Aw[, jjj_A]
    Q_z[h_ind, (J + 1):(2 * J), j] <- Bw[, jjj_B]
    s_z[,, j] <- var(cbind(R[jj, ], H[jj, ]))
    theta <- theta + Q_z[, , j] %*% c(Rbar[, j], Hbar[, j])
    vcv <- vcv + (1 / length(jj)) * (Q_z[,, j] %*% s_z[,, j] %*% t(Q_z[,, j]))
  }
  rho <- theta[r_ind]
  names(rho) <- colnames(A)
  psi <- theta[h_ind]
  names(psi) <- colnames(B)
  rownames(vcv) <- c(colnames(A), colnames(B))
  colnames(vcv) <- c(colnames(A), colnames(B))
  
  return(list(A = A, B = B, Hbar = Hbar, Rbar = Rbar, rho = rho,
              psi = psi, vcov = vcv))
}


factiv_cmd_fit <- function(y, d, z) {
  K <- dim(d)[2]
  N <- dim(d)[1]
  J <- 2 ^ K
  if (K != dim(z)[2]) stop("d/z dims do not match")
  
  dz_vals <- rep(list(c(1, 0)), 2 * K)
  ps_grid <- expand.grid(rep(list(c("a", "n", "c")), K))
  dz_d_grid <- expand.grid(dz_vals)[, 1:K]
  dz_z_grid <- expand.grid(dz_vals)[, (K + 1):(2 * K)]
  dz_z_grid_str <- do.call(paste0, dz_z_grid)
  dz_d_grid_str <- do.call(paste0, dz_d_grid)
  ps_type <- 2 + dz_z_grid - dz_d_grid + 2 * dz_d_grid * dz_z_grid
  ps_dict <- list("a", c("n", "c"), "n", c("a", "c"))
  
  A <- matrix(1, nrow = nrow(dz_d_grid), ncol = nrow(ps_grid))
  for (k in 1:K) {
    k_mat <- sapply(ps_dict[ps_type[,k]], function(x) 1 * (ps_grid[,k] %in% x))
    A <- A * t(k_mat)
  }
  colnames(A) <- do.call(paste0, ps_grid)
  rownames(A) <- paste0(dz_d_grid_str, "_", dz_z_grid_str)
  drop_mom <- grep(paste0(c(rep(0, times = K), "_"), collapse = ""),
                    rownames(A))
  drop_par <- grep(paste0(rep("n", times = K), collapse = ""), colnames(A))
  AA <- A[-drop_mom, -drop_par]
  Aw <- solve(crossprod(AA)) %*% t(AA)
  
  B <- matrix(0, nrow = nrow(dz_d_grid), ncol = nrow(dz_d_grid))
  rownames(B) <- rownames(A)
  hold <- list(c("n", "c"), c("a", "c"))
  for (j in 1:nrow(unique(dz_d_grid))) {
    this_str <- unique(dz_d_grid_str)[j]
    this_strata <- expand.grid(hold[unlist(unique(dz_d_grid)[j, ]) + 1])
    s_names <- do.call(paste0, this_strata)
    grab_rows <- dz_d_grid_str == this_str
    B[grab_rows, grab_rows] <- A[grab_rows, s_names]
    colnames(B)[grab_rows] <- paste0(this_str, "_", s_names)
  }
  Bw <- solve(crossprod(B)) %*% t(B)

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
  R <- R[, -J]
  s_z <- array(NA, dim = c(2 * J - 1, 2 * J - 1, J))
  Hbar <- matrix(NA, nrow = J, ncol = J)
  Rbar <- matrix(NA, nrow = J - 1, ncol = J)
  tot_p <- ncol(AA) + ncol(B)
  Q_z <- array(0, dim = c(tot_p, 2 * J - 1, J))
  r_ind <- 1:ncol(AA)
  h_ind <- (ncol(AA) + 1) : tot_p
  vcv <- matrix(0, ncol = tot_p, nrow = tot_p)
  theta <- rep(0, times = tot_p)
  HR_var <- matrix(0, nrow = ncol(R) * J + ncol(H) * J,
                   ncol = ncol(R) * J + ncol(H) * J)
  for (j in 1:J) {
    this_z <- z_grid_str[j]
    jj <- which(z_str == this_z)
    Hbar[, j] <- colMeans(H[jj, ])
    Rbar[, j] <- colMeans(R[jj, ])
    jjj_A <- grep(paste0("_", this_z), colnames(Aw))
    jjj_B <- grep(paste0("_", this_z), colnames(Bw))
    Q_z[r_ind, 1:(J - 1), j] <- Aw[, jjj_A]
    Q_z[h_ind, J:(2 * J - 1), j] <- Bw[, jjj_B]
    s_z[,,j] <- var(cbind(R[jj, ], H[jj, ]))
    HR_var[c(jjj_A, jjj_B + ncol(Aw)), c(jjj_A, jjj_B + ncol(Aw))] <- s_z[,, j]
    theta <- theta + Q_z[, , j] %*% c(Rbar[, j], Hbar[, j])
    vcv <- vcv + (1 / length(jj)) * (Q_z[,, j] %*% s_z[,, j] %*% t(Q_z[,, j]))
  }
  rho <- theta[r_ind]
  names(rho) <- colnames(AA)
  psi <- theta[h_ind]
  names(psi) <- colnames(B)

  if (qr(HR_var)$rank == ncol(HR_var)) {
      Ar <- nrow(AA)
  Aw_opt <- solve(crossprod(AA, solve(HR_var[1:Ar, 1:Ar]) %*% AA)) %*%
    t(AA) %*% solve(HR_var[1:Ar,1:Ar])
  Bw_opt <- solve(crossprod(B, solve(HR_var[-(1:Ar),-(1:Ar)]) %*% B)) %*%
      t(B) %*% solve(HR_var[-(1:Ar),-(1:Ar)]) 
  vcv2 <- matrix(0, ncol = tot_p, nrow = tot_p)
  theta2 <- rep(0, times = tot_p)
  for (j in 1:J) {
    this_z <- z_grid_str[j]
    jj <- which(z_str == this_z)
    jjj_A <- grep(paste0("_", this_z), colnames(Aw))
    jjj_B <- grep(paste0("_", this_z), colnames(Bw))
    Q_z[r_ind, 1:(J - 1), j] <- Aw_opt[, jjj_A]
    Q_z[h_ind, J:(2 * J - 1), j] <- Bw_opt[, jjj_B]
    theta2 <- theta2 + Q_z[, , j] %*% c(Rbar[, j], Hbar[, j])
    vcv2 <- vcv2 + (1 / length(jj)) * (Q_z[,, j] %*% s_z[,, j] %*% t(Q_z[,, j]))
  }
  rho2 <- theta2[r_ind]
  names(rho2) <- colnames(AA)
  psi2 <- theta2[h_ind]
  names(psi2) <- colnames(B)
  } else {
    warning("singular weight matrix with cmd, using lm...")
    rho2 <- rho
    psi2 <- psi
    vcv2 <- vcv
  }
  rownames(vcv2) <- c(colnames(AA), colnames(B))
  colnames(vcv2) <- c(colnames(AA), colnames(B))

  return(list(A = A, B = B, Hbar = Hbar, Rbar = Rbar, rho = rho2,
              psi = psi2, vcov = vcv2))
}

psi_to_tau <- function(psi, rho, K, vcv) {
  cons <- expand.grid(rep(list(c(1, -1)), times = K)) / 2 ^ (K - 1)
  rho_c <- paste0(rep("c", K), collapse = "")
  psi_c <- grep(rho_c, names(psi), value = TRUE)
  num <- c(t(cons) %*% psi[psi_c])
  den <- rho[rho_c]
  mains <- num / den
  if (sum(is.na(mains)) == 0) {
    num_var <- diag(t(cons) %*% vcv[psi_c, psi_c] %*% as.matrix(cons))
    den_var <- vcv[rho_c, rho_c]
    nd_cov <- c(t(cons) %*% vcv[psi_c, rho_c])
    mains_var <- den ^ (-2) * num_var +
      (num ^ 2 / den ^ 4) * den_var -
      2 * (num / den ^ 3) * nd_cov
    mains_se <- sqrt(mains_var)
  } else {
    warning("supercomplier effects not identified")
    mains_se <- rep(NA, times = K)
  }
  return(list(tau = mains, tau_se = mains_se))
}

#' @export
print.iv_factorial <- function(x, ...) {
  cat("\nCall:\n")
  print(x$call)

  cat("\nMain effects:\n")
  print(x$tau)
  cat("\n")
  invisible(x)
}

#' @export
summary.iv_factorial <- function(object, ...) {
  rdf <- object$df.residual
  tval <- object$scafe_est / object$scafe_se
  pval <- 2 * pt(abs(tval), rdf, lower.tail = FALSE)
  out <- object[c("call", "terms", "vcov")]
  out$coefficients <- cbind(object$scafe_est, object$scafe_se, tval, pval)
  out$c_prob <- object$rho[length(object$rho)]
  c_pos <- sum(!is.na(object$rho)) - 1
  out$c_prob_se <- sqrt(out$vcov[c_pos, c_pos])
  class(out) <- "summary.iv_factorial"
  out
}

#' @export
print.summary.iv_factorial <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\n Main effects among supercompliers:\n")
  stats::printCoefmat(x$coefficients, digits = digits)
  cat("\nEstimated prob. of supercompliers: ",
      formatC(x$c_prob, digits), "\tSE = ", formatC(x$c_prob_se, digits))
  cat("\n")
  invisible(x)
}
