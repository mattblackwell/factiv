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
##'   two components separated by the `|` symbol, with the first
##'   component containing the K binary treatment variables and the
##'   second component containing the K binary instruments associated
##'   with each treatment variable. The order of the variables in the
##'   formula must match.
##' @param data A data.frame on which to apply the `formula`.
##' @param subset subset of the data to pass to estimation.
##' @param method character indiciating if the estimator should be
##' `"lm"` using the least squares approach (default) or
##' `"cmd"` to estimate via efficent minimum distance estimator.
##' @param level the confidence level required.
##' @return A list of class `iv_factorial` that contains the following
##'   components:
##' \item{rho}{vector of estimated compliance class
##'   probabilities.}
##' \item{psi}{vector of the estimated conditional mean of the outcome
##'   within the compliance classes.}
##' \item{vcov}{estimated asymptotic variance matrix of the combined
##'   `rho`  and `psi` parameters.}
##' \item{pcafe_est}{vector of estimated main effects of each factor among
##'   perfect compliers.}
##' \item{pcafe_se}{vector of estimated standard errors for the
##'   estimated effects in `tau`.}
##' \item{pcafe_cis}{a matrix of confidence intervals for the PCAFE
##' estimates.}
##' \item{level}{the confidence level of the returned confience
##' intervals.}
##' @author Matt Blackwell
##' @references Matthew Blackwell (2017) Instrumental Variable Methods
##'   for Conditional Effects and Causal Interaction in Voter
##'   Mobilization Experiments, Journal of the American Statistical
##'   Association, 112:518, 590-599,
##'   \doi{10.1080/01621459.2016.1246363}
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
  effs <- psi_to_tau(out$psi, out$rho, K, out$vcov, colnames(D))
  out$pcafe_est <- effs$pcafe_est
  out$pcafe_se <- effs$pcafe_se
  out$mcafe_est <- effs$mcafe_est
  out$mcafe_se <- effs$mcafe_se
  out$level <- level
  alpha <- (1 - level) / 2
  qq <- abs(qnorm(alpha))
  out$pcafe_cis <- matrix(NA, nrow = length(out$pcafe_est), ncol = 2)
  out$pcafe_cis[, 1] <- out$pcafe_est - qq * out$pcafe_se
  out$pcafe_cis[, 2] <- out$pcafe_est + qq * out$pcafe_se
  out$mcafe_cis <- matrix(NA, nrow = length(out$mcafe_est), ncol = 2)
  out$mcafe_cis[, 1] <- out$mcafe_est - qq * out$mcafe_se
  out$mcafe_cis[, 2] <- out$mcafe_est + qq * out$mcafe_se
  rownames(out$pcafe_cis) <- names(out$pcafe_est)
  colnames(out$pcafe_cis) <- c("ci_lower", "ci_upper")
  rownames(out$mcafe_cis) <- names(out$mcafe_est)
  colnames(out$mcafe_cis) <- c("ci_lower", "ci_upper")

  out$pcafe_se[out$pcafe_se == 0] <- NA
  out$mcafe_se[out$mcafe_se == 0] <- NA
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

psi_to_tau <- function(psi, rho, K, vcv, var_names) {
  J <- 2 ^ K - 1
  L <- 2 ^ K

  ## reference grids
  ps_grid <- expand.grid(rep(list(c("a", "n", "c")), K))
  z_grid <- expand.grid(rep(list(c(1, 0)), times = K))
  z_grid_str <- do.call(paste0, z_grid)

  ## creating contrast matrices
  g <- expand.grid(rep(list(c(1, -1)), times = K))
  rownames(g) <- z_grid_str
  g_m_psi <- matrix(0, nrow = length(rho) + length(psi), ncol = J)
  rownames(g_m_psi) <- rownames(vcv)
  g_s_psi <- g_m_rho <- g_s_rho <- g_m_psi
  psi_d_grid <- sapply(strsplit(rownames(vcv), "_"),
                        function(x) x[[1]])
  psi_ps_grid <- sapply(strsplit(names(psi), "_"),
                        function(x) x[[2]])
  num_c <- rowSums(ps_grid == "c")
  names(num_c) <- do.call(paste0, ps_grid)
  psi_adj <- c(rep(0, times = length(rho)),
               num_c[psi_ps_grid])

  g_s_rho[paste0(rep("c", K), collapse = ""),] <- 1
  scomps_psi <- grep(paste0(c("_", rep("c", K)), collapse = ""),
                     rownames(vcv))
  comb_list <- all_subsets(1:K)
  eff_labs <- character(J)
  for (j in 1:J) {
    this_comb <- comb_list[[j]]
    k <- length(this_comb)
    eff_labs[j] <- paste0(var_names[this_comb], collapse = ":")

    this_contr <- g[, this_comb, drop = FALSE]
    this_contr <- apply(this_contr, 1, prod)
    mcomps <- rowSums(ps_grid[, this_comb, drop = FALSE] == rep("c", k)) == k
    mcomps <- apply(ps_grid[mcomps,, drop = FALSE], 1, paste0, collapse = "")
    mcomps_psi <- grep(paste0("_(", paste0(mcomps, collapse = "|"), ")"),
                       rownames(vcv))
    g_m_rho[mcomps, j] <- 1
    g_m_psi[mcomps_psi, j] <- this_contr[psi_d_grid[mcomps_psi]]
    g_s_psi[scomps_psi, j] <- this_contr[psi_d_grid[scomps_psi]]
  }
  colnames(g_m_psi) <- colnames(g_s_psi) <- eff_labs
  colnames(g_m_rho) <- colnames(g_s_rho) <- eff_labs
  g_m_psi <- g_m_psi / 2 ^ (psi_adj - 1)
  g_s_psi <- g_s_psi / 2 ^ (K - 1)

  theta <- c(rho, psi)
  m_num <- c(t(g_m_psi) %*% theta)
  m_den <- c(t(g_m_rho) %*% theta)
  s_num <- c(t(g_s_psi) %*% theta)
  s_den <- c(t(g_s_rho) %*% theta)

  mcafe_est <- m_num / m_den
  pcafe_est <- s_num / s_den
  m_num_var <- diag(t(g_m_psi) %*% vcv %*% g_m_psi)
  s_num_var <- diag(t(g_s_psi) %*% vcv %*% g_s_psi)
  m_den_var <- diag(t(g_m_rho) %*% vcv %*% g_m_rho)
  s_den_var <- diag(t(g_s_rho) %*% vcv %*% g_s_rho)
  m_cov <- diag(t(g_m_psi) %*% vcv %*% g_m_rho)
  s_cov <- diag(t(g_s_psi) %*% vcv %*% g_s_rho)
  mcafe_var <- m_den ^ (-2) * m_num_var +
    (m_num ^ 2 / m_den ^ 4) * m_den_var -
    2 * (m_num / m_den ^ 3) * m_cov
  pcafe_var <- s_den ^ (-2) * s_num_var +
    (s_num ^ 2 / s_den ^ 4) * s_den_var -
    2 * (s_num / s_den ^ 3) * s_cov
  mcafe_se <- sqrt(mcafe_var)
  pcafe_se <- sqrt(pcafe_var)

  names(pcafe_est) <- names(pcafe_se) <- eff_labs
  names(mcafe_est) <- names(mcafe_se) <- eff_labs

  ## last effect is really an pcafe
  mcafe_est <- mcafe_est[-J]
  mcafe_se <- mcafe_se[-J]
  return(list(mcafe_est = mcafe_est, mcafe_se = mcafe_se,
              pcafe_est = pcafe_est, pcafe_se = pcafe_se))
}

all_subsets <- function(x) {
  unlist(lapply(x, function(z) as.list(as.data.frame(combn(x, z)))),
         recursive = FALSE)
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
  tval <- object$pcafe_est / object$pcafe_se
  pval <- 2 * pt(abs(tval), rdf, lower.tail = FALSE)
  out <- object[c("call", "terms", "vcov")]
  out$coefficients <- cbind(object$pcafe_est, object$pcafe_se, tval, pval)
  out$c_prob <- object$rho[length(object$rho)]
  c_pos <- sum(!is.na(object$rho)) - 1
  out$c_prob_se <- sqrt(out$vcov[c_pos, c_pos])
  class(out) <- "summary.iv_factorial"
  out
}

#' @export
print.summary.iv_factorial <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\n Main effects among perfect compliers:\n")
  stats::printCoefmat(x$coefficients, digits = digits)
  cat("\nEstimated prob. of perfect compliers: ",
      formatC(x$c_prob, digits), "\tSE = ", formatC(x$c_prob_se, digits))
  cat("\n")
  invisible(x)
}


##' Tidy summarizes information about the components of a model.
##'
##'
##' @title Tidy an iv_factorial object
##' @param x An `iv_factorial` object produced by a call to [factiv::iv_factorial()]
##' @param conf.int Logical indicating whether or not to include a
##' confidence interval in the tidied output. Defaults to `FALSE`
##' @param conf.level The confidence level to use for the confidence
##' interval if `conf.int = TRUE`. Must be strictly greater than 0
##' and less than 1. Defaults to 0.95, which corresponds to a 95
##' percent confidence interval.
##' @param ... Additional arguments. Not used. Needed to match generic
##' signature only.
##' @return A [tibble::tibble()] with columns:
##' \item{term}{The name of the effect term.}
##' \item{estimand}{Which complier effect being estimated.}
##' \item{estimate}{The estimated value of the effect.}
##' \item{std.error}{The estimated standard error of the effect.}
##' \item{conf.low}{Lower bound of the confidence interval for the
##' estimate.}
##' \item{conf.high}{Upper bound of the confidence interval for the
##' estimate.} 
##' @author Matt Blackwell
##' @export
tidy.iv_factorial <- function(x, conf.int = FALSE, conf.level = 0.95, ...) {

  tms <- c(names(x$mcafe_est), names(x$pcafe_est))
  estimands <- c(rep("MCAFE", length(x$mcafe_est)),
                 rep("PCAFE", length(x$pcafe_est)))
  ests <- c(x$mcafe_est, x$pcafe_est)
  ses <- c(x$mcafe_se, x$pcafe_se)

  ret <- tibble::tibble(term = tms,
                        estimand = estimands,
                        estimate = ests,
                        std.error = ses)
  if (conf.int) {
    alpha <- (1 - conf.level) / 2
    qq <- abs(qnorm(alpha))
    ret <- dplyr::mutate(ret,
                         conf.low = estimate - qq * std.error,
                         conf.high = estimate + qq * std.error)
  }
  return(ret)
}

#' @importFrom utils globalVariables
globalVariables(
  c(
    "estimate",
    "std.error"
  )
)
