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
##' @param method character indiciating if the GMM should be
##' \code{"onestep"} using the identity weighting matrix (default) or
##' \code{"iterative"} to estimate the GMM weighting matrix. 
##' @param max_iter maximum number of iterations to conduct when
##'   estimating the interative GMM. 
##' @param tol criteria for stopping the iterative GMM procedure. 
##' @return A list of class \code{iv_factorial} that contains the following
##'   components: 
##' \item{rho}{vector of estimated compliance class
##'   probabilities.}
##' \item{psi}{vector of the estimated conditional mean of the outcome
##'   within the compliance classes.}
##' \item{init}{list of two vectors, with the initial starting values
##'   for the GMM estimation based on a just-identified estimation.}
##' \item{vcov}{estimated asymptotic variance matrix of the combined
##'   \code{rho} and \code{psi} parameters.}
##' \item{tau}{vector of estimated main effects of each factor among
##'   supercompliers.}
##' \item{tau_se}{vector of estimated standard errors for the
##'   estimated effects in \code{tau}.}
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

iv_factorial <- function(formula, data, subset, method = "onestep", max_iter = 15,
                         tol = 10e-7) {
  cl <- match.call(expand.dots = TRUE)
  mf <- match.call(expand.dots = FALSE)
  stopifnot(method %in% c("onestep", "iterative"))
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
  out <- iv_gmm_fit(Y, D, Z, method, max_iter = max_iter, tol = tol)
  out$call <- cl
  out$df.residual <- nrow(D) - ncol(out$vcov)
  return(out)
}



##' Fits a GMM for IV in 2^K Factorial Experiments
##'
##'
##' @title GMM for Estimation of 2^K Factorial Design
##' @param y vector of outcomes.
##' @param d K-column matrix of binary treatments.
##' @param z K-column matrix of binary instruments.
##' @param method character indiciating if the GMM should be
##' \code{"onestep"} using the identity weighting matrix (default) or
##' \code{"iterative"} to estimate the GMM weighting matrix. 
##' @param max_iter maximum number of iterations to conduct when
##'   estimating the iterative GMM. 
##' @param tol criteria for stopping the iterative GMM procedure. 
##' @return A list of class \code{iv_factorial} that contains the following
##'   components: 
##' \item{rho}{vector of estimated compliance class
##'   probabilities.}
##' \item{psi}{vector of the estimated conditional mean of the outcome
##'   within the compliance classes.}
##' \item{init}{list of two vectors, with the initial starting values
##'   for the GMM estimation based on a just-identified estimation.}
##' \item{vcov}{estimated asymptotic variance matrix of the combined
##'   \code{rho} and \code{psi} parameters.}
##' \item{tau}{vector of estimated main effects of each factor among
##'   supercompliers.}
##' \item{tau_se}{vector of estimated standard errors for the
##'   estimated effects in \code{tau}.}
##' @author Matt Blackwell
##' @importFrom stats optim

iv_gmm_fit <- function(y, d, z, method, max_iter = 15, tol = 10e-7) {  
  N <- length(y)
  K <- dim(d)[2]
  inits <- iv_init(y, d, z)
  A <- inits$A
  B <- inits$B
  dat <- cbind(d, z, y)
  W0 <- diag(nrow = length(inits$A_valid) + length(inits$B_valid))
  start <- c(inits$rho[!is.na(inits$rho)][-1],
             inits$psi[!is.na(inits$psi)])
  first <- optim(par = start,
                 fn = iv_g_loss, gr = iv_grad,
                 x = dat, W = W0, Z = inits$Ztilde, D = inits$Dtilde,
                 A = A, B = B, A_valid = inits$A_valid, B_valid = inits$B_valid,
                 method = "BFGS")
  res <- first$par
  g_out <- iv_g(res, x = dat, W = W0, Z = inits$Ztilde, D = inits$Dtilde,
               A = A, B = B, A_valid = inits$A_valid, B_valid = inits$B_valid)
  ghat <- g_out$moments
  delt <- 1000
  count <- 0
  
  if (method == "onestep" | qr(ghat)$rank < ncol(ghat)) {
    if (qr(ghat)$rank < ncol(ghat) & method == "iterative") {
      warning("iterative estimator has singular weight matrix, returning single step")
    }
    count <- max_iter + 1
    Lhat <- (1 / N) * crossprod(ghat)
    Ohat <- W0
  }
  
  while (count <= max_iter & delt > tol) {
    g_out <- iv_g(res, x = dat, W = W0, Z = inits$Ztilde, D = inits$Dtilde,
                  A = A, B = B, A_valid = inits$A_valid, B_valid = inits$B_valid)
    ghat <- g_out$moments
    Lhat <- (1 / N) * crossprod(ghat)
    qr_ghat <- qr(Lhat)
    if (qr_ghat$rank < ncol(ghat)) {
      warning("iterative estimator has singular weight matrix, returning single step")
      res <- first$par
      count <- max_iter + 1
    } else {
      Ohat <- solve(qr_ghat)
      
      second <- optim(par = start,
                      fn = iv_g_loss, gr = iv_grad,
                      x = dat, W = Ohat, Z = inits$Ztilde, D = inits$Dtilde,
                      A = A, B = B, A_valid = inits$A_valid, B_valid = inits$B_valid,
                      method = "BFGS")
      delt <- sum(abs(res - second$par) / res)
      res <- second$par
      count <- count + 1
    }
  }
  if (count == max_iter & delt > tol & method == "iterative") {
    warning("Reached maximum iterations without converging...\n")
  }
  Ghat <- (1 / N) * iv_g(res, x = dat, W = W0, Z = inits$Ztilde,
                         D = inits$Dtilde, A = A, B = B,
                         A_valid = inits$A_valid, B_valid = inits$B_valid)$Ghat
  res_var <- matrix(NA, nrow = ncol(Ghat), ncol = ncol(Ghat))
  if (qr(Ghat)$rank < ncol(Ghat)) {
    drop_ranks <- sapply(1:ncol(Ghat), function(x) qr(Ghat[, -x])$rank)
    bad_cols <- which(drop_ranks == max(drop_ranks))[-1]
    Ghat <- Ghat[, -bad_cols]
    Ahat <- crossprod(Ghat, Ohat %*% Ghat)
    Bhat <- crossprod(Ghat, (Ohat %*% Lhat %*% Ohat) %*% Ghat)    
    res_var[-bad_cols, -bad_cols] <- (solve(Ahat) %*% Bhat %*% solve(Ahat)) / N
  } else {
    Ahat <- crossprod(Ghat, Ohat %*% Ghat)
    Bhat <- crossprod(Ghat, (Ohat %*% Lhat %*% Ohat) %*% Ghat)    
    res_var <- (solve(Ahat) %*% Bhat %*% solve(Ahat)) / N
  }
  
  rho_params <- res[1:(sum(!is.na(inits$rho)) - 1)]
  
  out <- list()
  out$rho <- inits$rho
  out$rho[colnames(A)[-1]] <- rho_params
  out$rho[colnames(A)[1]] <- 1 - sum(rho_params)
  out$psi <- inits$psi
  out$psi[colnames(B)] <- res[colnames(B)]
  colnames(res_var) <- names(res)
  rownames(res_var) <- names(res)
  out$vcov <- res_var
  out$init <- inits
  
  effs <- psi_to_tau(out$psi, out$rho, K, res_var)
  out$tau <- effs$tau
  out$tau_se <- effs$tau_se
  names(out$tau) <- names(out$tau_se) <- colnames(d)
  out$tau_se[out$tau_se == 0] <- NA
  class(out) <- "iv_factorial"
  return(out)
}


iv_init <- function(y, d, z) {
  K <- dim(d)[2]
  N <- dim(d)[1]
  if (K != dim(z)[2]) stop("d/z dims do not match")
  
  dz_vals <- rep(list(c(1, 0)), 2 * K)
  ps_grid <- expand.grid(rep(list(c("a", "n", "c")), K))
  d_grid <- expand.grid(dz_vals)[, 1:K]
  z_grid <- expand.grid(dz_vals)[, (K + 1):(2 * K)]
  R <- nrow(z_grid)

  d_grid_str <- do.call(paste0, d_grid)
  z_grid_str <- do.call(paste0, z_grid)
  d_str <- apply(d, 1, paste0, collapse = "")
  z_str <- apply(z, 1, paste0, collapse = "")
  Dtilde <- matrix(0, nrow = N, ncol = R)
  Ztilde <- matrix(0, nrow = N, ncol = R)
  for (r in 1:R) {
    Dtilde[d_str == d_grid_str[r], r] <- 1
    Ztilde[z_str == z_grid_str[r], r] <- 1
  }
  
  ps_type <- 2 + z_grid - d_grid + 2 * d_grid * z_grid
  ps_dict <- list("a", c("n", "c"), "n", c("a", "c"))
  A <- matrix(1, nrow = nrow(d_grid), ncol = nrow(ps_grid))
  for (k in 1:K) {
    k_mat <- sapply(ps_dict[ps_type[,k]], function(x) 1 * (ps_grid[,k] %in% x))
    A <- A * t(k_mat)
  }
  colnames(A) <- do.call(paste0, ps_grid)
  rownames(A) <- paste0(do.call(paste0, d_grid), "_", do.call(paste0, z_grid))
  
  B <- matrix(0, nrow = nrow(d_grid), ncol = nrow(d_grid))
  rownames(B) <- rownames(A)

  hold <- list(c("n", "c"), c("a", "c"))
  psi_grid <- matrix(NA, nrow = 0, ncol = 2 * K)
  for (j in 1:nrow(unique(d_grid))) {
    this_str <- unique(d_grid_str)[j]
    this_strata <- expand.grid(hold[unlist(unique(d_grid)[j,]) + 1])
    s_names <- do.call(paste0, this_strata)
    B[d_grid_str == this_str, d_grid_str == this_str] <- A[d_grid_str == this_str, s_names]
    colnames(B)[d_grid_str == this_str] <- paste0(this_str, "_", s_names)
  }


  rho <- rep(NA, times = nrow(ps_grid))
  names(rho) <- colnames(A)
  psi <- rep(NA, times = ncol(B))
  names(psi) <- colnames(B)

  ZD <- colSums(Ztilde * Dtilde)
  f_dz <- ZD / colSums(Ztilde)
  rho_tmp <- solve(crossprod(A)) %*% crossprod(A, f_dz)
  ## sanity check - find any marginal comliance types that have 0 prob
  ## we do this because the rho_tmp sometimes "splits" the 0 prob into
  ## a small pos and neg number across two joint strata
  bad_rho <- rep(0, times = length(rho_tmp))
  for (k in 1:K) {
    marg_c <- sum(rho_tmp[ps_grid[,k] == "c"]) < (2 * .Machine$double.eps)
    bad_rho[ps_grid[,k] == "c"] <- bad_rho[ps_grid[,k] == "c"] + marg_c
    marg_a <- sum(rho_tmp[ps_grid[,k] == "a"]) < (2 * .Machine$double.eps)
    bad_rho[ps_grid[,k] == "a"] <- bad_rho[ps_grid[,k] == "a"] + marg_a
    marg_n <- sum(rho_tmp[ps_grid[,k] == "n"]) < (2 * .Machine$double.eps)
    bad_rho[ps_grid[,k] == "n"] <- bad_rho[ps_grid[,k] == "n"] + marg_n
  }
  bad_rho <- bad_rho + (rho_tmp <= 0)
  B_valid <- which(ZD >= 0)
  A_valid <- which(ZD >= 0 & f_dz <= 1)
  ##A_valid <- A_valid[duplicated(z_grid[A_valid,, drop = FALSE])]
  rho[bad_rho == 0] <- lm.fit(x = A[, bad_rho == 0],
                              y = f_dz)$coefficients
  rho_valid <- which(!is.na(rho))
  ##rho[-rho_valid] <- NA
  ##rho[rho_valid] <- rho[rho_valid] / sum(rho[rho_valid])
  ##A <- A[A_valid, rho_valid, drop = FALSE]

  S_dz <- colSums(Ztilde * y * Dtilde) / colSums(Ztilde)
  brho_names <- sapply(strsplit(colnames(B), "_"), function(x) x[2])
  brho_pos <- match(brho_names, names(rho))
  psi_valid <- which(brho_pos %in% rho_valid)
  ## Bp <- B[, psi_valid, drop = FALSE] %*% diag(rho[brho_pos[psi_valid]])
  ## colnames(Bp) <- colnames(B[, psi_valid])
  psi[psi_valid] <- lm.fit(x = B[, psi_valid], y = S_dz)$coefficients
  psi_valid <- which(!is.na(psi))
  B <- B[, psi_valid, drop = FALSE]

  ## dropping one row of A per combination of Z since they sum to 1
  ## we don't do this above because it creates inversion issues
  ## when trying to get the initial rho estimates
  A_valid <- A_valid[duplicated(z_grid[A_valid,, drop = FALSE])]
  A <- A[,rho_valid, drop = FALSE]
  return(list(A = A, B = B, Ztilde = Ztilde, Dtilde = Dtilde,
              A_valid = A_valid, B_valid = B_valid,
              rho = rho, psi = psi))
}

iv_g <- function(theta, x, W, Z, D, A, B, A_valid, B_valid) {

  N <- nrow(Z)
  K <- (dim(x)[2] - 1) / 2 ## probably shouldn't hard code this data
  dz_vals <- rep(list(c(1, 0)), 2 * K)
  d_grid <- expand.grid(dz_vals)[, 1:K]
  z_grid <- expand.grid(dz_vals)[, (K + 1):(2 * K)]

  ## structure
  ## drops <- seq(from = 2^K, to = ncol(Z), by = 2^K)
  
  rho <- c(1 - sum(theta[colnames(A)[-1]]), theta[colnames(A)[-1]])
  names(rho)[1] <- colnames(A)[1]
  psi <- theta[colnames(B)]
  y <- x[, 5]
  A <- A[A_valid,]
  B <- B[B_valid,]
  Arho <- matrix(A %*% rho, nrow = N, ncol = nrow(A), byrow = TRUE)
  brho_names <- sapply(strsplit(colnames(B), "_"), function(x) x[2])
  brho_pos <- match(brho_names, colnames(A))
  ## Bp <- B %*% diag(rho[brho_pos])
  Bpsi <- matrix(B %*% psi, nrow = N, ncol = nrow(B), byrow = TRUE)

  ## drops are redundant orthogonality conditions
  moments <- cbind(Z[, A_valid] * (D[, A_valid] - Arho),
                   Z[, B_valid] * (y * D[, B_valid] - Bpsi))
  
  
  ghats <- colMeans(moments)
  loss <- crossprod(ghats, W %*% ghats)

  Zbar <- colSums(Z)
  rho_psi <- matrix(0, nrow = nrow(B), ncol = ncol(A) - 1)
  ## for (j in 2:ncol(A)) {
  ##   Bscreen <- diag(1 * (brho_pos == j))
  ##   Bscreen[1,1] <- -1
  ##   rho_psi[,j-1] <- -Zbar[B_valid] * (B %*% Bscreen %*% psi)
  ## }
  A1 <- A[,-1] + -1 * A[,1]
  rho_grad <- cbind(-Zbar[A_valid] * A1,
                    matrix(0, nrow = nrow(A), ncol = ncol(B)))
  psi_grad <- cbind(rho_psi, -Zbar[B_valid] * B)
  Ghat <- rbind(rho_grad, psi_grad)
  Qgrad <- crossprod(Ghat, W %*% ghats)
  
  return(list(moments = moments, loss = loss, grad = Qgrad, Ghat = Ghat))
}

iv_grad <- function(...) {
  iv_g(...)$grad
}

iv_g_loss <- function(...) {
  iv_g(...)$loss
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
  tval <- object$tau / object$tau_se
  pval <- 2 * pt(abs(tval), rdf, lower.tail = FALSE)
  out <- object[c("call", "terms", "vcov")]
  out$coefficients <- cbind(object$tau, object$tau_se, tval, pval)
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
