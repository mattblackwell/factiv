##' Complier covariate profiles
##'
##' Calculates averages of covariates by compliance group in a
##' 2^K factorial setting.
##'
##'  
##' @param formula one-sided formula to indicate the treatment
##'   assignment, treatment uptake, and covariates. The right-hand
##'   side of the formula should have three components separated by
##'   the `|` symbol, with the first component containing the K binary
##'   treatment variables (treatment assignment), the second component
##'   containing the K binary instruments associated with each
##'   treatment variable (treatment uptake), and the third giving the
##'   covariates to be included in the profile. The order of the
##'   variables in the first two parts of the formula must match.
##' @param data a data.frame on which to apply the `formula`.
##' @param subset subset of the data to pass to estimation.
##' @return A list with two objects:
##' \item{raw_table}{a data.frame whose rows represent the covariates and
##' whose columns represent the different compliance groups. Each
##' entry is the estimated mean of the covariate for that compliance
##' group. }
##' \item{std_table}{a data.frame similarly structured to raw_table
##' but with the standardized difference between the compilance group
##' means and the overall means in place of the raw means.}
##'
##' @references
##'
##' Matthew Blackwell and Nicole Pashley (2021) "Noncompliance in
##'   Factorial Experiments." Working paper.
##'
##' @author Matthew Blackwell
##' @examples
##'
##' data(newhaven)
##'
##' cov_prof <- compliance_profile(~ inperson + phone | inperson_rand
##'   + phone_rand | age + maj_party + turnout_96, data = newhaven)
##'
##' cov_prof
##' @export
compliance_profile <- function(formula, data, subset) {
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
  if (inherits(try(terms(formula), silent = TRUE), "try-error")) {
    stop("cannot use dot '.' in formulas")
  }

  mt_d <- terms(formula, data = data, rhs = 1)
  attr(mt_d, "intercept") <- 0
  mt_z <- terms(formula, data = data, rhs = 2)
  attr(mt_z, "intercept") <- 0
  mt_x <- terms(formula, data = data, rhs = 3)
  attr(mt_x, "intercept") <- 0

  ## add to mf call
  mf$formula <- formula

  #  finally evaluate model.frame, create data matrix
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms") # terms object

  d <- model.matrix(mt_d, mf)
  z <- model.matrix(mt_z, mf)
  x <- model.matrix(mt_x, mf)
  out <- comp_cov_fit(x, d, z)
  sds <- apply(x, 2, sd, na.rm = TRUE)
  out_std <- out
  out_std[, -c(1,2)] <- (out[, -c(1, 2)] - out[, 2]) / sds
  return(list(raw_table = out, std_table = out_std))
}

comp_cov_fit <- function(x, d, z) {
  K <- dim(z)[2]
  N <- nrow(x)
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
  delta <- rep(0, times = n_eff)
    
  for (j in 1:J) {
    jj <- which(z_str == z_grid_str[j])
    Rbar <- colMeans(R[jj, ])
    this_A <- Aw[, which(zz_str_grid == z_grid_str[j])]
    delta <- delta + this_A %*% Rbar
  }

  out <- matrix(NA, ncol = n_eff + 2, nrow = ncol(x))
  out <- as.data.frame(out)
  out[, 1] <- colnames(x)


  out[, 2] <- colMeans(x, na.rm = TRUE)
  for (s in seq_len(ncol(x))) {
    R <- matrix(0, nrow = N, ncol = J)
    for (j in 1:J) {
      R[d_str == d_grid_str[j], j] <- x[d_str == d_grid_str[j], s]
    }
    theta <- rep(0, times = n_eff)
    
    for (j in 1:J) {
      jj <- which(z_str == z_grid_str[j])
      Rbar <- colMeans(R[jj, ])
      this_A <- Aw[, which(zz_str_grid == z_grid_str[j])]
      theta <- theta + this_A %*% Rbar
    }
    out[s, -c(1, 2)] <- theta / delta
  }
  out <- as.data.frame(out)
  names(out) <- c("term", "overall", eff_names)
  return(out)
}
