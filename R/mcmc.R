#' Metropolis-Hastings sampler for a single weight class
#'
#' Draws from the posterior distribution of the bias parameter and its
#' hyperparameters for a set of two-urn matching observations that share a
#' common weight. This is the hierarchical extension of the inferential
#' problem of Puza and Bonfrer (2018), and is most useful when several
#' matching experiments are pooled.
#'
#' @param B Integer; number of burn-in iterations.
#' @param J Integer; number of sampling iterations retained after burn-in.
#' @param N1vec,N2vec Numeric vectors; `N1` and `N2` for each observation.
#' @param N0vec Numeric vector; number of common items `N` for each
#'   observation.
#' @param M1vec,M2vec Numeric vectors; sample sizes `m1` and `m2` for each
#'   observation.
#' @param mvec Numeric vector; observed numbers of matches.
#' @param Cvec Integer vector; weight-class index for each observation.
#' @param eta,tau Numeric; shape and rate hyperparameters of the gamma prior
#'   on the precision `lambda`.
#' @param mu0,sig0 Numeric; mean and scale hyperparameters of the normal prior
#'   on `mu`.
#' @param muinit,rinit Numeric; starting values for `mu` and the latent
#'   weights.
#' @param del Numeric; random-walk proposal standard deviation(s).
#' @param u,v Numeric; parameters of the internal link function mapping the
#'   latent scale to the positive weight `w`.
#' @param hf Integer; if `0`, use the log-scale form of the log full
#'   conditional, otherwise the link-function form (the default).
#'
#' @return A list with components `muv` (sampled `mu`), `lamv` (sampled
#'   `lambda`), `rm` (matrix of sampled latent weights, one column per class)
#'   and `rar` (acceptance rate per class).
#'
#' @references
#' Puza, B. and Bonfrer, A. (2018). A series of two-urn biased sampling
#' problems. *Communications in Statistics - Theory and Methods*, 47(1),
#' 80-91. \doi{10.1080/03610926.2017.1300282}
#'
#' @seealso [restable()] and [estprop()] to summarise and diagnose the output.
#'
#' @examples
#' # Short illustrative run; use many more iterations for real inference.
#' set.seed(4172)
#' res <- MHALG(
#'   B = 5, J = 20,
#'   N1vec = 37, N2vec = 45, N0vec = 16,
#'   M1vec = 12, M2vec = 8, mvec = 3,
#'   Cvec = 1, del = 0.1
#' )
#' str(res)
#'
#' @importFrom stats rgamma rnorm runif density quantile
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
MHALG <- function(B = 5, J = 20,
                  N1vec = 150, N2vec = 200, N0vec = 50,
                  M1vec = 40, M2vec = 60, mvec = 10,
                  Cvec = 1, eta = 0.001, tau = 1000, mu0 = 0,
                  sig0 = 0.001, muinit = 0, rinit = 0, del = 0.1,
                  u = exp(20), v = 1, hf = 1) {
  n <- length(mvec)
  K <- max(Cvec)
  sig02 <- sig0^2
  mu <- muinit
  lam <- NA
  r <- rep(rinit, K)
  rct <- rep(0, K)
  rm <- matrix(0, nrow = 1 + B + J, ncol = K)
  rm[1, ] <- r
  muv <- lamv <- rep(0, 1 + B + J)
  muv[1] <- mu
  lamv[1] <- NA
  pb <- txtProgressBar(min = 0, max = B + J, style = 3)
  for (j in seq_len(B + J)) {
    setTxtProgressBar(pb, j)
    lam <- rgamma(1, eta + 0.5 * K, tau + 0.5 * sum((r - mu)^2))
    c0 <- K / (K + 1 / (lam * sig02))
    mustar <- (1 - c0) * mu0 + c0 * mean(r)
    sigstar <- sqrt(c0 / (K * lam))
    mu <- rnorm(1, mustar, sigstar)
    for (k in seq_len(K)) {
      rk <- r[k]
      rkp <- rnorm(1, rk, del[k])
      rp <- r
      rp[k] <- rkp
      hfun <- if (hf == 0) HFUN0 else HFUN
      hvalp <- hfun(
        n = n, k = k, r = rp, mu = mu, lam = lam,
        Cvec = Cvec, N1vec = N1vec, N2vec = N2vec,
        N0vec = N0vec, M1vec = M1vec, M2vec = M2vec,
        mvec = mvec, u = u, v = v
      )
      hval <- hfun(
        n = n, k = k, r = r, mu = mu, lam = lam,
        Cvec = Cvec, N1vec = N1vec, N2vec = N2vec,
        N0vec = N0vec, M1vec = M1vec, M2vec = M2vec,
        mvec = mvec, u = u, v = v
      )
      pk <- exp(hvalp - hval)
      if (runif(1) <= pk) {
        r[k] <- rkp
        if (j > B) rct[k] <- rct[k] + 1
      }
    }
    muv[1 + j] <- mu
    lamv[1 + j] <- lam
    rm[1 + j, ] <- r
  }
  list(muv = muv, lamv = lamv, rm = rm, rar = rct / J)
}

#' Metropolis-Hastings sampler with linear predictors
#'
#' Extends [MHALG()] so that the latent weight of each class depends on
#' covariates through a design matrix `X` and regression coefficients with a
#' multivariate normal prior.
#'
#' @param X Numeric matrix; design matrix linking classes to coefficients.
#' @param B0 Numeric vector; prior mean of the regression coefficients.
#' @param V0 Numeric matrix; prior covariance of the regression coefficients.
#' @inheritParams MHALG
#'
#' @return A list with components `betv` (sampled regression coefficients),
#'   `lamv`, `rm` and `rar`, as in [MHALG()].
#'
#' @references
#' Puza, B. and Bonfrer, A. (2018). A series of two-urn biased sampling
#' problems. *Communications in Statistics - Theory and Methods*, 47(1),
#' 80-91. \doi{10.1080/03610926.2017.1300282}
#'
#' @examples
#' set.seed(6215)
#' X <- cbind(1, c(1, 0, 1))
#' res <- RMMHALG(
#'   B = 2, J = 5,
#'   N1vec = c(150, 150, 150), N2vec = c(200, 200, 200),
#'   N0vec = c(50, 50, 50), M1vec = c(40, 40, 40), M2vec = c(60, 60, 60),
#'   mvec = c(10, 11, 9), Cvec = c(1, 2, 3),
#'   X = X, B0 = c(0, 0), V0 = diag(2)
#' )
#' str(res)
#'
#' @importFrom mvtnorm rmvnorm
#' @export
RMMHALG <- function(B = 5, J = 20, N1vec = 150, N2vec = 200, N0vec = 50,
                    M1vec = 40, M2vec = 60, mvec = 10,
                    Cvec = 1, X, B0, V0,
                    eta = 0.001, tau = 1000, mu0 = 0, sig0 = 0.001,
                    muinit = 0, rinit = 0, del = 0.1) {
  n <- length(mvec)
  K <- max(Cvec)
  bet <- matrix(B0, nrow = 1)
  betv <- bet
  invV0 <- solve(V0)
  tX <- t(X)
  tXX <- tX %*% X
  r <- rep(rinit, K)
  rct <- rep(0, K)
  if (length(del) == 1) del <- rep(del, K)
  rm <- matrix(0, nrow = J + B, ncol = K)
  lamv <- numeric()
  pb <- txtProgressBar(min = 0, max = J + B, style = 3)
  for (j in seq_len(B + J)) {
    setTxtProgressBar(pb, j)
    lam <- rgamma(1, eta + 0.5 * K, tau + 0.5 * sum((r - (X %*% t(bet)))^2))
    sigstar <- solve(invV0 + tXX * lam)
    mustar <- sigstar %*% (invV0 %*% B0 + lam * tX %*% r)
    bet <- mvtnorm::rmvnorm(1, mustar, sigstar)
    for (k in seq_len(K)) {
      rk <- r[k]
      rkp <- rnorm(1, rk, del[k])
      rp <- r
      rp[k] <- rkp
      mu <- (X %*% t(bet))[k]
      qk <- HFUN(
        n = n, k = k, r = rp, mu = mu, lam = lam,
        Cvec = Cvec, N1vec = N1vec, N2vec = N2vec,
        N0vec = N0vec, M1vec = M1vec, M2vec = M2vec, mvec = mvec
      ) -
        HFUN(
          n = n, k = k, r = r, mu = mu, lam = lam,
          Cvec = Cvec, N1vec = N1vec, N2vec = N2vec,
          N0vec = N0vec, M1vec = M1vec, M2vec = M2vec, mvec = mvec
        )
      pk <- exp(qk)
      if (runif(1) <= pk) {
        r[k] <- rkp
        if (j > B) rct[k] <- rct[k] + 1
      }
    }
    betv <- rbind(betv, bet)
    lamv <- c(lamv, lam)
    rm[j, ] <- r
  }
  list(betv = betv, lamv = lamv, rm = rm, rar = rct / J)
}

#' Metropolis-Hastings sampler with shared weights within classes
#'
#' Variant of [RMMHALG()] that infers the weight classes from the design
#' matrix `X` and constrains observations mapping to the same linear predictor
#' to share a common weight.
#'
#' @inheritParams RMMHALG
#' @param n Integer; number of observations (defaults to `length(mvec)`).
#'
#' @return A list with components `betv`, `lamv`, `rm` and `rar`, as in
#'   [RMMHALG()].
#'
#' @references
#' Puza, B. and Bonfrer, A. (2018). A series of two-urn biased sampling
#' problems. *Communications in Statistics - Theory and Methods*, 47(1),
#' 80-91. \doi{10.1080/03610926.2017.1300282}
#'
#' @examples
#' set.seed(8391)
#' X <- cbind(1, c(1, 0, 1))
#' res <- SMRMMHALG(
#'   B = 2, J = 5,
#'   N1vec = c(150, 150, 150), N2vec = c(200, 200, 200),
#'   N0vec = c(50, 50, 50), M1vec = c(40, 40, 40), M2vec = c(60, 60, 60),
#'   mvec = c(10, 11, 9),
#'   X = X, B0 = c(0, 0), V0 = diag(2)
#' )
#' str(res)
#'
#' @export
SMRMMHALG <- function(B = 5, J = 20, N1vec = 150, N2vec = 200, N0vec = 50,
                      M1vec = 40, M2vec = 60, mvec = 10,
                      X, B0, V0, eta = 0.001, tau = 1000, mu0 = 0, sig0 = 0.001,
                      muinit = 0, rinit = 0, del = 0.1, n = length(mvec)) {
  n <- length(mvec)
  DX <- ncol(X)
  BV <- matrix(log(1 + seq_len(DX)), nrow = 1)
  xb <- unique(X %*% t(BV))
  K <- length(xb)
  Cvec <- rep(NA, n)
  for (j in seq_len(n)) {
    for (k in seq_len(K)) {
      if ((X %*% t(BV))[j] == xb[k]) Cvec[j] <- k
    }
  }
  bet <- matrix(B0, nrow = 1)
  betv <- bet
  invV0 <- solve(V0)
  tX <- t(X)
  tXX <- tX %*% X
  r <- rep(0, K)
  rct <- rep(0, K)
  if (length(del) == 1) {
    del <- rep(del, K)
  } else {
    if (length(del) > K) {
      warning(sprintf("Only using first %s from del vector", K))
    }
    if (length(del) < K) stop(sprintf("Increase length of del vec to %s", K))
    del <- del[seq_len(K)]
  }
  rm <- matrix(0, nrow = J + B, ncol = K)
  lamv <- numeric()
  pb <- txtProgressBar(min = 0, max = J + B, style = 3)
  for (j in seq_len(B + J)) {
    setTxtProgressBar(pb, j)
    rwh <- rep(NA, n)
    for (k in seq_len(K)) rwh[which(Cvec == k)] <- r[k]
    rwh <- matrix(rwh, nrow = 1)
    lam <- rgamma(
      1, eta + 0.5 * K,
      tau + 0.5 * sum((unique(t(rwh) - (X %*% t(bet)))^2))
    )
    sigstar <- solve(invV0 + tXX * lam)
    mustar <- sigstar %*% (invV0 %*% B0 + lam * tX %*% t(rwh))
    bet <- mvtnorm::rmvnorm(1, mustar, sigstar)
    for (k in seq_len(K)) {
      rk <- r[k]
      rkp <- rnorm(1, rk, del[k])
      rp <- r
      rp[k] <- rkp
      mu <- (X %*% t(bet))[k]
      qk <- HFUN(
        n = n, k = k, r = rp, mu = mu, lam = lam,
        Cvec = Cvec, N1vec = N1vec, N2vec = N2vec,
        N0vec = N0vec, M1vec = M1vec, M2vec = M2vec, mvec = mvec
      ) -
        HFUN(
          n = n, k = k, r = r, mu = mu, lam = lam,
          Cvec = Cvec, N1vec = N1vec, N2vec = N2vec,
          N0vec = N0vec, M1vec = M1vec, M2vec = M2vec, mvec = mvec
        )
      pk <- exp(qk)
      if (runif(1) <= pk) {
        r[k] <- rkp
        if (j > B) rct[k] <- rct[k] + 1
      }
    }
    betv <- rbind(betv, bet)
    lamv <- c(lamv, lam)
    rm[j, ] <- r
  }
  list(betv = betv, lamv = lamv, rm = rm, rar = rct / J)
}

# Link function mapping the latent scale to a positive weight (odds-ratio form).
gFUN <- function(rval = 2.7, u = exp(20), v = 1) {
  u / (1 + u * exp(-v * rval))
}

# Inverse of the link function gFUN().
gINV <- function(wval = 15, u = exp(20), v = 1) {
  -log((1 / wval) - 1 / u) / v
}

# Log full conditional (up to a constant) for class k, on the log-weight scale.
HFUN0 <- function(k, r, mu, lam, Cvec, N1vec, N2vec, N0vec, M1vec, M2vec,
                  mvec, n = length(mvec), u = exp(20), v = 1) {
  rk <- r[k]
  term <- 0
  for (i in which(Cvec == k)) {
    .l <- log(max(1.0e-300, PROBM(
      m = mvec[i], N1 = N1vec[i], N2 = N2vec[i], N = N0vec[i],
      m1 = M1vec[i], m2 = M2vec[i], w = exp(rk)
    )))
    term <- term + .l
  }
  -0.5 * lam * (rk - mu)^2 + term
}

# Log full conditional (up to a constant) for class k, using the link gFUN().
HFUN <- function(n = length(mvec), k, r, mu, lam, Cvec, N1vec, N2vec, N0vec,
                 M1vec, M2vec, mvec, u = exp(20), v = 1) {
  rk <- r[k]
  term <- 0
  for (i in which(Cvec == k)) {
    .l <- log(max(1.0e-300, PROBM(
      m = mvec[i], N1 = N1vec[i], N2 = N2vec[i], N = N0vec[i],
      m1 = M1vec[i], m2 = M2vec[i], w = gFUN(rval = rk, u = u, v = v)
    )))
    term <- term + .l
  }
  term - 0.5 * lam * (rk - mu)^2
}
