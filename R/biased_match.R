#' Probability mass function for the number of matching items
#'
#' Given two urns containing overlapping sets of items, this function
#' computes the probability mass at a particular number of matches `m` when
#' samples are drawn without replacement from each urn with bias.  The
#' implementation uses Wallenius' and multivariate non‑central hypergeometric
#' distributions provided by the \pkg{BiasedUrn} package.
#'
#' Suppose urn 1 contains `N1` distinct items, urn 2 contains `N2` distinct
#' items and `N` items are common to both urns.  A sample of size `m1` is
#' drawn from urn 1 and a sample of size `m2` is drawn from urn 2.  Items
#' that are common to both urns are sampled with weight `w` relative to
#' non‑overlapping items.  The return value is the probability that the
#' number of matches between the two samples is exactly `m`.
#'
#' @param m Integer; the number of matches for which the probability is
#'   required.
#' @param N1 Integer; the number of distinct items in urn 1.
#' @param N2 Integer; the number of distinct items in urn 2.
#' @param N Integer; the number of items common to both urns.
#' @param m1 Integer; the size of the sample drawn from urn 1.
#' @param m2 Integer; the size of the sample drawn from urn 2.
#' @param w Numeric; the weight of a common item relative to a non‑common
#'   item.  A value greater than 1 indicates that common items are more likely
#'   to be selected.
#'
#' @return A single numeric value giving the probability mass at `m`.
#'
#' @details
#' The algorithm enumerates the possible values of `y1`, the number of
#' common items drawn from urn 1, and sums the probability of drawing
#' `m` matches conditional on `y1` multiplied by the probability of drawing
#' `y1` common items.  The distribution of `y1` is Wallenius' non‑central
#' hypergeometric distribution computed by
#' [BiasedUrn::dWNCHypergeo()], while the conditional distribution of
#' `m` given `y1` is computed using the multivariate Wallenius
#' distribution [BiasedUrn::dMWNCHypergeo()].
#'
#' @references
#' Puza, B. & Bonfrer, A. (2018). *A series of two‑urn biased sampling
#' problems*. Communications in Statistics – Theory and Methods, 47(1), 80–91.
#' The functions in this package implement the probability model and inference
#' routines described in that paper.
#'
#' @examples
#' # Probability of observing exactly 3 matches when both urns overlap by 16
#' # items and samples of size 12 and 8 are drawn with weight 2 for common
#' # items.
#' PROBM(m = 3, N1 = 37, N2 = 45, N = 16, m1 = 12, m2 = 8, w = 2)
#'
#' @importFrom BiasedUrn dWNCHypergeo dMWNCHypergeo
#' @export
PROBM <- function(m = 3, N1 = 37, N2 = 45, N = 16, m1 = 12, m2 = 8, w = 2) {
  # Range of possible matches
  a <- max(0, m1 + m2 - N1 - N2 + N)
  b <- min(N, m1, m2)
  f <- 0
  if (m >= a && m <= b) {
    # Vector of possible counts of common items drawn from urn 1
    y1vec <- max(m, m1 - N1 + N):min(N, m1, N2 + m - m2)
    leny1vec <- length(y1vec)
    fmgiveny1vec <- rep(NA, leny1vec)
    fy1vec <- BiasedUrn::dWNCHypergeo(y1vec, N, N1 - N, m1, w)
    # Conditional distribution of m given y1
    for (j in seq_len(leny1vec)) {
      fmgiveny1vec[j] <- fmgiveny1fun(m, y1vec[j], N2, N, m2, w)
    }
    f <- sum(fy1vec * fmgiveny1vec)
  }
  f
}

# Internal helper: conditional distribution of matches given y1
#
# Computes the sum over q of the multivariate Wallenius' probability for
# matching and non‑matching items.  This function is for internal use and
# not exported.
#
# @noRd
fmgiveny1fun <- function(m, y1, N2, N, m2, w) {
  A <- max(0, m2 - m - N2 + N)
  B <- min(N - y1, m2 - m)
  tot <- 0
  for (q in A:B) {
    tot <- tot + BiasedUrn::dMWNCHypergeo(
      c(m, q, m2 - m - q),
      c(y1, N - y1, N2 - N),
      m2,
      c(w, w, 1)
    )
  }
  tot
}

#' Distribution, mean and variance of the number of matches
#'
#' Computes the entire probability distribution of the number of matches for
#' given urn parameters together with the expected value and variance.
#'
#' @param N1 Integer; number of distinct items in urn 1.
#' @param N2 Integer; number of distinct items in urn 2.
#' @param N Integer; number of items common to both urns.
#' @param m1 Integer; sample size drawn from urn 1.
#' @param m2 Integer; sample size drawn from urn 2.
#' @param w Numeric; weight of common items relative to non‑common items.
#'
#' @return A list with components:
#'   \describe{
#'     \item{mvec}{Vector of possible numbers of matches.}
#'     \item{fmvec}{Vector of probability masses corresponding to `mvec`.}
#'     \item{sumfmvec}{The sum of the probability masses (numerical check).}
#'     \item{Em}{Expected number of matches.}
#'     \item{Vm}{Variance of the number of matches.}
#'     \item{SDm}{Standard deviation of the number of matches.}
#'   }
#'
#' @examples
#' # Compute full distribution when both urns overlap by 16 items and
#' # sample sizes are 12 and 8 with weight 2.
#' dist <- DISTM(N1 = 37, N2 = 45, N = 16, m1 = 12, m2 = 8, w = 2)
#' sum(dist$fmvec)  # Should be 1
#' dist$Em          # Expected number of matches
#'
#' @export
DISTM <- function(N1 = 37, N2 = 45, N = 16, m1 = 12, m2 = 8, w = 2) {
  a <- max(0, m1 + m2 - N1 - N2 + N)
  b <- min(N, m1, m2)
  mvec <- a:b
  fmvec <- rep(NA, length(mvec))
  for (i in seq_along(mvec)) {
    m <- mvec[i]
    y1vec <- max(m, m1 - N1 + N):min(N, m1, N2 + m - m2)
    leny1vec <- length(y1vec)
    fmgiveny1vec <- rep(NA, leny1vec)
    fy1vec <- BiasedUrn::dWNCHypergeo(y1vec, N, N1 - N, m1, w)
    for (j in seq_len(leny1vec)) {
      fmgiveny1vec[j] <- fmgiveny1fun(m, y1vec[j], N2, N, m2, w)
    }
    fmvec[i] <- sum(fy1vec * fmgiveny1vec)
  }
  Em <- sum(mvec * fmvec)
  Vm <- sum(mvec^2 * fmvec) - Em^2
  list(
    mvec = mvec,
    fmvec = fmvec,
    sumfmvec = sum(fmvec),
    Em = Em,
    Vm = Vm,
    SDm = sqrt(Vm)
  )
}

#' Metropolis–Hastings algorithm for a single weight class
#'
#' Implements a Metropolis–Hastings (MH) algorithm for estimating the
#' parameters of the biased matching model when all observations share a
#' common weight parameter.  This corresponds to the "Problem 4" algorithm
#' described by Puza and Bonfrer (2018).  The algorithm samples
#' from the posterior distributions of the latent weight parameters and
#' hyperparameters \eqn{\mu} and \eqn{\lambda} under conjugate priors.
#'
#' @param B Integer; number of burn‑in iterations.
#' @param J Integer; number of sampling iterations after burn‑in.
#' @param N1vec Numeric vector; values of `N1` for each observation.
#' @param N2vec Numeric vector; values of `N2` for each observation.
#' @param N0vec Numeric vector; values of `N` (number of common items) for each observation.
#' @param M1vec Numeric vector; values of `m1` (sample sizes from urn 1) for each observation.
#' @param M2vec Numeric vector; values of `m2` (sample sizes from urn 2) for each observation.
#' @param mvec Numeric vector; observed numbers of matches.
#' @param Cvec Integer vector; indicator assigning each observation to a weight class (should all be 1 for a single class).
#' @param eta,tau,mu0,sig0 Numeric scalars; hyperparameters for priors on \eqn{\lambda} and \eqn{\mu}.
#' @param muinit,rinit Numeric; starting values for \eqn{\mu} and the weight parameters.
#' @param del Numeric vector of tuning constants for proposal variances.
#' @param u,v Numeric; parameters controlling the link function for the weight parameter (`gFUN`).
#' @param hf Integer; if 0, use an older version of the helper function, otherwise the default.
#'
#' @return A list with components:
#'   \describe{
#'     \item{muv}{A vector of sampled \eqn{\mu} values of length `B+J+1`.}
#'     \item{lamv}{A vector of sampled \eqn{\lambda} values of length `B+J+1`.}
#'     \item{rm}{A matrix of sampled weight parameters with dimensions `(B+J+1) × K`,
#'       where `K = max(Cvec)`.  Each column corresponds to a weight class.}
#'     \item{rar}{A vector of acceptance rates for each weight class.}
#'   }
#'
#' @details
#' The algorithm iteratively samples \eqn{\lambda} from its gamma full
#' conditional, \eqn{\mu} from its normal full conditional and proposes new
#' values of the weight parameters via a random‑walk normal proposal.  The
#' proposed weight is transformed to a bias parameter using the link
#' function [gFUN()], and the acceptance probability is computed
#' using [HFUN()] (or [HFUN0()] when `hf = 0`).
#'
#' This function is primarily intended for internal use and for reproducing
#' the results of Puza and Bonfrer (2018).  Users may prefer to call
#' [SMRMMHALG()] when multiple weight classes share common weights.
#'
#' @examples
#' # A short MH run on synthetic data (for illustration only; increase J for
#' # real inference)
#' set.seed(123)
#' res <- MHALG(B = 2, J = 10,
#'              N1vec = 37, N2vec = 45, N0vec = 16,
#'              M1vec = 12, M2vec = 8, mvec = 3,
#'              Cvec = 1, del = 0.1)
#' head(res$muv)
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
  # preallocate
  rm <- matrix(0, nrow = 1 + B + J, ncol = K)
  rm[1, ] <- r
  muv <- lamv <- rep(0, 1 + B + J)
  muv[1] <- mu
  lamv[1] <- NA
  pb <- txtProgressBar(min = 0, max = B + J, style = 3)
  for (j in seq_len(B + J)) {
    setTxtProgressBar(pb, j)
    # Sample lambda
    lam <- rgamma(1, eta + 0.5 * K, tau + 0.5 * sum((r - mu)^2))
    # Sample mu
    c0 <- K / (K + 1 / (lam * sig02))
    mustar <- (1 - c0) * mu0 + c0 * mean(r)
    sigstar <- sqrt(c0 / (K * lam))
    mu <- rnorm(1, mustar, sigstar)
    # Update r
    for (k in seq_len(K)) {
      rk <- r[k]
      rkp <- rnorm(1, rk, del[k])
      rp <- r
      rp[k] <- rkp
      # Select helper function
      if (hf == 0) {
        hfun <- HFUN0
      } else {
        hfun <- HFUN
      }
      hvalp <- hfun(n = n, k = k, r = rp, mu = mu, lam = lam,
                    Cvec = Cvec, N1vec = N1vec, N2vec = N2vec,
                    N0vec = N0vec, M1vec = M1vec, M2vec = M2vec,
                    mvec = mvec, u = u, v = v)
      hval <- hfun(n = n, k = k, r = r, mu = mu, lam = lam,
                    Cvec = Cvec, N1vec = N1vec, N2vec = N2vec,
                    N0vec = N0vec, M1vec = M1vec, M2vec = M2vec,
                    mvec = mvec, u = u, v = v)
      qk <- hvalp - hval
      pk <- exp(qk)
      if (runif(1) <= pk) {
        r[k] <- rkp
        if (j > B) rct[k] <- rct[k] + 1
      }
    }
    # Store results
    muv[1 + j] <- mu
    lamv[1 + j] <- lam
    rm[1 + j, ] <- r
  }
  list(muv = muv, lamv = lamv, rm = rm, rar = rct / J)
}

#' Metropolis–Hastings algorithm with linear predictors
#'
#' Extends [MHALG()] by allowing the weight parameters to depend on a
#' design matrix `X` via a multivariate normal prior.  Each observation
#' belongs to one of `K` weight classes, and the latent log‑weights are
#' modelled as `X %*% beta`.  See Puza and Bonfrer (2018) for details.
#'
#' @param X Numeric matrix; design matrix linking observations to regression
#'   coefficients.
#' @param B0 Numeric vector; prior mean of regression coefficients.
#' @param V0 Numeric matrix; prior covariance of regression coefficients.
#' @inheritParams MHALG
#'
#' @return A list similar to that returned by [MHALG()] with additional
#'   components: `betv` (matrix of sampled regression coefficients), and
#'   `lamv`, `rm`, `rar` as before.
#'
#' @examples
#' # Simulate a small example with a design matrix of two covariates
#' set.seed(1)
#' X <- cbind(1, c(1, 0, 1))
#' res <- RMMHALG(B = 2, J = 5,
#'                N1vec = c(150, 150, 150), N2vec = c(200, 200, 200),
#'                N0vec = c(50, 50, 50), M1vec = c(40, 40, 40), M2vec = c(60, 60, 60),
#'                mvec = c(10, 11, 9), Cvec = c(1, 2, 3),
#'                X = X, B0 = c(0, 0), V0 = diag(2))
#' head(res$betv)
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
  bet <- B0
  betv <- matrix(bet, nrow = 1)
  invV0 <- solve(V0)
  tX <- t(X)
  tXX <- tX %*% X
  sig02 <- sig0^2
  mu <- muinit
  lam <- NA
  r <- rep(rinit, K)
  rct <- rep(0, K)
  if (length(del) == 1) del <- rep(del, K)
  rm <- matrix(0, nrow = J + B, ncol = K)
  lamv <- numeric()
  muv <- numeric(J + B)
  pb <- txtProgressBar(min = 0, max = J + B, style = 3)
  for (j in seq_len(B + J)) {
    setTxtProgressBar(pb, j)
    # Sample lambda
    lam <- rgamma(1, eta + 0.5 * K, tau + 0.5 * sum((t(r) - (X %*% t(bet)))^2))
    # Sample beta
    sigstar <- solve(invV0 + tXX * lam)
    mustar <- sigstar %*% (invV0 %*% B0 + lam * tX %*% r)
    bet <- mvtnorm::rmvnorm(1, mustar, sigstar)
    # Update r
    for (k in seq_len(K)) {
      rk <- r[k]
      rkp <- rnorm(1, rk, del[k])
      rp <- r
      rp[k] <- rkp
      # compute mu for this class
      mu <- (X %*% t(bet))[k]
      qk <- HFUN(n = n, k = k, r = rp, mu = mu, lam = lam,
                 Cvec = Cvec, N1vec = N1vec, N2vec = N2vec,
                 N0vec = N0vec, M1vec = M1vec, M2vec = M2vec,
                 mvec = mvec) -
        HFUN(n = n, k = k, r = r, mu = mu, lam = lam,
             Cvec = Cvec, N1vec = N1vec, N2vec = N2vec,
             N0vec = N0vec, M1vec = M1vec, M2vec = M2vec,
             mvec = mvec)
      pk <- exp(qk)
      if (runif(1) <= pk) {
        r[k] <- rkp
        if (j > B) rct[k] <- rct[k] + 1
      }
    }
    # collect results
    betv <- rbind(betv, bet)
    lamv <- c(lamv, lam)
    rm[j, ] <- r
  }
  list(betv = betv, lamv = lamv, rm = rm, rar = rct / J)
}

#' Metropolis–Hastings algorithm with common weights
#'
#' Implements the Metropolis–Hastings algorithm for the case where observations
#' sharing the same linear predictor are constrained to have equal weights.  This
#' variant infers a set of unique weight classes based on the design matrix `X`
#' and enforces common weights within each class.  See Puza and Bonfrer (2018)
#' for details.
#'
#' @inheritParams RMMHALG
#' @return A list with components `betv`, `lamv`, `rm` and `rar` analogous to
#'   those returned by [RMMHALG()].
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
  bet <- B0
  betv <- matrix(bet, nrow = 1)
  invV0 <- solve(V0)
  tX <- t(X)
  tXX <- tX %*% X
  r <- rep(0, K)
  rct <- rep(0, K)
  if (length(del) == 1) {
    del <- rep(del, K)
  } else {
    if (length(del) > K) warning(sprintf("Only using first %s from del vector", K))
    if (length(del) < K) stop(sprintf("Increase length of del vec to %s", K))
    del <- del[seq_len(K)]
  }
  rm <- matrix(0, nrow = J + B, ncol = K)
  lamv <- numeric()
  muv <- numeric(J + B)
  pb <- txtProgressBar(min = 0, max = J + B, style = 3)
  for (j in seq_len(B + J)) {
    setTxtProgressBar(pb, j)
    # replicate r for each observation
    rwh <- rep(NA, n)
    for (k in seq_len(K)) rwh[which(Cvec == k)] <- r[k]
    rwh <- matrix(rwh, nrow = 1)
    lam <- rgamma(1, eta + 0.5 * K,
                  tau + 0.5 * sum((unique(t(rwh) - (X %*% t(bet)))^2)))
    sigstar <- solve(invV0 + tXX * lam)
    mustar <- sigstar %*% (invV0 %*% B0 + lam * tX %*% t(rwh))
    bet <- mvtnorm::rmvnorm(1, mustar, sigstar)
    # update r
    for (k in seq_len(K)) {
      rk <- r[k]
      rkp <- rnorm(1, rk, del[k])
      rp <- r
      rp[k] <- rkp
      mu <- (X %*% t(bet))[k]
      qk <- HFUN(n = n, k = k, r = rp, mu = mu, lam = lam,
                 Cvec = Cvec, N1vec = N1vec, N2vec = N2vec,
                 N0vec = N0vec, M1vec = M1vec, M2vec = M2vec,
                 mvec = mvec) -
        HFUN(n = n, k = k, r = r, mu = mu, lam = lam,
             Cvec = Cvec, N1vec = N1vec, N2vec = N2vec,
             N0vec = N0vec, M1vec = M1vec, M2vec = M2vec,
             mvec = mvec)
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

#' Summarise MH algorithm output
#'
#' Provides point estimates, modes and credible intervals from the output of
#' [MHALG()], [RMMHALG()] or [SMRMMHALG()].  If `sim = TRUE` the function
#' also compares the estimates to supplied true parameter values; otherwise
#' only sample summaries are returned.
#'
#' @param res List; result returned by one of the MH algorithms.
#' @param Bn Integer; number of burn‑in iterations to discard.
#' @param Jn Integer; number of posterior samples to retain for summarisation.
#' @param rm.idx Integer vector; which columns of `res$rm` correspond to the
#'   weight classes of interest.  If `NULL`, all classes are used.
#' @param sim Logical; if `TRUE`, compute a comparison with known parameters.
#'
#' @return A `data.frame` with rows corresponding to the posterior mean,
#'   mode and quantiles (0.025 and 0.975) of each parameter.  Columns
#'   correspond to \eqn{\mu}, \eqn{\lambda}, \eqn{\sigma} and the weight
#'   parameters \eqn{r_k}.
#'
#' @export
restable <- function(res, Bn = 2, Jn = length(res$muv) - Bn,
                     rm.idx = 1:dim(res$rm)[2], sim = FALSE) {
  if (is.null(rm.idx)) K <- 1 else K <- length(rm.idx)
  Jn <- (Bn + Jn)
  den <- density(res$muv[Bn:Jn], adjust = 1.5)
  imode <- which.max(den$y)
  mumode <- den$x[imode]
  muci <- quantile(res$muv[Bn:Jn], c(0.025, 0.975))
  den <- density(res$lamv[Bn:Jn], adjust = 1.5)
  imode <- which.max(den$y)
  lammode <- den$x[imode]
  lamci <- quantile(res$lamv[Bn:Jn], c(0.025, 0.975))
  den <- density(1 / sqrt(res$lamv[Bn:Jn]), adjust = 1.5)
  imode <- which.max(den$y)
  sigmamode <- den$x[imode]
  sigmaci <- quantile(1 / sqrt(res$lamv)[Bn:Jn], c(0.025, 0.975))
  rci <- vector("list", K)
  rmodevec <- c(mumode, lammode, sigmamode)
  for (k in seq_len(K)) {
    den <- density(res$rm[Bn:Jn, rm.idx[k]], adjust = 1.5)
    imode <- which.max(den$y)
    rmodevec <- c(rmodevec, den$x[imode])
    rci[[k]] <- quantile(res$rm[Bn:Jn, rm.idx[k]], c(0.025, 0.975))
  }
  if (sim) {
    stop("Simulation comparison is not implemented in this simplified package")
  } else {
    if (K > 1) .rmm <- colMeans(res$rm[Bn:Jn, rm.idx]) else .rmm <- mean(res$rm[Bn:Jn])
    estparmat <- c(mean(res$muv[Bn:Jn]), mean(res$lamv[Bn:Jn]), mean(1 / sqrt(res$lamv[Bn:Jn])), .rmm)
  }
  estparmat <- rbind(estparmat, rmodevec)
  .vec <- c(muci[1], lamci[1], sigmaci[1])
  for (k in seq_len(K)) .vec <- c(.vec, rci[[k]][1])
  estparmat <- rbind(estparmat, .vec)
  .vec <- c(muci[2], lamci[2], sigmaci[2])
  for (k in seq_len(K)) .vec <- c(.vec, rci[[k]][2])
  estparmat <- rbind(estparmat, .vec)
  rn <- c("Mean", "Mode", "Q025", "Q975")
  estparmat <- data.frame(estparmat, row.names = rn)
  colnames(estparmat) <- c("mu", "lambda", "sigma", paste0("r", 1:K))
  estparmat
}

#' Diagnostic tests for convergence and serial correlation
#'
#' Applies a suite of diagnostic tests to the posterior draws returned by
#' [MHALG()], [RMMHALG()] or [SMRMMHALG()].  The tests include the
#' rank von Neumann serial correlation test, the AR1 Yule–Walker test,
#' Anderson–Darling, Cramér–von Mises and Lilliefors normality tests from
#' the \pkg{nortest} package and the Geweke convergence diagnostic from
#' \pkg{coda}.
#'
#' @param resg List; result returned by one of the MH algorithms.
#' @param Bn Integer; number of burn‑in iterations to discard.
#' @param Jn Integer; total number of iterations used in the diagnostics.
#' @param s Integer; starting index for thinning the chain.
#' @param step Integer; step size for thinning.
#'
#' @return A `data.frame` where each row corresponds to a diagnostic test and
#'   each column corresponds to a parameter (\eqn{\mu}, \eqn{\lambda}, \eqn{r_k}).
#'   Entries contain p‑values from the tests.
#'
#' @importFrom EnvStats serialCorrelationTest
#' @importFrom nortest ad.test cvm.test lillie.test
#' @importFrom coda geweke.diag
#' @importFrom stats pnorm
#' @export
estprop <- function(resg, Bn = 500, Jn = 2000, s = 4, step = 10) {
  T <- nrow(resg$rm)
  K <- ncol(resg$rm)
  if (Bn + Jn > T) stop(sprintf("Bn+Jn too large: T = %s", T))
  tg1 <- resg$muv[-seq_len(Bn)][seq(s, Jn, step)]
  tg2 <- resg$lamv[-seq_len(Bn)][seq(s, Jn, step)]
  if (K > 1) {
    tg3 <- resg$rm[-seq_len(Bn), ][seq(s, Jn, step), , drop = FALSE]
    .t1 <- apply(tg3, 2, function(x) EnvStats::serialCorrelationTest(x, test = "rank.von.Neumann")$p.value)
    .t2 <- apply(tg3, 2, function(x) EnvStats::serialCorrelationTest(x, test = "AR1.yw")$p.value)
    .t3 <- apply(tg3, 2, function(x) nortest::ad.test(x)$p.value)
    .t4 <- apply(tg3, 2, function(x) nortest::cvm.test(x)$p.value)
    .t5 <- apply(tg3, 2, function(x) nortest::lillie.test(x)$p.value)
    .t6 <- apply(tg3, 2, function(x) {
      z <- coda::geweke.diag(x)$z
      min((1 - pnorm(z)) * 2, (1 - pnorm(-z)) * 2)
    })
  } else {
    tg3 <- resg$rm[-seq_len(Bn)][seq(s, Jn, step)]
    .t1 <- EnvStats::serialCorrelationTest(tg3, test = "rank.von.Neumann")$p.value
    .t2 <- EnvStats::serialCorrelationTest(tg3, test = "AR1.yw")$p.value
    .t3 <- nortest::ad.test(tg3)$p.value
    .t4 <- nortest::cvm.test(tg3)$p.value
    .t5 <- nortest::lillie.test(tg3)$p.value
    z <- coda::geweke.diag(tg3)$z
    .t6 <- min((1 - pnorm(z)) * 2, (1 - pnorm(-z)) * 2)
  }
  s1 <- c(EnvStats::serialCorrelationTest(tg1, test = "rank.von.Neumann")$p.value,
          EnvStats::serialCorrelationTest(tg2, test = "rank.von.Neumann")$p.value,
          .t1)
  s2 <- c(EnvStats::serialCorrelationTest(tg1, test = "AR1.yw")$p.value,
          EnvStats::serialCorrelationTest(tg2, test = "AR1.yw")$p.value,
          .t2)
  s3 <- c(nortest::ad.test(tg1)$p.value, nortest::ad.test(tg2)$p.value, .t3)
  s4 <- c(nortest::cvm.test(tg1)$p.value, nortest::cvm.test(tg2)$p.value, .t4)
  s5 <- c(nortest::lillie.test(tg1)$p.value, nortest::lillie.test(tg2)$p.value, .t5)
  s6 <- c({z <- coda::geweke.diag(tg1)$z; min((1 - pnorm(z)) * 2, (1 - pnorm(-z)) * 2)},
          {z <- coda::geweke.diag(tg2)$z; min((1 - pnorm(z)) * 2, (1 - pnorm(-z)) * 2)},
          .t6)
  rn <- c("SerialCorrelationTest[rank.von.Neumann]",
          "SerialCorrelationTest[AR1.yw]",
          "NormalityTest[Anderson-Darling]",
          "NormalityTest[Cramér-von Mises]",
          "NormalityTest[Lilliefors]",
          "ConvergenceTest[Geweke]")
  estconv <- rbind(s1, s2, s3, s4, s5, s6)
  estconv <- data.frame(estconv, row.names = rn)
  colnames(estconv) <- c("mu", "lambda", paste0("r", 1:K))
  estconv
}

#' Approximate variance formula
#'
#' Computes the approximate variance of the number of matches when `w = 1`.
#' This formula is derived under the assumption of unbiased sampling and is
#' provided for comparison with the exact variance returned by [DISTM()].
#'
#' @param N1 Integer; number of distinct items in urn 1.
#' @param N2 Integer; number of distinct items in urn 2.
#' @param N Integer; number of items common to both urns.
#' @param m1 Integer; sample size drawn from urn 1.
#' @param m2 Integer; sample size drawn from urn 2.
#'
#' @return Numeric; the approximate variance.
#'
#' @examples
#' # Approximate variance for unbiased sampling
#' Vmf(N1 = 37, N2 = 45, N = 16, m1 = 12, m2 = 8)
#'
#' @export
Vmf <- function(N1 = 37, N2 = 45, N = 16, m1 = 12, m2 = 8) {
  c1 <- m1 / N1 * N
  c2 <- c1 * (N1 - N) / N1 * (N1 - m1) / (N1 - 1)
  c3 <- c1^2 + c2
  m2 * (N2 - m2) / (N2^2 * (N2 - 1)) * (N2 * c1 - c3) + m2^2 / N2^2 * c2
}

# Helper functions ---------------------------------------------------------

# Link function mapping log‑weights to odds ratio
gFUN <- function(rval = 2.7, u = exp(20), v = 1) {
  u / (1 + u * exp(-v * rval))
}

# Inverse of the link function
gINV <- function(wval = 15, u = exp(20), v = 1) {
  -log((1 / wval) - 1 / u) / v
}

# Original HFUN (for backward compatibility)
HFUN0 <- function(k, r, mu, lam, Cvec, N1vec, N2vec, N0vec, M1vec, M2vec, mvec, n = length(mvec)) {
  rk <- r[k]
  term <- 0
  for (i in which(Cvec == k)) {
    .l <- log(max(1.0e-300, PROBM(m = mvec[i], N1 = N1vec[i], N2 = N2vec[i], N = N0vec[i],
                                 m1 = M1vec[i], m2 = M2vec[i], w = exp(rk))))
    term <- term + .l
  }
  -0.5 * lam * (rk - mu)^2 + term
}

# Fixed HFUN using link function
HFUN <- function(n = length(mvec), k, r, mu, lam, Cvec, N1vec, N2vec, N0vec, M1vec, M2vec, mvec, u = exp(20), v = 1) {
  rk <- r[k]
  term <- 0
  for (i in which(Cvec == k)) {
    .l <- log(max(1.0e-300, PROBM(m = mvec[i], N1 = N1vec[i], N2 = N2vec[i], N = N0vec[i],
                                  m1 = M1vec[i], m2 = M2vec[i], w = gFUN(rval = rk, u = u, v = v))))
    term <- term + .l
  }
  term - 0.5 * lam * (rk - mu)^2
}
