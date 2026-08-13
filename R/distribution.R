#' Probability mass for the number of matching items
#'
#' Computes the probability that exactly `m` items match across two samples
#' drawn without replacement from two overlapping urns, under the biased
#' selection scheme of Puza and Bonfrer (2018).
#'
#' Urn 1 contains `N1` distinct items and urn 2 contains `N2` distinct items,
#' with `N` items common to both urns. A sample of size `m1` is drawn from
#' urn 1 and a sample of size `m2` from urn 2. Each item common to both urns
#' is given selection weight `w` while every other item has weight 1, so that
#' `w > 1` makes common items more likely to be selected. A match is an item
#' selected in both samples.
#'
#' @param m Integer; number of matches for which the probability is required.
#' @param N1 Integer; number of distinct items in urn 1.
#' @param N2 Integer; number of distinct items in urn 2.
#' @param N Integer; number of items common to both urns.
#' @param m1 Integer; sample size drawn from urn 1.
#' @param m2 Integer; sample size drawn from urn 2.
#' @param w Numeric (> 0); weight of a common item relative to a non-common
#'   item. `w = 1` gives unbiased simple random sampling without replacement.
#'
#' @return A single numeric probability.
#'
#' @details
#' The implementation follows Theorem 2 of Puza and Bonfrer (2018). It sums
#' over `y1`, the number of common items drawn from urn 1, weighting the
#' conditional match probability given `y1` by the probability of `y1`. The
#' distribution of `y1` is Wallenius' univariate non-central hypergeometric
#' distribution ([BiasedUrn::dWNCHypergeo()]); the conditional distribution of
#' matches given `y1` uses the bivariate Wallenius' distribution
#' ([BiasedUrn::dMWNCHypergeo()]). When `w = 1` the result reduces to the
#' hypergeometric mixture of Theorem 1.
#'
#' @references
#' Puza, B. and Bonfrer, A. (2018). A series of two-urn biased sampling
#' problems. *Communications in Statistics - Theory and Methods*, 47(1),
#' 80-91. \doi{10.1080/03610926.2017.1300282}
#'
#' @seealso [DISTM()] for the full distribution and [inferW()] for inference
#'   on `w`.
#'
#' @examples
#' # Probability of exactly 3 matches (Example 1 of the paper, w = 2)
#' PROBM(m = 3, N1 = 37, N2 = 45, N = 16, m1 = 12, m2 = 8, w = 2)
#'
#' @importFrom BiasedUrn dWNCHypergeo dMWNCHypergeo
#' @export
PROBM <- function(m = 3, N1 = 37, N2 = 45, N = 16, m1 = 12, m2 = 8, w = 2) {
  a <- max(0, m1 + m2 - N1 - N2 + N)
  b <- min(N, m1, m2)
  f <- 0
  if (m >= a && m <= b) {
    y1vec <- max(m, m1 - N1 + N):min(N, m1, N2 + m - m2)
    leny1vec <- length(y1vec)
    fmgiveny1vec <- rep(NA, leny1vec)
    fy1vec <- BiasedUrn::dWNCHypergeo(y1vec, N, N1 - N, m1, w)
    for (j in seq_len(leny1vec)) {
      fmgiveny1vec[j] <- fmgiveny1fun(m, y1vec[j], N2, N, m2, w)
    }
    f <- sum(fy1vec * fmgiveny1vec)
  }
  f
}

# Conditional probability of m matches given y1 common items drawn from urn 1.
# Sums the bivariate Wallenius' probability over the number q of common items
# drawn from urn 2 that are not among the y1 drawn from urn 1.
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
#' Computes the complete probability distribution of the number of matches
#' together with its mean and variance, for the biased two-urn sampling model.
#'
#' @inheritParams PROBM
#'
#' @return A list with components:
#' \describe{
#'   \item{mvec}{Integer vector of possible match counts.}
#'   \item{fmvec}{Probability masses corresponding to `mvec`.}
#'   \item{sumfmvec}{Sum of `fmvec`; a numerical check that should equal 1.}
#'   \item{Em}{Expected number of matches.}
#'   \item{Vm}{Variance of the number of matches.}
#'   \item{SDm}{Standard deviation of the number of matches.}
#' }
#'
#' @references
#' Puza, B. and Bonfrer, A. (2018). A series of two-urn biased sampling
#' problems. *Communications in Statistics - Theory and Methods*, 47(1),
#' 80-91. \doi{10.1080/03610926.2017.1300282}
#'
#' @seealso [PROBM()] for a single probability, [Vmf()] for the closed-form
#'   unbiased variance.
#'
#' @examples
#' # Full distribution for Example 1 of the paper (w = 2)
#' d <- DISTM(N1 = 37, N2 = 45, N = 16, m1 = 12, m2 = 8, w = 2)
#' d$sumfmvec # ~ 1
#' d$Em # expected number of matches
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

#' Closed-form variance of the number of matches under unbiased sampling
#'
#' Evaluates the closed-form variance of the number of matches when `w = 1`
#' (simple random sampling without replacement), given in Theorem 1(b) of
#' Puza and Bonfrer (2018). It is provided for comparison with the exact
#' variance returned by [DISTM()].
#'
#' @inheritParams PROBM
#'
#' @return A single numeric variance.
#'
#' @references
#' Puza, B. and Bonfrer, A. (2018). A series of two-urn biased sampling
#' problems. *Communications in Statistics - Theory and Methods*, 47(1),
#' 80-91. \doi{10.1080/03610926.2017.1300282}
#'
#' @examples
#' Vmf(N1 = 37, N2 = 45, N = 16, m1 = 12, m2 = 8)
#' # Matches the exact variance from DISTM() when w = 1:
#' DISTM(N1 = 37, N2 = 45, N = 16, m1 = 12, m2 = 8, w = 1)$Vm
#'
#' @export
Vmf <- function(N1 = 37, N2 = 45, N = 16, m1 = 12, m2 = 8) {
  c1 <- m1 / N1 * N
  c2 <- c1 * (N1 - N) / N1 * (N1 - m1) / (N1 - 1)
  c3 <- c1^2 + c2
  m2 * (N2 - m2) / (N2^2 * (N2 - 1)) * (N2 * c1 - c3) + m2^2 / N2^2 * c2
}
