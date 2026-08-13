#' Likelihood and posterior curves for the bias parameter
#'
#' Evaluates, over a grid of weight values, the likelihood of the bias
#' parameter `w` given an observed number of matches `m`, together with the
#' kernel of the posterior density under the prior `f(w)` proportional to
#' `1 / w`. Both curves are returned standardised to a maximum of 1, matching
#' the figures in Puza and Bonfrer (2018).
#'
#' @inheritParams PROBM
#' @param wgrid Numeric vector of positive weight values at which to evaluate
#'   the curves.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{w}{The weight grid (as supplied).}
#'   \item{likelihood}{Likelihood values standardised to a maximum of 1.}
#'   \item{posterior}{Posterior-kernel values (likelihood / `w`) standardised
#'     to a maximum of 1.}
#' }
#'
#' @references
#' Puza, B. and Bonfrer, A. (2018). A series of two-urn biased sampling
#' problems. *Communications in Statistics - Theory and Methods*, 47(1),
#' 80-91. \doi{10.1080/03610926.2017.1300282}
#'
#' @seealso [inferW()] for point estimates and the test of `w = 1`.
#'
#' @examples
#' # Reproduce the likelihood/posterior curves of Problem 3 (marketing data)
#' curves <- likW(
#'   m = 13, N1 = 292, N2 = 350, N = 192, m1 = 90, m2 = 59,
#'   wgrid = seq(0.05, 8, by = 0.05)
#' )
#' plot(curves$w, curves$likelihood, type = "l",
#'   xlab = "w", ylab = "standardised value")
#' lines(curves$w, curves$posterior, lty = 2)
#'
#' @export
likW <- function(m, N1, N2, N, m1, m2, wgrid = seq(0.05, 6, by = 0.05)) {
  if (any(wgrid <= 0)) {
    stop("`wgrid` must contain only positive values.", call. = FALSE)
  }
  lik <- vapply(
    wgrid,
    function(w) PROBM(m = m, N1 = N1, N2 = N2, N = N, m1 = m1, m2 = m2, w = w),
    numeric(1)
  )
  post <- lik / wgrid
  data.frame(
    w = wgrid,
    likelihood = lik / max(lik),
    posterior = post / max(post)
  )
}

#' Inference on the bias parameter from an observed number of matches
#'
#' "Inverts" the biased sampling model of [PROBM()]: given the fixed urn
#' constants and a single observed number of matches `m`, it returns the
#' maximum likelihood estimate and posterior mode of the bias parameter `w`,
#' and a p-value for testing no bias against a positive bias. This implements
#' the inferential problems (Problems 1-3) of Puza and Bonfrer (2018).
#'
#' @inheritParams PROBM
#' @param prior Character; the prior on `w` used for the posterior mode.
#'   `"jeffreys"` (the default) uses `f(w)` proportional to `1 / w`, as in the
#'   paper; `"flat"` uses a uniform prior, so the posterior mode equals the
#'   MLE.
#' @param wmax Numeric; upper bound of the search region for `w`. If the
#'   likelihood is still increasing at `wmax`, the MLE is reported as
#'   infinite.
#' @param alternative Character; the alternative hypothesis for the p-value.
#'   `"greater"` (default) tests `w = 1` against `w > 1` (bias towards common
#'   items); `"less"` tests against `w < 1`.
#'
#' @return A list with components:
#' \describe{
#'   \item{mle}{Maximum likelihood estimate of `w` (possibly `Inf`).}
#'   \item{posterior_mode}{Posterior mode of `w` under `prior`.}
#'   \item{p_value}{P-value for testing `w = 1`, computed from the null
#'     (`w = 1`) distribution of the number of matches.}
#'   \item{alternative}{The alternative hypothesis used.}
#'   \item{m}{The observed number of matches.}
#' }
#'
#' @details
#' The likelihood is `L(w) = PROBM(m, ..., w)` and the posterior kernel under
#' the Jeffreys-type prior is `L(w) / w`. Both are maximised numerically on
#' `(0, wmax]`. As noted in the paper, when the observed `m` is large enough
#' the likelihood is strictly increasing and the MLE is infinite; this is
#' detected and reported as `Inf`, while the posterior mode remains finite.
#' The p-value is the null-hypothesis (`w = 1`) tail probability of the number
#' of matches, obtained from [DISTM()].
#'
#' @references
#' Puza, B. and Bonfrer, A. (2018). A series of two-urn biased sampling
#' problems. *Communications in Statistics - Theory and Methods*, 47(1),
#' 80-91. \doi{10.1080/03610926.2017.1300282}
#'
#' @seealso [likW()] for the full likelihood and posterior curves.
#'
#' @examples
#' # Problem 3 of the paper (masked retailer promotion data)
#' inferW(m = 13, N1 = 292, N2 = 350, N = 192, m1 = 90, m2 = 59)
#'
#' # Problem 2 of the paper (gender bias in committee selection, B = 2)
#' inferW(m = 2, N1 = 52, N2 = 52, N = 22, m1 = 9, m2 = 6)
#'
#' @importFrom stats optimize
#' @export
inferW <- function(m, N1, N2, N, m1, m2,
                   prior = c("jeffreys", "flat"),
                   wmax = 200,
                   alternative = c("greater", "less")) {
  prior <- match.arg(prior)
  alternative <- match.arg(alternative)
  eps <- 1e-8

  loglik <- function(w) {
    log(PROBM(m = m, N1 = N1, N2 = N2, N = N, m1 = m1, m2 = m2, w = w))
  }

  # Maximum likelihood estimate. If the log-likelihood is still increasing at
  # the upper search bound, the MLE is infinite.
  opt_mle <- stats::optimize(loglik, c(eps, wmax), maximum = TRUE)
  if (loglik(wmax) >= loglik(wmax * (1 - 1e-4))) {
    mle <- Inf
  } else {
    mle <- opt_mle$maximum
  }

  # Posterior mode.
  if (prior == "jeffreys") {
    logpost <- function(w) loglik(w) - log(w)
  } else {
    logpost <- loglik
  }
  posterior_mode <- stats::optimize(logpost, c(eps, wmax), maximum = TRUE)$maximum

  # P-value from the null (w = 1) distribution of the number of matches.
  d0 <- DISTM(N1 = N1, N2 = N2, N = N, m1 = m1, m2 = m2, w = 1)
  if (alternative == "greater") {
    p_value <- sum(d0$fmvec[d0$mvec >= m])
  } else {
    p_value <- sum(d0$fmvec[d0$mvec <= m])
  }

  list(
    mle = mle,
    posterior_mode = posterior_mode,
    p_value = p_value,
    alternative = alternative,
    m = m
  )
}
