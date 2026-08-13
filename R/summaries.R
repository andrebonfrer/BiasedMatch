#' Summarise posterior draws from an MH sampler
#'
#' Produces posterior means, modes and 95% credible intervals for the
#' hyperparameters and weight classes returned by [MHALG()], [RMMHALG()] or
#' [SMRMMHALG()].
#'
#' @param res List; output from one of the MH samplers.
#' @param Bn Integer; number of initial iterations to discard as burn-in.
#' @param Jn Integer; number of retained iterations to summarise.
#' @param rm.idx Integer vector; columns of `res$rm` (weight classes) to
#'   summarise.
#' @param sim Logical; reserved for comparison against known parameters. Not
#'   implemented and must be `FALSE`.
#'
#' @return A data frame with rows `Mean`, `Mode`, `Q025` and `Q975`, and
#'   columns for `mu`, `lambda`, `sigma` and each weight class `r1`, `r2`, ...
#'
#' @seealso [MHALG()], [estprop()].
#'
#' @examples
#' set.seed(5093)
#' res <- MHALG(
#'   B = 5, J = 40,
#'   N1vec = 37, N2vec = 45, N0vec = 16,
#'   M1vec = 12, M2vec = 8, mvec = 3, Cvec = 1, del = 0.1
#' )
#' restable(res, Bn = 5)
#'
#' @importFrom stats density quantile
#' @export
restable <- function(res, Bn = 2, Jn = length(res$muv) - Bn,
                     rm.idx = 1:dim(res$rm)[2], sim = FALSE) {
  if (is.null(rm.idx)) K <- 1 else K <- length(rm.idx)
  Jn <- (Bn + Jn)
  den <- density(res$muv[Bn:Jn], adjust = 1.5)
  mumode <- den$x[which.max(den$y)]
  muci <- quantile(res$muv[Bn:Jn], c(0.025, 0.975))
  den <- density(res$lamv[Bn:Jn], adjust = 1.5)
  lammode <- den$x[which.max(den$y)]
  lamci <- quantile(res$lamv[Bn:Jn], c(0.025, 0.975))
  den <- density(1 / sqrt(res$lamv[Bn:Jn]), adjust = 1.5)
  sigmamode <- den$x[which.max(den$y)]
  sigmaci <- quantile(1 / sqrt(res$lamv)[Bn:Jn], c(0.025, 0.975))
  rci <- vector("list", K)
  rmodevec <- c(mumode, lammode, sigmamode)
  for (k in seq_len(K)) {
    den <- density(res$rm[Bn:Jn, rm.idx[k]], adjust = 1.5)
    rmodevec <- c(rmodevec, den$x[which.max(den$y)])
    rci[[k]] <- quantile(res$rm[Bn:Jn, rm.idx[k]], c(0.025, 0.975))
  }
  if (sim) {
    stop("Simulation comparison is not implemented.", call. = FALSE)
  }
  if (K > 1) {
    .rmm <- colMeans(res$rm[Bn:Jn, rm.idx])
  } else {
    .rmm <- mean(res$rm[Bn:Jn])
  }
  estparmat <- c(
    mean(res$muv[Bn:Jn]), mean(res$lamv[Bn:Jn]),
    mean(1 / sqrt(res$lamv[Bn:Jn])), .rmm
  )
  estparmat <- rbind(estparmat, rmodevec)
  .vec <- c(muci[1], lamci[1], sigmaci[1])
  for (k in seq_len(K)) .vec <- c(.vec, rci[[k]][1])
  estparmat <- rbind(estparmat, .vec)
  .vec <- c(muci[2], lamci[2], sigmaci[2])
  for (k in seq_len(K)) .vec <- c(.vec, rci[[k]][2])
  estparmat <- rbind(estparmat, .vec)
  estparmat <- data.frame(estparmat, row.names = c("Mean", "Mode", "Q025", "Q975"))
  colnames(estparmat) <- c("mu", "lambda", "sigma", paste0("r", 1:K))
  estparmat
}

#' Convergence and serial-correlation diagnostics for MH draws
#'
#' Applies a suite of diagnostic tests to the thinned posterior draws returned
#' by [MHALG()], [RMMHALG()] or [SMRMMHALG()]: rank von Neumann and AR(1)
#' Yule-Walker serial-correlation tests (\pkg{EnvStats}), the Anderson-Darling,
#' Cramer-von Mises and Lilliefors normality tests (\pkg{nortest}), and the
#' Geweke convergence diagnostic (\pkg{coda}).
#'
#' @param resg List; output from one of the MH samplers.
#' @param Bn Integer; number of initial iterations to discard as burn-in.
#' @param Jn Integer; number of post-burn-in iterations to use.
#' @param s Integer; starting index for thinning.
#' @param step Integer; thinning step size.
#'
#' @return A data frame of p-values with one row per diagnostic test and one
#'   column per parameter (`mu`, `lambda`, and each weight class).
#'
#' @seealso [MHALG()], [restable()].
#'
#' @examples
#' set.seed(7204)
#' res <- MHALG(
#'   B = 10, J = 200,
#'   N1vec = 37, N2vec = 45, N0vec = 16,
#'   M1vec = 12, M2vec = 8, mvec = 3, Cvec = 1, del = 0.1
#' )
#' estprop(res, Bn = 10, Jn = 150, s = 1, step = 2)
#'
#' @importFrom EnvStats serialCorrelationTest
#' @importFrom nortest ad.test cvm.test lillie.test
#' @importFrom coda geweke.diag
#' @importFrom stats pnorm
#' @export
estprop <- function(resg, Bn = 500, Jn = 2000, s = 4, step = 10) {
  T <- nrow(resg$rm)
  K <- ncol(resg$rm)
  if (Bn + Jn > T) stop(sprintf("Bn+Jn too large: T = %s", T), call. = FALSE)
  geweke_p <- function(x) {
    z <- coda::geweke.diag(x)$z
    min((1 - pnorm(z)) * 2, (1 - pnorm(-z)) * 2)
  }
  tg1 <- resg$muv[-seq_len(Bn)][seq(s, Jn, step)]
  tg2 <- resg$lamv[-seq_len(Bn)][seq(s, Jn, step)]
  if (K > 1) {
    tg3 <- resg$rm[-seq_len(Bn), ][seq(s, Jn, step), , drop = FALSE]
    .t1 <- apply(tg3, 2, function(x) {
      EnvStats::serialCorrelationTest(x, test = "rank.von.Neumann")$p.value
    })
    .t2 <- apply(tg3, 2, function(x) {
      EnvStats::serialCorrelationTest(x, test = "AR1.yw")$p.value
    })
    .t3 <- apply(tg3, 2, function(x) nortest::ad.test(x)$p.value)
    .t4 <- apply(tg3, 2, function(x) nortest::cvm.test(x)$p.value)
    .t5 <- apply(tg3, 2, function(x) nortest::lillie.test(x)$p.value)
    .t6 <- apply(tg3, 2, geweke_p)
  } else {
    tg3 <- resg$rm[-seq_len(Bn)][seq(s, Jn, step)]
    .t1 <- EnvStats::serialCorrelationTest(tg3, test = "rank.von.Neumann")$p.value
    .t2 <- EnvStats::serialCorrelationTest(tg3, test = "AR1.yw")$p.value
    .t3 <- nortest::ad.test(tg3)$p.value
    .t4 <- nortest::cvm.test(tg3)$p.value
    .t5 <- nortest::lillie.test(tg3)$p.value
    .t6 <- geweke_p(tg3)
  }
  s1 <- c(
    EnvStats::serialCorrelationTest(tg1, test = "rank.von.Neumann")$p.value,
    EnvStats::serialCorrelationTest(tg2, test = "rank.von.Neumann")$p.value,
    .t1
  )
  s2 <- c(
    EnvStats::serialCorrelationTest(tg1, test = "AR1.yw")$p.value,
    EnvStats::serialCorrelationTest(tg2, test = "AR1.yw")$p.value,
    .t2
  )
  s3 <- c(nortest::ad.test(tg1)$p.value, nortest::ad.test(tg2)$p.value, .t3)
  s4 <- c(nortest::cvm.test(tg1)$p.value, nortest::cvm.test(tg2)$p.value, .t4)
  s5 <- c(nortest::lillie.test(tg1)$p.value, nortest::lillie.test(tg2)$p.value, .t5)
  s6 <- c(geweke_p(tg1), geweke_p(tg2), .t6)
  rn <- c(
    "SerialCorrelationTest[rank.von.Neumann]",
    "SerialCorrelationTest[AR1.yw]",
    "NormalityTest[Anderson-Darling]",
    "NormalityTest[Cramer-von Mises]",
    "NormalityTest[Lilliefors]",
    "ConvergenceTest[Geweke]"
  )
  estconv <- rbind(s1, s2, s3, s4, s5, s6)
  estconv <- data.frame(estconv, row.names = rn)
  colnames(estconv) <- c("mu", "lambda", paste0("r", 1:K))
  estconv
}
