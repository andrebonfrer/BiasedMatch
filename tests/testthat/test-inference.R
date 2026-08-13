test_that("inferW reproduces Problem 3 (marketing data)", {
  res <- inferW(m = 13, N1 = 292, N2 = 350, N = 192, m1 = 90, m2 = 59)
  expect_equal(round(res$mle, 4), 1.5265)
  expect_equal(round(res$posterior_mode, 4), 1.3029)
  expect_equal(round(res$p_value, 4), 0.1737)
})

test_that("inferW reproduces Problem 2 (committee gender bias, B = 2)", {
  res <- inferW(m = 2, N1 = 52, N2 = 52, N = 22, m1 = 9, m2 = 6)
  expect_equal(round(res$mle, 2), 16.99)
  expect_equal(round(res$posterior_mode, 2), 2.29)
  expect_equal(round(res$p_value, 3), 0.061)
})

test_that("inferW reports an infinite MLE for large B", {
  # Problem 2 with B = 4: MLE is infinite, posterior mode finite (Table 2)
  res <- inferW(m = 4, N1 = 52, N2 = 52, N = 22, m1 = 9, m2 = 6)
  expect_equal(res$mle, Inf)
  expect_equal(round(res$posterior_mode, 2), 9.11)
  expect_equal(round(res$p_value, 6), 0.000177)
})

test_that("likW returns standardised curves on the supplied grid", {
  wgrid <- seq(0.1, 4, by = 0.1)
  curves <- likW(m = 13, N1 = 292, N2 = 350, N = 192, m1 = 90, m2 = 59, wgrid = wgrid)
  expect_named(curves, c("w", "likelihood", "posterior"))
  expect_equal(curves$w, wgrid)
  expect_equal(max(curves$likelihood), 1)
  expect_equal(max(curves$posterior), 1)
})

test_that("likW rejects non-positive weights", {
  expect_snapshot(
    error = TRUE,
    likW(m = 1, N1 = 10, N2 = 10, N = 3, m1 = 3, m2 = 3, wgrid = c(-1, 1))
  )
})
