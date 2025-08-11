# Run tests

# tests/testthat/test-probability.R
test_that("PROBM handles no overlap correctly", {
  # when there are no common items, only m=0 has probability 1
  expect_equal(PROBM(m = 0, N1 = 10, N2 = 10, N = 0, m1 = 3, m2 = 3, w = 1), 1)
  expect_equal(PROBM(m = 1, N1 = 10, N2 = 10, N = 0, m1 = 3, m2 = 3, w = 1), 0)
})

test_that("Distribution sums to 1", {
  dist <- DISTM(N1 = 10, N2 = 10, N = 3, m1 = 3, m2 = 3, w = 2)
  expect_equal(sum(dist$fmvec), 1, tolerance = 1e-6)
})

test_that("Link functions are inverses", {
  val <- BiasedMatch:::gFUN(2.7)  # gFUN and gINV are internal, so use :::
  expect_equal(BiasedMatch:::gINV(val), 2.7, tolerance = 1e-6)
})
