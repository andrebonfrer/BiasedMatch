test_that("PROBM handles the no-overlap case", {
  expect_equal(PROBM(m = 0, N1 = 10, N2 = 10, N = 0, m1 = 3, m2 = 3, w = 1), 1)
  expect_equal(PROBM(m = 1, N1 = 10, N2 = 10, N = 0, m1 = 3, m2 = 3, w = 1), 0)
})

test_that("PROBM matches paper values and is zero outside the support", {
  expect_equal(
    round(PROBM(m = 5, N1 = 37, N2 = 45, N = 16, m1 = 12, m2 = 8, w = 1), 4),
    0.0003
  )
  expect_equal(PROBM(m = 20, N1 = 37, N2 = 45, N = 16, m1 = 12, m2 = 8, w = 1), 0)
})

test_that("DISTM reproduces Table 1 of Puza and Bonfrer (2018)", {
  # w = 1 (Table 1a)
  d1 <- DISTM(N1 = 37, N2 = 45, N = 16, m1 = 12, m2 = 8, w = 1)
  expect_equal(
    round(d1$fmvec[1:6], 4),
    c(0.3573, 0.4137, 0.1836, 0.0404, 0.0047, 0.0003)
  )
  expect_equal(round(d1$Em, 4), 0.9225)

  # w = 2 (Table 1b)
  d2 <- DISTM(N1 = 37, N2 = 45, N = 16, m1 = 12, m2 = 8, w = 2)
  expect_equal(
    round(d2$fmvec[1:7], 4),
    c(0.1182, 0.3164, 0.3309, 0.1750, 0.0507, 0.0081, 0.0007)
  )
  expect_equal(round(d2$Em, 4), 1.7511)
})

test_that("DISTM probabilities sum to one", {
  expect_equal(DISTM(N1 = 10, N2 = 10, N = 3, m1 = 3, m2 = 3, w = 2)$sumfmvec, 1,
    tolerance = 1e-8
  )
  expect_equal(DISTM(N1 = 292, N2 = 350, N = 192, m1 = 90, m2 = 59, w = 1.5)$sumfmvec,
    1,
    tolerance = 1e-8
  )
})

test_that("Vmf matches the exact variance from DISTM when w = 1", {
  d <- DISTM(N1 = 37, N2 = 45, N = 16, m1 = 12, m2 = 8, w = 1)
  expect_equal(Vmf(N1 = 37, N2 = 45, N = 16, m1 = 12, m2 = 8), d$Vm,
    tolerance = 1e-8
  )
})
