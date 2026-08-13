test_that("link function gFUN and gINV are inverses", {
  val <- BiasedMatch:::gFUN(2.7)
  expect_equal(BiasedMatch:::gINV(val), 2.7, tolerance = 1e-6)
})

test_that("MHALG returns draws with the expected shape", {
  set.seed(3182)
  res <- MHALG(
    B = 3, J = 10,
    N1vec = 37, N2vec = 45, N0vec = 16,
    M1vec = 12, M2vec = 8, mvec = 3, Cvec = 1, del = 0.1
  )
  expect_named(res, c("muv", "lamv", "rm", "rar"))
  expect_length(res$muv, 14)
  expect_equal(dim(res$rm), c(14, 1))
  expect_true(res$rar >= 0 && res$rar <= 1)
})

test_that("restable summarises MHALG output", {
  set.seed(9930)
  res <- MHALG(
    B = 3, J = 30,
    N1vec = 37, N2vec = 45, N0vec = 16,
    M1vec = 12, M2vec = 8, mvec = 3, Cvec = 1, del = 0.1
  )
  tab <- restable(res, Bn = 3)
  expect_equal(rownames(tab), c("Mean", "Mode", "Q025", "Q975"))
  expect_true(all(c("mu", "lambda", "sigma", "r1") %in% colnames(tab)))
})
