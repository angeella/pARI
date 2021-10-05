context("Simulate Data")

m <- 200
power <- 0.8
rho <- 0
n <- 50
pi0 <- 0.8
alpha <- 0.05

test_that("output dimensions are as expected for two-sample tests", {
  out <- simulateData(pi0 = pi0,m = m,n = n, rho = rho, set.seed = 1234, power = power, alpha = alpha)
  
  expect_equal(dim(out)[1], m)
  expect_equal(dim(out)[2], n)

})

test_that("output dimensions are as expected for one-sample tests", {
  sim <- gaussianSamples(m, rho, n, pi0, SNR = 1, prob = 1)
  
  expect_equal(ncol(sim$X), n)
  expect_equal(nrow(sim$X), m)
  expect_length(sim$H, m)
  
  tbl <- table(sim$categ)
  expect_length(tbl, 1)
  
  tbl <- table(sim$H)
  expect_length(tbl, 2)
})


