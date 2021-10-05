context("One sample t-test")
library(pARI)

#----    input checks    ----


X <- matrix(rnorm(8*10), nrow = 8, byrow = T)


test_that("evaluate the correct one sample t-test", {
  
  #Two.sided
  out <- oneSample(X)
  Test <- round(out$Test,6)
  Test_test_two <- sapply(c(1:nrow(X)), function(x) round(t.test(X[x,])$statistic,6))
  names(Test_test_two) <- NULL
  P_test_two <- sapply(c(1:nrow(X)), function(x) round(t.test(X[x,])$p.value,6))
  names(P_test_two) <- NULL
  
  #Greater
  out <- oneSample(X, alternative = "greater")
  Test <- round(out$Test,6)
  Test_test_g <- sapply(c(1:nrow(X)), function(x) round(t.test(X[x,], alternative = "greater")$statistic,6))
  names(Test_test_g) <- NULL
  P_test_g <- sapply(c(1:nrow(X)), function(x) round(t.test(X[x,], alternative = "greater")$p.value,6))
  names(P_test_g) <- NULL
  
  #Lower
  out <- oneSample(X, alternative = "lower")
  Test <- round(out$Test,6)
  Test_test_l <- sapply(c(1:nrow(X)), function(x) round(t.test(X[x,], alternative = "less")$statistic,6))
  names(Test_test_l) <- NULL
  P_test_l <- sapply(c(1:nrow(X)), function(x) round(t.test(X[x,], alternative = "less")$p.value,6))
  names(P_test_l) <- NULL
  
  expect_equal(round(oneSample(X = X, alternative = "two.sided")$Test,6),
               Test_test_two)
  expect_equal(round(oneSample(X = X, alternative = "greater")$Test,6),
               Test_test_g)
  expect_equal(round(oneSample(X = X, alternative = "lower")$Test,6),
               Test_test_l)

})

test_that("evaluate the correct one sample t-test", {
  
  #Two.sided
  out <- oneSample(X)
  Test <- round(out$Test,6)
  Test_test_two <- sapply(c(1:nrow(X)), function(x) round(t.test(X[x,])$statistic,6))
  names(Test_test_two) <- NULL
  P_test_two <- sapply(c(1:nrow(X)), function(x) round(t.test(X[x,])$p.value,6))
  names(P_test_two) <- NULL
  
  #Greater
  out <- oneSample(X, alternative = "greater")
  Test <- round(out$Test,6)
  Test_test_g <- sapply(c(1:nrow(X)), function(x) round(t.test(X[x,], alternative = "greater")$statistic,6))
  names(Test_test_g) <- NULL
  P_test_g <- sapply(c(1:nrow(X)), function(x) round(t.test(X[x,], alternative = "greater")$p.value,6))
  names(P_test_g) <- NULL
  
  #Lower
  out <- oneSample(X, alternative = "lower")
  Test <- round(out$Test,6)
  Test_test_l <- sapply(c(1:nrow(X)), function(x) round(t.test(X[x,], alternative = "less")$statistic,6))
  names(Test_test_l) <- NULL
  P_test_l <- sapply(c(1:nrow(X)), function(x) round(t.test(X[x,], alternative = "less")$p.value,6))
  names(P_test_l) <- NULL
  
  expect_equal(round(oneSample(X = X, alternative = "two.sided")$pv,6),
               P_test_two)
  expect_equal(round(oneSample(X = X, alternative = "greater")$pv,6),
               P_test_g)
  expect_equal(round(oneSample(X = X, alternative = "lower")$pv,6),
               P_test_l)
  
})