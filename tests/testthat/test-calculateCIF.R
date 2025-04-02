test_that("calclateCIF works with a toy data", {
  n_test <- 10
  t_test <- 1:n_test
  epsilon_test <- rep(1, n_test)
  epsilon_test[9] <- 0
  epsilon_test[10] <- 2
  df_test <- data.frame(id = 1:n_test, t_test = t_test, epsilon_test = epsilon_test)
  tested <- calculateCIF2(t_test, epsilon_test)
  out <- cif(Event(t_test,epsilon_test) ~ +1, data=df_test, cause=1)
  expected <- out$mu
  expect_equal(expected, tested)
})


test_that("calclateCIF works with a toy data", {
  n_test <- 10
  t_test <- 1:n_test
  epsilon_test <- rep(1, n_test)
  epsilon_test[10] <- 0
  epsilon_test[9] <- 2
  df_test <- data.frame(id = 1:n_test, t_test = t_test, epsilon_test = epsilon_test)
  tested <- calculateCIF2(t_test, epsilon_test)
  out <- cif(Event(t_test,epsilon_test) ~ +1, data=df_test, cause=1)
  expected <- out$mu
  expect_equal(expected, tested)
})

test_that("calclateCIF works with a toy data", {
  n_test <- 10
  t_test <- 1:n_test
  epsilon_test <- rep(1, n_test)
  epsilon_test[1] <- 0
  epsilon_test[2] <- 2
  df_test <- data.frame(id = 1:n_test, t_test = t_test, epsilon_test = epsilon_test)
  tested <- calculateCIF2(t_test, epsilon_test)
  out <- cif(Event(t_test,epsilon_test) ~ +1, data=df_test, cause=1)
  expected <- out$mu
  expect_equal(expected, tested)
})

test_that("calclateCIF works with a toy data", {
  n_test <- 10
  t_test <- 1:n_test
  epsilon_test <- rep(1, n_test)
  epsilon_test[2] <- 0
  epsilon_test[1] <- 2
  df_test <- data.frame(id = 1:n_test, t_test = t_test, epsilon_test = epsilon_test)
  tested <- calculateCIF2(t_test, epsilon_test)
  out <- cif(Event(t_test,epsilon_test) ~ +1, data=df_test, cause=1)
  expected <- out$mu
  expect_equal(expected, tested)
})

test_that("calclateCIF works with a stratified data", {
  n_test <- 10
  t_test <- 1:n_test
  strata_test <- (t_test>5)
  epsilon_test <- rep(1, n_test)
  epsilon_test[10] <- 0
  epsilon_test[9] <- 2
  df_test <- data.frame(id = 1:n_test, t_test = t_test, epsilon_test = epsilon_test, strata_test = strata_test)
  tested <- calculateCIF4(t_test, epsilon_test, strata_test)
  out <- cif(Event(t_test,epsilon_test) ~ strata(strata_test), data=df_test, cause=1)
  expected <- out$mu
  expect_equal(expected, tested)
})
