test_that("estimateCIFcpp works with a toy data", {
  n_test <- 10
  t_test <- 1:n_test
  epsilon_test <- rep(1, n_test)
  epsilon_test[9] <- 0
  epsilon_test[10] <- 2
  df_test <- data.frame(id = 1:n_test, t_test = t_test, epsilon_test = epsilon_test)
  tested <- estimateCIFcpp(t = df_test$t_test,
                           epsilon = df_test$epsilon_test,
                           code_event1 = 1,
                           code_event2 = 2,
                           code_censoring = 0)
  expected <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.8, 0.8)
  expect_equal(expected, tested$CIF1)
})

