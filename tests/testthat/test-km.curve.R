
test_that("km.curve yileds the same outputs as survfit", {
  library(survival)
  library(Rcpp)
  df_test <- createTestData(20, 2, first_zero=TRUE, last_zero=FALSE, subset_present=FALSE, logical_strata=FALSE, na_strata=FALSE)
  e <- survfit(Surv(t, d)~1, df_test, weight=w, na.action=na.omit, conf.type = "none")
  t <- km.curve(Surv(t, d)~1, df_test, weight="w", na.action=na.omit, conf.type = "none")
  expected <- as.numeric(c(e$time, round(e$surv,digit=5), e$n, e$n.risk, e$n.event, e$n.censor, round(e$std.err,digit=5), e$strata))
  tested <- as.numeric(c(t$time, round(t$surv,digit=5), t$n, t$n.risk, t$n.event, t$n.censor, round(t$std.err,digit=5), t$strata))
  expect_equal(expected, tested)
})

test_that("km.curve yileds the same outputs as survfit", {
  library(survival)
  library(Rcpp)
  df_test <- createTestData(20, 1, first_zero=FALSE, last_zero=TRUE, subset_present=FALSE, logical_strata=TRUE, na_strata=FALSE)
  e <- survfit(Surv(t, d)~strata, df_test, weight=w, conf.type = "none")
  t <- km.curve(Surv(t, d)~strata, df_test, weight="w", conf.type = "none")
  expected <- as.numeric(c(e$time, round(e$surv,digit=5), e$n, e$n.risk, e$n.event, e$n.censor, round(e$std.err,digit=5), e$strata))
  tested <- as.numeric(c(t$time, round(t$surv,digit=5), t$n, t$n.risk, t$n.event, t$n.censor, round(t$std.err,digit=5), t$strata))
  expect_equal(expected, tested)
})

test_that("km.curve yileds the same outputs as survfit", {
  library(survival)
  library(Rcpp)
  df_test <- createTestData(200, 2, first_zero=TRUE, last_zero=TRUE, subset_present=FALSE, logical_strata=FALSE, na_strata=FALSE)
  e <- survfit(Surv(t, d)~strata, df_test, weight=w, conf.type = "log-log")
  t <- km.curve(Surv(t, d)~strata, df_test, weight="w", conf.type = "log-log")
  e$lower <- sapply(e$lower, function(x) ifelse(is.nan(x), NA, x))
  e$upper <- sapply(e$upper, function(x) ifelse(is.nan(x), NA, x))
  e$lower <- sapply(e$lower, function(x) ifelse(x>=1, NA, x))
  e$upper <- sapply(e$upper, function(x) ifelse(x>=1, NA, x))
  e$lower <- sapply(e$lower, function(x) ifelse(x<=0, NA, x))
  e$upper <- sapply(e$upper, function(x) ifelse(x<=0, NA, x))
  t$lower <- sapply(t$lower, function(x) ifelse(is.nan(x), NA, x))
  t$upper <- sapply(t$upper, function(x) ifelse(is.nan(x), NA, x))
  t$lower <- sapply(t$lower, function(x) ifelse(x>=1, NA, x))
  t$upper <- sapply(t$upper, function(x) ifelse(x>=1, NA, x))
  t$lower <- sapply(t$lower, function(x) ifelse(x<=0, NA, x))
  t$upper <- sapply(t$upper, function(x) ifelse(x<=0, NA, x))
  expected <- as.numeric(c(e$time, round(e$surv,digit=5), e$n, e$n.risk, e$n.event, e$n.censor, round(e$std.err,digit=5), round(e$lower,digit=5), round(e$upper,digit=5), e$strata))
  tested <- as.numeric(c(t$time, round(t$surv,digit=5), t$n, t$n.risk, t$n.event, t$n.censor, round(t$std.err,digit=5), round(t$lower,digit=5), round(t$upper,digit=5), t$strata))
  expect_equal(expected, tested)
})

test_that("km.curve yileds the same outputs as survfit", {
  library(survival)
  library(Rcpp)
  df_test <- createTestData(200, 2, first_zero=TRUE, last_zero=TRUE, subset_present=FALSE, logical_strata=FALSE, na_strata=FALSE)
  e <- survfit(Surv(t, d)~strata, df_test, weight=w, conf.type = "log")
  t <- km.curve(Surv(t, d)~strata, df_test, weight="w", conf.type = "log")
  e$lower <- sapply(e$lower, function(x) ifelse(is.nan(x), NA, x))
  e$upper <- sapply(e$upper, function(x) ifelse(is.nan(x), NA, x))
  e$lower <- sapply(e$lower, function(x) ifelse(x>=1, NA, x))
  e$upper <- sapply(e$upper, function(x) ifelse(x>=1, NA, x))
  e$lower <- sapply(e$lower, function(x) ifelse(x<=0, NA, x))
  e$upper <- sapply(e$upper, function(x) ifelse(x<=0, NA, x))
  t$lower <- sapply(t$lower, function(x) ifelse(is.nan(x), NA, x))
  t$upper <- sapply(t$upper, function(x) ifelse(is.nan(x), NA, x))
  t$lower <- sapply(t$lower, function(x) ifelse(x>=1, NA, x))
  t$upper <- sapply(t$upper, function(x) ifelse(x>=1, NA, x))
  t$lower <- sapply(t$lower, function(x) ifelse(x<=0, NA, x))
  t$upper <- sapply(t$upper, function(x) ifelse(x<=0, NA, x))

  expected <- as.numeric(c(e$time, round(e$surv,digit=5), e$n, e$n.risk, e$n.event, e$n.censor, round(e$std.err,digit=5), round(e$lower,digit=5), round(e$upper,digit=5), e$strata))
  tested <- as.numeric(c(t$time, round(t$surv,digit=5), t$n, t$n.risk, t$n.event, t$n.censor, round(t$std.err,digit=5), round(t$lower,digit=5), round(t$upper,digit=5), t$strata))
  expect_equal(expected, tested)
})

test_that("km.curve yileds the same outputs as survfit", {
  library(survival)
  library(Rcpp)
  df_test <- createTestData(200, 2, first_zero=TRUE, last_zero=TRUE, subset_present=FALSE, logical_strata=FALSE, na_strata=FALSE)
  e <- survfit(Surv(t, d)~strata, df_test, weight=w, conf.type = "a")
  t <- km.curve(Surv(t, d)~strata, df_test, weight="w", conf.type = "a")
  e$lower <- sapply(e$lower, function(x) ifelse(is.nan(x), NA, x))
  e$upper <- sapply(e$upper, function(x) ifelse(is.nan(x), NA, x))
  e$lower <- sapply(e$lower, function(x) ifelse(x>=1, NA, x))
  e$upper <- sapply(e$upper, function(x) ifelse(x>=1, NA, x))
  e$lower <- sapply(e$lower, function(x) ifelse(x<=0, NA, x))
  e$upper <- sapply(e$upper, function(x) ifelse(x<=0, NA, x))
  t$lower <- sapply(t$lower, function(x) ifelse(is.nan(x), NA, x))
  t$upper <- sapply(t$upper, function(x) ifelse(is.nan(x), NA, x))
  t$lower <- sapply(t$lower, function(x) ifelse(x>=1, NA, x))
  t$upper <- sapply(t$upper, function(x) ifelse(x>=1, NA, x))
  t$lower <- sapply(t$lower, function(x) ifelse(x<=0, NA, x))
  t$upper <- sapply(t$upper, function(x) ifelse(x<=0, NA, x))

  expected <- as.numeric(c(e$time, round(e$surv,digit=5), e$n, e$n.risk, e$n.event, e$n.censor, round(e$std.err,digit=5), round(e$lower,digit=5), round(e$upper,digit=5), e$strata))
  tested <- as.numeric(c(t$time, round(t$surv,digit=5), t$n, t$n.risk, t$n.event, t$n.censor, round(t$std.err,digit=5), round(t$lower,digit=5), round(t$upper,digit=5), t$strata))
  expect_equal(expected, tested)
})

test_that("km.curve yileds the same outputs as survfit", {
  library(survival)
  library(Rcpp)
  df_test <- createTestData(200, 2, first_zero=TRUE, last_zero=TRUE, subset_present=FALSE, logical_strata=FALSE, na_strata=FALSE)
  e <- survfit(Surv(t, d)~strata, df_test, weight=w, conf.type = "plain")
  t <- km.curve(Surv(t, d)~strata, df_test, weight="w", conf.type = "plain")
  e$lower <- sapply(e$lower, function(x) ifelse(is.nan(x), NA, x))
  e$upper <- sapply(e$upper, function(x) ifelse(is.nan(x), NA, x))
  e$lower <- sapply(e$lower, function(x) ifelse(x>=1, NA, x))
  e$upper <- sapply(e$upper, function(x) ifelse(x>=1, NA, x))
  e$lower <- sapply(e$lower, function(x) ifelse(x<=0, NA, x))
  e$upper <- sapply(e$upper, function(x) ifelse(x<=0, NA, x))
  t$lower <- sapply(t$lower, function(x) ifelse(is.nan(x), NA, x))
  t$upper <- sapply(t$upper, function(x) ifelse(is.nan(x), NA, x))
  t$lower <- sapply(t$lower, function(x) ifelse(x>=1, NA, x))
  t$upper <- sapply(t$upper, function(x) ifelse(x>=1, NA, x))
  t$lower <- sapply(t$lower, function(x) ifelse(x<=0, NA, x))
  t$upper <- sapply(t$upper, function(x) ifelse(x<=0, NA, x))

  expected <- as.numeric(c(e$time, round(e$surv,digit=5), e$n, e$n.risk, e$n.event, e$n.censor, round(e$std.err,digit=5), round(e$lower,digit=5), round(e$upper,digit=5), e$strata))
  tested <- as.numeric(c(t$time, round(t$surv,digit=5), t$n, t$n.risk, t$n.event, t$n.censor, round(t$std.err,digit=5), round(t$lower,digit=5), round(t$upper,digit=5), t$strata))
  expect_equal(expected, tested)
})

# Error due to exclusion when na_strata=TRUE and strata not used
test_that("km.curve yileds the same outputs as survfit", {
  library(survival)
  library(Rcpp)
  df_test <- createTestData(20, 1, first_zero=TRUE, last_zero=FALSE, subset_present=FALSE, logical_strata=FALSE, na_strata=TRUE)
  e <- survfit(Surv(t, d)~1, df_test, weight=w, na.action=na.omit, conf.type = "none")
  t <- km.curve(Surv(t, d)~1, df_test, weight="w", na.action=na.omit, conf.type = "none")
  expected <- as.numeric(c(e$n.censor))
  tested <- as.numeric(c(t$n.censor))
  expect_equal(expected, tested)
})
