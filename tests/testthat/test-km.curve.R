test_that("Surv yileds the same outputs as Surv of survfit", {
  library(survival)
  data(diabetes)
  expected <- survival::Surv(diabetes$time, diabetes$status)
  tested <- Surv_(diabetes$time, diabetes$status)
  expect_equal(expected, tested)
})

test_that("Surv yileds the same outputs as Surv of survfit wiht NA", {
  library(survival)
  data(diabetes)
  diabetes$time[1] <- NA
  diabetes$status[1] <- NA
  expected <- survival::Surv(diabetes$time, diabetes$status)
  tested <- Surv_(diabetes$time, diabetes$status)
  expect_equal(expected, tested)
})

test_that("Surv yileds the same outputs as Surv of survfit with a factor", {
  library(survival)
  data(diabetes)
  f <- as.factor(diabetes$status)
  expected <- survival::Surv(diabetes$time, f)
  tested <- Surv_(diabetes$time, f)

  ex <- survival::Surv(diabetes$time, diabetes$status)
  te <- Surv_(diabetes$time, diabetes$status)
  expect_equal(expected[,2], tested[,2])
})

#test_that("Surv yileds the same outputs as Surv of survfit when combined with formula", {
#  library(survival)
#  data(diabetes)
#  formula <- Surv_(diabetes$time, diabetes$status)~1
#  #cl <- match.call()
#  mf <- match.call(expand.dots = TRUE)[1:3]
#  #mf <- match.call()
#  special <- c("strata", "cluster", "offset")
#  out_terms <- terms(formula, special, data = diabetes)
#  mf$formula <- out_terms
#  mf[[1]] <- as.name("model.frame")
#  mf <- eval(mf, parent.frame())
#  Y <- model.extract(mf, "response")
#  expected <- diabetes$status
#  tested <- Y[,2]
#  expect_equal(expected, tested)
#})

test_that("km.curve yileds the same outputs as survfit", {
  library(survival)
  library(ggsurvfit)
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
  library(ggsurvfit)
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
  library(ggsurvfit)
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
  library(ggsurvfit)
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
  library(ggsurvfit)
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
  library(ggsurvfit)
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
  library(ggsurvfit)
  library(Rcpp)
  df_test <- createTestData(20, 1, first_zero=TRUE, last_zero=FALSE, subset_present=FALSE, logical_strata=FALSE, na_strata=TRUE)
  e <- survfit(Surv(t, d)~1, df_test, weight=w, na.action=na.omit, conf.type = "none")
  t <- km.curve(Surv(t, d)~1, df_test, weight="w", na.action=na.omit, conf.type = "none")
  expected <- as.numeric(c(e$n.censor))
  tested <- as.numeric(c(t$n.censor))
  expect_equal(expected, tested)
})

test_that("createAnalysisDataset produced expected variables with missing data", {
  data(diabetes.complications)
  nuisance.model <- t~fruitq1

  other.variables.analyzed <- NULL
  all_vars <- c(all.vars(nuisance.model), other.variables.analyzed)
  expected <- diabetes.complications[, all_vars, drop = FALSE]
  expected <- expected[-(1), ]
  expected <- expected$t
  diabetes.complications$t[1] <- NA
  tested <- createAnalysisDataset(formula=nuisance.model, data=diabetes.complications, other.variables.analyzed=NULL, subset.condition=NULL, na.action=na.omit)
  tested <- tested$t
  expect_equal(expected, tested)
})

test_that("createAnalysisDataset produced expected a subset dataset", {
  data(diabetes.complications)
  nuisance.model <- t~fruitq1

  other.variables.analyzed <- NULL
  all_vars <- c(all.vars(nuisance.model), other.variables.analyzed)
  expected <- diabetes.complications[, all_vars, drop = FALSE]
  subset.condition="(diabetes.complications$fruitq1 == 1)"
  expected <- subset(expected, eval(parse(text = subset.condition)))
  tested <- createAnalysisDataset(formula=nuisance.model, data=diabetes.complications, other.variables.analyzed=NULL, subset.condition="(diabetes.complications$fruitq1 == 1)", na.action=na.omit)
  expect_equal(expected, tested)
})

test_that("createAnalysisDataset produced expected a subset dataset of men", {
  data(diabetes.complications)
  nuisance.model <- t~fruitq1

  other.variables.analyzed <- "sex"
  all_vars <- c(all.vars(nuisance.model), other.variables.analyzed)
  expected <- diabetes.complications[, all_vars, drop = FALSE]
  subset.condition="(diabetes.complications$sex == 1)"
  expected <- subset(expected, eval(parse(text = subset.condition)))
  tested <- createAnalysisDataset(formula=nuisance.model, data=diabetes.complications, other.variables.analyzed="sex", subset.condition="(diabetes.complications$sex == 1)", na.action=na.omit)
  expect_equal(expected, tested)
})
