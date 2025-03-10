<<<<<<< HEAD
test_that("polyreg produced expected coefficients and var for retinopathy dataset", {
  data(retinopathy)
  model1 <- "Event(t,epsilon) ~ sex"
  model1 <- as.formula(model1)
  model2 <- "Event(t,epsilon) ~ age+sex+bmi+hba1c+drug_oha+drug_insulin"
  model2 <- as.formula(model2)
  model3 <- "Event(t,epsilon) ~ age+sex+bmi+hba1c+diabetes_duration+drug_oha+drug_insulin+sbp+ldl+hdl+tg+current_smoker+alcohol_drinker+ltpa"
  model3 <- as.formula(model3)
  result <- polyreg(nuisance.model = model1, exposure = 'fruitq1',
                    cens.model = Event(t,epsilon==0)~strata(strata), data = retinopathy, outcome.type = 'C',                  effect.measure1='RR', effect.measure2='RR', time.point=8)
  tested_coefficient <- round(result$out_coefficient,digit=5)
  tested_cov <- round(result$out_cov[1,],digit=5)
  tested <- as.vector(cbind(tested_coefficient,tested_cov))
  expected <- c(-1.58510,  0.46386,  0.33934, -3.53254, -1.25689, -0.01445, 0.01379, -0.01142, -0.00588, 0.00565, -0.00482, -0.00309)
  expect_equal(expected, tested)
})

test_that("polyreg produced expected coefficients and var for bmt dataset", {
  data(bmt)
  result <- polyreg(nuisance.model = Event(time, cause)~age+tcell, exposure = 'platelet', cens.model = Event(time,cause==0)~+1, data = bmt, effect.measure1='RR', effect.measure2='RR', time.point=24, outcome.type='COMPETINGRISK')
  tested_coefficient <- round(result$out_coefficient,digit=2)
  tested_cov <- round(result$out_cov[1,],digit=2)
  tested <- as.vector(cbind(tested_coefficient,tested_cov))
  expected <- c(-0.45, 1.01, -1.32, -0.48, -1.68, 0.41, 0.63, 0.12, 0.02, 0.00, -0.01, -0.01, 0.02, 0.00, -0.01, -0.01)
  #  expected <- c(-0.45, 1.01, -1.32, -0.48, -1.68, 0.41, 0.63, 0.12, 0.02, 0.00, -0.01, -0.01, 0.02, 0.01, -0.01, -0.01)
  expect_equal(expected, tested)
})

#test_that("predict.polyreg produced expected failure probabilities of first 2 obs bmt dataset", {
#  data(bmt)
#  bmt$cause1 <- as.numeric((bmt$cause>0))
#  result <- polyreg(nuisance.model = Event(time, cause1)~age+tcell, exposure = 'platelet', cens.model = Event(time,cause==0)~+1, data = bmt, effect.measure1='RR', effect.measure2='RR', time.point=24, outcome.type='SURVIVAL')
#  tmp1 <- predict.polyreg(formula = Event(time, cause)~age+tcell, exposure = 'platelet', data = bmt, coefficient = result$out_coefficient, effect.measure1='RR', effect.measure2='RR', outcome.type='SURVIVAL')
#  tmp2 <- round(tmp1[1,],digit=5)
#  tmp3 <- round(tmp1[2,],digit=5)
#  tested <- as.vector(cbind(tmp2,tmp3))
#  expected <- c(0.64673, 0.50754, 0.68552, 0.53798)
#  expect_equal(expected, tested)
#})

#test_that("predict.polyreg produced expected of cumulative incidence probabilities of first 2 obs for diabetes dataset", {
#data(diabetes)
#model1 <- "Event(t,epsilon) ~ sex"
#model1 <- as.formula(model1)
#model2 <- "Event(t,epsilon) ~ age+sex+bmi+hba1c+drug_oha+drug_insulin"
#model2 <- as.formula(model2)
#model3 <- "Event(t,epsilon) ~ age+sex+bmi+hba1c+diabetes_duration+drug_oha+drug_insulin+sbp+ldl+hdl+tg+current_smoker+alcohol_drinker+ltpa"
#model3 <- as.formula(model3)
#result <- polyreg(nuisance.model = model3, exposure = 'fruitq1',
#                  cens.model = Event(t,epsilon==0)~strata(strata), data = jdcs, outcome.type = 'C',
#                  effect.measure1='RR', effect.measure2='RR', time.point=8)
#}

#test_that("predict.polyreg produced expected of cumulative incidence probabilities of first 2 obs for bmt dataset", {
#  data(bmt)
#  result <- polyreg(nuisance.model = Event(time, cause)~age+tcell, exposure = 'platelet', cens.model = Event(time,cause==0)~+1, data = bmt, effect.measure1='RR', effect.measure2='RR', time.point=24, outcome.type='COMPETINGRISK')
#  tmp1 <- predict.polyreg(formula = Event(time, cause)~age+tcell, exposure = 'platelet', data = bmt, coefficient = result$out_coefficient, effect.measure1='RR', effect.measure2='RR', outcome.type='COMPETINGRISK')
#  tmp2 <- round(tmp1[1,],digit=5)
#  tmp3 <- round(tmp1[2,],digit=5)
#  tested <- as.vector(cbind(tmp2,tmp3))
#  expected <- c(0.47241, 0.17850, 0.29166, 0.20126, 0.52685, 0.17510, 0.32527, 0.19743)
#  expected <- c(0.47241, 0.29166, 0.17850, 0.20126, 0.52685, 0.32527, 0.17510, 0.19743)
=======
test_that("calculateIPCW produced expected IP weights in diabetes.complications", {
  data(diabetes.complications)
  output <- calculateIPCW(formula=Event(t,epsilon)~+1, data=diabetes.complications, code.censoring=0, strata_name='strata', specific.time=8)
  tested <- round(output[1:6],digit=3)
  expected <- c(1.656, 1.656, 0, 1.656, 1.656, 1.004)
  expect_equal(expected, tested)
})

test_that("polyreg produced expected coefficients and variance covariance matrix from competing risks data in diabetes.complications", {
  data(diabetes.complications)
  output <- polyreg(nuisance.model = Event(t,epsilon)~+1, exposure = 'fruitq1', data = diabetes.complications, effect.measure1='RR', effect.measure2='RR', time.point=8, outcome.type='C')
  tested_coefficient <- round(output$coefficient,digit=3)
  tested_cov <- round(output$cov[1,],digit=3)
  tested <- as.vector(cbind(tested_coefficient,tested_cov))
  expected <- c(-1.383, 0.300, -3.991, 0.076, 0.007, -0.005, -0.001, 0.005)
  expect_equal(expected, tested)
})

#test_that("polyreg produced expected coefficients and variance covariance matrix from competing risks data in diabetes.complications", {
#  library(dplyr)
#  library(boot)
#  data(prostate)
#  prostate <- prostate %>% mutate(epsilon=ifelse(status=="alive",0,
#                                                ifelse(status=="dead - prostatic ca",1,
#                                                ifelse(status=="dead - other ca",1,
#                                                ifelse(status=="dead - heart or vascular",2,
#                                                ifelse(status=="dead - cerebrovascular",2,2)
#                                  )))))
#  prostate$epsilon <- as.numeric(prostate$epsilon)
#  prostate$a <- as.numeric((prostate$rx=="placebo"))
#  prostate$t <- prostate$dtime/12
#  output <- polyreg(nuisance.model = Event(t,epsilon) ~ +1, exposure = 'a', strata='stage', data = prostate,
#                    effect.measure1='RR', effect.measure2='RR', time.point=1:5, outcome.type='PROPORTIONAL', boot.R=20)
#  tested_coefficient <- round(output$coefficient,digit=3)
#  tested_cov <- round(output$cov[1,],digit=3)
#  tested <- as.vector(cbind(tested_coefficient,tested_cov))
#  expected <- c(-1.383, 0.300, -3.991, 0.076, 0.007, -0.005, -0.001, 0.005)
>>>>>>> 26f9e354ed9810cbf9faaf4dfdcfe203945e81d3
#  expect_equal(expected, tested)
#})


<<<<<<< HEAD

=======
test_that("polyreg produced expected coefficients and variance covariance matrix when stratified IPCW used", {
  data(diabetes.complications)
  output <- polyreg(nuisance.model = Event(t,epsilon)~+1, exposure = 'fruitq1', data = diabetes.complications, strata = 'strata', effect.measure1='RR', effect.measure2='RR', time.point=8, outcome.type='C')
  tested_coefficient <- round(output$coefficient,digit=3)
  tested_cov <- round(output$cov[1,],digit=3)
  tested <- as.vector(cbind(tested_coefficient,tested_cov))
  #  expected <- c(-1.383, 0.300, -3.991, 0.076, 0.007, -0.005, -0.001, 0.005)
  expected <- c(-1.383, 0.300, -3.988, 0.078, 0.007, -0.005, -0.001, 0.005)
  expect_equal(expected, tested)
})

test_that("polyreg produced expected coefficients and variance covariance matrix under coding other than the default", {
  data(diabetes.complications)
  diabetes.complications$epsilon1 <- diabetes.complications$epsilon + 1
  output <- polyreg(nuisance.model = Event(t,epsilon1)~+1, exposure = 'fruitq1', data = diabetes.complications,
                    code.event1=2, code.event2=3, code.censoring=1, code.exposure.ref = 1,
                    effect.measure1='RR', effect.measure2='RR', time.point=8, outcome.type='C')
  tested_coefficient <- round(output$coefficient,digit=3)
  tested_cov <- round(output$cov[1,],digit=3)
  tested <- as.vector(cbind(tested_coefficient,tested_cov))
  expected <- c(-1.383, -0.300, -3.991, -0.076, 0.016, -0.012, 0.004, -0.004)
  expect_equal(expected, tested)
})

test_that("polyreg produced expected coefficients and variance covariance matrix from survival data in diabetes.complications", {
  data(diabetes.complications)
  diabetes.complications$d <- as.numeric(diabetes.complications$epsilon>0)
  output <- polyreg(nuisance.model = Event(t,d)~+1, exposure = 'fruitq1', data = diabetes.complications, effect.measure1='RR', effect.measure2='RR', time.point=8, outcome.type='S')
  tested_coefficient <- round(output$coefficient,digit=3)
  tested_cov <- round(output$cov[1,],digit=3)
  tested <- as.vector(cbind(tested_coefficient,tested_cov))
  expected <- c(-0.901, 0.252, 0.003, -0.003)
  expect_equal(expected, tested)
})

test_that("polyreg produced expected coefficients and variance covariance matrix when binomial regression is applied to diabetes.complications", {
  data(diabetes.complications)
  diabetes.complications$d <- as.numeric(diabetes.complications$epsilon>0)
  output <- polyreg(nuisance.model = d~+1, exposure = 'fruitq1', data = diabetes.complications, effect.measure1='RR', effect.measure2='RR', time.point=8, outcome.type='B')
  tested_coefficient <- round(output$coefficient,digit=3)
  tested_cov <- round(output$cov[1,],digit=3)
  tested <- as.vector(cbind(tested_coefficient,tested_cov))
  expected <- c(-0.877, 0.249, 0.003, -0.003)
  expect_equal(expected, tested)
})

#test_that("polyreg produced expected coefficients and variance covariance matrix from bmt dataset", {
#  library(mets)
#  library(nleqslv)
#  data(bmt)
#  output <- polyreg(nuisance.model = Event(time, cause)~age+tcell, exposure = 'platelet', cens.model = Event(time,cause==0)~+1, data = bmt, effect.measure1='RR', effect.measure2='RR', time.point=24, outcome.type='COMPETINGRISK')
#  tested_coefficient <- round(output$coefficient,digit=2)
#  tested_cov <- round(output$cov[1,],digit=2)
#  tested <- as.vector(cbind(tested_coefficient,tested_cov))
#  expected <- c(-0.45, 1.01, -1.32, -0.48, -1.68, 0.41, 0.63, 0.12, 0.02, 0.01, -0.02, -0.01, -0.00, 0.05, -0.02, 0.01)
#  expect_equal(expected, tested)
#})
>>>>>>> 26f9e354ed9810cbf9faaf4dfdcfe203945e81d3
