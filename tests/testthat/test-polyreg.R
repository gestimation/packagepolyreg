test_that("polyreg produced expected coefficients and var for bmt dataset", {
  data(bmt)
  result <- polyreg(nuisance.model = Event(time, cause)~age+tcell, exposure = 'platelet', cens.model = Event(time,cause==0)~+1, data = bmt, effect.measure1='RR', effect.measure2='RR', time.point=24, outcome.type='COMPETINGRISK')
  tested_coefficient <- round(result$out_coefficient,digit=2)
  tested_cov <- round(result$out_cov[1,],digit=2)
  tested <- as.vector(cbind(tested_coefficient,tested_cov))
  expected <- c(-2.81, 1.01, -1.32, -0.48, -2.65, 0.41, 0.63, 0.12, 0.09, -0.03, 0.00, -0.02, 0.05, -0.02, 0.00, -0.01)
  expect_equal(expected, tested)
})

test_that("predict.polyreg produced expected failure probabilities of first 2 obs bmt dataset", {
  data(bmt)
  bmt$cause1 <- as.numeric((bmt$cause>0))
  result <- polyreg(nuisance.model = Event(time, cause1)~age+tcell, exposure = 'platelet', cens.model = Event(time,cause==0)~+1, data = bmt, effect.measure1='RR', effect.measure2='RR', time.point=24, outcome.type='SURVIVAL')
  tmp1 <- predict.polyreg(formula = Event(time, cause)~age+tcell, exposure = 'platelet', data = bmt, coefficient = result$out_coefficient, effect.measure1='RR', effect.measure2='RR', outcome.type='SURVIVAL')
  tmp2 <- round(tmp1[1,],digit=5)
  tmp3 <- round(tmp1[2,],digit=5)
  tested <- as.vector(cbind(tmp2,tmp3))
  expected <- c(0.42830, 0.33612, 0.46791, 0.36720)
  expect_equal(expected, tested)
})


test_that("predict.polyreg produced expected of cumulative incidence probabilities of first 2 obs for bmt dataset", {
  data(bmt)
  result <- polyreg(nuisance.model = Event(time, cause)~age+tcell, exposure = 'platelet', cens.model = Event(time,cause==0)~+1, data = bmt, effect.measure1='RR', effect.measure2='RR', time.point=24, outcome.type='COMPETINGRISK')
  tmp1 <- predict.polyreg(formula = Event(time, cause)~age+tcell, exposure = 'platelet', data = bmt, coefficient = result$out_coefficient, effect.measure1='RR', effect.measure2='RR', outcome.type='COMPETINGRISK')
  tmp2 <- round(tmp1[1,],digit=5)
  tmp3 <- round(tmp1[2,],digit=5)
  tested <- as.vector(cbind(tmp2,tmp3))
  expected <- c(0.22138, 0.16741, 0.13667, 0.18875, 0.25956, 0.17265, 0.16025, 0.19467)
  expect_equal(expected, tested)
})
