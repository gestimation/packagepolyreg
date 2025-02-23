test_that("polyreg produced expected coefficients for bmt dataset", {
  data(bmt)
  result <- polyreg(nuisance.model = Event(time, cause)~age+tcell, exposure = 'platelet', cens.model = Event(time,cause==0)~+1, data = bmt, effect.measure1='RR', effect.measure2='RR', time.point=24, outcome.type='COMPETINGRISK')
  tested_coefficient <- round(result$out_coefficient,digit=2)
  expected_coefficient <- c(-2.81, 1.01, -1.32, -0.48, -2.65, 0.41, 0.63, 0.12)
  expect_equal(expected_coefficient, tested_coefficient)
})

test_that("polyreg produced expected first row of covariance matrix for bmt dataset", {
  data(bmt)
  result <- polyreg(nuisance.model = Event(time, cause)~age+tcell, exposure = 'platelet', cens.model = Event(time,cause==0)~+1, data = bmt, effect.measure1='RR', effect.measure2='RR', time.point=24, outcome.type='COMPETINGRISK')
  tested_cov <- round(result$out_cov[1,],digit=2)
  expected_cov <- c(0.09, -0.03, 0.00, -0.02, 0.05, -0.02, 0.00, -0.01)
  expect_equal(expected_cov, tested_cov)
})

