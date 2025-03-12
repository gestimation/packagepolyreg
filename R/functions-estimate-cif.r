# TO DO
# consider covariate
# use mets or survival package
# return should be adjusted for ggcuminc

estimateCIF <-  function(formula, data, code.event1, code.event2, code.censoring) {
  # extract variables from model frame
  mf <- model.frame(formula, data)
  response <- model.response(mf)
  t <- response[, 1]
  epsilon <- response[, 2]

  # combine all events
  epsilon_all <- ifelse(epsilon == code.censoring, 0, 1)

  # make df
  data <- data.frame(t = t, epsilon_all = epsilon_all, epsilon = epsilon,
                     id = 1:length(t))

  # sort by t
  sorted_data <- data[order(data$t), ]
  sorted_t <- sorted_data$t
  sorted_epsilon_all <- sorted_data$epsilon_all
  sorted_epsilon <- sorted_data$epsilon
  sorted_id <- sorted_data$id

  # calc risk set for all
  t_matrix <- matrix(rep(sorted_t, each = length(sorted_t)),
                    nrow = length(sorted_t), byrow = TRUE)
  atrisk <- t(t_matrix) >= t_matrix
  n_atrisk <- rowSums(atrisk)
  print(n_atrisk)

  s <- 1 - sorted_epsilon_all / n_atrisk
  log_s <- log(s)
  km <- exp(cumsum(log_s))

  # CIF_j(i) = sum{S(t_{i-1}) * [d_j(i) / Y(i)] }
  epsilon1 <- ifelse(sorted_epsilon == code.event1, 1, 0)
  epsilon2 <- ifelse(sorted_epsilon == code.event2, 1, 0)

  # lag km (km_lag[i] = km[i-1])
  km_lag <- c(1, km[-length(km)])

  # calculate dF1, dF2
  dF1 <- km_lag * (epsilon1 / n_atrisk)
  dF2 <- km_lag * (epsilon2 / n_atrisk)
  print(dF1)

  # cumsum
  CIF1 <- cumsum(dF1)
  CIF2 <- cumsum(dF2)

  result <- data.frame(
    t    = sorted_t,
    KM   = km,
    CIF1 = CIF1,
    CIF2 = CIF2
  )

  return(result)
}
