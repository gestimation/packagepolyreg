calc_ipw <- function(formula, data, specific.time) {
  cl <- match.call()
  mf <- match.call(expand.dots = TRUE)[1:3]
  special <- c("strata", "cluster", "offset")
  Terms <- terms(formula, special, data = data)
  mf$formula <- Terms
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  Y <- model.extract(mf, "response")
  if (!inherits(Y, c("Event", "Surv"))) 
    stop("Expected a 'Surv' or 'Event'-object")
  if (ncol(Y) == 2) {
    t <- Y[, 1]  # time variable
    epsilon <- Y[, 2]  # status variable
    if (any(t<0)) 
      stop("Expected non-negative time variable")
    if (!all(epsilon %in% c(0, 1, 2)))
      stop("Expected only 0, 1 or 2, with 0 representing censoring")
  } else {
    stop("Expected only right censored data")
  }
  if (!nrow(data) == nrow(mf)) 
    stop("Variables contain NA values")

  resC <- phreg(formula, data)
  cens.strata <- resC$strata[order(resC$ord)]
  cens.nstrata <- resC$nstrata
  if (resC$p > 0) kmt <- FALSE
  kmt <- TRUE
  out_predict <- suppressWarnings(predict(resC, newdata = data, type = "survival", times = t, individual.time = TRUE, se = FALSE, km = kmt, tminus = TRUE))
  out_predict_ <- suppressWarnings(predict(resC, newdata = data, type = "survival", times = specific.time, individual.time = TRUE, se = FALSE, km = kmt, tminus = TRUE))
  km <- out_predict$surv
  km_ <- out_predict_$surv[1]
  tmp1 <- ifelse(km > 0, 1 / km, 0)
  tmp2 <- ifelse(km_ > 0, 1 / km_, 0)
  tmp3 <- (t <= specific.time) * (epsilon==0) * tmp1
  tmp4 <- (t > specific.time) * tmp2
  ip_wt <- tmp3 + tmp4
  if (any(is.na(ip_wt)))
    stop("Inverse probability weights contain NA values")
  return(ip_wt)
}

kaplan_meier <- function(t, d){
  n = length(t)
  data = data.frame(t = t, d = d, id = 1:n)
  sorted_data = data[order(data$t), ]
  sorted_t = sorted_data$t
  sorted_d = sorted_data$d
  sorted_id = sorted_data$id
  t_matrix = matrix(rep(sorted_t, each = n), nrow = n, byrow = TRUE)
  atrisk = t(t_matrix) >= t_matrix
  n_atrisk = rowSums(atrisk)
  s = 1 - sorted_d / n_atrisk
  log_s = log(s)
  km = exp(cumsum(log_s))
  data = data.frame(id=sorted_id, km=km)
  sorted_data = data[order(data$id), ]
  km = sorted_data$km
  return(km)
}

d_nelsonaalen <- function(t, d) {
  atrisk <- outer(t, t, ">=")
  n_atrisk <- colSums(atrisk)
  na <- d / n_atrisk
  return(na)
}
