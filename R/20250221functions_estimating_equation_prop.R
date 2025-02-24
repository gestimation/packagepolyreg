estimating_equation_ipcw <- function(
    formula,
    data,
    exposure,
    ip_weight,
    alpha_beta,
    estimand,
    optim.method,
    prob.bound,
    initial_pred = NULL
) {
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
    if (any(t<0))
      stop("Expected non-negative time variable")
    epsilon <- Y[, 2]  # status variable
    if (!all(epsilon %in% c(0, 1, 2)))
      stop("Expected only 0, 1 or 2, with 0 representing censoring")
    time.name <- all.vars(formula)[1]
    status.name <- all.vars(formula)[2]
  } else {
    stop("Expected only right censored data")
  }
  if (!is.null(offsetpos <- attributes(Terms)$specials$offset)) {
    ts <- survival::untangle.specials(Terms, "offset")
    if (length(ts$vars) > 0) {
      Terms <- Terms[-ts$terms]
      offset <- mf[[ts$vars]]
    } else {
      offset <- rep(0, nrow(mf))
    }
  } else {
    offset <- rep(0, nrow(mf))
  }

  y_0 <- ifelse(epsilon == 0 | t > estimand$time.point, 1, 0)
  y_1 <- ifelse(epsilon == 1 & t <= estimand$time.point, 1, 0)
  y_2 <- ifelse(epsilon == 2 & t <= estimand$time.point, 1, 0)

  y_0_ <- ifelse(epsilon == 0, 1, 0)
  y_1_ <- ifelse(epsilon == 1, 1, 0)
  y_2_ <- ifelse(epsilon == 2, 1, 0)

  a_ <- as.factor(data[[exposure]])
  a <- model.matrix(~ a_)[, 2]
  epsilon0 <- epsilon[a == 0]
  epsilon1 <- epsilon[a == 1]
  x_l <- model.matrix(Terms, mf)

  n_para_1 <- ncol(x_l)
  n_para_2 <- n_para_1 + 1
  n_para_3 <- n_para_1 + 2
  n_para_4 <- 2*n_para_1 + 1
  n_para_5 <- 2*n_para_1 + 2
  one <- rep(1, nrow(x_l))

  pred <- calc_pred(alpha_beta,x_l,offset,epsilon0,epsilon1,one,n_para_1,estimand,optim.method,prob.bound,initial_pred)

  ey_1 <- pred[,3]*a + pred[,1]*(one - a)
  ey_2 <- pred[,4]*a + pred[,2]*(one - a)

  v11 <- ey_1 * (1 - ey_1)
  v12 <- -ey_1 * ey_2
  v22 <- ey_2 * (1 - ey_2)
  denom <- v11*v22 - v12*v12

  w11 <- v22 / denom
  w12 <- -v12 / denom
  w22 <- v11 / denom

  wy_1 <- ip_weight * y_1
  wy_2 <- ip_weight * y_2

  wy_1ey_1 <- w11*(wy_1 - ey_1) + w12*(wy_2 - ey_2)
  wy_2ey_2 <- w12*(wy_1 - ey_1) + w22*(wy_2 - ey_2)

  # d = ( [x_l||a, 0] ; [0, x_l||a] )
  x_la <- cbind(x_l, a)
  zero <- matrix(0, nrow=nrow(x_la), ncol=ncol(x_la))

  tmp1 <- cbind(x_la, zero)
  tmp2 <- cbind(zero, x_la)
  d    <- rbind(tmp1, tmp2)

  residual <- c(wy_1ey_1, wy_2ey_2)
  ret <- as.vector(t(d) %*% residual / nrow(x_l))

  n_col_d <- ncol(d)
  score_mat <- matrix(NA, nrow=nrow(d), ncol=n_col_d)
  for (j in seq_len(n_col_d)) {
    score_mat[,j] <- d[,j]*residual
  }

  out <- list(
    ret   = ret,
    score = score_mat,
    ey_1 = ey_1,
    ey_2 = ey_2,
    w11 = w11,
    w12 = w12,
    w22 = w22,
    pred = pred,
    t = t,
    y_0 = y_0,
    y_1 = y_1,
    y_2 = y_2,
    y_0_ = y_0_,
    y_1_ = y_1_,
    y_2_ = y_2_,
    a = a,
    x_l = x_l,
    ip_weight = ip_weight,
    n_para_1 = n_para_1,
    n_para_2 = n_para_2,
    n_para_3 = n_para_3,
    n_para_4 = n_para_4,
    n_para_5 = n_para_5
  )
  return(out)
}

estimating_equation_partial <- function(
    formula,
    data,
    exposure,
    ip_weight,
    alpha_beta_partial,
    alpha_beta1_current,
    alpha_beta2_current,
    optimizing_event_1,
    estimand,
    optim.method,
    prob.bound,
    initial_pred = NULL
) {
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
    if (any(t<0))
      stop("Expected non-negative time variable")
    epsilon <- Y[, 2]  # status variable
    if (!all(epsilon %in% c(0, 1, 2)))
      stop("Expected only 0, 1 or 2, with 0 representing censoring")
    time.name <- all.vars(formula)[1]
    status.name <- all.vars(formula)[2]
  } else {
    stop("Expected only right censored data")
  }
  if (!is.null(offsetpos <- attributes(Terms)$specials$offset)) {
    ts <- survival::untangle.specials(Terms, "offset")
    Terms <- Terms[-ts$terms]
    offset <- m[[ts$vars]]
  } else {
    offset <- rep(0, length(t))
  }

  y_0 <- ifelse(epsilon == 0 | t > estimand$time.point, 1, 0)
  y_1 <- ifelse(epsilon == 1 & t <= estimand$time.point, 1, 0)
  y_2 <- ifelse(epsilon == 2 & t <= estimand$time.point, 1, 0)

  y_0_ <- ifelse(epsilon == 0, 1, 0)
  y_1_ <- ifelse(epsilon == 1, 1, 0)
  y_2_ <- ifelse(epsilon == 2, 1, 0)

  a_ <- as.factor(data[[exposure]])
  a <- model.matrix(~ a_)[, 2]
  epsilon0 <- epsilon[a == 0]
  epsilon1 <- epsilon[a == 1]
  x_l <- model.matrix(Terms, mf)

  n_para_1 <- ncol(x_l)
  n_para_2 <- n_para_1 + 1
  n_para_3 <- n_para_1 + 2
  n_para_4 <- 2*n_para_1 + 1
  n_para_5 <- 2*n_para_1 + 2
  one <- rep(1, nrow(x_l))
  alpha_beta <- NULL

  if (optimizing_event_1 == TRUE) {
    alpha_beta[1:n_para_2] <- alpha_beta_partial
    alpha_beta[n_para_3:n_para_5] <- alpha_beta2_current
  } else {
    alpha_beta[1:n_para_2] <- alpha_beta1_current
    alpha_beta[n_para_3:n_para_5] <- alpha_beta_partial
  }

  pred <- calc_pred(alpha_beta,x_l,offset,epsilon0,epsilon1,one,n_para_1,estimand,optim.method,prob.bound,initial_pred)
  # pred: n?~4 => columns: p_10, p_20, p_11, p_21

  ey_1 <- pred[,3]*a + pred[,1]*(one - a)
  ey_2 <- pred[,4]*a + pred[,2]*(one - a)

  v11 <- ey_1 * (1 - ey_1)
  v12 <- -ey_1 * ey_2
  v22 <- ey_2 * (1 - ey_2)
  denom <- v11*v22 - v12*v12

  w11 <- v22 / denom
  w12 <- -v12 / denom
  w22 <- v11 / denom

  wy_1 <- ip_weight * y_1
  wy_2 <- ip_weight * y_2

  wy_1ey_1 <- w11*(wy_1 - ey_1) + w12*(wy_2 - ey_2)
  wy_2ey_2 <- w12*(wy_1 - ey_1) + w22*(wy_2 - ey_2)

  # d = ( [x_l||a, 0] ; [0, x_l||a] )
  x_la <- cbind(x_l, a)
  zero <- matrix(0, nrow=nrow(x_la), ncol=ncol(x_la))

  tmp1 <- cbind(x_la, zero)
  tmp2 <- cbind(zero, x_la)
  d    <- rbind(tmp1, tmp2)

  # residual = (wy_1ey_1, wy_2ey_2)
  residual <- c(wy_1ey_1, wy_2ey_2)
  ret <- as.vector(t(d) %*% residual / nrow(x_l))

  if (optimizing_event_1 == TRUE) {
    ret <- ret[1:n_para_2]
  } else {
    ret <- ret[n_para_3:n_para_5]
  }

  n_col_d <- ncol(d)
  score_mat <- matrix(NA, nrow=nrow(d), ncol=n_col_d)  # Score used for sandwich variance
  for (j in seq_len(n_col_d)) {
    score_mat[,j] <- d[,j]*residual
  }

  out <- list(
    ret   = ret,
    score = score_mat,
    ey_1 = ey_1,
    ey_2 = ey_2,
    w11 = w11,
    w12 = w12,
    w22 = w22,
    pred = pred,
    t = t,
    y_0 = y_0,
    y_1 = y_1,
    y_2 = y_2,
    y_0_ = y_0_,
    y_1_ = y_1_,
    y_2_ = y_2_,
    a = a,
    x_l = x_l,
    ip_weight = ip_weight,
    n_para_1 = n_para_1,
    n_para_2 = n_para_2,
    n_para_3 = n_para_3,
    n_para_4 = n_para_4,
    n_para_5 = n_para_5
  )
  return(out)
}

calc_cov <- function(objget_results, estimand, prob.bound)
  {
  w11 <- objget_results$w11
  w12 <- objget_results$w12
  w22 <- objget_results$w22
  score <- objget_results$score
  pred <- objget_results$pred
  t <- objget_results$t
  y_0_ <- objget_results$y_0_
  y_1 <- objget_results$y_1
  y_2 <- objget_results$y_2
  a <- objget_results$a
  x_l <- objget_results$x_l
  n_para_2 <- objget_results$n_para_2
  n_para_3 <- objget_results$n_para_3
  n_para_5 <- objget_results$n_para_5
  n <- length(t)
  one <- rep(1, n)

  ey_1		= pred[,3]*a + pred[,1]*(one-a);
  ey_2		= pred[,4]*a + pred[,2]*(one-a);

  censoring_dna	<- d_nelsonaalen(t,y_0_)
  censoring_martingale <- diag(y_0_) - (outer(t, t, ">") * censoring_dna)
  censoring_km <- kaplan_meier(t, y_0_)
  censoring_km[censoring_km == 0] <- 1e-5
  censoring_mkm <- censoring_martingale / censoring_km
  y_12 <- (y_1 + y_2 > 0)
  survival_km <- kaplan_meier(t, y_12)
  wy_1 <- w11 * (y_1 - ey_1) + w12 * (y_2 - ey_2)
  wy_2 <- w12 * (y_1 - ey_1) + w22 * (y_2 - ey_2)
  x_la <- cbind(x_l, a)
  AB1 <- score[1:n, 1:n_para_2]
  AB2 <- score[(n + 1):(2 * n), n_para_3:n_para_5]
  for (i_para in 1:n_para_2) {
    tmp0 <- x_la[, i_para]
    use <- (t <= estimand$time.point)
    tmp1 <- colSums((use * tmp0) * (outer(t, t, ">=") * wy_1))
    tmp1 <- tmp1 / survival_km / n
    integrand1 <- tmp1 * censoring_mkm
    tmp2 <- colSums((use * tmp0) * (outer(t, t, ">=") * wy_2))
    tmp2 <- tmp2 / survival_km / n
    integrand2 <- tmp2 * censoring_mkm
    for (i_score in 1:n) {
      integral1 <- cumsum(integrand1[i_score, ])
      AB1[i_score, i_para] <- AB1[i_score, i_para] + integral1[n]
      integral2 <- cumsum(integrand2[i_score, ])
      AB2[i_score, i_para] <- AB2[i_score, i_para] + integral2[n]
    }
  }

  total_score <- cbind(AB1, AB2)
  cov_AB <- crossprod(total_score,total_score) / n

  calc_d_out <- calc_d(pred, a, x_l, estimand, prob.bound)
  hesse_d11 <- crossprod(x_la, w11 * calc_d_out$d_11) / n
  hesse_d12 <- crossprod(x_la, w12 * calc_d_out$d_12) / n
  hesse_d22 <- crossprod(x_la, w22 * calc_d_out$d_22) / n

  hesse_d1 <- cbind(hesse_d11, hesse_d12)
  hesse_d2 <- cbind(hesse_d12, hesse_d22)
  hesse <- rbind(hesse_d1, hesse_d2)

#  cov_estimated_ <- solve(hesse) %*% cov_AB %*% t(solve(hesse)) / n
  influence_function <- total_score %*% t(solve(hesse))
  cov_estimated <- crossprod(influence_function,influence_function) / n / n

  return(list(cov_estimated = cov_estimated, score_function = total_score, influence_function = influence_function))
}

calc_d <- function(pred, a, x_l, estimand, prob.bound) {
  pred_1 <- ifelse(pred[, 1] == 0, prob.bound, ifelse(pred[, 1] == 1, 1 - prob.bound, pred[, 1]))
  pred_2 <- ifelse(pred[, 2] == 0, prob.bound, ifelse(pred[, 2] == 1, 1 - prob.bound, pred[, 2]))
  pred_3 <- ifelse(pred[, 3] == 0, prob.bound, ifelse(pred[, 3] == 1, 1 - prob.bound, pred[, 3]))
  pred_4 <- ifelse(pred[, 4] == 0, prob.bound, ifelse(pred[, 4] == 1, 1 - prob.bound, pred[, 4]))
  a11 <- a12 <- a22 <- NULL
  n <- length(a)

  if (estimand$effect.measure1 == 'RR') {
    a11 <- a * (1 / pred_3) + (1 - a) * (1 / pred_1)
  } else if (estimand$effect.measure1 == 'OR') {
    a11 <- a * (1 / pred_3 + 1 / (1 - pred_3)) + (1 - a) * (1 / pred_1 + 1 / (1 - pred_1))
  } else if (estimand$effect.measure1 == 'SHR') {
    tmp31 <- -1 / (1 - pred_3)
    tmp32 <- log(1 - pred_3)
    tmp11 <- -1 / (1 - pred_1)
    tmp12 <- log(1 - pred_1)
    a11 <- a * (tmp31 / tmp32) + (1 - a) * (tmp11 / tmp12)
  } else {
    stop("Invalid effect_measure. Must be RR, OR or SHR.")
  }

  if (estimand$effect.measure2 == 'RR') {
    a22 <- a * (1 / pred_4) + (1 - a) * (1 / pred_2)
  } else if (estimand$effect.measure2 == 'OR') {
    a22 <- a * (1 / pred_4 + 1 / (1 - pred_4)) + (1 - a) * (1 / pred_2 + 1 / (1 - pred_2))
  } else if (estimand$effect.measure2 == 'SHR') {
    tmp41 <- -1 / (1 - pred_4)
    tmp42 <- log(1 - pred_4)
    tmp21 <- -1 / (1 - pred_2)
    tmp22 <- log(1 - pred_2)
    a22 <- a * (tmp41 / tmp42) + (1 - a) * (tmp21 / tmp22)
  } else {
    stop("Invalid effect_measure. Must be RR, OR or SHR.")
  }
  a12 <- matrix(0, nrow = n, ncol = 1)

  d_ey_d_beta_11 <- a22 / (a11 * a22 - a12 * a12)
  d_ey_d_beta_12 <- -a12 / (a11 * a22 - a12 * a12)
  d_ey_d_beta_22 <- a11 / (a11 * a22 - a12 * a12)
  c12 <- a * (1 / (1 - pred_3 - pred_4)) + (1 - a) * (1 / (1 - pred_1 - pred_2))
  c11 <- a11 + c12
  c22 <- a22 + c12
  d_ey_d_alpha_11 <- c22 / (c11 * c22 - c12 * c12)
  d_ey_d_alpha_12 <- -c12 / (c11 * c22 - c12 * c12)
  d_ey_d_alpha_22 <- c11 / (c11 * c22 - c12 * c12)
  d_11 <- cbind((d_ey_d_alpha_11 * x_l), (d_ey_d_beta_11 * a))
  d_12 <- cbind((d_ey_d_alpha_12 * x_l), (d_ey_d_beta_12 * a))
  d_22 <- cbind((d_ey_d_alpha_22 * x_l), (d_ey_d_beta_22 * a))
  return(list(d_11 = d_11, d_12 = d_12, d_22 = d_22))
}

estimating_equation_survival <- function(
    formula,
    data,
    exposure,
    ip_weight,
    alpha_beta,
    estimand, prob.bound, initial_pred = NULL
) {
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
    if (any(t<0))
      stop("Expected non-negative time variable")
    epsilon <- Y[, 2]  # status variable
    #    if (!all(epsilon %in% c(0, 1, 2)))
    #      stop("Expected only 0, 1 or 2, with 0 representing censoring")
  } else {
    stop("Expected only right censored data")
  }
  if (!is.null(offsetpos <- attributes(Terms)$specials$offset)) {
    ts <- survival::untangle.specials(Terms, "offset")
    Terms <- Terms[-ts$terms]
    offset <- m[[ts$vars]]
  } else {
    offset <- rep(0, length(t))
  }

  y_0 <- ifelse(epsilon == 0 | t > estimand$time.point, 1, 0)
  y_1 <- ifelse(epsilon == 1 & t <= estimand$time.point, 1, 0)
  y_0_ <- ifelse(epsilon == 0, 1, 0)
  y_1_ <- ifelse(epsilon == 1, 1, 0)

  a_ <- as.factor(data[[exposure]])
  a <- model.matrix(~ a_)[, 2]
  epsilon0 <- epsilon[a == 0]
  epsilon1 <- epsilon[a == 1]
  x_l <- model.matrix(Terms, mf)

  n_para_1 <- ncol(x_l)
  n_para_2 <- n_para_1 + 1
  n_para_3 <- n_para_1 + 2
  one <- rep(1, nrow(x_l))

  pred <- calc_pred_survival(alpha_beta, x_l, offset, estimand)
  ey_1 <- pred[,2]*a + pred[,1]*(one - a)
  w11 <- 1 / (ey_1 * (1 - ey_1))
  wy_1 <- ip_weight * y_1
  wy_1ey_1 <- w11*(wy_1 - ey_1)
  d <- cbind(x_l, a)
  residual <- wy_1ey_1
  ret <- as.vector(t(d) %*% residual / nrow(x_l))

  n_col_d <- ncol(d)
  score_mat <- matrix(NA, nrow=nrow(d), ncol=n_col_d)
  for (j in seq_len(n_col_d)) {
    score_mat[,j] <- d[,j]*residual
  }

  out <- list(
    ret   = ret,
    score = score_mat,
    ey_1 = ey_1,
    w11 = w11,
    pred = pred,
    t = t,
    y_0 = y_0,
    y_1 = y_1,
    y_0_ = y_0_,
    y_1_ = y_1_,
    a = a,
    x_l = x_l,
    ip_weight = ip_weight,
    n_para_1 = n_para_1,
    n_para_2 = n_para_2,
    n_para_3 = n_para_3
  )
  return(out)
}

calc_cov_survival <- function(objget_results, estimand, prob.bound)
{
  w11 <- objget_results$w11
  score <- objget_results$score
  pred <- objget_results$pred
  t <- objget_results$t
  y_0_ <- objget_results$y_0_
  y_1 <- objget_results$y_1
  a <- objget_results$a
  x_l <- objget_results$x_l
  n_para_2 <- objget_results$n_para_2
  n <- length(t)
  one <- rep(1, n)

  ey_1		= pred[,2]*a + pred[,1]*(one-a);
  censoring_dna	<- d_nelsonaalen(t,y_0_)
  censoring_martingale <- diag(y_0_) - (outer(t, t, ">") * censoring_dna)
  censoring_km <- kaplan_meier(t, y_0_)
  censoring_km[censoring_km == 0] <- 1e-5
  censoring_mkm <- censoring_martingale / censoring_km
  y_12 <- (y_1 > 0)
  survival_km <- kaplan_meier(t, y_12)
  wy_1 <- w11 * (y_1 - ey_1)
  x_la <- cbind(x_l, a)
  AB1 <- score[1:n, 1:n_para_2]
  for (i_para in 1:n_para_2) {
    tmp0 <- x_la[, i_para]
    use <- (t <= estimand$time.point)
    tmp1 <- colSums((use * tmp0) * (outer(t, t, ">=") * wy_1))
    tmp1 <- tmp1 / survival_km / n
    integrand1 <- tmp1 * censoring_mkm
    for (i_score in 1:n) {
      integral1 <- cumsum(integrand1[i_score, ])
      AB1[i_score, i_para] <- AB1[i_score, i_para] + integral1[n]
    }
  }

  total_score <- AB1
  cov_AB <- crossprod(total_score,total_score) / n

  calc_d_out <- calc_d_survival(pred=pred, a=a, x_l=x_l, estimand=estimand, prob.bound=prob.bound)
  hesse <- crossprod(x_la, w11 * calc_d_out) / n
  #  cov_estimated_ <- solve(hesse) %*% cov_AB %*% t(solve(hesse)) / n
  influence_function <- total_score %*% t(solve(hesse))
  cov_estimated <- crossprod(influence_function,influence_function) / n / n

  return(list(cov_estimated = cov_estimated, score_function = total_score, influence_function = influence_function))
}

calc_d_survival <- function(pred, a, x_l, estimand, prob.bound) {
  pred_1 <- ifelse(pred[, 1] == 0, prob.bound, ifelse(pred[, 1] == 1, 1 - prob.bound, pred[, 1]))
  pred_2 <- ifelse(pred[, 2] == 0, prob.bound, ifelse(pred[, 2] == 1, 1 - prob.bound, pred[, 2]))
  a11 <- a12 <- a22 <- NULL
  n <- length(a)

  if (estimand$effect.measure1 == 'RR') {
    a11 <- a * (1 / pred_2) + (1 - a) * (1 / pred_1)
  } else if (estimand$effect.measure1 == 'OR') {
    a11 <- a * (1 / pred_2 + 1 / (1 - pred_2)) + (1 - a) * (1 / pred_1 + 1 / (1 - pred_1))
  } else if (estimand$effect.measure1 == 'SHR') {
    tmp31 <- -1 / (1 - pred_2)
    tmp32 <- log(1 - pred_2)
    tmp11 <- -1 / (1 - pred_1)
    tmp12 <- log(1 - pred_1)
    a11 <- a * (tmp31 / tmp32) + (1 - a) * (tmp11 / tmp12)
  } else {
    stop("Invalid effect_measure. Must be RR, OR or SHR.")
  }
  d_ey_d_beta_11 <- 1 / a11
  d_ey_d_alpha_11 <- 1 / a11
  d_11 <- cbind((d_ey_d_alpha_11 * x_l), (d_ey_d_beta_11 * a))
  return(d_11)
}

estimating_equation_proportional <- function(
    formula,
    data,
    exposure,
    ip_wt_matrix,
    alpha_beta,
    estimand,
    time.point,
    optim.method,
    prob.bound,
    initial_pred = NULL
) {
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
    if (any(t<0))
      stop("Expected non-negative time variable")
    epsilon <- Y[, 2]  # status variable
    if (!all(epsilon %in% c(0, 1, 2)))
      stop("Expected only 0, 1 or 2, with 0 representing censoring")
    time.name <- all.vars(formula)[1]
    status.name <- all.vars(formula)[2]
  } else {
    stop("Expected only right censored data")
  }
  if (!is.null(offsetpos <- attributes(Terms)$specials$offset)) {
    ts <- survival::untangle.specials(Terms, "offset")
    Terms <- Terms[-ts$terms]
    offset <- m[[ts$vars]]
  } else {
    offset <- rep(0, length(t))
  }

  y_0_ <- ifelse(epsilon == 0, 1, 0) # ?ł??؂??ϐ?
  y_1_ <- ifelse(epsilon == 1, 1, 0)
  y_2_ <- ifelse(epsilon == 2, 1, 0)

  a_ <- as.factor(data[[exposure]])
  a <- model.matrix(~ a_)[, 2]
  epsilon0 <- epsilon[a == 0]
  epsilon1 <- epsilon[a == 1]
  x_l <- model.matrix(Terms, mf)

  n_para_1 <- ncol(x_l)
  n_para_2 <- n_para_1 + 1
  n_para_3 <- n_para_1 + 2
  n_para_4 <- 2*n_para_1 + 1
  n_para_5 <- 2*n_para_1 + 2
  n_para_6 <- length(time.point)*(n_para_5-2) + 2
  one <- rep(1, nrow(x_l))

  score_beta <- NULL
  score_alpha1 <- 0
  score_alpha2 <- 0
  i_time <- 0
  alpha_beta_i <- rep(NA, n_para_5)
  for (specific.time in time.point) {
    i_time <- i_time + 1
    i_para <- n_para_1*(i_time-1)+1
    alpha_beta_i[1:n_para_1]        <- alpha_beta[i_para:(i_para+n_para_1-1)]
    alpha_beta_i[n_para_2]          <- alpha_beta[n_para_6/2]
    alpha_beta_i[n_para_3:n_para_4] <- alpha_beta[(n_para_6/2+i_para):(n_para_6/2+i_para+n_para_1-1)]
    alpha_beta_i[n_para_5]          <- alpha_beta[n_para_6]

    y_0 <- ifelse(epsilon == 0 | t > specific.time, 1, 0)
    y_1 <- ifelse(epsilon == 1 & t <= specific.time, 1, 0)
    y_2 <- ifelse(epsilon == 2 & t <= specific.time, 1, 0)

    pred <- calc_pred(alpha_beta_i,x_l,offset,epsilon0,epsilon1,one,n_para_1,estimand,optim.method,prob.bound,initial_pred)
    # pred: n?~4 => columns: p_10, p_20, p_11, p_21

    ey_1 <- pred[,3]*a + pred[,1]*(one - a)
    ey_2 <- pred[,4]*a + pred[,2]*(one - a)

    v11 <- ey_1 * (1 - ey_1)
    v12 <- -ey_1 * ey_2
    v22 <- ey_2 * (1 - ey_2)
    denom <- v11*v22 - v12*v12

    w11 <- v22 / denom
    w12 <- -v12 / denom
    w22 <- v11 / denom

    wy_1 <- ip_wt_matrix[,i_time] * y_1
    wy_2 <- ip_wt_matrix[,i_time] * y_2

    wy_1ey_1 <- w11*(wy_1 - ey_1) + w12*(wy_2 - ey_2)
    wy_2ey_2 <- w12*(wy_1 - ey_1) + w22*(wy_2 - ey_2)

    # d = ( [x_l||a, 0] ; [0, x_l||a] )
    x_la <- cbind(x_l, a)
    zero <- matrix(0, nrow=nrow(x_la), ncol=ncol(x_la))
    tmp1 <- cbind(x_la, zero)
    tmp2 <- cbind(zero, x_la)
    d    <- rbind(tmp1, tmp2)

    residual <- c(wy_1ey_1, wy_2ey_2)
    subscore <- as.matrix(t(d) %*% residual / nrow(x_l))
    tmp1 <- t(subscore[1:n_para_1,])
    tmp2 <- t(subscore[n_para_3:n_para_4,])
    score_beta <- cbind(score_beta, tmp1, tmp2)
    score_alpha1 <- score_alpha1 + subscore[n_para_2,]
    score_alpha2 <- score_alpha2 + subscore[n_para_5,]
  }
  score <- cbind(score_beta, score_alpha1, score_alpha2)
  out <- list(ret = score)
  return(out)
}

