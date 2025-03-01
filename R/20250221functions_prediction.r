calc_pred <- function(
    alpha_beta_tmp,
    x_a, x_l, offset,
    epsilon,
    estimand,
    optim.method,
    prob.bound,
    initial_pred = NULL
) {
  i_parameter <- rep(NA, 7)
  i_parameter <- calc_i_parameter(i_parameter,x_l,x_a)
  alpha_1 <- alpha_beta_tmp[1:i_parameter[1]]
  beta_1  <- alpha_beta_tmp[i_parameter[2]:i_parameter[3]]
  alpha_2 <- alpha_beta_tmp[i_parameter[4]:i_parameter[5]]
  beta_2  <- alpha_beta_tmp[i_parameter[6]:i_parameter[7]]
  alpha_tmp_1_vals <- x_l %*% as.matrix(alpha_1) + offset
  alpha_tmp_2_vals <- x_l %*% as.matrix(alpha_2) + offset
  beta_tmp_1_vals  <- x_a %*% as.matrix(beta_1)
  beta_tmp_2_vals  <- x_a %*% as.matrix(beta_2)

  pred_list <- vector("list", nrow(x_l))

  # Initialize previous prediction values
  prev_pred <- NULL
  # Loop through each observation in the sorted x_l
  for (i_x in seq_len(nrow(x_l))) {
    # Skip the calculation if the current observation is the same as the previous one
    if (!is.null(prev_pred) && all(x_l[i_x, ] == x_l[i_x-1, ])) {
      pred_list[[i_x]] <- prev_pred
      next
    }

    # Use the previous prediction value of observation i_x if initial_pred is not NULL
    if (!is.null(initial_pred)) {
      log_p0 <- log(initial_pred[i_x, ])
    } else {
      n_epsilon <- length(epsilon)
      p0 <- c(
        sum(epsilon == 1) / n_epsilon + prob.bound,
        sum(epsilon == 2) / n_epsilon + prob.bound,
        sum(epsilon == 1) / n_epsilon + prob.bound,
        sum(epsilon == 2) / n_epsilon + prob.bound
      )
      log_p0 <- log(p0)
    }

    if (optim.method$inner.optim.method == 'optim' | optim.method$inner.optim.method == 'BFGS') {
      eq_fn <- function(lp) {
        estimating_equation_pred(
          log_p          = lp,
          alpha_tmp_1    = alpha_tmp_1_vals[i_x],
          beta_tmp_1     = beta_tmp_1_vals[i_x],
          alpha_tmp_2    = alpha_tmp_2_vals[i_x],
          beta_tmp_2     = beta_tmp_2_vals[i_x],
          estimand       = estimand,
          optim.method   = optim.method,
          prob.bound     = prob.bound
        )
      }
      sol <- optim(par = log_p0, fn = eq_fn, method = "BFGS", control = list(maxit = optim.method$optim.parameter8, reltol = optim.method$optim.parameter7))
      pred_list[[i_x]] <- exp(sol$par)
      prev_pred <- pred_list[[i_x]]
    } else if (optim.method$inner.optim.method == 'SANN') {
      eq_fn <- function(lp) {
        estimating_equation_pred(
          log_p          = lp,
          alpha_tmp_1    = alpha_tmp_1_vals[i_x],
          beta_tmp_1     = beta_tmp_1_vals[i_x],
          alpha_tmp_2    = alpha_tmp_2_vals[i_x],
          beta_tmp_2     = beta_tmp_2_vals[i_x],
          estimand       = estimand,
          optim.method   = optim.method,
          prob.bound     = prob.bound
        )
      }
      sol <- optim(par = log_p0, fn = eq_fn, method = "SANN", control = list(maxit = optim.method$optim.parameter8, reltol = optim.method$optim.parameter7))
      pred_list[[i_x]] <- exp(sol$par)
      prev_pred <- pred_list[[i_x]]
    } else if (optim.method$inner.optim.method == 'multiroot'){
      eq_fn <- function(lp) {
        estimating_equation_pred(
          log_p          = lp,
          alpha_tmp_1    = alpha_tmp_1_vals[i_x],
          beta_tmp_1     = beta_tmp_1_vals[i_x],
          alpha_tmp_2    = alpha_tmp_2_vals[i_x],
          beta_tmp_2     = beta_tmp_2_vals[i_x],
          estimand       = estimand,
          optim.method   = optim.method,
          prob.bound     = prob.bound
        )
      }
      sol <- multiroot(eq_fn, start = log_p0)
      pred_list[[i_x]] <- exp(sol$root)
      prev_pred <- pred_list[[i_x]]
    }
  }
#  pred_observed <- matrix(NA,nrow(x_l),2)
#  for (i_x in seq_len(nrow(x_l))) {
#    tmp1 <- pred_list[[i_x]]
#    tmp2 <- cbind(1-sum(x_a[i_x,]),x_a[i_x,])
#    pred_observed[i_x,1] <- tmp2 %*% tmp1[1:(length(tmp1)/2)]
#    pred_observed[i_x,2] <- tmp2 %*% tmp1[(length(tmp1)/2+1):length(tmp1)]
#  }
  out <- do.call(rbind, pred_list)
  return(out)
  #  return(pred_observed)
}

estimating_equation_pred <- function(
    log_p,      # c(log_p10, log_p20, log_p11, log_p21) > c(log_p10, log_p11, log_p20, log_p21)
    alpha_tmp_1, beta_tmp_1,
    alpha_tmp_2, beta_tmp_2,
    estimand, optim.method, prob.bound
) {
  objective_function <- function(log_p) {
    ret <- numeric(length(log_p))
    clog_p <- numeric(length(log_p))

    clamp_log_p <- function(lp) {
      val <- exp(lp)
      if (val < prob.bound) {
        return(log(prob.bound))
      } else if ((1 - val) < prob.bound) {
        return(log(1 - prob.bound))
      } else {
        return(lp)
      }
    }
    for (i in seq_len(length(log_p))) {
      clog_p[i] <- clamp_log_p(log_p[i])
    }

    log_rr1 <- log_p[2:(length(clog_p)/2)]/clog_p[1]
    log_rr2 <- log_p[(length(clog_p)/2+2):length(clog_p)]/clog_p[length(clog_p)/2+1]
    plor_1n <- sum(clog_p[1:(length(clog_p)/2)])
    plor_2n <- sum(clog_p[(length(clog_p)/2+1):length(clog_p)])
    p <- exp(clog_p)
    p01 <- 1-sum(p[1:(length(clog_p)/2)])
    p02 <- 1-sum(p[(length(clog_p)/2+1):length(clog_p)])
    log_p01 <- log(p01)
    log_p02 <- log(p02)
    plor_d <- log_p01 + log_p02

    if ((1 - sum(p[1:(length(clog_p)/2)]) < prob.bound) |
      (1 - sum(p[(length(clog_p)/2+1):length(clog_p)]) < prob.bound)) {
        plor_d <- log(prob.bound)
        #      lp0102 <- log(prob.bound)
      } else {
        p01 <- 1-sum(p[1:(length(clog_p)/2)])
        p02 <- 1-sum(p[(length(clog_p)/2+1):length(clog_p)])
        log_p01 <- log(p01)
        log_p02 <- log(p02)
        plor_d <- log_p01 + log_p02
      }
    plor_1 <- plor_1n - plor_d
    plor_2 <- plor_2n - plor_d

    log_p1 <- (log_p[1])
    log_p2 <- (log_p[2])
    log_p3 <- (log_p[3])
    log_p4 <- (log_p[4])

if (estimand$effect.measure1 == 'RR') {
  #ret[1] <- alpha_tmp_1 - plor_1
  #ret[2] <- beta_tmp_1  - log_rr1
        ret[1] <- alpha_tmp_1 - log_p1 - log_p2 + plor_d
        ret[2] <- beta_tmp_1  - log_p2 + log_p1
  #      ret[1] <- alpha_tmp_1 - log_p1 - log_p3 + lp0102
  #      ret[2] <- beta_tmp_1  - log_p3 + log_p1
} else if (estimand$effect.measure1 == 'OR') {
  ret[1] <- alpha_tmp_1 - plor_1
  ret[2] <- beta_tmp_1  - log_p3 + log_p1 + log(1 - exp_lp3) - log(1 - exp_lp1)
} else if (estimand$effect.measure1 == 'SHR') {
  ret[1] <- alpha_tmp_1 - plor_1
  ret[2] <- exp(beta_tmp_1) - ( log(1 - exp_lp3) / log(1 - exp_lp1) )
} else {
  stop("Invalid effect_measure. Must be RR, OR or SHR.")
}

if (estimand$effect.measure2 == 'RR') {
  #ret[3] <- alpha_tmp_2 - plor_2
  #ret[4] <- beta_tmp_2  - log_rr2
  ret[3] <- alpha_tmp_1 - log_p3 - log_p4 + plor_d
  ret[4] <- beta_tmp_1  - log_p4 + log_p3
} else if (estimand$effect.measure2 == 'OR') {
  ret[3] <- alpha_tmp_2 - plor_2
  ret[4] <- beta_tmp_2  - log_p4 + log_p2 + log(1 - exp_lp4) - log(1 - exp_lp2)
} else if (estimand$effect.measure2 == 'SHR') {
  ret[3] <- alpha_tmp_2 - plor_2
  ret[4] <- exp(beta_tmp_2) - ( log(1 - exp_lp4) / log(1 - exp_lp2) )
} else {
  stop("Invalid effect_measure. Must be RR, OR or SHR.")
}
if (optim.method$inner.optim.method == 'multiroot'){
  return(ret)
} else {
  return(sum(ret^2))
}
  }
  return(objective_function(log_p))
}

calc_pred_survival <- function(alpha_beta, x_a, x_l, offset, estimand) {
  if (estimand$effect.measure1 == 'RR') {
    one <- rep(1, nrow(x_l))
    n_para_1 <- ncol(x_l)
    n_para_2 <- ncol(x_l) + 1
    alpha_beta_ <- as.matrix(as.vector(alpha_beta))
    if (ncol(alpha_beta_ == 1)) {
      alpha_beta_ <- t(alpha_beta_)
    }
    phi <- x_l %*% alpha_beta_[, 1:n_para_1] + offset
    theta <- one * alpha_beta_[, n_para_2]
    expphi <- exp(phi)
    exptheta <- exp(theta)
    if (all(phi == 0)) {
      p_10 <- one / (one + exptheta)
      p_11 <- exptheta * p_10
    } else {
      denomi_1 <- -(exptheta + one) * expphi
      denomi_2 <- sqrt(exp(2 * phi) * (exptheta + one) * (exptheta + one) + 4 * exp(theta + phi) * (one - expphi))
      denomi <- denomi_1 + denomi_2
      numera <- 2 * exptheta * (one - expphi)
      p_10 <- denomi / numera
      p_11 <- exptheta * p_10
    }
  } else if (estimand$effect.measure1 == 'OR') {
    one <- rep(1, nrow(x_l))
    n_para_1 <- ncol(x_l)
    n_para_2 <- ncol(x_l) + 1
    phi <- x_l %*% as.matrix(alpha_beta)[, 1:n_para_1] + offset
    theta <- one * as.matrix(alpha_beta)[, n_para_2]
    sqrt1 <- sqrt(exp(-theta - phi))
    sqrt2 <- sqrt(exp(theta - phi))
    if (all(phi == theta)) {
      p_10 <- 0.5 * one
      p_11 <- one/(one + sqrt1)
    } else {
      p_10 <- one/(one + sqrt2)
      p_11 <- one/(one + sqrt1)
    }
  } else {
    stop("Invalid effect_measure. Must be RR or OR.")
  }
  return(cbind(p_10, p_11))
}

