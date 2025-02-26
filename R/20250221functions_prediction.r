calc_pred <- function(
    alpha_beta_tmp,
    x_em, x_l, offset,
    epsilon0, epsilon1,
    one,
    n_para_1,
    estimand,
    optim.method,
    prob.bound,
    initial_pred = NULL
) {
  i_parameter <- rep(NA, 7)
  i_parameter[1] <- ncol(x_l)
  i_parameter[2] <- i_parameter[1] + 1
  i_parameter[3] <- i_parameter[1] + ncol(x_em)
  i_parameter[4] <- i_parameter[1] + ncol(x_em) + 1
  i_parameter[5] <- 2 * i_parameter[1] + ncol(x_em)
  i_parameter[6] <- 2 * i_parameter[1] + ncol(x_em) + 1
  i_parameter[7] <- 2 * i_parameter[1] + 2 * ncol(x_em)
  alpha_1 <- alpha_beta_tmp[1:i_parameter[1]]
  beta_1  <- alpha_beta_tmp[i_parameter[2]:i_parameter[3]]
  alpha_2 <- alpha_beta_tmp[i_parameter[4]:i_parameter[5]]
  beta_2  <- alpha_beta_tmp[i_parameter[6]:i_parameter[7]]
  alpha_tmp_1_vals <- x_l %*% as.matrix(alpha_1) + offset
  alpha_tmp_2_vals <- x_l %*% as.matrix(alpha_2) + offset
  beta_tmp_1_vals  <- x_em %*% as.matrix(beta_1)
  beta_tmp_2_vals  <- x_em %*% as.matrix(beta_2)

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
      n_epsilon0 <- length(epsilon0)
      n_epsilon1 <- length(epsilon1)
      p0 <- c(
        sum(epsilon0 == 1) / n_epsilon0 + prob.bound,
        sum(epsilon0 == 2) / n_epsilon0 + prob.bound,
        sum(epsilon1 == 1) / n_epsilon1 + prob.bound,
        sum(epsilon1 == 2) / n_epsilon1 + prob.bound
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
  out <- do.call(rbind, pred_list)
  return(out)
}

estimating_equation_no_log_pred <- function(
    pred,
    alpha_tmp_1, beta_tmp_1,
    alpha_tmp_2, beta_tmp_2,
    estimand, optim.method, prob.bound
) {
  objective_function <- function(pred) {
    ret <- numeric(4)

    clamp_p <- function(p) {
      if (p < prob.bound) {
        return(prob.bound)
      } else if ((1 - p) < prob.bound) {
        return(1 - prob.bound)
      } else {
        return(p)
      }
    }
    pred1 <- clamp_p(pred[1])
    pred2 <- clamp_p(pred[2])
    pred3 <- clamp_p(pred[3])
    pred4 <- clamp_p(pred[4])

    if ((1 - pred[1] - pred[2] < prob.bound) |
        (1 - pred[3] - pred[4] < prob.bound)) {
      p0102 <- prob.bound
    } else {
      p01   <- abs(1 - pred1 - pred2)
      p02   <- abs(1 - pred3 - pred4)
      p0102 <- p01*p02
    }

    if (estimand$effect.measure1 == 'RR') {
      ret[1] <- exp(alpha_tmp_1) - pred1*pred3/p0102
      ret[2] <- exp(beta_tmp_1)  - pred3/pred1
    } else if (estimand$effect.measure1 == 'OR') {
      ret[1] <- exp(alpha_tmp_1) - pred1*pred3/p0102
      ret[2] <- exp(beta_tmp_1)  - pred3/(1 - pred3)/pred1*(1 - pred1)
    } else if (estimand$effect.measure1 == 'SHR') {
      ret[1] <- exp(alpha_tmp_1) - (pred1*pred3/p0102)
      ret[2] <- exp(beta_tmp_1) - ( log(1 - pred3) / log(1 - pred1) )
    } else {
      stop("Invalid effect_measure. Must be RR, OR or SHR.")
    }

    if (estimand$effect.measure2 == 'RR') {
      ret[3] <- exp(alpha_tmp_2) - (pred2*pred4/p0102)
      ret[4] <- exp(beta_tmp_2)  - (pred4/pred2)
    } else if (estimand$effect.measure2 == 'OR') {
      ret[3] <- exp(alpha_tmp_2) - (pred2*pred4/p0102)
      ret[4] <- exp(beta_tmp_2)  - pred4/(1 - pred4)/pred2*(1 - pred2)
    } else if (estimand$effect.measure2 == 'SHR') {
      ret[3] <- exp(alpha_tmp_2) - (pred2*pred4/p0102)
      ret[4] <- exp(beta_tmp_2) - ( log(1 - pred4) / log(1 - pred2) )
    } else {
      stop("Invalid effect_measure. Must be RR, OR or SHR.")
    }
    if (optim.method$inner.optim.method == 'multiroot'){
      return(ret)
    } else {
      return(sum(ret^2))
    }
  }
  return(objective_function(pred))
}

estimating_equation_pred_old <- function(
    log_p,
    alpha_tmp_1, beta_tmp_1,
    alpha_tmp_2, beta_tmp_2,
    estimand, optim.method, prob.bound
) {
  objective_function <- function(log_p) {
    ret <- numeric(4)

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
    log_p1 <- clamp_log_p(log_p[1])
    log_p2 <- clamp_log_p(log_p[2])
    log_p3 <- clamp_log_p(log_p[3])
    log_p4 <- clamp_log_p(log_p[4])

    exp_lp1 <- exp(log_p1)
    exp_lp2 <- exp(log_p2)
    exp_lp3 <- exp(log_p3)
    exp_lp4 <- exp(log_p4)

    if ((1 - exp(log_p[1]) - exp(log_p[2]) < prob.bound) |
        (1 - exp(log_p[3]) - exp(log_p[4]) < prob.bound)) {
      lp0102 <- log(prob.bound)
    } else {
      lp01   <- log(abs(1 - exp_lp1 - exp_lp2))
      lp02   <- log(abs(1 - exp_lp3 - exp_lp4))
      lp0102 <- lp01 + lp02
    }

    if (estimand$effect.measure1 == 'RR') {
      ret[1] <- alpha_tmp_1 - log_p1 - log_p3 + lp0102
      ret[2] <- beta_tmp_1  - log_p3 + log_p1
    } else if (estimand$effect.measure1 == 'OR') {
      ret[1] <- alpha_tmp_1 - log_p1 - log_p3 + lp0102
      ret[2] <- beta_tmp_1  - log_p3 + log_p1 +
        log(1 - exp_lp3) - log(1 - exp_lp1)
    } else if (estimand$effect.measure1 == 'SHR') {
      ret[1] <- alpha_tmp_1 - log_p1 - log_p3 + lp0102
      ret[2] <- exp(beta_tmp_1) - ( log(1 - exp_lp3) / log(1 - exp_lp1) )
    } else {
      stop("Invalid effect_measure. Must be RR, OR or SHR.")
    }

    if (estimand$effect.measure2 == 'RR') {
      ret[3] <- alpha_tmp_2 - log_p2 - log_p4 + lp0102
      ret[4] <- beta_tmp_2  - log_p4 + log_p2
    } else if (estimand$effect.measure2 == 'OR') {
      ret[3] <- alpha_tmp_2 - log_p2 - log_p4 + lp0102
      ret[4] <- beta_tmp_2  - log_p4 + log_p2 +
        log(1 - exp_lp4) - log(1 - exp_lp2)
    } else if (estimand$effect.measure2 == 'SHR') {
      ret[3] <- alpha_tmp_2 - log_p2 - log_p4 + lp0102
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



estimating_equation_pred <- function(
    log_p,                 # c(log_p10, log_p20, log_p11, log_p21)
    alpha_tmp_1, beta_tmp_1,
    alpha_tmp_2, beta_tmp_2,
    estimand, optim.method, prob.bound
) {
  objective_function <- function(log_p) {
    ret <- numeric(4)

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
    log_p1 <- clamp_log_p(log_p[1])
    log_p2 <- clamp_log_p(log_p[2])
    log_p3 <- clamp_log_p(log_p[3])
    log_p4 <- clamp_log_p(log_p[4])

    exp_lp1 <- exp(log_p1)
    exp_lp2 <- exp(log_p2)
    exp_lp3 <- exp(log_p3)
    exp_lp4 <- exp(log_p4)

    if ((1 - exp(log_p[1]) - exp(log_p[2]) < prob.bound) |
        (1 - exp(log_p[3]) - exp(log_p[4]) < prob.bound)) {
      lp0102 <- log(prob.bound)
    } else {
      lp01   <- log(abs(1 - exp_lp1 - exp_lp2))
      lp02   <- log(abs(1 - exp_lp3 - exp_lp4))
      lp0102 <- lp01 + lp02
    }

    if (estimand$effect.measure1 == 'RR') {
      ret[1] <- alpha_tmp_1 - log_p1 - log_p3 + lp0102
      ret[2] <- beta_tmp_1  - log_p3 + log_p1
    } else if (estimand$effect.measure1 == 'OR') {
      ret[1] <- alpha_tmp_1 - log_p1 - log_p3 + lp0102
      ret[2] <- beta_tmp_1  - log_p3 + log_p1 +
        log(1 - exp_lp3) - log(1 - exp_lp1)
    } else if (estimand$effect.measure1 == 'SHR') {
      ret[1] <- alpha_tmp_1 - log_p1 - log_p3 + lp0102
      ret[2] <- exp(beta_tmp_1) - ( log(1 - exp_lp3) / log(1 - exp_lp1) )
    } else {
      stop("Invalid effect_measure. Must be RR, OR or SHR.")
    }

    if (estimand$effect.measure2 == 'RR') {
      ret[3] <- alpha_tmp_2 - log_p2 - log_p4 + lp0102
      ret[4] <- beta_tmp_2  - log_p4 + log_p2
    } else if (estimand$effect.measure2 == 'OR') {
      ret[3] <- alpha_tmp_2 - log_p2 - log_p4 + lp0102
      ret[4] <- beta_tmp_2  - log_p4 + log_p2 +
        log(1 - exp_lp4) - log(1 - exp_lp2)
    } else if (estimand$effect.measure2 == 'SHR') {
      ret[3] <- alpha_tmp_2 - log_p2 - log_p4 + lp0102
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

calc_pred_survival <- function(alpha_beta, x_l, offset, estimand) {
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

