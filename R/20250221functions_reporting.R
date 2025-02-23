reporting_competing_risk <- function(nuisance.model, exposure, estimand, alpha_beta_estimated, cov_estimated, 
                                     objget_results, iteration, diff_val, sol, conf.level, outer.optim.method) {
  alpha <- 1 - conf.level
  critical_value <- qnorm(1 - alpha / 2)
  n_para_1 <- objget_results$n_para_1
  n_para_2 <- objget_results$n_para_2
  n_para_3 <- objget_results$n_para_3
  n_para_4 <- objget_results$n_para_4
  n_para_5 <- objget_results$n_para_5
  rr_1 <- exp(alpha_beta_estimated[n_para_2])
  rr_2 <- exp(alpha_beta_estimated[n_para_5])
  se_1 <- sqrt(cov_estimated[n_para_2, n_para_2])
  se_2 <- sqrt(cov_estimated[n_para_5, n_para_5])
  rr_1_l <- exp(alpha_beta_estimated[n_para_2] - critical_value * se_1)
  rr_1_u <- exp(alpha_beta_estimated[n_para_2] + critical_value * se_1)
  rr_2_l <- exp(alpha_beta_estimated[n_para_5] - critical_value * se_2)
  rr_2_u <- exp(alpha_beta_estimated[n_para_5] + critical_value * se_2)
  z_1 <- abs(alpha_beta_estimated[n_para_2]) / se_1
  z_2 <- abs(alpha_beta_estimated[n_para_5]) / se_2
  p_1 <- 2 * (1 - pnorm(z_1))
  p_2 <- 2 * (1 - pnorm(z_2))
  
  coef_1 <- alpha_beta_estimated[1:n_para_2]
  coef_2 <- alpha_beta_estimated[n_para_3:n_para_5]
  coef_se_1 <- sqrt(diag(cov_estimated)[1:n_para_2])
  coef_se_2 <- sqrt(diag(cov_estimated)[n_para_3:n_para_5])
  coef_1_l <- alpha_beta_estimated[1:n_para_2] - critical_value * coef_se_1
  coef_1_u <- alpha_beta_estimated[1:n_para_2] + critical_value * coef_se_1
  coef_2_l <- alpha_beta_estimated[n_para_3:n_para_5] - critical_value * coef_se_2
  coef_2_u <- alpha_beta_estimated[n_para_3:n_para_5] + critical_value * coef_se_2
  coef_z_1 <- abs(alpha_beta_estimated[1:n_para_2]) / coef_se_1
  coef_z_2 <- abs(alpha_beta_estimated[n_para_3:n_para_5]) / coef_se_2
  coef_p_1 <- 2 * (1 - pnorm(coef_z_1))
  coef_p_2 <- 2 * (1 - pnorm(coef_z_2))
  
  t <- objget_results$t
  y_0_ <- objget_results$y_0_
  y_1_ <- objget_results$y_1_
  y_2_ <- objget_results$y_2_
  y_0 <- objget_results$y_0
  y_1 <- objget_results$y_1
  y_2 <- objget_results$y_2
  
  text1 <- c("Intercept", attr(terms(nuisance.model), "term.labels"), exposure)
  text2 <- paste(estimand$effect.measure1, "of", exposure, "at", estimand$time.point)
  text3 <- paste(estimand$effect.measure2, "of", exposure, "at", estimand$time.point)
  text4 <- paste(sum(y_1), "events in", length(y_1), "observations")
  text5 <- paste(sum(y_2), "events in", length(y_2), "observations")
  text6 <- paste(median(t), "[", min(t), ",", max(t), "]")
  
  ti_1 <- data.frame(
    term = text1,
    estimate = coef_1,  
    std.error = coef_se_1,
    p.value = coef_p_1,  
    conf.low = coef_1_l,  
    conf.high = coef_1_u
  )
  ti_2 <- data.frame(
    term = text1,
    estimate = coef_2,  
    std.error = coef_se_2,
    p.value = coef_p_2,  
    conf.low = coef_2_l,  
    conf.high = coef_2_u
  )
  if (outer.optim.method == 'nleqslv' || outer.optim.method == 'Newton' || outer.optim.method == 'Broyden' || outer.optim.method == 'partial'){
    gl_1 <- data.frame(
      effect.measure = text2, 
      n.events = text4,
      median.follow.up = text6,
      #    n.parameter = length(alpha_beta_estimated[1:n_para_2]),
      #    min.parameter = min(alpha_beta_estimated[1:n_para_2]),
      #    max.parameter = max(alpha_beta_estimated[1:n_para_2]),
      n.loop.iteration = iteration,  
      #    max.estimate.difference = diff_val, 
      n.optimization.iteration = sol$iter, 
      max.function.value = max(sol$fvec),
      optimization.message = sol$message
    )
    gl_2 <- data.frame(
      effect.measure = text3, 
      n.events = text5,
      median.follow.up = text6,
      #    n.parameter = length(alpha_beta_estimated[n_para_3:n_para_5]),
      #    min.parameter = min(alpha_beta_estimated[n_para_3:n_para_5]),
      #    max.parameter = max(alpha_beta_estimated[n_para_3:n_para_5]),
      n.loop.iteration = "-",  
      #    max.estimate.difference = diff_val, 
      n.optimization.iteration = "-", 
      max.function.value = "-",
      optimization.message = "-"
    )
  } else if (outer.optim.method == 'optim' || outer.optim.method == 'BFGS' || outer.optim.method == 'SANN') {
    gl_1 <- data.frame(
      effect.measure = text2, 
      n.events = text4,
      median.follow.up = text6,
      #n.parameter = length(alpha_beta_estimated[1:n_para_2]),
      #min.parameter = min(alpha_beta_estimated[1:n_para_2]),
      #max.parameter = max(alpha_beta_estimated[1:n_para_2]),
      n.loop.iteration = iteration,  
      #max.estimate.difference = diff_val, 
      n.optimization.iteration = as.vector(sol$counts)[2], 
      max.function.value = max(sol$value),
      optimization.message = if (!is.null(sol$message)) sol$message else "None"
    )
    gl_2 <- data.frame(
      effect.measure = text3, 
      n.events = text5,
      median.follow.up = text6,
      #n.parameter = length(alpha_beta_estimated[n_para_3:n_para_5]),
      #min.parameter = min(alpha_beta_estimated[n_para_3:n_para_5]),
      #max.parameter = max(alpha_beta_estimated[n_para_3:n_para_5]),
      n.loop.iteration = "-",  
      #max.estimate.difference = diff_val, 
      n.optimization.iteration = "-", 
      max.function.value = "-",
      optimization.message = "-"
    )
  } else {
    gl_1 <- data.frame(
      effect.measure = text2, 
      n.events = text4,
      median.follow.up = text6,
      #n.parameter = length(alpha_beta_estimated[1:n_para_2]),
      #min.parameter = min(alpha_beta_estimated[1:n_para_2]),
      #max.parameter = max(alpha_beta_estimated[1:n_para_2]),
      n.loop.iteration = iteration, 
      #max.estimate.difference = diff_val,
      n.optimization.iteration = sol$iter, 
      max.function.value = max(sol$f.root),
      optimization.message = "-"
    )
    gl_2 <- data.frame(
      effect.measure = text3, 
      n.events = text5,
      median.follow.up = text6,
      #    n.parameter = length(alpha_beta_estimated[n_para_3:n_para_5]),
      #    min.parameter = min(alpha_beta_estimated[n_para_3:n_para_5]),
      #    max.parameter = max(alpha_beta_estimated[n_para_3:n_para_5]),
      n.loop.iteration = "-", 
      #    max.estimate.difference = diff_val,
      n.optimization.iteration = "-", 
      max.function.value = "-",
      optimization.message = "-"
    )
  }
  tg1 <- list(tidy = ti_1,glance = gl_1)
  tg2 <- list(tidy = ti_2,glance = gl_2)
  class(tg1) <- "modelsummary_list"
  class(tg2) <- "modelsummary_list"
  tg <- list()
  tg[[1]] <- tg1
  tg[[2]] <- tg2
  names(tg) <- c("event 1", "event 2")
  return(tg)
}

reporting_survival <- function(nuisance.model, exposure, estimand, alpha_beta_estimated, cov_estimated, 
                               objget_results, iteration, diff_val, sol, conf.level, outer.optim.method) {
  alpha <- 1 - conf.level
  critical_value <- qnorm(1 - alpha / 2)
  n_para_1 <- objget_results$n_para_1
  n_para_2 <- objget_results$n_para_2
  n_para_3 <- objget_results$n_para_3
  rr_1 <- exp(alpha_beta_estimated[n_para_2])
  se_1 <- sqrt(cov_estimated[n_para_2, n_para_2])
  rr_1_l <- exp(alpha_beta_estimated[n_para_2] - critical_value * se_1)
  rr_1_u <- exp(alpha_beta_estimated[n_para_2] + critical_value * se_1)
  z_1 <- abs(alpha_beta_estimated[n_para_2]) / se_1
  p_1 <- 2 * (1 - pnorm(z_1))
  
  coef_1 <- alpha_beta_estimated[1:n_para_2]
  coef_se_1 <- sqrt(diag(cov_estimated)[1:n_para_2])
  coef_1_l <- alpha_beta_estimated[1:n_para_2] - critical_value * coef_se_1
  coef_1_u <- alpha_beta_estimated[1:n_para_2] + critical_value * coef_se_1
  coef_z_1 <- abs(alpha_beta_estimated[1:n_para_2]) / coef_se_1
  coef_p_1 <- 2 * (1 - pnorm(coef_z_1))
  
  t <- objget_results$t
  y_0_ <- objget_results$y_0_
  y_1_ <- objget_results$y_1_
  y_0 <- objget_results$y_0
  y_1 <- objget_results$y_1
  
  text1 <- c("Intercept", attr(terms(nuisance.model), "term.labels"), exposure)
  text2 <- paste(estimand$effect.measure1, "of", exposure, "at", estimand$time.point)
  text4 <- paste(sum(y_1), "events in", length(y_1), "observations")
  text6 <- paste(median(t), "[", min(t), ",", max(t), "]")
  
  ti_1 <- data.frame(
    term = text1,
    estimate = coef_1,  
    std.error = coef_se_1,
    p.value = coef_p_1,  
    conf.low = coef_1_l,  
    conf.high = coef_1_u
  )

  if (outer.optim.method == 'nleqslv' || outer.optim.method == 'Newton' || outer.optim.method == 'Broyden' || outer.optim.method == 'partial'){
    gl_1 <- data.frame(
      effect.measure = text2, 
      n.events = text4,
      median.follow.up = text6,
      #    n.parameter = length(alpha_beta_estimated[1:n_para_2]),
      #    min.parameter = min(alpha_beta_estimated[1:n_para_2]),
      #    max.parameter = max(alpha_beta_estimated[1:n_para_2]),
      n.loop.iteration = iteration,  
      #    max.estimate.difference = diff_val, 
      #    max.function.value = max(sol$fvec),
      n.optimization.iteration = sol$iter, 
      optimization.message = sol$message
    )
  } else if (outer.optim.method == 'optim' || outer.optim.method == 'BFGS' || outer.optim.method == 'SANN') {
    gl_1 <- data.frame(
      effect.measure = text2, 
      n.events = text4,
      median.follow.up = text6,
      #n.parameter = length(alpha_beta_estimated[1:n_para_2]),
      #min.parameter = min(alpha_beta_estimated[1:n_para_2]),
      #max.parameter = max(alpha_beta_estimated[1:n_para_2]),
      n.loop.iteration = iteration,  
      #max.estimate.difference = diff_val, 
      #max.function.value = max(sol$value),
      n.solver.iteration = as.vector(sol$counts)[2], 
      optimization.message = if (!is.null(sol$message)) sol$message else "None"
    )
  } else {
    gl_1 <- data.frame(
      effect.measure = text2, 
      n.events = text4,
      median.follow.up = text6,
      #n.parameter = length(alpha_beta_estimated[1:n_para_2]),
      #min.parameter = min(alpha_beta_estimated[1:n_para_2]),
      #max.parameter = max(alpha_beta_estimated[1:n_para_2]),
      n.loop.iteration = iteration, 
      #max.estimate.difference = diff_val,
      #max.function.value = max(sol$f.root),
      n.optimization.iteration = sol$iter, 
      optimization.message = "-"
    )
  }
  
  tg1 <- list(tidy = ti_1,glance = gl_1)
  class(tg1) <- "modelsummary_list"
  tg <- list()
  tg[[1]] <- tg1
  names(tg) <- c("event 1 (no competing risk)")
  return(tg)
}

reporting_prediction <- function(formula,data,exposure,alpha_beta,cov,outcome.type,estimand,optim.method,prob.bound,initial_pred = NULL) 
{
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
  
  a_ <- as.factor(data[[exposure]])
  a <- model.matrix(~ a_)[, 2]
  epsilon0 <- epsilon[a == 0]
  epsilon1 <- epsilon[a == 1]
  x_l <- model.matrix(Terms, mf)
  n_para_1 <- ncol(x_l)
  one <- rep(1, nrow(x_l))
  if (outcome.type == 'SURVIVAL') {
    pred <- calc_pred_survival(alpha_beta,x_l,offset,estimand)
  } else {
    pred <- calc_pred(alpha_beta,x_l,offset,epsilon0,epsilon1,one,n_para_1,estimand,optim.method,prob.bound,initial_pred)
  }
  return(pred)
}
