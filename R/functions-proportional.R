solveEstimatingEquationP <- function(
    nuisance.model,
    exposure,
    strata,
    sorted_data,
    estimand = estimand,
    optim.method = optim.method,
    data.initial.values,
    prob.bound = 1e-5
) {
  outcome.type <- 'PROPORTIONAL'
  i_time <- 0
  ip.weight.matrix <- matrix(NA, nrow=nrow(sorted_data), ncol=length(estimand$time.point))
  for (specific.time in estimand$time.point) {
    i_time <- i_time + 1
    ip.weight.matrix[,i_time] <- calculateIPCW(formula = nuisance.model, data = sorted_data, code.censoring=estimand$code.censoring, strata=strata, specific.time = specific.time)
  }

  makeObjectiveFunction <- function() {
    out_ipcw <- list()
    initial.CIFs <- NULL
    estimating_equation_p <- function(p) {
      out_ipcw <- estimating_equation_proportional(
        formula = nuisance.model,
        data = sorted_data,
        exposure = exposure,
        ip.weight.matrix = ip.weight.matrix,
        alpha_beta = p,
        estimand = estimand,
        time.point = estimand$time.point,
        optim.method = optim.method,
        prob.bound = prob.bound,
        initial.CIFs = initial.CIFs)
      out_ipcw <<- out_ipcw
      return(out_ipcw$ret)
    }
    setInitialCIFs <- function(new.CIFs) {
      initial.CIFs <<- new.CIFs
    }
    getResults <- function() {
      out_ipcw
    }
    list(
      estimating_equation_p = estimating_equation_p,
      setInitialCIFs = setInitialCIFs,
      getResults = getResults
    )
  }

  iteration <- 0
  max_param_diff  <- Inf
  obj <- makeObjectiveFunction()
  sol_list <- list()
  diff_list <- list()

  current_params <- data.initial.values
  while ((iteration < optim.method$optim.parameter2) & (max_param_diff > optim.method$optim.parameter1)) {
    iteration <- iteration + 1
    if (optim.method$outer.optim.method == "nleqslv" | optim.method$outer.optim.method == "Broyden"){
      sol <- nleqslv(current_params, obj$estimating_equation_p, method="Broyden", control=list(maxit=optim.method$optim.parameter5, allowSingular=FALSE))
      new_params <- sol$x
    } else if (optim.method$optim.method$outer.optim.method == "Newton"){
      sol <- nleqslv(current_params, obj$estimating_equation_p, method="Newton", control=list(maxit=optim.method$optim.parameter5, allowSingular=FALSE))
      new_params <- sol$x
    } else if (optim.method$outer.optim.method == "multiroot") {
      sol <- multiroot(obj$estimating_equation_p, start = current_params, maxiter=optim.method$optim.parameter5, rtol=optim.method$optim.parameter4)
      new_params <- sol$root
    } else if (optim.method$outer.optim.method == "optim" | optim.method$outer.optim.method == "SANN"){
      sol <- optim(par = current_params,
                   fn = function(params) {
                     sum(obj$estimating_equation_p(params)^2)
                   },
                   method = "SANN",  control = list(maxit=optim.method$optim.parameter5, reltol=optim.method$optim.parameter4)
      )
      new_params <- sol$par
    } else if (optim.method$outer.optim.method == "BFGS"){
      sol <- optim(par = current_params,
                   fn = function(params) {
                     sum(obj$estimating_equation_p(params)^2)
                   },
                   method = "BFGS",  control = list(maxit=optim.method$optim.parameter5, reltol=optim.method$optim.parameter4)
      )
      new_params <- sol$par
    }
    if (any(abs(new_params) > optim.method$optim.parameter3)) {
      stop("Estimates are either too large or too small, and convergence might not be achieved.")
    }
    param_diff <- abs(new_params - current_params)
    max_param_diff   <- max(param_diff)
    current_params <- new_params

    obj$setInitialCIFs(obj$getResults()$potential.CIFs)
    sol_list[[iteration]] <- sol
    diff_list[[iteration]] <- max_param_diff
  }
  return(current_params)
}

reportProportional <- function(nuisance.model,
                               exposure,
                               estimand,
                               coef,
                               coef_se,
                               p_value,
                               conf_lower,
                               conf_upper,
                               out_getResults,
                               iteration,
                               param_diff,
                               sol,
                               outer.optim.method) {
  i_parameter <- out_getResults$i_parameter

  coef1 <- list(coef = coef[1], coef_se = coef_se[1], conf_lower = conf_lower[1], conf_upper = conf_upper[1], p_value = p_value[1])
  coef2 <- list(coef = coef[2], coef_se = coef_se[2], conf_lower = conf_lower[2], conf_upper = conf_upper[2], p_value = p_value[2])

  summary_text <- function(events, observations, exposed) {
    paste(events, "events in", observations, if (exposed) "exposed observations" else "unexposed observations")
  }

  text_values <- list(
    exposure_text   = paste(exposure, "( ref =", estimand$code.exposure.ref, ")"),
    effect1_text    = paste(estimand$effect.measure1, "of", exposure, "over time"),
    effect2_text    = paste(estimand$effect.measure2, "of", exposure, "over time"),
    median_followup = paste(round(median(out_getResults$t),digit=2), "[", round(min(out_getResults$t),digit=3), ",", round(max(out_getResults$t),digit=2), "]")
  )

  tidy_df <- function(coef, text) {
    data.frame(
      term = text,
      estimate = coef$coef,
      std.error = coef$coef_se,
      p.value = coef$p_value,
      conf.low = coef$conf_lower,
      conf.high = coef$conf_upper
    )
  }

  glance_df <- function(effect_text, events, exposed_events, unexposed_events, param_len, message) {
    data.frame(
      effect.measure = effect_text,
      n.events = events,
      n.events.exposed = exposed_events,
      n.events.unexposed = unexposed_events,
      median.follow.up = text_values$median_followup,
      n.parameters = param_len,
      optimization.message = message
    )
  }

  getMessage <- function() {
    if (outer.optim.method %in% c("nleqslv", "Newton", "Broyden")) return(sol$message)
    if (outer.optim.method %in% c("optim", "BFGS", "SANN")) return(ifelse(!is.null(sol$message), sol$message, "None"))
    return("-")
  }

  tg <- list(
    event1 = list(
      tidy = tidy_df(coef1, text_values$exposure_text),
      glance = glance_df(text_values$effect1_text,
                         sum(out_getResults$y_1),
                         summary_text(sum(out_getResults$y_1 * out_getResults$x_a), sum(out_getResults$x_a), TRUE),
                         summary_text(sum(out_getResults$y_1 * (1 - out_getResults$x_a)), length(out_getResults$y_1) - sum(out_getResults$x_a), FALSE),
                         i_parameter[3],
                         getMessage())
    ),
    event2 = list(
      tidy = tidy_df(coef2, text_values$exposure_text),
      glance = glance_df(text_values$effect2_text,
                         sum(out_getResults$y_2), summary_text(sum(out_getResults$y_2 * out_getResults$x_a), sum(out_getResults$x_a), TRUE),
                         summary_text(sum(out_getResults$y_2 * (1 - out_getResults$x_a)), length(out_getResults$y_2) - sum(out_getResults$x_a), FALSE),
                         i_parameter[7] - i_parameter[4] + 1,
                         "-")
    )
  )

  class(tg$event1) <- class(tg$event2) <- "modelsummary_list"
  return(tg)
}
