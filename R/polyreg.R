#' @title Direct polynomial regression for survival and competing risks analysis
#'
#' @param nuisance.model formula Model formula representing outcome, exposure and covariates
#' @param exposure character Column name representing the exposure (1 = exposed, 0 = not exposed).
#' @param cens.model formula Model formula representing censoring and covariates
#' @param data data.frame Input dataset containing survival data.
#' @param effect.measure1 character Specifies the effect measure for event (RR, OR, SHR).
#' @param effect.measure2 character Specifies the effect measure for competing risk (RR, OR, SHR).
#' @param time.point numeric The time point(s) for analysis.
#' @param outcome.type character Specifies the type of outcome (COMPETINGRISK or SURVIVAL).
#' @param conf.level numeric The level for confidence intervals.
#' @param outer.optim.method character Specifies the method of optimization (nleqslv, optim, multiroot).
#' @param inner.optim.method character Specifies the method of optimization (nleqslv, optim, multiroot).
#' @param optim.parameter1 numeric Convergence threshold for outer loop. Defaults to 1e-5.
#' @param optim.parameter2 integer Maximum number of iterations. Defaults to 100.
#' @param optim.parameter3 numeric
#' @param optim.parameter4 numeric Constraint range for parameters. Defaults to 20.
#' @param optim.parameter5 numeric Convergence threshold for optim in outer loop. Defaults to 1e-5.
#' @param optim.parameter6 integer Maximum number of iterations for nleqslv/optim in outer loop. Defaults to 200.
#' @param optim.parameter7 numeric Convergence threshold for optim in inner loop. Defaults to 1e-7.
#' @param optim.parameter8 integer Maximum number of iterations for optim in inner loop. Defaults to 200.
#' @param data.initlal.values data.frame Optional dataset containing initial values. Defaults to NULL.
#' @param covariate.normalization logical Indicates whether covariates are normalized (TRUE = normalize, FALSE = otherwise). Defaults to TRUE.
#' @param sort.data logical Indicates whether data are initially sorted to reduce computation steps (TRUE = sort, FALSE = otherwise). Defaults to TRUE.
#' @param prob.bound numeric Small threshold for clamping probabilities. Defaults to 1e-5.
#'
#' @return A list of results from direct polynomial regression. summary meets requirement of msummary function.
#' @export
#'
#' @examples
#' data(bmt)
#' result <- polyreg(nuisance.model = Event(time, cause)~age+tcell, exposure = 'platelet',
#' cens.model = Event(time,cause==0)~+1, data = bmt, effect.measure1='RR', effect.measure2='RR', time.point=24, outcome.type='COMPETINGRISK')
#' msummary(result$out_summary, statistic = c("conf.int"), exponentiate = TRUE)
polyreg <- function(
    nuisance.model,
    exposure,
    cens.model,
    data,
    effect.measure1 = 'RR',
    effect.measure2 = 'RR',
    time.point,
    outcome.type = 'COMPETINGRISK',
    conf.level = 0.95,
    outer.optim.method = 'nleqslv',
    inner.optim.method = 'optim',
    optim.parameter1 = 1e-5,
    optim.parameter2 = 100,
    optim.parameter3 = 1e-5, # Not used now
    optim.parameter4 = 20,
    optim.parameter5 = 1e-5,
    optim.parameter6 = 200,
    optim.parameter7 = 1e-7,
    optim.parameter8 = 200,
    data.initlal.values = NULL,
    covariate.normalization = TRUE,
    sort.data = TRUE,
    prob.bound = 1e-5
) {

  #######################################################################################
  # 1. Pre-processing (function: normalize_covariate, sort_by_covariate)
  #######################################################################################
  if (outcome.type %in% c("COMPETINGRISK", "C", "CR", 'COMPETING RISK', 'COMPETINGRISKS', 'COMPETING RISKS', 'Competingrisk', 'Competing risk', 'Competingrisks', 'Competing risks', 'competingrisk', 'competing risk', 'competingrisks', 'competing risks')) {
    outcome.type <- 'COMPETINGRISK'
    if (length(time.point)>1) {time.point <- max(time.point)}
  } else if (outcome.type %in% c("SURVIVAL", "S", "Survival", "Survival")) {
    outcome.type <- 'SURVIVAL'
    if (length(time.point)>1) {time.point <- max(time.point)}
  } else if (outcome.type %in% c("PROPORTIONAL", "P", "Proportional", "proportional")) {
    outcome.type <- 'PROPORTIONAL'
    if (length(time.point)==1) {outcome.type <- 'COMPETINGRISK'}
  } else {
    stop("Invalid input for 'outcome.type', Choose 'COMPETINGRISK', 'SURVIVAL', or 'PROPORTIONAL'.")
  }
  if (conf.level <= 0 || conf.level >= 1)
    stop("Confidence level must be between 0 and 1")

  estimand <- list(effect.measure1=effect.measure1, effect.measure2=effect.measure2, time.point=time.point)
  optim.method <- list(
    outer.optim.method = outer.optim.method,
    inner.optim.method = inner.optim.method,
    optim.parameter1 = optim.parameter1,
    optim.parameter2 = optim.parameter2,
    optim.parameter3 = optim.parameter3,
    optim.parameter4 = optim.parameter4,
    optim.parameter5 = optim.parameter5,
    optim.parameter6 = optim.parameter6,
    optim.parameter7 = optim.parameter7,
    optim.parameter8 = optim.parameter8
  )

  out_normalize_covariate <- normalize_covariate(nuisance.model, data, covariate.normalization, outcome.type)
  normalized_data <- out_normalize_covariate$normalized_data
  sorted_data <- sort_by_covariate(nuisance.model, normalized_data, sort.data, out_normalize_covariate$n_covariate)

  #######################################################################################
  # 2. Pre-processing and Calculating initial values alpha_beta_0 (function: calc_initial_values)
  #######################################################################################
  if (outcome.type == 'PROPORTIONAL') {
    n_para_1 <- out_normalize_covariate$n_covariate+1
    n_para_2 <- out_normalize_covariate$n_covariate+2
    n_para_3 <- out_normalize_covariate$n_covariate+3
    n_para_4 <- 2*out_normalize_covariate$n_covariate+3
    n_para_5 <- 2*out_normalize_covariate$n_covariate+4
    n_para_6 <- length(time.point)*(n_para_5-2) + 2
    alpha_beta_0 <- rep(NA, n_para_6) # initial parameters for event 1 and 2 over time points
    i_time <- 0
    sum1 <- 0
    sum2 <- 0
    for (specific.time in time.point) {
      specific.time <- as.numeric(specific.time)
      out_calc_initial_values <- calc_initial_values(
        formula = nuisance.model,
        data = sorted_data,
        exposure = exposure,
        data.initlal.values = data.initlal.values,
        estimand = estimand,
        specific.time = specific.time,
        outcome.type = outcome.type,
        optim.method = optim.method,
        prob.bound = prob.bound
      )
      i_time <- i_time + 1
      i_para <- n_para_1*(i_time-1)+1
      tmp1 <- out_calc_initial_values$init_vals[1:n_para_1]
      tmp2 <- out_calc_initial_values$init_vals[n_para_3:n_para_4]
      alpha_beta_0[i_para:(i_para+n_para_1-1)]                          <- tmp1
      alpha_beta_0[(n_para_6/2+i_para):(n_para_6/2+i_para+n_para_1-1)]  <- tmp2
      sum1 <- sum1 + out_calc_initial_values$init_vals[n_para_2]
      sum2 <- sum2 + out_calc_initial_values$init_vals[n_para_5]
    }
    alpha_beta_0[n_para_6/2]  <- sum1/length(time.point)
    alpha_beta_0[n_para_6]    <- sum2/length(time.point)
  } else {
    out_calc_initial_values <- calc_initial_values(
      formula = nuisance.model,
      data = sorted_data,
      exposure = exposure,
      data.initlal.values = data.initlal.values,
      estimand = estimand,
      specific.time = estimand$time.point,
      outcome.type = outcome.type,
      optim.method = optim.method,
      prob.bound = prob.bound
    )
    alpha_beta_0 <- out_calc_initial_values$init_vals
  }

  #######################################################################################
  # 3. Calculating IPCW (function: IP_WEIGHT_1)
  #######################################################################################
  if (outcome.type == 'PROPORTIONAL') {
    i_time <- 0
    ip_wt_matrix <- matrix(NA, nrow=out_normalize_covariate$n, ncol=length(time.point))
    for (specific.time in time.point) {
      i_time <- i_time + 1
      ip_wt_matrix[,i_time] <- calc_ipw(formula = cens.model, data = sorted_data, specific.time = specific.time)
    }
  } else {
    ip_wt <- calc_ipw(formula = cens.model, data = sorted_data, specific.time = estimand$time.point)
  }

  #######################################################################################
  # 4. Parameter estimation (functions: estimating_equation_ipcw, _partial, _survival)
  #######################################################################################

  makeObjectiveFunction <- function() {
    out_ipcw <- list()
    initial_pred <- NULL
    estimating_equation_i <- function(p) {
      out_ipcw <- estimating_equation_ipcw(
        formula = nuisance.model,
        data = sorted_data,
        exposure = exposure,
        ip_weight = ip_wt,
        alpha_beta = p,
        estimand = estimand,
        optim.method = optim.method,
        prob.bound = prob.bound,
        initial_pred = initial_pred)
      out_ipcw <<- out_ipcw
      return(out_ipcw$ret)
    }
    estimating_equation_s <- function(p) {
      out_ipcw <- estimating_equation_survival(
        formula = nuisance.model,
        data = sorted_data,
        exposure = exposure,
        ip_weight = ip_wt,
        alpha_beta = p,
        estimand = estimand,
        prob.bound = prob.bound,
        initial_pred = initial_pred)
      out_ipcw <<- out_ipcw
      return(out_ipcw$ret)
    }
    estimating_equation_pa <- function(p) {
      out_ipcw <- estimating_equation_partial(
        formula = nuisance.model,
        data = sorted_data,
        exposure = exposure,
        ip_weight = ip_wt,
        alpha_beta_partial = p,
        alpha_beta1_current = current_params1,
        alpha_beta2_current = current_params2,
        optimizing_event_1 = optimizing_event_1,
        estimand = estimand,
        optim.method = optim.method,
        prob.bound = prob.bound,
        initial_pred = initial_pred)
      out_ipcw <<- out_ipcw
      return(out_ipcw$ret)
    }
    estimating_equation_p <- function(p) {
      out_ipcw <- estimating_equation_proportional(
        formula = nuisance.model,
        data = sorted_data,
        exposure = exposure,
        ip_wt_matrix = ip_wt_matrix,
        alpha_beta = p,
        estimand = estimand,
        time.point = time.point,
        optim.method = optim.method,
        prob.bound = prob.bound,
        initial_pred = initial_pred)
      out_ipcw <<- out_ipcw
      return(out_ipcw$ret)
    }
    set_initial_pred <- function(new_pred) {
      initial_pred <<- new_pred
    }
    get_results <- function() {
      out_ipcw
    }
    list(
      estimating_equation_i = estimating_equation_i,
      estimating_equation_s = estimating_equation_s,
      estimating_equation_p = estimating_equation_p,
      estimating_equation_pa = estimating_equation_pa,
      set_initial_pred = set_initial_pred,
      get_results = get_results
    )
  }

  iteration <- 0
  diff_val  <- Inf
  obj <- makeObjectiveFunction()
  sol_list <- list()
  diff_list <- list()

  if (outcome.type == 'COMPETINGRISK' && !outer.optim.method == "partial") {
    current_params <- alpha_beta_0
    while ((iteration < optim.parameter2) && (diff_val > optim.parameter1)) {
      iteration <- iteration + 1
      if (outer.optim.method == "nleqslv" || outer.optim.method == "Broyden"){
        sol <- nleqslv(current_params, obj$estimating_equation_i, method="Broyden", control=list(maxit=optim.parameter6, allowSingular=FALSE))
        new_params <- sol$x
      } else if (outer.optim.method == "Newton"){
        sol <- nleqslv(current_params, obj$estimating_equation_i, method="Newton", control=list(maxit=optim.parameter6, allowSingular=FALSE))
        new_params <- sol$x
      } else if (outer.optim.method == "multiroot") {
        sol <- multiroot(obj$estimating_equation_i, start = current_params, maxiter = optim.parameter6, rtol = optim.parameter5)
        new_params <- sol$root
      } else if (outer.optim.method == "optim" || outer.optim.method == "SANN"){
        sol <- optim(par = current_params,
                     fn = function(params) {
                       sum(obj$estimating_equation_i(params)^2)
                     },
                     method = "SANN",  control = list(maxit = optim.parameter6, reltol = optim.parameter5)
        )
        new_params <- sol$par
      } else if (outer.optim.method == "BFGS"){
        sol <- optim(par = current_params,
                     fn = function(params) {
                       sum(obj$estimating_equation_i(params)^2)
                     },
                     method = "BFGS",  control = list(maxit = optim.parameter6, reltol = optim.parameter5)
        )
        new_params <- sol$par
      }
      if (any(abs(new_params) > optim.parameter4)) {
        stop("Estimates are either too large or too small, and convergence might not be achieved.")
      }
      param_diff <- abs(new_params - current_params)
      diff_val   <- max(param_diff)
      current_params <- new_params

      obj$set_initial_pred(obj$get_results()$pred)
      sol_list[[iteration]] <- sol
      diff_list[[iteration]] <- diff_val
    }
  } else if (outer.optim.method == "partial") {
    n_para_2  <- length(alpha_beta_0)/2
    n_para_3  <- length(alpha_beta_0)/2 + 1
    n_para_5  <- length(alpha_beta_0)
    current_params1 <- alpha_beta_0[1:n_para_2]
    current_params2 <- alpha_beta_0[n_para_3:n_para_5]
    while ((iteration < optim.parameter2) && (diff_val > optim.parameter1)) {
      iteration <- iteration + 1
      optimizing_event_1 <- TRUE # Optimize alpha_beta[1:n_para_2], parameters of event j=1
      sol <- nleqslv(current_params1, obj$estimating_equation_pa, method="Broyden", control=list(maxit=optim.parameter6, allowSingular=FALSE))
      new_params1 <- sol$x
      param_diff1 <- abs(new_params1 - current_params1)
      current_params1 <- new_params1
      optimizing_event_1 <- FALSE # Optimize alpha_beta[n_para_3:n_para_5], parameters of event j=2
      sol <- nleqslv(current_params2, obj$estimating_equation_pa, method="Broyden", control=list(maxit=optim.parameter6, allowSingular=FALSE))
      new_params2 <- sol$x
      param_diff2 <- abs(new_params2 - current_params2)
      current_params2 <- new_params2
      if (any(abs(new_params1) > optim.parameter4 | abs(new_params2) > optim.parameter4)) {
        stop("Estimates are either too large or too small, and convergence might not be achieved.")
      }
      diff_val   <- max(param_diff1, param_diff2)
      obj$set_initial_pred(obj$get_results()$pred)
      sol_list[[iteration]] <- sol
      diff_list[[iteration]] <- diff_val
    }
    current_params <- cbind(current_params1, current_params2)
  } else if (outcome.type == 'SURVIVAL') {
    current_params <- alpha_beta_0
    while ((iteration < optim.parameter2) && (diff_val > optim.parameter1)) {
      iteration <- iteration + 1
      if (outer.optim.method == "nleqslv" || outer.optim.method == "Broyden"){
        sol <- nleqslv(current_params, obj$estimating_equation_s, method="Broyden", control=list(maxit=optim.parameter6, allowSingular=FALSE))
        new_params <- sol$x
      } else if (outer.optim.method == "Newton"){
        sol <- nleqslv(current_params, obj$estimating_equation_s, method="Newton", control=list(maxit=optim.parameter6, allowSingular=FALSE))
        new_params <- sol$x
      } else if (outer.optim.method == "multiroot") {
        sol <- multiroot(obj$estimating_equation_s, start = current_params, maxiter = optim.parameter6, rtol = optim.parameter5)
        new_params <- sol$root
      } else if (outer.optim.method == "optim" || outer.optim.method == "SANN"){
        sol <- optim(par = current_params,
                     fn = function(params) {
                       sum(obj$estimating_equation_s(params)^2)
                     },
                     method = "SANN",  control = list(maxit = optim.parameter6, reltol = optim.parameter5)
        )
        new_params <- sol$par
      } else if (outer.optim.method == "BFGS"){
        sol <- optim(par = current_params,
                     fn = function(params) {
                       sum(obj$estimating_equation_s(params)^2)
                     },
                     method = "BFGS",  control = list(maxit = optim.parameter6, reltol = optim.parameter5)
        )
        new_params <- sol$par
      }
      if (any(abs(new_params) > optim.parameter4)) {
        stop("Estimates are either too large or too small, and convergence might not be achieved.")
      }
      param_diff <- abs(new_params - current_params)
      diff_val   <- max(param_diff)
      current_params <- new_params

      obj$set_initial_pred(obj$get_results()$pred)
      sol_list[[iteration]] <- sol
      diff_list[[iteration]] <- diff_val
    }
  } else if (outcome.type == 'PROPORTIONAL') {
    current_params <- alpha_beta_0
    while ((iteration < optim.parameter2) && (diff_val > optim.parameter1)) {
      iteration <- iteration + 1
      if (outer.optim.method == "nleqslv" || outer.optim.method == "Broyden"){
        sol <- nleqslv(current_params, obj$estimating_equation_p, method="Broyden", control=list(maxit=optim.parameter6, allowSingular=FALSE))
        new_params <- sol$x
      } else if (outer.optim.method == "Newton"){
        sol <- nleqslv(current_params, obj$estimating_equation_p, method="Newton", control=list(maxit=optim.parameter6, allowSingular=FALSE))
        new_params <- sol$x
      } else if (outer.optim.method == "multiroot") {
        sol <- multiroot(obj$estimating_equation_p, start = current_params, maxiter = optim.parameter6, rtol = optim.parameter5)
        new_params <- sol$root
      } else if (outer.optim.method == "optim" || outer.optim.method == "SANN"){
        sol <- optim(par = current_params,
                     fn = function(params) {
                       sum(obj$estimating_equation_p(params)^2)
                     },
                     method = "SANN",  control = list(maxit = optim.parameter6, reltol = optim.parameter5)
        )
        new_params <- sol$par
      } else if (outer.optim.method == "BFGS"){
        sol <- optim(par = current_params,
                     fn = function(params) {
                       sum(obj$estimating_equation_p(params)^2)
                     },
                     method = "BFGS",  control = list(maxit = optim.parameter6, reltol = optim.parameter5)
        )
        new_params <- sol$par
      }
      if (any(abs(new_params) > optim.parameter4)) {
        stop("Estimates are either too large or too small, and convergence might not be achieved.")
      }
      param_diff <- abs(new_params - current_params)
      diff_val   <- max(param_diff)
      current_params <- new_params

      obj$set_initial_pred(obj$get_results()$pred)
      sol_list[[iteration]] <- sol
      diff_list[[iteration]] <- diff_val
    }
    alpha_beta_estimated <- current_params
    return(alpha_beta_estimated)
  }

  #######################################################################################
  # 5. Calculating variance (functions: calc_cov, calc_cov_survival)
  #######################################################################################
  objget_results <- obj$get_results()
  normalize_estimate <- function(out_calc_cov, covariate.normalization, current_params, out_normalize_covariate) {
    if (covariate.normalization == FALSE) {
      alpha_beta_estimated <- current_params
      cov_estimated <- out_calc_cov$cov_estimated
    } else {
      adj <- 1 / as.vector(out_normalize_covariate$range)
      alpha_beta_estimated <- adj * current_params
      adj_matrix <- diag(adj, length(adj), length(adj))
      cov_estimated <- adj_matrix %*% out_calc_cov$cov_estimated %*% adj_matrix
    }
    return(list(alpha_beta_estimated = alpha_beta_estimated, cov_estimated = cov_estimated))
  }

  # Main calculation process for each outcome.type
  if (outcome.type == 'COMPETINGRISK' && outer.optim.method == 'partial') {
    out_calc_cov <- calc_cov(objget_results, estimand, prob.bound)
    out_normalize_estimate <- normalize_estimate(out_calc_cov, covariate.normalization, cbind(current_params1, current_params2), out_normalize_covariate)
  } else if (outcome.type == 'COMPETINGRISK') {
    out_calc_cov <- calc_cov(objget_results, estimand, prob.bound)
    out_normalize_estimate <- normalize_estimate(out_calc_cov, covariate.normalization, current_params, out_normalize_covariate)
  } else if (outcome.type == 'SURVIVAL') {
    out_calc_cov <- calc_cov_survival(objget_results, estimand, prob.bound)
    out_normalize_estimate <- normalize_estimate(out_calc_cov, covariate.normalization, current_params, out_normalize_covariate)
  }

  alpha_beta_estimated <- out_normalize_estimate$alpha_beta_estimated
  cov_estimated <- out_normalize_estimate$cov_estimated

  #######################################################################################
  # 6. Output (functions: reporting_survival, _competing_risk, _prediction)
  #######################################################################################
  report_results <- list(
    SURVIVAL = reporting_survival,
    COMPETINGRISK = reporting_competing_risk
  )
  if (outcome.type %in% names(report_results)) {
    out_summary <- report_results[[outcome.type]](
      nuisance.model, exposure, estimand, alpha_beta_estimated,
      cov_estimated, objget_results, iteration, diff_val, sol,
      conf.level, optim.method$outer.optim.method
    )
    out_prediction <- reporting_prediction(nuisance.model,data,exposure,
                                           alpha_beta_estimated,cov_estimated,outcome.type,estimand,optim.method,prob.bound)
  }
  out_reporting <- list(out_summary = out_summary, out_prediction = out_prediction, out_coefficient=alpha_beta_estimated, out_cov=cov_estimated)
  return(out_reporting)
}

