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


calc_initial_values <- function(
    formula, data, exposure, data.initlal.values, estimand, specific.time, outcome.type, optim.method, prob.bound
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
    check_input(epsilon, outcome.type, optim.method)
  } else {
    stop("Invalid outcome variables. Must be survival or competing risks outcome")
  }
  if (!is.null(offsetpos <- attributes(Terms)$specials$offset)) {
    ts <- survival::untangle.specials(Terms, "offset")
    Terms <- Terms[-ts$terms]
    offset <- m[[ts$vars]]
  } else {
    offset <- rep(0, length(t))
  }
  if (!nrow(data) == nrow(mf))
    stop("Variables contain NA values")
  if (any(is.na(data[[exposure]])))
    stop("Variables contain NA values")
  y_0 <- ifelse(epsilon == 0 | t > specific.time, 1, 0)
  y_1 <- ifelse(epsilon == 1 & t <= specific.time, 1, 0)
  y_2 <- ifelse(epsilon == 2 & t <= specific.time, 1, 0)

  y_0_ <- ifelse(epsilon == 0, 1, 0) # 打ち切り変数
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

  tmp <- cbind(x_l, a)
  zero <- matrix(0, nrow=nrow(tmp), ncol=ncol(tmp))
  tmp1 <- cbind(tmp, zero)
  tmp2 <- cbind(zero, tmp)
  x_la <- rbind(tmp1, tmp2) # d used for calculation of score

  if (!is.null(data.initlal.values)) {
    x_l <- model.matrix(Terms, mf)
    if (!(1+ncol(x_l))*2 == length(data.initlal.values))
      stop("Invalid initial value dataset. Must contain the same number of initial values as parameters")
    out <- list(init_vals = data.initlal.values,
                y_0=y_0, y_1=y_1, y_2=y_2,
                y_0_=y_0_, y_1_=y_1_, y_2_=y_2_,
                a=a, t=t, epsilon=epsilon, epsilon0=epsilon0, epsilon1=epsilon1, x_l=x_l, x_la=x_la, offset=offset,
                n_para_1=n_para_1, n_para_2=n_para_2, n_para_3=n_para_3, n_para_4=n_para_4, n_para_5=n_para_5, one=one)
    return(out)
  }

  binarize_if_continuous <- function(x) {
    if (outcome.type == 'SURVIVAL') {
      if (is.numeric(x) && length(unique(x)) > 2) {  # Check if numeric & not binary
        l <- as.numeric((x >= median(x)) == TRUE)
        out <- calc_initial_1_survival(t, epsilon, a, l, estimand, specific.time, prob.bound)
        return(out)
      } else if (length(unique(x)) == 2) {
        l <- x
        out <- calc_initial_1_survival(t, epsilon, a, l, estimand, specific.time, prob.bound)
        return(out)
      }
    } else {
      if (is.numeric(x) && length(unique(x)) > 2) {  # Check if numeric & not binary
        l <- as.numeric((x >= median(x)) == TRUE)
        out <- calc_initial_1_competing_risk(t, epsilon, a, l, estimand, specific.time, prob.bound)
        return(out)
      } else if (length(unique(x)) == 2) {
        l <- x
        out <- calc_initial_1_competing_risk(t, epsilon, a, l, estimand, specific.time, prob.bound)
        return(out)
      }
    }
  }

  if (all(epsilon %in% c(0, 1))) {
    if (n_para_1>1) {
      out_bic_1 <- t(binarize_if_continuous(x_l[,2])[1,1:2])
      if (n_para_1>2) {
        for (i_para in 3:n_para_1) {
          out_bic_i <- binarize_if_continuous(x_l[,i_para])
          out_bic_1 <- cbind(out_bic_1, t(out_bic_i[1,2]))
        }
        out_bic_1 <- cbind(out_bic_1, t(out_bic_i[1,3]))
      } else {
        out_bic_1 <- cbind(out_bic_1, t(binarize_if_continuous(x_l[,2])[1,3]))
      }
      init_vals <- out_bic_1
    } else {
      init_vals <- calc_initial_2_survival(t, epsilon, a, estimand, specific.time, prob.bound)
    }
  } else {
    if (n_para_1>1) {
      out_bic_1 <- t(binarize_if_continuous(x_l[,2])[1,1:2])
      out_bic_2 <- t(binarize_if_continuous(x_l[,2])[1,4:5])
      if (n_para_1>2) {
        for (i_para in 3:n_para_1) {
          out_bic_i <- binarize_if_continuous(x_l[,i_para])
          out_bic_1 <- cbind(out_bic_1, t(out_bic_i[1,2]))
          out_bic_2 <- cbind(out_bic_2, t(out_bic_i[1,5]))
        }
        out_bic_1 <- cbind(out_bic_1, t(out_bic_i[1,3]))
        out_bic_2 <- cbind(out_bic_2, t(out_bic_i[1,6]))
      } else {
        out_bic_1 <- cbind(out_bic_1, t(binarize_if_continuous(x_l[,2])[1,3]))
        out_bic_2 <- cbind(out_bic_2, t(binarize_if_continuous(x_l[,2])[1,6]))
      }
      init_vals <- cbind(out_bic_1, out_bic_2)
    } else {
      init_vals <- calc_initial_2_competing_risk(t, epsilon, a, estimand, specific.time, prob.bound)
    }
  }
  out <- list(init_vals = init_vals,
              y_0=y_0, y_1=y_1, y_2=y_2,
              y_0_=y_0_, y_1_=y_1_, y_2_=y_2_,
              a=a, t=t, epsilon=epsilon, epsilon0=epsilon0, epsilon1=epsilon1, x_l=x_l, x_la=x_la, offset=offset,
              n_para_1=n_para_1, n_para_2=n_para_2, n_para_3=n_para_3, n_para_4=n_para_4, n_para_5=n_para_5, one=one)
  return(out)
}

check_input <- function(epsilon, outcome.type, optim.method) {
  if (!all(epsilon %in% c(0, 1, 2))) {
    stop("Invalid event code. Must be 0, 1 or 2, with 0 representing censoring")
  } else if (outcome.type == 'COMPETINGRISK' && all(epsilon %in% c(0, 1))) {
    stop("Invalid event code. Expected both event code 1 and 2 with competing risks outcome")
  } else if (outcome.type == 'SURVIVAL' && !all(epsilon %in% c(0, 1))) {
    stop("Invalid event code. Must be 0 or 1, with 0 representing censoring")
  }
  if (outcome.type == "COMPETINGRISK" && !optim.method$outer.optim.method %in% c("nleqslv","Newton","Broyden","optim","BFGS","SANN","multiroot","partial")) {
    stop("Invalid input for 'optimization'. Choose 'nleqslv', 'Newton', 'Broyden', 'optim', 'BFGS', 'SANN', 'multiroot' or 'partial'.")
  }
  if (outcome.type == "SURVIVAL" && optim.method$outer.optim.method == "partial") {
    stop("Invalid input for 'optimization'. Choose 'nleqslv', 'Newton', 'Broyden', 'optim', 'BFGS', 'SANN' or 'multiroot'.")
  }
  if (!optim.method$inner.optim.method %in% c("optim","BFGS","SANN","multiroot")) {
    stop("Invalid input for 'optimization'. Choose 'optim', 'BFGS', 'SANN' or 'multiroot'.")
  }
}

calc_initial_1_competing_risk <- function(t, epsilon, a, l, estimand, specific.time, prob.bound) {
  epsilon00 <- epsilon[a == 0 & l == 0]
  epsilon10 <- epsilon[a == 1 & l == 0]
  epsilon01 <- epsilon[a == 0 & l == 1]
  epsilon11 <- epsilon[a == 1 & l == 1]
  t00 <- t[a == 0 & l == 0]
  t10 <- t[a == 1 & l == 0]
  t01 <- t[a == 0 & l == 1]
  t11 <- t[a == 1 & l == 1]

  p_100 <- (sum(epsilon00==1 & t00<=specific.time)/length(epsilon00)) + prob.bound
  p_200 <- (sum(epsilon00==2 & t00<=specific.time)/length(epsilon00)) + prob.bound
  p_000 <- 1 - p_100 - p_200
  p_110 <- (sum(epsilon10==1 & t10<=specific.time)/length(epsilon10)) + prob.bound
  p_210 <- (sum(epsilon10==2 & t10<=specific.time)/length(epsilon10)) + prob.bound
  p_010 <- 1 - p_110 - p_210
  p_101 <- (sum(epsilon01==1 & t01<=specific.time)/length(epsilon01)) + prob.bound
  p_201 <- (sum(epsilon01==2 & t01<=specific.time)/length(epsilon01)) + prob.bound
  p_001 <- 1 - p_101 - p_201
  p_111 <- (sum(epsilon11==1 & t11<=specific.time)/length(epsilon11)) + prob.bound
  p_211 <- (sum(epsilon11==2 & t11<=specific.time)/length(epsilon11)) + prob.bound
  p_011 <- 1 - p_111 - p_211

  if (p_100 == 0 | is.na(p_100)) { stop("Compelete separation detected in initial value search") }
  if (p_200 == 0 | is.na(p_200)) { stop("Compelete separation detected in initial value search") }
  if (p_000 == 0 | is.na(p_000)) { stop("Compelete separation detected in initial value search") }
  if (p_110 == 0 | is.na(p_100)) { stop("Compelete separation detected in initial value search") }
  if (p_210 == 0 | is.na(p_200)) { stop("Compelete separation detected in initial value search") }
  if (p_010 == 0 | is.na(p_000)) { stop("Compelete separation detected in initial value search") }
  if (p_101 == 0 | is.na(p_101)) { stop("Compelete separation detected in initial value search") }
  if (p_201 == 0 | is.na(p_201)) { stop("Compelete separation detected in initial value search") }
  if (p_001 == 0 | is.na(p_001)) { stop("Compelete separation detected in initial value search") }
  if (p_111 == 0 | is.na(p_101)) { stop("Compelete separation detected in initial value search") }
  if (p_211 == 0 | is.na(p_201)) { stop("Compelete separation detected in initial value search") }
  if (p_011 == 0 | is.na(p_001)) { stop("Compelete separation detected in initial value search") }

  alpha_10 <- log( (p_100*p_110/(p_000*p_010) ))
  alpha_20 <- log( (p_200*p_210/(p_000*p_010) ))
  alpha_11 <- log( (p_101*p_111/(p_000*p_011) )) - alpha_10
  alpha_21 <- log( (p_201*p_211/(p_000*p_011) )) - alpha_20

  epsilon0 <- epsilon[a == 0]
  epsilon1 <- epsilon[a == 1]
  t0 <- t[a == 0]
  t1 <- t[a == 1]

  p_10 <- (sum(epsilon0==1 & t0<=specific.time)/length(epsilon0)) + prob.bound
  p_20 <- (sum(epsilon0==2 & t0<=specific.time)/length(epsilon0)) + prob.bound
  p_00 <- 1 - p_10 - p_20
  p_11 <- (sum(epsilon1==1 & t1<=specific.time)/length(epsilon1)) + prob.bound
  p_21 <- (sum(epsilon1==2 & t1<=specific.time)/length(epsilon1)) + prob.bound
  p_01 <- 1 - p_11 - p_21

  if (estimand$effect.measure1 == 'RR') {
    beta_1 <- log(p_11/p_10)
  } else if (estimand$effect.measure1 == 'OR') {
    beta_1 <- log( (p_11/(1-p_11)) / (p_10/(1-p_10)) )
  } else if (estimand$effect.measure1 == 'SHR') {
    beta_1 <- log( log(1 - p_11) / log(1 - p_10) )
  } else {
    stop("Invalid effect measure code. Must be RR, OR or SHR.")
  }
  if (estimand$effect.measure2 == 'RR') {
    beta_2 <- log(p_21/p_20)
  } else if (estimand$effect.measure2 == 'OR') {
    beta_2 <- log( (p_21/(1-p_21)) / (p_20/(1-p_20)) )
  } else if (estimand$effect.measure2 == 'SHR') {
    beta_2 <- log( log(1 - p_21) / log(1 - p_20) )
  } else {
    stop("Invalid effect measure code. Must be RR, OR or SHR.")
  }
  alpha_beta_univariable=cbind(alpha_10, alpha_11, beta_1, alpha_20, alpha_21, beta_2)
  return(alpha_beta_univariable)
}

calc_initial_2_competing_risk <- function(t, epsilon, a, estimand, specific.time, prob.bound) {
  epsilon0 <- epsilon[a == 0]
  epsilon1 <- epsilon[a == 1]
  t0 <- t[a == 0]
  t1 <- t[a == 1]

  p_10 <- (sum(epsilon0==1 & t0<=specific.time)/length(epsilon0)) + prob.bound
  p_20 <- (sum(epsilon0==2 & t0<=specific.time)/length(epsilon0)) + prob.bound
  p_00 <- 1 - p_10 - p_20
  p_11 <- (sum(epsilon1==1 & t1<=specific.time)/length(epsilon1)) + prob.bound
  p_21 <- (sum(epsilon1==2 & t1<=specific.time)/length(epsilon1)) + prob.bound
  p_01 <- 1 - p_11 - p_21

  if (estimand$effect.measure1 == 'RR') {
    beta_1 <- log(p_11/p_10)
  } else if (estimand$effect.measure1 == 'OR') {
    beta_1 <- log( (p_11/(1-p_11)) / (p_10/(1-p_10)) )
  } else if (estimand$effect.measure1 == 'SHR') {
    beta_1 <- log( log(1 - p_11) / log(1 - p_10) )
  } else {
    stop("Invalid effect measure code. Must be RR, OR or SHR.")
  }
  if (estimand$effect.measure2 == 'RR') {
    beta_2 <- log(p_21/p_20)
  } else if (estimand$effect.measure2 == 'OR') {
    beta_2 <- log( (p_21/(1-p_21)) / (p_20/(1-p_20)) )
  } else if (estimand$effect.measure2 == 'SHR') {
    beta_2 <- log( log(1 - p_21) / log(1 - p_20) )
  } else {
    stop("Invalid effect measure code. Must be RR, OR or SHR.")
  }
  alpha_1 <- log( (p_10*p_11/(p_00*p_01) ))
  alpha_2 <- log( (p_20*p_21/(p_00*p_01) ))
  alpha_beta_nocovariates <- cbind(alpha_1, beta_1, alpha_2, beta_2)
  return(alpha_beta_nocovariates)
}

calc_initial_1_survival <- function(t, epsilon, a, l, estimand, specific.time, prob.bound) {
  epsilon00 <- epsilon[a == 0 & l == 0]
  epsilon10 <- epsilon[a == 1 & l == 0]
  epsilon01 <- epsilon[a == 0 & l == 1]
  epsilon11 <- epsilon[a == 1 & l == 1]
  t00 <- t[a == 0 & l == 0]
  t10 <- t[a == 1 & l == 0]
  t01 <- t[a == 0 & l == 1]
  t11 <- t[a == 1 & l == 1]

  p_100 <- (sum(epsilon00==1 & t00<=specific.time)/length(epsilon00)) + prob.bound
  p_000 <- 1 - p_100
  p_110 <- (sum(epsilon10==1 & t10<=specific.time)/length(epsilon10)) + prob.bound
  p_010 <- 1 - p_110
  p_101 <- (sum(epsilon01==1 & t01<=specific.time)/length(epsilon01)) + prob.bound
  p_001 <- 1 - p_101
  p_111 <- (sum(epsilon11==1 & t11<=specific.time)/length(epsilon11)) + prob.bound
  p_011 <- 1 - p_111

  if (p_100 == 0 | is.na(p_100)) { stop("Compelete separation detected in initial value search") }
  if (p_000 == 0 | is.na(p_000)) { stop("Compelete separation detected in initial value search") }
  if (p_110 == 0 | is.na(p_100)) { stop("Compelete separation detected in initial value search") }
  if (p_010 == 0 | is.na(p_000)) { stop("Compelete separation detected in initial value search") }
  if (p_101 == 0 | is.na(p_101)) { stop("Compelete separation detected in initial value search") }
  if (p_001 == 0 | is.na(p_001)) { stop("Compelete separation detected in initial value search") }
  if (p_111 == 0 | is.na(p_101)) { stop("Compelete separation detected in initial value search") }
  if (p_011 == 0 | is.na(p_001)) { stop("Compelete separation detected in initial value search") }

  alpha_10 <- log( (p_100*p_110/(p_000*p_010) ))
  alpha_11 <- log( (p_101*p_111/(p_000*p_011) )) - alpha_10

  epsilon0 <- epsilon[a == 0]
  epsilon1 <- epsilon[a == 1]
  t0 <- t[a == 0]
  t1 <- t[a == 1]
  p_10 <- (sum(epsilon0==1 & t0<=specific.time)/length(epsilon0)) + prob.bound
  p_00 <- 1 - p_10
  p_11 <- (sum(epsilon1==1 & t1<=specific.time)/length(epsilon1)) + prob.bound
  p_01 <- 1 - p_11

  if (estimand$effect.measure1 == 'RR') {
    beta_1 <- log(p_11/p_10)
  } else if (estimand$effect.measure1 == 'OR') {
    beta_1 <- log( (p_11/(1-p_11)) / (p_10/(1-p_10)) )
  } else if (estimand$effect.measure1 == 'SHR') {
    beta_1 <- log( log(1 - p_11) / log(1 - p_10) )
  } else {
    stop("Invalid effect measure code. Must be RR, OR or SHR.")
  }
  alpha_beta_univariable=cbind(alpha_10, alpha_11, beta_1)
  return(alpha_beta_univariable)
}

calc_initial_2_survival <- function(t, epsilon, a, estimand, specific.time, prob.bound) {
  epsilon0 <- epsilon[a == 0]
  epsilon1 <- epsilon[a == 1]
  t0 <- t[a == 0]
  t1 <- t[a == 1]
  p_10 <- (sum(epsilon0==1 & t0<=specific.time)/length(epsilon0)) + prob.bound
  p_00 <- 1 - p_10
  p_11 <- (sum(epsilon1==1 & t1<=specific.time)/length(epsilon1)) + prob.bound
  p_01 <- 1 - p_11

  if (estimand$effect.measure1 == 'RR') {
    beta_1 <- log(p_11/p_10)
  } else if (estimand$effect.measure1 == 'OR') {
    beta_1 <- log( (p_11/(1-p_11)) / (p_10/(1-p_10)) )
  } else if (estimand$effect.measure1 == 'SHR') {
    beta_1 <- log( log(1 - p_11) / log(1 - p_10) )
  } else {
    stop("Invalid effect measure code. Must be RR, OR or SHR.")
  }
  alpha_1 <- log( (p_10*p_11/(p_00*p_01) ))
  alpha_beta_nocovariates <- cbind(alpha_1, beta_1)
  return(alpha_beta_nocovariates)
}

normalize_covariate <- function(formula, data, covariate.normalization, outcome.type) {
  mf <- model.frame(formula, data)
  Y <- model.extract(mf, "response")
  response_term <- formula[[2]]
  if (inherits(mf[[1]], "Surv") || inherits(mf[[1]], "Event")) {
    # If it's a Surv object, extract the covariates excluding time and event
    response_vars <- all.vars(response_term)
    covariate_cols <- setdiff(all.vars(formula), response_vars)  # Remove time and event
  } else {
    # If it's a regular response, exclude just the response variable
    covariate_cols <- all.vars(formula)[-1]  # Exclude the response variable
  }

  normalized_data <- data
  range_vector <- 1
  if (covariate.normalization == TRUE && length(covariate_cols)>0) {
    normalize <- function(x) {
      return ((x - min(x)) / (max(x) - min(x)))
    }
    for (col in covariate_cols) {
      x <- normalized_data[[col]]
      range <- max(x)-min(x)
      normalized_data[[col]] <- (x - min(x)) / range
      range_vector <- cbind(range_vector,range)
    }
    if (outcome.type == 'PROPORTIONAL') {
      range_vector <- cbind(range_vector)
    } else if (outcome.type == 'SURVIVAL') {
      range_vector <- cbind(range_vector,1)
    } else {
      range_vector <- cbind(range_vector,1,range_vector,1)
    }
  } else {
    if (outcome.type == 'PROPORTIONAL') {
      range_vector <- rep(1, (length(covariate_cols)+1))
    } else if (outcome.type == 'SURVIVAL') {
      range_vector <- rep(1, (length(covariate_cols)+2))
    } else {
      range_vector <- rep(1, (2*length(covariate_cols)+4))
    }
  }
  n_covariate <- length(covariate_cols)
  out <- list(normalized_data=normalized_data, range=range_vector, n_covariate=n_covariate, n=nrow(data))
  return(out)
}

sort_by_covariate <- function(formula, data, sort.data, n_covariate) {
  if (sort.data == TRUE && n_covariate>0) {
    terms_obj <- terms(formula)
    covariate_names <- attr(terms_obj, "term.labels")
    missing_vars <- setdiff(covariate_names, names(data))
    if (length(missing_vars) > 0) {
      stop("The following covariates are missing: ", paste(missing_vars, collapse = ", "))
    }
    sorted_data <- data[do.call(order, data[covariate_names]), , drop = FALSE]
    return(sorted_data)
  } else {
    return(data)
  }
}

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

  pred <- calc_pred(alpha_beta,x_l,offset,epsilon0,epsilon1,one,n_para_1,estimand,optim.method,prob.bound,initial_pred)
  # pred: n×4 => columns: p_10, p_20, p_11, p_21

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
  # pred: n×4 => columns: p_10, p_20, p_11, p_21

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

  y_0_ <- ifelse(epsilon == 0, 1, 0) # 打ち切り変数
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
    # pred: n×4 => columns: p_10, p_20, p_11, p_21

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

calc_pred <- function(
    alpha_beta_tmp,
    x_l, offset,
    epsilon0, epsilon1,
    one,
    n_para_1,
    estimand,
    optim.method,
    prob.bound,
    initial_pred = NULL
) {
  idx_alpha1 <- 1:n_para_1
  idx_beta1  <- n_para_1 + 1
  idx_alpha2 <- (n_para_1 + 2):(2 * n_para_1 + 1)
  idx_beta2  <- 2 * n_para_1 + 2

  alpha_1 <- alpha_beta_tmp[idx_alpha1]
  beta_1  <- alpha_beta_tmp[idx_beta1]
  alpha_2 <- alpha_beta_tmp[idx_alpha2]
  beta_2  <- alpha_beta_tmp[idx_beta2]

  alpha_tmp_1_vals <- x_l[, 1:n_para_1] %*% as.matrix(alpha_1) + offset
  alpha_tmp_2_vals <- x_l[, 1:n_para_1] %*% as.matrix(alpha_2) + offset
  beta_tmp_1_val  <- beta_1
  beta_tmp_2_val  <- beta_2

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
      log_p <- log(initial_pred[i_x, ])
    } else {
      n_epsilon0 <- length(epsilon0)
      n_epsilon1 <- length(epsilon1)
      p0 <- c(
        sum(epsilon0 == 1) / n_epsilon0 + prob.bound,
        sum(epsilon0 == 2) / n_epsilon0 + prob.bound,
        sum(epsilon1 == 1) / n_epsilon1 + prob.bound,
        sum(epsilon1 == 2) / n_epsilon1 + prob.bound
      )
      log_p <- log(p0)
    }

    if (optim.method$inner.optim.method == 'optim' || optim.method$inner.optim.method == 'BFGS') {
      eq_fn <- function(lp) {
        estimating_equation_pred(
          log_p          = lp,
          alpha_tmp_1    = alpha_tmp_1_vals[i_x],
          beta_tmp_1     = beta_tmp_1_val,
          alpha_tmp_2    = alpha_tmp_2_vals[i_x],
          beta_tmp_2     = beta_tmp_2_val,
          estimand = estimand,
          optim.method = optim.method,
          prob.bound = prob.bound
        )
      }
      sol <- optim(par = log_p, fn = eq_fn, method = "BFGS", control = list(maxit = optim.method$optim.parameter8, reltol = optim.method$optim.parameter7))
      pred_list[[i_x]] <- exp(sol$par)
      prev_pred <- pred_list[[i_x]]
    } else if (optim.method$inner.optim.method == 'SANN') {
      eq_fn <- function(lp) {
        estimating_equation_pred(
          log_p          = lp,
          alpha_tmp_1    = alpha_tmp_1_vals[i_x],
          beta_tmp_1     = beta_tmp_1_val,
          alpha_tmp_2    = alpha_tmp_2_vals[i_x],
          beta_tmp_2     = beta_tmp_2_val,
          estimand = estimand,
          optim.method = optim.method,
          prob.bound = prob.bound
        )
      }
      sol <- optim(par = log_p, fn = eq_fn, method = "SANN", control = list(maxit = optim.method$optim.parameter8, reltol = optim.method$optim.parameter7))
      pred_list[[i_x]] <- exp(sol$par)
      prev_pred <- pred_list[[i_x]]
    } else if (optim.method$inner.optim.method == 'multiroot'){
      eq_fn <- function(lp) {
        estimating_equation_pred(
          log_p          = lp,
          alpha_tmp_1    = alpha_tmp_1_vals[i_x],
          beta_tmp_1     = beta_tmp_1_val,
          alpha_tmp_2    = alpha_tmp_2_vals[i_x],
          beta_tmp_2     = beta_tmp_2_val,
          estimand = estimand,
          optim.method = optim.method,
          prob.bound = prob.bound
        )
      }
      sol <- multiroot(eq_fn, start = log_p)
      pred_list[[i_x]] <- exp(sol$root)
      prev_pred <- pred_list[[i_x]]
    }
  }
  pred <- do.call(rbind, pred_list)
  return(pred)
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

    if ((1 - exp(log_p[1]) - exp(log_p[2]) < prob.bound) ||
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
