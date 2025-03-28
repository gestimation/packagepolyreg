#' Title
#'
#' @param formula formula Model formula representing outcome and strata
#' @param data data.frame Input dataset containing survival data.
#' @param weights character Column name representing the weights. The weights must be nonnegative and it is strongly recommended that they be strictly positive, since zero weights are ambiguous, compared to use of the subset argument.
#' @param subset character Specifies a condition for subsetting the data. Defaults to NULL.
#' @param code.event integer Specifies the code of event. Defaults to 1.
#' @param code.censoring integer Specifies the code of censoring. Defaults to 0.
#' @param na.action character Specifies a missing-data filter function, applied to the model frame, after any subset argument has been used. Defaults to na.pass.
#' @param conf.int numeric The level for a two-sided confidence interval on the survival probabilities. Defaults to 0.95.
#' @param error character either the string "greenwood" for the Greenwood formula or "tsiatis" for the Tsiatis formula. Defaults to "greenwood".
#' @param conf.type character Specifies transformation used to construct the confidence interval on the probabilities. Defaults to "arcsine-square root".
#'
#' @returns
#' @export km.curve
#'
#' @examples
km.curve <- function(formula,
                     data,
                     weights = NULL,
                     subset = NULL,
                     code.event = 1,
                     code.censoring = 0,
                     na.action = na.pass,
                     conf.int = 0.95,
                     error = "greenwood",
                     conf.type = "arcsine-square root"
) {
  checkDependentPackages()
  out_readSurv <- readSurv(formula, data, weights, code.event, code.censoring, subset, na.action)
  out_calculateKaplanMeier_rcpp <- calculateKaplanMeier_rcpp(out_readSurv$t, out_readSurv$d, out_readSurv$w, out_readSurv$strata, error)
  ci <- calculateConfidenceInterval(out_calculateKaplanMeier_rcpp, conf.int, conf.type, conf.lower)
  if (all(out_readSurv$strata == 1)) {
    survfit_object <- list(
      time = out_calculateKaplanMeier_rcpp$time,
      surv = out_calculateKaplanMeier_rcpp$surv,
      n = out_calculateKaplanMeier_rcpp$n,
      n.risk = out_calculateKaplanMeier_rcpp$n.risk,
      n.event = out_calculateKaplanMeier_rcpp$n.event,
      n.censor = out_calculateKaplanMeier_rcpp$n.censor,
      std.err = out_calculateKaplanMeier_rcpp$std.err,
      upper = ci$upper,
      lower = ci$lower,
      conf.type = conf.type,
      call = match.call(),
      type = "kaplan-meier",
      method = "Kaplan-Meier"
    )
  } else {
    survfit_object <- list(
      time = out_calculateKaplanMeier_rcpp$time,
      surv = out_calculateKaplanMeier_rcpp$surv,
      n = out_calculateKaplanMeier_rcpp$n,
      n.risk = out_calculateKaplanMeier_rcpp$n.risk,
      n.event = out_calculateKaplanMeier_rcpp$n.event,
      n.censor = out_calculateKaplanMeier_rcpp$n.censor,
      std.err = out_calculateKaplanMeier_rcpp$std.err,
      upper = ci$upper,
      lower = ci$lower,
      conf.type = conf.type,
      strata = out_calculateKaplanMeier_rcpp$strata,
      call = match.call(),
      type = "kaplan-meier",
      method = "Kaplan-Meier"
    )
  }
  class(survfit_object) <- c("survfit")
  return(survfit_object)
}

calculateConfidenceInterval <- function(survfit_object, conf.int, conf.type, conf.lower) {
  if (conf.int <= 0 | conf.int >= 1)
    stop("Confidence level must be between 0 and 1")
  alpha <- 1 - conf.int
  critical_value <- qnorm(1 - alpha / 2)
  if (is.null(conf.type) | conf.type == "none") {
    lower <- NULL
    upper <- NULL
  } else if (conf.type == "arcsine-square root" | conf.type == "arcsin" | conf.type == "a") {
    se <- survfit_object$surv*survfit_object$std.err/2/sqrt(survfit_object$surv * (1 - survfit_object$surv))
    lower <- sin(pmax(asin(sqrt(survfit_object$surv)) - critical_value*se, 0))^2
    upper <- sin(pmin(asin(sqrt(survfit_object$surv)) + critical_value*se, pi/2))^2
  } else if (conf.type == "plain" | conf.type == "p" | conf.type == "linear") {
    lower <- pmax(survfit_object$surv - critical_value*survfit_object$surv*survfit_object$std.err, 0)
    upper <- pmin(survfit_object$surv + critical_value*survfit_object$surv*survfit_object$std.err, 1)
  } else if (conf.type == "log") {
    se <- survfit_object$std.err
    lower <- survfit_object$surv * exp(-critical_value*se)
    upper <- pmin(survfit_object$surv * exp(critical_value*se), 1)
  } else if (conf.type == "log-log") {
    se <- survfit_object$std.err / log(survfit_object$surv)
    lower <- survfit_object$surv^exp(-critical_value*se)
    upper <- survfit_object$surv^exp(critical_value*se)
  } else if (conf.type == "logit") {
    se <- survfit_object$std.err/(1 - survfit_object$surv)
    lower <- survfit_object$surv / (survfit_object$surv + (1 - survfit_object$surv)*exp(critical_value*se))
    upper <- survfit_object$surv / (survfit_object$surv + (1 - survfit_object$surv)*exp(-critical_value*se))
  }
  lower <- sapply(lower, function(x) ifelse(is.nan(x), NA, x))
  upper <- sapply(upper, function(x) ifelse(is.nan(x), NA, x))
  lower <- sapply(lower, function(x) ifelse(x>=1, 1, x))
  upper <- sapply(upper, function(x) ifelse(x>=1, 1, x))
  lower <- sapply(lower, function(x) ifelse(x<=0, 0, x))
  upper <- sapply(upper, function(x) ifelse(x<=0, 0, x))
  return(list(upper=upper, lower=lower))
}

createAnalysisDataset <- function(formula, data, weights, subset_condition, na.action) {
  if (!is.null(subset_condition)) {
    analysis_dataset <- subset(data, eval(parse(text = subset_condition)))
  } else {
    analysis_dataset <- data
  }
  all_vars <- all.vars(formula)
  all_vars <- c(all_vars, weights)
  analysis_dataset <- analysis_dataset[, all_vars, drop = FALSE]
  return(na.action(analysis_dataset))
}

checkDependentPackages <- function() {
  if (requireNamespace("mets", quietly = TRUE) & requireNamespace("Rcpp", quietly = TRUE)) {
    suppressWarnings(library(mets))
    suppressWarnings(library(Rcpp))
  } else {
    stop("Required packages 'mets' and/or 'Rcpp' are not installed.")
  }
}

readSurv <- function(formula, data, weights, code.event, code.censoring, subset_condition, na.action) {
  data <- createAnalysisDataset(formula, data, weights, subset_condition, na.action)
  cl <- match.call()
  if (missing(formula))
    stop("A formula argument is required")
  mf <- match.call(expand.dots = TRUE)[1:3]
  special <- c("strata", "offset", "cluster")
  out_terms <- terms(formula, special, data = data)
  if (!is.null(attr(out_terms, "specials")$strata))
    stop("strata() cannot appear in formula")
  if (!is.null(attr(out_terms, "specials")$offset))
    stop("offset() cannot appear in formula")
  if (!is.null(attr(out_terms, "specials")$cluster))
    stop("cluster() cannot appear in formula")
  mf$formula <- out_terms
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  Y <- model.extract(mf, "response")
  if (!inherits(Y, c("Event", "Surv"))) {
    stop("Expected a 'Surv' or 'Event'-object when outcome.type is COMPETINGRISK, SURVIVAL or PROPORTIONAL. ")
  } else {
    t <- as.numeric(Y[, 1])
    if (any(t<0)) {
      stop("Invalid time variable. Expected non-negative values. ")
    }
    if (any(is.na(t)))
      stop("Time variable contains NA values")
    if (any(is.na(Y[, 2]))) {
      stop("Event variable contains NA values")
    } else if (!all(Y[, 2] %in% c(code.event, code.censoring))) {
      stop("Invalid event codes. Must be 0 or 1, with 0 representing censoring, if event codes are not specified. ")
    } else {
      d <- ifelse(Y[, 2] == code.censoring, 0, 1)
    }
  }
  if (is.na(all.vars(out_terms)[3])) {
    strata <- rep(1, nrow(data))
  } else {
    strata_name <- all.vars(out_terms)[3]
    if (any(is.na(data[[strata_name]]))) {
      stop("Strata variable contains NA values")
    } else{
      strata <- as.integer(as.factor(data[[strata_name]]))
    }
  }
  if (is.null(weights)) {
    w <- rep(1, nrow(data))
  } else {
    w <- data[[weights]]
    if (!is.numeric(w))
      stop("weights must be numeric")
    if (any(!is.finite(w)))
      stop("weights must be finite")
    if (any(w < 0))
      stop("weights must be non-negative")
    if (any(is.na(w)))
      stop("weights contain NA values")
  }
  return(list(t = t, d = d, strata = strata, w=w))
}

##############################################################################################################
library(survival)
createTestData <- function(n, w, first_zero=FALSE, last_zero=FALSE, subset_present=FALSE, logical_strata=FALSE, na_strata=FALSE) {
  one <- rep(1, n)
  t <- c(1:(n/2), 1:(n/2))
  epsilon <- rep(1, n)
  epsilon[2] <- 2
  epsilon[3] <- 2
  if (first_zero==TRUE) {
    epsilon[1] <- 0
    epsilon[n/2+1] <- 0
  }
  if (last_zero==TRUE) {
    epsilon[n/2] <- 0
    epsilon[n] <- 0
  }
  w <- rep(w, n)
  if (logical_strata==TRUE) {
    strata <- (t %% 2 == 0)
  } else {
    strata <- as.factor((t %% 2 == 0))
  }
  if (na_strata==TRUE) {
    strata[1] <- NA
  }
  subset <- rep(1, n)
  if (subset_present==TRUE) {
    subset[1] <- 0
  }
  d <- as.numeric(epsilon>0)
  return(data.frame(id = 1:n, t = t, epsilon = epsilon, d = d, w = w, strata = strata, subset=subset))
}


##############################################################################################################
library(Rcpp)

cppFunction('
Rcpp::List calculateKaplanMeier_rcpp(Rcpp::NumericVector t, Rcpp::IntegerVector d,
                                    Rcpp::NumericVector w = Rcpp::NumericVector::create(),
                                    Rcpp::IntegerVector strata = Rcpp::IntegerVector::create(),
                                    Rcpp::CharacterVector error = Rcpp::CharacterVector::create("greenwood")) {

//  Rcpp::Rcout << "Method: " << error[0] << std::endl;

  Rcpp::List km_list;

  std::vector<double> combined_times;
  std::vector<double> combined_surv;
  std::vector<int> combined_n_risk;
  std::vector<int> combined_n_event;
  std::vector<int> combined_n_censor;
  std::vector<double> combined_std_err;
  std::vector<int> combined_n_stratum;
  std::vector<int> combined_u_stratum;

  if ((strata.size() == 0 || Rcpp::unique(strata).size() == 1) && (w.size() == 0 || (Rcpp::unique(w).size() == 1 && w[0] == 1))) {

    Rcpp::NumericVector unique_times = Rcpp::unique(t);
    std::sort(unique_times.begin(), unique_times.end());

    int u = unique_times.size();
    Rcpp::NumericVector km(u);
    Rcpp::NumericVector km_i(u);
    Rcpp::IntegerVector weighted_n_risk(u);
    Rcpp::IntegerVector weighted_n_event(u);
    Rcpp::IntegerVector weighted_n_censor(u);
    Rcpp::NumericVector std_err(u);

    int n_stratum = t.size();
    int n = t.size();
    for (int i = 0; i < u; ++i) {
      double time = unique_times[i];
      double weighted_n_i = 0;
      for (int j = 0; j < n; ++j) {
        if (t[j] >= time) {
          weighted_n_i++;
        }
      }
      weighted_n_risk[i] = weighted_n_i;

      double weighted_events = 0;
      for (int j = 0; j < n; ++j) {
        if (t[j] == time && d[j] == 1) {
          weighted_events++;
        }
      }

      if (weighted_n_risk[i] > 0) {
        km_i[i] = 1 - weighted_events / weighted_n_risk[i];
        weighted_n_event[i] += weighted_events;
      } else {
        km_i[i] = 1;
      }

      if (i > 0) {
        weighted_n_censor[i-1] = weighted_n_risk[i-1] - weighted_n_risk[i] - weighted_n_event[i-1];
      }
      if (i == u-1) {
        weighted_n_censor[i] = weighted_n_risk[i] - weighted_n_event[i];
      }

      if (i == 0) {
        km[i] = km_i[i];
      } else {
        km[i] = km[i - 1] * km_i[i];
      }

      double sum_se = 0;
      for (int j = 0; j <= i; ++j) {
        double n_i = weighted_n_risk[j];
        double d_i = 0;
        for (int k = 0; k < n; ++k) {
          if (t[k] == unique_times[j] && d[k] == 1) {
            d_i ++;
          }
        }
        if (n_i > d_i) {
          if (error[0] == "tsiatis") {
            sum_se += (d_i / (n_i * n_i));
          } else {
            sum_se += (d_i / (n_i * (n_i - d_i)));
          }
        } else {
          sum_se = std::numeric_limits<double>::infinity();
        }
      }
      std_err[i] = sqrt(sum_se);
    }

    combined_times.insert(combined_times.end(), unique_times.begin(), unique_times.end());
    combined_surv.insert(combined_surv.end(), km.begin(), km.end());
    combined_n_risk.insert(combined_n_risk.end(), weighted_n_risk.begin(), weighted_n_risk.end());
    combined_n_event.insert(combined_n_event.end(), weighted_n_event.begin(), weighted_n_event.end());
    combined_n_censor.insert(combined_n_censor.end(), weighted_n_censor.begin(), weighted_n_censor.end());
    combined_std_err.insert(combined_std_err.end(), std_err.begin(), std_err.end());
    combined_n_stratum.insert(combined_n_stratum.end(), n_stratum);

  } else if (strata.size() == 0 || Rcpp::unique(strata).size() == 1) {

    Rcpp::NumericVector unique_times = Rcpp::unique(t);
    std::sort(unique_times.begin(), unique_times.end());

    int u = unique_times.size();
    Rcpp::NumericVector km(u);
    Rcpp::NumericVector km_i(u);
    Rcpp::IntegerVector weighted_n_risk(u);
    Rcpp::IntegerVector weighted_n_event(u);
    Rcpp::IntegerVector weighted_n_censor(u);
    Rcpp::NumericVector std_err(u);

    int n_stratum = t.size();
    int n = t.size();
    for (int i = 0; i < u; ++i) {
      double time = unique_times[i];
      double weighted_n_i = 0;
      double unweighted_n_i = 0;
      for (int j = 0; j < n; ++j) {
        if (t[j] >= time) {
          weighted_n_i += w[j];
        }
      }
      weighted_n_risk[i] = weighted_n_i;

      double weighted_events = 0;
      for (int j = 0; j < n; ++j) {
        if (t[j] == time && d[j] == 1) {
          weighted_events += w[j];
        }
      }

      if (weighted_n_risk[i] > 0) {
        km_i[i] = 1 - weighted_events / weighted_n_risk[i];
        weighted_n_event[i] += weighted_events;
      } else {
        km_i[i] = 1;
      }

      if (i > 0) {
        weighted_n_censor[i-1] = weighted_n_risk[i-1] - weighted_n_risk[i] - weighted_n_event[i-1];
      }
      if (i == u-1) {
        weighted_n_censor[i] = weighted_n_risk[i] - weighted_n_event[i];
      }

      if (i == 0) {
        km[i] = km_i[i];
      } else {
        km[i] = km[i - 1] * km_i[i];
      }

      double sum_se = 0;
      for (int j = 0; j <= i; ++j) {
        double n_i = weighted_n_risk[j];
        double d_i = 0;
        for (int k = 0; k < n; ++k) {
          if (t[k] == unique_times[j] && d[k] == 1) {
            d_i += w[k];
          }
        }
        if (n_i > d_i) {
          if (error[0] == "tsiatis") {
            sum_se += (d_i / (n_i * n_i));
          } else {
            sum_se += (d_i / (n_i * (n_i - d_i)));
          }
        } else {
          sum_se = std::numeric_limits<double>::infinity();
        }
      }
      std_err[i] = sqrt(sum_se);
    }

    combined_times.insert(combined_times.end(), unique_times.begin(), unique_times.end());
    combined_surv.insert(combined_surv.end(), km.begin(), km.end());
    combined_n_risk.insert(combined_n_risk.end(), weighted_n_risk.begin(), weighted_n_risk.end());
    combined_n_event.insert(combined_n_event.end(), weighted_n_event.begin(), weighted_n_event.end());
    combined_n_censor.insert(combined_n_censor.end(), weighted_n_censor.begin(), weighted_n_censor.end());
    combined_std_err.insert(combined_std_err.end(), std_err.begin(), std_err.end());
    combined_n_stratum.insert(combined_n_stratum.end(), n_stratum);

  } else if (w.size() == 0 || (Rcpp::unique(w).size() == 1 && w[0] == 1)) {

    Rcpp::IntegerVector strata_vec = strata;

    // Loop over each strata level and calculate the Kaplan-Meier estimate
    for (int i = 0; i < Rcpp::max(strata_vec); ++i) {
      // Select the data subset corresponding to the current strata level
      Rcpp::LogicalVector strata_condition = (strata_vec == (i + 1));
      Rcpp::NumericVector t_selected = t[strata_condition];
      Rcpp::IntegerVector d_selected = d[strata_condition];

      // Calculate weighted_n for this strata (sum of weights in the selected dataset)
      int n_stratum = t_selected.size();
      double weighted_n_stratum = 0;
      int n_event_stratum = 0;
      for (int j = 0; j < n_stratum; ++j) {
        weighted_n_stratum += d_selected[j];
        n_event_stratum += d_selected[j];
      }

      // Calculate Kaplan-Meier for the selected subset
      Rcpp::NumericVector unique_times = Rcpp::unique(t_selected);
      std::sort(unique_times.begin(), unique_times.end());
      int u_stratum = unique_times.size();

      Rcpp::NumericVector km(u_stratum);
      Rcpp::NumericVector km_i(u_stratum);
      Rcpp::IntegerVector weighted_n_risk(u_stratum);
      Rcpp::IntegerVector weighted_n_event(u_stratum);
      Rcpp::IntegerVector weighted_n_censor(u_stratum);
      Rcpp::NumericVector std_err(u_stratum);

      for (int j = 0; j < u_stratum; ++j) {
        double time = unique_times[j];
        double weighted_n_sub = 0;
        for (int k = 0; k < t_selected.size(); ++k) {
          if (t_selected[k] >= time) {
            weighted_n_sub++;
          }
        }
        weighted_n_risk[j] = weighted_n_sub;

        double weighted_events = 0;
        for (int k = 0; k < t_selected.size(); ++k) {
          if (t_selected[k] == time && d_selected[k] == 1) {
            weighted_events++;
          }
        }

        if (weighted_n_risk[j] > 0) {
          km_i[j] = 1 - weighted_events / weighted_n_risk[j];
          weighted_n_event[j] += weighted_events;
        } else {
          km_i[j] = 1;
        }

        if (j > 0) {
          weighted_n_censor[j-1] = weighted_n_risk[j-1] - weighted_n_risk[j] - weighted_n_event[j-1];
        }
        if (j == u_stratum-1) {
          weighted_n_censor[j] = weighted_n_risk[j] - weighted_n_event[j];
        }

        km[j] = (j == 0) ? km_i[j] : km[j - 1] * km_i[j];

        double sum_se = 0;
        for (int k = 0; k <= j; ++k) {
          double n_i = weighted_n_risk[k];
          double d_i = 0;
          for (int m = 0; m < t_selected.size(); ++m) {
            if (t_selected[m] == unique_times[k] && d_selected[m] == 1) {
              d_i += d_selected[m];
            }
          }
          if (n_i > d_i) {
            if (error[0] == "tsiatis") {
              sum_se += (d_i / (n_i * n_i));
            } else {
              sum_se += (d_i / (n_i * (n_i - d_i)));
            }
          } else {
            sum_se = std::numeric_limits<double>::infinity();
          }
        }
        std_err[j] = sqrt(sum_se);
      }

      combined_times.insert(combined_times.end(), unique_times.begin(), unique_times.end());
      combined_surv.insert(combined_surv.end(), km.begin(), km.end());
      combined_n_risk.insert(combined_n_risk.end(), weighted_n_risk.begin(), weighted_n_risk.end());
      combined_n_event.insert(combined_n_event.end(), weighted_n_event.begin(), weighted_n_event.end());
      combined_n_censor.insert(combined_n_censor.end(), weighted_n_censor.begin(), weighted_n_censor.end());
      combined_std_err.insert(combined_std_err.end(), std_err.begin(), std_err.end());
      combined_n_stratum.insert(combined_n_stratum.end(), n_stratum);
      combined_u_stratum.insert(combined_u_stratum.end(), u_stratum);
    }
  } else {

    Rcpp::IntegerVector strata_vec = strata;

    // Loop over each strata level and calculate the Kaplan-Meier estimate
    for (int i = 0; i < Rcpp::max(strata_vec); ++i) {
      // Select the data subset corresponding to the current strata level
      Rcpp::LogicalVector strata_condition = (strata_vec == (i + 1));  // 1-based indexing for strata
      Rcpp::NumericVector t_selected = t[strata_condition];
      Rcpp::IntegerVector d_selected = d[strata_condition];
      Rcpp::NumericVector w_selected = w[strata_condition];

      // Calculate weighted_n for this strata (sum of weights in the selected dataset)
      int n_stratum = t_selected.size();
      double weighted_n_stratum = 0;
      int n_event_stratum = 0;
      for (int j = 0; j < n_stratum; ++j) {
        weighted_n_stratum += w_selected[j];
        n_event_stratum += d_selected[j];
      }

      // Calculate Kaplan-Meier for the selected subset
      Rcpp::NumericVector unique_times = Rcpp::unique(t_selected);
      std::sort(unique_times.begin(), unique_times.end());
      int u_stratum = unique_times.size();

      Rcpp::NumericVector km(u_stratum);
      Rcpp::NumericVector km_i(u_stratum);
      Rcpp::IntegerVector weighted_n_risk(u_stratum);
      Rcpp::IntegerVector weighted_n_event(u_stratum);
      Rcpp::IntegerVector weighted_n_censor(u_stratum);
      Rcpp::NumericVector std_err(u_stratum);

      for (int j = 0; j < u_stratum; ++j) {
        double time = unique_times[j];

        double weighted_n_sub = 0;
        for (int k = 0; k < t_selected.size(); ++k) {
          if (t_selected[k] >= time) {
            weighted_n_sub += w_selected[k];
          }
        }
        weighted_n_risk[j] = weighted_n_sub;

        double weighted_events = 0;
        for (int k = 0; k < t_selected.size(); ++k) {
          if (t_selected[k] == time && d_selected[k] == 1) {
            weighted_events += w_selected[k];
          }
        }

        if (weighted_n_risk[j] > 0) {
          km_i[j] = 1 - weighted_events / weighted_n_risk[j];
          weighted_n_event[j] += weighted_events;
        } else {
          km_i[j] = 1;
        }

        if (j > 0) {
          weighted_n_censor[j-1] = weighted_n_risk[j-1]-weighted_n_risk[j] - weighted_n_event[j-1];
        }
        if (j == u_stratum-1) {
          weighted_n_censor[j] = weighted_n_risk[j] - weighted_n_event[j];
        }

        km[j] = (j == 0) ? km_i[j] : km[j - 1] * km_i[j];

        double sum_se = 0;
        for (int k = 0; k <= j; ++k) {
          double n_i = weighted_n_risk[k];
          double d_i = 0;
          for (int m = 0; m < t_selected.size(); ++m) {
            if (t_selected[m] == unique_times[k] && d_selected[m] == 1) {
              d_i += w_selected[m];
            }
          }
          if (n_i > d_i) {
            if (error[0] == "tsiatis") {
              sum_se += (d_i / (n_i * n_i));
            } else {
              sum_se += (d_i / (n_i * (n_i - d_i)));
            }
          } else {
            sum_se = std::numeric_limits<double>::infinity();
          }
        }
        std_err[j] = sqrt(sum_se);
      }

      combined_times.insert(combined_times.end(), unique_times.begin(), unique_times.end());
      combined_surv.insert(combined_surv.end(), km.begin(), km.end());
      combined_n_risk.insert(combined_n_risk.end(), weighted_n_risk.begin(), weighted_n_risk.end());
      combined_n_event.insert(combined_n_event.end(), weighted_n_event.begin(), weighted_n_event.end());
      combined_n_censor.insert(combined_n_censor.end(), weighted_n_censor.begin(), weighted_n_censor.end());
      combined_std_err.insert(combined_std_err.end(), std_err.begin(), std_err.end());
      combined_n_stratum.insert(combined_n_stratum.end(), n_stratum);
      combined_u_stratum.insert(combined_u_stratum.end(), u_stratum);
    }
  }

  Rcpp::NumericVector all_times = Rcpp::wrap(combined_times);
  Rcpp::NumericVector all_surv = Rcpp::wrap(combined_surv);
  Rcpp::IntegerVector all_n_risk = Rcpp::wrap(combined_n_risk);
  Rcpp::IntegerVector all_n_event = Rcpp::wrap(combined_n_event);
  Rcpp::IntegerVector all_n_censor = Rcpp::wrap(combined_n_censor);
  Rcpp::NumericVector all_std_err = Rcpp::wrap(combined_std_err);
  Rcpp::IntegerVector all_n_stratum = Rcpp::wrap(combined_n_stratum);
  Rcpp::IntegerVector all_u_stratum = Rcpp::wrap(combined_u_stratum);

  // Create final combined output for all strata
  km_list = Rcpp::List::create(
    Rcpp::_["time"] = all_times,
    Rcpp::_["surv"] = all_surv,
    Rcpp::_["n.risk"] = all_n_risk,
    Rcpp::_["n"] = combined_n_stratum,
    Rcpp::_["n.event"] = all_n_event,
    Rcpp::_["n.censor"] = all_n_censor,
    Rcpp::_["unweighted.n"] = all_n_stratum,
    Rcpp::_["std.err"] = all_std_err,
    Rcpp::_["high"] = R_NilValue,
    Rcpp::_["low"] = R_NilValue,
    Rcpp::_["conf.type"] = "log-log",
    Rcpp::_["strata"] = all_u_stratum,
    Rcpp::_["type"] = "kaplan-meier",
    Rcpp::_["method"] = "Kaplan-Meier"
  );
  return km_list;
}
')












cppFunction('
Rcpp::List calculateKaplanMeier_rcpp_old2(Rcpp::NumericVector t, Rcpp::IntegerVector d,
                                    Rcpp::NumericVector w = Rcpp::NumericVector::create(),
                                    Rcpp::IntegerVector strata = Rcpp::IntegerVector::create(),
                                    Rcpp::CharacterVector error = Rcpp::CharacterVector::create("greenwood")) {

//  Rcpp::Rcout << "Method: " << error[0] << std::endl;

  Rcpp::List km_list;

  std::vector<double> combined_times;
  std::vector<double> combined_surv;
  std::vector<int> combined_n_risk;
  std::vector<int> combined_n_event;
  std::vector<int> combined_n_censor;
//  std::vector<int> combined_unweighted_n_risk;
  std::vector<double> combined_std_err;
  std::vector<int> combined_n_stratum;
  std::vector<int> combined_u_stratum;

//  double weighted_n_event = 0;
//  double weighted_n_censor = 0;
//  double weighted_n = 0;
//  int n_event = 0;
//  int n_censor = 0;
//  int n = 0;

  if ((strata.size() == 0 || Rcpp::unique(strata).size() == 1) && (w.size() == 0 || (Rcpp::unique(w).size() == 1 && w[0] == 1))) {
  Rcpp::Rcout << "Method1: " << error[0] << std::endl;

    Rcpp::NumericVector unique_times = Rcpp::unique(t);
    std::sort(unique_times.begin(), unique_times.end());

    int u = unique_times.size();
    Rcpp::NumericVector km(u);
    Rcpp::NumericVector km_i(u);
    Rcpp::IntegerVector weighted_n_risk(u);
    Rcpp::IntegerVector weighted_n_event(u);
    Rcpp::IntegerVector weighted_n_censor(u);
//    Rcpp::IntegerVector unweighted_n_risk(u);
    Rcpp::NumericVector std_err(u);

    int n = w.size();
//    for (int j = 0; j < n; ++j) {
//      weighted_n += w[j];
//      n_event += d[j];
//    }

    for (int i = 0; i < u; ++i) {
      double time = unique_times[i];

      double weighted_n_i = 0;
      double unweighted_n_i = 0;
      for (int j = 0; j < n; ++j) {
        if (t[j] >= time) {
//          weighted_n_i += w[j];
          weighted_n_i++;
          unweighted_n_i++;
        }
      }

      weighted_n_risk[i] = weighted_n_i;
//      unweighted_n_risk[i] = unweighted_n_i;

      double weighted_events = 0;
      for (int j = 0; j < n; ++j) {
        if (t[j] == time && d[j] == 1) {
//          weighted_events += w[j];
          weighted_events++;
        }
      }

      if (weighted_n_risk[i] > 0) {
        km_i[i] = 1 - weighted_events / weighted_n_risk[i];
//        weighted_n_event += weighted_events;
        weighted_n_event[i] += weighted_events;
      } else {
        km_i[i] = 1;
      }
      if (i < u-1) {
        weighted_n_censor[i] = (weighted_n_risk[i]-weighted_n_risk[i+1]) - weighted_n_event[i];
      } else {
        weighted_n_censor[i] = weighted_n_risk[i] - weighted_n_event[i];
      }

      if (i == 0) {
        km[i] = km_i[i];
      } else {
        km[i] = km[i - 1] * km_i[i];
      }

      double sum_se = 0;
      for (int j = 0; j <= i; ++j) {
        double n_i = weighted_n_risk[j];
        double d_i = 0;
        for (int k = 0; k < n; ++k) {
          if (t[k] == unique_times[j] && d[k] == 1) {
//            d_i += w[k];
            d_i ++;
          }
        }
        if (n_i > d_i) {
          if (error[0] == "tsiatis") {
            sum_se += (d_i / (n_i * n_i));
          } else {
            sum_se += (d_i / (n_i * (n_i - d_i)));
          }
        } else {
          sum_se = std::numeric_limits<double>::infinity();
        }
      }
      std_err[i] = sqrt(sum_se);
    }

//    n_censor = n - n_event;
//    weighted_n_censor = weighted_n - weighted_n_event;

    combined_times.insert(combined_times.end(), unique_times.begin(), unique_times.end());
    combined_surv.insert(combined_surv.end(), km.begin(), km.end());
    combined_n_risk.insert(combined_n_risk.end(), weighted_n_risk.begin(), weighted_n_risk.end());
//    combined_unweighted_n_risk.insert(combined_unweighted_n_risk.end(), unweighted_n_risk.begin(), unweighted_n_risk.end());
    combined_n_event.insert(combined_n_event.end(), weighted_n_event.begin(), weighted_n_event.end());
//    combined_unweighted_n_event.insert(combined_unweighted_n_event.end(), unweighted_n_event.begin(), unweighted_n_event.end());
    combined_n_censor.insert(combined_n_censor.end(), weighted_n_censor.begin(), weighted_n_censor.end());
//    combined_unweighted_n_censor.insert(combined_unweighted_n_censor.end(), unweighted_n_censor.begin(), unweighted_n_censor.end());
    combined_std_err.insert(combined_std_err.end(), std_err.begin(), std_err.end());

  } else if (strata.size() == 0 || Rcpp::unique(strata).size() == 1) {
  Rcpp::Rcout << "Method2: " << error[0] << std::endl;

    Rcpp::NumericVector unique_times = Rcpp::unique(t);
    std::sort(unique_times.begin(), unique_times.end());

    int u = unique_times.size();
    Rcpp::NumericVector km(u);
    Rcpp::NumericVector km_i(u);
    Rcpp::IntegerVector weighted_n_risk(u);
    Rcpp::IntegerVector weighted_n_event(u);
    Rcpp::IntegerVector weighted_n_censor(u);
//    Rcpp::IntegerVector unweighted_n_risk(u);
    Rcpp::NumericVector std_err(u);

    int n = w.size();
//    for (int j = 0; j < n; ++j) {
//      weighted_n += w[j];
//      n_event += d[j];
//    }

    for (int i = 0; i < u; ++i) {
      double time = unique_times[i];

      double weighted_n_i = 0;
      double unweighted_n_i = 0;
      for (int j = 0; j < n; ++j) {
        if (t[j] >= time) {
          weighted_n_i += w[j];
          unweighted_n_i++;
        }
      }

      weighted_n_risk[i] = weighted_n_i;
//      unweighted_n_risk[i] = unweighted_n_i;

      double weighted_events = 0;
      for (int j = 0; j < n; ++j) {
        if (t[j] == time && d[j] == 1) {
          weighted_events += w[j];
        }
      }

      if (weighted_n_risk[i] > 0) {
        km_i[i] = 1 - weighted_events / weighted_n_risk[i];
//        weighted_n_event += weighted_events;
        weighted_n_event[i] += weighted_events;
      } else {
        km_i[i] = 1;
      }

      if (i < u-1) {
        weighted_n_censor[i] = (weighted_n_risk[i]-weighted_n_risk[i+1]) - weighted_n_event[i];
      } else {
        weighted_n_censor[i] = weighted_n_risk[i] - weighted_n_event[i];
      }

      if (i == 0) {
        km[i] = km_i[i];
      } else {
        km[i] = km[i - 1] * km_i[i];
      }

      double sum_se = 0;
      for (int j = 0; j <= i; ++j) {
        double n_i = weighted_n_risk[j];
        double d_i = 0;
        for (int k = 0; k < n; ++k) {
          if (t[k] == unique_times[j] && d[k] == 1) {
            d_i += w[k];
          }
        }
        if (n_i > d_i) {
          if (error[0] == "tsiatis") {
            sum_se += (d_i / (n_i * n_i));
          } else {
            sum_se += (d_i / (n_i * (n_i - d_i)));
          }
        } else {
          sum_se = std::numeric_limits<double>::infinity();
        }
      }
      std_err[i] = sqrt(sum_se);
    }

//    n_censor = n - n_event;
//    weighted_n_censor = weighted_n - weighted_n_event;

    combined_times.insert(combined_times.end(), unique_times.begin(), unique_times.end());
    combined_surv.insert(combined_surv.end(), km.begin(), km.end());
    combined_n_risk.insert(combined_n_risk.end(), weighted_n_risk.begin(), weighted_n_risk.end());
//    combined_unweighted_n_risk.insert(combined_unweighted_n_risk.end(), unweighted_n_risk.begin(), unweighted_n_risk.end());
    combined_n_event.insert(combined_n_event.end(), weighted_n_event.begin(), weighted_n_event.end());
//    combined_unweighted_n_event.insert(combined_unweighted_n_event.end(), unweighted_n_event.begin(), unweighted_n_event.end());
    combined_n_censor.insert(combined_n_censor.end(), weighted_n_censor.begin(), weighted_n_censor.end());
//    combined_unweighted_n_censor.insert(combined_unweighted_n_censor.end(), unweighted_n_censor.begin(), unweighted_n_censor.end());
    combined_std_err.insert(combined_std_err.end(), std_err.begin(), std_err.end());

  } else if (w.size() == 0 || (Rcpp::unique(w).size() == 1 && w[0] == 1)) {
  Rcpp::Rcout << "Method3: " << error[0] << std::endl;

    Rcpp::IntegerVector strata_vec = strata;

    // Loop over each strata level and calculate the Kaplan-Meier estimate
    for (int i = 0; i < Rcpp::max(strata_vec); ++i) {
      // Select the data subset corresponding to the current strata level
      Rcpp::LogicalVector strata_condition = (strata_vec == (i + 1));  // 1-based indexing for strata
      Rcpp::NumericVector t_selected = t[strata_condition];
      Rcpp::IntegerVector d_selected = d[strata_condition];
      Rcpp::NumericVector w_selected = w[strata_condition];

      // Calculate weighted_n for this strata (sum of weights in the selected dataset)
      int n_stratum = t_selected.size();
      double weighted_n_stratum = 0;
      int n_event_stratum = 0;
      for (int j = 0; j < n_stratum; ++j) {
//        weighted_n_stratum += w_selected[j];
        weighted_n_stratum += d_selected[j];
        n_event_stratum += d_selected[j];
      }

      // Calculate Kaplan-Meier for the selected subset
      Rcpp::NumericVector unique_times = Rcpp::unique(t_selected);
      std::sort(unique_times.begin(), unique_times.end());
      int u_stratum = unique_times.size();

      Rcpp::NumericVector km(u_stratum);
      Rcpp::NumericVector km_i(u_stratum);
      Rcpp::IntegerVector weighted_n_risk(u_stratum);
      Rcpp::IntegerVector weighted_n_event(u_stratum);
      Rcpp::IntegerVector weighted_n_censor(u_stratum);
//      Rcpp::IntegerVector unweighted_n_risk(u_stratum);
      Rcpp::NumericVector std_err(u_stratum);

      for (int j = 0; j < u_stratum; ++j) {
        double time = unique_times[j];

        double weighted_n_sub = 0;
//        double unweighted_n_sub = 0;
        for (int k = 0; k < t_selected.size(); ++k) {
          if (t_selected[k] >= time) {
//            weighted_n_sub += w_selected[k];
            weighted_n_sub += d_selected[k];
//            unweighted_n_sub++;
          }
        }

        weighted_n_risk[j] = weighted_n_sub;
//        unweighted_n_risk[j] = unweighted_n_sub;

        double weighted_events = 0;
        for (int k = 0; k < t_selected.size(); ++k) {
          if (t_selected[k] == time && d_selected[k] == 1) {
//            weighted_events += w_selected[k];
            weighted_events += d_selected[k];
          }
        }

        if (weighted_n_risk[j] > 0) {
          km_i[j] = 1 - weighted_events / weighted_n_risk[j];
          weighted_n_event[j] += weighted_events;
        } else {
          km_i[j] = 1;
        }
        if (j < u_stratum-1) {
          weighted_n_censor[j] = (weighted_n_risk[j]-weighted_n_risk[j+1]) - weighted_n_event[j];
        } else {
          weighted_n_censor[j] = weighted_n_risk[j] - weighted_n_event[j];
        }

        km[j] = (j == 0) ? km_i[j] : km[j - 1] * km_i[j];

        double sum_se = 0;
        for (int k = 0; k <= j; ++k) {
          double n_i = weighted_n_risk[k];
          double d_i = 0;
          for (int m = 0; m < t_selected.size(); ++m) {
            if (t_selected[m] == unique_times[k] && d_selected[m] == 1) {
//              d_i += w_selected[m];
              d_i += d_selected[m];
            }
          }
          if (n_i > d_i) {
            if (error[0] == "tsiatis") {
              sum_se += (d_i / (n_i * n_i));
            } else {
              sum_se += (d_i / (n_i * (n_i - d_i)));
            }
          } else {
            sum_se = std::numeric_limits<double>::infinity();
          }
        }
        std_err[j] = sqrt(sum_se);
      }
//      n += n_stratum;
//      weighted_n += weighted_n_stratum;
//      weighted_n_event += weighted_n_event_stratum;
//      n_event += n_event_stratum;

      combined_times.insert(combined_times.end(), unique_times.begin(), unique_times.end());
      combined_surv.insert(combined_surv.end(), km.begin(), km.end());
      combined_n_risk.insert(combined_n_risk.end(), weighted_n_risk.begin(), weighted_n_risk.end());
//      combined_unweighted_n_risk.insert(combined_unweighted_n_risk.end(), unweighted_n_risk.begin(), unweighted_n_risk.end());
      combined_n_event.insert(combined_n_event.end(), weighted_n_event.begin(), weighted_n_event.end());
//      combined_unweighted_n_event.insert(combined_unweighted_n_event.end(), unweighted_n_event.begin(), unweighted_n_event.end());
      combined_n_censor.insert(combined_n_censor.end(), weighted_n_censor.begin(), weighted_n_censor.end());
//      combined_unweighted_n_censor.insert(combined_unweighted_n_censor.end(), unweighted_n_censor.begin(), unweighted_n_censor.end());
      combined_std_err.insert(combined_std_err.end(), std_err.begin(), std_err.end());
      combined_n_stratum.insert(combined_n_stratum.end(), n_stratum);
      combined_u_stratum.insert(combined_u_stratum.end(), u_stratum);
    }
  } else {
    Rcpp::Rcout << "Method4: " << error[0] << std::endl;

    Rcpp::IntegerVector strata_vec = strata;

    // Loop over each strata level and calculate the Kaplan-Meier estimate
    for (int i = 0; i < Rcpp::max(strata_vec); ++i) {
      // Select the data subset corresponding to the current strata level
      Rcpp::LogicalVector strata_condition = (strata_vec == (i + 1));  // 1-based indexing for strata
      Rcpp::NumericVector t_selected = t[strata_condition];
      Rcpp::IntegerVector d_selected = d[strata_condition];
      Rcpp::NumericVector w_selected = w[strata_condition];

      // Calculate weighted_n for this strata (sum of weights in the selected dataset)
      int n_stratum = t_selected.size();
      double weighted_n_stratum = 0;
      int n_event_stratum = 0;
      for (int j = 0; j < n_stratum; ++j) {
        weighted_n_stratum += w_selected[j];
        n_event_stratum += d_selected[j];
      }

      // Calculate Kaplan-Meier for the selected subset
      Rcpp::NumericVector unique_times = Rcpp::unique(t_selected);
      std::sort(unique_times.begin(), unique_times.end());
      int u_stratum = unique_times.size();

      Rcpp::NumericVector km(u_stratum);
      Rcpp::NumericVector km_i(u_stratum);
      Rcpp::IntegerVector weighted_n_risk(u_stratum);
      Rcpp::IntegerVector weighted_n_event(u_stratum);
      Rcpp::IntegerVector weighted_n_censor(u_stratum);
//      Rcpp::IntegerVector unweighted_n_risk(u_stratum);
      Rcpp::NumericVector std_err(u_stratum);

      for (int j = 0; j < u_stratum; ++j) {
        double time = unique_times[j];

        double weighted_n_sub = 0;
//        double unweighted_n_sub = 0;
        for (int k = 0; k < t_selected.size(); ++k) {
          if (t_selected[k] >= time) {
            weighted_n_sub += w_selected[k];
//            unweighted_n_sub++;
          }
        }

        weighted_n_risk[j] = weighted_n_sub;
//        unweighted_n_risk[j] = unweighted_n_sub;

        double weighted_events = 0;
        for (int k = 0; k < t_selected.size(); ++k) {
          if (t_selected[k] == time && d_selected[k] == 1) {
            weighted_events += w_selected[k];
          }
        }

        if (weighted_n_risk[j] > 0) {
          km_i[j] = 1 - weighted_events / weighted_n_risk[j];
          weighted_n_event[j] += weighted_events;
        } else {
          km_i[j] = 1;
        }

        if (j < u_stratum-1) {
          weighted_n_censor[j] = (weighted_n_risk[j]-weighted_n_risk[j+1]) - weighted_n_event[j];
        } else {
          weighted_n_censor[j] = weighted_n_risk[j] - weighted_n_event[j];
        }

        km[j] = (j == 0) ? km_i[j] : km[j - 1] * km_i[j];

        double sum_se = 0;
        for (int k = 0; k <= j; ++k) {
          double n_i = weighted_n_risk[k];
          double d_i = 0;
          for (int m = 0; m < t_selected.size(); ++m) {
            if (t_selected[m] == unique_times[k] && d_selected[m] == 1) {
              d_i += w_selected[m];
            }
          }
          if (n_i > d_i) {
            if (error[0] == "tsiatis") {
              sum_se += (d_i / (n_i * n_i));
            } else {
              sum_se += (d_i / (n_i * (n_i - d_i)));
            }
          } else {
            sum_se = std::numeric_limits<double>::infinity();
          }
        }
        std_err[j] = sqrt(sum_se);
      }
//      n += n_stratum;
//      weighted_n += weighted_n_stratum;
//      weighted_n_event += weighted_n_event_stratum;
//      n_event += n_event_stratum;

      combined_times.insert(combined_times.end(), unique_times.begin(), unique_times.end());
      combined_surv.insert(combined_surv.end(), km.begin(), km.end());
      combined_n_risk.insert(combined_n_risk.end(), weighted_n_risk.begin(), weighted_n_risk.end());
//      combined_unweighted_n_risk.insert(combined_unweighted_n_risk.end(), unweighted_n_risk.begin(), unweighted_n_risk.end());
      combined_n_event.insert(combined_n_event.end(), weighted_n_event.begin(), weighted_n_event.end());
//      combined_unweighted_n_event.insert(combined_unweighted_n_event.end(), unweighted_n_event.begin(), unweighted_n_event.end());
      combined_n_censor.insert(combined_n_censor.end(), weighted_n_censor.begin(), weighted_n_censor.end());
//      combined_unweighted_n_censor.insert(combined_unweighted_n_censor.end(), unweighted_n_censor.begin(), unweighted_n_censor.end());
      combined_std_err.insert(combined_std_err.end(), std_err.begin(), std_err.end());
      combined_n_stratum.insert(combined_n_stratum.end(), n_stratum);
      combined_u_stratum.insert(combined_u_stratum.end(), u_stratum);
    }
  }

//  n_censor = n - n_event;
//  weighted_n_censor = weighted_n - weighted_n_event;

  Rcpp::NumericVector all_times = Rcpp::wrap(combined_times);
  Rcpp::NumericVector all_surv = Rcpp::wrap(combined_surv);
  Rcpp::IntegerVector all_n_risk = Rcpp::wrap(combined_n_risk);
//  Rcpp::IntegerVector all_unweighted_n_risk = Rcpp::wrap(combined_unweighted_n_risk);
  Rcpp::IntegerVector all_n_event = Rcpp::wrap(combined_n_event);
//  Rcpp::IntegerVector all_unweighted_n_event = Rcpp::wrap(combined_unweighted_n_event);
  Rcpp::IntegerVector all_n_censor = Rcpp::wrap(combined_n_censor);
//  Rcpp::IntegerVector all_unweighted_n_censor = Rcpp::wrap(combined_unweighted_n_censor);
  Rcpp::NumericVector all_std_err = Rcpp::wrap(combined_std_err);
  Rcpp::IntegerVector all_n_stratum = Rcpp::wrap(combined_n_stratum);
  Rcpp::IntegerVector all_u_stratum = Rcpp::wrap(combined_u_stratum);

  // Create final combined output for all strata
  km_list = Rcpp::List::create(
    Rcpp::_["time"] = all_times,
    Rcpp::_["surv"] = all_surv,
    Rcpp::_["n.risk"] = all_n_risk,
//    Rcpp::_["unweighted.n.risk"] = all_unweighted_n_risk,
    Rcpp::_["n"] = combined_n_stratum,
    Rcpp::_["n.event"] = all_n_event,
    Rcpp::_["n.censor"] = all_n_censor,
//    Rcpp::_["unweighted.n"] = all_n_stratum,
//    Rcpp::_["unweighted.n.event"] = all_unweighted_n_event,
//    Rcpp::_["unweighted.n.censor"] = all_unweighted_n_censor,
    Rcpp::_["std.err"] = all_std_err,
    Rcpp::_["high"] = R_NilValue,
    Rcpp::_["low"] = R_NilValue,
    Rcpp::_["conf.type"] = R_NilValue,
    Rcpp::_["strata"] = all_u_stratum,
    Rcpp::_["type"] = "kaplan-meier",
    Rcpp::_["method"] = "Kaplan-Meier"
  );
  return km_list;
}
')

















cppFunction('
Rcpp::List calculateKaplanMeier_rcpp_old(Rcpp::NumericVector t, Rcpp::IntegerVector d, Rcpp::NumericVector w,
                                     Rcpp::IntegerVector strata = Rcpp::IntegerVector::create(),
                                     Rcpp::CharacterVector error = Rcpp::CharacterVector::create("greenwood")) {

//  Rcpp::Rcout << "Method: " << error[0] << std::endl;

  Rcpp::List km_list;

  std::vector<double> combined_times;
  std::vector<double> combined_surv;
  std::vector<int> combined_n_risk;
  std::vector<int> combined_n_event;
  std::vector<int> combined_n_censor;
//  std::vector<int> combined_unweighted_n_risk;
  std::vector<double> combined_std_err;
  std::vector<int> combined_n_stratum;
  std::vector<int> combined_u_stratum;

//  double weighted_n_event = 0;
//  double weighted_n_censor = 0;
//  double weighted_n = 0;
//  int n_event = 0;
//  int n_censor = 0;
//  int n = 0;

  if (strata.size() == 0 || Rcpp::unique(strata).size() == 1) {

    Rcpp::NumericVector unique_times = Rcpp::unique(t);
    std::sort(unique_times.begin(), unique_times.end());

    int u = unique_times.size();
    Rcpp::NumericVector km(u);
    Rcpp::NumericVector km_i(u);
    Rcpp::IntegerVector weighted_n_risk(u);
    Rcpp::IntegerVector weighted_n_event(u);
    Rcpp::IntegerVector weighted_n_censor(u);
//    Rcpp::IntegerVector unweighted_n_risk(u);
    Rcpp::NumericVector std_err(u);

    int n = w.size();
//    for (int j = 0; j < n; ++j) {
//      weighted_n += w[j];
//      n_event += d[j];
//    }

    for (int i = 0; i < u; ++i) {
      double time = unique_times[i];

      double weighted_n_i = 0;
      double unweighted_n_i = 0;
      for (int j = 0; j < n; ++j) {
        if (t[j] >= time) {
          weighted_n_i += w[j];
          unweighted_n_i++;
        }
      }

      weighted_n_risk[i] = weighted_n_i;
//      unweighted_n_risk[i] = unweighted_n_i;

      double weighted_events = 0;
      for (int j = 0; j < n; ++j) {
        if (t[j] == time && d[j] == 1) {
          weighted_events += w[j];
        }
      }

      if (weighted_n_risk[i] > 0) {
        km_i[i] = 1 - weighted_events / weighted_n_risk[i];
//        weighted_n_event += weighted_events;
        weighted_n_event[i] += weighted_events;
      } else {
        km_i[i] = 1;
      }
      if (i < u-1) {
        weighted_n_censor[i] = (weighted_n_risk[i]-weighted_n_risk[i+1]) - weighted_n_event[i];
      } else {
        weighted_n_censor[i] = weighted_n_risk[i] - weighted_n_event[i];
      }

      if (i == 0) {
        km[i] = km_i[i];
      } else {
        km[i] = km[i - 1] * km_i[i];
      }

      double sum_se = 0;
      for (int j = 0; j <= i; ++j) {
        double n_i = weighted_n_risk[j];
        double d_i = 0;
        for (int k = 0; k < n; ++k) {
          if (t[k] == unique_times[j] && d[k] == 1) {
            d_i += w[k];
          }
        }
        if (n_i > d_i) {
          if (error[0] == "tsiatis") {
            sum_se += (d_i / (n_i * n_i));
          } else {
            sum_se += (d_i / (n_i * (n_i - d_i)));
          }
        } else {
          sum_se = std::numeric_limits<double>::infinity();
        }
      }
      std_err[i] = sqrt(sum_se);
    }

//    n_censor = n - n_event;
//    weighted_n_censor = weighted_n - weighted_n_event;

    combined_times.insert(combined_times.end(), unique_times.begin(), unique_times.end());
    combined_surv.insert(combined_surv.end(), km.begin(), km.end());
    combined_n_risk.insert(combined_n_risk.end(), weighted_n_risk.begin(), weighted_n_risk.end());
//    combined_unweighted_n_risk.insert(combined_unweighted_n_risk.end(), unweighted_n_risk.begin(), unweighted_n_risk.end());
    combined_n_event.insert(combined_n_event.end(), weighted_n_event.begin(), weighted_n_event.end());
//    combined_unweighted_n_event.insert(combined_unweighted_n_event.end(), unweighted_n_event.begin(), unweighted_n_event.end());
    combined_n_censor.insert(combined_n_censor.end(), weighted_n_censor.begin(), weighted_n_censor.end());
//    combined_unweighted_n_censor.insert(combined_unweighted_n_censor.end(), unweighted_n_censor.begin(), unweighted_n_censor.end());
    combined_std_err.insert(combined_std_err.end(), std_err.begin(), std_err.end());

  } else {
    Rcpp::IntegerVector strata_vec = strata;

    // Loop over each strata level and calculate the Kaplan-Meier estimate
    for (int i = 0; i < Rcpp::max(strata_vec); ++i) {
      // Select the data subset corresponding to the current strata level
      Rcpp::LogicalVector strata_condition = (strata_vec == (i + 1));  // 1-based indexing for strata
      Rcpp::NumericVector t_selected = t[strata_condition];
      Rcpp::IntegerVector d_selected = d[strata_condition];
      Rcpp::NumericVector w_selected = w[strata_condition];

      // Calculate weighted_n for this strata (sum of weights in the selected dataset)
      int n_stratum = t_selected.size();
      double weighted_n_stratum = 0;
      int n_event_stratum = 0;
      for (int j = 0; j < n_stratum; ++j) {
        weighted_n_stratum += w_selected[j];
        n_event_stratum += d_selected[j];
      }

      // Calculate Kaplan-Meier for the selected subset
      Rcpp::NumericVector unique_times = Rcpp::unique(t_selected);
      std::sort(unique_times.begin(), unique_times.end());
      int u_stratum = unique_times.size();

      Rcpp::NumericVector km(u_stratum);
      Rcpp::NumericVector km_i(u_stratum);
      Rcpp::IntegerVector weighted_n_risk(u_stratum);
      Rcpp::IntegerVector weighted_n_event(u_stratum);
      Rcpp::IntegerVector weighted_n_censor(u_stratum);
//      Rcpp::IntegerVector unweighted_n_risk(u_stratum);
      Rcpp::NumericVector std_err(u_stratum);

      for (int j = 0; j < u_stratum; ++j) {
        double time = unique_times[j];

        double weighted_n_sub = 0;
//        double unweighted_n_sub = 0;
        for (int k = 0; k < t_selected.size(); ++k) {
          if (t_selected[k] >= time) {
            weighted_n_sub += w_selected[k];
//            unweighted_n_sub++;
          }
        }

        weighted_n_risk[j] = weighted_n_sub;
//        unweighted_n_risk[j] = unweighted_n_sub;

        double weighted_events = 0;
        for (int k = 0; k < t_selected.size(); ++k) {
          if (t_selected[k] == time && d_selected[k] == 1) {
            weighted_events += w_selected[k];
          }
        }

        if (weighted_n_risk[j] > 0) {
          km_i[j] = 1 - weighted_events / weighted_n_risk[j];
          weighted_n_event[j] += weighted_events;
        } else {
          km_i[j] = 1;
        }
        if (j < u_stratum-1) {
          weighted_n_censor[j] = (weighted_n_risk[j]-weighted_n_risk[j+1]) - weighted_n_event[j];
        } else {
          weighted_n_censor[j] = weighted_n_risk[j] - weighted_n_event[j];
        }

        km[j] = (j == 0) ? km_i[j] : km[j - 1] * km_i[j];

        double sum_se = 0;
        for (int k = 0; k <= j; ++k) {
          double n_i = weighted_n_risk[k];
          double d_i = 0;
          for (int m = 0; m < t_selected.size(); ++m) {
            if (t_selected[m] == unique_times[k] && d_selected[m] == 1) {
              d_i += w_selected[m];
            }
          }
          if (n_i > d_i) {
            if (error[0] == "tsiatis") {
              sum_se += (d_i / (n_i * n_i));
            } else {
              sum_se += (d_i / (n_i * (n_i - d_i)));
            }
          } else {
            sum_se = std::numeric_limits<double>::infinity();
          }
        }
        std_err[j] = sqrt(sum_se);
      }
//      n += n_stratum;
//      weighted_n += weighted_n_stratum;
//      weighted_n_event += weighted_n_event_stratum;
//      n_event += n_event_stratum;

      combined_times.insert(combined_times.end(), unique_times.begin(), unique_times.end());
      combined_surv.insert(combined_surv.end(), km.begin(), km.end());
      combined_n_risk.insert(combined_n_risk.end(), weighted_n_risk.begin(), weighted_n_risk.end());
//      combined_unweighted_n_risk.insert(combined_unweighted_n_risk.end(), unweighted_n_risk.begin(), unweighted_n_risk.end());
      combined_n_event.insert(combined_n_event.end(), weighted_n_event.begin(), weighted_n_event.end());
//      combined_unweighted_n_event.insert(combined_unweighted_n_event.end(), unweighted_n_event.begin(), unweighted_n_event.end());
      combined_n_censor.insert(combined_n_censor.end(), weighted_n_censor.begin(), weighted_n_censor.end());
//      combined_unweighted_n_censor.insert(combined_unweighted_n_censor.end(), unweighted_n_censor.begin(), unweighted_n_censor.end());
      combined_std_err.insert(combined_std_err.end(), std_err.begin(), std_err.end());
      combined_n_stratum.insert(combined_n_stratum.end(), n_stratum);
      combined_u_stratum.insert(combined_u_stratum.end(), u_stratum);
    }
  }

//  n_censor = n - n_event;
//  weighted_n_censor = weighted_n - weighted_n_event;

  Rcpp::NumericVector all_times = Rcpp::wrap(combined_times);
  Rcpp::NumericVector all_surv = Rcpp::wrap(combined_surv);
  Rcpp::IntegerVector all_n_risk = Rcpp::wrap(combined_n_risk);
//  Rcpp::IntegerVector all_unweighted_n_risk = Rcpp::wrap(combined_unweighted_n_risk);
  Rcpp::IntegerVector all_n_event = Rcpp::wrap(combined_n_event);
//  Rcpp::IntegerVector all_unweighted_n_event = Rcpp::wrap(combined_unweighted_n_event);
  Rcpp::IntegerVector all_n_censor = Rcpp::wrap(combined_n_censor);
//  Rcpp::IntegerVector all_unweighted_n_censor = Rcpp::wrap(combined_unweighted_n_censor);
  Rcpp::NumericVector all_std_err = Rcpp::wrap(combined_std_err);
  Rcpp::IntegerVector all_n_stratum = Rcpp::wrap(combined_n_stratum);
  Rcpp::IntegerVector all_u_stratum = Rcpp::wrap(combined_u_stratum);

  // Create final combined output for all strata
  km_list = Rcpp::List::create(
    Rcpp::_["time"] = all_times,
    Rcpp::_["surv"] = all_surv,
    Rcpp::_["n.risk"] = all_n_risk,
//    Rcpp::_["unweighted.n.risk"] = all_unweighted_n_risk,
    Rcpp::_["n"] = combined_n_stratum,
    Rcpp::_["n.event"] = all_n_event,
    Rcpp::_["n.censor"] = all_n_censor,
    Rcpp::_["unweighted.n"] = all_n_stratum,
//    Rcpp::_["unweighted.n.event"] = all_unweighted_n_event,
//    Rcpp::_["unweighted.n.censor"] = all_unweighted_n_censor,
    Rcpp::_["std.err"] = all_std_err,
    Rcpp::_["high"] = R_NilValue,
    Rcpp::_["low"] = R_NilValue,
    Rcpp::_["conf.type"] = "log-log",
    Rcpp::_["strata"] = all_u_stratum,
    Rcpp::_["type"] = "kaplan-meier",
    Rcpp::_["method"] = "Kaplan-Meier"
  );
  return km_list;
}
')



