km.curve <- function(formula,
                     data,
                     weights=NULL,
                     subset=NULL,
                     conf.int=.95,
                     error=NULL,
                     conf.type=NULL,
                     conf.lower=NULL
) {
  cl <- match.call()
  mf <- match.call(expand.dots = TRUE)[1:3]
  special <- NULL
  out_terms <- terms(formula, special, data = data)
  mf$formula <- out_terms
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  Y <- model.extract(mf, "response")
  t <- as.numeric(Y[, 1])  # time variable
  d <- as.integer(Y[, 2])  # status variable
  strata_name <- all.vars(out_terms)[3]
  strata <- as.integer(as.factor(data[[strata_name]]))
  #  w <- data[[weights]]
  w <- rep(1, nrow(data))
  subset <- as.integer(rep(1, nrow(data)))
#  strata <- rep(1, nrow(data))
  out_calculateKaplanMeier_rcpp <- calculateKaplanMeier_rcpp(t, d, w, strata, subset)
  km_object <- list(
    time = out_calculateKaplanMeier_rcpp$time,
    surv = out_calculateKaplanMeier_rcpp$surv,
    n = out_calculateKaplanMeier_rcpp$n,
    n.risk = out_calculateKaplanMeier_rcpp$n.risk,
    n.event = out_calculateKaplanMeier_rcpp$n.event,
    n.censor = out_calculateKaplanMeier_rcpp$n.censor,
    std.err = out_calculateKaplanMeier_rcpp$std.err,
    high = NULL,
    low = NULL,
    conf.type = conf.type,
    strata = NULL,
    call = match.call(),
    type = "kaplan-meier",
    method = "Kaplan-Meier"
  )
  return(km_object)
}

library(survival)
library(Rcpp)

cppFunction('
Rcpp::List calculateKaplanMeier_rcpp_old(Rcpp::NumericVector t, Rcpp::IntegerVector d, Rcpp::NumericVector w,
                                     Rcpp::IntegerVector strata = Rcpp::IntegerVector::create(),
                                     Rcpp::IntegerVector subset = Rcpp::IntegerVector::create()) {

  int n = t.size();
  Rcpp::List km_list;

  // Handle weight, default to 1 if not provided
  if (w.size() == 0) {
    w = Rcpp::rep(1.0, n);
  }

  // If strata is not provided or has only one level, calculate Kaplan-Meier without strata
  if (strata.size() == 0 || Rcpp::unique(strata).size() == 1) {
    // Perform standard Kaplan-Meier calculation
    Rcpp::NumericVector unique_times = Rcpp::unique(t);
    std::sort(unique_times.begin(), unique_times.end());
    int u = unique_times.size();

    Rcpp::NumericVector km(u);
    Rcpp::NumericVector km_i(u);
    Rcpp::IntegerVector weighted_n_risk(u);
    Rcpp::IntegerVector unweighted_n_risk(u);
    Rcpp::NumericVector std_err(u);  // Standard error vector
    double weighted_n_event = 0;
    double weighted_n_censor = 0;

    for (int i = 0; i < u; ++i) {
      double time = unique_times[i];

      double weighted_n = 0;
      double unweighted_n = 0;
      for (int j = 0; j < n; ++j) {
        if (t[j] >= time) {
          weighted_n += w[j];
          unweighted_n++;
        }
      }
      weighted_n_risk[i] = weighted_n;
      unweighted_n_risk[i] = unweighted_n;

      double weighted_events = 0;
      for (int j = 0; j < n; ++j) {
        if (t[j] == time && d[j] == 1) {
          weighted_events += w[j];
        }
      }

      if (weighted_n_risk[i] > 0) {
        km_i[i] = 1 - weighted_events / weighted_n_risk[i];
        weighted_n_event += weighted_events;
      } else {
        km_i[i] = 1;
      }

      if (i == 0) {
        km[i] = km_i[i];
      } else {
        km[i] = km[i - 1] * km_i[i];
      }

      double sum_se = 0;
      for (int j = 0; j <= i; ++j) {  // Summing over all times <= current time (t_j)
        double n_i = weighted_n_risk[j];
        double d_i = 0;
        for (int k = 0; k < n; ++k) {
          if (t[k] == unique_times[j] && d[k] == 1) {
            d_i += w[k];
          }
        }
        if (n_i > d_i) {
          sum_se += d_i / (n_i * (n_i - d_i));
        }
      }
      std_err[i] = sqrt(sum_se);
    }

    weighted_n_censor = std::accumulate(w.begin(), w.end(), 0.0) - weighted_n_event;

    Rcpp::List km_object = Rcpp::List::create(
      Rcpp::_["time"] = unique_times,
      Rcpp::_["surv"] = km,
      Rcpp::_["n"] = n,
      Rcpp::_["n.risk"] = weighted_n_risk,
      Rcpp::_["unweighted.n.risk"] = unweighted_n_risk,
      Rcpp::_["n.event"] = weighted_n_event,
      Rcpp::_["n.censor"] = weighted_n_censor,
      Rcpp::_["std.err"] = std_err,
      Rcpp::_["high"] = R_NilValue,
      Rcpp::_["low"] = R_NilValue,
      Rcpp::_["conf.type"] = R_NilValue,
      Rcpp::_["strata"] = R_NilValue,
      Rcpp::_["type"] = "kaplan-meier",
      Rcpp::_["method"] = "Kaplan-Meier"
    );

    km_list.push_back(km_object);

  } else {
    // If strata is provided, perform stratified analysis
    Rcpp::IntegerVector strata_vec = strata;

    // Loop over each strata level
    for (int i = 0; i < Rcpp::max(strata_vec); ++i) {
      // Select the data subset corresponding to the current strata level
      Rcpp::LogicalVector strata_condition = (strata_vec == (i + 1));  // 1-based indexing for strata
      Rcpp::LogicalVector subset_condition = Rcpp::LogicalVector(n, true);

      if (subset.size() > 0) {
        // Use the subset condition if provided
        Rcpp::IntegerVector subset_vec = subset;
        subset_condition = (subset_vec == 1);
      }

      Rcpp::LogicalVector final_condition = strata_condition & subset_condition;

      Rcpp::NumericVector t_selected = t[final_condition];
      Rcpp::IntegerVector d_selected = d[final_condition];
      Rcpp::NumericVector w_selected = w[final_condition];

      // Calculate Kaplan-Meier for the selected subset (same code as before)
      Rcpp::NumericVector unique_times = Rcpp::unique(t_selected);
      std::sort(unique_times.begin(), unique_times.end());
      int u = unique_times.size();

      Rcpp::NumericVector km(u);
      Rcpp::NumericVector km_i(u);
      Rcpp::IntegerVector weighted_n_risk(u);
      Rcpp::IntegerVector unweighted_n_risk(u);
      Rcpp::NumericVector std_err(u);
      double weighted_n_event = 0;
      double weighted_n_censor = 0;

      for (int j = 0; j < u; ++j) {
        double time = unique_times[j];

        double weighted_n = 0;
        double unweighted_n = 0;
        for (int k = 0; k < t_selected.size(); ++k) {
          if (t_selected[k] >= time) {
            weighted_n += w_selected[k];
            unweighted_n++;
          }
        }
        weighted_n_risk[j] = weighted_n;
        unweighted_n_risk[j] = unweighted_n;

        double weighted_events = 0;
        for (int k = 0; k < t_selected.size(); ++k) {
          if (t_selected[k] == time && d_selected[k] == 1) {
            weighted_events += w_selected[k];
          }
        }

        if (weighted_n_risk[j] > 0) {
          km_i[j] = 1 - weighted_events / weighted_n_risk[j];
          weighted_n_event += weighted_events;
        } else {
          km_i[j] = 1;
        }

        km[j] = (j == 0) ? km_i[j] : km[j - 1] * km_i[j];

        // Calculate the standard error using the correct Greenwood formula
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
            sum_se += (d_i / (n_i * (n_i - d_i)));
          }
        }
        std_err[j] = sqrt(sum_se);  // Apply the standard error formula
      }

      weighted_n_censor = std::accumulate(w_selected.begin(), w_selected.end(), 0.0) - weighted_n_event;

      // Create the km_object for this strata
      km_list.push_back(Rcpp::List::create(
        Rcpp::_["time"] = unique_times,
        Rcpp::_["surv"] = km,
        Rcpp::_["n"] = t_selected.size(),
        Rcpp::_["n.risk"] = weighted_n_risk,
        Rcpp::_["unweighted.n.risk"] = unweighted_n_risk,
        Rcpp::_["n.event"] = weighted_n_event,
        Rcpp::_["n.censor"] = weighted_n_censor,
        Rcpp::_["std.err"] = std_err,
        Rcpp::_["high"] = R_NilValue,
        Rcpp::_["low"] = R_NilValue,
        Rcpp::_["conf.type"] = R_NilValue,
        Rcpp::_["strata"] = i + 1,
        Rcpp::_["type"] = "kaplan-meier",
        Rcpp::_["method"] = "Kaplan-Meier"
      ));
    }
  }

  return km_list;
}
')


cppFunction('
Rcpp::List calculateKaplanMeier_rcpp_old2(Rcpp::NumericVector t, Rcpp::IntegerVector d, Rcpp::NumericVector w,
                                     Rcpp::IntegerVector strata = Rcpp::IntegerVector::create(),
                                     Rcpp::IntegerVector subset = Rcpp::IntegerVector::create()) {

  int n = t.size();
  Rcpp::List km_list;

  // Handle weight, default to 1 if not provided
  if (w.size() == 0) {
    w = Rcpp::rep(1.0, n);
  }

  // Variables to accumulate results across all strata
  std::vector<double> combined_times;
  std::vector<double> combined_surv;
  std::vector<int> combined_n_risk;
  std::vector<int> combined_unweighted_n_risk;
  std::vector<double> combined_std_err;

  double weighted_n_event = 0;
  double weighted_n_censor = 0;
  double weighted_n = 0;

  // If strata is not provided or has only one level, calculate Kaplan-Meier without strata
  if (strata.size() == 0 || Rcpp::unique(strata).size() == 1) {

    Rcpp::NumericVector unique_times = Rcpp::unique(t);
    std::sort(unique_times.begin(), unique_times.end());
    int u = unique_times.size();

    Rcpp::NumericVector km(u);
    Rcpp::NumericVector km_i(u);
    Rcpp::IntegerVector weighted_n_risk(u);
    Rcpp::IntegerVector unweighted_n_risk(u);
    Rcpp::NumericVector std_err(u);

    // Calculate weighted_n (sum of weights in the selected dataset)
    for (int j = 0; j < n; ++j) {
      weighted_n += w[j];  // Sum of weights for all subjects
    }

    for (int i = 0; i < u; ++i) {
      double time = unique_times[i];

      double weighted_n_stratum = 0;
      double unweighted_n_stratum = 0;
      for (int j = 0; j < n; ++j) {
        if (t[j] >= time) {
          weighted_n_stratum += w[j];
          unweighted_n_stratum++;
        }
      }

      weighted_n_risk[i] = weighted_n_stratum;
      unweighted_n_risk[i] = unweighted_n_stratum;

      double weighted_events = 0;
      for (int j = 0; j < n; ++j) {
        if (t[j] == time && d[j] == 1) {
          weighted_events += w[j];
        }
      }

      if (weighted_n_risk[i] > 0) {
        km_i[i] = 1 - weighted_events / weighted_n_risk[i];
        weighted_n_event += weighted_events;
      } else {
        km_i[i] = 1;
      }

      if (i == 0) {
        km[i] = km_i[i];
      } else {
        km[i] = km[i - 1] * km_i[i];
      }

      // Calculate weighted_n_censor as the difference between weighted_n and weighted_n_event
      weighted_n_censor = weighted_n - weighted_n_event;

      // Calculate the standard error using the correct Greenwood formula
      double sum_se = 0;
      for (int j = 0; j <= i; ++j) {  // Summing over all times <= current time (t_j)
        double n_i = weighted_n_risk[j];
        double d_i = 0;
        for (int k = 0; k < n; ++k) {
          if (t[k] == unique_times[j] && d[k] == 1) {
            d_i += w[k];
          }
        }
        if (n_i > d_i) {
          sum_se += (d_i / (n_i * (n_i - d_i)));
        }
      }
      std_err[i] = sqrt(sum_se);
    }

    // Add results to combined vectors
    combined_times.insert(combined_times.end(), unique_times.begin(), unique_times.end());
    combined_surv.insert(combined_surv.end(), km.begin(), km.end());
    combined_n_risk.insert(combined_n_risk.end(), weighted_n_risk.begin(), weighted_n_risk.end());
    combined_unweighted_n_risk.insert(combined_unweighted_n_risk.end(), unweighted_n_risk.begin(), unweighted_n_risk.end());
    combined_std_err.insert(combined_std_err.end(), std_err.begin(), std_err.end());

  } else {
    Rcpp::IntegerVector strata_vec = strata;
    double weighted_n_event = 0;

    // Loop over each strata level and calculate the Kaplan-Meier estimate
    for (int i = 0; i < Rcpp::max(strata_vec); ++i) {

      double weighted_n_stratum = 0;
      double weighted_n_event_stratum = 0;
      int n_event_stratum = 0;
      int n_stratum;

      // Select the data subset corresponding to the current strata level
      Rcpp::LogicalVector strata_condition = (strata_vec == (i + 1));  // 1-based indexing for strata
      Rcpp::LogicalVector subset_condition = Rcpp::LogicalVector(n, true);

      if (subset.size() > 0) {
        // Use the subset condition if provided
        Rcpp::IntegerVector subset_vec = subset;
        subset_condition = (subset_vec == 1);
      }

      Rcpp::LogicalVector final_condition = strata_condition & subset_condition;

      Rcpp::NumericVector t_selected = t[final_condition];
      Rcpp::IntegerVector d_selected = d[final_condition];
      Rcpp::NumericVector w_selected = w[final_condition];

      // Calculate weighted_n for this strata (sum of weights in the selected dataset)
      n_stratum = t_selected.size();
      for (int j = 0; j < n_stratum; ++j) {
        weighted_n_stratum += w_selected[j];  // Sum of weights for the strata
      }

      // Calculate Kaplan-Meier for the selected subset
      Rcpp::NumericVector unique_times = Rcpp::unique(t_selected);
      std::sort(unique_times.begin(), unique_times.end());
      int u = unique_times.size();

      Rcpp::NumericVector km(u);
      Rcpp::NumericVector km_i(u);
      Rcpp::IntegerVector weighted_n_risk(u);
      Rcpp::IntegerVector unweighted_n_risk(u);
      Rcpp::NumericVector std_err(u);

      for (int j = 0; j < u; ++j) {
        double time = unique_times[j];

        double weighted_n_sub = 0;
        double unweighted_n_sub = 0;
        for (int k = 0; k < t_selected.size(); ++k) {
          if (t_selected[k] >= time) {
            weighted_n_sub += w_selected[k];
            unweighted_n_sub++;
          }
        }

        weighted_n_risk[j] = weighted_n_sub;
        unweighted_n_risk[j] = unweighted_n_sub;

        double weighted_events = 0;
        for (int k = 0; k < t_selected.size(); ++k) {
          if (t_selected[k] == time && d_selected[k] == 1) {
            weighted_events += w_selected[k];
          }
        }

        if (weighted_n_risk[j] > 0) {
          km_i[j] = 1 - weighted_events / weighted_n_risk[j];
          weighted_n_event_stratum += weighted_events;
        } else {
          km_i[j] = 1;
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
            sum_se += (d_i / (n_i * (n_i - d_i)));
          }
        }
        std_err[j] = sqrt(sum_se);
      }

      weighted_n_event += weighted_n_event_stratum;
      weighted_n_censor = weighted_n - weighted_n_event;

      // Add results to combined vectors
      combined_times.insert(combined_times.end(), unique_times.begin(), unique_times.end());
      combined_surv.insert(combined_surv.end(), km.begin(), km.end());
      combined_n_risk.insert(combined_n_risk.end(), weighted_n_risk.begin(), weighted_n_risk.end());
      combined_unweighted_n_risk.insert(combined_unweighted_n_risk.end(), unweighted_n_risk.begin(), unweighted_n_risk.end());
      combined_std_err.insert(combined_std_err.end(), std_err.begin(), std_err.end());
    }
  }

  Rcpp::NumericVector all_times = Rcpp::wrap(combined_times);
  Rcpp::NumericVector all_surv = Rcpp::wrap(combined_surv);
  Rcpp::IntegerVector all_n_risk = Rcpp::wrap(combined_n_risk);
  Rcpp::IntegerVector all_unweighted_n_risk = Rcpp::wrap(combined_unweighted_n_risk);
  Rcpp::NumericVector all_std_err = Rcpp::wrap(combined_std_err);

  // Create final combined output for all strata
  km_list = Rcpp::List::create(
    Rcpp::_["time"] = all_times,
    Rcpp::_["surv"] = all_surv,
    Rcpp::_["n.risk"] = all_n_risk,
    Rcpp::_["unweighted.n.risk"] = all_unweighted_n_risk,
    Rcpp::_["n"] = weighted_n,
    Rcpp::_["n.event"] = weighted_n_event,
    Rcpp::_["n.censor"] = weighted_n_censor,
    Rcpp::_["unweighted.n"] = n,
    Rcpp::_["unweighted.n.event"] = weighted_n_event,
    Rcpp::_["unweighted.n.censor"] = weighted_n_censor,
    Rcpp::_["std.err"] = all_std_err,
    Rcpp::_["high"] = R_NilValue,
    Rcpp::_["low"] = R_NilValue,
    Rcpp::_["conf.type"] = "log-log",
    Rcpp::_["strata"] = R_NilValue,
    Rcpp::_["type"] = "kaplan-meier",
    Rcpp::_["method"] = "Kaplan-Meier"
  );
  return km_list;
}
')









cppFunction('
Rcpp::List calculateKaplanMeier_rcpp(Rcpp::NumericVector t, Rcpp::IntegerVector d, Rcpp::NumericVector w,
                                     Rcpp::IntegerVector strata = Rcpp::IntegerVector::create(),
                                     Rcpp::IntegerVector subset = Rcpp::IntegerVector::create()) {

  int n_ = t.size();
  Rcpp::List km_list;

  // Handle weight, default to 1 if not provided
  if (w.size() == 0) {
    w = Rcpp::rep(1.0, n_);
  }

  // Variables to accumulate results across all strata
  std::vector<double> combined_times;
  std::vector<double> combined_surv;
  std::vector<int> combined_n_risk;
  std::vector<int> combined_unweighted_n_risk;
  std::vector<double> combined_std_err;

  double weighted_n_event = 0;
  double weighted_n_censor = 0;
  double weighted_n = 0;
  int n_event = 0;
  int n_censor = 0;
  int n = 0;

  if (strata.size() == 0 || Rcpp::unique(strata).size() == 1) {

    Rcpp::LogicalVector subset_condition = Rcpp::LogicalVector(n, true);
    if (subset.size() > 0) {
      Rcpp::IntegerVector subset_vec = subset;
      subset_condition = (subset_vec == 1);
    }
    Rcpp::NumericVector t_selected = t[subset_condition];
    Rcpp::IntegerVector d_selected = d[subset_condition];
    Rcpp::NumericVector w_selected = w[subset_condition];

    Rcpp::NumericVector unique_times = Rcpp::unique(t_selected);
    std::sort(unique_times.begin(), unique_times.end());
    int u = unique_times.size();

    Rcpp::NumericVector km(u);
    Rcpp::NumericVector km_i(u);
    Rcpp::IntegerVector weighted_n_risk(u);
    Rcpp::IntegerVector unweighted_n_risk(u);
    Rcpp::NumericVector std_err(u);

    // Calculate weighted_n (sum of weights in the selected dataset)
    n = w_selected.size();
    for (int j = 0; j < n; ++j) {
      weighted_n += w_selected[j];
      n_event += d_selected[j];
    }

    for (int i = 0; i < u; ++i) {
      double time = unique_times[i];

      double weighted_n_i = 0;
      double unweighted_n_i = 0;
      for (int j = 0; j < n; ++j) {
        if (t[j] >= time) {
          weighted_n_i += w_selected[j];
          unweighted_n_i++;
        }
      }

      weighted_n_risk[i] = weighted_n_i;
      unweighted_n_risk[i] = unweighted_n_i;

      double weighted_events = 0;
      for (int j = 0; j < n; ++j) {
        if (t_selected[j] == time && d_selected[j] == 1) {
          weighted_events += w_selected[j];
        }
      }

      if (weighted_n_risk[i] > 0) {
        km_i[i] = 1 - weighted_events / weighted_n_risk[i];
        weighted_n_event += weighted_events;
      } else {
        km_i[i] = 1;
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
          if (t[k] == unique_times[j] && d_selected[k] == 1) {
            d_i += w_selected[k];
          }
        }
        if (n_i > d_i) {
          sum_se += (d_i / (n_i * (n_i - d_i)));
        }
      }
      std_err[i] = sqrt(sum_se);
    }

    n_censor = n - n_event;
    weighted_n_censor = weighted_n - weighted_n_event;

    combined_times.insert(combined_times.end(), unique_times.begin(), unique_times.end());
    combined_surv.insert(combined_surv.end(), km.begin(), km.end());
    combined_n_risk.insert(combined_n_risk.end(), weighted_n_risk.begin(), weighted_n_risk.end());
    combined_unweighted_n_risk.insert(combined_unweighted_n_risk.end(), unweighted_n_risk.begin(), unweighted_n_risk.end());
    combined_std_err.insert(combined_std_err.end(), std_err.begin(), std_err.end());

  } else {
    Rcpp::IntegerVector strata_vec = strata;

    // Loop over each strata level and calculate the Kaplan-Meier estimate
    for (int i = 0; i < Rcpp::max(strata_vec); ++i) {
      // Select the data subset corresponding to the current strata level
      Rcpp::LogicalVector strata_condition = (strata_vec == (i + 1));  // 1-based indexing for strata
      Rcpp::LogicalVector subset_condition = Rcpp::LogicalVector(n, true);

      if (subset.size() > 0) {
        Rcpp::IntegerVector subset_vec = subset;
        subset_condition = (subset_vec == 1);
      }

      Rcpp::LogicalVector final_condition = strata_condition & subset_condition;

      Rcpp::NumericVector t_selected = t[final_condition];
      Rcpp::IntegerVector d_selected = d[final_condition];
      Rcpp::NumericVector w_selected = w[final_condition];

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
      int u = unique_times.size();

      Rcpp::NumericVector km(u);
      Rcpp::NumericVector km_i(u);
      Rcpp::IntegerVector weighted_n_risk(u);
      Rcpp::IntegerVector unweighted_n_risk(u);
      Rcpp::NumericVector std_err(u);
      double weighted_n_event_stratum = 0;

      for (int j = 0; j < u; ++j) {
        double time = unique_times[j];

        double weighted_n_sub = 0;
        double unweighted_n_sub = 0;
        for (int k = 0; k < t_selected.size(); ++k) {
          if (t_selected[k] >= time) {
            weighted_n_sub += w_selected[k];
            unweighted_n_sub++;
          }
        }

        weighted_n_risk[j] = weighted_n_sub;
        unweighted_n_risk[j] = unweighted_n_sub;

        double weighted_events = 0;
        for (int k = 0; k < t_selected.size(); ++k) {
          if (t_selected[k] == time && d_selected[k] == 1) {
            weighted_events += w_selected[k];
          }
        }

        if (weighted_n_risk[j] > 0) {
          km_i[j] = 1 - weighted_events / weighted_n_risk[j];
          weighted_n_event += weighted_events;
        } else {
          km_i[j] = 1;
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
            sum_se += (d_i / (n_i * (n_i - d_i)));
          }
        }
        std_err[j] = sqrt(sum_se);
      }

      n += n_stratum;
      weighted_n += weighted_n_stratum;
      weighted_n_event += weighted_n_event_stratum;
      n_event += n_event_stratum;

      combined_times.insert(combined_times.end(), unique_times.begin(), unique_times.end());
      combined_surv.insert(combined_surv.end(), km.begin(), km.end());
      combined_n_risk.insert(combined_n_risk.end(), weighted_n_risk.begin(), weighted_n_risk.end());
      combined_unweighted_n_risk.insert(combined_unweighted_n_risk.end(), unweighted_n_risk.begin(), unweighted_n_risk.end());
      combined_std_err.insert(combined_std_err.end(), std_err.begin(), std_err.end());
    }
  }

  n_censor = n - n_event;
  weighted_n_censor = weighted_n - weighted_n_event;

  Rcpp::NumericVector all_times = Rcpp::wrap(combined_times);
  Rcpp::NumericVector all_surv = Rcpp::wrap(combined_surv);
  Rcpp::IntegerVector all_n_risk = Rcpp::wrap(combined_n_risk);
  Rcpp::IntegerVector all_unweighted_n_risk = Rcpp::wrap(combined_unweighted_n_risk);
  Rcpp::NumericVector all_std_err = Rcpp::wrap(combined_std_err);

  // Create final combined output for all strata
  km_list = Rcpp::List::create(
    Rcpp::_["time"] = all_times,
    Rcpp::_["surv"] = all_surv,
    Rcpp::_["n.risk"] = all_n_risk,
    Rcpp::_["unweighted.n.risk"] = all_unweighted_n_risk,
    Rcpp::_["n"] = weighted_n,
    Rcpp::_["n.event"] = weighted_n_event,
    Rcpp::_["n.censor"] = weighted_n_censor,
    Rcpp::_["unweighted.n"] = n,
    Rcpp::_["unweighted.n.event"] = n_event,
    Rcpp::_["unweighted.n.censor"] = n_censor,
    Rcpp::_["std.err"] = all_std_err,
    Rcpp::_["high"] = R_NilValue,
    Rcpp::_["low"] = R_NilValue,
    Rcpp::_["conf.type"] = "log-log",
    Rcpp::_["strata"] = R_NilValue,
    Rcpp::_["type"] = "kaplan-meier",
    Rcpp::_["method"] = "Kaplan-Meier"
  );
  return km_list;
}
')


n_test <- 20
one <- rep(1, n_test)
w_test <- rep(2, n_test)
subset_test <- rep(1, n_test)
subset_test[1] <- 0
t_test <- c(1:(n_test/2),1:(n_test/2))
strata_test <- as.factor((t_test %% 2 == 0))
epsilon_test <- rep(1, n_test)
epsilon_test[1] <- 0
epsilon_test[2] <- 0
epsilon_test[n_test/2-1] <- 0
epsilon_test[n_test-1] <- 0
epsilon_test[n_test/2] <- 0
epsilon_test[n_test/2+1] <- 0
epsilon_test[n_test] <- 0
d_test <- as.numeric(epsilon_test>0)
df_test <- data.frame(id = 1:n_test, t_test = t_test, epsilon_test = epsilon_test, d_test = d_test, strata_test = strata_test)

tested <- calculateKaplanMeier_rcpp(t_test, d_test, w_test, strata_test, subset_test)
#tested <- calculateKaplanMeier_rcpp(t_test, d_test, w_test, one, subset_test)
print(tested)


#tested <- km.curve(Surv(t_test, d_test)~strata_test, df_test)

