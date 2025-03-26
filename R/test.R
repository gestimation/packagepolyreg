library(Rcpp)
library(RcppParallel)




cppFunction('
Rcpp::List calculateKaplanMeier_rcpp1(Rcpp::NumericVector t, Rcpp::IntegerVector d, Rcpp::NumericVector w,
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

    // Precompute weighted events and risk for all unique times
    for (int i = 0; i < u; ++i) {
      double time = unique_times[i];
      double weighted_n = 0;
      double unweighted_n = 0;
      double weighted_events = 0;

      // Efficiently compute weighted_n, unweighted_n, and weighted_events in one loop
      for (int j = 0; j < n; ++j) {
        if (t[j] >= time) {
          weighted_n += w[j];
          unweighted_n++;
        }
        if (t[j] == time && d[j] == 1) {
          weighted_events += w[j];
        }
      }

      weighted_n_risk[i] = weighted_n;
      unweighted_n_risk[i] = unweighted_n;

      if (weighted_n_risk[i] > 0) {
        km_i[i] = 1 - weighted_events / weighted_n_risk[i];
        weighted_n_event += weighted_events;
      } else {
        km_i[i] = 1;
      }

      km[i] = (i == 0) ? km_i[i] : km[i - 1] * km_i[i];

      // Efficient standard error calculation using precomputed values
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
          sum_se += (d_i / (n_i * (n_i - d_i)));
        }
      }
      std_err[i] = sqrt(sum_se);
    }

    weighted_n_censor = std::accumulate(w.begin(), w.end(), 0.0) - weighted_n_event;

    km_list.push_back(Rcpp::List::create(
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
    ));
  } else {
    // Stratified analysis
    for (int i = 0; i < Rcpp::max(strata); ++i) {
      Rcpp::LogicalVector strata_condition = (strata == (i + 1));  // 1-based indexing for strata
      Rcpp::LogicalVector subset_condition = Rcpp::LogicalVector(n, true);

      if (subset.size() > 0) {
        Rcpp::IntegerVector subset_vec = subset;
        subset_condition = (subset_vec == 1);
      }

      Rcpp::LogicalVector final_condition = strata_condition & subset_condition;
      Rcpp::NumericVector t_selected = t[final_condition];
      Rcpp::IntegerVector d_selected = d[final_condition];
      Rcpp::NumericVector w_selected = w[final_condition];

      Rcpp::NumericVector unique_times = Rcpp::unique(t_selected);
      std::sort(unique_times.begin(), unique_times.end());
      int u = unique_times.size();

      Rcpp::NumericVector km(u);
      Rcpp::NumericVector km_i(u);
      Rcpp::IntegerVector weighted_n_risk(u);
      Rcpp::IntegerVector unweighted_n_risk(u);
      Rcpp::NumericVector std_err(u);  // Standard error vector
      double weighted_n_event = 0;

      // Precompute weighted events and risk for the selected strata subset
      for (int j = 0; j < u; ++j) {
        double time = unique_times[j];
        double weighted_n = 0;
        double unweighted_n = 0;
        double weighted_events = 0;

        for (int k = 0; k < t_selected.size(); ++k) {
          if (t_selected[k] >= time) {
            weighted_n += w_selected[k];
            unweighted_n++;
          }
          if (t_selected[k] == time && d_selected[k] == 1) {
            weighted_events += w_selected[k];
          }
        }

        weighted_n_risk[j] = weighted_n;
        unweighted_n_risk[j] = unweighted_n;

        if (weighted_n_risk[j] > 0) {
          km_i[j] = 1 - weighted_events / weighted_n_risk[j];
          weighted_n_event += weighted_events;
        } else {
          km_i[j] = 1;
        }

        km[j] = (j == 0) ? km_i[j] : km[j - 1] * km_i[j];

        // Efficient standard error calculation
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

      double weighted_n_censor = std::accumulate(w_selected.begin(), w_selected.end(), 0.0) - weighted_n_event;

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
Rcpp::List calculateKaplanMeier_rcpp2(Rcpp::NumericVector t, Rcpp::IntegerVector d, Rcpp::NumericVector w,
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

    // Precompute weighted events and risk for all unique times
    for (int i = 0; i < u; ++i) {
      double time = unique_times[i];
      double weighted_n = 0;
      double unweighted_n = 0;
      double weighted_events = 0;

      // Efficiently compute weighted_n, unweighted_n, and weighted_events in one loop
      for (int j = 0; j < n; ++j) {
        if (t[j] >= time) {
          weighted_n += w[j];
          unweighted_n++;
        }
        if (t[j] == time && d[j] == 1) {
          weighted_events += w[j];
        }
      }

      weighted_n_risk[i] = weighted_n;
      unweighted_n_risk[i] = unweighted_n;

      if (weighted_n_risk[i] > 0) {
        km_i[i] = 1 - weighted_events / weighted_n_risk[i];
        weighted_n_event += weighted_events;
      } else {
        km_i[i] = 1;
      }

      km[i] = (i == 0) ? km_i[i] : km[i - 1] * km_i[i];

      // Efficient standard error calculation using precomputed values
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
          sum_se += (d_i / (n_i * (n_i - d_i)));
        }
      }
      std_err[i] = sqrt(sum_se);
    }

    weighted_n_censor = std::accumulate(w.begin(), w.end(), 0.0) - weighted_n_event;

    km_list.push_back(Rcpp::List::create(
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
    ));
  } else {
    // Stratified analysis
    for (int i = 0; i < Rcpp::max(strata); ++i) {
      Rcpp::LogicalVector strata_condition = (strata == (i + 1));  // 1-based indexing for strata
      Rcpp::LogicalVector subset_condition = Rcpp::LogicalVector(n, true);

      if (subset.size() > 0) {
        Rcpp::IntegerVector subset_vec = subset;
        subset_condition = (subset_vec == 1);
      }

      Rcpp::LogicalVector final_condition = strata_condition & subset_condition;
      Rcpp::NumericVector t_selected = t[final_condition];
      Rcpp::IntegerVector d_selected = d[final_condition];
      Rcpp::NumericVector w_selected = w[final_condition];

      Rcpp::NumericVector unique_times = Rcpp::unique(t_selected);
      std::sort(unique_times.begin(), unique_times.end());
      int u = unique_times.size();

      Rcpp::NumericVector km(u);
      Rcpp::NumericVector km_i(u);
      Rcpp::IntegerVector weighted_n_risk(u);
      Rcpp::IntegerVector unweighted_n_risk(u);
      Rcpp::NumericVector std_err(u);  // Standard error vector
      double weighted_n_event = 0;

      // Precompute weighted events and risk for the selected strata subset
      for (int j = 0; j < u; ++j) {
        double time = unique_times[j];
        double weighted_n = 0;
        double unweighted_n = 0;
        double weighted_events = 0;

        for (int k = 0; k < t_selected.size(); ++k) {
          if (t_selected[k] >= time) {
            weighted_n += w_selected[k];
            unweighted_n++;
          }
          if (t_selected[k] == time && d_selected[k] == 1) {
            weighted_events += w_selected[k];
          }
        }

        weighted_n_risk[j] = weighted_n;
        unweighted_n_risk[j] = unweighted_n;

        if (weighted_n_risk[j] > 0) {
          km_i[j] = 1 - weighted_events / weighted_n_risk[j];
          weighted_n_event += weighted_events;
        } else {
          km_i[j] = 1;
        }

        km[j] = (j == 0) ? km_i[j] : km[j - 1] * km_i[j];

        // Efficient standard error calculation
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

      double weighted_n_censor = std::accumulate(w_selected.begin(), w_selected.end(), 0.0) - weighted_n_event;

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
Rcpp::List calculateKaplanMeier_rcpp3(Rcpp::NumericVector t, Rcpp::IntegerVector d, Rcpp::NumericVector w,
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
      Rcpp::_["std.err"] = std_err,  // Add standard error to the result
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
      Rcpp::NumericVector std_err(u);  // Standard error vector
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
        Rcpp::_["std.err"] = std_err,  // Add standard error to the result
        Rcpp::_["high"] = R_NilValue,
        Rcpp::_["low"] = R_NilValue,
        Rcpp::_["conf.type"] = R_NilValue,
        Rcpp::_["strata"] = i + 1, // Add strata level
        Rcpp::_["type"] = "kaplan-meier",
        Rcpp::_["method"] = "Kaplan-Meier"
      ));
    }
  }

  return km_list;
}
')


n_test <- 20
w_test <- rep(1, n_test)
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

tested <- calculateKaplanMeier_rcpp3(t_test, as.numeric(epsilon_test>0), w_test)
print(tested)


#library(microbenchmark)
#microbenchmark(survfit(Surv(t_test, d_test) ~ strata(strata_test), data = df_test),
#             calculateKaplanMeier_rcpp(t_test, as.numeric(epsilon_test>0), strata_test),
#            calculateKaplanMeier_rcpp3(t_test, as.numeric(epsilon_test>0), strata_test),
#           times = 20)

