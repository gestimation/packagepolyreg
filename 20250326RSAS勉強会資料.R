###########################################################################################
# calculateKaplanMeier()の出力をsurvfit()の出力形式に合わせて修正
# survfit$survは全観測値ではなく全時点の推定値を返すことに注意

  calculateKaplanMeier_unique_times <- function(t, d, w){
    n <- length(t)
    data <- data.frame(t = t, d = d, w = w, id = 1:n)
    sorted_data <- data[order(data$t), ]
    sorted_t <- sorted_data$t
    sorted_d <- sorted_data$d
    sorted_w <- sorted_data$w
    sorted_id <- sorted_data$id

    unique_times <- unique(sorted_t)
    u <- length(unique_times)
    km <- numeric(u)
    n.risk <- numeric(u)
    km_i <- numeric(u)

    for (i in 1:u) {
      time <- unique_times[i]
      risk_set <- (sorted_t >= time)
      n.risk[i] <- sum(risk_set)
      events_at_time <- sum(sorted_t == time & sorted_d == 1)
      weighted_events <- sum(sorted_w[sorted_t == time & sorted_d == 1])
      if (n.risk[i] > 0) {
        km_i[i] <- 1 - weighted_events / n.risk[i]
      } else {
        km_i[i] <- 1
      }
      if (i == 1) {
        km[i] <- km_i[i]
      } else {
        km[i] <- km[i-1] * km_i[i]
      }
    }
    km_object <- list(
      time = unique_times,
      surv = km,
      n = n,
      n.risk = n.risk,
      n.event = sum(d),
      n.censor = n-sum(d),
      std.err = NULL,
      high = NULL,
      low = NULL,
      conf.type = NULL,
      strata = NULL,
      call = match.call(),
      type = "kaplan-meier",
      method = "Kaplan-Meier"
    )
    return(km_object)
  }

###########################################################################################
# 層別・重み付きKaplan-Meier推定量に対応

  calculateKaplanMeier_strata <- function(t, d, strata = NULL, subset = NULL, weight = NULL){
    km_list <- NULL
    ifelse(is.null(weight), w <- rep(1,length(t)), w <- weight)
    if (is.null(strata)) {
      km_object <- calculateKaplanMeier_unique_times(t, d, w)
      class(km_object) <- "survfit"
      return(km_object)
    } else {
      for (level in levels(as.factor(strata))) {
        ifelse(is.null(subset), index <- (strata == level), index <- (strata == level & (subset == 1)))
        t_selected <- t[index]
        d_selected <- d[index]
        w_selected <- w[index]
        km_selected <- calculateKaplanMeier_unique_times(t_selected, d_selected, w_selected)
        km_selected$strata <- level
        km_list[[level]] <- km_selected
      }
    }
    combined_km_object <- list(
      time = unique(unlist(lapply(km_list, function(x) x$time))),
      surv = unlist(lapply(km_list, function(x) x$surv)),
      n = sum(sapply(km_list, function(x) x$n)),
      n.risk = unlist(lapply(km_list, function(x) x$n.risk)),
      n.event = unlist(lapply(km_list, function(x) x$n.event)),
      n.censor = unlist(lapply(km_list, function(x) x$n.censor)),
      std.err = unlist(lapply(km_list, function(x) x$std.err)),
      high = unlist(lapply(km_list, function(x) x$high)),
      low = unlist(lapply(km_list, function(x) x$low)),
      conf.type = NULL,
      strata = sapply(km_list, function(x) x$strata),
      call = match.call(),
      type = "kaplan-meier",
      method = "Kaplan-Meier"
    )
    class(combined_km_object) <- "survfit"
    return(combined_km_object)
  }


###########################################################################################
# テスト用のデータセットを作成

library(survival)
n_test <- 20
w_test <- rep(1, n_test)
t_test <- c(1:(n_test/2),1:(n_test/2))
strata_test <- as.factor((t_test %% 2 == 0))
epsilon_test <- rep(1, n_test)
epsilon_test[1] <- 1
epsilon_test[2] <- 1
epsilon_test[n_test/2-1] <- 1
epsilon_test[n_test-1] <- 1
epsilon_test[n_test/2] <- 1
epsilon_test[n_test/2+1] <- 1
epsilon_test[n_test] <- 0
d_test <- as.numeric(epsilon_test>0)
df_test <- data.frame(id = 1:n_test, t_test = t_test, epsilon_test = epsilon_test, d_test = d_test, strata_test = strata_test)
expected <- survfit(Surv(t_test, d_test) ~ strata_test, data = df_test)
print(expected$surv)
print(expected$time)
print(expected$n.risk)
print(expected$strata)

tested <- calculateKaplanMeier_strata(t_test, as.numeric(epsilon_test>0), strata_test)
print(tested$surv)
print(tested$time)
print(tested$n.risk)

###########################################################################################
expected <- survfit(Surv(t_test, d_test) ~ +1, data = df_test)
print(expected$surv)
print(expected$time)
print(expected$n.risk)
print(epsilon_test)

tested <- calculateKaplanMeier_unique_times(t_test, as.numeric(epsilon_test>0), w_test)
tested <- calculateKaplanMeier_strata(t_test, as.numeric(epsilon_test>0))
print(tested$n.event)
print(tested$surv)
print(tested$time)

###########################################################################################
# データセットが大きくなるにつれsurvfit()の計算時間の方が速くなる

library(microbenchmark)
microbenchmark(survfit(Surv(t_test, d_test) ~ strata(strata_test), data = df_test),
               calculateKaplanMeier_strata(t_test, d_test, strata_test), times = 100)

###########################################################################################
# Rccp関数を用いてC++で実装すると速くなる

library(Rcpp)


cppFunction('
Rcpp::List calculateKaplanMeier_rcpp(Rcpp::NumericVector t, Rcpp::IntegerVector d, Rcpp::NumericVector w,
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
        // SE(S_t) = S_t * sqrt(sum(d_i / (n_i * (n_i - d_i))))
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
t_test <- 1:n_test
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

tested <- calculateKaplanMeier_rcpp(t_test, as.numeric(epsilon_test>0), w_test, strata_test)
print(tested)

expected <- survfit(Surv(t_test, d_test) ~ strata_test, data = df_test, error="tsiatis")
print(expected$std.err)
print(expected$surv)
expected <- survfit(Surv(t_test, d_test) ~ +1, data = df_test, error="greenwood")
print(expected$std.err)


library(microbenchmark)
microbenchmark(survfit(Surv(t_test, d_test) ~ strata(strata_test), data = df_test),
              calculateKaplanMeier_rcpp(t_test, as.numeric(epsilon_test>0), strata_test), times = 100)














cppFunction('
Rcpp::List calculateKaplanMeier_rcpp(Rcpp::NumericVector t, Rcpp::IntegerVector d, Rcpp::NumericVector w,
                                     Rcpp::IntegerVector strata = Rcpp::IntegerVector::create(),
                                     Rcpp::IntegerVector subset = Rcpp::IntegerVector::create()) {

  int n = t.size();
  Rcpp::List km_list;

  // Handle weight, default to 1 if not provided
  if (w.size() == 0) {
    w = Rcpp::rep(1.0, n);
  }

  // Determine if stratification is needed
  bool stratified = (strata.size() > 0 && Rcpp::unique(strata).size() > 1);

  if (!stratified) {
    // If strata is not provided or has only one level, calculate Kaplan-Meier without strata
    Rcpp::NumericVector unique_times = Rcpp::unique(t);
    std::sort(unique_times.begin(), unique_times.end());
    int u = unique_times.size();

    // Apply subset condition regardless of whether strata is defined
    Rcpp::LogicalVector subset_condition = Rcpp::LogicalVector(n, true);
    if (subset.size() > 0) {
      // Use the subset condition if provided
      Rcpp::IntegerVector subset_vec = subset;
      subset_condition = (subset_vec == 1);
    }

    // Apply subset condition to the data
    Rcpp::NumericVector t_selected = t[subset_condition];
    Rcpp::IntegerVector d_selected = d[subset_condition];
    Rcpp::NumericVector w_selected = w[subset_condition];

    // Preallocate memory for the result vectors
    Rcpp::NumericVector km(u);
    Rcpp::NumericVector std_err(u);
    Rcpp::IntegerVector weighted_n_risk(u);
    Rcpp::IntegerVector unweighted_n_risk(u);
    double weighted_n_event = 0;

    int idx = 0;
    for (int i = 0; i < u; ++i) {
      double time = unique_times[i];
      double weighted_n = 0;
      double unweighted_n = 0;
      double weighted_events = 0;

      // Find the relevant entries for this unique time
      while (idx < t_selected.size() && t_selected[idx] <= time) {
        if (d_selected[idx] == 1) {
          weighted_events += w_selected[idx];  // Event weight
        }
        weighted_n += w_selected[idx];  // Risk set size (weighted)
        unweighted_n++;  // Risk set size (unweighted)
        idx++;
      }

      weighted_n_risk[i] = weighted_n;
      unweighted_n_risk[i] = unweighted_n;

      if (weighted_n_risk[i] > 0) {
        km[i] = (i == 0) ? 1 - weighted_events / weighted_n_risk[i] : km[i - 1] * (1 - weighted_events / weighted_n_risk[i]);
        weighted_n_event += weighted_events;
      } else {
        km[i] = 1;
      }

            double sum_se = 0.0;
            for (int j = 0; j <= i; ++j) {
              if (weighted_n_risk[j] > 0) {
                sum_se += (d_selected[j] / (weighted_n_risk[j] * (weighted_n_risk[j] - d_selected[j])));
              }
            }
            std_err[i] = km[i];
    }

    double weighted_n_censor = std::accumulate(w_selected.begin(), w_selected.end(), 0.0) - weighted_n_event;

    Rcpp::List km_object = Rcpp::List::create(
      Rcpp::_["time"] = unique_times,
      Rcpp::_["surv"] = km,
      Rcpp::_["n"] = t_selected.size(),
      Rcpp::_["n.risk"] = weighted_n_risk,
      Rcpp::_["unweighted.n.risk"] = unweighted_n_risk,
      Rcpp::_["n.event"] = weighted_n_event,
      Rcpp::_["n.censor"] = weighted_n_censor,
      Rcpp::_["std.err"] = std_err,  // Standard error added
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
    int max_strata = Rcpp::max(strata);

    // Loop over each strata level
    for (int i = 0; i < max_strata; ++i) {
      // Select the data subset corresponding to the current strata level
      Rcpp::LogicalVector strata_condition = (strata == (i + 1));  // 1-based indexing for strata

      // Apply the subset condition
      Rcpp::LogicalVector subset_condition = Rcpp::LogicalVector(n, true);
      if (subset.size() > 0) {
        // Use the subset condition if provided
        Rcpp::IntegerVector subset_vec = subset;
        subset_condition = (subset_vec == 1);
      }

      Rcpp::LogicalVector final_condition = strata_condition & subset_condition;

      // Subset the data using the condition
      Rcpp::NumericVector t_selected = t[final_condition];
      Rcpp::IntegerVector d_selected = d[final_condition];
      Rcpp::NumericVector w_selected = w[final_condition];

      // Get unique times and sort them
      Rcpp::NumericVector unique_times = Rcpp::unique(t_selected);
      std::sort(unique_times.begin(), unique_times.end());
      int u = unique_times.size();

      // Preallocate memory for the result vectors
      Rcpp::NumericVector km(u);
      Rcpp::NumericVector std_err(u);
      Rcpp::IntegerVector weighted_n_risk(u);
      Rcpp::IntegerVector unweighted_n_risk(u);
      double weighted_n_event = 0;

      int idx = 0;
      for (int j = 0; j < u; ++j) {
        double time = unique_times[j];
        double weighted_n = 0;
        double unweighted_n = 0;
        double weighted_events = 0;

        // Find the relevant entries for this unique time
        while (idx < t_selected.size() && t_selected[idx] <= time) {
          if (d_selected[idx] == 1) {
            weighted_events += w_selected[idx];  // Event weight
          }
          weighted_n += w_selected[idx];  // Risk set size (weighted)
          unweighted_n++;  // Risk set size (unweighted)
          idx++;
        }

        weighted_n_risk[j] = weighted_n;
        unweighted_n_risk[j] = unweighted_n;

        if (weighted_n_risk[j] > 0) {
          km[j] = (j == 0) ? 1 - weighted_events / weighted_n_risk[j] : km[j - 1] * (1 - weighted_events / weighted_n_risk[j]);
          weighted_n_event += weighted_events;
        } else {
          km[j] = 1;
        }

            double sum_se = 0.0;
            for (int k = 0; k <= j; ++k) {
              if (weighted_n_risk[k] > 0) {
                sum_se += (d_selected[k] / (weighted_n_risk[k] * (weighted_n_risk[k] - d_selected[k])));
              }
            }
            std_err[j] = km[j];
      }

      double weighted_n_censor = std::accumulate(w_selected.begin(), w_selected.end(), 0.0) - weighted_n_event;

      Rcpp::List km_object = Rcpp::List::create(
        Rcpp::_["time"] = unique_times,
        Rcpp::_["surv"] = km,
        Rcpp::_["n"] = t_selected.size(),
        Rcpp::_["n.risk"] = weighted_n_risk,
        Rcpp::_["unweighted.n.risk"] = unweighted_n_risk,
        Rcpp::_["n.event"] = weighted_n_event,
        Rcpp::_["n.censor"] = weighted_n_censor,
        Rcpp::_["std.err"] = std_err,  // Standard error added
        Rcpp::_["high"] = R_NilValue,
        Rcpp::_["low"] = R_NilValue,
        Rcpp::_["conf.type"] = R_NilValue,
        Rcpp::_["strata"] = i + 1,  // Add strata level
        Rcpp::_["type"] = "kaplan-meier",
        Rcpp::_["method"] = "Kaplan-Meier"
      );

      km_list.push_back(km_object);
    }
  }

  return km_list;
}
')
tested <- calculateKaplanMeier_rcpp(t_test, as.numeric(epsilon_test>0), w_test)
print(tested)





##################################################################################################################################
df_test <- createTestData(20, 2, first_zero=TRUE, last_zero=FALSE, subset_present=FALSE, logical_strata=FALSE, na_strata=TRUE)
print(df_test)
expected <- survfit(Surv(t, d)~strata, df_test, weights=w)
print(expected$surv)
print(expected$time)
print(expected$n)
print(expected$n.risk)
print(expected$n.event)
print(expected$n.censor)
print(expected$strata)
print(expected$std.err)
print(expected$low)
print(expected$high)

tested <- km.curve(Surv(t, d)~strata, df_test, conf.type="log-log", na.action=na.omit, weights=w)
print(tested$surv)
print(tested$time)
print(tested$n)
print(tested$n.risk)
print(tested$n.event)
print(tested$n.censor)
print(tested$strata)
print(tested$std.err)
print(tested$low)
print(tested$high)


##################################################################################################################################
library(microbenchmark)
df_test <- createTestData(1000, 1, first_zero=TRUE, last_zero=FALSE, subset_present=FALSE, logical_strata=FALSE, na_strata=TRUE)
microbenchmark(survfit(Surv(t, d)~strata, df_test),
               km.curve(Surv(t, d)~strata, df_test, na.action=na.omit),
               km.curve(Surv(t, d)~strata, df_test, conf.type=NULL, na.action=na.omit),
               times = 200)



###########################################################################################
# 層別なしKaplan-Meier推定量

cppFunction('
Rcpp::List calculateKaplanMeier_rcpp(Rcpp::NumericVector t, Rcpp::IntegerVector d, Rcpp::NumericVector w) {
  int n = t.size();

  Rcpp::IntegerVector indices = Rcpp::seq(0, n-1); // Create index vector
  std::vector<std::tuple<double, int, double, int>> sorted_data(n);
  for (int i = 0; i < n; ++i) {
    sorted_data[i] = std::make_tuple(t[i], d[i], w[i], indices[i]);
  }
  std::sort(sorted_data.begin(), sorted_data.end(), [](const std::tuple<double, int, double, int>& a,
                                                       const std::tuple<double, int, double, int>& b) {
    return std::get<0>(a) < std::get<0>(b);  // Sort by time (first element)
  });

  // Reassign the sorted values back to t, d, w, and indices, based on sorted order
  for (int i = 0; i < n; ++i) {
    t[i] = std::get<0>(sorted_data[i]);
    d[i] = std::get<1>(sorted_data[i]);
    w[i] = std::get<2>(sorted_data[i]);
    indices[i] = std::get<3>(sorted_data[i]);
  }

  // Get unique times
  Rcpp::NumericVector unique_times = Rcpp::unique(t);
  std::sort(unique_times.begin(), unique_times.end());
  int u = unique_times.size();

  // Initialize vectors to store the survival probabilities and at-risk set sizes
  Rcpp::NumericVector km(u);
  Rcpp::NumericVector km_i(u);
  Rcpp::IntegerVector weighted_n_risk(u);
  Rcpp::IntegerVector unweighted_n_risk(u);
  double weighted_n_event = 0;
  double weighted_n_censor = 0;

  for (int i = 0; i < u; ++i) {
    double time = unique_times[i];

    double weighted_n = 0;
    double unweighted_n = 0;
    for (int j = 0; j < n; ++j) {
      if (t[j] >= time) {
        weighted_n += w[j];  // Weight the individuals in the risk set
        unweighted_n++;
      }
    }
    weighted_n_risk[i] = weighted_n;
    unweighted_n_risk[i] = unweighted_n;

    double weighted_events = 0;
    for (int j = 0; j < n; ++j) {
      if (t[j] == time && d[j] == 1) {
        weighted_events += w[j];  // Weight the event occurrence
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
    Rcpp::_["std.err"] = R_NilValue,  // NULL equivalent in C++
    Rcpp::_["high"] = R_NilValue,     // NULL equivalent in C++
    Rcpp::_["low"] = R_NilValue,      // NULL equivalent in C++
    Rcpp::_["conf.type"] = R_NilValue, // NULL equivalent in C++
    Rcpp::_["strata"] = R_NilValue,   // NULL equivalent in C++
    // Rcpp::_["call"] = Rcpp::wrap(Rcpp::call("match.call")()),
    Rcpp::_["type"] = "kaplan-meier",
    Rcpp::_["method"] = "Kaplan-Meier"
  );
  return km_object;
}
')

