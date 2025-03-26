#include <Rcpp.h>
#include <RcppParallel.h>
#include <algorithm>
#include <numeric>

// Define the parallel loop for Kaplan-Meier calculation
struct KaplanMeierWorker : public RcppParallel::Worker {
    // Input vectors
    Rcpp::NumericVector t;
    Rcpp::IntegerVector d;
    Rcpp::NumericVector w;
    int n;
    Rcpp::NumericVector unique_times;
    
    // Output vectors
    Rcpp::NumericVector km;
    Rcpp::NumericVector km_i;
    Rcpp::IntegerVector weighted_n_risk;
    Rcpp::IntegerVector unweighted_n_risk;
    Rcpp::NumericVector std_err;

    // Constructor
    KaplanMeierWorker(Rcpp::NumericVector t, Rcpp::IntegerVector d, Rcpp::NumericVector w, 
                      Rcpp::NumericVector unique_times, Rcpp::NumericVector km, Rcpp::NumericVector km_i, 
                      Rcpp::IntegerVector weighted_n_risk, Rcpp::IntegerVector unweighted_n_risk, 
                      Rcpp::NumericVector std_err, int n) 
        : t(t), d(d), w(w), unique_times(unique_times), km(km), km_i(km_i), 
          weighted_n_risk(weighted_n_risk), unweighted_n_risk(unweighted_n_risk), std_err(std_err), n(n) {}

    // This is the function that will be called for each "chunk" of the data
    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; ++i) {
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
            } else {
                km_i[i] = 1;
            }

            km[i] = (i == 0) ? km_i[i] : km[i - 1] * km_i[i];

            // Standard error calculation using Greenwood formula
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
    }
};

// Main function for Kaplan-Meier estimation
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
        Rcpp::NumericVector unique_times = Rcpp::unique(t);
        std::sort(unique_times.begin(), unique_times.end());
        int u = unique_times.size();

        Rcpp::NumericVector km(u);
        Rcpp::NumericVector km_i(u);
        Rcpp::IntegerVector weighted_n_risk(u);
        Rcpp::IntegerVector unweighted_n_risk(u);
        Rcpp::NumericVector std_err(u);  // Standard error vector

        // Create the worker object for parallelization
        KaplanMeierWorker worker(t, d, w, unique_times, km, km_i, weighted_n_risk, unweighted_n_risk, std_err, n);

        // Use RcppParallel to run the calculation in parallel
        RcppParallel::parallelFor(0, u, worker);

        // Calculating weighted and unweighted number of events at each unique time
        Rcpp::NumericVector n_event(u);
        for (int i = 0; i < u; ++i) {
            double weighted_events = 0;
            for (int j = 0; j < n; ++j) {
                if (t[j] == unique_times[i] && d[j] == 1) {
                    weighted_events += w[j];
                }
            }
            n_event[i] = weighted_events;
        }

        Rcpp::List km_object = Rcpp::List::create(
            Rcpp::_["time"] = unique_times,
            Rcpp::_["surv"] = km,
            Rcpp::_["n"] = n,
            Rcpp::_["n.risk"] = weighted_n_risk,
            Rcpp::_["unweighted.n.risk"] = unweighted_n_risk,
            Rcpp::_["n.event"] = n_event,
            Rcpp::_["n.censor"] = weighted_n_risk - n_event,
            Rcpp::_["std.err"] = std_err,
            Rcpp::_["type"] = "kaplan-meier",
            Rcpp::_["method"] = "Kaplan-Meier"
        );

        km_list.push_back(km_object);
    } else {
        // Perform stratified analysis (similar to before)
        // Add your stratified calculation code here...
    }

    return km_list;
}
