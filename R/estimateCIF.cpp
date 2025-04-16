// estimateCIF.cpp
#include <Rcpp.h>
#include <algorithm>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame estimateCIFcpp(NumericVector t,
                         NumericVector epsilon,
                         double code_event1,
                         double code_event2,
                         double code_censoring) {
  int n = t.size();

  // インデックス（0～n-1）を作成し、t の昇順に並び替えるための順序を求める
  std::vector<int> idx(n);
  for (int i = 0; i < n; i++) {
    idx[i] = i;
  }
  std::sort(idx.begin(), idx.end(), [&](int i, int j) {
    return t[i] < t[j];
  });

  // ソート済みの各ベクトルを作成
  NumericVector sorted_t(n);
  NumericVector sorted_epsilon(n);
  NumericVector sorted_epsilon_all(n);
  IntegerVector sorted_id(n);

  for (int i = 0; i < n; i++) {
    int ii = idx[i];
    sorted_t[i] = t[ii];
    sorted_epsilon[i] = epsilon[ii];
    sorted_id[i] = ii + 1; // R の id は 1 から始まる
    // 打ち切りなら 0、それ以外なら 1
    sorted_epsilon_all[i] = (epsilon[ii] == code_censoring) ? 0 : 1;
  }

  // リスクセットのサイズを計算
  NumericVector n_atrisk(n);
  for (int i = 0; i < n; i++) {
    n_atrisk[i] = n - i;
  }

  // Kaplan-Meier 推定量の計算
  NumericVector km(n);
  double cum_log = 0.0;
  for (int i = 0; i < n; i++) {
    double s = 1.0 - sorted_epsilon_all[i] / n_atrisk[i];
    cum_log += std::log(s);
    km[i] = std::exp(cum_log);
  }

  // イベントごとのインジケータ作成
  NumericVector epsilon1(n);
  NumericVector epsilon2(n);
  for (int i = 0; i < n; i++) {
    epsilon1[i] = (sorted_epsilon[i] == code_event1) ? 1.0 : 0.0;
    epsilon2[i] = (sorted_epsilon[i] == code_event2) ? 1.0 : 0.0;
  }

  // lagged km を計算：最初は 1、それ以降は前の km の値
  NumericVector km_lag(n);
  km_lag[0] = 1.0;
  for (int i = 1; i < n; i++) {
    km_lag[i] = km[i - 1];
  }

  // dF1, dF2 の計算： dF = km_lag * (event_ind / n_atrisk)
  NumericVector dF1(n);
  NumericVector dF2(n);
  for (int i = 0; i < n; i++) {
    dF1[i] = km_lag[i] * (epsilon1[i] / n_atrisk[i]);
    dF2[i] = km_lag[i] * (epsilon2[i] / n_atrisk[i]);
  }

  // CIF の累積和計算
  NumericVector CIF1(n);
  NumericVector CIF2(n);
  double sum1 = 0.0, sum2 = 0.0;
  for (int i = 0; i < n; i++) {
    sum1 += dF1[i];
    sum2 += dF2[i];
    CIF1[i] = sum1;
    CIF2[i] = sum2;
  }

  return DataFrame::create(
    Named("t") = sorted_t,
    Named("KM") = km,
    Named("CIF1") = CIF1,
    Named("CIF2") = CIF2
  );
}
