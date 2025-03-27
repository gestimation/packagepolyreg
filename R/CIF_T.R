
################################################################################
#  CIF_H()におけるCIFを計算するコード
#  n_1 <- length(event1_time)
#  CIF1_value <- rep(NA, n_1)
#  cum <- 0
#  for(i in 1:n_1){
#    index <- sum(as.numeric((s_time-event1_time[i])<0))
#    if(index!=0){
#      cum <- cum + n.event1[i]*s_surv[index]/s_atrisk[index]
#    }
#    CIF1_value[i] <- cum
#  }

calculateCIF1 <- function(X, Y) {
  event_surv <- as.numeric(Y>0)
  km_surv <- summary(survfit(Surv(X, event_surv) ~ 1))
  s_time <- km_surv$time # time of s_hat
  s_surv <- km_surv$surv # value of s_hat
  s_atrisk <- km_surv$n.risk # Y_l
  event_one <- as.numeric(Y==1)
  km_event1 <- summary(survfit(Surv(X, event_one) ~ 1))
  event1_time <- km_event1$time # time of d_1l
  n.event1 <- km_event1$n.event # value of d_1l

  n_1 <- length(event1_time)
  CIF1_value <- rep(NA, n_1)
  cum <- 0
  for(i in 1:n_1){
    index <- sum(as.numeric((s_time-event1_time[i])<0))
    if(index!=0){
      cum <- cum + n.event1[i]*s_surv[index]/s_atrisk[index]
    }
    CIF1_value[i] <- cum
  }
  return(CIF1_value)
}

calculateCIF2 <- function(X, Y) {
  event_surv <- as.numeric(Y>0)
  km_surv <- summary(survfit(Surv(X, event_surv) ~ 1))
  s_time <- km_surv$time # time of s_hat
  s_surv <- km_surv$surv # value of s_hat
  s_atrisk <- km_surv$n.risk # Y_l
  event_one <- as.numeric(Y==1)
  km_event1 <- summary(survfit(Surv(X, event_one) ~ 1))
  event1_time <- km_event1$time # time of d_1l
  n.event1 <- km_event1$n.event # value of d_1l

  n_1 <- length(event1_time)
  CIF1_value <- rep(NA, n_1)
  cum <- 0
  for(i in 1:n_1){
    index <- sum(as.numeric((s_time-event1_time[i])<0))
#    index <- sum(as.numeric((s_time-event1_time[i])<=0))
    if (index==0) {
      cum <- cum + n.event1[i]/s_atrisk[index+1]
    } else {
      cum <- cum + n.event1[i]*s_surv[index]/s_atrisk[index+1]
    }
    CIF1_value[i] <- cum
  }
  return(CIF1_value)
}

calculateCIF3 <- function(event1_time, s_time, n.event1, s_surv, s_atrisk) {
  n_1 <- length(event1_time)
  CIF1_value <- rep(NA, n_1)
  cum <- 0
  for(i in 1:n_1){
    index <- sum(as.numeric((s_time-event1_time[i])<0))
    if (index==0) {
      cum <- cum + n.event1[i]/s_atrisk[index+1]
    } else {
      cum <- cum + n.event1[i]*s_surv[index]/s_atrisk[index+1]
    }
    CIF1_value[i] <- cum
  }
  return(CIF1_value)
}

calculateCIF4 <- function(X, Y, Z) {
  event_surv <- as.numeric(Y>0)
  km_surv <- summary(survfit(Surv(X, event_surv) ~ Z))
  s_time <- km_surv$time # time of s_hat
  s_surv <- km_surv$surv # value of s_hat
  s_atrisk <- km_surv$n.risk # Y_l
  s_strata <- km_surv$strata # Stratified levels of Z
  event_one <- as.numeric(Y==1)
  km_event1 <- summary(survfit(Surv(X, event_one) ~ Z))
  event1_time <- km_event1$time # time of d_1l
  n.event1 <- km_event1$n.event # value of d_1l
  event1_strata <- km_event1$strata # Stratified levels of Z

  CIF1_value <- rep(NA, length(event1_time))
  for (level in levels(s_strata)) { # Loop over the levels of Z
    index_strata_1 <- which(s_strata == level) # Identify elements of s_time, s_surv and s_atrisk in the stratum
    stratified_s_time <- s_time[index_strata_1]
    stratified_s_surv <- s_surv[index_strata_1]
    stratified_s_atrisk <- s_atrisk[index_strata_1]
    index_strata_2 <- which(event1_strata == level) # Identify elements of event1_time in the stratum

    cum <- 0
    for (index in index_strata_2) { # Loop over event1_time within a stratum
      if (event1_strata[index] == level) {
        index_km <- sum(as.numeric((stratified_s_time-event1_time[index])<0))
        if (index_km == 0) {
          cum <- cum + n.event1[index] / stratified_s_atrisk[1]
        } else {
          cum <- cum + n.event1[index] * stratified_s_surv[index_km] / stratified_s_atrisk[index_km + 1]
        }
        CIF1_value[index] <- cum
      }
    }
  }
  return(CIF1_value)
}

createAtRiskMatrix <- function(t) {
  atrisk <- outer(t, t, "<=")
  return(atrisk)
}

calculateKaplanMeier <- function(t, d){
  n <- length(t)
  data <- data.frame(t = t, d = d, id = 1:n)
  sorted_data <- data[order(data$t), ]
  sorted_t <- sorted_data$t
  sorted_d <- sorted_data$d
  sorted_id <- sorted_data$id
  atrisk <- createAtRiskMatrix(sorted_t)
  n_atrisk <- rowSums(atrisk)
  s <- 1 - sorted_d / n_atrisk
  km <- cumprod(s)
  data <- data.frame(id=sorted_id, km=km)
  sorted_data <- data[order(data$id), ]
  km <- sorted_data$km
  return(km)
}


#cppFunction('
#NumericVector calculateKaplanMeier_new(NumericVector t, NumericVector d) {
#  int n = t.size();
#  NumericVector km(n);
#  IntegerVector indices = seq(0, n-1);
#  std::sort(indices.begin(), indices.end(), [&](int i, int j) {
#    return t[i] < t[j];
#  });

#  NumericVector sorted_t(n);
#  NumericVector sorted_d(n);
#  IntegerVector sorted_id(n);
#  for (int i = 0; i < n; i++) {
#    sorted_t[i] = t[indices[i]];
#    sorted_d[i] = d[indices[i]];
#    sorted_id[i] = indices[i] + 1; // 1-based indexing
#  }

#  int n_atrisk = n;
#  double km_i = 1.0;
#  for (int i = 0; i < n; i++) {
#    if (i > 0 && sorted_t[i] != sorted_t[i-1]) {
#      n_atrisk = n - i;
#    }
#    km_i *= (1 - sorted_d[i] / (double)n_atrisk);
#    km[i] = km_i;
#  }
#  return km;
#}
#')

calculateKaplanMeier1 <- function(t, d){
  n <- length(t)
  data <- data.frame(t = t, d = d, id = 1:n)
  sorted_data <- data[order(data$t), ]
  sorted_t <- sorted_data$t
  sorted_d <- sorted_data$d
  sorted_id <- sorted_data$id
  n_atrisk <- sapply(sorted_t, function(s) sum(sorted_t >= s))
  s <- 1 - sorted_d / n_atrisk
  km <- cumprod(s)
  data <- data.frame(id=sorted_id, km=km)
  sorted_data <- data[order(data$id), ]
  km <- sorted_data$km
  return(km)
}


calculateKaplanMeier2 <- function(t, d){
  n <- length(t)
  data <- data.frame(t = t, d = d, id = 1:n)
  sorted_data <- data[order(data$t), ]
  sorted_t <- sorted_data$t
  sorted_d <- sorted_data$d
  sorted_id <- sorted_data$id
  n_atrisk <- sapply(sorted_t, function(s) sum(sorted_t >= s))
  #  atrisk <- createAtRiskMatrix(sorted_t)
  #  n_atrisk <- rowSums(atrisk)
  km <- cumprod(1 - sorted_d / n_atrisk)
  data <- data.frame(id=sorted_id, km=km)
  sorted_data <- data[order(data$id), ]
  km <- sorted_data$km
  return(km)
}


calculateKaplanMeier3 <- function(t, d){
  unique_t <- sort(unique(t))
  n_risk <- sapply(unique_t, function(s) sum(t >= s))
  n_event <- sapply(unique_t, function(s) sum(d[t == s]))
  km <- cumprod(1 - n_event / n_risk)
  return(km)
}

calculateKaplanMeier4 <- function(t, d){
  #unique_t <- sort(unique(t))
  unique_t <- sort(unique(t), method = "quick")
  n_risk <- sapply(unique_t, function(s) sum(t >= s))
  n_event <- sapply(unique_t, function(s) sum(d[t == s]))
  km <- cumprod(1 - n_event / n_risk)
  return(km)
}

calculateKaplanMeier5 <- function(t, d){
  unique_t <- t %>%
    unique() %>%
    sort()
  n_risk <- sapply(unique_t, function(s) sum(t >= s))
  n_event <- sapply(unique_t, function(s) sum(d[t == s]))
  km <- cumprod(1 - n_event / n_risk)
  return(km)
}

calculateKaplanMeier6 <- function(t, d){
  unique_t <- sort_unique(t)
  n_risk <- sapply(unique_t, function(s) sum(t >= s))
  n_event <- sapply(unique_t, function(s) sum(d[t == s]))
  km <- cumprod(1 - n_event / n_risk)
  return(km)
}


CIF_H <- function(
    nuisance.model,
    exposure,
    #strata = NULL,
    data,
    code.event1 = 1,
    code.event2 = 2,
    code.censoring = 0
){
  # 曝露変数
  expo <- data[[exposure]]

  # 左辺
  lhs <- terms(nuisance.model)
  # 時間・イベントの変数名
  time_var <- all.vars(lhs)[1]
  event_var <- all.vars(lhs)[2]
  # 時間とイベントコードの情報
  time <- data[[time_var]]
  event_all_ <- data[[event_var]]
  #codeから0,1,2に変更
  event_all <- ifelse(event_all_ == code.event1, 1,
                      ifelse(event_all_ == code.event2, 2,
                             ifelse(event_all_ == code.censoring, 0, NA)))


  # event1, event2を両方1に（s_hatのため）
  event_surv <- ifelse(event_all == 2, 1, event_all)

  # event1を1, event2を0に（event1のCIFのため）
  event_one <- ifelse(event_all == 1, 1, 0)

  # event2を1, event1を0に（event2のCIFのため）
  event_two <- ifelse(event_all == 2, 1, 0)

  # s_hatのKaplan-Meier推定
  km_surv <- summary(survfit(Surv(time, event_surv) ~ 1))
  s_time <- km_surv$time # time of s_hat
  s_surv <- km_surv$surv # value of s_hat
  s_atrisk <- km_surv$n.risk # Y_l

  # event1のKaplan-Meier推定
  km_event1 <- summary(survfit(Surv(time, event_one) ~ 1))
  event1_time <- km_event1$time # time of d_1l
  n.event1 <- km_event1$n.event # value of d_1l

  # event2のKaplan-Meier推定
  km_event2 <- summary(survfit(Surv(time, event_two) ~ 1))
  event2_time <- km_event2$time # time of d_2l
  n.event2 <- km_event2$n.event # value of d_2l

  # event1のCIF
  # 横軸はevent1_time
  n_1 <- length(event1_time)
  CIF1_value <- rep(NA, n_1)
  cum <- 0
  for(i in 1:n_1){
    index <- sum(as.numeric((s_time-event1_time[i])<0))
    if(index!=0){
      cum <- cum + n.event1[i]*s_surv[index]/s_atrisk[index]
    }
    CIF1_value[i] <- cum
  }
  CIF1 <- cbind(event1_time, CIF1_value)
  # 右端を追加
  CIF1 <- rbind(CIF1, t(as.matrix(c(s_time[length(s_time)], CIF1_value[n_1]))))

  # event2のCIF
  # 横軸はevent2_time
  n_2 <- length(event2_time)
  CIF2_value <- rep(NA, n_2)
  cum <- 0
  for(i in 1:n_2){
    index <- sum(as.numeric((s_time-event2_time[i])<0))
    if(index!=0){
      cum <- cum + n.event2[i]*s_surv[index]/s_atrisk[index]
    }
    CIF2_value[i] <- cum
  }
  CIF2 <- cbind(event2_time, CIF2_value)
  # 右端を追加
  CIF2 <- rbind(CIF2, t(as.matrix(c(s_time[length(s_time)], CIF2_value[n_2]))))

  # 重ねてプロット
  plot(CIF1[,1], CIF1[,2], type = "s", main = "Cumulative Incidence Function",
       xlab = "Time", ylab = "Cumulative Incidence", col = "blue", lwd = 2, ylim = c(0, max(CIF1[,2], CIF2[,2])))
  lines(CIF2[,1], CIF2[,2], type = "s", col = "red", lwd = 2)
  legend("topleft", legend = c("Event 1", "Event 2"), col = c("blue", "red"), lwd = 2, cex = 0.6)
}
