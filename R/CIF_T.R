
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

