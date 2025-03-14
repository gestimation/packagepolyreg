# author: Honda
# date: 2025/03/11
#library(survival)


calcCIF_H <- function(
    nuisance.model,
    exposure,
    #strata = NULL,
    data,
    code.event1 = 1,
    code.event2 = 2,
    code.censoring = 0,
    code.exposure.ref= 0,
    devideByExposure = FALSE,
    drawGraph = TRUE
){
  # 曝露変数
  expo_ <- data[[exposure]]
  #codeから0,1に変更
  expo <- ifelse(expo_ == code.exposure.ref, 0, 1)

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


  if(!devideByExposure){
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
    # 右端を追加するのをやめた
    #CIF1 <- rbind(CIF1, t(as.matrix(c(s_time[length(s_time)], CIF1_value[n_1]))))

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
    # 右端を追加するのをやめた
    #CIF2 <- rbind(CIF2, t(as.matrix(c(s_time[length(s_time)], CIF2_value[n_2]))))

    if(drawGraph){
      # 重ねてプロット
      plot(CIF1[,1], CIF1[,2], type = "s", main = "Cumulative Incidence Function",
         xlab = "Time", ylab = "Cumulative Incidence", col = "blue", lwd = 2, ylim = c(0, max(CIF1[,2], CIF2[,2])))
      lines(CIF2[,1], CIF2[,2], type = "s", col = "red", lwd = 2)
      legend("topleft", legend = c("Event 1", "Event 2"), col = c("blue", "red"), lwd = 2, cex = 0.5)
    }
    return(list(CIF1 = CIF1, CIF2 = CIF2))
  }
  else{
    # 曝露群別にデータ分け
    time_0 <- time[expo == 0]
    time_1 <- time[expo == 1]

    event_all_0 <- event_all[expo == 0]
    event_all_1 <- event_all[expo == 1]

    event_surv_0 <- event_surv[expo == 0]
    event_surv_1 <- event_surv[expo == 1]

    event_one_0 <- event_one[expo == 0]
    event_one_1 <- event_one[expo == 1]

    event_two_0 <- event_two[expo == 0]
    event_two_1 <- event_two[expo == 1]

    # s_hatのKaplan-Meier推定(expo == 0)
    km_surv_0 <- summary(survfit(Surv(time_0, event_surv_0) ~ 1))
    s_time_0 <- km_surv_0$time # time of s_hat
    s_surv_0 <- km_surv_0$surv # value of s_hat
    s_atrisk_0 <- km_surv_0$n.risk # Y_l

    # s_hatのKaplan-Meier推定(expo == 1)
    km_surv_1 <- summary(survfit(Surv(time_1, event_surv_1) ~ 1))
    s_time_1 <- km_surv_1$time # time of s_hat
    s_surv_1 <- km_surv_1$surv # value of s_hat
    s_atrisk_1 <- km_surv_1$n.risk # Y_l

    # event1のKaplan-Meier推定(expo == 0)
    km_event1_0 <- summary(survfit(Surv(time_0, event_one_0) ~ 1))
    event1_time_0 <- km_event1_0$time # time of d_1l
    n.event1_0 <- km_event1_0$n.event # value of d_1l

    # event1のKaplan-Meier推定(expo == 1)
    km_event1_1 <- summary(survfit(Surv(time_1, event_one_1) ~ 1))
    event1_time_1 <- km_event1_1$time # time of d_1l
    n.event1_1 <- km_event1_1$n.event # value of d_1l

    # event2のKaplan-Meier推定(expo == 0)
    km_event2_0 <- summary(survfit(Surv(time_0, event_two_0) ~ 1))
    event2_time_0 <- km_event2_0$time # time of d_2l
    n.event2_0 <- km_event2_0$n.event # value of d_2l

    # event2のKaplan-Meier推定(expo == 1)
    km_event2_1 <- summary(survfit(Surv(time_1, event_two_1) ~ 1))
    event2_time_1 <- km_event2_1$time # time of d_2l
    n.event2_1 <- km_event2_1$n.event # value of d_2l

    # event1のCIF(expo == 0)
    # 横軸はevent1_time
    n_1_0 <- length(event1_time_0)
    CIF1_value_0 <- rep(NA, n_1_0)
    cum <- 0
    for(i in 1:n_1_0){
      index <- sum(as.numeric((s_time_0-event1_time_0[i])<0))
      if(index!=0){
        cum <- cum + n.event1_0[i]*s_surv_0[index]/s_atrisk_0[index]
      }
      CIF1_value_0[i] <- cum
    }
    CIF1_0 <- cbind(event1_time_0, CIF1_value_0)
    # 右端を追加するのをやめた
    #CIF1_0 <- rbind(CIF1_0, t(as.matrix(c(s_time_0[length(s_time_0)], CIF1_value_0[n_1_0]))))

    # event1のCIF(expo == 1)
    # 横軸はevent1_time
    n_1_1 <- length(event1_time_1)
    CIF1_value_1 <- rep(NA, n_1_1)
    cum <- 0
    for(i in 1:n_1_1){
      index <- sum(as.numeric((s_time_1-event1_time_1[i])<0))
      if(index!=0){
        cum <- cum + n.event1_1[i]*s_surv_1[index]/s_atrisk_1[index]
      }
      CIF1_value_1[i] <- cum
    }
    CIF1_1 <- cbind(event1_time_1, CIF1_value_1)
    # 右端を追加するのをやめた
    #CIF1_1 <- rbind(CIF1_1, t(as.matrix(c(s_time_1[length(s_time_1)], CIF1_value_1[n_1_1]))))

    # event2のCIF(expo == 0)
    # 横軸はevent2_time
    n_2_0 <- length(event2_time_0)
    CIF2_value_0 <- rep(NA, n_2_0)
    cum <- 0
    for(i in 1:n_2_0){
      index <- sum(as.numeric((s_time_0-event2_time_0[i])<0))
      if(index!=0){
        cum <- cum + n.event2_0[i]*s_surv_0[index]/s_atrisk_0[index]
      }
      CIF2_value_0[i] <- cum
    }
    CIF2_0 <- cbind(event2_time_0, CIF2_value_0)
    # 右端を追加するのをやめた
    #CIF2_0 <- rbind(CIF2_0, t(as.matrix(c(s_time_0[length(s_time_0)], CIF2_value_0[n_2_0]))))

    # event2のCIF(expo == 1)
    # 横軸はevent2_time
    n_2_1 <- length(event2_time_1)
    CIF2_value_1 <- rep(NA, n_2_1)
    cum <- 0
    for(i in 1:n_2_1){
      index <- sum(as.numeric((s_time_1-event2_time_1[i])<0))
      if(index!=0){
        cum <- cum + n.event2_1[i]*s_surv_1[index]/s_atrisk_1[index]
      }
      CIF2_value_1[i] <- cum
    }
    CIF2_1 <- cbind(event2_time_1, CIF2_value_1)
    # 右端を追加するのをやめた
    #CIF2_1 <- rbind(CIF2_1, t(as.matrix(c(s_time_1[length(s_time_1)], CIF2_value_1[n_2_1]))))

    if(drawGraph){
      ylim_max <- max(CIF1_0[,2], CIF1_1[,2], CIF2_0[,2], CIF2_1[,2])

      plot(CIF1_0[,1], CIF1_0[,2], type = "s", main = "Cumulative Incidence Function",
           xlab = "Time", ylab = "Cumulative Incidence", col = "blue", lwd = 2, lty = 1, ylim = c(0, ylim_max))

      lines(CIF1_1[,1], CIF1_1[,2], type = "s", col = "blue", lwd = 2, lty = 2)
      lines(CIF2_0[,1], CIF2_0[,2], type = "s", col = "red", lwd = 2, lty = 1)
      lines(CIF2_1[,1], CIF2_1[,2], type = "s", col = "red", lwd = 2, lty = 2)

      legend("topleft",
             legend = c("Event 1, expo=0", "Event 1, expo=1", "Event 2, expo=0", "Event 2, expo=1"),
             col = c("blue", "blue", "red", "red"),
             lwd = 2, lty = c(1, 2, 1, 2), cex = 0.5)
    }
  }
  return(list(CIF1_unexposured = CIF1_0, CIF1_exposured = CIF1_1,
              CIF2_unexposured = CIF2_0, CIF2_exposured = CIF2_1))
}


# test code below
path = "C:/PROJECT/polyregData/20250210diabetes_dataset.csv"
jdcs = read.csv(path)

model1 <- "Event(t,epsilon) ~ sex"
model1 <- as.formula(model1)
model2 <- "Event(t,epsilon) ~ age+sex+bmi+hba1c+drug_oha+drug_insulin"
model2 <- as.formula(model2)
model3 <- "Event(t,epsilon) ~ age+sex+bmi+hba1c+diabetes_duration+drug_oha+drug_insulin+sbp+ldl+hdl+tg+current_smoker+alcohol_drinker+ltpa"
model3 <- as.formula(model3)

#CIF_H(nuisance.model = model3, exposure = 'fruitq1',
#      data = jdcs, devideByExposure = TRUE, drawGraph = TRUE)
CIF_H(nuisance.model = model3, exposure = 'fruitq1',
      data = jdcs, devideByExposure = FALSE, drawGraph = TRUE)
