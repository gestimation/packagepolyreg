# author: Honda
# date: 2025/03/11
#library(survival)


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


# test code below
path = "C:/PROJECT/polyregData/20250210diabetes_dataset.csv"
jdcs = read.csv(path)

model1 <- "Event(t,epsilon) ~ sex"
model1 <- as.formula(model1)
model2 <- "Event(t,epsilon) ~ age+sex+bmi+hba1c+drug_oha+drug_insulin"
model2 <- as.formula(model2)
model3 <- "Event(t,epsilon) ~ age+sex+bmi+hba1c+diabetes_duration+drug_oha+drug_insulin+sbp+ldl+hdl+tg+current_smoker+alcohol_drinker+ltpa"
model3 <- as.formula(model3)

CIF_H(nuisance.model = model3, exposure = 'fruitq1',data = jdcs)

