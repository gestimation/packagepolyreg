# author: Honda
# date: 2025/03/11
#library(survival)


calcCIF <- function(
    nuisance.model,
    exposure,
    #strata = NULL,
    data
    #,code.event1 = 1,
    #code.event2 = 2,
    #code.censoring = 0,
    #code.exposure.ref = 0
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
  event_all <- data[[event_var]]

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
  s_atrisk <- km_surv$n.atrisk # Y_l
  s_data <- cbind(s_time, s_surv, s_atrisk)

  # event1のKaplan-Meier推定
  km_event1 <- summary(survfit(Surv(time, event_one) ~ 1))
  event1_time <- km_event1$time # time of d_1l
  n.event1 <- km_event1$n.event # value of d_1l
  event1_data <- cbind(event1_time, n.event1)

  # event2のKaplan-Meier推定
  km_event2 <- summary(survfit(Surv(time, event_two) ~ 1))
  event2_time <- km_event2$time # time of d_2l
  n.event2 <- km_event2$n.event # value of d_2l
  event2_data <- cbind(event2_time, n.event2)

  #print(dim(s_data))
  #print(dim(event1_data))
  #print(dim(event1_data))

  # event1のCIF
}

path = "C:/PROJECT/polyregData/20250210diabetes_dataset.csv"
jdcs = read.csv(path)

model1 <- "Event(t,epsilon) ~ sex"
model1 <- as.formula(model1)
model2 <- "Event(t,epsilon) ~ age+sex+bmi+hba1c+drug_oha+drug_insulin"
model2 <- as.formula(model2)
model3 <- "Event(t,epsilon) ~ age+sex+bmi+hba1c+diabetes_duration+drug_oha+drug_insulin+sbp+ldl+hdl+tg+current_smoker+alcohol_drinker+ltpa"
model3 <- as.formula(model3)


calcCIF(nuisance.model = model3, exposure = 'fruitq1',data = jdcs)

