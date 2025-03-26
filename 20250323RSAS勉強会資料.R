# SAS/R勉強会の参加者にやっておいてほしいこと

###########################################################################################
# GitHubの個人アカウント作成（必須）
# 詳細は以下のURL参照
# https://docs.github.com/ja/get-started/start-your-journey/creating-an-account-on-github

###########################################################################################
# RiffomonasプロジェクトのYouTube動画をみる（時間があれば）
# https://www.youtube.com/watch?v=XjolVT16YNw

###########################################################################################
# Rの行列と関数の扱いに慣れる(重要)

# 例題1. 行列

a <- 1:10
b <- 1:10
print(a)          # 1から10までの値を持つベクトルができた
print(a*b)        # aとbを掛け算すると1から100までの値を持つベクトルに
print(length(a))  # aの長さは10

A <- as.matrix(a) # ベクトルaを行列Aに変換した
B <- as.matrix(b)
print(dim(A))     # 関数dimで次元を確認すると10×1行列

print(A*b)        # 行列Aとベクトルbの「要素ごとの積」
print(A*B)        # 行列Aと行列Bの「要素ごとの積」

print(A %*% B)    # 行列の積は%*%で計算する（右の式は行列の次元があわないから計算できない）
print(t(A) %*% B) # 行列Aの転置行列と行列Bの積


# 例題2. ベクトル・行列の特殊計算

print(cumsum(a))  # 関数cumsumでベクトルaの累積の和を計算
print(cumprod(a)) # 関数cumprodでベクトルaの累積の積を計算
?cumsum           # 知らない関数がでてきたら

C <- matrix(1:9, nrow = 3, byrow = TRUE)
print(C)
print(C[2,2])     # 行列の要素
print(sum(C))     # 行列の和
print(rowSums(C)) # 行列の行の和
print(colSums(C)) # 行列の列の和


# 例題3. 与えられた生存時間データからアットリスク数を計算する関数を自作（実行して結果を確かめよ）

calculateNumberAtRisk <- function(time) {
  atrisk <- outer(time, time, "<=")          # なぜこの関数でアットリスク数が得られるかは時間があるときに考えよ
  n_atrisk <- rowSums(atrisk)
  return(n_atrisk)
}

lambda <- 1                                  # ハザードを設定
N <- 10
survival_time <- rexp(N, rate = lambda)      # 指数分布から生存時間を発生
print(survival_time)
n_atrisk <- calculateNumberAtRisk(survival_time)
print(n_atrisk)                              # アットリスク数が正しいか確認すること


# 例題4. これらを組み合わせると生存時間データからのKaplan-Meier推定量を計算できる（実行して結果を確かめよ）

calculateKaplanMeier <- function(t, d){
  n <- length(t)
  data <- data.frame(t = t, d = d, id = 1:n)
  sorted_data <- data[order(data$t), ]
  sorted_t <- sorted_data$t
  sorted_d <- sorted_data$d
  sorted_id <- sorted_data$id
  atrisk <- outer(sorted_t, sorted_t, "<=") # なぜこの関数でアットリスク数が得られるかは時間があるときに考えよ
  n_atrisk <- rowSums(atrisk)
  s <- 1 - sorted_d / n_atrisk
  km <- cumprod(s)
  data <- data.frame(id=sorted_id, km=km)
  sorted_data <- data[order(data$id), ]
  km <- sorted_data$km
  return(km)
}

survival_time <- rexp(N, rate = lambda)     # 指数分布から生存時間を発生
censoring_time <- rexp(N, rate = 0.5)       # 指数分布から打ち切り時間を発生
observed_time <- min(survival_time,censoring_time)
event <- as.numeric(survival_time<=censoring_time)

km <- calculateKaplanMeier(observed_time, event)
print(km)


###########################################################################################
# 補足. データフレーム
# データフレームは行列と異なり, 数値とテキストを混在することができ, 各行・各列がラベルを持っている

example <- data.frame(age = c(45, 32, 20, 25), sex = c("M", "F", "F", "M"), id = 1:4)
print(example)
print(as.matrix(example))
print(as.numeric(as.matrix(example)))

# データフレームXの列iはX[i]で参照できるが, ラベル（変数名）xということがわかればX[["x"]]やX$xと書くこともできる
print(example$age)
print(example[1])
print(example[["age"]])

# データフレームの行数や列数はnrow関数やncol関数で取得できる
print(nrow(example))
print(ncol(example))

# データフレームの行名や列名はrownames関数やcolnames関数で取得できる
print(rownames(example))
print(colnames(example))

# データフレームの行名や列名を変更するときはrownames関数やcolnames関数を使って変更する
colnames(example) <- c("age_in_year", "sex_in_MF", "id")
print(example)

# order(): 数値またはテキストから順位情報を読み取る関数
print(order(example$age))
print(order(example$sex))

# ある要素の順位情報に従ってデータフレームの順番を変えたいときは, order関数またはsort関数が利用できる
sorted_data <- example[order(example$age), ]
print(sorted_data)


###########################################################################################
# CIF_H（本田さんの関数）
# packagepolyregを読み込むことでCIF_Hが使用できる

CIF_H(nuisance.model = Event(t,epsilon) ~ +1, exposure = 'fruitq1', data = diabetes.complications)


###########################################################################################
# cif（パッケージmetsの関数）

out <- cif(Event(t,epsilon) ~ +1, data=diabetes.complications, cause=1)
par(mfrow=c(1,2))
bplot(out,se=TRUE)


###########################################################################################
# 関数について情報を得るにはgetAnywhere関数と?を用いる
getAnywhere(cif)
?cif

getAnywhere(CIF_H)
?CIF_H


###########################################################################################
# テスト用のデータセットを作成

n_test <- 10
t_test <- 1:n_test
strata_test <- (t_test>5)
epsilon_test <- rep(1, n_test)
epsilon_test[9] <- 0
epsilon_test[10] <- 2
df_test <- data.frame(id = 1:n_test, t_test = t_test, epsilon_test = epsilon_test, strata_test = strata_test)
print(df_test)


###########################################################################################
# CIF_Hのテスト結果

CIF_H(nuisance.model = Event(t_test,epsilon_test) ~ +1, exposure = "t_test", data = df_test)


###########################################################################################
# cifのテスト結果
# CIF_H内部の計算結果を調べなければならないことがわかる

out <- cif(Event(t_test,epsilon_test) ~ +1, data=df_test, cause=1)
#par(mfrow=c(1,2))
#bplot(out,se=TRUE)
print(cbind(out$mu, out$times))


###########################################################################################
# CIF_Hの補足1. フォーミュラ型
# フォーミュラ型は, 「~」の左側と右側に変数名を記述するR特有の表現
# 多くの場合, 左辺はアウトカム, 右辺は共変量
# データフレームのどの変数を用いるかを指定するときに有用

example1 <- y ~ x1+x2
class(example1)
print(example1[1])
print(example1[2])
print(example1[3])

# CIF_Hでは, フォーミュラ（nuisance.model）から変数名を読み取っている
# Surv(, )またはEvent(, ): 生存時間アウトカムを表すための関数（パッケージsurvivalまたはmetsに依存するので注意）
# terms(): フォーミュラから変数名を読み取るための関数
# all.vars(): terms()の出力から変数名を読み取るための関数
# ~+1: 定数項

nuisance.model <- Event(t_test, epsilon_test) ~ +1
data <- df_test

lhs <- terms(nuisance.model)
time_var <- all.vars(lhs)[1] # 時間の変数名
event_var <- all.vars(lhs)[2] # イベントの変数名
time <- data[[time_var]] # 時間の情報
event_all_ <- data[[event_var]] # イベントコードの情報
print(time)
print(event_all_)

# 別の読み取り方もある
# 以下のコードではフォーミュラ（nuisance.model）とデータフレーム（df_test）を
# モデルフレーム型（example2）に変換し, データを取得している
# モデルフレームはデータそのもので, モデルフレームを扱うために用意された関数が利用できることも特長
# 一方でフォーミュラは変数名の情報にすぎず, 実在するデータと結びついていない
# model.frame(): フォーミュラとデータフレームを結びつけるための関数
# model.extract(): モデルフレームからデータを読み取るための関数

nuisance.model <- Event(t_test, epsilon_test) ~ +1
special <- NULL
out_terms <- terms(nuisance.model, special, data = df_test)
example2 <- match.call()
example2$formula <- out_terms
example2[[1]] <- as.name("model.frame")
example2 <- eval(example2, parent.frame())
y <- model.extract(example2, "response")
time <- y[, 1] # 時間の情報
event_all_ <- y[, 2] # イベントコードの情報

class(example2)
print(y)
print(time)
print(event_all_)

# model.matrix(): モデルフレームから回帰モデルのデザイン行列を作るための関数
# デザイン行列が特異行列にならないように自動的にダミー変数を作成してくれるが
# デザイン行列のコーディングによって, 回帰係数の解釈が異なるため, 仕様を理解する必要がある
# コーディングは, 対比（contrasts.arg）によって指定できる
# デフォルト: contr.treatment（切片項は1, ダミー変数は0と1でコーディング）

stratified.model <- Event(t_test, epsilon_test) ~ strata_test
special <- NULL
out_terms <- terms(stratified.model, special, data = df_test)
example3 <- match.call()
example3$formula <- out_terms
example3[[1]] <- as.name("model.frame")
example3 <- eval(example3, parent.frame())
design_matrix_t <- model.matrix(out_terms, example3, contrasts.arg = list(strata_test = "contr.treatment"))
design_matrix_s <- model.matrix(out_terms, example3, contrasts.arg = list(strata_test = "contr.sum"))
design_matrix_p <- model.matrix(out_terms, example3, contrasts.arg = list(strata_test = "contr.poly"))
design_matrix_S <- model.matrix(out_terms, example3, contrasts.arg = list(strata_test = "contr.SAS"))
design_matrix_h <- model.matrix(out_terms, example3, contrasts.arg = list(strata_test = "contr.helmert"))
print(design_matrix_S)

# フォーミュラから変数xを減らすときは, マイナスを用いて-xと指定する
# 以下のフォーミュラは切片項を含まないデザイン行列に対応する
# このデザイン行列を用いて線型モデルを当てはめると, 共変量の水準ごとの平均が推定される

no.intercept.model <- Event(t_test, epsilon_test) ~ -1+strata_test
special <- NULL
out_terms <- terms(no.intercept.model, special, data = df_test)
example4 <- match.call()
example4$formula <- out_terms
example4[[1]] <- as.name("model.frame")
example4 <- eval(example4, parent.frame())
design_matrix_n <- model.matrix(out_terms, example4)
print(design_matrix_n)


# CIF_Hの補足2. 外部関数の利用
# ifelse(): TRUEまたはFALSEを判定し, 判定に応じて別の値を返す
# survfit(): Kaplan-Meier推定量を計算（パッケージsurvivalに依存するので注意）

event_surv <- ifelse(event_all_ == 2, 1, event_all_)
km_surv <- summary(survfit(Surv(time, event_surv) ~ 1))
print(km_surv)

# CIF_Hの補足3. インデックスを用いた行列演算
# Rで行われる計算のほとんどは行列演算であり
# プログラミング上も, ベクトルや行列の線型演算として書いた方が, 可読性・速度が高い
# 一方で, べき乗や微分積分など線型ではない計算をプログラミングするときは
# forループなどを用いる必要がある
# ベクトルXの要素iはX[i]により参照できる
# 以下のように, iについてループを回すと, ベクトルの要素ごとの計算ができる

n_1 <- length(t_test)
cum <- 0
for(i in 1:n_1){
  cum <- cum + t_test[i]
}
print(1:n_1)
print(cum)

# forループを書くときには, 最低限でもループの最初・最後について, 挙動を確認しておくべきである
# 上のような簡単なアルゴリズムですら, たとえばcumが定義されていないと, 最初の計算がうまくいかない

n_1 <- length(t_test)
cum <- NULL
for(i in 1:n_1){
  cum <- cum + t_test[i]
}
print(1:n_1)
print(cum)

# 最後のループも要確認である. 以下の例では, インデックスの終点がベクトルの範囲を超えてしまっている

n_1 <- length(t_test)+1
cum <- 0
for(i in 1:n_1){
  cum <- cum + t_test[i]
}
print(1:n_1)
print(cum)


###########################################################################################
# CIF_H内部の計算結果

CIF1_value <- calculateCIF1(t_test, epsilon_test)
print(CIF1_value)

CIF1_value <- calculateCIF2(t_test, epsilon_test)
print(CIF1_value)

###########################################################################################
# CIF_HとsurvfitをcalculateKaplanMeierに置き換えた結果

s_time <- t_test
s_surv <- calculateKaplanMeier(t_test, as.numeric(epsilon_test>0))
s_atrisk <- rowSums(createAtRiskMatrix(t_test))
n.event1 <- as.numeric(epsilon_test==1)
event1_time <- t_test
CIF1_value <- calculateCIF3(event1_time, s_time, n.event1, s_surv, s_atrisk)
print(CIF1_value)

###########################################################################################
# 計算スピードを測る
# tic(): 計算のスタートを決める関数
# toc(): 計算のゴールを決める関数

library(tictoc)
n_test <- 10000
t_test <- 1:n_test
epsilon_test <- rep(1, n_test)
epsilon_test[1] <- 0
epsilon_test[2] <- 2
epsilon_test[10] <- 2
df_test <- data.frame(id = 1:n_test, t_test = t_test, epsilon_test = epsilon_test)

tic()
CIF1_value <- calculateCIF2(t_test, epsilon_test)
toc()

tic()
s_time <- t_test
s_surv <- calculateKaplanMeier(t_test, as.numeric(epsilon_test>0))
s_atrisk <- rowSums(createAtRiskMatrix(t_test))
n.event1 <- as.numeric(epsilon_test==1)
event1_time <- t_test
CIF1_value <- calculateCIF3(event1_time, s_time, n.event1, s_surv, s_atrisk)
toc()


###########################################################################################
# calculateCIF2はKaplan-Meier推定量をsurvfit関数（内部でCを使用）で計算しているため速い
# calculateKaplanMeierはRcpp関数を用いて高速化できる

cppFunction('
NumericVector calculateKaplanMeier_new(NumericVector t, NumericVector d) {
  int n = t.size();
  NumericVector km(n);
  IntegerVector indices = seq(0, n-1);
  std::sort(indices.begin(), indices.end(), [&](int i, int j) {
    return t[i] < t[j];
  });

  NumericVector sorted_t(n);
  NumericVector sorted_d(n);
  IntegerVector sorted_id(n);
  for (int i = 0; i < n; i++) {
    sorted_t[i] = t[indices[i]];
    sorted_d[i] = d[indices[i]];
    sorted_id[i] = indices[i] + 1;
  }

  int n_atrisk = n;
  double km_i = 1.0;
  for (int i = 0; i < n; i++) {
    if (i > 0 && sorted_t[i] != sorted_t[i-1]) {
      n_atrisk = n - i;
    }
    km_i *= (1 - sorted_d[i] / (double)n_atrisk);
    km[i] = km_i;
  }
  return km;
}
')

library(Rcpp)
library(microbenchmark)
n_test <- 30000
t_test <- 1:n_test
epsilon_test <- rep(1, n_test)
epsilon_test[1] <- 0
epsilon_test[2] <- 2
epsilon_test[10] <- 2
df_test <- data.frame(id = 1:n_test, t_test = t_test, epsilon_test = epsilon_test)
microbenchmark(
  calculateKaplanMeier(t_test, as.numeric(epsilon_test>0)),
  calculateKaplanMeier_new(t_test, as.numeric(epsilon_test>0))
)


