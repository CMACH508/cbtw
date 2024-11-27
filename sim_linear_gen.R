library(MASS)
library(ggplot2)

path <- "./sim_data/"
feature_seq <- c(12, 24)
filenum <- 1000
set_mu <- 3
set_sigma <- 3
ATE <- 10
cov_wt <- 2
intercep <- 5
sample_seq <- c(100, 150, 200, 250, 300)
noise_sigma <- 1

for (feature_num in feature_seq) {
  meshX <- matrix(rep(1 : feature_num, feature_num),feature_num)
  meshY <- matrix(rep(1 : feature_num, each = feature_num),feature_num)
  data_mu <- rep(set_mu, feature_num)
  data_sigma <- set_sigma * 0.5 ^ abs(meshX - meshY)
  wt <- c(ATE, (1 : feature_num) / cov_wt, intercep)
  logis_wt <- rep(c(1, -1), feature_num / 2)
  for (sample_num in sample_seq) {
    filedir = paste0(path, feature_num, "/nocenter", "/samplenum_", sample_num, "/")
    dir.create(filedir, showWarnings = F, recursive = T)
    for (k in 1 : filenum) {
      data_cov <- mvrnorm(sample_num,  data_mu, data_sigma)
      log_exp <- exp(data_cov %*% logis_wt + rnorm(sample_num))
      post_p <- 1 / (1 + log_exp)
      is_treat <- rbinom(post_p, 1, post_p)
      outcome <- cbind(is_treat, data_cov, 1) %*% wt + rnorm(sample_num, sd = noise_sigma)
      is_treat[is_treat == 1] <- "no"
      is_treat[is_treat == 0] <- "yes"
      data <- data.frame(outcome, is_treat, data_cov)
      write.csv(data, file = paste0(filedir, k, ".csv"), row.names =  F)
    }
  }
}
