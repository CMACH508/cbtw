source("lbmse.R")
library(parallel)

path <- "./sim_data/"
trtname <- "is_treat"
outname <- "outcome"
file_num <- 1000
feature_seq <- c(24)
method <- c("optim", "overlap", "IPW", "none")
sample_seq <- c(100, 150, 200, 250, 300)
whether_linear <- c("/nolinear", "/nolinear")
reserve_cores <- 0

cores <- detectCores() - reserve_cores
for (feature_num in feature_seq) {
  for (whether_li in whether_linear) {
    for (sample_num in sample_seq) {
      filedir = paste0(path, feature_num, whether_li, "/samplenum_", sample_num, "/")
      data_list <- mclapply(as.list(paste0(filedir, 1 : file_num, ".csv")),
                            read.csv,
                            stringsAsFactors = T,
                            mc.cores = cores)
      result <- mclapply(data_list,
                         ATE_smd_est,
                         trtname = trtname,
                         outname = outname,
                         method = method,
                         mc.cores = cores)
      list_result <- do.call(Map, c(f = rbind, result))
      write.csv(list_result[[1]], file = paste0(filedir, "LBMSE_smd.csv"), row.names = F)
      write.csv(list_result[[2]], file = paste0(filedir, "LBMSE_ATE.csv"), row.names = F)
    }
  }
}
