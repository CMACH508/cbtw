source("lbmse.R")
library(parallel)
library(plyr)

path <- "./real_data/lalonde/"
trtname <- "is_treat"
outname <- "re78"
file_num <- 1000
method <- c("overlap", "IPW", "none", "optim")
sample_seq <- c(100, 150, 200, 250, 300)
reserve_cores <- 0

cores <- detectCores() - reserve_cores
for (sample_num in sample_seq) {
  filedir <- paste0(path, sample_num, "/")
  filenames <- head(list.files(filedir), file_num)
  data_list <- mclapply(as.list(paste0(filedir, filenames)),
                        read.csv,
                        stringsAsFactors = T,
                        mc.cores = cores)
  result <- mclapply(data_list,
                     ATE_smd_est,
                     trtname = trtname,
                     outname = outname,
                     method = method,
                     mc.cores = cores)
  list_result <- do.call(Map, c(f = rbind.fill, result))
  write.csv(list_result[[1]], file = paste0(filedir, "LBMSE_smd.csv"), row.names = F)
  write.csv(list_result[[2]], file = paste0(filedir, "LBMSE_ATE.csv"), row.names = F)
}