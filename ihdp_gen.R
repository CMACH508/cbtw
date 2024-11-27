library(parallel)

ihdp_gen <- function(data_list, data, sample_num, filedir){
  while (TRUE) {
    subdata <- data[sample(1: nrow(data), sample_num), ]
    ident <- sapply(lapply(subdata, unique), length)
    if (! any(ident == 1)) {
      write.csv(subdata, file = paste0(filedir, as.numeric(Sys.time()), ".csv"), row.names =  F)
      break
    }
  }
}

path <- "./real_data/"
outputpath <- "./real_data/ihdp/"
file <- "ihdp.csv"
sample_seq <- c(100, 150, 200, 250, 300)
reps <- 1000
cores <- detectCores()

for (sample_num in sample_seq) {
  filedir <- paste0(outputpath, sample_num, "/")
  data <- read.csv(paste0(path, file), stringsAsFactors = T)
  # control_data <- read.csv(paste0(path, control_file), stringsAsFactors = T)
  data_list <- rep(list(0), reps)
  dir.create(filedir, showWarnings = F, recursive = T)
  mclapply(data_list, ihdp_gen, data, sample_num, filedir, mc.cores = cores)
}
