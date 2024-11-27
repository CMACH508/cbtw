library(parallel)

lalonde_gen <- function(data, treat_data, control_data, sample_num, filedir){
  while (TRUE) {
    subdata_con <- control_data[sample(1: nrow(control_data), sample_num / 2), ]
    subdata_tre <- treat_data[sample(1: nrow(treat_data), sample_num / 2), ]
    subdata <- rbind(subdata_con, subdata_tre)
    ident <- sapply(lapply(subdata, unique), length)
    if (! any(ident == 1)) {
          write.csv(subdata, file = paste0(filedir, as.numeric(Sys.time()), ".csv"), row.names =  F)
          break
    }
  }
}

path <- "./real_data/"
outputpath <- "./real_data/lalonde/"
treat_file <- "la_treat.csv"
control_file <- "la_control.csv"
sample_seq <- c(100, 150, 200, 250, 300)
reps <- 1000
cores <- detectCores()

for (sample_num in sample_seq) {
  filedir <- paste0(outputpath, sample_num, "/")
  treat_data <- read.csv(paste0(path, treat_file), stringsAsFactors = T)
  control_data <- read.csv(paste0(path, control_file), stringsAsFactors = T)
  data_list <- rep(list(0), reps)
  dir.create(filedir, showWarnings = F, recursive = T)
  mclapply(data_list, lalonde_gen, treat_data, control_data, sample_num, filedir, mc.cores = cores)
}

