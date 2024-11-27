library(Rsolnp)
library(leaps)
library(reshape2)
library(ggplot2)

best_reg <- function(data, trtname, outname){
  covname <- names(data)[! (names(data) %in% c(trtname, outname))]
  formula <- paste0(trtname, "~", paste(covname, collapse = "+"))
  covM <- model.matrix(formula(formula), data)
  covM <- data.frame(covM[, - 1, drop = F])
  isbig <- F
  if (length(data) > 51){
    isbig <- T
  }
  reg_exha <- summary(regsubsets(as.formula(paste0(outname, "~.")), 
                                 data, 
                                 nvmax = length(data) - 1,
                                 really.big = isbig))
  mod_exha <- reg_exha$which[which.max(reg_exha$adjr2), ]
  covM <- covM[names(covM) %in% names(mod_exha[mod_exha])]
  data <- cbind(data[, c(outname, trtname)], covM)
  return(data)
}

OWweight <- function(data, trtname, outname){
  covname <- names(data)[! (names(data) %in% c(trtname, outname))]
  formula <- paste0(trtname, "~", paste(covname, collapse = "+"))
  fit <- glm(formula = formula, data = data, family = binomial(link = "logit"))
  e <- fit$fitted.values
  e <- cbind(e, 1 - e)
  z <- as.numeric(data[, trtname])
  w <- rep(0, nrow(data))
  w[z == 1] <- e[z == 1, 1]
  w[z == 2] <- e[z == 2, 2]
  return(list(w, z, formula))
}

IPW <- function(data, trtname, outname){
  covname <- names(data)[! (names(data) %in% c(trtname, outname))]
  formula <- paste0(trtname, "~", paste(covname, collapse = "+"))
  fit <- glm(formula = formula, data = data, family = binomial(link = "logit"))
  e <- fit$fitted.values
  e <- cbind(1 / (1 - e), 1 / e)
  z <- as.numeric(data[, trtname])
  w <- rep(0, nrow(data))
  w[z == 1] <- e[z == 1, 1]
  w[z == 2] <- e[z == 2, 2]
  return(list(w, z, formula))
}

lam_reg <- function(data, trtname, outname, covM){
  reg_list <- summary(lm(paste0(outname, "~."), data))
  coef_beta <- reg_list$coefficients[- 1 , 1]
  coef_beta <- coef_beta[colnames(covM)]
  coef_beta[is.na(coef_beta)] <- 0
  coef_sigma <- reg_list$sigma
  return(list(coef_sigma, coef_beta))
}

objfun <- function(w, covM1, covM2, z, coef_sigma, coef_beta){
  w1 <- w[z == 1]
  w2 <- w[z == 2]
  m1 <- apply(covM1, 2, weighted.mean, w1)
  m2 <- apply(covM2, 2, weighted.mean, w2)
  optmean <- sum((coef_beta * (m1 - m2)) ^ 2)
  optweight <- sum(w1 ^ 2) + sum(w2 ^ 2)
  return(optmean + coef_sigma ^ 2 * optweight / length(coef_beta))
}

eqn <- function(w, covM1, covM2, z, coef_sigma, coef_beta){
  w1 <- w[z == 1]
  w2 <- w[z == 2]
  z1 <- sum(w1)
  z2 <- sum(w2)
  return(c(z1, z2))
}

lbmseweight <- function(data, trtname, outname){
  covname <- names(data)[! (names(data) %in% c(trtname, outname))]
  formula <- paste0(trtname, "~", paste(covname, collapse = "+"))
  covM <- model.matrix(formula(formula), data)
  covM <- covM[, - 1, drop = F]
  opt_par <- lam_reg(data, trtname, outname, covM)
  z <- as.numeric(data[, trtname])
  covM1 <- covM[z == 1, , drop = F]
  covM2 <- covM[z == 2, , drop = F]
  w0 <- rep(0, nrow(data))
  w0[z == 1] <- 1 / sum(z == 1)
  w0[z == 2] <- 1 / sum(z == 2)
  opt <- solnp(pars = w0, 
               fun = objfun, 
               eqfun = eqn, 
               eqB = c(1, 1), 
               LB = rep(0, nrow(data)), 
               covM1 = covM1,
               covM2 = covM2,
               z = z, 
               coef_sigma = opt_par[[1]], 
               coef_beta = opt_par[[2]], 
               control = list(trace = 0))
  w <- opt$pars
  return(list(w, z, formula))
}

crude <- function(data, trtname, outname){
  covname <- names(data)[! (names(data) %in% c(trtname, outname))]
  formula <- paste0(trtname, "~", paste(covname, collapse = "+"))
  z <- as.numeric(data[, trtname])
  w <- rep(1, nrow(data))
  return(list(w, z, formula))
}

stdiff <- function(data, w, z, formula){
  covM <- model.matrix(formula(formula), data)
  covM <- covM[, - 1, drop = F]
  covM1 <- covM[z == 1, , drop = F]
  covM2 <- covM[z == 2, , drop = F]
  w1 <- w[z == 1]
  w2 <- w[z == 2]
  m1 <- apply(covM1, 2, weighted.mean, w1)
  m2 <- apply(covM2, 2, weighted.mean, w2)
  var1 <- apply(covM1, 2, var)
  var2 <- apply(covM2, 2, var)
  md <- abs((m1 - m2) / sqrt((var1 + var2) / 2))
  return(md)
}

ATE_est <- function(data, w, z, outname){
  outcome1 <- data[, outname][z == 1]
  outcome2 <- data[, outname][z == 2]
  w1 <- w[z == 1]
  w2 <- w[z == 2]
  outm1 <- weighted.mean(outcome1, w1)
  outm2 <- weighted.mean(outcome2, w2)
  ATE_md <- outm1 - outm2
  return(ATE_md)
}

draw_smd <- function(smd){
  plot_df <- melt(data = smd,
                  id.vars = "meth", 
                  variable.name = "covname", 
                  value.name = "smd")
  pt <- ggplot(plot_df, aes(x = smd, y = covname)) + 
    theme_bw() +
    geom_point(aes(color = meth, shape = meth)) +
    scale_x_continuous(limits = c(0, max(plot_df$smd, 1)), 
                       breaks = c(0, 0.1, 1)) +
    scale_y_discrete(limits = rev(levels(plot_df$covname))) +
    xlab("Standardized Mean Differences") + 
    ylab("") +
    theme(axis.text.x = element_text(size = 13), 
          axis.text.y = element_text(size = 10),
          axis.title.x = element_text(size = 16), 
          legend.text = element_text(size = 13),
          legend.title = element_blank()) + 
    geom_vline(xintercept = 0.1, linetype = "dotted", lwd = 0.5)
  print(pt)
}

ATE_smd_est <- function(data, trtname, outname, method = "optim", use_exha = F, plot_smd = F){
  if (use_exha){
    data <- best_reg(data, trtname, outname)
  }
  smd <- data.frame()
  ATE <- data.frame()
  for (meth in method) {
    if (meth == "overlap"){
      obj <- OWweight(data, trtname, outname)
    }
    else if (meth == "IPW"){
      obj <- IPW(data, trtname, outname)
    }
    else if (meth == "none"){
      obj <- crude(data, trtname, outname)
    }
    else if (meth == "optim"){
      obj <- lbmseweight(data, trtname, outname)
    }
    else{
      stop("weight option not recognized", "\n")
    }
    md <- stdiff(data, obj[[1]], obj[[2]], obj[[3]])
    ATE_md <- ATE_est(data, obj[[1]], obj[[2]], outname)
    smd <- rbind(smd, cbind(meth, data.frame(t(md))))
    ATE <- rbind(ATE, cbind(meth, data.frame(ATE_md)))
  }
  if (plot_smd) {
    draw_smd(smd)
  }
  return(list(smd, ATE))
}