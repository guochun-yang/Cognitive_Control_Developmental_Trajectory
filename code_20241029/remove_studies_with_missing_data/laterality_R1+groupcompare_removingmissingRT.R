rm(list=ls())

basedir = 'C:/Users/Guochun Yang/OneDrive - University of Iowa/WithLi/Paper/code' ######
scriptdir = paste(basedir,'code_20241029',sep = '/') #####
funcdir = paste(scriptdir,"func",sep = '/')
funcfiles <- list.files(funcdir, pattern = "\\.R$", full.names = TRUE)
lapply(funcfiles, source)

library(metafor)
library(stringr)
library(rms)
library(mgcv)
library(MuMIn)
library(visreg)
library(readxl)
library(fastDummies)

## load data and calculate the averaged gs and variances

datadir = shortPathName(paste(basedir, "/meta_data/laterality/",sep = ''))
savedir = datadir

if (!dir.exists(savedir)) {
  dir.create(savedir, recursive = TRUE)
}

filename = "laterality_g_all.csv"

AICs = matrix(, nrow = 3, ncol = 4)
BICs = matrix(, nrow = 3, ncol = 4)
QMps = matrix(, nrow = 3, ncol = 4)
QMs = matrix(, nrow = 3, ncol = 4)
modelweights = matrix(, nrow = 3, ncol = 4)

bs<- reses_lin <- reses_qua <- reses_cub <- reses_log <- reses_sqrt <- c()

data = read.table(
  paste(datadir, filename, sep = ''),
  header = TRUE,
  sep = ',',
  dec = '.',
  na.strings = 'NA'
)

# calculate variance based on effect size and sample size
data$laterality_abs = abs(data$laterality_pmn)
data$laterality_nmp = -1 * data$laterality_pmn
data$variance = 1/data$n
data$weight = data$n #1/data$variance

data$LR = sign(data$laterality_pmn)
data$LR12 = (data$LR+3)/2
data$age = as.numeric(data$age)
agefull = max(data$age)+min(data$age)


## load SDM table
sdmdir = paste(basedir, "meta_data",sep = '/')
sdmtable = read.table(
  paste(sdmdir, '/sdm_table.txt', sep = ''),
  header = TRUE,
  sep = '',
  dec = '.'
)


## load the covariates
supptable = read_excel(
  paste(scriptdir, 'covariates.xlsx', sep = '/'),
  sheet = 'covariates'
)

# add corresponding covariates
idx2s<-c()
for (istudy in 1:length(data$study)) {
  # print(istudy)
  idx <- which(supptable$AuthoYear == data$study[istudy])
  data$TaskCode[istudy] <- supptable$TaskCode[idx]
  data$HandnessCode[istudy] <- supptable$HandnessCode[idx]
  data$ContrastCode[istudy] <- supptable$ContrastCode[idx]
  data$ErrorTrialCode[istudy] <- supptable$ErrorTrialCode[idx]
  data$DesignCode[istudy] <- supptable$DesignCode[idx]
  data$SRC[istudy] <- supptable$SRC0[idx]
  data$SRCcode[istudy] <- supptable$SRCcode[idx]
  
  #add age from sdmtable
  idx2 <- which(sdmtable$study == data$study[istudy])
  idx2s<-c(idx2s,idx2)
  data$age[istudy] <- sdmtable$avgAge[idx2]
  data$age0[istudy] <- sdmtable$avgAge[idx2]
}


for (i in 1:3) { #1:absolute value; 2: raw R-L value ; 3:raw L-R #Note: Only the 3 is reported
  if (i == 1) {
    data$laterality = data$laterality_abs
    name = 'absolute_'
  } else if (i == 2) {
    data$laterality = data$laterality_pmn
    name = 'raw_R-L'
  } else {
    data$laterality = data$laterality_nmp
    name = 'raw_L-R'
  }
  
  # delete outliers
  outup = mean(data$laterality) + 3*sd(data$laterality)
  outdown = mean(data$laterality) - 3*sd(data$laterality)
  if (length(which(data$laterality>outup | data$laterality < outdown)) > 0) {
    data = data[-which(data$laterality>outup | data$laterality < outdown),]
  }

  ## add dummy columns
  data <- dummy_cols(data, 
                     select_columns = c("TaskCode","HandnessCode","ContrastCode","ErrorTrialCode","DesignCode"))
  
  ## z-score these columns
  idx_col <- which(colnames(data) %in% c('TaskCode_1','TaskCode_2','TaskCode_3','TaskCode_4',
                                         'HandnessCode_1','HandnessCode_2','HandnessCode_3',
                                         'ContrastCode_1','ContrastCode_2','ContrastCode_3',
                                         'ErrorTrialCode_1','ErrorTrialCode_2','ErrorTrialCode_3',
                                         'DesignCode_1','DesignCode_0','SRCcode'))
  for (icol in 1:length(idx_col)) {
    data[,idx_col[icol]] = as.vector(scale(data[,idx_col[icol]]))
  }
  
  # filter rows
  data <- subset(data,!is.na(SRC))

  ### do two t-test across groups ###
  # between youth and young adults
  data_group12 <- data[which(data$age < 60),]
  data_group23 <- data[which(data$age >= 18),]
  
  for (istudy in 1:length(data_group12$study)) {
    if (data_group12$age[istudy] < 18) {
      data_group12$group[istudy] = 1 
    } else {
      data_group12$group[istudy] = 0  
    }
  }
  for (istudy in 1:length(data_group23$study)) {
    if (data_group23$age[istudy] < 60) {
      data_group23$group[istudy] = 0 
    } else {
      data_group23$group[istudy] = 1  
    }
  }
  res.lin.g12 <- 
    rma(laterality, 
        variance, 
        mods = ~ group + TaskCode_1 + TaskCode_2 + TaskCode_3 + HandnessCode_1 + HandnessCode_2 + ContrastCode_1 + ContrastCode_2 + ErrorTrialCode_1 + DesignCode_1 + SRC,
        data = data_group12)
  p_Rtail_g12 <- 1-pnorm(res.lin.g12$zval[2])
  b12_ub = res.lin.g12$beta[2] + 1.645*res.lin.g12$se[2]
  b12_lb = res.lin.g12$beta[2] - 1.645*res.lin.g12$se[2]
  res.lin.g23 <- 
    rma(laterality, 
        variance, 
        mods = ~ group + TaskCode_1 + TaskCode_2 + TaskCode_3 + HandnessCode_1 + HandnessCode_2 + ContrastCode_1 + ContrastCode_2 + ErrorTrialCode_1 + DesignCode_1 + SRC,
        data = data_group23)  
  p_Rtail_g23 <- 1-pnorm(res.lin.g23$zval[2])
  b23_ub = res.lin.g23$beta[2] + 1.645*res.lin.g23$se[2]
  b23_lb = res.lin.g23$beta[2] - 1.645*res.lin.g23$se[2]
  
  # models
  b <- gam(laterality~s(age)+ TaskCode_1 + TaskCode_2 + TaskCode_3 + HandnessCode_1 + HandnessCode_2 + ContrastCode_1 + ContrastCode_2 + ErrorTrialCode_1 + DesignCode_1 + SRC,data=data, weights = weight, method = 'REML')
  bs[[i]] = b
  res.lin <- 
    rma(laterality, 
        variance, 
        mods = ~ age
        + TaskCode_1 + TaskCode_2 + TaskCode_3 + HandnessCode_1 + HandnessCode_2 + ContrastCode_1 + ContrastCode_2 + ErrorTrialCode_1 + DesignCode_1 + SRC,
        data = data)
  res.qua <- 
    rma(laterality, 
        variance, 
        mods = ~ age + I(age^2)+ TaskCode_1 + TaskCode_2 + TaskCode_3 + HandnessCode_1 + HandnessCode_2 + ContrastCode_1 + ContrastCode_2 + ErrorTrialCode_1 + DesignCode_1 + SRC,
        data = data)
  res.cub <-
    rma(laterality,
        variance,
        mods = ~ age + I(age^2) + I(age^3)+ TaskCode_1 + TaskCode_2 + TaskCode_3 + HandnessCode_1 + HandnessCode_2 + ContrastCode_1 + ContrastCode_2 + ErrorTrialCode_1 + DesignCode_1 + SRC,
        data = data)
  res.log <-
    rma(laterality,
        variance,
        mods = ~ I(log(age)) + I((log(age))^2)+ TaskCode_1 + TaskCode_2 + TaskCode_3 + HandnessCode_1 + HandnessCode_2 + ContrastCode_1 + ContrastCode_2 + ErrorTrialCode_1 + DesignCode_1 + SRC,
        data = data)

  res.sqrt <-
    rma(laterality,
        variance,
        mods = ~ I(sqrt(age)) + age+ TaskCode_1 + TaskCode_2 + TaskCode_3 + HandnessCode_1 + HandnessCode_2 + ContrastCode_1 + ContrastCode_2 + ErrorTrialCode_1 + DesignCode_1 + SRC, 
        data = data)

  reses_qua[[i]] = res.qua
  reses_cub[[i]] = res.cub
  reses_log[[i]] = res.log
  reses_sqrt[[i]] = res.sqrt
  
  AICs[i, 1] = AIC(res.qua)
  AICs[i, 2] = AIC(res.cub)
  AICs[i, 3] = AIC(res.log)
  AICs[i, 4] = AIC(res.sqrt)
  
  BICs[i, 1] = BIC(res.qua)
  BICs[i, 2] = BIC(res.cub)
  BICs[i, 3] = BIC(res.log)
  BICs[i, 4] = BIC(res.sqrt)
  
  QMps[i, 1] = res.qua$QMp
  QMps[i, 2] = res.cub$QMp
  QMps[i, 3] = res.log$QMp
  QMps[i, 4] = res.sqrt$QMp
  
  QMs[i, 1] = res.qua$QM
  QMs[i, 2] = res.cub$QM
  QMs[i, 3] = res.log$QM
  QMs[i, 4] = res.sqrt$QM
  
  mod.comparison <- AIC(res.qua, res.cub, res.log, res.sqrt)
  modelweights[i,] <- Weights(mod.comparison)
}

minAIC <- minBIC <- c()
for (i in 1:3) {
  minAIC[i] = which.min(AICs[i,])
  minBIC[i] = which.min(AICs[i,])
}

library(MASS)

modelweights2 <- format(round(modelweights, 3), nsmall = 3)
  

### draw plot
xs <- seq(8, 75, length=500)

for (i in 3) {
  if (i == 1) {
    name = 'absolute_'
    ylabname = 'laterality(absolute value)'
  } else if (i == 2) {
    name = 'raw_R-L_'
    ylabname = 'laterality(R-L)'
  } else {
    name = 'raw_L-R_'
    ylabname = 'laterality(L-R)'  
  }
  
  mtx_cov_xs <- matrix(0,nrow = length(xs), ncol = length(reses_qua[[i]]$beta)-3)
  if (minAIC[i] == 1) {
    modelbest = reses_qua[[i]]
    title = "QUADRATIC"
    sav <-
      predict(modelbest, newmods = unname(cbind(poly(
        xs, degree = 2, raw = TRUE
      ),mtx_cov_xs)))
    peak_rma = -0.5*modelbest$b[2,1]/modelbest$b[3,1]
  } else if (minAIC[i] == 2) {
    modelbest = reses_cub[[i]]
    title = "CUBIC"
    sav <-
      predict(modelbest, newmods = unname(cbind(poly(
        log(xs), degree = 2, raw = TRUE
      ),mtx_cov_xs)))
    if (modelbest$b[4,1] > 0) {
      peak_rma = (-2*modelbest$b[3,1] - sqrt(4*(modelbest$b[3,1])^2 - 12*modelbest$b[2,1]*modelbest$b[4,1]))/(6*modelbest$b[4,1])
    } else {
      peak_rma = (-2*modelbest$b[3,1] + sqrt(4*(modelbest$b[3,1])^2 - 12*modelbest$b[2,1]*modelbest$b[4,1]))/(6*modelbest$b[4,1])
    }
  } else if (minAIC[i] == 3) {
    modelbest = reses_log[[i]]
    title = "LOG"
    sav <-
      predict(modelbest, newmods = unname(cbind(poly(
        log(xs), degree = 2, raw = TRUE
      ),mtx_cov_xs)))
    peak_rma = exp(-modelbest$b[2,1]/(2*modelbest$b[3,1]))
  } else if (minAIC[i] == 4) {
    modelbest = reses_sqrt[[i]]
    title = "SQRT"
    sav <-
      predict(modelbest, newmods = unname(cbind(poly(
        sqrt(xs), degree = 2, raw = TRUE
      ),mtx_cov_xs)))
    peak_rma = 0.25*(modelbest$b[2,1]/modelbest$b[3,1])^2
  }

  filename = 'log'
    
  
  ylim0 <- c()
  maxvar = max(data$variance)
  ylim0[1] = min(floor((data$laterality-maxvar)*100)/100)
  ylim0[2] = max(ceiling((data$laterality+maxvar)*100)/100)

  modpred <- data.frame(xs,sav$pred,sav$se)
  modpred$upper <- modpred$sav.pred + modpred$sav.se
  modpred$lower <- modpred$sav.pred - modpred$sav.se
  
  
  ### gam prediction
  b <- bs[[i]]
  xs0 <- b$model$age
  preds0 <- predict.gam(b,se.fit = TRUE)
  
  preds <- predict(b, newdata = data.frame(
    age = xs,
    TaskCode_1 = 0,
    TaskCode_2 = 0,
    TaskCode_3 = 0,
    HandnessCode_1 = 0,
    HandnessCode_2 = 0,
    ContrastCode_1 = 0,
    ContrastCode_2 = 0,
    ErrorTrialCode_1 = 0,
    ErrorTrialCode_2 = 0,
    DesignCode_1 = 0,
    SRC = 0
  ),se.fit = TRUE)
  
  peak_gam <- xs[which.min(preds$fit)]
  if (peak_gam > 60) {
    peak_gam <- NaN
  }

  gampred <- data.frame(xs,preds$fit,preds$se.fit)
  gampred$upper <- gampred$preds.fit + gampred$preds.se.fit
  gampred$lower <- gampred$preds.fit - gampred$preds.se.fit

  
  ### plot all at once
  # tiff(
  #   file = paste(savedir, name, title, "_nopeak_rmmissing.tiff", sep = ''),
  #   width = 6,
  #   height = 5,
  #   units = "in",
  #   res = 300
  # )

  # windows()
  par(cex.lab=2)
  par(cex.axis=1.8)
  par(mai=c(1,1.2,0.2,0.2))
  par(mgp=c(3.5,1,0))
  plt = regplot(
    reses_qua[[i]], # this is only to make sure all dots are with the original scale (we can also use reses_cub but not log)
    mod = 2,
    lcol = 'red',
    # pred = FALSE,
    ci=FALSE,
    pred = FALSE,
    xvals = xs,
    las = 1,
    digits = 1,
    bty = "l",
    psize = .20 / sqrt(modelbest$vi),
    xlab = "Age (years)",ylab = "Laterality(L-R)",
    xlim = c(8, 75),
    ylim = c(ylim0[1],ylim0[2])
  )
  lim0 <- par("usr")

  par(new=TRUE)
  # windows()
  lines(modpred$xs, modpred$sav.pred, type = 'l', lwd = 3, col = 'red')  # or type = 'n' if you just want to set up axes
  transparent_red <- rgb(1, 0, 0, alpha = 0.2)  # For example, gray with 50% transparency
  polygon(c(modpred$xs, rev(modpred$xs)), c(modpred$upper, rev(modpred$lower)), col = transparent_red, border = NA)
  lines(modpred$xs, modpred$sav.pred)  # Add main line after to ensure it's on top


  par(new=TRUE)
  transparent_color <- rgb(0, 0, 1, alpha = 0.2)  # For example, gray with 50% transparency
  lines(gampred$xs, gampred$preds.fit, type = 'l', lwd = 3, col = 'blue')  # or type = 'n' if you just want to set up axes
  polygon(c(gampred$xs, rev(gampred$xs)), c(gampred$upper, rev(gampred$lower)), col = transparent_color, border = NA)
  lines(gampred$xs, gampred$preds.fit)  # Add main line after to ensure it's on top
  # dev.off()
}


psb <- array(,dim = c(2,2)) #2(L,R)x(peakGAM,peakRMA)
for (ipeak in 1:2) {
  # for (i in 1:length(bs)) {
  b <- bs[[3]]
  bsum <- summary(b)
  # if (bsum$s.pv > 0.05) {next} #module --ignore_cache load "stack/2022.1"
  datab <- b$model
  datab$var <- 1/datab$`(weights)`
  datab$weight <- datab$`(weights)`
  datab$sqrtage <- sqrt(datab$age)
  predictb <- predict(b, newdata = data.frame(
    age = datab$age,
    TaskCode_1 = 0,
    TaskCode_2 = 0,
    TaskCode_3 = 0,
    HandnessCode_1 = 0,
    HandnessCode_2 = 0,
    ContrastCode_1 = 0,
    ContrastCode_2 = 0,
    ErrorTrialCode_1 = 0,
    ErrorTrialCode_2 = 0,
    DesignCode_1 = 0,
    SRC = 0,
    SRCcode = 0
  ))
  if (ipeak == 1) {
    idx <- which(predictb==min(predictb))
  } else {
    peak <- peak_rma
    idx <- which.min(abs(datab$age-peak))
  }
  
  if (length(idx) == 1) {
    databL <- datab[1:idx,]
    databR <- datab[idx:length(predictb),]
  } else {
    databL <- datab[1:idx[1],]
    databR <- datab[idx[-1]:length(predictb),]
  }
  
  # fit linear model and see if the slope is significant
  step = 0.1
  for (iana in 1:2) { #left and right
    if (iana == 1) {databX <- databL} else {databX <- databR}
    bX <- gam(laterality~age + TaskCode_1 + TaskCode_2 + TaskCode_3 + HandnessCode_1 + HandnessCode_2 + ContrastCode_1 + ContrastCode_2 + ErrorTrialCode_1 + DesignCode_1, data=databX, weights = weight, method = 'REML')
    bXsum <- summary(bX)
    psb[iana,ipeak] <- bXsum$p.pv[2]
  }
  # }
}

psb2 = psb/2
