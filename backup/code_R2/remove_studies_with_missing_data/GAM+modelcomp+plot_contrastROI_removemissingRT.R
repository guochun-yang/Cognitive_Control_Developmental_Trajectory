rm(list=ls())
library(metafor)
library(stringr)
library(performance)
library(car)
library(visreg)
library(MuMIn)
library(mgcv)
library(MASS)
library("readxl")
library(fastDummies)
# library(mice)
source('H:/meta_development/Results_R1/SDM/scripts_NHB/naturalSortFunctions.R')
source('H:/meta_development/Results_R1/SDM/scripts_NHB/pmtx_adjust.R')

## load data and calculate the averaged betas and variances
datadir = shortPathName("H:/meta_development/Results_R1/SDM/mean/analysis_Adult-ChildOld/") 
sdmdir = 'H:/meta_development/Results_R1/SDM/mean/' #sdm_good

filenames = list.files(path = paste(datadir,'analysis_MyLinearModel/extracts/backup',sep = ''),
                       pattern = "^multivoxel_extract_MyLinearModel_adult_childolder_good_z_voxelCorrected_p_0.00100_1_blob.*\\.txt$")
filenames = natural_sort(filenames,'blob','.txt')

savedir = shortPathName("H:/meta_development/Results_R2/SDM/analysis_Adult-ChildOld/plot/blob_adult-chilolder_removemissingRT/")
if (!dir.exists(savedir)) {
  dir.create(savedir, recursive = TRUE)
}

nfile = length(filenames)

## load SDM table
sdmtable = read.table(
  paste(sdmdir, '/sdm_table.txt', sep = ''),
  header = TRUE,
  sep = '',
  dec = '.'
)

## load the covariates
supptable = read_excel(
  'H:/meta_development/Results_R2/SDM/scripts_NHB/covariates.xlsx',
  sheet = 'covariates'
)

logliks = matrix(, nrow = nfile, ncol = 5)
AICs = matrix(, nrow = nfile, ncol = 5)
BICs = matrix(, nrow = nfile, ncol = 5)
QMps = matrix(, nrow = nfile, ncol = 5)
modelweights = matrix(, nrow = nfile, ncol = 5)

models = list()


Fs_gam <- ps_gam <- data2s <- bs <- devexp_gam <- reses_qua <- reses_cub <- reses_log <- reses_sqrt <- reses_revlog <- reses_revsqrt <- reses_lin <- c()
Fs <- ps <- df1s <- df2s <- array(,dim = c(nfile))
QMss <- zss <- array(, dim = c(nfile, 5)) #ROI, 4, #perm
for (i in 1:nfile) {
  data = read.table(
    paste(datadir, 'analysis_MyLinearModel/extracts/backup', filenames[i], sep = '/'),
    header = FALSE,
    sep = '',
    dec = '.',
    skip = 13
  )
  
  data$study = sub("_I000.*", "", data$V1)
  data$betas = data$V2
  
  data$variances = data$V3
  # data$weight = 1/data$variances
  
  #screen data through sdm table
  idx_incl <- which(data$study %in% sdmtable$study)
  data_good <- data[idx_incl,] 
  
  beta = aggregate(data_good$betas, list(data_good$study), FUN = mean)
  variance = aggregate(data_good$variances, list(data_good$study), FUN = mean)
  data2 <- data.frame(beta, variance$x)
  
  
  colnames(data2)[1] <- 'study'
  colnames(data2)[2] <- 'beta'
  colnames(data2)[3] <- 'variance'
  
  # add corresponding covariates
  for (istudy in 1:length(data2$study)) {
    # print(istudy)
    idx <- which(supptable$AuthoYear == data2$study[istudy])
    data2$varianceAgerange[istudy] <- as.numeric(supptable$Variance[idx])
    # data2$varianceAdjusted[istudy] <- data2$variance[istudy] * as.numeric(supptable$Variance[idx]) #adjust the variance by multiply the two variance sources
    data2$TaskCode[istudy] <- supptable$TaskCode[idx]
    data2$HandnessCode[istudy] <- supptable$HandnessCode[idx]
    data2$ContrastCode[istudy] <- supptable$ContrastCode[idx]
    data2$ErrorTrialCode[istudy] <- supptable$ErrorTrialCode[idx]
    data2$DesignCode[istudy] <- supptable$DesignCode[idx]
    data2$SRC[istudy] <- supptable$SRC0[idx]
    data2$SRCcode[istudy] <- supptable$SRCcode[idx]
    
    #add age from sdmtable
    idx2 <- which(sdmtable$study == data2$study[istudy])
    data2$age[istudy] <- sdmtable$avgAge[idx2]
    data2$age0[istudy] <- sdmtable$avgAge[idx2]
  }

  
  ## add dummy columns
  data2 <- dummy_cols(data2, 
                      select_columns = c("TaskCode","HandnessCode","ContrastCode","ErrorTrialCode","DesignCode"))
  
  ## z-score these columns
  idx_col <- which(colnames(data2) %in% c('TaskCode_1','TaskCode_2','TaskCode_3','TaskCode_4',
                                          'HandnessCode_1','HandnessCode_2','HandnessCode_3',
                                          'ContrastCode_1','ContrastCode_2','ContrastCode_3',
                                          'ErrorTrialCode_1','ErrorTrialCode_2','ErrorTrialCode_3',
                                          'DesignCode_1','DesignCode_0','SRCcode'))
  
  for (icol in 1:length(idx_col)) {
    data2[,idx_col[icol]] = as.vector(scale(data2[,idx_col[icol]]))
    # data2[,idx_col[icol]] = as.vector(data2[,idx_col[icol]])
  }
  
  # delete outliers
  outup = mean(data2$beta) + 3*sd(data2$beta)
  outdown = mean(data2$beta) - 3*sd(data2$beta)
  data2 = data2[-which(data2$beta>outup | data2$beta < outdown),]
  
  # sort data2 with the age
  data2 <- data2[order(data2$age), ]
  
  ## select analysis type
  data2$var = data2$variance
  data2$weight = 1/data2$variance
  
  # randomly choose the age from the noramal distribution of (age, varirance) for 1000 times
  options(digits = 7) #7 is default
  
  ages_perm <- array(, dim = c(dim(data2)[1]))
  
  data2$age <- data2$age0
  
  
  # filter rows
  data2 <- subset(data2,!is.na(SRC))
  data2$SRC = as.vector(scale(data2$SRC))
  data2s[[i]] = data2
  
  
  b <- gam(beta~s(age) + TaskCode_1 + TaskCode_2 + TaskCode_3 + HandnessCode_1 + HandnessCode_2 + ContrastCode_1 + ContrastCode_2 + ErrorTrialCode_1 + ErrorTrialCode_2 + DesignCode_1 + SRC,data=data2, weights = weight, method = 'REML')
  
  # # multiple imputation with mice
  # imputed_data <- mice(data2, m = 5, method = c("", "", "", "", "", "", "", "", "", "pmm","", "", "", "", "", "", "", "", "", "","", "", "", "", "", "", "", "", "", ""), seed = 123)
  # models <- with(imputed_data, gam(beta ~ s(age) + TaskCode_1 + TaskCode_2 + TaskCode_3 + 
  #                                    HandnessCode_1 + HandnessCode_2 + ContrastCode_1 + 
  #                                    ContrastCode_2 + ErrorTrialCode_1 + ErrorTrialCode_2 + 
  #                                    DesignCode_1 + SRC, weights = weight, 
  #                                  method = 'REML'))
  # pooled_results <- pool(models)
  # summary(pooled_results)
  
  sumb <- summary(b)
  F <- sumb$s.table[1,'F']
  Fs[i] <- F
  p <- sumb$s.table[1,'p-value']
  ps[i] <- p
  df1 <- sumb$edf
  df2 <- b$df.residual
  df1s[i] = df1
  df2s[i] = df2
  
  bs[[i]] = b
  # sumb<-summary(b)
  ps_gam[i] <- p #sumb$s.table[1,'p-value']
  Fs_gam[i] <- F #sumb$s.table[1,'F']
  devexp_gam[i] <- sumb$dev.expl
  
  
  
  step = 0.1 #avoid not converging
  
  zvals <- pvals <- bvals <- c()
  
  res.base <- 
    rma(beta, 
        var, 
        mods = ~ TaskCode_1 + TaskCode_2 + TaskCode_3 + HandnessCode_1 + HandnessCode_2 + ContrastCode_1 + ContrastCode_2 + ErrorTrialCode_1 + DesignCode_1 + SRC,
        data = data2,
        control=list(stepadj=step))
  if (sum(data2$ErrorTrialCode_2) == 0) {
    res.qua <- 
      rma(beta, 
          var, 
          mods = ~ age + I(age^2)
          + TaskCode_1 + TaskCode_2 + TaskCode_3 + HandnessCode_1 + HandnessCode_2 + ContrastCode_1 + ContrastCode_2 + ErrorTrialCode_1 + DesignCode_1 + SRC,
          data = data2,
          control=list(stepadj=step))
    res.cub <-
      rma(beta,
          var,
          mods = ~ age + I(age^2) + I(age^3)
          + TaskCode_1 + TaskCode_2 + TaskCode_3 + HandnessCode_1 + HandnessCode_2 + ContrastCode_1 + ContrastCode_2 + ErrorTrialCode_1 + DesignCode_1 + SRC, 
          data = data2,
          control=list(stepadj=step))
    res.log <-
      rma(beta,
          var,
          mods = ~ I(log(age)) + I((log(age))^2)
          + TaskCode_1 + TaskCode_2 + TaskCode_3 + HandnessCode_1 + HandnessCode_2 + ContrastCode_1 + ContrastCode_2 + ErrorTrialCode_1 + DesignCode_1 + SRC, 
          data = data2,
          control=list(stepadj=step))
    res.sqrt <-
      rma(beta,
          var,
          mods = ~ I(sqrt(age)) + age
          + TaskCode_1 + TaskCode_2 + TaskCode_3 + HandnessCode_1 + HandnessCode_2 + ContrastCode_1 + ContrastCode_2 + ErrorTrialCode_1 + DesignCode_1 + SRC, 
          data = data2,
          control=list(stepadj=step))
    res.lin <-
      rma(beta,
          var,
          mods = ~ age
          + TaskCode_1 + TaskCode_2 + TaskCode_3 + HandnessCode_1 + HandnessCode_2 + ContrastCode_1 + ContrastCode_2 + ErrorTrialCode_1 + DesignCode_1 + SRC, 
          data = data2,
          control=list(stepadj=step))
  } else {
    res.qua <- 
      rma(beta, 
          var, 
          mods = ~ age + I(age^2)
          + TaskCode_1 + TaskCode_2 + TaskCode_3 + HandnessCode_1 + HandnessCode_2 + ContrastCode_1 + ContrastCode_2 + ErrorTrialCode_1 + ErrorTrialCode_2 + DesignCode_1 + SRC,
          data = data2,
          control=list(stepadj=step))
    res.cub <-
      rma(beta,
          var,
          mods = ~ age + I(age^2) + I(age^3)
          + TaskCode_1 + TaskCode_2 + TaskCode_3 + HandnessCode_1 + HandnessCode_2 + ContrastCode_1 + ContrastCode_2 + ErrorTrialCode_1 + ErrorTrialCode_2 + DesignCode_1 + SRC, 
          data = data2,
          control=list(stepadj=step))
    res.log <-
      rma(beta,
          var,
          mods = ~ I(log(age)) + I((log(age))^2)
          + TaskCode_1 + TaskCode_2 + TaskCode_3 + HandnessCode_1 + HandnessCode_2 + ContrastCode_1 + ContrastCode_2 + ErrorTrialCode_1 + ErrorTrialCode_2 + DesignCode_1 + SRC, 
          data = data2,
          control=list(stepadj=step))
    res.sqrt <-
      rma(beta,
          var,
          mods = ~ I(sqrt(age)) + age
          + TaskCode_1 + TaskCode_2 + TaskCode_3 + HandnessCode_1 + HandnessCode_2 + ContrastCode_1 + ContrastCode_2 + ErrorTrialCode_1 + ErrorTrialCode_2 + DesignCode_1 + SRC, 
          data = data2,
          control=list(stepadj=step))
    res.lin <-
      rma(beta,
          var,
          mods = ~ age
          + TaskCode_1 + TaskCode_2 + TaskCode_3 + HandnessCode_1 + HandnessCode_2 + ContrastCode_1 + ContrastCode_2 + ErrorTrialCode_1 + ErrorTrialCode_2 + DesignCode_1 + SRC, 
          data = data2,
          control=list(stepadj=step))
    # }
    QMss[i,1] <- res.qua$QM
    QMss[i,2] <- res.cub$QM
    QMss[i,3] <- res.log$QM
    QMss[i,4] <- res.sqrt$QM
    QMss[i,5] <- res.lin$QM
    zss[i,1] <- res.qua$zval[3]
    zss[i,2] <- res.cub$zval[4]
    zss[i,3] <- res.log$zval[3]
    zss[i,4] <- res.sqrt$zval[2]
    zss[i,5] <- res.lin$zval[2]
    
    
    reses_qua[[i]] = res.qua
    reses_cub[[i]] = res.cub
    reses_log[[i]] = res.log
    reses_sqrt[[i]] = res.sqrt
    reses_lin[[i]] = res.lin
    
    logliks[i, 1] = as.numeric(logLik.rma(res.qua))
    logliks[i, 2] = as.numeric(logLik.rma(res.cub))
    logliks[i, 3] = as.numeric(logLik.rma(res.log))
    logliks[i, 4] = as.numeric(logLik.rma(res.sqrt))
    logliks[i, 5] = as.numeric(logLik.rma(res.lin))
    
    AICs[i, 1] = AIC(res.qua)
    AICs[i, 2] = AIC(res.cub)
    AICs[i, 3] = AIC(res.log)
    AICs[i, 4] = AIC(res.sqrt)
    AICs[i, 5] = AIC(res.lin)
    
    BICs[i, 1] = BIC(res.qua)
    BICs[i, 2] = BIC(res.cub)
    BICs[i, 3] = BIC(res.log)
    BICs[i, 4] = BIC(res.sqrt)
    BICs[i, 5] = BIC(res.lin)
    
    QMps[i, 1] = res.qua$QMp
    QMps[i, 2] = res.cub$QMp
    QMps[i, 3] = res.log$QMp
    QMps[i, 4] = res.sqrt$QMp
    QMps[i, 5] = res.lin$QMp
    
    
    mod.comparison <- AICc(res.qua, res.cub, res.log, res.sqrt)
    modelweights[i,1:4] <- Weights(mod.comparison)
  }
}

# # check results
p_gam_fdr <- p.adjust(ps_gam, method = 'fdr', n = length(ps_gam))

p_values <- p_zs <- zval_meanage <- zval_randage <- array(,dim = c(nfile,5))
F_meanage <- F_randage <- array(,dim = c(nfile,1))
b_sqrt <- bl_sqrt <- bu_sqrt <- array(,dim = c(nfile,1))
for (i in 1:nfile) {
    
  # calculate the modulator effect (including age and others)
  p_value <- pchisq(mean(QMss[i,1],na.rm=TRUE), reses_qua[[i]]$QMdf[1], lower.tail = FALSE)
  p_values[i,1] <- p_value
  p_value <- pchisq(mean(QMss[i,2],na.rm=TRUE), reses_cub[[i]]$QMdf[1], lower.tail = FALSE)
  p_values[i,2] <- p_value
  p_value <- pchisq(mean(QMss[i,3],na.rm=TRUE), reses_log[[i]]$QMdf[1], lower.tail = FALSE)
  p_values[i,3] <- p_value
  p_value <- pchisq(mean(QMss[i,4],na.rm=TRUE), reses_sqrt[[i]]$QMdf[1], lower.tail = FALSE)
  p_values[i,4] <- p_value
  p_value <- pchisq(mean(QMss[i,5],na.rm=TRUE), reses_lin[[i]]$QMdf[1], lower.tail = FALSE)
  p_values[i,5] <- p_value
  
  #calculate the age modulator effect in each model
  p_z <- pnorm(mean(zss[i,1],na.rm=TRUE), lower.tail = TRUE)
  p_zs[i,1] <- p_z
  p_z <- pnorm(mean(zss[i,2],na.rm=TRUE), lower.tail = FALSE)
  p_zs[i,2] <- p_z
  p_z <- pnorm(mean(zss[i,3],na.rm=TRUE), lower.tail = TRUE)
  p_zs[i,3] <- p_z
  p_z <- pnorm(mean(zss[i,4],na.rm=TRUE), lower.tail = FALSE)
  p_zs[i,4] <- p_z
  p_z <- pnorm(mean(zss[i,5],na.rm=TRUE), lower.tail = FALSE)
  p_zs[i,5] <- p_z
  
  zval_meanage[i,1] <- reses_qua[[i]]$zval[3]
  zval_randage[i,1] <- mean(zss[i,1],na.rm=TRUE)
  zval_meanage[i,2] <- reses_cub[[i]]$zval[4]
  zval_randage[i,2] <- mean(zss[i,2],na.rm=TRUE)
  zval_meanage[i,3] <- reses_log[[i]]$zval[3]
  zval_randage[i,3] <- mean(zss[i,3],na.rm=TRUE)
  zval_meanage[i,4] <- reses_sqrt[[i]]$zval[2]
  zval_randage[i,4] <- mean(zss[i,4],na.rm=TRUE)
  zval_meanage[i,5] <- reses_lin[[i]]$zval[2]
  zval_randage[i,5] <- mean(zss[i,5],na.rm=TRUE)
  b_sqrt[i,1] <- reses_sqrt[[i]]$beta[2]
  bl_sqrt[i,1] <- reses_sqrt[[i]]$ci.lb[2]
  bu_sqrt[i,1] <- reses_sqrt[[i]]$ci.ub[2]
  
  F_meanage[i,1] <- Fs_gam[i]
  F_randage[i,1] <- mean(Fs[i],na.rm=TRUE)
}


options(scipen=999)
# pss

minAIC <- minBIC <- maxLogLik <- c()
for (i in 1:nfile) {
  minAIC[i] = which.min(AICs[i,1:4])
  minBIC[i] = which.min(AICs[i,1:4])
  maxLogLik[i] = which.max(logliks[i,1:4])
}


modelweights2 <- format(round(modelweights, 3), nsmall = 3)
write.matrix(modelweights2, file = paste(savedir,"modelweights.csv"), sep = ',')

# get effect sizes and confident interval
ds_sqrt <- ds_qua <- r_sqs <- array(dim = c(nfile,1),)
cis_sqrt <- cis_qua <- array(dim = c(nfile,2),)
for (i in 1:9) {
  ds_sqrt[i] <- coef(reses_sqrt[[i]])[2]
  ds_qua[i] <- coef(reses_qua[[i]])[3]
  cis_sqrt[i,1] <- reses_sqrt[[i]]$ci.lb[2]
  cis_sqrt[i,2] <- reses_sqrt[[i]]$ci.ub[2]
  cis_qua[i,1] <- reses_qua[[i]]$ci.lb[3]
  cis_qua[i,2] <- reses_qua[[i]]$ci.ub[3]
  b <- summary(bs[[i]])
  r_sqs[i] <- b$r.sq #explaned ratio
}


###  plot the figures  ###
xs <- seq(8, 75, length=500)

minmax = minAIC
peak_rmas <- peak_gams <- c()
for (i in 1:nfile) {

  mtx_cov_xs <- matrix(0,nrow = length(xs), ncol = length(reses_qua[[i]]$beta)-3)
  if (minmax[i] == 1) {
    modelbest = reses_qua[[i]]
    title = "QUADRATIC"
    sav <-
      predict(modelbest, newmods = unname(cbind(poly(
        xs, degree = 2, raw = TRUE
      ),mtx_cov_xs)))
    peak_rma = -0.5*modelbest$b[2,1]/modelbest$b[3,1]
    if (modelbest$pval[3]/2 > 0.05) {
      sav$pred <- NaN
      peak_rma <- NaN
    }
  } else if (minmax[i] == 2) {
    modelbest = reses_cub[[i]]
    title = "CUBIC"
    sav <-
      predict(modelbest, newmods = unname(cbind(poly(
        xs, degree = 3, raw = TRUE
      ),mtx_cov_xs)))
    if (modelbest$b[4,1] > 0) {
      peak_rma = (-2*modelbest$b[3,1] - sqrt(4*(modelbest$b[3,1])^2 - 12*modelbest$b[2,1]*modelbest$b[4,1]))/(6*modelbest$b[4,1])
    } else {
      peak_rma = (-2*modelbest$b[3,1] + sqrt(4*(modelbest$b[3,1])^2 - 12*modelbest$b[2,1]*modelbest$b[4,1]))/(6*modelbest$b[4,1])
    }
    if (modelbest$pval[4]/2 > 0.05) {
      sav$pred <- NaN
      peak_rma <- NaN
    }
  } else if (minmax[i] == 3) {
    modelbest = reses_log[[i]]
    title = "LOG"
    sav <-
      predict(modelbest, newmods = unname(cbind(poly(
        log(xs), degree = 2, raw = TRUE
      ),mtx_cov_xs)))
    peak_rma = exp(-modelbest$b[2,1]/(2*modelbest$b[3,1]))
    if (modelbest$pval[3]/2 > 0.05) {
      sav$pred <- NaN
      peak_rma <- NaN
    }
  } else if (minmax[i] == 4) {
    modelbest = reses_sqrt[[i]]
    title = "SQRT"
    sav <-
      predict(modelbest, newmods = unname(cbind(poly(
        sqrt(xs), degree = 2, raw = TRUE
      ),mtx_cov_xs)))
    peak_rma = 0.25*(modelbest$b[2,1]/modelbest$b[3,1])^2
    if (modelbest$pval[2]/2 > 0.05) {
      sav$pred <- NaN
      peak_rma <- NaN
    }
  } else if (minmax[i] == 5) {
    modelbest = reses_lin[[i]]
    title = "LIN"
    sav <-
      predict(modelbest, newmods = unname(cbind(poly(
        xs, degree = 1, raw = TRUE
      ),mtx_cov_xs)))
    peak_rma = NaN
    # if (modelbest$pval[2]/2 > 0.05 && p_gam_fdr[i] < 0.05) {
    #   sav$pred <- NaN
    # }
  }
  peak_rmas[i] = peak_rma

  
  ylim0 = c(-0.2,1.0)
  ylim1 = c(-0.18,1.02)
  
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
  
  peak_gam <- xs[which.max(preds$fit)]
  if (peak_gam > 60 || p_gam_fdr[i] > 0.05) {
    peak_gam <- NaN
  }
  peak_gams[i] = peak_gam #round(peak_gam,1)
  
  gampred <- data.frame(xs,preds$fit,preds$se.fit)
  gampred$upper <- gampred$preds.fit + gampred$preds.se.fit
  gampred$lower <- gampred$preds.fit - gampred$preds.se.fit
  
  if (p_gam_fdr[i] > 0.05) {
    gampred$preds.fit <- gampred$upper <- gampred$lower <- NaN
  }
  

### plot all at once
# filename = str_match(filenames[i], "_p_0.00100_10_\\s*(.*?)\\s*.txt")
filename = str_match(filenames[i], "extract_\\s*(.*?)\\s*.txt")
tiff(
  file = paste(savedir, filename[, 2], "_sameylim.tiff", sep = ''),
  width = 6,
  height = 5,
  units = "in",
  res = 300
)

# windows()
par(cex.lab=2)
par(cex.axis=2)
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
  xlab = "Age (years)",ylab = "Effect size",
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
# lines(xs0, preds0$fit, type = 'l', lwd = 3, col = 'blue')  # or type = 'n' if you just want to set up axes
# polygon(c(xs0, rev(xs0)), c(preds0$fit + preds0$se.fit, rev(preds0$fit - preds0$se.fit)), col = transparent_color, border = NA)
# lines(xs0, preds0$fit)  # Add main line after to ensure it's on top

if (i %in% c(1,4,5,6,7,8,9)) {
  par(new=TRUE)
  lines(c(peak_rma,peak_rma),ylim0,lty = 2,col = 'red',lwd = 2)
}
if (i %in% c(3)) {
  par(new=TRUE)
  lines(c(peak_rma,peak_rma),ylim1,lty = 2,col = 'red',lwd = 2) # in fig 3 the two peaks overlap
}
par(new=TRUE)
lines(c(peak_gam,peak_gam),ylim0,lty = 2,col = 'blue',lwd = 2)
dev.off()
}
r_peak <- cor.test(peak_rmas,peak_gams)



## test the significance of peaks with the Simonsohn (2018)'s two-line approach. 
# Simonsohn, U. (2018) Two Lines: A Valid Alternative to the Invalid Testing of U-Shaped Relationships With Quadratic Regressions. Advances in Methods and Practices in Psychological Science, 1(4), 538-555.
psb <- array(,dim = c(length(bs),2,2)) #9x2(L,R)x(peakGAM,peakRMA)
for (ipeak in 1:2) {
  for (i in 1:length(bs)) {
    b <- bs[[i]]
    bsum <- summary(b)
    if (bsum$s.pv > 0.05) {psb[i,iana,ipeak] = NA; next}
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
      SRC = 0
    ))
    if (ipeak == 1) {
      idx <- which(predictb==max(predictb))
    } else {
      peak <- peak_rmas[i]
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
      bX <- gam(beta~age + TaskCode_1 + TaskCode_2 + TaskCode_3 + HandnessCode_1 + HandnessCode_2 + ContrastCode_1 + ContrastCode_2 + ErrorTrialCode_1 + ErrorTrialCode_2 + DesignCode_1 + SRC,data=databX, weights = weight, method = 'REML')
      bXsum <- summary(bX)
      psb[i,iana,ipeak] <- bXsum$p.pv[2]
    }
  }
}
psb_fdr <- array(dim = dim(psb))
for (ipeak in 1:2) {
  psb_fdr[,,ipeak] <- t(apply(psb[,,ipeak]/2, 1, function(x) p.adjust(x, method = "fdr")))
}