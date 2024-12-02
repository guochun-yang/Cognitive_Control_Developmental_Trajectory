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
source('./naturalSortFunctions.R')
source('./pmtx_adjust.R')

## load data and calculate the averaged betas and variances
datadir = shortPathName("../data/") 
sdmdir = '../data/'

filenames = list.files(path = paste(datadir,'extracted_data_contrastanalysis',sep = ''),
                       pattern = "^multivoxel_extract_MyLinearModel_adult_childolder_good_z_voxelCorrected_p_0.00100_1_blob.*\\.txt$")
filenames = natural_sort(filenames,'blob','.txt')

savedir = shortPathName("../plot/blob_adult-chilolder_good+design+SRC/")
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

## load the supplementary table to read the covariates
supptable = read_excel(
  '../data/Table.xlsx',
  sheet = 'Sheet2',
  skip = 1
)


logliks = matrix(, nrow = nfile, ncol = 5)
AICs = matrix(, nrow = nfile, ncol = 5)
BICs = matrix(, nrow = nfile, ncol = 5)
QMps = matrix(, nrow = nfile, ncol = 5)
modelweights = matrix(, nrow = nfile, ncol = 5)

models = list()


Fs_gam <- ps_gam <- data2s <- bs <- reses_qua <- reses_cub <- reses_log <- reses_sqrt <- reses_revlog <- reses_revsqrt <- reses_lin <- c()
Fs <- ps <- df1s <- df2s <- array(,dim = c(nfile))
QMss <- zss <- array(, dim = c(nfile, 5)) #ROI, 4, #perm
for (i in 1:nfile) {
  data = read.table(
    paste(datadir, 'extracted_data_contrastanalysis', filenames[i], sep = '/'),
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
    data2$SRC[istudy] <- supptable$SRC[idx]
    
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
                                          'DesignCode_1','DesignCode_0'))
  
  for (icol in 1:length(idx_col)) {
    data2[,idx_col[icol]] = as.vector(scale(data2[,idx_col[icol]]))
  }
  
  # delete outliers
  outup = mean(data2$beta) + 3*sd(data2$beta)
  outdown = mean(data2$beta) - 3*sd(data2$beta)
  data2 = data2[-which(data2$beta>outup | data2$beta < outdown),]
  data2s[[i]] = data2
  
  # sort data2 with the age
  data2 <- data2[order(data2$age), ]
  
  ## select analysis type
  data2$var = data2$variance
  data2$weight = 1/data2$variance
  
  # randomly choose the age from the noramal distribution of (age, varirance) for 1000 times
  options(digits = 7) #7 is default
  
  ages_perm <- array(, dim = c(dim(data2)[1]))
  
  data2$age <- data2$age0
  
  
  # # as.numeric(data2$SRC)
  # # filter rows
  # data3 <- subset(data2,!is.na(SRC))
  # data3$SRC = as.vector(scale(data3$SRC))
  
  
  # models
  # b1 <- gam(beta~s(age) + TaskCode_1 + TaskCode_2 + TaskCode_3 + HandnessCode_1 + HandnessCode_2 + ContrastCode_1 + ContrastCode_2 + ErrorTrialCode_1 + ErrorTrialCode_2 + DesignCode_1,data=data2, weights = weight, method = 'REML')
  b <- gam(beta~s(age) + TaskCode_1 + TaskCode_2 + TaskCode_3 + HandnessCode_1 + HandnessCode_2 + ContrastCode_1 + ContrastCode_2 + ErrorTrialCode_1 + ErrorTrialCode_2 + DesignCode_1 + SRC,data=data2, weights = weight, method = 'REML')
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
  
  
  
  step = 0.1 #avoid not converging
  
  zvals <- pvals <- bvals <- c()
  
  data2$age <- data2$age0
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