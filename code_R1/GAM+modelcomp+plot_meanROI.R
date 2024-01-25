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
source('H:/meta_development/Results_R1/SDM/scripts/naturalSortFunctions.R')

## load data and calculate the averaged betas and variances
datadir = shortPathName("H:/meta_development/Results_R1/SDM/mean/analysis_Adult-ChildOld/") 
sdmdir = 'H:/meta_development/Results_R1/SDM/mean/' #sdm_good

filenames = list.files(path = paste(datadir,'analysis_MyLinearModel/extracts',sep = ''),
                       pattern = "^multivoxel_extract_mean_good_z_voxelCorrected_p_0.00100_10_blob.*\\.txt$")
filenames = natural_sort(filenames,'blob','.txt')

savedir = shortPathName("H:/meta_development/Results_R1/SDM/mean/analysis_Adult-ChildOld/plot/blob_mean_good/")
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
  'C:/Users/Guochun Yang/OneDrive - University of Iowa/WithLi/Paper/data_share/R1/SRC元分析文献汇总_lizh_ygc.xlsx',
  sheet = 'TableS1-整合',
  skip = 1
)

logliks = matrix(, nrow = nfile, ncol = 5)
AICs = matrix(, nrow = nfile, ncol = 5)
BICs = matrix(, nrow = nfile, ncol = 5)
QMps = matrix(, nrow = nfile, ncol = 5)
QMs = matrix(, nrow = nfile, ncol = 5)
modelweights = matrix(, nrow = nfile, ncol = 5)

models = list()

# nperm = 500
# iana=2
# addcov = iana-1

Fs_gam <- ps_gam <- data2s <- bs <- reses_qua <- reses_cub <- reses_log <- reses_sqrt <- reses_revlog <- reses_revsqrt <- reses_lin <- c()
Fs <- ps <- df1s <- df2s <- array(,dim = c(nfile))
QMss <- zss <- array(, dim = c(nfile, 5)) #ROI, 4, #perm
for (i in 1:nfile) {
  data = read.table(
    paste(datadir, 'analysis_MyLinearModel/extracts', filenames[i], sep = '/'),
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
    
    #add age from sdmtable
    idx2 <- which(sdmtable$study == data2$study[istudy])
    data2$age[istudy] <- sdmtable$avgAge[idx2]
    data2$age0[istudy] <- sdmtable$avgAge[idx2]
  }
  
  ## add dummy columns
  data2 <- dummy_cols(data2, 
                      select_columns = c("TaskCode","HandnessCode","ContrastCode","ErrorTrialCode"))
  
  ## z-score these columns
  idx_col <- which(colnames(data2) %in% c('TaskCode_1','TaskCode_2','TaskCode_3','TaskCode_4',
                                          'HandnessCode_1','HandnessCode_2','HandnessCode_3',
                                          'ContrastCode_1','ContrastCode_2','ContrastCode_3',
                                          'ErrorTrialCode_1','ErrorTrialCode_2','ErrorTrialCode_3'))
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
  
    # if (iperm == nperm+1) {
      data2$age <- data2$age0
    # } else{
    #   for (istudy in 1:length(data2$study)) {
    #     repeat {
    #       randage <- rnorm(n = 1, mean = data2$age0[istudy], sd = data2$varianceAgerange[istudy])
    #     # data2$age[istudy] <- runif(n = 1, min = data2$age0[istudy] - 2*data2$varianceAgerange[istudy], max = data2$age0[istudy] + 2*data2$varianceAgerange[istudy])
    #       # if (randage > data2$age0[istudy] + data2$varianceAgerange[istudy] || randage < data2$age0[istudy] - data2$varianceAgerange[istudy]) {
    #       if (randage < data2$age0[istudy] + data2$varianceAgerange[istudy] && randage > data2$age0[istudy] - data2$varianceAgerange[istudy]) {
    #         break
    #       }
    #     }
    #     data2$age[istudy] <- randage
    #   }
    # }
    # ages_perm[,iperm] <- data2$age
    
    
    
    # models
    # b <- gam(beta~s(age) + TaskCode_1 + TaskCode_2 + TaskCode_3,data=data2, weights = weight, method = 'REML')
    b <- gam(beta~s(age) + TaskCode_1 + TaskCode_2 + TaskCode_3 + HandnessCode_1 + HandnessCode_2 + ContrastCode_1 + ContrastCode_2 + ErrorTrialCode_1 + ErrorTrialCode_2,data=data2, weights = weight, method = 'REML')
    # b2 <- gam(beta~s(age) + TaskCode_1 + TaskCode_2 + TaskCode_3, data=data2, weights = weight, method = 'REML')
    
    # b2 <- gam(beta~s(age) + s(TaskCode_1) + s(TaskCode_2) + s(TaskCode_3) + s(HandnessCode_1) + s(HandnessCode_2) + s(ContrastCode_1) + s(ContrastCode_2) + s(ErrorTrialCode_1) + s(ErrorTrialCode_2),data=data2, weights = weight, method = 'REML')
    sumb <- summary(b)
    F <- sumb$s.table[1,'F']
    Fs[i] <- F
    p <- sumb$s.table[1,'p-value']
    ps[i] <- p
    df1 <- sumb$edf
    df2 <- b$df.residual
    df1s[i] = df1
    df2s[i] = df2
    #later we will use the mean of F and the df1/df2 to generate the mean p value
    # pf(q = F, df1 = df1, df2 = df2, lower.tail = FALSE)
    
    # if (iperm == nperm+1) {
      # data2$age <- data2$age0
      # b <- gam(beta~s(age) + TaskCode_1 + TaskCode_2 + TaskCode_3 + HandnessCode_1 + HandnessCode_2 + ContrastCode_1 + ContrastCode_2 + ErrorTrialCode_1 + ErrorTrialCode_2,data=data2, weights = weight, method = 'REML')
      bs[[i]] = b
      # sumb<-summary(b)
      ps_gam[i] <- p #sumb$s.table[1,'p-value']
      Fs_gam[i] <- F #sumb$s.table[1,'F']
    # }
  
  # mean(Fs[1,1:500],na.rm=TRUE)
  # exp(mean(log(Fs[1,1:500]),na.rm=TRUE))
  # Fs[1,501]
  # mean(ps[1,1:500],na.rm=TRUE)
  # ps[1,501]

  
  
  
  step = 0.1 #avoid not converging
  # res.qua <- rma(beta, variance, mods = ~ age + I(age^2), data = data2)
  
  zvals <- pvals <- bvals <- c()
  # for (iperm in 1:nperm+1) {
  #   print(iperm)
    # if (iperm == nperm+1) {
      data2$age <- data2$age0
    # } else{
    #   for (istudy in 1:length(data2$study)) {
    #     data2$age[istudy] <- rnorm(n = 1, mean = data2$age0[istudy], sd = data2$varianceAgerange[istudy])
    #   }
    # }
    # if (addcov == 0) {
    #   res.qua <- 
    #     rma(beta, 
    #         var, 
    #         mods = ~ age + I(age^2),
    #         data = data2,
    #         control=list(stepadj=step))
    #   res.cub <-
    #     rma(beta,
    #         var,
    #         mods = ~ age + I(age^2) + I(age^3),
    #         data = data2,
    #         control=list(stepadj=step))
    #   res.log <-
    #     rma(beta,
    #         var,
    #         mods = ~ I(log(age)) + I((log(age))^2),
    #         data = data2,
    #         control=list(stepadj=step))
    #   res.sqrt <-
    #     rma(beta,
    #         var,
    #         mods = ~ I(sqrt(age)) + age,
    #         data = data2,
    #         control=list(stepadj=step))
    #   res.lin <-
    #     rma(beta,
    #         var,
    #         mods = ~ age,
    #         data = data2,
    #         control=list(stepadj=step))
    #   res.softplus <-
    #     rma(beta,
    #         var,
    #         mods = ~ I(log(1+exp(age))),
    #         data = data2,
    #         control=list(stepadj=step))
    # } else {
      if (sum(data2$ErrorTrialCode_2) == 0) {
        res.qua <- 
          rma(beta, 
              var, 
              mods = ~ age + I(age^2)
              + TaskCode_1 + TaskCode_2 + TaskCode_3 + HandnessCode_1 + HandnessCode_2 + ContrastCode_1 + ContrastCode_2 + ErrorTrialCode_1,
              data = data2,
              control=list(stepadj=step))
        res.cub <-
          rma(beta,
              var,
              mods = ~ age + I(age^2) + I(age^3)
              + TaskCode_1 + TaskCode_2 + TaskCode_3 + HandnessCode_1 + HandnessCode_2 + ContrastCode_1 + ContrastCode_2 + ErrorTrialCode_1, 
              data = data2,
              control=list(stepadj=step))
        res.log <-
          rma(beta,
              var,
              mods = ~ I(log(age)) + I((log(age))^2)
              + TaskCode_1 + TaskCode_2 + TaskCode_3 + HandnessCode_1 + HandnessCode_2 + ContrastCode_1 + ContrastCode_2 + ErrorTrialCode_1, 
              data = data2,
              control=list(stepadj=step))
        res.sqrt <-
          rma(beta,
              var,
              mods = ~ I(sqrt(age)) + age
              + TaskCode_1 + TaskCode_2 + TaskCode_3 + HandnessCode_1 + HandnessCode_2 + ContrastCode_1 + ContrastCode_2 + ErrorTrialCode_1, 
              data = data2,
              control=list(stepadj=step))
        res.lin <-
          rma(beta,
              var,
              mods = ~ age
              + TaskCode_1 + TaskCode_2 + TaskCode_3 + HandnessCode_1 + HandnessCode_2 + ContrastCode_1 + ContrastCode_2 + ErrorTrialCode_1, 
              data = data2,
              control=list(stepadj=step))
      } else {
        res.qua <- 
          rma(beta, 
              var, 
              mods = ~ age + I(age^2)
              + TaskCode_1 + TaskCode_2 + TaskCode_3 + HandnessCode_1 + HandnessCode_2 + ContrastCode_1 + ContrastCode_2 + ErrorTrialCode_1 + ErrorTrialCode_2,
              data = data2,
              control=list(stepadj=step))
        res.cub <-
          rma(beta,
              var,
              mods = ~ age + I(age^2) + I(age^3)
              + TaskCode_1 + TaskCode_2 + TaskCode_3 + HandnessCode_1 + HandnessCode_2 + ContrastCode_1 + ContrastCode_2 + ErrorTrialCode_1 + ErrorTrialCode_2, 
              data = data2,
              control=list(stepadj=step))
        res.log <-
          rma(beta,
              var,
              mods = ~ I(log(age)) + I((log(age))^2)
              + TaskCode_1 + TaskCode_2 + TaskCode_3 + HandnessCode_1 + HandnessCode_2 + ContrastCode_1 + ContrastCode_2 + ErrorTrialCode_1 + ErrorTrialCode_2, 
              data = data2,
              control=list(stepadj=step))
        res.sqrt <-
          rma(beta,
              var,
              mods = ~ I(sqrt(age)) + age
              + TaskCode_1 + TaskCode_2 + TaskCode_3 + HandnessCode_1 + HandnessCode_2 + ContrastCode_1 + ContrastCode_2 + ErrorTrialCode_1 + ErrorTrialCode_2, 
              data = data2,
              control=list(stepadj=step))
        res.lin <-
          rma(beta,
              var,
              mods = ~ age
              + TaskCode_1 + TaskCode_2 + TaskCode_3 + HandnessCode_1 + HandnessCode_2 + ContrastCode_1 + ContrastCode_2 + ErrorTrialCode_1 + ErrorTrialCode_2, 
              data = data2,
              control=list(stepadj=step))
      }
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
      # later we can calculate the p of mean QMs
      # p_value <- pchisq(res.qua$QM, res.qua$QMdf, lower.tail = FALSE)
    # }
    
    # if (iperm == nperm+1) {
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
      
      QMs[i, 1] = res.qua$QM
      QMs[i, 2] = res.cub$QM
      QMs[i, 3] = res.log$QM
      QMs[i, 4] = res.sqrt$QM
      QMs[i, 5] = res.lin$QM
      
      mod.comparison <- AICc(res.qua, res.cub, res.log, res.sqrt, res.lin)
      modelweights[i,] <- Weights(mod.comparison)
    # } 
  }
# }

# # check results
p_gam_fdr <- p.adjust(ps_gam, method = 'fdr', n = length(ps_gam))

p_values <- p_zs <- zval_meanage <- zval_randage <- array(,dim = c(nfile,5))
F_meanage <- F_randage <- array(,dim = c(nfile,1))
for (i in 1:nfile) {
  
  # calculate the modulator effect (including age and others) on average across the nperm values
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
  
  #calculate the age modulator effect in each model on average across the nperm values
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

# ## plot the distribution, mean of distribution and the value based on mean age
# source('H:/meta_development/Results_R1/SDM/scripts/plot_distribution.R')
# savedir = shortPathName("H:/meta_development/Results_R1/SDM/mean/analysis_Adult-ChildOld/plot/testmeanage/")
# if (!dir.exists(savedir)) {
#   dir.create(savedir, recursive = TRUE)
# }
# for (i in 1:nfile) {
#   filename = str_match(filenames[i], "extract_\\s*(.*?)\\s*.txt")
#   tiff(file = paste(savedir, filename[, 2], "_gamF.tiff", sep = ''),
#        width = 6,height = 5,units = "in",res = 300)
#   plot_distribution(Fs[i,1:nperm],F_meanage[i,1])
#   dev.off()
#   tiff(file = paste(savedir, filename[, 2], "_sqrtZ.tiff", sep = ''),
#        width = 6,height = 5,units = "in",res = 300)
#   plot_distribution(zss[i,4,1:nperm], zval_meanage[i,4]) #sqrt, due to it is the only one significant in 2 regions
#   dev.off()
#   tiff(file = paste(savedir, filename[, 2], "_linZ.tiff", sep = ''),
#        width = 6,height = 5,units = "in",res = 300)
#   plot_distribution(zss[i,5,1:nperm], zval_meanage[i,5]) #lin, due to it is the simplest
#   dev.off()
# }


# rs <- c()
# for (imodel in 1:4) {
#   rs[imodel] <- cor(zval_meanage[,imodel],zval_randage[,imodel])
# }
# 
# # apply(p_values, 1, min)
# apply(p_zs, 1, min)
# # check if the mean age results represent the mean of random age results
# 
# }
# }
# pss = cbind(ps00,ps01,ps10,ps11)

options(scipen=999)
# pss

minAIC <- minBIC <- maxLogLik <- c()
for (i in 1:nfile) {
  minAIC[i] = which.min(AICs[i,])
  minBIC[i] = which.min(AICs[i,])
  maxLogLik[i] = which.max(logliks[i,])
}


modelweights2 <- format(round(modelweights, 3), nsmall = 3)
write.matrix(modelweights2, file = paste(savedir,"modelweights.csv"), sep = ',')


###  plot the figures  ###
xs <- seq(8, 75, length=500)

minmax = maxLogLik
peak_rmas <- peak_gams <- c()
for (i in 1:nfile) {
  # the covariate matrix
  # data2 <- as.data.frame(reses_qua[[i]]$X)
  data2 <- data2s[[i]]
  
  mtx_cov_xs <- matrix(0,nrow = length(xs), ncol = length(reses_qua[[i]]$beta)-3)
  if (minmax[i] == 1) {
    modelbest = reses_qua[[i]]
    title = "QUADRATIC"
    sav <-
      predict(modelbest, newmods = unname(cbind(poly(
        xs, degree = 2, raw = TRUE
      ),mtx_cov_xs)))
    peak_rma = -0.5*modelbest$b[2,1]/modelbest$b[3,1]
    if (modelbest$pval[3] > 0.05) {
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
    if (modelbest$pval[4] > 0.05) {
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
    if (modelbest$pval[3] > 0.05) {
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
    if (modelbest$pval[2] > 0.05) {
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
    # if (modelbest$pval[2] > 0.05 && p_gam_fdr[i] < 0.05) {
    #   sav$pred <- NaN
    # }
  }
  peak_rmas[i] = peak_rma

  
  
  # ylim0 <- c()
  # maxvar = max(data2$variance)
  # ylim0[1] = min(floor((data2$beta-maxvar)*100)/100)
  # ylim0[2] = max(ceiling((data2$beta+maxvar)*100)/100)
  ylim0 = c(-0.1,0.7)
  
  modpred <- data.frame(xs,sav$pred,sav$se)
  modpred$upper <- modpred$sav.pred + modpred$sav.se
  modpred$lower <- modpred$sav.pred - modpred$sav.se
  
  
  ### gam prediction
  b <- bs[[i]]
  xs0 <- b$model$age
  preds0 <- predict.gam(b,se.fit = TRUE)
  
  # mtx_b <- data.frame(model.matrix(b)) # if outliers removed the data points are different with data2, so here use the model matrix to get the covariate values
  
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
    ErrorTrialCode_2 = 0
  ),se.fit = TRUE)
  
  peak_gam <- xs[which.max(preds$fit)]
  if (peak_gam > 60 || p_gam_fdr[i] > 0.05) {
    peak_gam <- NaN
  }
  peak_gams[i] = round(peak_gam,1)
  
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
    file = paste(savedir, filename[, 2], "_sameylim_0.7.tiff", sep = ''),
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
  
  # par(new=TRUE)
  # # windows()
  # lines(modpred$xs, modpred$sav.pred, type = 'l', lwd = 3, col = 'red')  # or type = 'n' if you just want to set up axes
  # transparent_red <- rgb(1, 0, 0, alpha = 0.2)  # For example, gray with 50% transparency
  # polygon(c(modpred$xs, rev(modpred$xs)), c(modpred$upper, rev(modpred$lower)), col = transparent_red, border = NA)
  # lines(modpred$xs, modpred$sav.pred)  # Add main line after to ensure it's on top
  
  
  par(new=TRUE)
  transparent_color <- rgb(0, 0, 1, alpha = 0.2)  # For example, gray with 50% transparency
  lines(gampred$xs, gampred$preds.fit, type = 'l', lwd = 3, col = 'blue')  # or type = 'n' if you just want to set up axes
  polygon(c(gampred$xs, rev(gampred$xs)), c(gampred$upper, rev(gampred$lower)), col = transparent_color, border = NA)
  lines(gampred$xs, gampred$preds.fit)  # Add main line after to ensure it's on top
  # lines(xs0, preds0$fit, type = 'l', lwd = 3, col = 'blue')  # or type = 'n' if you just want to set up axes
  # polygon(c(xs0, rev(xs0)), c(preds0$fit + preds0$se.fit, rev(preds0$fit - preds0$se.fit)), col = transparent_color, border = NA)
  # lines(xs0, preds0$fit)  # Add main line after to ensure it's on top
  
  # par(new=TRUE)
  # lines(c(peak_rma,peak_rma),ylim0,lty = 2,col = 'red',lwd = 2)
  if (i > 3) {
    par(new=TRUE)
    lines(c(peak_gam,peak_gam),ylim0,lty = 2,col = 'blue',lwd = 2)    
  }
  dev.off()
}



## test the significance of peaks with the Simonsohn (2018)'s two-line approach. 
# Simonsohn, U. (2018) Two Lines: A Valid Alternative to the Invalid Testing of U-Shaped Relationships With Quadratic Regressions. Advances in Methods and Practices in Psychological Science, 1(4), 538-555.
psb <- array(,dim = c(length(bs),2,1)) #9x2(L,R)x(peakGAM #,peakRMA)
for (ipeak in 1) {
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
      ErrorTrialCode_2 = 0
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
      bX <- gam(beta~age + TaskCode_1 + TaskCode_2 + TaskCode_3 + HandnessCode_1 + HandnessCode_2 + ContrastCode_1 + ContrastCode_2 + ErrorTrialCode_1 + ErrorTrialCode_2,data=databX, weights = weight, method = 'REML')
      bXsum <- summary(bX)
      psb[i,iana,ipeak] <- bXsum$p.pv[2]
    }
  }
}
psb_fdr <- array(dim = dim(psb))
for (ipeak in 1) {
  psb_fdr[,,ipeak] <- t(apply(psb[,,ipeak]/2, 1, function(x) p.adjust(x, method = "fdr")))
}
