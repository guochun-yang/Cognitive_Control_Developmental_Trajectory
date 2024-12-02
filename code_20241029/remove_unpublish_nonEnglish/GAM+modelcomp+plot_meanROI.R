rm(list=ls())

basedir = 'C:/Users/Guochun Yang/OneDrive - University of Iowa/WithLi/Paper/code' ######
scriptdir = paste(basedir,'code_20241029',sep = '/') #####
funcdir = paste(scriptdir,"func",sep = '/')
funcfiles <- list.files(funcdir, pattern = "\\.R$", full.names = TRUE)
lapply(funcfiles, source)

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


## load data and calculate the averaged betas and variances
datadir = shortPathName(paste(basedir, "/meta_data/analysis_Adult-ChildOld/",sep = '')) 
sdmdir = paste(basedir, "meta_data",sep = '/') 

filenames = list.files(path = datadir,
                       pattern = "^multivoxel_extract_mean_good_z_voxelCorrected_p_0.00100_10_blob.*\\.txt$")
filenames = natural_sort(filenames,'blob','.txt')

savedir = shortPathName(paste(datadir,"plot/blob_mean_rmunpublish/",sep='/'))

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
  paste(scriptdir, 'covariates.xlsx', sep = '/'),
  sheet = 'covariates'
)


Fs_gam <- ps_gam <- data2s <- bs <- reses_lin <- c()
Fs <- ps <- df1s <- df2s <- array(,dim = c(nfile))
for (i in 1:nfile) {
  data = read.table(
    paste(datadir, filenames[i], sep = '/'),
    header = FALSE,
    sep = '',
    dec = '.',
    skip = 13
  )
  
  data$study = sub("_I000.*", "", data$V1)
  data$betas = data$V2
  
  data$variances = data$V3
  
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
    idx <- which(supptable$AuthoYear == data2$study[istudy])
    data2$varianceAgerange[istudy] <- as.numeric(supptable$Variance[idx])
    data2$TaskCode[istudy] <- supptable$TaskCode[idx]
    data2$HandnessCode[istudy] <- supptable$HandnessCode[idx]
    data2$ContrastCode[istudy] <- supptable$ContrastCode[idx]
    data2$ErrorTrialCode[istudy] <- supptable$ErrorTrialCode[idx]
    data2$DesignCode[istudy] <- supptable$DesignCode[idx]
    data2$pubLanguage[istudy] <- supptable$pubLanguage[idx]
    data2$SRC[istudy] <- supptable$SRC_imp_median[idx]
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
  }
  
  # delete outliers
  outup = mean(data2$beta) + 3*sd(data2$beta)
  outdown = mean(data2$beta) - 3*sd(data2$beta)
  data2 = data2[-which(data2$beta>outup | data2$beta < outdown),]
  data2s[[i]] = data2
  
  # remove unpublished/non-English ones
  data2 = data2[which(data2$pubLanguage==1),]
  
  # sort data2 with the age
  data2 <- data2[order(data2$age), ]
  
  ## select analysis type
  data2$var = data2$variance
  data2$weight = 1/data2$variance
  
  # randomly choose the age from the noramal distribution of (age, varirance) for 1000 times
  options(digits = 7) #7 is default

  data2$age <- data2$age0

  # models
  b <- gam(beta~s(age) + TaskCode_1 + TaskCode_2 + TaskCode_3 + HandnessCode_1 + HandnessCode_2 + ContrastCode_1 + ContrastCode_2 + ErrorTrialCode_1 + ErrorTrialCode_2 + DesignCode_1 + SRC + SRCcode,data=data2, weights = weight, method = 'REML')
  
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
  
  bs[[i]] = b
  ps_gam[i] <- p
  Fs_gam[i] <- F
  
  step = 0.1 #avoid not converging
  data2$age <- data2$age0
  if (sum(data2$ErrorTrialCode_2) == 0) {
    res.lin <-
      rma(beta,
          var,
          mods = ~ age
          + TaskCode_1 + TaskCode_2 + TaskCode_3 + HandnessCode_1 + HandnessCode_2 + ContrastCode_1 + ContrastCode_2 + ErrorTrialCode_1 + DesignCode_1 + SRC + SRCcode,
          data = data2,
          control=list(stepadj=step))
  } else {
    res.lin <-
      rma(beta,
          var,
          mods = ~ age
          + TaskCode_1 + TaskCode_2 + TaskCode_3 + HandnessCode_1 + HandnessCode_2 + ContrastCode_1 + ContrastCode_2 + ErrorTrialCode_1 + ErrorTrialCode_2 + DesignCode_1 + SRC + SRCcode,
          data = data2,
          control=list(stepadj=step))
  }
  reses_lin[[i]] = res.lin
}

# # check results
p_gam_fdr <- p.adjust(ps_gam, method = 'fdr', n = length(ps_gam))


options(scipen=999)

###  plot the figures  ###
xs <- seq(8, 75, length=500)

# minmax = maxLogLik
peak_gams <- c()
for (i in 1:nfile) {
  # the covariate matrix
  data2 <- data2s[[i]]
  
  ylim0 = c(-0.1,0.7)
  
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
    ErrorTrialCode_2 = 0,
    DesignCode_1 = 0,
    SRC = 0,
    SRCcode = 0
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
    reses_lin[[i]],
    mod = 2,
    lcol = 'red',
    # pred = FALSE,
    ci=FALSE,
    pred = FALSE,
    xvals = xs,
    las = 1,
    digits = 1,
    bty = "l",
    psize = .20 / sqrt(reses_lin[[i]]$vi),
    xlab = "Age (years)",ylab = "Effect size",
    xlim = c(8, 75),
    ylim = c(ylim0[1],ylim0[2])
  )
  lim0 <- par("usr")
  
  par(new=TRUE)
  transparent_color <- rgb(0, 0, 1, alpha = 0.2)  # For example, gray with 50% transparency
  lines(gampred$xs, gampred$preds.fit, type = 'l', lwd = 3, col = 'blue')  # or type = 'n' if you just want to set up axes
  polygon(c(gampred$xs, rev(gampred$xs)), c(gampred$upper, rev(gampred$lower)), col = transparent_color, border = NA)
  lines(gampred$xs, gampred$preds.fit)  # Add main line after to ensure it's on top
  if (i > 2) {
    par(new=TRUE)
    lines(c(peak_gam,peak_gam),ylim0,lty = 2,col = 'blue',lwd = 2)    
  }
  dev.off()
}


# 
# ## test the significance of peaks with the Simonsohn (2018)'s two-line approach. 
# # Simonsohn, U. (2018) Two Lines: A Valid Alternative to the Invalid Testing of U-Shaped Relationships With Quadratic Regressions. Advances in Methods and Practices in Psychological Science, 1(4), 538-555.
# psb <- array(,dim = c(length(bs),2,1)) #9x2(L,R)x(peakGAM #,peakRMA)
# for (ipeak in 1:2) {
#   for (i in 1:length(bs)) {
#     b <- bs[[i]]
#     bsum <- summary(b)
#     if (bsum$s.pv > 0.05) {psb[i,iana,ipeak] = NA; next}
#     datab <- b$model
#     datab$var <- 1/datab$`(weights)`
#     datab$weight <- datab$`(weights)`
#     datab$sqrtage <- sqrt(datab$age)
#     predictb <- predict(b, newdata = data.frame(
#       age = datab$age,
#       TaskCode_1 = 0,
#       TaskCode_2 = 0,
#       TaskCode_3 = 0,
#       HandnessCode_1 = 0,
#       HandnessCode_2 = 0,
#       ContrastCode_1 = 0,
#       ContrastCode_2 = 0,
#       ErrorTrialCode_1 = 0,
#       ErrorTrialCode_2 = 0,
#       DesignCode_1 = 0,
#       SRC = 0,
#       SRCcode = 0
#     ))
#     if (ipeak == 1) {
#       idx <- which(predictb==max(predictb))
#     } else {
#       peak <- peak_rmas[i]
#       idx <- which.min(abs(datab$age-peak))
#     }
#     
#     if (length(idx) == 1) {
#       databL <- datab[1:idx,]
#       databR <- datab[idx:length(predictb),]
#     } else {
#       databL <- datab[1:idx[1],]
#       databR <- datab[idx[-1]:length(predictb),]
#     }
#     
#     # fit linear model and see if the slope is significant
#     step = 0.1
#     for (iana in 1:2) { #left and right
#       if (iana == 1) {databX <- databL} else {databX <- databR}
#       bX <- gam(beta~age + TaskCode_1 + TaskCode_2 + TaskCode_3 + HandnessCode_1 + HandnessCode_2 + ContrastCode_1 + ContrastCode_2 + ErrorTrialCode_1 + ErrorTrialCode_2 + DesignCode_1 + SRC + SRCcode,data=databX, weights = weight, method = 'REML')
#       bXsum <- summary(bX)
#       psb[i,iana,ipeak] <- bXsum$p.pv[2]
#     }
#   }
# }
# psb_fdr <- array(dim = dim(psb))
# for (ipeak in 1) {
#   psb_fdr[,,ipeak] <- t(apply(psb[,,ipeak]/2, 1, function(x) p.adjust(x, method = "fdr")))
# }
