rm(list=ls())
library(metafor)
library(performance)
library(car)
library(visreg)
library(MuMIn)
library(mgcv)
library(MASS)
library("readxl")
library(fastDummies)
library(lme4)
library(gammit)
source('H:/meta_development/Results_R1/SDM/scripts/naturalSortFunctions.R')
source('H:/meta_development/Results_R1/SDM/scripts/plot_distribution.R')

## load data and calculate the averaged betas and variances
datadir = shortPathName("H:/meta_development/Results_R1/SDM/mean/analysis_Adult-ChildOld/") 
sdmdir = 'H:/meta_development/Results_R1/SDM/mean/' #sdm_good

filenames = list.files(path = paste(datadir,'analysis_MyLinearModel/extracts',sep = ''),
                       pattern = "^multivoxel_extract_mean_good_z_voxelCorrected_p_0.00100_10_blob.*\\.txt$")
filenames = natural_sort(filenames,'blob','.txt')

nperm = 1000

savedir = shortPathName(paste0("H:/meta_development/Results_R1/SDM/mean/analysis_Adult-ChildOld/plot/testmeanage/",nperm,"/"))
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


# Fs_gam <- ps_gam <- data2s <- bs <- reses_qua <- reses_cub <- reses_log <- reses_sqrt <- reses_revlog <- reses_revsqrt <- reses_lin <- c()
# Fs <- ps <- df1s <- df2s <- array(,dim = c(nfile,nperm+1))
# QMss <- zss <- array(, dim = c(nfile, 5, nperm+1)) #ROI, 4, #perm
rmeans <- rmins <- c()
for (i in 1:nfile) {
# for (i in c(1:6,8)) {
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
  
  #screen data through sdm table
  idx_incl <- which(data$study %in% sdmtable$study)
  data_good <- data[idx_incl,] 
  
  beta = aggregate(data_good$betas, list(data_good$study), FUN = mean)
  variance = aggregate(data_good$variances, list(data_good$study), FUN = mean)
  data2 <- data.frame(beta, variance$x)
  
  
  colnames(data2)[1] <- 'study'
  colnames(data2)[2] <- 'beta'
  colnames(data2)[3] <- 'variance'
  
  
  # delete outliers
  outup = mean(data2$beta) + 3*sd(data2$beta)
  outdown = mean(data2$beta) - 3*sd(data2$beta)
  data2 = data2[-which(data2$beta>outup | data2$beta < outdown),]
  
  
  # add variables
  #note data2$varianceAgerange is standard deviation, whereas data2$variance is deviation 
  for (istudy in 1:length(data2$study)) {
    # print(istudy)
    idx <- which(supptable$AuthoYear == data2$study[istudy])
    data2$varianceAgerange[istudy] <- as.numeric(supptable$Variance[idx])
    data2$nsub[istudy] <- as.numeric(supptable$nsub[idx])
    data2$upbound[istudy] <- as.numeric(supptable$upbound[idx])
    data2$lowbound[istudy] <- as.numeric(supptable$lowbound[idx])
    
    #add age from sdmtable
    idx2 <- which(sdmtable$study == data2$study[istudy])
    data2$age[istudy] <- sdmtable$avgAge[idx2]
    data2$age0[istudy] <- sdmtable$avgAge[idx2]
  }
  
  
  # sort data2 with the age
  data2 <- data2[order(data2$age), ]

  ## select analysis type
  data2$var = data2$variance
  data2$weight = 1/data2$variance
  
  nstudy = length(data2$study)
  
  # simulate the ages and effect sizes for each subject from each study based on noramal distribution of (age0, varianceAgerange) and (beta, variance) for 1000 times
  preds_sim <- preds_sim_se <- array(,dim = c(nstudy,nperm))
  
  for (iperm in 1:(nperm+1)) {
    print(iperm)
    if (iperm < (nperm+1)) {
      studies <- ages <- betas <- c()
      for (istudy in 1:length(data2$study)) {
        studies <- c(studies,rep(istudy,data2$nsub[istudy],replace=TRUE))
        
        if (! is.nan(data2$upbound[istudy])) {
          upbound <- data2$upbound[istudy]
          lowbound <- data2$lowbound[istudy]
        } else {
          upbound <- data2$age0[istudy] + 3*data2$varianceAgerange[istudy]
          lowbound <- data2$age0[istudy] - 3*data2$varianceAgerange[istudy]
          if (data2$age0[istudy] < 18) {
            upbound <- min(c(upbound,18))
            lowbound <- max(c(lowbound,0))
          } else if (data2$age0[istudy] < 52) {
            upbound <- min(c(upbound,60))
            lowbound <- max(c(lowbound,18))
          } else {
            upbound <- min(c(upbound,100))
            lowbound <- max(c(lowbound,60))
          }
        }
        upbound_beta <- data2$beta[istudy] + 3*sqrt(data2$variance[istudy])
        lowbound_beta <- data2$beta[istudy] - 3*sqrt(data2$variance[istudy])
        
        simages <- rnorm(data2$nsub[istudy], mean = data2$age0[istudy], sd = data2$varianceAgerange[istudy])
        bounded_ages <- pmax(pmin(simages, upbound), lowbound)
        adjusted_simages <- data2$age0[istudy] + (bounded_ages - mean(bounded_ages))
        ages <- c(ages,adjusted_simages)
        
        simbetas <- rnorm(data2$nsub[istudy], mean = data2$beta[istudy], sd = sqrt(data2$variance[istudy]))
        bounded_betass <- pmax(pmin(simbetas, upbound_beta), lowbound_beta)
        adjusted_simbetas <- data2$beta[istudy] + (bounded_betass - mean(bounded_betass))
        betas <- c(betas,adjusted_simbetas)
      }
      data_sim <- data.frame(studies,ages,betas)
      # data_sim$sqrtages <- sqrt(data_sim$ages)
      
      # models
      b_sim <- gam(betas~s(ages) + s(studies, bs = "re") + s(ages, studies, bs = "re"),data=data_sim, method = 'REML')
      sumb <- summary(b_sim)

      pred_sim <- predict(b_sim, newdata = data.frame(ages = data2$age0, studies = c(1:nstudy)),se.fit = TRUE)
      preds_sim[,iperm] <- pred_sim$fit
      preds_sim_se[,iperm] <- pred_sim$se.fit
      # pred_sim <- predict_gamm(b_sim,newdata=data.frame(ages = data2$age0, studies = c(1:nstudy)),re_form=NA,se.fit = TRUE)
      # preds_sim[,iperm] <- pred_sim$prediction.fit
      # preds_sim_se[,iperm] <- pred_sim$prediction.se.fit
    } else if (iperm == nperm+1) {
      b <- gam(beta~s(age0),data=data2, weights = weight, method = 'REML') # meta-regression
      preds <- predict(b, se.fit = TRUE)
      preds$upper <- preds$fit + preds$se.fit
      preds$lower <- preds$fit - preds$se.fit
    }
  }
  
  rs <- c()
  for (iperm in 1:nperm) {
    r <- cor(preds_sim[,iperm],preds$fit)
    rs <- c(rs,r)
  }
  rmin <- min(rs)
  rmins <- c(rmins,rmin)
  rmean <- cor(rowMeans(preds_sim),preds$fit)
  rmeans <- c(rmeans,rmean)
  

  filename <- gsub(paste0(".*", "extract_", "(.*?)", ".txt", ".*"), "\\1", filenames[i])
  tiff(
    file = paste(savedir, filename, "_simdata",nperm, ".tiff", sep = ''),
    width = 6,
    height = 5,
    units = "in",
    res = 300
  )
  
  simpreds <- c()
  simpreds$mean <- rowMeans(preds_sim,na.rm = TRUE)
  simpreds$se <- rowMeans(preds_sim_se,na.rm = TRUE)
  simpreds$upper <- simpreds$mean + simpreds$se
  simpreds$lower <- simpreds$mean - simpreds$se
  
  ylim0 = c(min(c(simpreds$lower,preds$lower))-0.01,max(c(simpreds$upper,preds$upper))+0.01)
  ylim0 = c(0,0.30)
  
  # windows()
  par(cex.lab=2, cex.axis=2, mar=c(5,5,4,2))
  plot(data2$age0,simpreds$mean,col='red',type = "n",bty='l', ylim = ylim0, xlab = "Age (years)",ylab = "Effect size")
  lines(data2$age0,simpreds$mean, type = 'l', lwd = 3, col = 'red')  # or type = 'n' if you just want to set up axes
  polygon(c(data2$age0, rev(data2$age0)), c(simpreds$upper, rev(simpreds$lower)), col = rgb(1, 0, 0, alpha = 0.3), border = NA)
  par(new=TRUE)
  box(lty=0)
  lines(data2$age0,preds$fit, type = 'l', lwd = 3, col = 'blue',ylim=ylim0)  # or type = 'n' if you just want to set up axes
  polygon(c(data2$age0, rev(data2$age0)), c(preds$upper, rev(preds$lower)), col = rgb(0, 0, 1, alpha = 0.3), border = NA)
  # if (i == 1) {
  #   legend("topright", legend=c("simulation", "meta-regression"), col=c("red", "blue"), lty=1, lwd=3)
  # }
  dev.off()
}

r_table <- data.frame(rmins,rmeans)
write.matrix(r_table, file = paste(savedir,"rs.csv"), sep = ',')