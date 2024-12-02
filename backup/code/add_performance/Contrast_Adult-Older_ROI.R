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
datadir = '../data/'
sdmdir = '../data/'

filenames = list.files(path = paste(datadir,'extracted_data_contrastanalysis',sep = ''),
                       pattern = "^multivoxel_extract_Adult_minus_Older_good_z_voxelCorrected_p_0.00100_10_pblob.*\\.txt$")

filenames = natural_sort(filenames,'pblob','.txt')

savedir = shortPathName("../plot/blob_adult-older_good+design+SRC/")
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


# Fs_gam <- ps_gam <- data2s <- bs <- reses_qua <- reses_cub <- reses_log <- reses_sqrt <- reses_revlog <- reses_revsqrt <- reses_lin <- c()
reses_lin <- c()
Fs <- ps <- df1s <- df2s <- array(,dim = c(nfile))
zss <- pss <- array(, dim = c(nfile, 1)) #ROI, 4, #perm
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
    # idx2 <- which(sdmtable$study == data2$study[istudy])
    # data2$SRC[istudy] <- sdmtable$performance[idx2]
    
    #add age from sdmtable
    idx2 <- which(sdmtable$study == data2$study[istudy])
    data2$age[istudy] <- sdmtable$avgAge[idx2]
    data2$age0[istudy] <- sdmtable$avgAge[idx2]
    data2$filter_AdultOlder[istudy] <- sdmtable$filter_AdultOlder[idx2]
    data2$Adult_Older[istudy] <- sdmtable$Adult_Older[idx2]
  }
  
  # idx_adultolder <- which(data2$filter_AdultOlder == 1 & !is.na(data2$SRC) &data2$performance<999)
  data2_bk <- data2
  idx_adultolder <- which(data2$filter_AdultOlder == 1 & !is.na(data2$SRC))
  data2 <- data2[idx_adultolder,]
  
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
  
  # # delete outliers
  # outup = mean(data2$beta) + 3*sd(data2$beta)
  # outdown = mean(data2$beta) - 3*sd(data2$beta)
  # data2 = data2[-which(data2$beta>outup | data2$beta < outdown),]
  # data2s[[i]] = data2
  
  # sort data2 with the age
  data2 <- data2[order(data2$age), ]
  
  ## select analysis type
  data2$var = data2$variance
  data2$weight = 1/data2$variance
  
  data2$age <- data2$age0
  
  
  step = 0.1 #avoid not converging
  
  zvals <- pvals <- bvals <- c()
  
  data2$age <- data2$age0
  if (sum(data2$ErrorTrialCode_2) == 0) {
    res.lin <-
      rma(beta,
          var,
          mods = ~ Adult_Older
          + TaskCode_1 + TaskCode_2 + TaskCode_3 + HandnessCode_1 + HandnessCode_2 + ContrastCode_1 + ContrastCode_2 + ErrorTrialCode_1 + DesignCode_1 + SRC, 
          data = data2,
          control=list(stepadj=step))
  } else {
    res.lin <-
      rma(beta,
          var,
          mods = ~ age
          + TaskCode_1 + TaskCode_2 + TaskCode_3 + HandnessCode_1 + HandnessCode_2 + ContrastCode_1 + ContrastCode_2 + ErrorTrialCode_1 + ErrorTrialCode_2 + DesignCode_1 + SRC, 
          data = data2,
          control=list(stepadj=step))
    }
    zss[i,1] <- res.lin$zval[2]
    pss[i,1] <- res.lin$pval[2]
    
    reses_lin[[i]] = res.lin
}

# zss = c(3.29 2.57 2.81 2.278 2.51 2.87 2.38 2.36)
# pss = c(0.0009921706 0.01015866 0.004966197 0.02282496 0.01203655 0.004071789 0.01739482 0.01831173)