rm(list=ls())
library(metafor)
library(stringr)
library(performance)
library(car)
library(visreg)
library(MuMIn)
library(mgcv)
library(MASS)

## load data and calculate the averaged betas and variances
datadir = shortPathName("./data_share/")
savedir = shortPathName("UNDEFINED")
  
filenames = list.files(path = paste(datadir, 'extracted_data_contrastanalysis', sep = '/'),
                       pattern = "^extract_M.*\\.txt$")
nfile = length(filenames)

AICs = matrix(, nrow = nfile, ncol = 4)
BICs = matrix(, nrow = nfile, ncol = 4)
QMps = matrix(, nrow = nfile, ncol = 4)
QMs = matrix(, nrow = nfile, ncol = 4)
modelweights = matrix(, nrow = nfile, ncol = 4)

models = list()

Fs_gam <- ps_gam <- data3s <- bs <- reses_qua <- reses_cub <- reses_log <- reses_sqrt <- reses_revlog <- reses_revsqrt <- c()

for (i in 1:nfile) {
  data = read.table(
    paste(datadir, 'extracted_data_contrastanalysis', filenames[i], sep = '/'),
    header = FALSE,
    sep = '',
    dec = '.',
    skip = 13
  )
  sdmtable = read.table(
    paste(datadir, '/sdm_table.txt', sep = ''),
    header = TRUE,
    sep = '',
    dec = '.'
  )
  
  data$study = sub("_I000.*", "", data$V1)
  data$betas = data$V2
  
  data$variances = data$V3
  data$weight = 1/data$variances
  
  beta = aggregate(data$betas, list(data$study), FUN = mean)
  variance = aggregate(data$variances, list(data$study), FUN = mean)
  data2 <- data.frame(beta, variance)
  
  colnames(data2)[2] <- 'beta'
  colnames(data2)[3] <- 'study'
  colnames(data2)[4] <- 'variance'
  
  data3 <- data2[order(data2$study), ]
  data3 <- sdmtable[order(sdmtable$study), ]
  
  data3$study2  = data2$study
  data3$age = data3$avgAge
  data3$beta = data2$beta
  data3$variance = data2$variance
  data3$weight = 1/data3$variance
  agefull = max(data3$age) + min(data3$age)
  
  # delete outliers
  outup = mean(data3$beta) + 3*sd(data3$beta)
  outdown = mean(data3$beta) - 3*sd(data3$beta)
  data3 = data3[-which(data3$beta>outup | data3$beta < outdown),]
  data3s[[i]] = data3
  
  # models
  b <- gam(beta~s(age),data=data3, weights = weight, method = 'REML')
  bs[[i]] = b
  x<-summary(b)
  ps_gam[i] <- x$s.table[4]
  Fs_gam[i] <- x$s.table[3]
  
  res.qua <- rma(beta, variance, mods = ~ age + I(age^2), data = data3)
  res.cub <-
    rma(beta,
        variance,
        mods = ~ age + I(age^2) + I(age^3),
        data = data3)
  res.log <-
    rma(beta,
        variance,
        mods = ~ I(log(age)) + I((log(age))^2),
        data = data3)
  res.sqrt <-
    rma(beta,
        variance,
        mods = ~ I(sqrt(age)) + age,
        data = data3)
  
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
  
  mod.comparison <- AICc(res.qua, res.cub, res.log, res.sqrt)
  modelweights[i,] <- Weights(mod.comparison)
}

minAIC <- minBIC <- c()
for (i in 1:nfile) {
  minAIC[i] = which.min(AICs[i,])
  minBIC[i] = which.min(AICs[i,])
}



modelweights2 <- format(round(modelweights, 3), nsmall = 3)
write.matrix(modelweights2, file = paste(savedir,"modelweights.csv"), sep = ',')


### draw plot
xs <- seq(8, 75, length=500)

peak_rmas <- peak_gams <- c()
for (i in 1:15) {
  if (minAIC[i] == 1) {
    modelbest = reses_qua[[i]]
    title = "QUADRATIC"
    sav <-
      predict(modelbest, newmods = unname(poly(
        xs, degree = 2, raw = TRUE
      )))
    peak_rma = -0.5*modelbest$b[2,1]/modelbest$b[3,1]
  } else if (minAIC[i] == 3) {
    modelbest = reses_sqrt[[i]]
    title = "LOG"
    sav <-
      predict(modelbest, newmods = unname(poly(
        log(xs), degree = 2, raw = TRUE
      )))
  } else if (minAIC[i] == 4) {
    modelbest = reses_sqrt[[i]]
    title = "SQRT"
    sav <-
      predict(modelbest, newmods = unname(poly(
        sqrt(xs), degree = 2, raw = TRUE
      )))
    peak_rma = 0.25*(modelbest$b[2,1]/modelbest$b[3,1])^2
  }
  peak_rmas[i] = peak_rma

  filename = str_match(filenames[i], "_p_0.00100_10_\\s*(.*?)\\s*.txt")
  tiff(
    file = paste(savedir, filename[, 2], "_sameylim.tiff", sep = ''),
    width = 6,
    height = 5,
    units = "in",
    res = 300
  )

  data3 <- data3s[[i]]
  ylim0 <- c()
  maxvar = max(data3$variance)
  # ylim0[1] = min(floor((data3$beta-maxvar)*100)/100)
  # ylim0[2] = max(ceiling((data3$beta+maxvar)*100)/100)
  ylim0 = c(-0.2,1.0)
  
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
    pred = sav,
    xvals = xs,
    las = 1,
    digits = 1,
    bty = "l",
    psize = .20 / sqrt(modelbest$vi),
    xlab = "Age",ylab = "effect size",
    xlim = c(8, 75),
    ylim = c(ylim0[1],ylim0[2])
  )
  # lim0 <- par("usr")

  par(new=TRUE)
  # add visreg
  b <- bs[[i]]
  # b <- gam(beta~s(age),data=data3, weights = weight, method = 'REML')
  preds <- predict.gam(b, newdata = data.frame(age = xs))
  peak_gam <- xs[which.max(preds)]
  peak_gams[i] = round(peak_gam,1)
  # visreg(b,axes=FALSE, pch = NA, xlab='',ylab='',yaxs='i',ylim = c(lim0[3],lim0[4]))
  visreg(b,"age",
         overlay = TRUE,partial = FALSE, rug = FALSE,
         line=list(lty=1, col="blue"), alpha=0,
         points=list(cex=1, pch=26),
         axes=FALSE, xlab='',ylab='',yaxs='i',
         ylim = c(ylim0[1],ylim0[2])
         )
  
  par(new=TRUE)
  lines(c(peak_rma,peak_rma),ylim0,lty = 2,col = 'red',lwd = 2)
  par(new=TRUE)
  lines(c(peak_gam,peak_gam),ylim0,lty = 2,col = 'blue',lwd = 2)
  dev.off()
}



