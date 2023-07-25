rm(list = ls())
library(metafor)
library(stringr)
library(rms)
library(mgcv)
library(MuMIn)
library(visreg)

## load data and calculate the averaged gs and variances
datadir = shortPathName("C:/Users/Guochun Yang/OneDrive - University of Iowa/WithLi/laterality/")
savedir = shortPathName("C:/Users/Guochun Yang/OneDrive - University of Iowa/WithLi/laterality/plot_weightn/")

filename = "laterality_g_all.csv"

AICs = matrix(, nrow = 3, ncol = 6)
BICs = matrix(, nrow = 3, ncol = 6)
QMps = matrix(, nrow = 3, ncol = 6)
QMs = matrix(, nrow = 3, ncol = 6)
modelweights = matrix(, nrow = 3, ncol = 6)


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
# data$variance = data$laterality_abs
# data$variance0 = 2/data$n + (1-((8*data$n - 9)^2*(data$n - 2))/((8*data$n - 12)^2*(data$n - 1))) * (data$laterality_abs)^2
data$variance = 1/data$n
data$weight = data$n #1/data$variance

data$LR = sign(data$laterality_pmn)
data$LR12 = (data$LR+3)/2

agefull = max(data$age)+min(data$age)

  
for (i in 1:3) { #1:absolute value; 2: raw R-L value ; 3:raw L-R
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

  # models
  b <- gam(laterality~s(age),data=data, weights = weight, method = 'REML')
  bs[[i]] = b
  res.lin <- 
    rma(laterality, 
        variance, 
        mods = ~ age,
        data = data)
  res.qua <- 
    rma(laterality, 
        variance, 
        mods = ~ age + I(age^2),# + (age + I(age^2)|LR), 
        data = data)
  res.cub <-
    rma(laterality,
        variance,
        mods = ~ age + I(age^2) + I(age^3), 
        data = data)
  res.log <-
    rma(laterality,
        variance,
        mods = ~ I(log(age)) + I((log(age))^2),
        data = data)
  
  data_type <- data[which(abs(data$Type) > 0),]
  res.log_int <-
    rma(laterality,
        variance,
        mods = ~ I(log(age)) + I((log(age))^2) + Type + I(log(age)):Type + I((log(age))^2):Type,
        data = data_type)
  res.sqrt <-
    rma(laterality,
        variance,
        mods = ~ I(sqrt(age)) + age, 
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

for (i in 1:3) {
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
    
  if (minAIC[i] == 1) {
    modelbest = reses_qua[[i]]
    title = "QUADRATIC"
    sav <-
      predict(modelbest, newmods = unname(poly(
        xs, degree = 2, raw = TRUE
      )))
  } else if (minAIC[i] == 2) {
    modelbest = reses_cub[[i]]
    title = "CUBIC"
    sav <-
      predict(modelbest, newmods = unname(poly(
        log(xs), degree = 2, raw = TRUE
      )))
  } else if (minAIC[i] == 3) {
    modelbest = reses_log[[i]]
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
  }

  filename = 'log'
    
  tiff(
    file = paste(savedir, name, title, ".tiff", sep = ''),
    width = 6,
    height = 5,
    units = "in",
    res = 300
  )
  
  ylim0 <- c()
  maxvar = max(data$variance)
  ylim0[1] = min(floor((data$laterality-maxvar)*100)/100)
  ylim0[2] = max(ceiling((data$laterality+maxvar)*100)/100)
  
  par(cex.lab=2)
  par(cex.axis=2)
  par(mai=c(1,1.2,0.2,0.2))
  # par(mgp=c(2,0.5,0),xaxs='i')
  # par(mgp=c(2,1,0),mgp=c(5,0.5,0))
  # par(mgp=c(3.5,0.5,0),yaxs='i')
  plt = regplot(
    reses_qua[[i]], # this is only to make sure all dots are with the original scale (we can also use reses_cub but not log)
    mod = 2,
    lcol = 'red',
    ci=FALSE,
    pred = sav,
    xvals = xs,
    las = 1,
    digits = 1,
    bty = "l",
    xlim = c(8, 75),
    psize = .20 / sqrt(modelbest$vi),
    xlab = ' ',
    ylab = ' '
    
    #main = title
  )
  title(ylab = ylabname, mgp = c(4, 1, 0))
  title(xlab = 'Age', mgp = c(2.5, 1, 0))
  # axis(1,mgp=c(3,0.5,0))
  
  par(new=TRUE)
  # add visreg
  b <- bs[[i]]
  preds <- predict.gam(b, newdata = data.frame(age = xs))
  if(i==2){
    peak_gam <- xs[which.max(preds)]
    peak_rma <- xs[which.max(sav$pred)]
  } else {
    peak_gam <- xs[which.min(preds)]
    peak_rma <- xs[which.min(sav$pred)]
  }

  # peak_gams[i] = round(peak_gam,1)
  # visreg(b,axes=FALSE, pch = NA, xlab='',ylab='',yaxs='i',ylim = c(lim0[3],lim0[4]))
  visreg(b,"age",
         overlay = TRUE,partial = FALSE, rug = FALSE,
         line=list(lty=1, col="blue"), alpha=0,
         points=list(cex=1, pch=26),
         axes=FALSE, xlab='',ylab='',yaxs='i',ylim = c(ylim0[1],ylim0[2]))
  
  par(new=TRUE)
  lines(c(peak_rma,peak_rma),ylim0,lty = 2,col = 'red',lwd = 2)
  par(new=TRUE)
  lines(c(peak_gam,peak_gam),ylim0,lty = 2,col = 'blue',lwd = 2)
  
  dev.off()
}