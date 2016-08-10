###################
# Shock Mitochondrial DNA Analysis
# Madeline Williams
# July 20, 2016

# HOUSEKEEPING (loading data, formatting, etc.)
# Load libraries 
library(reshape2)
library(ggplot2)
library(stringr)
library(GGally)
library(survival)


# Load necessary functions
rz.transform<-function(y) {
  rankY=rank(y, ties.method="average", na.last="keep")
  rzT=qnorm(rankY/(length(na.exclude(rankY))+1))  
  rzT
}

panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y, use="complete.obs")
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor*1.0, col=c("gray60", "black")[(abs(r)>0.5)+1])
}

panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2],0,1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}



####################
# 4WC mice are next in the analysis

# Load the cleaned 4wc data for p-values on the correlation coefficients
fwc <- read.csv("http://melania.jax.org/files/projects1/186_JAC_4way_Cross/phenotypes/FWC.phenotypes.cleaned.csv")

# Change the sexs from false to F
fwc$sex[fwc$sex == "FALSE"] <- "F"

# Make Object for mice that lived past 1 year to do analysis on
fwc.surv <- subset(fwc, fwc$Lifespan >= 365)

# Compute the p-values for the correlation coefficients 
cor.test(fwc.surv$Lifespan,fwc.surv$IGF1.24wk)
cor.test(fwc.surv$Lifespan,fwc.surv$IGF1.52wk)
cor.test(fwc.surv$Lifespan,fwc.surv$IGF1.76wk)
# All p-values for the the correlations are significant (below 0.05)
# Anova p-values will be the same as these because there are no covariates in the data

# Run an ANOVA regression on all the phenotypes in the fwc.surv frame
# as predictors of lifespan
sig.fwc.predictors <- data.frame()
for(i in c(28,31,34)){
  anovaresult.fwc <- anova(lm(fwc.surv$Lifespan~ fwc.surv$pgm 
                               + fwc.surv[,i])) #lm(lifespan~Strains+phenotype)
  
  insert.fwc <- cbind(colnames(fwc.surv[i]), anovaresult.fwc[2,5])
  sig.fwc.predictors=rbind(sig.fwc.predictors, insert.fwc)
  
}
colnames(sig.fwc.predictors) <- c("phenotype", "| p-val for phenotype predicting lifespan")
print(sig.fwc.predictors)


##################
# PLxKK cross IGF1 Analysis

# Load data
plkk <- read.csv("http://melania.jax.org/files/projects1/187_JAC_PLxKK/phenotypes/KKPL.phenotypes.final.csv")

# Make Object for mice that lived past 1 year to do analysis on
plkk.surv <- subset(plkk, plkk$Lifespan >= 365)

# Make correlation coefficients 
cor.test(plkk.surv$Lifespan, plkk.surv$IGF1.7wk)
cor.test(plkk.surv$Lifespan, plkk.surv$IGF1.16wk)
cor.test(plkk.surv$Lifespan, plkk.surv$IGF1.24wk)
cor.test(plkk.surv$Lifespan, plkk.surv$IGF1.52wk)
cor.test(plkk.surv$Lifespan, plkk.surv$IGF1.76wk)

# Run an ANOVA regression on all the phenotypes in the plkk.surv frame
# as predictors of lifespan
sig.plkk.predictors <- data.frame()
for(i in c(24,28,32,36,40)){
  anovaresult.plkk <- anova(lm(plkk.surv$Lifespan~ plkk.surv$Sex + plkk.surv$Strains 
                          + plkk.surv[,i])) #lm(lifespan~Sex+Strains+phenotype)
  
  insert.plkk <- cbind(colnames(plkk.surv[i]), anovaresult.plkk[3,5])
  sig.plkk.predictors=rbind(sig.plkk.predictors, insert.plkk)
  
}
colnames(sig.plkk.predictors) <- c("phenotype", "| p-val for phenotype predicting lifespan")
print(sig.plkk.predictors)



