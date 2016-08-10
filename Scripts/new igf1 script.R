#################
# IGF-1 Statistical Analysis (mtDNA analysis included)
# Madeline Williams
# August 1st, 2016
# Last modified August 5, 2016

# Load libraries
library(stringr)
library(reshape2)
library(survival)
library(ggplot2)
library(GGally)

# Load functions
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
#
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2],0,1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}

# Load data files
lifespan_data <- read.csv("/Users/s-willim/IGF1_project-master/Aging/Data/Data_for_Maddy.Elliot/JAC_long_lifespan.csv")
mtDNA.data <- read.csv("/Users/s-willim/Desktop/shock_DO_inventory2.csv")
last.igf1.data <- read.csv("/Users/s-willim/Downloads/JAC_longitudinal_IGF1.csv")

# Reshape the data from long form to extra short form
long.new.igf1 <- melt(last.igf1.data,
                       id.vars= c("Mouse.ID", "Sex", "Generation", "Age.at.exp.date"), 
                       measure.vars=c("BW", "Glucose", "IGF1..ng.mL."), 
                       variable.name="phenotype", 
                       value.name="count")

extra.short.new <- dcast(long.new.igf1,
                          Mouse.ID + Sex + Generation ~ Age.at.exp.date + phenotype,
                          value.var = "count", fun.aggregate = NULL)

# Remove the ear tag information from the mouse id column
mouse.id.ear.marks <- str_split_fixed(lifespan_data$Mouse.ID, " ", n = 2)
lifespan_data <- cbind(mouse.id.ear.marks, lifespan_data)
colnames(lifespan_data)[1] <- "mouse.id"

# Match the IGF.1 data frame to the length of the lifespan data to get rid of mice that 
# don't have matching data
a <- extra.short.new$Mouse.ID
b <- lifespan_data$mouse.id
c <- which(a %in% b)
d <- which(b %in% a)
extra.short.new <- extra.short.new[c,]
lifespan_data <- lifespan_data[d,]
rm(a,b,c,d)

# Get rid of the duplicated row
extra.short.new <- extra.short.new[!duplicated(extra.short.new$Mouse.ID),]

# Now bind lifespan and event data to the phenotype data
extra.short.new <- cbind(extra.short.new, lifespan_data$lifespan, lifespan_data$event)

# Rename the columns of extra.short.new
colnames(extra.short.new) <- c("Mouse.ID", "Sex", "Generation", "6.BW", 
                                "6.Glu", "6.IGF1", "12.BW", "12.Glu", 
                                "12.IGF1", "18.BW", "18.Glu", 
                                "18.IGF1", "Lifespan", "event")

# Change the typo in the body weight catergory
extra.short.new$`6.BW`[extra.short.new$`6.BW` == 3273] <- 32.73 

# Make the change over time phenotypes and add to delta.phenos matrix
delta.phenos <- cbind((as.numeric(extra.short.new[,7])-as.numeric(extra.short.new[,4])), (as.numeric(extra.short.new[,10])-as.numeric(extra.short.new[,7])), 
                      (as.numeric(extra.short.new[,10])-as.numeric(extra.short.new[,4])), (as.numeric(extra.short.new[,8])-as.numeric(extra.short.new[,5])),
                      (as.numeric(extra.short.new[,11])-as.numeric(extra.short.new[,8])), (as.numeric(extra.short.new[,11])-as.numeric(extra.short.new[,5])),
                      (as.numeric(extra.short.new[,9])-as.numeric(extra.short.new[,6])), (as.numeric(extra.short.new[,12])-as.numeric(extra.short.new[,9])),
                      (as.numeric(extra.short.new[,12])-as.numeric(extra.short.new[,6])))

# Rename columns
colnames(delta.phenos) <- c("BW 12-6", "BW 18-12", "BW 18-6","Glu 12-6", "Glu 18-12", "Glu 18-6",
                            "IGF1 12-6", "IGF1 18-12", "IGF1 18-6")

# Bind the regular phenotypes to the change over time phenotypes
extra.short.new <- cbind(extra.short.new, delta.phenos)

# Fix the mouse.ids to have a period instead of a dash 
extra.short.new$Mouse.ID <- gsub("-", ".", extra.short.new$Mouse.ID)

# Rename the phenotype data from extra.short.new to shock phenos
shock.phenos <- extra.short.new

# RZ transform all of the phenotypes
rz.bw6 <- rz.transform(shock.phenos[,4])
rz.glu6 <- rz.transform(shock.phenos[,5])
rz.igf6 <- rz.transform(shock.phenos[,6])
rz.bw12 <- rz.transform(shock.phenos[,7])
rz.glu12 <- rz.transform(shock.phenos[,8])
rz.igf12 <- rz.transform(shock.phenos[,9])
rz.bw18 <- rz.transform(shock.phenos[,10])
rz.glu18 <- rz.transform(shock.phenos[,11])
rz.igf18 <- rz.transform(shock.phenos[,12])
rz.lifespan <- rz.transform(shock.phenos[,13])
rz.bw12.6 <- rz.transform(shock.phenos[,14])
rz.bw18.12 <- rz.transform(shock.phenos[,15])
rz.bw18.6 <- rz.transform(shock.phenos[,16])
rz.glu12.6 <- rz.transform(shock.phenos[,17])
rz.glu18.12 <- rz.transform(shock.phenos[,18])
rz.glu18.6 <- rz.transform(shock.phenos[,19])
rz.igf12.6 <- rz.transform(shock.phenos[,20])
rz.igf18.12 <- rz.transform(shock.phenos[,21])
rz.igf18.6 <- rz.transform(shock.phenos[,22])

# Bind rz phenos to the shock.phenos matrix
shock.phenos <- cbind(shock.phenos, rz.bw6, rz.bw12, rz.bw18, rz.glu6, rz.glu12, rz.glu18,
                      rz.igf6, rz.igf12, rz.igf18, rz.bw12.6, rz.bw18.12, rz.bw18.6,
                      rz.glu12.6, rz.glu18.12, rz.glu18.6, rz.igf12.6, rz.igf18.12,
                      rz.igf18.6, rz.lifespan)

# Now that that's out of the way, move on to statistical analysis.
# Before anything else, need to run correlation/anova/coxph analyses.

# Make a threshold so that mice that lived less than 365 days can be separated from those
# that did live longer than a year.
threshold <- 365
# Make separate data sets for the mice that lived past 365 days and those that didn't 
surv_more_than_year <- subset(shock.phenos, shock.phenos$Lifespan >= threshold)
surv_less_than_year <- subset(shock.phenos, shock.phenos$Lifespan < threshold)

# Make a pairwise scan for each of the following sets of phenos:

# 6 month phenos vs. lifespan
pairs(surv_more_than_year[,c(24,27,30,42)],
      use = "pairwise.complete.obs", 
      upper.panel = panel.cor, 
      diag.panel = panel.hist)

# 12 month vs. lifespan
pairs(surv_more_than_year[,c(25,28,31,42)],
      use = "pairwise.complete.obs", 
      upper.panel = panel.cor, 
      diag.panel = panel.hist)

# 18 month vs. lifespan
pairs(surv_more_than_year[,c(26,29,32,42)],
      use = "pairwise.complete.obs", 
      upper.panel = panel.cor, 
      diag.panel = panel.hist)

# BW phenos vs. lifespan
pairs(surv_more_than_year[,c(24:26,33:35,42)],
      use = "pairwise.complete.obs", 
      upper.panel = panel.cor, 
      diag.panel = panel.hist)

# Glu vs. lifespan
pairs(surv_more_than_year[,c(27:29,36:38,42)],
      use = "pairwise.complete.obs", 
      upper.panel = panel.cor, 
      diag.panel = panel.hist)

# IGF1 vs. lifespan
pairs(surv_more_than_year[,c(30:32,39:41,42)],
      use = "pairwise.complete.obs", 
      upper.panel = panel.cor, 
      diag.panel = panel.hist)

#### Statistical Tests ####

# Check if there is a difference between the 6 month phenotypes between the two groups of mice
# (those that lived past one year and those that didn't)
t.test(surv_less_than_year$rz.igf6, surv_more_than_year$rz.igf6, paired = FALSE) #  p = 0.7934
t.test(surv_less_than_year$rz.glu6, surv_more_than_year$rz.glu6, paired = FALSE) #  p = 0.7757
t.test(surv_less_than_year$rz.bw6, surv_more_than_year$rz.bw6, paired = FALSE) #  p = 0.1055

# Now assess the significance of the individual correlation coefficients that each IGF-1
# phenotype has with lifespan.
cor.test(surv_more_than_year$rz.igf6, surv_more_than_year$rz.lifespan) # p = 0.009047**
cor.test(surv_more_than_year$rz.igf12, surv_more_than_year$rz.lifespan) # p = 0.172
cor.test(surv_more_than_year$rz.igf18, surv_more_than_year$rz.lifespan) # p = 0.05688
cor.test(surv_more_than_year$rz.igf12.6, surv_more_than_year$rz.lifespan) # p = 0.02617*
cor.test(surv_more_than_year$rz.igf18.12, surv_more_than_year$rz.lifespan) # p = 0.3742
cor.test(surv_more_than_year$rz.igf18.6, surv_more_than_year$rz.lifespan) # p = 0.007651**

# Assess if any of the phenotypes (IGF-1, BW, or Glucose) are significant predictors of lifespan
sig.lifespan.predictors <- data.frame()
for(i in c(24:41)){
  anovaresult <- anova(lm(surv_more_than_year$rz.lifespan~ surv_more_than_year$Sex + surv_more_than_year$Generation 
                          + surv_more_than_year[,i]))
                          #lm(lifespan~Sex+Gen+phenotype)
  if(anovaresult[3,5] < 0.05) { # Tells function to compile all p-values that are significant
  insert <- cbind(colnames(surv_more_than_year[i]), anovaresult[3,5])
  sig.lifespan.predictors=rbind(sig.lifespan.predictors, insert)
  }
}
colnames(sig.lifespan.predictors) <- c("phenotype", "| p-val for phenotype predicting lifespan")
print(sig.lifespan.predictors)
# IGF-1 at 6 months is a significant predictor(p= 0.0169*), but want to look at cox
# ph regressions before moving on to other stats tests


# Run Cox PH regressions for all phenos as predictors of lifespan
for(i in c(24:41)){
  coxph.tests = coxph(Surv(Lifespan) ~ Sex + Generation + surv_more_than_year[,i],
                      data = surv_more_than_year)
  # (lifespan ~ sex + gen. + pheno)
  assign(paste0("coxph.",colnames(surv_more_than_year)[i]), coxph.tests)
}


# See if any of the igf-1 phenos are significant predictors of lifespan
summary(coxph.rz.igf6) # p = 0.00494 **
summary(coxph.rz.igf12) # p = 0.0262 *
summary(coxph.rz.igf18) # p = 0.4068 
summary(coxph.rz.igf12.6) # p = 0.1791
summary(coxph.rz.igf18.12) # p = 0.701
summary(coxph.rz.igf18.6) # p =  0.0119 *


# Make change over time graphs to better understand the phenotypes
age_at_test <- factor(last.igf1.data$Age.at.exp.date)
levels(age_at_test)[levels(age_at_test)=="6"] <- "6.Months"
levels(age_at_test)[levels(age_at_test)=="12"] <- "12.Months"
levels(age_at_test)[levels(age_at_test)=="18"] <- "18.Months"

# Start with body weight and compile data from orignal igf-1 file
bw_data <- cbind(last.igf1.data[,c("Sex", "BW")], age_at_test)

# Change outlier
bw_data[1020,2] = 32.73

# Make an object for mean bw by age and sex
mean_bw_by_age <- dcast(bw_data, age_at_test ~ Sex, fun.aggregate = mean, na.rm = TRUE, value.var= "BW")

# Convert to correct format for line graph 
bw_linegraph_data <- melt(mean_bw_by_age, id.vars= c("age_at_test"), measure.vars= c("F", "M"), value.name = "BW")

# Change the sex columns name to Sex
colnames(bw_linegraph_data)[2] <- "Sex"

# Compute the standard error
bw_se <- sd(na.omit(bw_data$BW))/sqrt(565)

# Plot line graph of bw change over time (by sex) with standard error lines
ggplot(data = bw_linegraph_data, aes(x= age_at_test, y= BW, colour= Sex, group= Sex)) +
  geom_line() + 
  geom_errorbar(data = bw_linegraph_data, aes(ymin=BW-bw_se, ymax=BW+bw_se, width=.05)) + 
  theme_bw() +
  ggtitle("Change of Body Weight over Time") +
  xlab("Age") +
  ylab("Body Weight")

# Repeat these steps for glucose and igf1 phenotypes.
# See code below for those steps

# Make Glu Change Over Time Graph
glu_data <- cbind(last.igf1.data[,c("Sex", "Glucose")], age_at_test)
mean_glu_by_age <- dcast(glu_data, age_at_test ~ Sex,
                         fun.aggregate = mean, na.rm = TRUE, value.var= "Glucose")
glu_se <- sd(glu_data$Glucose, na.rm = TRUE)/sqrt(565)
glu_linegraph_data <- melt(mean_glu_by_age, id.vars= c("age_at_test"), measure.vars= c("F", "M"), value.name = "Glucose")
ggplot(data = glu_linegraph_data, aes(x= age_at_test, y= Glucose, colour= variable, group= variable))+
  geom_line() +
  geom_errorbar(data = glu_linegraph_data, aes(ymin=Glucose-glu_se, ymax=Glucose+glu_se, width=.05)) + 
  theme_bw()


# Make IGF1 Change Over Time Graph
igf1_challenge_data <- cbind(last.igf1.data[,c("Sex", "IGF1..ng.mL.")], age_at_test)


mean_igf1_by_age <- dcast(igf1_challenge_data, age_at_test ~ Sex, 
                          fun.aggregate = mean, na.rm = TRUE, value.var= "IGF1..ng.mL.")
igf1_se <- sd(na.omit(igf1_challenge_data$IGF1..ng.mL.))/sqrt(565)
igf1_linegraph_data <- melt(mean_igf1_by_age, id.vars= c("age_at_test"), measure.vars= c("F", "M"), value.name = "IGF1..ng.mL.")
# Change the name of sex column back to sex
colnames(igf1_linegraph_data)[2] <- "Sex"
ggplot(data = igf1_linegraph_data, aes(x= age_at_test, y= IGF1..ng.mL., colour= Sex, group= Sex))+
  geom_line() +
  geom_errorbar(data = igf1_linegraph_data, aes(ymin=IGF1..ng.mL.-igf1_se,ymax=IGF1..ng.mL.+igf1_se, width=.05)) + 
  theme_bw() +
  ggtitle("Change of IGF1 over Time") +
  xlab("Age") +
  ylab("IGF1 Levels")

#### Mitochondrial Analysis ####

# Match the rows in the haplotype data to the phenotype data frame
a <- surv_more_than_year$Mouse.ID
b <- mtDNA.data$JCMS.Mouse.ID
c <- which(a %in% b)
d <- which(b %in% a)
surv_more_than_year <- surv_more_than_year[c,]
mtDNA.data <- mtDNA.data[d,]
rm(a,b,c,d)
# Leaves 540 mice

# Now bind the mtDNA genotypes to the phenotype data frame and rename it to mtDNA
surv_more_than_year <- cbind(mtDNA.data$chrM, surv_more_than_year)
colnames(surv_more_than_year)[1] <- "mtDNA"

# Run an anova regression to see which phenotypes mtDNA is a good predictor of
sig.mtDNA.shock.phenos <- data.frame()
for(i in 25:43){
  anovaresult2 <- anova(lm(surv_more_than_year[,i] ~ surv_more_than_year$Sex + 
                           surv_more_than_year$Generation + as.factor(surv_more_than_year$mtDNA)))
  #lm(phenotype ~ Sex + Gen. + mtDNA geno)
 # if(anovaresult2[3,5] < 0.05){ # compiles all phenotypes that mtDNA sig. predicts
  insert2 <- cbind(colnames(surv_more_than_year)[i],anovaresult2[1,5], anovaresult2[2,5], anovaresult2[3,5])
  sig.mtDNA.shock.phenos = rbind(sig.mtDNA.shock.phenos, insert2)
 #  }
}
colnames(sig.mtDNA.shock.phenos) <- c("phenotype", "pval for Sex", "pval for Gen", "pval for mtDNA")
print(sig.mtDNA.shock.phenos)
# IGF1 at 6 Months is significantly predicted by mtDNA haplotype
# Want to see what the distribution is on that with summary function

# Make subsets with only mice of the same mtDNA haplotype for individual pairwise scans
mtDNA.ABCD.mice <- subset(surv_more_than_year, surv_more_than_year$mtDNA == "ABCD")

mtDNA.E.mice <- subset(surv_more_than_year, surv_more_than_year$mtDNA == "E")

mtDNA.F.mice <- subset(surv_more_than_year, surv_more_than_year$mtDNA == "F")

mtDNA.G.mice <- subset(surv_more_than_year, surv_more_than_year$mtDNA == "G")

mtDNA.H.mice <- subset(surv_more_than_year, surv_more_than_year$mtDNA == "H")


# Look at scatterplot matrices to understand IGF-1's relationship to lifespan
# in the individual mtDNA haplotype groups
pairs(mtDNA.ABCD.mice[,c(31:33,40:43)],
      use = "pairwise.complete.obs", 
      upper.panel = panel.cor, 
      diag.panel = panel.hist)
pairs(mtDNA.E.mice[,c(31:33,40:43)],
      use = "pairwise.complete.obs", 
      upper.panel = panel.cor, 
      diag.panel = panel.hist)
pairs(mtDNA.F.mice[,c(31:33,40:43)],
      use = "pairwise.complete.obs", 
      upper.panel = panel.cor, 
      diag.panel = panel.hist)
pairs(mtDNA.G.mice[,c(31:33,40:43)],
      use = "pairwise.complete.obs", 
      upper.panel = panel.cor, 
      diag.panel = panel.hist)
pairs(mtDNA.H.mice[,c(31:33,40:43)],
      use = "pairwise.complete.obs", 
      upper.panel = panel.cor, 
      diag.panel = panel.hist)

# Now I want to see if there is an interaction between igf1 at 6 months and mtDNA
mtDNA.interaction <- data.frame()
for(i in c(26:43)){
  anovaresult.inter <- anova(lm(surv_more_than_year$rz.lifespan ~ surv_more_than_year$Sex + surv_more_than_year$Generation +
                                surv_more_than_year$mtDNA*surv_more_than_year[,i]))
  #lm(lifespan ~ Sex + Gen. + mtDNA + phenotype)
  insert.inter <- cbind((colnames(surv_more_than_year)[i]), anovaresult.inter[1,5], 
                       anovaresult.inter[2,5], anovaresult.inter[3,5], anovaresult.inter[4,5],
                       anovaresult.inter[5,5])
  mtDNA.interaction = rbind(mtDNA.interaction, insert.inter)
}
print(mtDNA.interaction)
# There is no interaction between IGF-1 at 6 Months and mtDNA haplotype (at least not a significant
# one becasue the p-value is 0.6209)

# Plot the survival curve with the separate mtDNA groups to visualize the relationship between
# lifespan and mtDNA haplotype
surv.by.mtDNA <- survfit(Surv(Lifespan, event == 1) ~ mtDNA, data = surv_more_than_year, conf.type = "log-log")
ggsurv(surv.by.mtDNA, plot.cens = FALSE, surv.col = c("red", "blue", "green", "orange", "black"),
       back.white = TRUE, xlab = "days", ylab = "proportion alive", 
       main = "Survival of Shock Mice by mtDNA")
# There doesn't seem to be any shocking trends here

# Use survdiff function to determine if there is a sig. difference between the lifespans of 
# the different mtDNA groups
age.mtDNA.survdiff = survdiff(Surv(Lifespan, event == 1) ~ mtDNA, data = surv_more_than_year)
age.mtDNA.survdiff 
# Appears that there is not a significant difference because the p-value is 0.917

# Now run coxph regressions again with mtDNA as a covariate
for(i in c(26:43)){
  coxph.tests.mtDNA = coxph(Surv(Lifespan) ~ Sex + Generation + mtDNA + surv_more_than_year[,i],
                      data = surv_more_than_year)
  assign(paste0("coxph.mtDNA.",colnames(surv_more_than_year)[i]), coxph.tests.mtDNA)
}

summary(coxph.mtDNA.rz.bw6) # p = 0.000338 ***
summary(coxph.mtDNA.rz.bw12) # p = 0.0187 *
summary(coxph.mtDNA.rz.bw18) # p = 0.7604
summary(coxph.mtDNA.rz.glu6) # p = 0.8579
summary(coxph.mtDNA.rz.glu12) # p = 0.2035  
summary(coxph.mtDNA.rz.glu18) # p = 0.00469 **
summary(coxph.mtDNA.rz.igf6) # p = 0.00141 **
summary(coxph.mtDNA.rz.igf12) # p = 0.0682 .
summary(coxph.mtDNA.rz.igf18) # p = 0.6056  
summary(coxph.mtDNA.rz.bw12.6) # p = < 2e-16 *** # Generations also have low p-values?? tied to gen
summary(coxph.mtDNA.rz.bw18.12) # p = 0.2781  
summary(coxph.mtDNA.rz.bw18.6) # p = 0.000457 ***
summary(coxph.mtDNA.rz.glu12.6) # p = 0.0130 *
summary(coxph.mtDNA.rz.glu18.12) # p = 0.3569
summary(coxph.mtDNA.rz.glu18.6) # p = 3.88e-05 ***
summary(coxph.mtDNA.rz.igf12.6) # p = 0.0572 .
summary(coxph.mtDNA.rz.igf18.12) # p = 0.474 
summary(coxph.mtDNA.rz.igf18.6) # p = 0.1307

short.scatterplot.data <- subset(surv_more_than_year, !(is.na(surv_more_than_year$mtDNA)))

# Now that I'm suspicious that lifespan and mtDNA actually are connected (sort of) I want to 
# see a scatterplot with the mtDNA groups separated by color with lifespan on the y axis
qplot(as.numeric(`6.IGF1`), as.numeric(Lifespan), data = short.scatterplot.data, colour = mtDNA) + 
      geom_smooth(aes(group = mtDNA), method = "lm", se = FALSE) +
      xlab("IGF-1 at 6 Months") + ylab("Lifespan") + theme_bw()
# F and G mtDNA phenos look fairly different from the rest


#### Sex by Sex Analysis of IGF-1 at 6 Months ####

males <- subset(surv_more_than_year, surv_more_than_year$Sex == "M")
females <- subset(surv_more_than_year, surv_more_than_year$Sex == "F")

# scatterplot of lifespan and igf-1 at 6 months in mice that survived past 1 yr
ggplot(surv_more_than_year, aes(x = as.numeric(`6.IGF1`), y = Lifespan, group = Sex, colour = Sex)) + 
  geom_point() +
  geom_smooth(aes(group = Sex, colour = Sex), method = "lm", se = FALSE) +
  theme_bw() +
  ggtitle("IGF-1 at 6 Months vs. Lifespan") +
  xlab("6 Month IGF-1 Levels") +
  ylab("Lifespan")

# Lifespan vs. igf-1 at 6 months in males
ggplot(males, aes(x = as.numeric(`6.IGF1`), y = Lifespan)) + geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  ggtitle("Male IGF-1 at 6 Months vs. Lifespan") +
  xlab("6 Month IGF-1 Levels") +
  ylab("Lifespan")

# lifespan vs. igf-1 at 6 months in females
ggplot(females, aes(x = as.numeric(`6.IGF1`), y = Lifespan)) + geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  ggtitle("Female IGF-1 at 6 Months vs. Lifespan") +
  xlab("6 Month IGF-1 Levels") +
  ylab("Lifespan")

# Need to check the correlation coefficients and p-values for the relationship between 
# lifespan and igf-1 at 6 months in both sexs
cor.test(as.numeric(males$`6.IGF1`), males$Lifespan) #  p = 0.8228
cor.test(as.numeric(females$`6.IGF1`), females$Lifespan) #  p = 0.03194*
# IGF-1 at 6 months in females is a sig. predictor of lifespan

# Now check if mtDNA is a predictor of any phenotypes in either males or females
# Males
male.mtDNA.predictors <- data.frame()
for(i in 25:43){
  anovaresult.males <- anova(lm(males[,i] ~ males$Generation +
                                  as.factor(males$mtDNA)))
  #lm(phenotype ~ Gen. + mtDNA geno)
  if(anovaresult.males[2,5] < 0.05){
    insert.males <- cbind(colnames(males)[i],anovaresult.males[1,5], anovaresult.males[2,5])
    male.mtDNA.predictors = rbind(male.mtDNA.predictors, insert.males)
  }
}
colnames(male.mtDNA.predictors) <- c("phenotype", "pval for Gen", "pval for mtDNA")
print(male.mtDNA.predictors)
# mtDNA in males predicts IGF-1 at 6 months

# Females
female.mtDNA.predictors <- data.frame()
for(i in 25:43){
  anovaresult.females <- anova(lm(females[,i] ~ females$Generation +
                                  as.factor(females$mtDNA)))
  #lm(phenotype ~ Gen. + mtDNA geno)

    insert.females <- cbind(colnames(females)[i],anovaresult.females[1,5], anovaresult.females[2,5])
    female.mtDNA.predictors = rbind(female.mtDNA.predictors, insert.females)

}
colnames(female.mtDNA.predictors) <- c("phenotype", "pval for Gen", "pval for mtDNA")
print(female.mtDNA.predictors)
# mtDNA predicts nothing in the females in this data


#### Survival Curves of both Sexes #### 


# Males first
median(as.numeric(males$`6.IGF1`)) #  Median is 576.5223
# now make a factor for mice above and below the median 6 month IGF-1 reading

# make two groups of male mice: above and below median IGF-1 levels at 6 months
six.mon.igf1 <- factor(males$`6.IGF1`)
levels(six.mon.igf1)[levels(six.mon.igf1) >= 576.5223] <- "High IGF-1"
levels(six.mon.igf1)[levels(six.mon.igf1) < 576.5223] <- "Low IGF-1"

# Make male survival curve
surv.of.males <- survfit(Surv(Lifespan, event == 1) ~ six.mon.igf1, data = males, conf.type = "log-log")
ggsurv(surv.of.males, plot.cens = FALSE, surv.col = c("black", "purple"),
       back.white = TRUE, xlab = "days", ylab = "proportion alive", 
       main = "Males")

male.igf1.survdiff = survdiff(Surv(Lifespan, event == 1) ~ six.mon.igf1, data = males)
male.igf1.survdiff
# p= 0.98 which means there is no difference in male survival when the group is split up 
# by high and low igf-1 levels at 6 months

qplot(mtDNA,as.numeric(`6.IGF1`), data = males) + 
  geom_smooth(method = "lm", se = FALSE) +
  xlab("mtDNA") + 
  ylab("IGF-1 at 6 Months") + 
  theme_bw()


# Females next
median(as.numeric(females$`6.IGF1`))

# make two groups of female mice: above and below median IGF-1 levels at 6 months
six.mon.igf1.fe <- factor(females$`6.IGF1`)
levels(six.mon.igf1.fe)[levels(six.mon.igf1.fe) >= 494.1141] <- "High IGF-1"
levels(six.mon.igf1.fe)[levels(six.mon.igf1.fe) < 494.1141] <- "Low IGF-1"

# Make survival curve
surv.of.females <- survfit(Surv(Lifespan, event == 1) ~ six.mon.igf1.fe, data = females, conf.type = "log-log")
ggsurv(surv.of.females, plot.cens = FALSE, surv.col = c("purple", "black"),
       back.white = TRUE, xlab = "days", ylab = "proportion alive", 
       main = "Females")

female.igf1.survdiff = survdiff(Surv(Lifespan, event == 1) ~ six.mon.igf1.fe, data = females)
female.igf1.survdiff
# p= 0.00529 which means there is a sig. difference in female survival when the
# group is split up by high and low igf-1 levels at 6 months

qplot(mtDNA,as.numeric(`6.IGF1`), data = females) + 
  geom_smooth(method = "lm", se = FALSE) +
  xlab("mtDNA") + 
  ylab("IGF-1 at 6 Months") + 
  theme_bw()


# Now just check if any phenotypes are significant predictors of lifespan in 
# the male or female populations

male.lifespan.predictors <- data.frame()
for(i in c(25:42)){
  anovaresult.lm <- anova(lm(males$rz.lifespan~ males$Generation 
                          + males[,i])) #lm(lifespan~Gen+phenotype)
  if(anovaresult.lm[2,5] < 0.05){
    insert.lm <- cbind(colnames(males[i]), anovaresult.lm[2,5])
    male.lifespan.predictors=rbind(male.lifespan.predictors, insert.lm)
  }
}
colnames(male.lifespan.predictors) <- c("phenotype", "| p-val for phenotype predicting lifespan")
print(male.lifespan.predictors)
# A few phenos like bw6/12 and igf1 at 18 months are predictors of male lifespan

# female predictors next
female.lifespan.predictors <- data.frame()
for(i in c(25:42)){
  anovaresult.lf <- anova(lm(females$rz.lifespan~ females$Generation 
                             + females[,i])) #lm(lifespan~Gen+phenotype)
  if(anovaresult.lf[2,5] < 0.05){
    insert.lf <- cbind(colnames(females[i]), anovaresult.lf[2,5])
    female.lifespan.predictors=rbind(female.lifespan.predictors, insert.lf)
  }
}
colnames(female.lifespan.predictors) <- c("phenotype", "| p-val for phenotype predicting lifespan")
print(female.lifespan.predictors)
# IGF-1 at 6 months and IGF-1 at 12 months are significant predictors of lifespan in females

