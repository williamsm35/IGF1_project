##############
#Year long DO mice survival analysis
# Maddy Williams
# June 22, 2016
##############


# Making data more easily usable
# Packages to load
library(ISwR)
library(ggplot2)
library(reshape2)
library(GGally)
library(survival)

# Read lifespan and phenotype data into R and format for transformation later
igf1_data <- read.csv("/Users/s-willim/IGF1_project-master/Aging/Data/Data_for_Maddy.Elliot/shock_igf1_glucose_long_JCMS2.csv")
lifespan_data <- read.csv("/Users/s-willim/IGF1_project-master/Aging/Data/Data_for_Maddy.Elliot/JAC_long_lifespan.csv")
long_igf1_data <- melt(igf1_data,
                       id.vars= c("Mouse.ID", "Sex", "Generation", "penName", "Test.Name", "Age.at.Exp.Date"), 
                       measure.vars=c("BW", "Glucose", "IGF.1", "IGF.1.control", "Age.at.Exp.Date"), 
                       variable.name="phenotype", 
                       value.name="count")
extra_short_igf1 <- dcast(long_igf1_data,
                          Mouse.ID + Sex + Generation ~ Test.Name + phenotype,
                          value.var="count")

# Remove mice who from lifespan data that don't have matching shock data
lifespan_data2 <- lifespan_data[-c(19,40,43,91,92,103,105,145,165,169,208,211,235,248,298,
                                   386,388,418,425,428,431,433,436,462,463,464,468,493,
                                   540,541,542,544,554,568), ]

# Make a data frame with phenos and lifespan
extra_short_igf1 <- cbind(extra_short_igf1, lifespan_data2$lifespan, lifespan_data2$event)

# Functions for removing or replacing outliers
extra_short_igf1[237,9] = NA # for mouse 251 igf1 at 12
extra_short_igf1[32,9] = NA # for mouse 33 igf1 at 12
extra_short_igf1[367,4] = 32.73 # for mouse 382 BW at 6
extra_short_igf1[207,6] = NA # for mouse 219 igf1 at 6

# Make separate data frame for mice that didn't survive past a year (365 days). 
# extra_short_igf1[,19] stands for the lifespan column
threshold <- 365
surv_less_than_year <- subset(extra_short_igf1, extra_short_igf1[, 19] < threshold)

# Remove unimportant phenos (igf1 control and age at experiment date)
surv_less_than_year <- surv_less_than_year[,-c(7, 8, 12, 13, 17, 18)]

# Rename columns for easier identification
colnames(surv_less_than_year) <- c("Mouse.ID","Sex","Generation","6.Mon.BW",
                                   "6.Mon.Glu","6.Mon.IGF1","12.Mon.BW","12.Mon.Glu",
                                   "12.Mon.IGF1","18.Mon.BW","18.Mon.Glu",
                                   "18.Mon.IGF1","Lifespan", "Event")

# Remove mice who did not have lifespans greater than or equal to a year (365 days). 
# extra_short_igf1[,19] stands for the lifespan column in the data.
year_long_surv <- subset(extra_short_igf1, extra_short_igf1[, 19] >= threshold)

# Remove unimportant phenos (igf1 control and age at experiment date)
year_long_surv <- year_long_surv[,-c(7, 8, 12, 13, 17, 18)]

# Rename columns for easier identification
colnames(year_long_surv) <- c("Mouse.ID","Sex","Generation","6.Mon.BW",
                                "6.Mon.Glu","6.Mon.IGF1","12.Mon.BW","12.Mon.Glu",
                                "12.Mon.IGF1","18.Mon.BW","18.Mon.Glu",
                                "18.Mon.IGF1","Lifespan", "Event")

# Comparisons between these two data sets:
# year long survivors: 531 mice; 270 females, 261 males.
# less than year survived: 34 mice; 14 females, 20 males.
# Want to move on to t-test of six month phenotypes to see if there is 
# a difference between the mice that survived past a year and those that didn't.


#######################
### T.Tests

# BW at 6 months t.test between mouse that survived past a year and those that didn't
# [,4] indicates that the fourth column is being used which in this case is 6 mon BW
t.test(na.omit(year_long_surv[,4]), 
       na.omit(surv_less_than_year[,4]), 
       var.equal = FALSE, conf.level = 0.95)
# P-value is .04753 which means that there is likely a significant difference.
# Confidence interval corroborates this story because it does not contain 0.

# Glu at 6 months t.test
# [,5] indicates that it is column five being analyzed (which is Glu at 6 mon)
t.test(na.omit(year_long_surv[,5]),
       na.omit(surv_less_than_year[,5]),
       var.equal = FALSE, conf.level = 0.95)
# P-value is .07895 which means there is likely not a significant difference.
# Confidence interval corroborates this because it includes 0.

# IGF1 at 6 months t.test
# [,6] indicates the IGF1 at 6 mon column in the data
t.test(na.omit(year_long_surv[,6]),
       na.omit(surv_less_than_year[,6]),
       var.equal = FALSE, conf.level = 0.95)
# P-value is .7027 which means that there is likely not a significant difference.
# Confidence interval corroborates because it includes 0.

###############
# Survival analysis using survival curves, log-rank tests, and coxph tests
# Removing igf1 control and age at exp. from original data frame
extra_short_igf1 <- extra_short_igf1[,-c(7, 8, 12, 13, 17, 18)]

# Renaming phenotypes
colnames(extra_short_igf1) <- c("Mouse.ID","Sex","Generation","6.Mon.BW",
                                "6.Mon.Glu","6.Mon.IGF1","12.Mon.BW","12.Mon.Glu",
                                "12.Mon.IGF1","18.Mon.BW","18.Mon.Glu",
                                "18.Mon.IGF1","lifespan","event")

# Now on to Kaplan Meier Curves
surv_by_sex <- survfit(Surv(lifespan, event == 1) ~ Sex, data = extra_short_igf1, conf.type = "log-log")
ggsurv(surv_by_sex, plot.cens = FALSE, surv.col = c("red", "blue"),
       back.white = TRUE, xlab = "days", ylab = "proportion alive", 
       main = "Survival of DO Mice by Sex")

# Significance Test for difference between lifespan of the two sexes
age.sex.survdiff = survdiff(Surv(lifespan, event == 1) ~ Sex, data=extra_short_igf1)
age.sex.survdiff 
# Survival is significantly different between the sexes (P-value = 0.0485) for all mice that 
# have shock data and lifespan data (when all mice were included, the pvalue was .0781)

# Now to make a new event for those that died before one year
extra_short_igf1$event[extra_short_igf1$lifespan<365] <- 3

# Want to repeat the Kaplan Meier Curve with the mice that died before 1 
# year censored (don't know if that's the right term, but only the mice
# with a 1 in the event column will be analyzed)
surv_by_sex_mod <- survfit(Surv(lifespan, event == 1) ~ Sex, data = extra_short_igf1, conf.type = "log-log")
ggsurv(surv_by_sex_mod, plot.cens = FALSE, surv.col = c("red", "blue"),
       back.white = TRUE, xlab = "days", ylab = "proportion alive", 
       main = "Survival of DO Mice by Sex")


# Cox proportional hazard regression
# And boxplots to make sense of what is happening in the hazard regressions
boxplot(lifespan~Generation, data = extra_short_igf1, xlab = "Gen", ylab = "lifespan") #(x=gen, y=lifespan)
boxplot(lifespan~Sex, data = extra_short_igf1, xlab = "Sex", ylab= "lifespan") #(x= sex, y= lifespan)
boxplot(extra_short_igf1[,4]~Generation, data = extra_short_igf1, xlab = "Gen", ylab = "bw_6") #(x= gen, y=bw at 6mon)
boxplot(extra_short_igf1[,5]~Generation, data = extra_short_igf1, xlab = "Gen", ylab = "glu_6") #(x= gen, y=glu at 6mon)
boxplot(extra_short_igf1[,6]~Generation, data = extra_short_igf1, xlab = "Gen", ylab = "igf1_6") #(x= gen, y=igf1 at 6mon)
boxplot(extra_short_igf1[,7]~Generation, data = extra_short_igf1, xlab = "Gen", ylab = "bw_12") #(x= gen, y=bw at 12mon)
boxplot(extra_short_igf1[,8]~Generation, data = extra_short_igf1, xlab = "Gen", ylab = "glu_12") #(x= gen, y=glu at 12mon)
boxplot(extra_short_igf1[,9]~Generation, data = extra_short_igf1, xlab = "Gen", ylab = "igf1_12") #(x= gen, y=igf1 at 12mon)
boxplot(extra_short_igf1[,10]~Generation, data = extra_short_igf1, xlab = "Gen", ylab = "bw_18") #(x= gen, y=bw at 18mon)
boxplot(extra_short_igf1[,11]~Generation, data = extra_short_igf1, xlab = "Gen", ylab = "glu_18") #(x= gen, y=glu at 18mon)
boxplot(extra_short_igf1[,12]~Generation, data = extra_short_igf1, xlab = "Gen", ylab = "igf1_18") #(x= gen, y=igf1 at 18mon)

# bw at 6 mon = [,4]
coxph_bw6 <- coxph(Surv(lifespan, event == 1) ~ Sex + Generation + extra_short_igf1[,4],
                    data = extra_short_igf1)
summary(coxph_bw6)
# Seems that bw at 6 mon and G7 are significant predictors of lifespan(p-values 
# for regression below 0.05)

# Cox proportional hazard: glu at 6 mon = [,5]
coxph_glu6 <- coxph(Surv(lifespan, event == 1) ~ Sex + Generation + extra_short_igf1[,5],
                  data = extra_short_igf1)
summary(coxph_glu6)
# Seems that G7 is a predictor of lifespan(p-values 
# for regression below 0.05)

# igf1 at 6 mon = [,6]
coxph_igf1_6 <- coxph(Surv(lifespan, event == 1) ~ Sex + Generation + extra_short_igf1[,6],
                    data = extra_short_igf1)
summary(coxph_igf1_6)
# Again G7 and igf1 at 6 mon are significant predictors of lifespan (p-values 
# for regression below 0.05)

# Bw at 12 mon = [,7]
coxph_bw12 <- coxph(Surv(lifespan, event == 1) ~ Sex + Generation + extra_short_igf1[,7],
                   data = extra_short_igf1, use = "pairwise.complete.obs")
summary(coxph_bw12)
# P-value of .0559 (not significant)

# Glu at 12 mon = [,8]
coxph_glu12 <- coxph(Surv(lifespan, event == 1) ~ Sex + Generation + extra_short_igf1[,8],
                    data = extra_short_igf1, use = "pairwise.complete.obs")
summary(coxph_glu12)
# P-value of 0.0912  (not significant)

# IGF1 at 12 mon = [,9]
coxph_igf12 <- coxph(Surv(lifespan, event == 1) ~ Sex + Generation + extra_short_igf1[,9],
                     data = extra_short_igf1)
summary(coxph_igf12)
# P-value of .379 (not significant)

# BW at 18 mon = [,10]
coxph_bw18 <- coxph(Surv(lifespan, event == 1) ~ Sex + Generation + extra_short_igf1[,10],
                     data = extra_short_igf1, use = "pairwise.complete.obs")
summary(coxph_bw18)
# P-value = .928 (not significant)

# Glu at 18 mon = [,11]
coxph_glu18 <- coxph(Surv(lifespan, event == 1) ~ Sex + Generation + extra_short_igf1[,11],
                     data = extra_short_igf1)
summary(coxph_glu18)
# p-value = 0.0592 (not significant)

# IGF1 at 18 mon = [,12]
coxph_igf18 <- coxph(Surv(lifespan, event == 1) ~ Sex + Generation + extra_short_igf1[,12],
                     data = extra_short_igf1)
summary(coxph_igf18)
# P-value = 0.5558 (not significant)


# Pairwise scans for mice that lived past one year and cleaned up data
# 6 month matrix
pairs(year_long_surv[,c(4,5,6,13)],
      use = "pairwise.complete.obs", 
      upper.panel = panel.cor, 
      diag.panel = panel.hist)
# 12 month matrix
pairs(year_long_surv[,c(7,8,9,13)],
      use = "pairwise.complete.obs", 
      upper.panel = panel.cor, 
      diag.panel = panel.hist)
# 18 month matrix
pairs(year_long_surv[,c(10,11,12,13)],
      use = "pairwise.complete.obs", 
      upper.panel = panel.cor, 
      diag.panel = panel.hist)

# BW matrix
pairs(year_long_surv[,c(4,7,10)],
      use = "pairwise.complete.obs", 
      upper.panel = panel.cor, 
      diag.panel = panel.hist)
# Glu matrix
pairs(year_long_surv[,c(5,8,11)],
      use = "pairwise.complete.obs", 
      upper.panel = panel.cor, 
      diag.panel = panel.hist)
# IGF1 matrix
pairs(year_long_surv[,c(6,9,12)],
      use = "pairwise.complete.obs", 
      upper.panel = panel.cor, 
      diag.panel = panel.hist)