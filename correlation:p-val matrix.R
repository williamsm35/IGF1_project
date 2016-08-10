##############
# Correlation/P-value Matrix Script for Shock Data
# Madeline Williams
# August 4, 2016
# Last modified August 5, 2016

# Load libraries
library(RColorBrewer)
options(stringsAsFactors = FALSE)

# Load Data
lifespan_data <- read.csv("/Users/s-willim/IGF1_project-master/Aging/Data/Data_for_Maddy.Elliot/JAC_long_lifespan.csv")
last.igf1.data <- read.csv("/Users/s-willim/Downloads/JAC_longitudinal_IGF1.csv")

# Split the ear mark data from the Mouse ID's in the lifespan data
mouse.id.ear.marks <- str_split_fixed(lifespan_data$Mouse.ID, " ", n = 2)
lifespan_data <- cbind(mouse.id.ear.marks, lifespan_data)

# Name the new column Mouse.ID
colnames(lifespan_data)[1] <- "Mouse.ID"

# Remove the ear marks and extra id columns from the lifespan data
lifespan_data <- lifespan_data[,-c(2:4)]

# Match the igf1 data to the mice in the lifespan data frame
a <- last.igf1.data$Mouse.ID
b <- lifespan_data$Mouse.ID
c <- which(a %in% b)
d <- which(b %in% a)
last.igf1.data <- last.igf1.data[c,]
lifespan_data <- lifespan_data[d,]
rm(a,b,c,d)

# Merge the lifespan data with the igf1 data now that they include the same mouse ids
m = merge(last.igf1.data, lifespan_data, by = "Mouse.ID", all = FALSE)

# Remove the irrelevant columns
m <- m[,-c(5:8,12,17:21)]

# Make the variables (IGF1, BW, Glucose) numeric variables
m$lifespan = as.numeric(m$lifespan)
m$BW = as.numeric(m$BW)
m$Glucose = as.numeric(m$Glucose)
m$IGF1..ng.mL. = as.numeric(m$IGF1..ng.mL.)

# Rename the IGF1 column to IGF1
colnames(m)[11] <- "IGF1"

# Keep only mice that lived >= 365 days.
m = m[m$lifespan >= 365,]

# rank z transform all of the phenotype data
for(i in 9:11) {
  
  sc = scale(m[,i])
  sc[abs(sc) > 5] = NA
  rz = rank(sc, na.last = FALSE)
  rz = qnorm(rz / (max(rz, na.rm = T) + 1))
  m[,i] = rz
  
}

# Make an object that includes the names of the different phenotypes
pheno = colnames(m)[9:11]
# Make an object that includes the specific times that the phenotype data was recorded
time = unique(m$Age.at.exp.date)
# Make a new matrix for the correlations
cor = matrix(0, nrow = 3, ncol = length(pheno), dimnames = 
               list(c(6,12,18), pheno))
# Make a new matrix for the p-values
pv = matrix(0, nrow = 3, ncol = length(pheno), dimnames = 
              list(c(6,12,18), pheno))

# Compile the p-values and correlations of each trait at each time as a predictor of lifespan
layout(matrix(1:length(cor), nrow(cor), ncol(cor)))
par(plt = c(0,1,0,1))
for(i in 1:length(pheno)) {
  
  spl = split(m[,c("Mouse.ID", "Sex.x", "Generation.x", "lifespan", pheno[i])], 
              m$Age.at.exp.date)
  for(j in 1:length(spl)) {
    
    mod1 = lm(lifespan ~ Sex.x + Generation.x, data = spl[[j]],
              na.action = na.exclude)
    spl[[j]][,4] = residuals(mod1)
    mod2 = lm(spl[[j]][,5] ~ Sex.x + Generation.x, data = spl[[j]],
              na.action = na.exclude)
    spl[[j]][,5] = residuals(mod2)
    
    t = cor.test(spl[[j]][,4], spl[[j]][,5])
    cor[j,i] = t$estimate
    pv[j,i]  = t$p.value
    plot(spl[[j]][,5], spl[[j]][,4], pch = 16)
    
  }
  
}

# Rename the columns in the cor matrix to reflect the phenotypes
colnames(cor) = c("Body Weight", "Glucose", "IGF1")
colnames(pv) = colnames(cor)

# Reorder the cor and pv data frames so that they match
ord = order(cor[3,])
cor = cor[,ord]
pv  = pv[,ord]
ord = colnames(cor)

# Make the graphic that shows the correlations and p-values with color
# brks = is a function that defines the range of the correlation color scheme.
# -160:160 means the range is -.16 to .16. Change accordingly.
png("Corplot_shock_phenos.png", width = 1500, height = 1000,
    res = 164)
layout(matrix(1:2, 1, 2), width = c(0.85, 0.15))
brks = -170:170/1000
col = rev(colorRampPalette(brewer.pal(n = 11, name = "Spectral"))(length(brks) - 1))
par(plt = c(0.15, 0.95, 0.05, 0.92))
image(1:nrow(cor), 1:ncol(cor), cor, breaks = brks, col = col, axes = F, 
      ann = F)
mtext(side = 3, line = 0.25, at = 1:3, text = paste(c(6, 12, 18), "Mo"))
mtext(side = 2, at = 1:ncol(cor), text = colnames(cor), las = 1)
for(i in 1:ncol(cor)) {
  text(x = 1:3, y = rep(i, 3), labels = format(pv[,i], digits = 3))
} # for(i)

par(plt = c(0, 0.5, 0.05, 0.92))
image(matrix(brks, nrow = 1), breaks = brks, col = col, axes = F,
      ann = F)
at = seq(1,length(brks), 50)
mtext(side = 4, line = 0.25, at = at / length(brks), text = brks[at], las = 1)
mtext(side = 3, line = 0.25, text = "Cor")
dev.off()

# Get rid of all of the original environment except for the ord object and the original data
rm(at, brks, col, cor, i , j, mod1, mod2, pheno, pv, rz, sc, spl, t, time)

# Now make a new data frame for the females
females <- subset(m, m$Sex.x == "F")

# Run through the same steps again, but for the females
pheno = colnames(females)[9:11]
time = unique(females$Age.at.exp.date)
cor = matrix(0, nrow = 3, ncol = length(pheno), dimnames = 
               list(c(6,12,18), pheno))
pv = matrix(0, nrow = 3, ncol = length(pheno), dimnames = 
              list(c(6,12,18), pheno))

layout(matrix(1:length(cor), nrow(cor), ncol(cor)))
par(plt = c(0,1,0,1))
for(i in 1:length(pheno)) {
  
  spl = split(females[,c("Mouse.ID", "Generation.x", "lifespan", pheno[i])], 
              females$Age.at.exp.date)
  for(j in 1:length(spl)) {
    
    mod1 = lm(lifespan ~ Generation.x, data = spl[[j]],
              na.action = na.exclude)
    spl[[j]][,3] = residuals(mod1)
    mod2 = lm(spl[[j]][,4] ~ Generation.x, data = spl[[j]],
              na.action = na.exclude)
    spl[[j]][,4] = residuals(mod2)
    
    t = cor.test(spl[[j]][,3], spl[[j]][,4])
    cor[j,i] = t$estimate
    pv[j,i]  = t$p.value
    plot(spl[[j]][,4], spl[[j]][,3], pch = 16)
    
  }
  
}

colnames(cor) = c("Body Weight", "Glucose", "IGF1")
colnames(pv) = colnames(cor)

cor = cor[,ord]
pv  = pv[,ord]

png("Corplot_females.png", width = 1500, height = 1000,
    res = 164)
layout(matrix(1:2, 1, 2), width = c(0.85, 0.15))
brks = -170:170/1000
col = rev(colorRampPalette(brewer.pal(n = 11, name = "Spectral"))(length(brks) - 1))
par(plt = c(0.15, 0.95, 0.05, 0.92))
image(1:nrow(cor), 1:ncol(cor), cor, breaks = brks, col = col, axes = F, 
      ann = F)
mtext(side = 3, line = 0.25, at = 1:3, text = paste(c(6, 12, 18), "Mo"))
mtext(side = 2, at = 1:ncol(cor), text = colnames(cor), las = 1)
for(i in 1:ncol(cor)) {
  text(x = 1:3, y = rep(i, 3), labels = format(pv[,i], digits = 3))
} # for(i)

par(plt = c(0, 0.5, 0.05, 0.92))
image(matrix(brks, nrow = 1), breaks = brks, col = col, axes = F,
      ann = F)
at = seq(1,length(brks), 50)
mtext(side = 4, line = 0.25, at = at / length(brks), text = brks[at], las = 1)
mtext(side = 3, line = 0.25, text = "Cor")
dev.off()


# Remove everything again
rm(at, brks, col, cor, i , j, mod1, mod2, pheno, pv, rz, spl, t, time)

# now make the male matrix
males <- subset(m, m$Sex.x == "M")

# Run through the same steps again, but for the females
pheno = colnames(males)[9:11]
time = unique(males$Age.at.exp.date)
cor = matrix(0, nrow = 3, ncol = length(pheno), dimnames = 
               list(c(6,12,18), pheno))
pv = matrix(0, nrow = 3, ncol = length(pheno), dimnames = 
              list(c(6,12,18), pheno))

layout(matrix(1:length(cor), nrow(cor), ncol(cor)))
par(plt = c(0,1,0,1))
for(i in 1:length(pheno)) {
  
  spl = split(males[,c("Mouse.ID", "Generation.x", "lifespan", pheno[i])], 
              males$Age.at.exp.date)
  for(j in 1:length(spl)) {
    
    mod1 = lm(lifespan ~ Generation.x, data = spl[[j]],
              na.action = na.exclude)
    spl[[j]][,3] = residuals(mod1)
    mod2 = lm(spl[[j]][,4] ~ Generation.x, data = spl[[j]],
              na.action = na.exclude)
    spl[[j]][,4] = residuals(mod2)
    
    t = cor.test(spl[[j]][,3], spl[[j]][,4])
    cor[j,i] = t$estimate
    pv[j,i]  = t$p.value
    plot(spl[[j]][,4], spl[[j]][,3], pch = 16)
    
  }
  
}

colnames(cor) = c("Body Weight", "Glucose", "IGF1")
colnames(pv) = colnames(cor)

cor = cor[,ord]
pv  = pv[,ord]

png("Corplot_males.png", width = 1500, height = 1000,
    res = 164)
layout(matrix(1:2, 1, 2), width = c(0.85, 0.15))
brks = -170:170/1000
col = rev(colorRampPalette(brewer.pal(n = 11, name = "Spectral"))(length(brks) - 1))
par(plt = c(0.15, 0.95, 0.05, 0.92))
image(1:nrow(cor), 1:ncol(cor), cor, breaks = brks, col = col, axes = F, 
      ann = F)
mtext(side = 3, line = 0.25, at = 1:3, text = paste(c(6, 12, 18), "Mo"))
mtext(side = 2, at = 1:ncol(cor), text = colnames(cor), las = 1)
for(i in 1:ncol(cor)) {
  text(x = 1:3, y = rep(i, 3), labels = format(pv[,i], digits = 3))
} # for(i)

par(plt = c(0, 0.5, 0.05, 0.92))
image(matrix(brks, nrow = 1), breaks = brks, col = col, axes = F,
      ann = F)
at = seq(1,length(brks), 50)
mtext(side = 4, line = 0.25, at = at / length(brks), text = brks[at], las = 1)
mtext(side = 3, line = 0.25, text = "Cor")
dev.off()