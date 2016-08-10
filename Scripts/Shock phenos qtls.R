################
# IGF-1 Shock Phenotype QTL analysis
# Madeline Williams
# August 5, 2016
# Last modified August 5, 2016

# Load libraries
library(DOQTL)

# Load data files
load(file = "/Users/s-willim/IGF1_project-master/Aging/Data/probs.Rdata")
lifespan_data <- read.csv("/Users/s-willim/IGF1_project-master/Aging/Data/Data_for_Maddy.Elliot/JAC_long_lifespan.csv")
last.igf1.data <- read.csv("/Users/s-willim/Downloads/JAC_longitudinal_IGF1.csv")

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

#### Reformat Data ####

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

# Change the column names
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

# Make the rownames the Mouse.IDs for the shock.phenos matrix
rownames(shock.phenos) <- shock.phenos$Mouse.ID

# Make two sets of data; one for those that did and didn't live past one year
threshold <- 365
# Make separate data sets for the mice that lived past 365 days and those that didn't 
surv_more_than_year <- subset(shock.phenos, shock.phenos$Lifespan >= threshold)
surv_less_than_year <- subset(shock.phenos, shock.phenos$Lifespan < threshold)

#### Pairwise Scans ####

# Make a pairwise scan for each of the following sets of phenos before making QTL plots:
# 6 month phenos vs. lifespan
pairs(surv_more_than_year[,c(24,27,30,42)],
      use = "pairwise.complete.obs", 
      upper.panel = panel.cor, 
      diag.panel = panel.hist)

# 12 month phenos vs. lifespan
pairs(surv_more_than_year[,c(25,28,31,42)],
      use = "pairwise.complete.obs", 
      upper.panel = panel.cor, 
      diag.panel = panel.hist)

# 18 month phenos vs. lifespan
pairs(surv_more_than_year[,c(26,29,32,42)],
      use = "pairwise.complete.obs", 
      upper.panel = panel.cor, 
      diag.panel = panel.hist)

# BW phenos vs. lifespan
pairs(surv_more_than_year[,c(24:26,33:35,42)],
      use = "pairwise.complete.obs", 
      upper.panel = panel.cor, 
      diag.panel = panel.hist)

# Glu phenos vs. lifespan
pairs(surv_more_than_year[,c(27:29,36:38,42)],
      use = "pairwise.complete.obs", 
      upper.panel = panel.cor, 
      diag.panel = panel.hist)

# IGF1 phenos vs. lifespan
pairs(surv_more_than_year[,c(30:32,39:41,42)],
      use = "pairwise.complete.obs", 
      upper.panel = panel.cor, 
      diag.panel = panel.hist)

#### QTL plotting ####

# Make sure to eliminate the non-matching mice from the probs data
a <- rownames(surv_more_than_year)
b <- rownames(probs)
c <- which(a %in% b)
d <- which(b %in% a)
surv_more_than_year <- surv_more_than_year[c,]
adj.probs <- probs[d,,]
rm(a,b,c,d)

# Make the covariate object/matrix
shock.addcovar = model.matrix(~Sex + Generation, data = surv_more_than_year)[,-1]
colnames(shock.addcovar)[1] = "Sex"

# Load snps
load(url("ftp://ftp.jax.org/MUGA/MM_snps.Rdata"))

# Make sure the rownames of the addcovar object match the Mouse.IDs.
# Otherwise DOQTL won't make your qtl objects
rownames(shock.addcovar) <- rownames(surv_more_than_year)

# Create a kinship object
K = kinship.probs(adj.probs, snps = MM_snps, bychr = TRUE)

# Create the qtl objects for the rz transformed phenotypes  
for (i in 24:42){
  qtl = scanone(pheno = surv_more_than_year, pheno.col = colnames(surv_more_than_year)[i], 
                probs = adj.probs, K = K, addcovar = shock.addcovar, snps = MM_snps)
  assign(paste0(colnames(surv_more_than_year)[i]), qtl)
}

# Make a file and then store the qtl plots in the file. File will show up in working directory
pdf(file = "shock qtls")
for (i in 24:42){
  plot(eval(parse(text=paste0(colnames(surv_more_than_year)[i]))),main=paste0(colnames(surv_more_than_year)[i]),
       sig.thr = thr, sig.col = c("red","blue","black"))
}
dev.off()
