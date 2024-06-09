library(dplyr)
library(ggpubr)

source("wdi-utils.R")

datafolder <-  "~/data/WDI/"  # insert data folder name here 
whocountries <- read.csv(paste0(datafolder,"who-countries.csv"), header=FALSE)
whocodes <- whocountries[,2]
mydat <- read.csv(paste0(datafolder,"WDI-5-vbles.csv"), header=TRUE)

# Create list of country-specific data sets
wlist <- list()
ccodes <- sort(unique(mydat$countrycode))
for (co in ccodes) {
  wlist[[co]] <- subset(mydat, countrycode==co)
}

######################

# DR Congo example
cor2_plus_plot("COD","investment","fem.sec.enroll", main="")

tmp <- cor2("COD","fem.sec.enroll","investment")
intersect(tmp$yr1, tmp$yr2)
range(tmp$yr1)
range(tmp$yr2)

###################

# Grid of histograms of correlations between pairs of variables
hisdat <- NULL 
system.time(for (i in 4:7) for (j in (i+1):8) {
  v1 <- names(mydat)[i]; v2 <- names(mydat)[j]
  cat(v1,'\t',v2,'\n')
  cors <- c()
  for (m in 1:216) {  # errors occur for some countries due to insufficient (overlapping) data
    trycor <- try(cor2(names(wlist)[m], v1, v2, min.overlap = 12))
    if (!inherits(trycor, "try-error")) hisdat <- rbind(hisdat, c(names(wlist)[m],v1,v2,trycor$cor))
  }
})
hisdat <- as.data.frame(hisdat)
nrow(hisdat) / 216 / choose(5,2)   # proportion of correlations that can be computed
names(hisdat) <- c("countrycode","v1", "v2", "correlation")
hisdat$correlation <- as.numeric(hisdat$correlation)
hisdat$v1 <- factor(hisdat$v1, names(mydat)[4:7)])
hisdat$v2 <- factor(hisdat$v2, names(mydat)[5:8])

gghistogram(hisdat, "correlation", facet.by=c('v1','v2'), fill="blue") 

