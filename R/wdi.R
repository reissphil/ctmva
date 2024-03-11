# Functions to be used below for World Development Indices analyses

# Obtain country-specific CT correlation between two variables
cor2 <- function(ccode, v1, v2, lag=0, min.overlap=0) {
  require(mgcv); require(ctmva)
  dat <- wlist[[ccode]]
  nm1 <- !is.na(dat[,v1])
  nm2 <- !is.na(dat[,v2])
  n1 <- sum(nm1)
  n2 <- sum(nm2)
  yr1 <- dat$year[nm1]+lag
  yr2 <- dat$year[nm2]
  if (min(n1,n2)<8) stop("Not enough data")
  
  bsb1 <- fda::create.bspline.basis(range(yr1), min(10, n1-2))
  B1 <- fda::eval.basis(yr1, bsb1)
  P1 <- fda::getbasispenalty(bsb1)
  mod1 <- gam(dat[nm1,v1]~B1-1, paraPen=list(B1=list(P1)), method="REML")
  fd1 <- fda::fd(mod1$coef, bsb1)
  
  bsb2 <- fda::create.bspline.basis(range(yr2), min(10, n2-2))
  if (min(bsb1$rangeval[2],bsb2$rangeval[2])-max(bsb1$rangeval[1],bsb2$rangeval[1]) < min.overlap)
    stop("Insufficient overlap of time ranges")
  B2 <- fda::eval.basis(yr2, bsb2)
  P2 <- fda::getbasispenalty(bsb2)
  mod2 <- gam(dat[nm2,v2]~B2-1, paraPen=list(B2=list(P2)), method="REML")
  fd2 <- fda::fd(mod2$coef, bsb2)
  
  res <- list(v1=v1, yr1=yr1, dat.v1=dat[nm1,v1], mod1=mod1, bsb1=bsb1, 
              v2=v2, yr2=yr2, dat.v2=dat[nm2,v2], mod2=mod2, bsb2=bsb2,
              fd1=fd1, fd2=fd2, cor=cor.ct(fd1, fd2))
  class(res) <- "cor2"
  res
}

# Plot method for "cor2" objects
plot.cor2 <- function(x, legend=TRUE, legpos="topright", main=NULL, xlim=NULL, ...) {
  require(viridis)
  cols <- viridis(90)[c(2,40)]
  allx <- c(x$yr1, x$yr2)
  ally <- c(x$dat.v1, x$dat.v2)
  n1 <- length(x$yr1); n2 <- length(x$yr2)
  par(mfrow=2:1)
  par(mar=c(2, 4, 2, 2) + 0.1)
  if (is.null(main)) main <- round(x$cor,3)
  if (is.null(xlim)) xlim <- c(min(c(x$yr1,x$yr2)), max(c(x$yr1,x$yr2)))
  
  grid1 <- seq(min(x$yr1),max(x$yr1),,200)
  X1 <- fda::eval.basis(grid1,x$bsb1)
  pred1 <- X1 %*% x$mod1$coef
  se1 <- sqrt(rowSums((X1 %*% x$mod1$Vc)*X1))
  yrange1 <- c(min(c(x$dat.v1,pred1-2*se1)), max(c(x$dat.v1,pred1+2*se1)))
  plot(x$yr1, x$dat.v1, xlab="Year", ylab=x$v1, xlim=xlim, ylim=yrange1, type="n", main=main, ...)
  polygon(c(grid1,rev(grid1)),c(pred1-2*se1,rev(pred1+2*se1)),col="lightgrey",border=NA)
  points(x$yr1, x$dat.v1, pch=16, col=cols[1], ...)
  lines(grid1, pred1, col=cols[1])
  
  grid2 <- seq(min(x$yr2),max(x$yr2),,200)
  X2 <- fda::eval.basis(grid2,x$bsb2)
  pred2 <- X2 %*% x$mod2$coef
  se2 <- sqrt(rowSums((X2 %*% x$mod2$Vc)*X2))
  yrange2 <- c(min(c(x$dat.v2,pred2-2*se2)), max(c(x$dat.v2,pred2+2*se2)))
  plot(x$yr2, x$dat.v2, xlab="Year", ylab=x$v2, xlim=xlim, ylim=yrange2, type="n", ...)
  polygon(c(grid2,rev(grid2)),c(pred2-2*se2,rev(pred2+2*se2)),col="lightgrey",border=NA)
  points(x$yr2, x$dat.v2, pch=16, col=cols[2], ...)
  lines(grid2, pred2, col=cols[2])
  
  par(mar=c(5, 4, 4, 2) + 0.1)
}

# Run cor2 and plot the result
cor2_plus_plot <- function(ccode, v1, v2, lag=0, min.overlap=0, main=NULL, ...) {
  co <- cor2(ccode, v1, v2, lag, min.overlap)
  if (is.null(main)) main <- paste(ccode,round(co$cor,3))
  plot(co, main=main, ...)
}

#########################
#########################

library(dplyr)
library(ggpubr)

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

# Angola example
cor2_plus_plot("AGO","investment","fem.sec.enroll", main="")

tmp <- cor2("AGO","fem.sec.enroll","investment")
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

