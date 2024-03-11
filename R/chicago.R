# Some analyses of the Chicago air pollution data from the gamair package

require(gamair)
require(mgcv)
require(fda)
require(ctmva)
require(viridis)

data(chicago)

bsb <- create.bspline.basis(range(chicago$time), 200) # create B-spline basis
B <- eval.basis(chicago$time, bsb) # evaluation matrix of the B-spline basis functions
P <- getbasispenalty(bsb) 

# Smooth four variables ignoring autocorrelation
pm10 <- gam(chicago$pm10median ~ B - 1, paraPen=list(B=list(P)), method="REML")
o3 <- gam(chicago$o3median ~ B - 1, paraPen=list(B=list(P)), method="REML")
so2 <- gam(chicago$so2median ~ B - 1, paraPen=list(B=list(P)), method="REML")
temp <- gam(chicago$tmpd ~ B - 1, paraPen=list(B=list(P)), method="REML")

# Smooth four variables with AR(1) residuals 
# (cf. Cairo temperature example in Wood's book)
REML1 <- REML2 <- REML3 <- REML4 <- rho <- 0.1 + 0:70/100
for (i in 1:length(rho)) {
  cat(rho[i], '\n')
  REML1[i] <- bam(chicago$pm10median ~ B - 1, paraPen=list(B=list(P)), rho=rho[i])$gcv.ubre
  REML2[i] <- bam(chicago$o3median ~ B - 1, paraPen=list(B=list(P)), rho=rho[i])$gcv.ubre
  REML3[i] <- bam(chicago$so2median ~ B - 1, paraPen=list(B=list(P)), rho=rho[i])$gcv.ubre
  REML4[i] <- bam(chicago$tmpd ~ B - 1, paraPen=list(B=list(P)), rho=rho[i])$gcv.ubre
}
par(mfrow=c(2,2))
plot(rho, REML1)
plot(rho, REML2)
plot(rho, REML3)
plot(rho, REML4)
pm10ar <- bam(chicago$pm10median ~ B - 1, paraPen=list(B=list(P)), rho=rho[which.min(REML1)])
o3ar <- bam(chicago$o3median ~ B - 1, paraPen=list(B=list(P)), rho=rho[which.min(REML2)])
so2ar <- bam(chicago$so2median ~ B - 1, paraPen=list(B=list(P, sp=exp(14))), rho=0.33) # for larger rho, sp blows up
tempar <- bam(chicago$tmpd ~ B - 1, paraPen=list(B=list(P)), rho=rho[which.min(REML4)])

# Compare fits with and without accounting for autocorrelation
par(mfrow=c(2,2))
matplot(chicago$time[!is.na(chicago$pm10median)], cbind(predict(pm10),predict(pm10ar)), type='l',lty=1, main="PM10")
matplot(chicago$time, cbind(predict(o3),predict(o3ar)), type='l',lty=1, main="Ozone")
matplot(chicago$time[!is.na(chicago$so2)], cbind(predict(so2),predict(so2ar)), type='l',lty=1, main="SO2")
matplot(chicago$time, cbind(predict(temp),predict(tempar)), type='l', lty=1, main="Temperature")

# Define and standardize a functional data object
chifd <- fd(coef=cbind(pm10ar$coef, o3ar$coef,so2ar$coef, tempar$coef),
            basis=bsb, 
            fdnames=list(args="day", vars=c("pm10","o3","so2","temp"), values="value"))
# Standardize the functional data object
schifd <- standardize.ct(chifd)

# Obtain boundaries between 1987 and 1988, ..., 1999 and 2000
yearbreaks <- -2557+365*0:14
yearbreaks[3:14] <- yearbreaks[3:14] + 1   # 1988 is leap year
yearbreaks[7:14] <- yearbreaks[7:14] + 1   # 1992
yearbreaks[11:14] <- yearbreaks[11:14] + 1 # 1996
yearbreaks[15] <- 2557
midyear <- (yearbreaks[-1] + yearbreaks[-15]) / 2

# Plot the four standardized functions 
par(mfrow=c(1,1))
plot(schifd, axes=FALSE, lty=1, xlab="", ylab="Standardized values")
for (j in yearbreaks) abline(v=j,col='grey')
axis(1, at=midyear, labels=1987:2000, lwd.ticks=0)
axis(2)
legend("topright", lty=1, legend=c("PM10", "Ozone", "SO2", "Temperature"), col=1:4)
box()

##################

# PCA
pp <- pca.ct(schifd)
100*pp$var/sum(pp$var)
round(pp$loadings, 3)
start <- min(chicago$time)
end <- max(chicago$time)
gridd <- seq(start, end, , 11200)
sc = eval.fd(gridd, pp$scorefd)
v800 <- c(rev(rocket(400)), rocket(400))

plot(sc[,1],sc[,2], col=rep(v800, 14), pch=15, cex=.5, xlab="First principal component", 
     ylab="Second principal component")
# Add arrows indicating time direction
spotz <- c(20, seq(790,11190,800))
spotz[6] <- spotz[6] - 26
for (j in 1:15) 
  arrows(sc[spotz[j],1], sc[spotz[j],2], sc[spotz[j]+1,1], sc[spotz[j]+1,2], 
         col=v800[spotz[j]%%800], lwd=2, length=.18)
spotz2 <- seq(400,10800,800)
spotz2[2] <- 1176
for (j in 1:14) 
  arrows(sc[spotz2[j],1], sc[spotz2[j],2], sc[spotz2[j]+1,1], sc[spotz2[j]+1,2], 
         col=v800[spotz2[j]%%800], lwd=2, length=.18)
text(sc[800*1:14,1],sc[800*1:14,2], labels =1987:2000)

##################

# 3-means clustering
set.seed(185)
km <- kmeans.ct(schifd,3)   

plot(km, axes = FALSE, xlab="", legend=c("PM10", "Ozone", "SO2", "Temperature"), 
     lty=1, ylab="Standardized values")
abline(h=0, col='grey')
axis(1, at=midyear, labels=1987:2000, lwd.ticks=0)
axis(1, at=yearbreaks, labels=rep("",15))
axis(2)
box()

# Display cluster transitions separately for each year
transits <- clusts <- list()
for (i in 1:14) {
  transits[[i]] <- km$transitions[km$transitions>yearbreaks[i] & km$transitions<yearbreaks[i+1]] - yearbreaks[i]
  wich <- which(km$transitions>yearbreaks[i] & km$transitions<yearbreaks[i+1])
  wich <- c(wich, max(wich)+1)
  clusts[[i]] <- km$cluster[wich]
}

yy <-c(1,14); xx <- c(0,400)
plot(xx,yy, axes=FALSE, type="n", xlab="", ylab="")
colvec <- viridis(3)[c(1,3,2)]  # scico(3, palette="hawaii")  
for (i in 1:14) {
  yearlength <- diff(yearbreaks)[i]
  ntran <- length(transits[[i]])
  tran <- c(transits[[i]], yearlength) * 365/yearlength
  segments(0, 15-i, tran[1], col=colvec[clusts[[i]][1]], lwd=14)
  for (j in 1:ntran) segments(tran[j], 15-i, tran[j+1], col=colvec[clusts[[i]][j+1]], lwd=14)
}
fda::axisIntervals()
axis(2, las=2, at=1:14, labels=2000:1987)
legend("topright", title="Cluster", col=colvec, lty=1, lwd=14, legend=1:3)
box()

# Silhouette width
set.seed(993)
sil <- list()
silvec <- c()
for (k in 2:9) {
  if (k==3) set.seed(185)  # to get same clustering as above
  sil[[k-1]] <- silhouette.ct(kmeans.ct(schifd,k), ngrid=6000)
  silvec[k-1] <- sil[[k-1]]$mean
}

par(mfrow=1:2)
plot(sil[[2]], axes=FALSE, xlab="")  # 3-means
axis(1, at=midyear, labels=1987:2000, lwd.ticks=0)
axis(1, at=yearbreaks, labels=rep("",15))
axis(2)
box()
plot(2:9, silvec, xlab="k", type="o", ylab="Mean silhouette width")

# Check sensitivity to starting value
kms <- list()
ntrans <- twvec <- c()
for (i in 1:50) {
  set.seed(300+i)
  kms[[i]] <- kmeans.ct(schifd,3) 
  ntrans[i] <- length(kms[[i]]$transitions)
  twvec[i] <- kms[[i]]$tot.withinisd
}
plot(ntrans, twvec)  # conclusion: tot.withinisd is essentially fixed for given # transitions (40,54)
1-max(twvec)/kms[[1]]$totisd
1-min(twvec)/kms[[1]]$totisd # conclusion: between/total ISD is about 80%,88% when ntrans is 40,54

# Check variation in cluster transition points among 54-transition solutions, 40-transition solutions
transmat54 <- NULL
for (i in which(ntrans==54)) transmat54 <- cbind(transmat54, kms[[i]]$transitions)
summary(ranges54 <- apply(transmat54, 1, function(v) diff(range(v))))
transmat40 <- NULL
for (i in which(ntrans==40)) transmat40 <- cbind(transmat40, kms[[i]]$transitions)
summary(ranges40 <- apply(transmat40, 1, function(v) diff(range(v))))
# conclusion: each transition occurs consistenly within a <1-day range

# Ordinary 3-means
rawchi <- chicago[,c(2,4,5,7)]
rawchi <- rawchi[apply(rawchi,1,function(v) !any(is.na(v))), ]
set.seed(175)
km3r <- kmeans(rawchi,3)
sum(diff(km3r$cluster)!=0)  # that's a lot of transitions!

##################

# LDA
ll <- lda.ct(schifd, yearbreaks[2:14], part.names=1987:2000)
round(ll$scaling, 4)

plot(ll, axes = FALSE, xlab="", which=1:2, col=mako(8)[c(2,6)])
axis(1, at=midyear, labels=1987:2000, lwd.ticks=0)
axis(2)
box()

ms <- with(ll, means %*% scaling)
cex. <- .7
par(mfrow=1:2)
plot(ms[,1], ms[,3], type="n", xlab="Linear discriminant 1", ylab="Linear discriminant 2")
text(ms[,1], ms[,3], 1987:2000, col=plasma(17)[1:14], cex=cex.)
plot(ll$means[,1], ll$means[,3], type="n", xlab="PM10", ylab="SO2")
text(ll$means[,1], ll$means[,3], 1987:2000, col=plasma(17)[1:14], cex=cex.)
