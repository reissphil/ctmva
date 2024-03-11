require(mgcv)
require(fda)
require(ctmva)
require(viridis)
require(ggplot2)
require(corrplot)

daybasis <- create.fourier.basis(c(0,365), nbasis=45)  

cwdat <- data.frame(tmp=1:35)
cwdat$temp <-t(CanadianWeather$dailyAv[,,"Temperature.C"])
cwdat$prec <- t(CanadianWeather$dailyAv[,,"log10precip"])

# Use intercept-only pffr to diagnose need for detrending
require(refund)
tmod<-pffr(temp~1, data=cwdat, bs.int=list(bs="cc", k=20, m=c(2, 1)))
pmod <- pffr(prec~1, data=cwdat, bs.int=list(bs="cc", k=20, m=c(2, 1)))
summary(tmod)$r.sq
summary(pmod)$r.sq

###############################

# Define fd objects, with lambda chosen by REML (see old versions of this file for GCV)

require(gamair)
data(canWeather)
cw <- CanWeather 
# Given coordinates are really degrees and minutes 
# (except Prince George)
# Function to fix this
min2dec <- Vectorize(function(coo) {
  minpart <- coo-floor(coo)
  floor(coo) + 5*minpart/3
})
CanadianWeather$coordinates[5, ] <-c(46.24,63.13) # fix Charlottetown
for (j in setdiff(1:35,c(5,28))) 
  CanadianWeather$coordinates[j, ] <- min2dec(CanadianWeather$coordinates[j, ])
cw$longitude <- rep(CanadianWeather$coordinates[,2], each=365)
cw$latitude <- rep(CanadianWeather$coordinates[,1], each=365)
cw$log10precip <- as.vector(CanadianWeather$dailyAv[,,"log10precip"])

cw$coast <- "Inland"
cw$coast[1:(365*6)] <- "Atlantic"
cw$coast[cw$place %in% c("Vancouver","Victoria","Pr. Rupert")] <- "Pacific"

############################

B <- eval.basis(day.5, daybasis) 
P <- getbasispenalty(daybasis) 
coef.temp <- coef.logprec <- matrix(NA, daybasis$nbasis, 35)
pred.temp <- pred.logprec <-c()
for (i in 1:35) {
  tempfit <- gam(CanadianWeather$dailyAv[,i,"Temperature.C"] ~ B - 1,
                 paraPen=list(B=list(P)), method="REML")
  precfit <- gam(CanadianWeather$dailyAv[,i,"log10precip"] ~ B - 1,
                 paraPen=list(B=list(P)), method="REML")
  coef.temp[,i] <- tempfit$coef
  coef.logprec[,i] <- precfit$coef
  pred.temp <- c(pred.temp, predict(tempfit))
  pred.logprec <- c(pred.logprec, predict(precfit))
}

cw$tempfit <- pred.temp
cw$precfit <- pred.logprec
fd_names <- list(time=day.5, reps=CanadianWeather$place, values="value")
tempfd <- fd(coef=coef.temp, basisobj=daybasis, fdnames=fd_names)
precfd <- fd(coef=coef.logprec, basisobj=daybasis, fdnames=fd_names)

#################################

# Plot raw and smoothed temperature and log precipitation
wplot <- function(y, colvar, ylab, cold=NULL) {
  ggplot(cw, aes(time,{{y}})) +
    geom_line(aes(time, {{y}}, group=place, color={{colvar}})) +
    {if (is.null(cold)) scale_color_gradientn(colors=rev(viridis(40)))} +
    {if (!is.null(cold)) scale_color_discrete(type=cold)} +
    scale_x_time(breaks = fda::monthMid+0.5, 
                 labels= c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                           "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
    labs(y=ylab, x="")
}

wplot(T, longitude, "Mean temperature")
wplot(tempfit, longitude, "Mean temperature")

aip <- c("blue4","darkolivegreen3","turquoise2")
wplot(log10precip, coast, "log10 of mean precipitation in mm", cold=aip)
wplot(precfit, coast, "log10 of mean precipitation in mm", cold=aip)

##################

# Map of the 35 stations, adapted from Scheipl et al., JCGS, 2015
require(maps)
mapcanada <- map(database="world", regions="can", resolution=1, plot=TRUE)

plot(mapcanada, type="l", xaxt="n", yaxt="n", ylab="", xlab="", bty="n",
     xlim = c(-141,-52), 
     ylim=c(40,78), 
     mar=c(0,0,0,0))
colvec <- rep(aip[2], 35)
colvec[1:6] <- aip[1]
colvec[CanadianWeather$place %in% c("Vancouver","Victoria","Pr. Rupert")] <- aip[3]
with(CanadianWeather, points(-coordinates[,2], coordinates[,1], 
                             pch = 16, col = colvec, cex = 1.2))
legend("bottomleft", pch=16, col=aip,
       legend=c("Atlantic", "Inland", "Pacific"))

################################################

# Correlation matrix for temperature
# Without versus with detrending
rawcor <- cor.ct(tempfd, common_trend = FALSE)
dtcor  <- cor.ct(tempfd, common_trend = TRUE)

par(mfrow=1:2)
corrplot(rawcor, method = 'square', tl.col="black", tl.cex = 0.6)
corrplot(dtcor, method = 'square', tl.col="black", tl.cex = 0.6)

# Correlation matrix for log precipitation
# Raw-data versus CT correlations, with reordering
par(mfrow=1:2)
rco <-cor(CanadianWeather$dailyAv[,,"log10precip"])
ood <- corrMatOrder(rco)
corrplot(rco[rev(ood),rev(ood)], tl.col="black", tl.cex = 0.6)
pco <- cor.ct(precfd)
ord <- corrMatOrder(pco)
corrplot(pco[ord,ord], tl.col="black", tl.cex = 0.6)
