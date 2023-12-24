saveplots <- FALSE
pdif <- function(filename, ...) if (saveplots) pdf(paste0("~/Dropbox/ctmva/output/", filename), ...)
done <- function() if (saveplots) dev.off()

require(ggplot2)
require(eegkit)
require(ctmva)
require(viridis)

data(eegdata)
data(eegcoord)

# Function to create functional data object for a single trial from the Begleiter EEG data
eeg1 <- function(sub, tri, nbasis=30) {
    require(fda)
    longdat <- subset(eegdata, subject==sub & trial==tri)
    widedat <- reshape(longdat, direction="wide", drop=c("subject","group","condition","trial"),
                       v.names="voltage",idvar="channel")
    bsb <- create.bspline.basis(c(0,255),nbasis)
    return(list(longdat = longdat, widedat = widedat,
                fdobj = Data2fd(argvals=0:255, 
                                y=t(as.matrix(widedat[,-1])), basisobj=bsb)))
}

with(eegdata, table(subject, trial))

fdo <- eeg1("co2a0000369", 0)
plot(ef <- fdo$fdobj)
library(corrplot)
corrplot(cor.ct(ef))
corrplot(cor.ct(ef, common_trend=TRUE))

p = pca.ct(ef)
sc = eval.fd(0:255, p$scorefd)
pairs(sc[,1:5], panel=lines)

cc = pca.ct(ef, common_trend=TRUE)
sc = eval.fd(0:255, cc$scorefd)
pairs(sc[,1:5], panel=lines)

pdif("eeggplot.pdf")
ggplot(fdo$longdat, aes(time, voltage)) + geom_line() + facet_wrap(vars(channel))
dev.off()

pve <- 100*cc$var/sum(cc$var)
pdif("eeg3load.pdf", 13)
par(mfrow=c(1,3))
cidx <- match(fdo$widedat[,1],rownames(eegcoord))
eegspace(eegcoord[cidx,4:5],cc$loadings[,1], colorlab="PC1 loadings", 
         mycolors=cividis(25), main=paste0(round(pve[1],0), "%"), mar=c(17,3,12,2), cex.main=2)
eegspace(eegcoord[cidx,4:5],cc$loadings[,2], colorlab="PC2 loadings", 
         mycolors=cividis(25), main=paste0(round(pve[2],0), "%"), mar=c(17,3,12,2), cex.main=2)
eegspace(eegcoord[cidx,4:5],cc$loadings[,3], colorlab="PC3 loadings", 
         mycolors=cividis(25), main=paste0(round(pve[3],0), "%"), mar=c(17,3,12,2), cex.main=2)
done()
