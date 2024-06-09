# Functions to be used for WDI analyses:
#
# wrangle - extracts desired WDI's for each country
# cor2 - obtains CT cor between two WDI's for a given country
# plot.cor2 - plot method for the output of cor2()
# cor2_plus_plot - a wrapper that runs cor2 and makes a plot 
# bootcor2 - computes bootstrap CI for CT cor of two curves

wrangle <- function(incodes, newnames=NULL) {
    wdisome <- subset(wdiall, Indicator.Code %in% incodes)
    wdil <- reshape(wdisome, direction="long", timevar="year", varying=5:67, sep="")
    wdil$coyear <- with(wdil, paste0(Country.Code, year))
    wdiw <- reshape(wdil, direction="wide", drop=c("Indicator.Name","id"), 
                    timevar="Indicator.Code", idvar="coyear", v.names="X")
    rownames(wdiw) <- wdiw$coyear
    wdiw <- wdiw %>% rename(countrycode=Country.Code, countryname=Country.Name) %>%
        mutate(coyear=NULL)
    for (j in 1:length(names(wdiw))) {
        name <- names(wdiw)[j]
        if (substr(name,1,2)=="X.") names(wdiw)[j] <- substr(name,3,nchar(name))
    }
    wdiw
}


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
    #print(allx); print(ally); print(n1); print(n2)
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

# This function computes bootstrap CI of CT correlation.
# It assumes two curves sampled at different time points
# ... such that the two sets of observations must be resampled separately
# m1 and m2 are each two-column matrices, with time in first column and variable of interest in second
bootcor2 <- function(m1, m2, lag=0, min.overlap=0, k.vec=10, nboot=999, conf=0.95) {
    nm1 <- !is.na(m1[,2])
    nm2 <- !is.na(m2[,2])
    n1 <- sum(nm1)
    n2 <- sum(nm2)
    t1 <- m1[nm1,1]+lag
    t2 <- m2[nm2,1]
    if (any(is.na(t1)) | any(is.na(t2))) stop("Missing time points")
    y1 <- m1[nm1,2]
    y2 <- m2[nm2,2]
    if (min(n1,n2)<8) stop("Not enough data")
    
    corvec <- c()
    b <- 1
    while (b <= nboot) {
        ind1 <- sort(sample.int(n1, n1, TRUE))
        ind2 <- sort(sample.int(n2, n2, TRUE))
        t1b <- t1[ind1]
        t2b <- t2[ind2]
        
        if (min(max(t1b),max(t2b))-max(min(t1b),min(t2b)) < min.overlap) {
            warning("Insufficient overlap of time ranges")
            next
        }
        
        y1b <- y1[ind1]
        y2b <- y2[ind2]
        
        modlist1 <- modlist2 <- bsb1 <- bsb2 <- list()
        AIC.vec <- c()
        for (k in 1:length(k.vec)) {
            bsb1[[k]] <- create.bspline.basis(range(t1b), nbasis = k.vec[k])  
            B1 <- eval.basis(t1b, bsb1[[k]])
            P1 <- getbasispenalty(bsb1[[k]])
            modlist1[[k]] <- gam(y1b ~ B1 - 1, paraPen = list(B1 = list(P1)), method = "REML")
            bsb2[[k]] <- create.bspline.basis(range(t2b), nbasis = k.vec[k])  
            B2 <- eval.basis(t2b, bsb2[[k]])
            P2 <- getbasispenalty(bsb2[[k]])
            modlist2[[k]] <- gam(y2b ~ B2 - 1, paraPen = list(B2 = list(P2)), method = "REML")
            AIC.vec[k] <- sum(AIC(modlist1[[k]]), AIC(modlist2[[k]]))
        }
        wmaic <- which.min(AIC.vec) 
        gam1 <- modlist1[[wmaic]]
        gam2 <- modlist2[[wmaic]]
        
        fd1 <- fd(coef = gam1$coef, basisobj = bsb1[[wmaic]])
        fd2 <- fd(coef = gam2$coef, basisobj = bsb2[[wmaic]])
        
        corvec[b] <- cor.ct(fd1, fd2)
        b <- b+1
    }
    list(cors=corvec, 
         lo=sort(corvec)[floor((nboot+1)*(1-conf)/2)], 
         up=sort(corvec)[ceiling((nboot+1)*(1+conf)/2)])
}

