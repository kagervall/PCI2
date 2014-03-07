##library(ggplot2)

plot.PCI <- function (x, y, z, col = c(1, 1), maxsize = 3, do.sqrt = TRUE,
    main = "", cex.main = 1, xlab = "", ylab = "", minnbubble = 8,
    xlim = NULL, ylim = NULL, axis1 = TRUE, xlabels = TRUE, xlimextra = 1, add = FALSE,
    yadj=.15, las = 1, allopen = TRUE)
    ## Code borrowed from bubble3 in r4ss
{
    az <- abs(z)
    az[az==0] <- .01
    if (do.sqrt)
        az <- sqrt(az)
    cex <- maxsize * az/max(az)
    z.col <- ifelse(z < 0, col[1], col[2])
    if (is.null(xlim)) {558
        xlim <- range(x)
        if (length(unique(x)) < minnbubble)
            xlim = xlim + c(-1, 1) * xlimextra
    }
    pch <- z
    pch[pch == 0] <- 1
    pch[pch > 0] <- 16
    pch[pch < 0] <- 1
    if (allopen)
        pch[!is.na(pch)] <- 1
    if (!add) {
        if (is.null(ylim))
            ylim <- range(y)
        plot(x, y, type = "n", xlim = xlim, ylim = ylim,
            cex.main = cex.main, xlab = xlab, ylab = ylab, axes = FALSE)
        text(x = 3, y = 1.8, labels=main, cex=1.2)
        xvec <- unique(x)
        ## if (axis1)
        ##     axis(1, at = floor(unique(x)), labels=xlabels)
        ## axis(2, las = las)
        box()
    }
    abline(h = 0, col="lightgray")
    lines(x, y, col = z.col)
    points(x, y, pch = pch, cex = cex, col = z.col)
##    yadj <- ylim[2] / cex
##    yadj <- .2
    text(x, y - yadj, cex=cex.main * .8,
         labels=sub(".00", "0", sub("0.", ".",sprintf("%0.2f", z))), col="black")
}


D <- function(rx,ry,d=1){
    if (sign(rx) == sign(ry)) return(0)
    abs(rx - ry) - 2 + d
}

PCI2 <- function(x, scale="bipolar", width=5, distance=1, power=1) {
    min.val <- -(width %/% 2)
    max.val <- (width %/% 2)
   k.levels <- min.val : max.val
   freq <- table(factor(x, levels=k.levels))
   n_t <- sum(freq)
   x.mod <- n_t %% 2
   totdist <- D(min.val, max.val, d=distance)^power * (n_t - x.mod)^2 / 2
    actdist <- 0
    for (k in which(sign(k.levels) > 0)) {
        for (h in which(sign(k.levels) < 0)) {
            actdist <- actdist + 2 * freq[k] * freq[h] *
                D(k.levels[k], k.levels[h], d=distance)^power
        }
    }
    return(list(PCI=actdist/totdist,actdist=actdist,
                totdist=totdist, mean=mean(x), freq=freq))
}

## =IF(MOD($B20;2)=0;
##      (3+$B$4-1)^$B$8*$B20*$B20/2;
##      (3+$B$4-1)^$B$8*($B20+1)*($B20-1)/2)

all.1 <- PCI2(mydata$D8_1)
all.2 <- PCI2(mydata$D8_2)
all.3 <- PCI2(mydata$D8_3)
all.4 <- PCI2(mydata$D8_4)
all.5 <- PCI2(mydata$D8_5)

keep.1 <- PCI2(mydata[mydata$cluster=="Keep",]$D8_1)
keep.2 <- PCI2(mydata[mydata$cluster=="Keep",]$D8_2)
keep.3 <- PCI2(mydata[mydata$cluster=="Keep",]$D8_3)
keep.4 <- PCI2(mydata[mydata$cluster=="Keep",]$D8_4)
keep.5 <- PCI2(mydata[mydata$cluster=="Keep",]$D8_5)
mixed.1 <- PCI2(mydata[mydata$cluster=="Mixed",]$D8_1)
mixed.2 <- PCI2(mydata[mydata$cluster=="Mixed",]$D8_2)
mixed.3 <- PCI2(mydata[mydata$cluster=="Mixed",]$D8_3)
mixed.4 <- PCI2(mydata[mydata$cluster=="Mixed",]$D8_4)
mixed.5 <- PCI2(mydata[mydata$cluster=="Mixed",]$D8_5)
release.1 <- PCI2(mydata[mydata$cluster=="Release",]$D8_1)
release.2 <- PCI2(mydata[mydata$cluster=="Release",]$D8_2)
release.3 <- PCI2(mydata[mydata$cluster=="Release",]$D8_3)
release.4 <- PCI2(mydata[mydata$cluster=="Release",]$D8_4)
release.5 <- PCI2(mydata[mydata$cluster=="Release",]$D8_5)
## acceptance <- data.frame(keep=c(keep.1$mean, mixed.1$mean,release.1$mean),
##                          rel25=c(keep.2$mean, mixed.2$mean,release.2$mean),
##                          rel50=c(keep.3$mean, mixed.3$mean,release.3$mean),
##                          rel75=c(keep.4$mean, mixed.4$mean,release.4$mean),
##                          rel100=c(keep.5$mean, mixed.5$mean,release.5$mean)
## )
## acceptance <- as.data.frame(t(acceptance))
## colnames(acceptance) <- c("Keep","Mixed","Release")
all <- data.frame(mean=c(all.1$mean,all.2$mean,all.3$mean,
                   all.4$mean,all.5$mean),
                   PCI=c(all.1$PCI,all.2$PCI,
                   all.3$PCI,all.4$PCI,all.5$PCI),
                   behavior=c(1:5),
                   clu="All")

keep <- data.frame(mean=c(keep.1$mean,keep.2$mean,keep.3$mean,
                   keep.4$mean,keep.5$mean),
                   PCI=c(keep.1$PCI,keep.2$PCI,
                   keep.3$PCI,keep.4$PCI,keep.5$PCI),
                   behavior=c(1:5),
                   clu="Keep")
mixed <- data.frame(mean=c(mixed.1$mean,mixed.2$mean,mixed.3$mean,
                    mixed.4$mean,mixed.5$mean),
                    PCI=c(mixed.1$PCI,mixed.2$PCI,
                    mixed.3$PCI,mixed.4$PCI,mixed.5$PCI),
                    behavior=c(1:5),
                    clu="Mixed")
release <- data.frame(mean=c(release.1$mean,release.2$mean,release.3$mean,
                      release.4$mean,release.5$mean),
                      PCI=c(release.1$PCI,release.2$PCI,
                      release.3$PCI,release.4$PCI,release.5$PCI),
                      behavior=c(1:5),
                      clu="Release")
acceptance <- rbind(all,keep,mixed,release)
acceptance$psize <- 4 + acceptance$PCI * 100

## qplot(behavior, mean, data=acceptance, color=clu, size=psize)

## qp <- ggplot(acceptance)
## qp <- qp + geom_point(aes(y=mean, x=behavior, color=clu, size=psize))
## print(qp)

 win.metafile(file = "norms-PCI.wmf")
##png(file = "norms-PCI.png", width=6.5, height=6.5, units="in", res=300)
##setEPS()
##postscript(file="norms-PCI.eps")
## pdf(file="norms-PCI.pdf")

lbls <- c(0,25,50,75,100)
##opar <- par(mfrow=c(2,2))
opar <- par(mfrow=c(2,2), omi=rep(0.75, 4), mar=rep(0, 4))
with(all, plot.PCI(behavior, mean, PCI,col=c("darkgray","darkgray"), ylim=c(-2,2),
                   maxsize=5, allopen=FALSE, main="(a) Total sample",
                   xlabels=lbls,  xlab = "Percent",
                   ylab="Acceptance", yadj=c(.37, .35, .36, .37, .41)))
axis(3, at = 1:5, labels=FALSE)
axis(2, las = labels=FALSE)

with(release, plot.PCI(behavior, mean, PCI,col=c("darkgray","darkgray"), ylim=c(-2,2),
                       maxsize=5, allopen=FALSE, main="(b) Release",
                       xlabels=lbls,  xlab = "Percent",
                       ylab="Acceptance", yadj=c(.14, .14, .35, .41, .32)))
axis(3, at = 1:5, labels=lbls)
axis(4, las = 1)

with(mixed, plot.PCI(behavior, mean, PCI,col=c("darkgray","darkgray"), ylim=c(-2,2),
                     maxsize=5, allopen=FALSE, main="(c) Mixed",
                     xlabels=lbls, xlab = "Percent",
                     ylab="Acceptance", yadj=c(.41, .32, .18, .32, .41)))
axis(1, at = 1:5, labels=lbls)
axis(2, las = 1)


with(keep, plot.PCI(behavior, mean, PCI,col=c("darkgray","darkgray"), ylim=c(-2,2), maxsize=5,
                    allopen=FALSE, main="(d) Keep",  xlabels=lbls,  xlab = "Percent",
                    ylab="Acceptance", yadj=c(.4, .4, .35, .18, .14)))
mtext("Percent release", 1, outer=TRUE, padj=2.7)
mtext("Acceptance", 2, outer=TRUE, padj=-2.7)
par(opar)
dev.off()


plot.plain <- function (x, y, main = "",
                        xlab = "", ylab = "",
                        xlim = NULL, ylim = NULL,
                        xlabels = TRUE, add = FALSE,
                        las = 1)
{
    if (is.null(ylim))
        ylim <- range(y)
    plot(x, y, type = "n", xlim = xlim, ylim = ylim,# main = main,
         xlab = xlab, ylab = ylab, axes = FALSE)
    text(x = 3, y = 1.8, labels=main, cex=1.5)
    xvec <- unique(x)
##    if (axis1)
    axis(1, at = floor(unique(x)), labels=xlabels)
    axis(2, las = las)
    box()
    abline(h = 0)
    lines(x, y, type="b", lwd=2)
}

## Plot simple graph without PCI bubbles
 win.metafile(file = "norms-noPCI.wmf")
##setEPS()
##postscript(file="norms-noPCI.eps")
## pdf(file="norms-noPCI.pdf")

opar <- par(mfrow=c(2,2), mar=c(4,4,.5,2) + .1 )

with(all, plot.plain(behavior, mean, ylim=c(-2,2),
                     main="(a) Total sample",
                     xlabels=lbls,  xlab = "Release (%)",
                     ylab="Acceptance"))
with(release, plot.plain(behavior, mean, ylim=c(-2,2),
                     main="(b) Release",
                     xlabels=lbls,  xlab = "Release (%)",
                     ylab="Acceptance"))
with(mixed, plot.plain(behavior, mean, ylim=c(-2,2),
                     main="(c) Mixed",
                     xlabels=lbls,  xlab = "Release (%)",
                     ylab="Acceptance"))
with(keep, plot.plain(behavior, mean, ylim=c(-2,2),
                     main="(d) Keep",
                     xlabels=lbls,  xlab = "Release (%)",
                     ylab="Acceptance"))


par(opar)
dev.off()
