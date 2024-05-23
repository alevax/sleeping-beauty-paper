# Functions for gsea analysis (From Mariano's atools -- without DESeq)
# TODO
# Add ledge

#' Gene Set Enrichment Analysis
#' 
#' This function performs gene eset enrichment analysis and plots
#' 
#' @param signature Numeric vector containing the gene expression signature
#' @param geneset Vector of character strings indicating the genes in the geneset, or named numeric vector for 2-tail analysis
#' @param score Number indicating the exponent score for GSEA
#' @param twoTails Logical, whether the ES should be computed as the different between each tail ES
#' @param pout Logical, whether a plot should be generated
#' @param per Integer indicating the number of permutations for p-value estimation
#' @param alternative Character string indicating the tail to test for statistical significance
#' @param colSig Vector indicating the colors for the signature, including hue for negative values, hue for positive values, mate value and gama
#' @param colHit Vector indicating the colors for the hits, including hue for negative values, hue for positive values, mate value and gama
#' @param ylim Optional numeric vector indicating the limits for the y-axis
#' @param axes Logical, whether axis should be ploted
#' @param xlab Character string indicating the label for the x-axis
#' @param ylab Character string indicating the label for the y-axis
#' @param lwd Number indicating the line width for the enrichment plot
#' @param maxhit Integer indicating the maximum number of hits to plot, 0 means all
#' @param ... Additional parameters to pass to plot function
#' @return List of results and plot
#' @export
gsea <- function(signature, geneset, score=1, twoTails=FALSE, pout=TRUE, per=0, nesnull=NULL, alternative=c("two.sided", "greater", "less"), colSig=c(.45, .15, .3, 1), colHit=c(.58, .05, .1, 2), ylim=NULL, axes=TRUE, xlab="Signature", ylab="ES", lwd=1, maxhit=0, ...) {
  alternative <- match.arg(alternative)
  if (is.null(names(geneset))) {
    tmp <- geneset
    geneset <- rep(1, length(geneset))
    names(geneset) <- tmp
  }
  if (!is.null(ncol(signature))) signature <- signature[, 1]
  if (prod(range(signature))>=0) { #One tail enrichment
    twoTails <- FALSE
    signature <- sort(abs(signature), decreasing=TRUE)
    es1 <- gsea.es(signature, which(names(signature) %in% names(geneset)), score)
    pos <- which.max(abs(es1))
    es <- es1[pos]
    # LEDGE
    if (es<0) ledgesig <- signature[pos:length(signature)]
    else ledgesig <- rev(signature[1:pos])
    ledge <- ledgesig[names(ledgesig) %in% names(geneset)] 
    nes <- list(nes=NULL, p.value=NULL)
    if (per>0 || !is.null(nesnull)) {
      if (is.null(nesnull)){
        nesnull <- sapply(1:per, function(i, signature, setsize, score) {
          es1 <- gsea.es(signature, sample(length(signature), setsize), score)
          es1[which.max(abs(es1))]
        }, signature=signature, setsize=length(which(names(signature) %in% names(geneset))), score=score)
      }
      nes <- aecdf(nesnull, symmetric=TRUE)(es, alternative)
    }
    if (pout) {
      if (is.null(ylim)) ylim <- c(min(0, min(es1)), max(0, max(es1)))
      if (length(ylim)==1) ylim <- c(-abs(ylim), abs(ylim))
      wd <- diff(ylim)
      base.bar <- ylim[2]+wd*.05
      top.bar <- ylim[2]+wd*.15
      base.sig <- ylim[1]-wd*.15
      top.sig <- ylim[1]-wd*.05
      plot(0, 0, type="n", xlim=c(0, length(signature)+1), ylim=c(base.sig, top.bar), axes=axes, xlab=xlab, ylab=ylab, ...)
      abline(h=0, col="grey")
      pos <- which.max(abs(es1))
      lines(c(pos, pos), c(0, es1[pos]), col="grey")    
      lines(es1, lwd=2)
      addDensityBar(which(names(signature) %in% names(geneset)), length(signature), col=colHit[1], base.bar=base.bar, top.bar=top.bar, mate=colHit[3], gama=colHit[4], lwd=lwd, maxhit=maxhit)
      addSigBar(signature, colSig[1], base.sig, top.sig, colSig[3], colSig[4])
      if (!is.null(nes$nes)) {
        axis(3, length(signature)*.7, paste("NES: ", round(nes$nes, 2), ", p:", sep=""), hadj=1, tick=FALSE)
        axis(3, length(signature)*.71, niceExponent(nes$p.value), hadj=0, tick=FALSE)
      }
    }
    return(list(es=es, nes=nes$nes, p.value=nes$p.value, ledge=ledge, ledgesig=ledgesig))
  }    
  if (prod(range(geneset))>=0) { #1-tail gsea on a +- signature
    signature <- sort(signature)
    es1 <- -gsea.es(signature, which(names(signature) %in% names(geneset)), score)
    # Check whether the analysis is 1 or 2 tails !!!!
    if (twoTails) {
      pos <- sort(c(which.min(es1), which.max(es1)))
      es <- diff(range(es1))
      # LEDGE
      ledgesig <- signature[c(1:pos[1], pos[2]:length(signature))]
      ledge <- ledgesig[names(ledgesig) %in% names(geneset)] 
      nes <- list(nes=NULL, p.value=NULL)
      if (per>0 || !is.null(nesnull)) {
        if (is.null(nesnull)){
          nesnull <- sapply(1:per, function(i, signature, setsize, score) {
            es1 <- -gsea.es(signature, sample(length(signature), setsize), score)
            diff(range(es1))
          }, signature=signature, setsize=length(which(names(signature) %in% names(geneset))), score=score)
        }
        nes <- aecdf(nesnull, symmetric=TRUE)(es, alternative)
      }
      if (pout) {
        if (is.null(ylim)) ylim <- c(min(0, min(es1)), max(0, max(es1)))
        if (length(ylim)==1) ylim <- c(-abs(ylim), abs(ylim))
        wd <- diff(ylim)
        base.bar <- ylim[2]+wd*.05
        top.bar <- ylim[2]+wd*.15
        base.sig <- ylim[1]-wd*.15
        top.sig <- ylim[1]-wd*.05
        plot(0, 0, type="n", xlim=c(0, length(signature)+1), ylim=c(base.sig, top.bar), axes=axes, xlab=xlab, ylab=ylab, ...)
        abline(h=0, col="grey")
        pos <- which.max(es1)
        lines(c(pos, pos), c(0, es1[pos]), col="grey")
        pos <- which.min(es1)
        lines(c(pos, pos), c(0, es1[pos]), col="grey")
        lines(es1, lwd=2)
        addDensityBar(which(names(signature) %in% names(geneset)), length(signature), col=colHit[1], base.bar=base.bar, top.bar=top.bar, mate=colHit[3], gama=colHit[4], lwd=lwd, maxhit=maxhit)
        addSigBar(signature, colSig[1:2], base.sig, top.sig, colSig[3], colSig[4])
        if (!is.null(nes$nes)) {
          axis(3, length(signature)*.7, paste("NES: ", round(nes$nes, 2), ", p:", sep=""), hadj=1, tick=FALSE)
          axis(3, length(signature)*.71, niceExponent(nes$p.value), hadj=0, tick=FALSE)
        }
      }
      return(list(es=es, nes=nes$nes, p.value=nes$p.value, ledge=ledge, ledgesig=ledgesig))
    }
    else {
      pos <- which.max(abs(es1))
      es <- es1[pos]
      # ledge
      if (es<0) ledgesig <- signature[1:pos]
      else ledgesig <- rev(signature[pos:length(signature)])
      ledge <- ledgesig[names(ledgesig) %in% names(geneset)]
      nes <- list(nes=NULL, p.value=NULL)
      if (per>0 || !is.null(nesnull)) {
        if (is.null(nesnull)){
          nesnull <- sapply(1:per, function(i, signature, setsize, score) {
            es1 <- -gsea.es(signature, sample(length(signature), setsize), score)
            es1[which.max(abs(es1))]
          }, signature=signature, setsize=length(which(names(signature) %in% names(geneset))), score=score)
        }
        nes <- aecdf(nesnull, symmetric=TRUE)(es, alternative)
      }
      if (pout) {
        if (is.null(ylim)) ylim <- c(min(0, min(es1)), max(0, max(es1)))
        if (length(ylim)==1) ylim <- c(-abs(ylim), abs(ylim))
        wd <- diff(ylim)
        base.bar <- ylim[2]+wd*.05
        top.bar <- ylim[2]+wd*.15
        base.sig <- ylim[1]-wd*.15
        top.sig <- ylim[1]-wd*.05
        plot(0, 0, type="n", xlim=c(0, length(signature)+1), ylim=c(base.sig, top.bar), axes=axes, xlab=xlab, ylab=ylab, ...)
        abline(h=0, col="grey")
        pos <- which.max(abs(es1))
        lines(c(pos, pos), c(0, es1[pos]), col="grey")
        lines(es1, lwd=2)
        addDensityBar(which(names(signature) %in% names(geneset)), length(signature), col=colHit[1], base.bar=base.bar, top.bar=top.bar, mate=colHit[3], gama=colHit[4], lwd=lwd, maxhit=maxhit)
        addSigBar(signature, colSig[1:2], base.sig, top.sig, colSig[3], colSig[4])
        if (!is.null(nes$nes)) {
          axis(3, length(signature)*.7, paste("NES: ", round(nes$nes, 2), ", p:", sep=""), hadj=1, tick=FALSE)
          axis(3, length(signature)*.71, niceExponent(nes$p.value), hadj=0, tick=FALSE)
        }
      }
      return(list(es=es, nes=nes$nes, p.value=nes$p.value, ledge=ledge, ledgesig=ledgesig))
    }
  }
  # 2-tails GSEA
  es1 <- gsea2.es(signature, geneset, score)
  pos <- which.max(abs(es1$es))
  es <- es1$es[pos]
  # ledge
  if (es>0) {
    tmp <- es1$rlist[1:pos]
    ledge <- tmp[names(tmp) %in% names(geneset)]
    ledge <- ledge * sign(geneset)[match(names(ledge), names(geneset))]
    ledgesig=1
    #        ledgesig <- signature[c(1:max(which(names(signature) %in% names(ledge[ledge<0]))), min(which(names(signature) %in% names(ledge[ledge>0]))):length(signature))]
  }
  else {
    tmp <- rev(es1$rlist[pos:length(es1$rlist)])
    ledge <- tmp[names(tmp) %in% names(geneset)]
    ledge <- ledge * sign(geneset)[match(names(ledge), names(geneset))]
    ledgesig=1
    #        ledgesig <- signature[c(1:max(which(names(signature) %in% names(ledge[ledge>0]))), min(which(names(signature) %in% names(ledge[ledge<0]))):length(signature))]
  }   
  nes <- list(nes=NULL, p.value=NULL)
  
  if (per>0 || !is.null(nesnull)) {
    if (is.null(nesnull)){
      nesnull <- sapply(1:per, function(i, signature, geneset, score) {
        names(geneset) <- sample(names(signature), length(geneset))
        es1 <- gsea2.es(signature, geneset, score)
        es1$es[which.max(abs(es1$es))]
      }, signature=signature, geneset=geneset, score=score)
    }
    nes <- aecdf(nesnull, symmetric=TRUE)(es, alternative)
  }
  if (pout) {
    signature <- sort(signature)
    es1 <- gsea.es(signature, which(names(signature) %in% names(geneset[geneset<0])), score)
    es2 <- gsea.es(signature, which(names(signature) %in% names(geneset[geneset>=0])), score)
    if (is.null(ylim)) ylim <- range(c(es1, es2))
    if (length(ylim)==1) ylim <- c(-abs(ylim), abs(ylim))
    wd <- diff(ylim)
    base.sig <- ylim[1]-wd*.15
    top.sig <- ylim[1]-wd*.05
    plot(0, 0, type="n", xlim=c(0, length(signature)+1), ylim=c(base.sig, ylim[2]+wd*.25), axes=axes, xlab=xlab, ylab=ylab, ...)
    abline(h=0, col="grey")
    pos <- which.max(abs(es1))
    if (colHit[1]>1) {
      lines(c(pos, pos), c(0, es1[pos]), col="grey")
      lines(es1, lwd=2)
    }
    else {
      lines(c(pos, pos), c(0, es1[pos]), col=hsv(colHit[1], .5, 1-colHit[3]))
      lines(es1, lwd=2, col=hsv(colHit[1], 1, 1-colHit[3]))
    }
    pos <- which.max(abs(es2))
    if (colHit[2]>1) {
      lines(c(pos, pos), c(0, es2[pos]), col="grey")
      lines(es2, lwd=2)
    }
    else {
      lines(c(pos, pos), c(0, es2[pos]), col=hsv(colHit[2], .5, 1-colHit[3]))
      lines(es2, lwd=2, col=hsv(colHit[2], 1, 1-colHit[3]))
    }
    base.bar <- ylim[2]+wd*.05
    top.bar <- ylim[2]+wd*.15
    addDensityBar(which(names(signature) %in% names(geneset[geneset<0])), length(signature), col=colHit[1], base.bar=base.bar, top.bar=top.bar, mate=colHit[3], gama=colHit[4], lwd=lwd, maxhit=maxhit)
    base.bar <- ylim[2]+wd*.15
    top.bar <- ylim[2]+wd*.25
    addDensityBar(which(names(signature) %in% names(geneset[geneset>=0])), length(signature), col=colHit[2], base.bar=base.bar, top.bar=top.bar, mate=colHit[3], gama=colHit[4], lwd=lwd, maxhit=maxhit)
    addSigBar(signature, colSig[1:2], base.sig, top.sig, colSig[3], colSig[4])
    if (!is.null(nes$nes)) {
      axis(3, length(signature)*.7, paste("NES: ", round(nes$nes, 2), ", p:", sep=""), hadj=1, tick=FALSE)
      axis(3, length(signature)*.71, niceExponent(nes$p.value), hadj=0, tick=FALSE)
    }
  }
  return(list(es=es, nes=nes$nes, p.value=nes$p.value, ledge=ledge, ledgesig=ledgesig))
}

gsea.es <- function(rlist, x, score) {
  nr <- sum(abs(rlist[x])^score)
  nh <- length(rlist)-length(x)
  es <- rep(-(1/nh),length(rlist))
  es[x] <- abs(rlist[x])^score/nr
  return(cumsum(es))
}

gsea2.es <- function(rlist, x, score) {
  x1 <- names(x)[x>=0]
  x2 <- names(x)[x<0]
  px1 <- match(x1, names(sort(rlist, decreasing=T)))
  px2 <- match(x2, names(sort(rlist, decreasing=F)))
  names(px1) <- x1
  names(px2) <- x2
  rlistr <- rank(rlist)
  rlistr <- rlistr[!(names(rlist) %in% names(x))]
  rlistr <- sort(c(rlistr, px1, px2))
  rlist <- rlist[match(names(rlistr), names(rlist))]
  x <- which(names(rlist) %in% names(x))
  nr <- sum(abs(rlist[x])^score)
  nh <- length(rlist)-length(x)
  es <- rep(-(1/nh),length(rlist))
  es[x] <- abs(rlist[x])^score/nr
  return(list(es=cumsum(es), rlist=rlist))
}

addDensityBar <- function(x, sigLength, col, base.bar, top.bar, mate=.2, gama=2, lwd=1, maxhit=0) {
  den <- density(x, adj=.1)
  if (maxhit>0 & maxhit<length(x)) {
    x <- x[unique(round(seq(1, length(x), length=maxhit)))]
  }
  den <- approx(den, xout=x)
  den$y <- den$y/max(den$y)
  den <- lapply(den, function(x, pos) x[pos], pos=order(den$y))
  if (col>1) col <- hsv(0, 0, (1-den$y^gama))
  else col <- hsv(col, den$y^gama, 1-mate*den$y)
  for (i in 1:length(x)) lines(c(den$x[i], den$x[i]), c(base.bar, top.bar), col=col[i], lwd=lwd)
  rect(0, base.bar, sigLength+1, top.bar)
}

addSigBar <- function(signature, col, base.sig, top.sig, mate=.1, gama, res=1000) {
  x <- signature/max(abs(signature))
  if (length(col)==1) col <- rep(col, 2)
  col1 <- rep(hsv(0, 0, 0), length(x))
  pos <- x<0
  if (col[1]>1) col1[pos] <- hsv(0, 0, (1-abs(x[pos])^gama))
  else col1[pos] <- hsv(col[1], abs(x[pos])^gama, 1-mate*abs(x[pos]))
  if (col[2]>1) col1[!pos] <- hsv(0, 0, (1-x[!pos])^gama)
  else col1[!pos] <- hsv(col[2], x[!pos]^gama, 1-mate*x[!pos])
  for (i in unique(round(seq(1, length(x), length=res)))) lines(c(i, i), c(base.sig, top.sig), col=col1[i])
  rect(0, base.sig, length(x)+1, top.sig)
}

#' Approximate empirical commulative distribution function
#'
#' This function generates an empirical null model that computes a normalized statistics and p-value
#' 
#' @param dnull Numerical vector representing the null model
#' @param symmetric Logical, whether the distribution should betreated as symmetric around zero and only one tail should be approximated
#' @param n Integer indicating the number of points to evaluate the empirical cummulative probability function
#' @return function with two parameters, \code{x} and \code{alternative}
#' @export

aecdf <- function(dnull, symmetric=FALSE, n=100) {
  dnull <- dnull[is.finite(dnull)]
  if (symmetric) {
    tmp <- sort(abs(dnull), decreasing=T)
    i <- 4
    n <- 4
    while(n<14) {
      i <- i+1
      n <- length(unique(tmp[1:i]))
      if (n==5) iq1 <- i
    }
    tl1 <- i
    iqr <- quantile(abs(dnull), c(.5, 1-iq1/length(dnull)))
    epd <- ecdf(abs(dnull))
    a <- list(x=knots(epd), y=epd(knots(epd)))
    fit <- lm(y~0+x, data=list(x=a$x[length(a$x)-(tl1:iq1)+1]-iqr[2], y=log(1-epd(iqr[2]))-log(1-a$y[length(a$x)-(tl1:iq1)+1])))
    val <- seq(0, iqr[2], length=n)
    pd <- approxfun(val, epd(val), method="linear", yleft=0, rule=2)
    dnull <- function(x, alternative=c("two.sided", "greater", "less")) {
      alternative <- match.arg(alternative)
      x1 <- abs(x)
      p <- exp(log(1-pd(iqr[2]))-predict(fit, list(x=x1-iqr[2])))
      p[!is.finite(p)] <- 1
      p <- p * (x1>iqr[2]) + (1-pd(x1)) * (x1<=iqr[2])
      nes <- qnorm(p/2, lower.tail=F)*sign(x)
      switch(alternative,
             two.sided={p <- p},
             greater={p <- p/2; p[x<0] <- 1-p[x<0]},
             less={p <- p/2; p[x>0] <- 1-p[x>0]}
      )
      names(nes) <- names(p) <- names(x)
      list(nes=nes, p.value=p)
    }
    return(dnull)
  }
  tmp <- sort(dnull, decreasing=FALSE)
  i <- 4
  n <- 4
  while(n<14) {
    i <- i+1
    n <- length(unique(tmp[1:i]))
    if (n==5) iq1 <- i
  }
  tl1 <- i
  tmp <- sort(dnull, decreasing=TRUE)
  i <- 4
  n <- 4
  while(n<14) {
    i <- i+1
    n <- length(unique(tmp[1:i]))
    if (n==5) iq2 <- i
  }
  tl2 <- i
  iqr <- quantile(dnull, c(iq1/length(dnull), .5, 1-iq2/length(dnull)))
  epd <- ecdf(dnull)
  a <- list(x=knots(epd), y=epd(knots(epd)))
  fit1 <- lm(y~0+x, data=list(x=a$x[iq1:tl1]-iqr[1], y=log(epd(iqr[1]))-log(a$y[iq1:tl1])))
  fit2 <- lm(y~0+x, data=list(x=a$x[length(a$x)-(tl2:iq2)+1]-iqr[3], y=log(1-epd(iqr[3]))-log(1-a$y[length(a$x)-(tl2:iq2)+1])))
  val <- seq(iqr[1], iqr[3], length=n)
  pd <- approxfun(val, epd(val), method="linear", rule=2)
  dnull <- function(x, alternative=c("two.sided", "greater", "less")) {
    alternative <- match.arg(alternative)
    p1 <- exp(log(pd(iqr[1]))-predict(fit1, list(x=x-iqr[1])))
    p2 <- exp(log(1-pd(iqr[3]))-predict(fit2, list(x=x-iqr[3])))
    p1[!is.finite(p1)] <- 1
    p2[!is.finite(p2)] <- 1
    p <- p1*(x<iqr[1]) + p2*(x>iqr[3]) + pd(x)*(x>=iqr[1] & x<iqr[2]) + (1-pd(x))*(x>=iqr[2] & x<=iqr[3])
    nes <- qnorm(p, lower.tail=F)*sign(x-iqr[2])
    switch(alternative,
           two.sided={p <- p*2},
           greater={p[x<iqr[2]] <- 1-p[x<iqr[2]]},
           less={p[x>=iqr[2]] <- 1-p[x>=iqr[2]]}
    )
    names(nes) <- names(p) <- names(x)
    list(nes=nes, p.value=p)
  }
  return(dnull)
}

#' Nice Exponential representations of scientific notation
#' 
#' This function generates a plotmath or latex representation of scientific notation
#' 
#' @param x Numeric vector
#' @param drop.1 Logical, whether 1 in 1 x type of representatons should be dropped
#' @param sub10 Either logical, "10", a non-negative integer or a length 2 integer vector, indicating if some expression should be formatted traditionally, when integer, all expression before the integer are simplified. when a 2 elements vector, all between the indicated range are simplified
#' @param digits Number of significant digits
#' @param lab.type Character string indicating how the result should look like, either plotmath or latex
#' @param lab.sep Character separator between mantissa and exponent
#' @return Vector of formated numbers
#' @export
niceExponent <- function(x, drop.1 = TRUE, sub10 = "10", digits = 2, digits.fuzz, lab.type = c("plotmath", "latex"), lab.sep = c("cdot", "times"))
{
  lab.type <- match.arg(lab.type)
  lab.sep <- match.arg(lab.sep)    
  eT <- floor(log10(abs(x)) + 10^-digits)
  mT <- signif(x / 10^eT, digits)
  ss <- vector("list", length(x))
  if(sub.10 <- !identical(sub10, FALSE)) {
    if(identical(sub10, TRUE))
      sub10 <- c(0,0)
    else if(identical(sub10, "10"))
      sub10 <- 0:1
    sub10 <- as.integer(sub10)
    noE <-
      if(length(sub10) == 1) {
        if(sub10 < 0)
          stop("'sub10' must not be negative if a single number")
        eT <= sub10
      } else if(length(sub10) == 2) {
        stopifnot(sub10[1] <= sub10[2])
        sub10[1] <= eT & eT <= sub10[2]
      } else stop("invalid 'sub10'")
    mT[noE] <- mT[noE] * 10^eT[noE]
  }
  if (lab.type == "plotmath") {
    for(i in seq(along = x))
      ss[[i]] <-
        if(x[i] == 0) quote(0)
    else if(sub.10 &&  noE[i]    ) substitute( A, list(A = mT[i]))
    else if(drop.1 && mT[i] ==  1) substitute( 10^E, list(E = eT[i]))
    else if(drop.1 && mT[i] == -1) substitute(-10^E, list(E = eT[i]))
    else substitute(A %*% 10^E, list(A = mT[i], E = eT[i]))
    do.call("expression", ss)
  }
  else { 
    mTf <- format(mT)
    eTf <- format(eT)
    for(i in seq(along = x))
      ss[[i]] <-
      if(x[i] == 0) ""
    else if(sub.10 &&  noE[i]    ) mTf[i]
    else if(drop.1 && mT[i] ==  1) sprintf("$10^{%s}$", eTf[i])
    else if(drop.1 && mT[i] == -1) sprintf("$-10^{%s}$",eTf[i])
    else sprintf("$%s \\%s 10^{%s}$", mTf[i], lab.sep,  eTf[i])
    ss  ## perhaps unlist(ss) ?
  }
}
