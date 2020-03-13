#!/usr/local/bin/Rscript

# functions.R -  R script
# Copyright (C) 2018  Sebastian Mackowiak
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.



## reverse function to logB('x') with B=base
## pfun(x,2) => 2^x
## useful to apply to lists with lapply(list,pfun,base=2)
pfun <- function(x,base){base^x}

## remove this later
.libPaths( c( "/scratch/local2/R-3.6.0/lib64/R/library" , .libPaths() ) )

zz=file("/dev/null", open = "wt")
#sink(zz,type="message")
require("data.table")

nn <-function(x){
	  return(format(x, big.mark=",", scientific=FALSE,trim=TRUE))
}


hm <- function(x,...){
	require(gplots)
	heatmap.2(x,Rowv=F,Colv=F,dendrogram="none",trace="none",...)
}

id2num <- function(x,lrange=1:4,size=48){
    abc <- letters[lrange];
    v <- substr(x,0,1);
    w<-as.integer(substr(x,2,3));
    return ((which(abc %in%v)-1)*size+w)
}
num2id <- function(x,lrange=1:4,size=48){abc <- letters[lrange];v <- ceiling(x/size) ;w<-x%%size; return (paste(abc[v],w,sep=""))}


## write table function
## rownames=F,colnames=F,quote=F,variable name will be written to file variable_name.csv
wt <- function(f,q=F,r=F,c=F){
    ## rownames=F,colnames=F,quote=F,variable name will be written to file variable_name.csv
    write.table(f,file=paste(deparse(substitute(f)),"csv",sep="."),quote=q,row.names=r,col.names=c,sep="\t")
}


#quantile normalization of matrix
mqq <-function(m){
    mo=apply(m,2,order)
    dm=dim(m)[2]
    rmm=matrix(rep(rowMeans(apply(m,2,sort)),dm),ncol=dm)
    qqm=matrix(rmm[order(mo)],ncol=dm,byrow=T)
    qqm
}

reordermatrix <- function(orig,reord,rown=FALSE,h=0){
	if(h==1){
		print("reordermatrix(orig,reord)")
		print("the function sorts the rows of matrix reord to be in the same order as rows of matrix orig")
		return(0)
	}
	x=rownames(reord)
	y=rownames(orig)
        if(rown){
            x=reord[,rown]
            y=orig[,rown]
        }else{
            print("rown = FALSE")
        }



	keepx=c()

	## we go over the rownames of l2 from 1:length(l2) and get the index in l1 that has the same name as l2
	## so we resort l1 to have the same order as l2
	for(i in 1:length(y)){
		rn=which(x%in%y[i])
		if(length(rn) == 1){
			keepx=c(keepx,rn)
		}else{
			print(paste("name",y[i],"not found in matrix 1\n"))
			keepx=c(keepx,0)
		}
	}
	return(keepx)
}

## open a new plotting window with size different from default 7,7
library(repr)
d1 <- function(w=10,h=10){
	options(repr.plot.width=w, repr.plot.height=h)
}


d <- function(w=10,h=10){
	dev.new(width=w, height=h)
}

## close device , shortcut for dev.off()
do <- function(x=dev.cur()){
	if(x != 1){
		dev.off(x)
	}else{
		print("only null device is open")
	}
}

## this function is a negative plot function just giving a white window by default
## add parameters to modify what is shown
nplot <- function(x=0,col="white",axes=F,main="",xlab="",ylab="",bty="n",...){
	plot(x,col=col,axes=axes,main=main,xlab=xlab,ylab=ylab,bty=bty,...)
}

rt<-function(f,col.names=F,row.names=F){
	sk=0
	if(col.names){
		sk=1
	}
	bf1 <- fread(f,skip=sk)
	if(row.names){
		m=as.matrix(bf1[,-1])
		rownames(m)=unlist(bf1[,1])
	}else{
		m=as.matrix(bf1)
	}
	if(col.names){
		colnames(m)=scan(f,nlines=1,what="char")
	}
	return(m)
}

rpm <-function(m,te=0,mult=0){
    if(te==0){
        te=dim(m)[1]
    }

    if(mult==0){
        mult=1e6
    }
    mult*(t(t(m[1:te,])/colSums(m[1:te,])))
}

tpm <- function(fpkm){
	exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

## get variances for matrix either on  rows=1 or columns=2
vars <- function(mat,j){apply(mat,j,var)}
## same for medians
medians <- function(mat,j){apply(mat,j,median)}
means <- function(mat,j){apply(mat,j,mean)}

## check which values are 0, same as is.na
is0=function(m){res=which(as.numeric(m)==0);if(length(res) ==0){0}else{res}}


## check the number of zeros in a vector
is0l <- function(m){length(which(as.numeric(m)==0))}

nis0l <- function(m){length(m)-length(which(as.numeric(m)==0))}

isxl<-function(m,minr=4){length(which(as.numeric(m)<minr))}


## intersects 2 vectors
inters <- function(x,y){intersect(x,y)}

## standardize
z <- function(x){ (x-mean(x))/sd(x)}

skewness <- function(x){  (1/length(x))*sum(z(x)^3)}

pskewness <- function(x){ (mean(x)-median(x)/sd(x) )}

kurtosis <- function(x) { (1/length(x))*sum(z(x)^4)}


## total SS
tss <- function(y){my=mean(y); sum((y-my)^2)}
## regression SS
ssr <- function(y,lmm){my=mean(y); sum((fitted(lmms)-my)^2)}
## residual SS
rss <- function(y,lmm){ sum((y-fitted(lmms))^2)}
## r squarred - tells how much of the variance is explained by the model
rs <- function(rss,tss){ 1-(rss/tss)}



getB <- function(countsM,mol=siM,thres=0,tl=1,xshift=0){
	totake=which(countsM > 0)
	take2=which(siM > thres)
	ftake=intersect(totake,take2)
	if(length(totake) == 0){ return(1) }
	if(tl==1){
		lmm=lm(log2(countsM[ftake])~log2(mol[ftake]))$coeff[[1]]

	}else{
		lmc=lm(countsM[ftake]~mol[ftake])
		lmm=lmc$coeff[[1]]
		if(xshift != 0){
			lmm=lmm+xshift*lmc$coeff[[2]]
		}
		return(lmm)
	}
	return(2^lmm)
}

getBs1 <- function(countsM,mol=siM,thres=0,tl=1,xshift=0){
	totake=which(countsM > 0)
	take2=which(siM > thres)
	ftake=intersect(totake,take2)

	if(tl == 1){
		lmm=lm(log2(countsM[ftake])~offset(log2(mol[ftake])))$coeff[[1]]
	}else{
		lmm=lm(countsM[ftake]~offset(mol[ftake]))$coeff[[1]]
		if(xshift != 0){
			lmm=lmm+xshift
		}
		return(lmm)
	}
	return(2^lmm)
}

## bin data values
binme <- function(u,v,br=1000){
	ttake <-which(is.finite(u)&is.finite(v))
	u=u[ttake]
	v=v[ttake]


	cc <- cut(u,br,labels=F)
	newv <- c()
	newu <- c()
	sdu <- c()
	sdv <- c()

	for(i in 1:br){
		index <- which(cc == i)
		if(length(index) > 0){
			newu <- c(newu,mean(u[index]))
			newv <- c(newv,mean(v[index]))
			if(length(index) > 1){
				sdu <- c(sdu,sd(u[index]))
				sdv <- c(sdv,sd(v[index]))
			}else{
				sdu <- c(sdu,1e-9)
				sdv <- c(sdv,1e-9)
			}
		}
	}
	return(cbind(newu,newv,sdu,sdv))
}


## this function is like cut but it puts equal numbers of data points into each bin
cut2 <- function(x,breaks=4){
	no=order(x)
	ne=floor(length(x)/breaks)
	rest=length(x)%%breaks
	out=rep(1:breaks,each=ne)
	if(rest > 0){
		out=c(out,rep(breaks+1,rest))
	}
	return(out[order(no)])
}



binme2 <- function(u,v,br=1000){
	cc <- cut2(u,br)
	newv <- c()
	newu <- c()
	sdu <- c()
	sdv <- c()
	for(i in 1:br){
		index <- which(cc == i)
		if(length(index) > 0){
			newu <- c(newu,mean(u[index]))
			newv <- c(newv,mean(v[index]))
			if(length(index) > 1){
				sdu <- c(sdu,sd(u[index]))
				sdv <- c(sdv,sd(v[index]))
			}else{
				sdu <- c(sdu,1e-9)
				sdv <- c(sdv,1e-9)
			}
		}
	}
	return(cbind(newu,newv,sdu,sdv))
}


## make correlation values
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{


	usr <- par("usr"); on.exit(par(usr))
	par(usr = c(0, 1, 0, 1))
	r <- abs(cor(x, y,method="spearman"))
	txt <- format(c(r, 0.123456789), digits = digits)[1]
	txt <- paste0(prefix, txt)
	if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
	text(0.5, 0.5, txt, cex = cex.cor * r)
}

panel.hist <- function(x,to=1,dens=0,...)
{
	usr <- par("usr"); on.exit(par(usr))
	par(usr = c(usr[1:2], 0, 1.5) )
	h <- hist(x, plot = FALSE)

	breaks <- h$breaks; nB <- length(breaks)
	y <- h$counts; y <- y/max(y)
	rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
	if(dens == 1){
		lines(density(x),col="red")
	}
}

mtick <- function(x,y,w=0.02,color="red"){
	for(i in 1:length(x)){
		lines(c(x[i],x[i]+w),c(y[i],y[i]),col=color)
	}
}

panel.bar <- function(x,...){
	usr <- par("usr"); on.exit(par(usr))
	z=0
	if(max(x) >=1){
		z=x
		x <- x/max(x)
	}


	par(usr = c(0,1, 0,1.5))
	## sort values
	px=sort(x)
	breaks <- seq(1,length(x)+1,1)/(length(x)+1); nB <- length(breaks)
	## normalise to height 1
	y <- px;
	rect(breaks[-nB], 0, breaks[-1], y, col = rgb(1-y+y/2,1-y,y,1), border=rgb(1-y+y/2,1-y,y,1),...)

	if(length(z) > 1){
		af=max(z)/4
		options("scipen"=-100, "digits"=1)
		text(0,c(0.25,0.5,0.75,1),labels=c(signif(af,2),signif(2*af,2),signif(3*af,2),signif(4*af,2)),pos=4)
		options("scipen"=0, "digits"=7)
	}else{
		options("scipen"=0, "digits"=7)
		text(0,c(0.25,0.5,0.75,1),labels=c(0.25,0.5,0.75,1),pos=4)
	}
	mtick(rep(0,4),c(0.25,0.5,0.75,1))

}

MsmoothScatter <- function(x, y = NULL, nbin = 128, bandwidth, colramp = colorRampPalette(c("white", blues9)),nrpoints = 100, pch = ".", cex = 1, col = "black", transformation = function(x) x^0.25, postPlotHook = box, xlab = NULL, ylab = NULL, xlim, ylim, xaxs = par("xaxs"), yaxs = par("yaxs"),to=1,dline=1, ...) {


	## take out points that are 0 in both variables
	if(to == 1){
		to=which(x == 0 & y==0)
		x=x[-to]
		y=y[-to]
	}

	if (!is.numeric(nrpoints) | (nrpoints < 0) | (length(nrpoints) != 
	1)) 
	stop("'nrpoints' should be numeric scalar with value >= 0.")
	xlabel <- if (!missing(x)) 
	deparse(substitute(x))
	ylabel <- if (!missing(y)) 
	deparse(substitute(y))
	xy <- xy.coords(x, y, xlabel, ylabel)
	xlab <- if (is.null(xlab)) 
	xy$xlab
	else xlab
	ylab <- if (is.null(ylab)) 
	xy$ylab
	else ylab
	x <- cbind(xy$x, xy$y)[is.finite(xy$x) & is.finite(xy$y),	, drop = FALSE]
	if (!missing(xlim)) {
		stopifnot(is.numeric(xlim), length(xlim) == 2, is.finite(xlim))
		x <- x[min(xlim) <= x[, 1] & x[, 1] <= max(xlim), ]
	}
	else {
		xlim <- range(x[, 1])
	}
	if (!missing(ylim)) {
		stopifnot(is.numeric(ylim), length(ylim) == 2, is.finite(ylim))
		x <- x[min(ylim) <= x[, 2] & x[, 2] <= max(ylim), ]
	}
	else {
		ylim <- range(x[, 2])
	}
	map <- grDevices:::.smoothScatterCalcDensity(x, nbin, bandwidth)
	xm <- map$x1
	ym <- map$x2
	dens <- map$fhat
	dens[] <- transformation(dens)
	image(xm, ym, z = dens, col = colramp(256), xlab = xlab,	ylab = ylab, xlim = xlim, ylim = ylim, xaxs = xaxs, yaxs = yaxs,	...)
	if (!is.null(postPlotHook)) 
	postPlotHook()
	if (nrpoints > 0) {
		nrpoints <- min(nrow(x), ceiling(nrpoints))
		stopifnot((nx <- length(xm)) == nrow(dens), (ny <- length(ym)) == ncol(dens))
		ixm <- 1L + as.integer((nx - 1) * (x[, 1] - xm[1])/(xm[nx] - xm[1]))
		iym <- 1L + as.integer((ny - 1) * (x[, 2] - ym[1])/(ym[ny] - ym[1]))
		sel <- order(dens[cbind(ixm, iym)])[seq_len(nrpoints)]
		points(x[sel, ], pch = pch, cex = cex, col = col)
	}
	if(dline==1){
		abline(0,1,col="red")
	}
}




pairs.default <- function (x, labels, panel = points, ..., horInd = 1:nc, verInd = 1:nc,	lower.panel = panel, upper.panel = panel, diag.panel = NULL,	text.panel = textPanel, label.pos = 0.5 + has.diag/3, line.main = 3,	cex.labels = NULL, font.labels = 1, row1attop = TRUE, gap = 1,	log = "") 
{
	if (doText <- missing(text.panel) || is.function(text.panel)) 
	textPanel <- function(x = 0.5, y = 0.5, txt, cex, font) text(x, y, txt, cex = cex, font = font)
	localAxis <- function(side, x, y, xpd, bg, col = NULL, main, oma, ...) {
		xpd <- NA
		if (side%%2L == 1L && xl[j]) 
		xpd <- FALSE
		if (side%%2L == 0L && yl[i]) 
		xpd <- FALSE
		if (side%%2L == 1L) 
		Axis(x, side = side, xpd = xpd, ...)
		else Axis(y, side = side, xpd = xpd, ...)
	}
	localPlot <- function(..., main, oma, font.main, cex.main) plot(...)
	localLowerPanel <- function(..., main, oma, font.main, cex.main) lower.panel(...)
	localUpperPanel <- function(..., main, oma, font.main, cex.main) upper.panel(...)
	localDiagPanel <- function(..., main, oma, font.main, cex.main) diag.panel(...)
	dots <- list(...)
	nmdots <- names(dots)
	if (!is.matrix(x)) {
		x <- as.data.frame(x)
		for (i in seq_along(names(x))) {
			if (is.factor(x[[i]]) || is.logical(x[[i]])) 
			x[[i]] <- as.numeric(x[[i]])
			if (!is.numeric(unclass(x[[i]]))) 
			stop("non-numeric argument to 'pairs'")
		}
	}
	else if (!is.numeric(x)) 
	stop("non-numeric argument to 'pairs'")
	panel <- match.fun(panel)
	if ((has.lower <- !is.null(lower.panel)) && !missing(lower.panel)) 
	lower.panel <- match.fun(lower.panel)
	if ((has.upper <- !is.null(upper.panel)) && !missing(upper.panel)) 
	upper.panel <- match.fun(upper.panel)
	if ((has.diag <- !is.null(diag.panel)) && !missing(diag.panel)) 
	diag.panel <- match.fun(diag.panel)
	if (row1attop) {
		tmp <- lower.panel
		lower.panel <- upper.panel
		upper.panel <- tmp
		tmp <- has.lower
		has.lower <- has.upper
		has.upper <- tmp
	}
	nc <- ncol(x)
	if (nc < 2L) 
	stop("only one column in the argument to 'pairs'")
	if (!all(horInd >= 1L && horInd <= nc)) 
	stop("invalid argument 'horInd'")
	if (!all(verInd >= 1L && verInd <= nc)) 
	stop("invalid argument 'verInd'")
	if (doText) {
		if (missing(labels)) {
			labels <- colnames(x)
			if (is.null(labels)) 
			labels <- paste("var", 1L:nc)
		}
		else if (is.null(labels)) 
		doText <- FALSE
	}
	oma <- if ("oma" %in% nmdots) 
	dots$oma
	main <- if ("main" %in% nmdots) 
	dots$main
	if (is.null(oma)) 
	oma <- c(4, 4, if (!is.null(main)) 6 else 4, 4)
	opar <- par(mfrow = c(length(horInd), length(verInd)), mar = rep.int(gap/2,			4), oma = oma)
	on.exit(par(opar))
	dev.hold()
	on.exit(dev.flush(), add = TRUE)
	xl <- yl <- logical(nc)
	if (is.numeric(log)) 
	xl[log] <- yl[log] <- TRUE
	else {
		xl[] <- grepl("x", log)
		yl[] <- grepl("y", log)
	}
	for (i in if (row1attop) 
	verInd
	else rev(verInd)) for (j in horInd) {
		l <- paste0(ifelse(xl[j], "x", ""), ifelse(yl[i], "y", ""))
		localPlot(x[, j], x[, i], xlab = "", ylab = "", axes = FALSE, type = "n", ..., log = l)
		if (i == j || (i < j && has.lower) || (i > j && has.upper)) {
			box()
			if (i == 1 && (!(j%%2L) || !has.upper || !has.lower)) 
			localAxis(1L + 2L * row1attop, x[, j], x[, i],						...)
			if (i == nc && (j%%2L || !has.upper || !has.lower)) 
			localAxis(3L - 2L * row1attop, x[, j], x[, i], ...)
			if (j == 1 && (!(i%%2L) || !has.upper || !has.lower)) 
			localAxis(2L, x[, j], x[, i], ...)
			if (j == nc && (i%%2L || !has.upper || !has.lower)) 
			localAxis(4L, x[, j], x[, i], ...)
			mfg <- par("mfg")
			if (i == j) {
				if (has.diag) 
				localDiagPanel(as.vector(x[, i]), ...)
				if (doText) {
					par(usr = c(0, 1, 0, 1))
					if (is.null(cex.labels)) {
						l.wid <- strwidth(labels, "user")
						cex.labels <- max(0.8, min(2, 0.9/max(l.wid)))
					}
					xlp <- if (xl[i]) 
					10^0.5
					else 0.5
					ylp <- if (yl[j]) 
					10^label.pos
					else label.pos
					text.panel(xlp, ylp, labels[i], cex = cex.labels,									font = font.labels)
				}
			}
			else if (i < j) 
			localLowerPanel(as.vector(x[, j]), as.vector(x[, i]), ...)
			else localUpperPanel(as.vector(x[, j]), as.vector(x[, i]), ...)
			if (any(par("mfg") != mfg)) 
			stop("the 'panel' function made a new plot")
		}
		else par(new = FALSE)
	}
	if (!is.null(main)) {
		font.main <- if ("font.main" %in% nmdots) 
		dots$font.main
		else par("font.main")
		cex.main <- if ("cex.main" %in% nmdots) 
		dots$cex.main
		else par("cex.main")
		mtext(main, 3, line.main, outer = TRUE, at = 0.5, cex = cex.main, font = font.main)
	}
	invisible(NULL)
}

dplot <- function(x,y,mt="plot",am=0,bm=0){
	xl1 <- min(x)
	xl2 <- max(x)

	yl1 <- min(y)
	yl2 <- max(y)

	amin=min(xl1,yl1)
	amax=max(xl2,yl2)

	if(am){
		amin=am
		amax=bm
	}


	move=0.7

	par(fig=c(0,move,0,move))
	plot(x,y,pch=16,col="red",xlim=c(amin,amax),ylim=c(amin,amax),main="")

	par(fig=c(0,move,move-0.25,1),new=T)
	plot(density(x),col="red",xlim=c(amin,amax),xlab="",ylab="",xaxt="none",main=mt)

	par(fig=c(move-0.16,1,0,move),new=T)
	d.b=density(y)
	plot(d.b$y,d.b$x,"l",col="blue",ylim=c(amin,amax),ylab="",xlab="",yaxt="none")

	par(fig=c(move-0.09,1,move-0.2,1),new=T)

	plot(density(x),"l",col="red",xlim=c(amin,amax),ylim=c(0,max(density(x)$y,density(y)$y)),xlab="",ylab="",main="")
	lines(density(y),col="blue")
}


scatterhist <- function(x, y, px=c(),py=c(),xl="", yl="",xmin=-10,xmax=10,maint="plot"){
	zones=matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
	layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))
	xhist = hist(x, plot=FALSE,breaks=50)
	yhist = hist(y, plot=FALSE,breaks=50)
	top = max(c(xhist$density, yhist$density))
	par(mar=c(5,4,1,1),cex=2)
	plot(x,y,xlim=c(xmin,xmax),ylim=c(xmin,xmax),pch=20,xlab=xl,ylab=yl)
	if(length(px) > 0 & length(py) > 0){
		points(px,py,col="blue",pch=20)
		lx=length(x)
		ffrac=length(px)/lx
		legend("topleft",c(paste("within 2sd:",lx-length(px)),paste("outside 2sd:",length(px)),paste("exp. outside 2sd\n 5% of",lx,"=",0.05*lx) ),pch=20,col=c("black","blue","white"))

	}


	par(mar=c(0,3,3,1))
	xhist = hist(x, breaks=50,freq=F,xlim=c(xmin,xmax),xaxt="n",yaxt="n",main=maint,col="blue")
	lines(density(x,na.rm=T),col="red",lwd=2)
	par(mar=c(3,0,1,1))

	plot(0,col="white",xlab="",ylab="",xlim=c(0,max(yhist$density)+0.1),ylim=c(xmin,xmax),xaxt="n",yaxt="n",main="",bty="n")		
	dx=density(y,na.rm=T)
	size=abs(abs(yhist$mids[2])-abs(yhist$mids[1]))	
	rect(rep(0,length(yhist$mids)),yhist$mids,yhist$density,yhist$mids+size,col="blue")
	lines(dx$y,dx$x,lwd=2,col="red")


	#par(oma=c(3,3,0,0))
	#mtext(xlab, side=1, line=1, outer=TRUE, adj=0,at=.8 * (mean(x) - min(x))/(max(x)-min(x)))
	#mtext(ylab, side=2, line=1, outer=TRUE, adj=0,at=(.8 * (mean(y) - min(y))/(max(y) - min(y))))
}

## get indices of matrix given the unlisted index ...
geti <- function(index,mat){
	col=ceil(index/dim(mat)[2]); 
	row=index-( (col-1)*dim(mat)[1]);
	ret=cbind(row,col);
	rownames(ret)=index;
	return(ret)
}

exact.p <- function(r,n){
	r= -abs(r)
	tstat <- r*sqrt((n-2)/(1-r^2))
	return (pt(tstat,n-2))
}


## return x and y values of the maximum of a density function object
get_max_d <- function(v){mi=which(v$y==max(v$y));c(v$x[mi],v$y[mi])}


## functions for upper triangle matrix indices if ordered by columns and 
## triangle numbers are the last elem of each column
## get column in upper.tri matrix given index
col.i <- function(x){ceiling(-0.5+0.5*sqrt(1+8*x))}
## alle decimalzahlen sind groesser als die Dreieckszahl und daher in der naechsten Spalte, als ceiling

## Fuer die Zeile nehmen wir wieder den index 
## row is, n is value of col.i, x is the index as before, n ist ja immer N-1, wenn N unsere Quadratmatrix ist, da die diagonale nicht zaehlt, 
## daher ist elems dann die Dreieckszahl in der Spalte vor der Spalte von x, daher ist die differenz gleich der Zeile von x
row.i <-function(x,n){
    elems=n*(n+1)/2-n
    x-elems
}


## beide funktionen zusammengefasst!
## so getting an index in upper triangle logic we return the row and column of the square matrix
uti <- function(x,diag=F){
    if(diag){
        matrix(c(row.i(x,col.i(x)),col.i(x)),ncol=2)
    }else{
        matrix(c(row.i(x,col.i(x)),col.i(x)+1),ncol=2)
    }
}

## having row and column of square matrix we get the index in upper triangle logic
## we just use the sum formula from gauss 
## we take one column less, get the triangle number and add the rows
utir <- function(x,y){ 0.5*(y-1)*(y)-(y-1)+x}
utirm=function(m){0.5*(m[,2]-1)*(m[,2])-(m[,2]-1)+m[,1]}



## get upper triangle matrix
ut <- function(x){x[upper.tri(x)]}
## get number of elements
utl <- function(x){n=dim(x)[1];n*(n+1)/2-n}


## calculate big cor matrix, cor calculated on columnwise comparisons
## taken from https://www.r-bloggers.com/bigcor-large-correlation-matrices-in-r/
bigcor <- function(MAT, nblocks = 10, verbose = TRUE, ...){
    require(ff, quietly = TRUE)
    require(Hmisc)
    NCOL <- ncol(MAT)
     
    ## test if ncol(x) %% nblocks gives remainder 0
    if (NCOL %% nblocks != 0) stop("Choose different 'nblocks' so that ncol(x) %% nblocks = 0!")
     
    ## preallocate square matrix of dimension
    ## ncol(x) in 'ff' single format
    corMAT <- ff(vmode = "single", dim = c(NCOL, NCOL))
    pMAT <- ff(vmode = "single", dim = c(NCOL, NCOL))
     
    ## split column numbers into 'nblocks' groups
    SPLIT <- split(1:NCOL, rep(1:nblocks, each = NCOL/nblocks))
     
    ## create all unique combinations of blocks
    COMBS <- expand.grid(1:length(SPLIT), 1:length(SPLIT))
    COMBS <- t(apply(COMBS, 1, sort))
    COMBS <- unique(COMBS)
     
    ## iterate through each block combination, calculate correlation matrix
    ## between blocks and store them in the preallocated matrix on both
    ## symmetric sides of the diagonal
    for (i in 1:nrow(COMBS)) {
        print(i)
        COMB <- COMBS[i, ]
        G1 <- SPLIT[[COMB[1]]]
        G2 <- SPLIT[[COMB[2]]]
        if (verbose) cat("Block", COMB[1], "with Block", COMB[2], "\n")
        flush.console()

        ## It's a bit inefficienct here since we dont need to calculate two full matrices but only one triangle of them 
        COR <- rcorr(MAT[, G1], MAT[, G2],...)
        xl=length(G1)
        yl=length(G2)
        tx=1:xl
        ty=(xl+1):(xl+yl)
        

        corMAT[G1, G2] <- COR[[1]][tx,ty]
        corMAT[G2, G1] <- t(COR[[1]][tx,ty])
        pMAT[G1, G2] <-   COR[[3]][tx,ty]
        pMAT[G2, G1] <-   t(COR[[3]][tx,ty])

        COR <- NULL
    }
    gc()
    return(list(corMAT,pMAT))
}

an <- function(x){as.numeric(x)}


## this function checks if two matrices have the same dimentions and same column order based on colname comparison
check_col_names <- function(m1,m2){
    if(dim(m1)[2] != dim(m2)[2]){
        print("matrices have different dimensions")
        return(0)
    }
    z2=vapply(colnames(m1),findx,c(1),m=colnames(m2))
    z1=vapply(colnames(m2),findx,c(1),m=colnames(m1))
    s12=sum(z1==z2)
    if(s12==dim(m1)[2]){
        return(1)
    }else{
        print("matrices have different colnames or different order of columns")
        return(0)
    }
}

## helper function to return indices in matrix m for given 
## index as x.
## we can use it with apply or vapply
## for each position in vector v we get the row index in matrix m
## use it as unlist(apply(as.matrix(names(if1)),1,findx,m=hl[,1]))
## or vapply(names(if1),findx,c(1),m=hl[,1])
findx<-function(m=matrix(c(1:3,6,4)),x=1){
    o=which(m %in%x)
    if(length(o) > 0){
        return(o)
    }
    return(0)
}


## transform varible name to a string
## name myVar is returned as "myVar"
get_name <- function(x){
	deparse(substitute(x))
}

## return variable content if it exists
## x is a string here
## if string is "myVar" and a variable names myVar exists
## then myVar is return -> so save it to some other variable then 
get_var <- function(x)
if (!is.null(r <- get0(x))){
	r
}


object.sizes <- function(){
	return(rev(sort(sapply(ls(envir=.GlobalEnv), function (object.name)
	object.size(get(object.name))))))
}

## ######################
## ######################
## ######################
## ######################
## pca plotting functions

## get total variance
mytvar <-function(x){
	return(signif(100*x$sdev^2/sum(x$sdev^2),3))
}

## nice scree plot
myscree <- function(x,pca=F,cs=F,varline=0){
	varex=x
	if(pca){
		varex=mytvar(x)
	}
	my=ceiling(max(varex)/10)*10
	if(cs){
		par(mar=c(5,4,4,4))
	}
	bp<-barplot(varex,ylab="% Variance",ylim=c(0,my),names.arg=1:length(varex),xlab="PC",col="cornflowerblue",las=1,border=NA)
	## if true add line of cumulative variance
	if(cs){
		sf<- 100/my
		sf2=(bp[2]-bp[1])/2
		lines(bp-sf2,cumsum(varex)/sf,type="s",lwd=1,col="brown")
		as <- seq(0,my,length.out=11)
		axis(4,as,seq(0,100,10),las=1)
		mtext("% cumulative Variance",4,line=2,col="brown")
		if(varline >0){
			abline(h=my*varline/100,col="grey")
			abline(v=which(cumsum(varex)>varline)[1]+(bp[2]-bp[1])/2,col="grey")

		}
	}
}

## nicer pca plot
mypcaplot <- function(tmat,pcplot=c(1,2),cols="cornflowerblue",use.text=0,varex="na",...){
	if(length(varex) ==1){
		varex=rep("na",dim(tmat)[2])
	}
	plot(tmat[,pcplot],pch=16,xlab=paste("PC",pcplot[1]," Var =",varex[pcplot[1]],"%"),ylab=paste("PC",pcplot[2]," Var =",varex[pcplot[2]],"%"),col=cols,...)
	if(use.text > 0){
		mytext=rownames(tmat)
		if(use.text==2){
			mytext=unlist(strsplit(rownames(tmat)," "))[seq(1,2*dim(tmat)[1],2)]
			text(tmat[,pcplot[1]],tmat[,pcplot[2]],mytext,pos=4,col=cols,cex=1.1)
		}else if(use.text==3){
			mytext=unlist(strsplit(rownames(tmat)," "))[seq(1,2*dim(tmat)[1],2)]
			mytext=unlist(strsplit(mytext,"dbd"))[seq(2,2*dim(tmat)[1],2)]
			text(tmat[,pcplot[1]],tmat[,pcplot[2]],mytext,col=cols,cex=1.1,pos=4)
		}else{
			text(tmat[,pcplot[1]],tmat[,pcplot[2]],mytext,pos=4,col=cols,cex=1.1)
		}
	}
}

## add braces
addbraces <- function(x){
	return(paste("(",x,")",sep=""))
}

rembraces <- function(x,start=2,end=0){
	lx=str_length(x)
	if(end == 0){
		end=lx-1
	}
	return(substr(x,start,end))
}


### BIC and AIC for kmeans object
BIC2 <- function(fit){
	m = ncol(fit$centers)
    n = length(fit$cluster)
	k = nrow(fit$centers)
    D = fit$tot.withinss
	return(data.frame(AIC = D + 2*m*k,BIC = D + log(n)*m*k,AICc = D + 2*m*k + 2*k*(k+1)/(n-k-1)))
}

