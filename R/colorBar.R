# a helper function to plot a color bar.
colorBar<- function(xl, yb, xr, yt, v.min, v.max, col.l, alpha=1, color.spec="hsv",title=""){
    nticks <- 11; 
    ticks  <- seq(v.min, v.max, len=nticks);
    tmp.v  <- seq(v.min, v.max, length.out=100)
    lut    <- col.l;
    #lut   <- colorRampPalette(c("white", "red"))(100);
    xs    <- (xr-xl)/(length(lut)-1);
    xl.old <- xl;
    # change from 1:(length(lut)-1) to seq_len(length(lut)-1))
    for (i in seq_len(length(lut)-1)) {
        rect(xl.old, yb, xl.old+xs, yt, col=lut[i], border=NA)
        xl.old <- xl.old + xs;
    }
    v.med <- (v.max+v.min)/2;
    xm  <- (xl+xr)/2;
    v.max <- round(v.max,2);
    v.min <- round(v.min,2);
    v.med <- round(v.med,2);
    x.d   <- xr-xl;
    y.d   <- yt-yb;
    text(xl, yt+x.d, title, cex=1);
    text(xr, yb-y.d, v.max, cex=0.6);
    text(xl, yb-y.d, v.min, cex=0.6);
    text(xm, yb-y.d, v.med, cex=0.6);
    segments(xl, yb-y.d/5, xr,  yb-y.d/5, col="black", lwd=0.8);
    segments(xl, yb-y.d/5, xl,  yb-y.d/2, col="black", lwd=0.8);
    segments(xm, yb-y.d/5, xm,  yb-y.d/2, col="black", lwd=0.8);
    segments(xr, yb-y.d/5, xr,  yb-y.d/2, col="black", lwd=0.8);
}