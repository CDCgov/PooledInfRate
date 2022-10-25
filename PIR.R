
# CDC Color
 cdc.blue <- "#005DAA"

"polygonBand" <-  function(x,lower,upper,f1=NULL,f2=NULL,fill=TRUE,col="lightgray",density,border=FALSE,...){
     if(!is.null(f1)) lower <- f1(x)
     if(!is.null(f2)) upper <- f2(x)
     if(fill) density <- -1
     if(!fill & missing(density)){
         density <- 10
         border <- TRUE
     }
     polygon(c(x,rev(x)), c(lower, rev(upper)), density=density, border=border,col=col,...)
}

"ellipse" <- function(a=1,b=1,h=0,k=0,from=0,to=2*pi,npts=1000){
  theta <- seq(from,to,length=npts)
  x <- cos(theta)
  y <- sin(theta)
  x <- x/a-h
  y <- y/b-k
  list(x=x,y=y)
}

"lemniscate" <- function(a=1,b=1,from=pi/2,to=3*pi/2,orientation=c("vertical","horizontal"),npts=1000){
   orientation <- match.arg(orientation)
   theta <- seq(from,to,length=npts)
   switch(orientation,
        vertical =  ans <- list(y=b*a*cos(theta)/(1+sin(theta)^2),x=a*sin(theta)*cos(theta)/(1+sin(theta)^2)),
        horizontal =  ans <- list(x=a*cos(theta)/(1+sin(theta)^2),y=b*a*sin(theta)*cos(theta)/(1+sin(theta)^2))
   )
   ans
}

# figure built with these xlim=c(-5,5), ylim=c(-20,20)
"drawMosquito" <- function(scale=1,bg=c("transparent","white"),lwd=3,col="black"){
    scale <- 1/scale # so the interpretation at the call is bigger scale = bigger image
    xlim <- scale*c(-5,5)
    ylim <- scale*c(-20,20)
    #BG <- match.arg(bg)
    par(bg=bg[1], pty="s")

    plot(lemniscate(1,10,pi/2,3*pi/2,orient="v"),lwd=lwd,xlim=xlim,ylim=ylim,axes=FALSE,xlab="",ylab="",type="l",col=col)
        lmn <- lemniscate(1,10,pi/2,3*pi/2,orient="v")
        polygonBand(lmn$x,lmn$x,lmn$y,col=col,border=FALSE)

      lines(lemniscate(1,4,from=0,to=pi/2,orient="v"),lwd=lwd,col=col)
        lmn <- lemniscate(1,4,from=0,to=pi/2,orient="v")
        polygonBand(lmn$x,lmn$x,lmn$y,col=col,border=FALSE)

      lines(lemniscate(1,4,from=3*pi/2,to=2*pi,orient="v"),lwd=lwd,col=col)
        lmn <- lemniscate(1,4,from=3*pi/2,to=2*pi,orient="v")
        polygonBand(lmn$x,lmn$x,lmn$y,col=col,border=FALSE)

      lines(ellipse(a=5,b=1.5,h=0,k=-4.7),lwd=lwd,col=col)
        ell <- ellipse(a=5,b=1.5,h=0,k=-4.7)
        polygonBand(ell$x,ell$x,ell$y,col=col,border=FALSE)

      segments(0,0,1,-5,lwd=lwd,col=col)
        segments(1,-5,0.5,-11,lwd=lwd,col=col)
        segments(0.5,-11,1.2,-18,lwd=lwd,col=col)
      segments(0,0,-1,-5,lwd=lwd,col=col)
        segments(-1,-5,-0.5,-11,lwd=lwd,col=col)
        segments(-0.5,-11,-1.2,-18,lwd=lwd,col=col)
      segments(0,0,1.4,-4,lwd=lwd,col=col)
        segments(1.4,-4,1,-10,lwd=lwd,col=col)
        segments(1,-10,1.3,-15,lwd=lwd,col=col)
      segments(0,0,-1.4,-4,lwd=lwd,col=col)
        segments(-1.4,-4,-1,-10,lwd=lwd,col=col)
        segments(-1,-10,-1.3,-15,lwd=lwd,col=col)
      segments(0.35,2.8,1,4,lwd=lwd,col=col)
        segments(1,4,1.15,8,lwd=lwd,col=col)
        segments(1.15,8,2,11,lwd=lwd,col=col)
      segments(-0.35,2.8,-1,4,lwd=lwd,col=col)
        segments(-1,4,-1.15,8,lwd=lwd,col=col)
        segments(-1.15,8,-2,11,lwd=lwd,col=col)
      # proboscis & antennas
      #segments(0,5.4,0,9.5,lwd=lwd,col=col)
      #  segments(0,5.4,0.45,8,lwd=lwd,col=col)
      #  segments(0,5.4,-0.45,8,lwd=lwd,col=col)
      # proboscis & antennas
        # updated antennas with input from John-Paul and Janet
      segments(0,5.4,0,10.5,lwd=lwd,col=col)
        segments(0,5.4,0.15,5.6,lwd=lwd,col=col)
        segments(0,5.4,-0.15,5.6,lwd=lwd,col=col)
      # wings
      segments(0.2,0.9,3,2.4,lwd=lwd,col=col)
        segments(0,0,2,-0.5,lwd=lwd,col=col)
        segments(2,-0.5,3.1,1.65,lwd=lwd,col=col)
        ell.right <- ellipse(a=5,b=2.5,h=-3,k=-2,from=-pi/3,to=pi/2) # wingtip
        lines(ell.right,lwd=lwd,col=col)
     segments(-0.2,0.9,-3,2.4,lwd=lwd,col=col)
        segments(0,0,-2,-0.5,lwd=lwd,col=col)
        segments(-2,-0.5,-3.1,1.65,lwd=lwd,col=col)
        ell.left <- ell.right
        ell.left$x <- -ell.right$x
        lines(ell.left,lwd=lwd,col=col) # wingtip
    invisible(scale)
}


"addTickRot" <- function(a=2,b=5,x.shift=0,y.shift=0,scale=1,theta=0,bg=c("transparent","white"),lwd=3,col="black",add=TRUE){
    scale <- 1/scale # so the interpretation at the call is bigger scale = bigger image
    xlim <- scale*c(-5,5)
    ylim <- scale*c(-20,20)
    x.scale <- xlim[2]
    y.scale <- ylim[2]

    #BG <- match.arg(bg)
    par(bg=bg[1])

    rotMat <- matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2,byrow=TRUE)

    lmn <- lemniscate(a,b,from=pi/2,to=3*pi/2,orientation="horizontal")
    xy.rot <- rotMat %*% rbind(lmn$x/x.scale,lmn$y/y.scale)
    lmn$x <- xy.rot[1,]
    lmn$y <- xy.rot[2,]
    lmn$x <- x.scale*lmn$x + x.shift
    lmn$y <- y.scale*lmn$y + y.shift

    if(add){
        lines(lmn,col=col)
    } else {
        plot(lmn,lwd=lwd,xlim=xlim,ylim=ylim,axes=FALSE,xlab="",ylab="",type="l",col=col)
    }
    polygonBand(lmn$x,lmn$x,lmn$y,col=col,border=FALSE)

    points(0+x.shift,0+y.shift,pch=16,col=col,cex=1.5)
    #
    #segments(-0.25+x.shift,-1+y.shift,-0.2+x.shift,-2.5+y.shift,lwd=3,col=cdc.blue)
    #segments(-0.2+x.shift,-2.5+y.shift,0.5+x.shift,-3.5+y.shift,lwd=3,col=cdc.blue)
    #segments(-0.25+x.shift,1+y.shift,-0.2+x.shift,2.5+y.shift,lwd=3,col=cdc.blue)
    #segments(-0.2+x.shift,2.5+y.shift,0.5+x.shift,3.5+y.shift,lwd=3,col=cdc.blue)
   #
   # segments(-0.5+x.shift,-2+y.shift,-0.3+x.shift,-5.5+y.shift,lwd=3,col=cdc.blue)
   # segments(-0.5+x.shift,2+y.shift,-0.3+x.shift,5.5+y.shift,lwd=3,col=cdc.blue)
   #
   # segments(-0.8+x.shift,-3.1+y.shift,-1.5+x.shift,-6.5+y.shift,lwd=3,col=cdc.blue)
   # segments(-0.8+x.shift,3.1+y.shift,-1.5+x.shift,6.5+y.shift,lwd=3,col=cdc.blue)
   #
   # segments(-1.3+x.shift,-3.6+y.shift,-2.5+x.shift,-6.25+y.shift,lwd=3,col=cdc.blue)
   # segments(-1.3+x.shift,3.6+y.shift,-2.5+x.shift,6.25+y.shift,lwd=3,col=cdc.blue)

   seg.points <- matrix(c(
       -0.25, -1, -0.2, -2.5,
       -0.2,-2.5,0.5,-3.5,
       -0.25,1,-0.2,2.5,
       -0.2,2.5,0.5,3.5,
       #
       -0.5,-2,-0.3,-5.5,
       -0.5,2,-0.3,5.5,
       #
       -0.8,-3.1,-1.5,-6.5,
       -0.8,3.1,-1.5,6.5,
       #
       -1.3,-3.6,-2.5,-6.25,
       -1.3,3.6,-2.5,6.25),
       nrow=2)

   seg.points[1,] <- seg.points[1,]/x.scale
   seg.points[2,] <- seg.points[2,]/y.scale

   seg.points.rot <- rotMat %*% seg.points

   seg.points.rot <- matrix(as.vector(seg.points.rot),ncol=4,byrow=TRUE)

    for(i in 1:nrow(seg.points.rot))
        segments(x.scale*seg.points.rot[i,1]+x.shift, y.scale*seg.points.rot[i,2]+y.shift,
                 x.scale*seg.points.rot[i,3]+x.shift, y.scale*seg.points.rot[i,4]+y.shift, lwd=3,col=col)

    invisible(scale)
}

"addFleaRot" <- function(a=2,b=5,x.shift=0,y.shift=0,scale=1,theta=0,bg=c("transparent","white"),lwd=3,col="black",add=TRUE){
    scale <- 1/scale # so the interpretation at the call is bigger scale = bigger image
    xlim <- scale*c(-5,5)
    ylim <- scale*c(-20,20)
    x.scale <- xlim[2]
    y.scale <- ylim[2]
    #BG <- match.arg(bg)
    par(bg=bg[1])

    rotMat <- matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2,byrow=TRUE)


    # front side
    lmn <- lemniscate(a,b,from=-pi/2,to=pi/2,orientation="horizontal")
    xy.rot <- rotMat %*% rbind(lmn$x/x.scale,lmn$y/y.scale)
    lmn$x <- xy.rot[1,]
    lmn$y <- xy.rot[2,]
    lmn$x <- x.scale*lmn$x + x.shift
    lmn$y <- y.scale*lmn$y + y.shift

    if(add){
        lines(lmn,col=col)
    } else {
        plot(lmn,lwd=lwd,xlim=xlim,ylim=ylim,axes=FALSE,xlab="",ylab="",type="l",col=col)
    }
    polygonBand(lmn$x,lmn$x,lmn$y,col=col,border=FALSE)

    # back side
    lmn <- lemniscate(a,b,from=pi/2,to=3*pi/2,orientation="horizontal")
    xy.rot <- rotMat %*% rbind((lmn$x+1.2*a)/x.scale,lmn$y/y.scale)
    lmn$x <- xy.rot[1,]
    lmn$y <- xy.rot[2,]
    lmn$x <- x.scale*(lmn$x) + x.shift
    lmn$y <- y.scale*lmn$y + y.shift

    lines(lmn,col=col)
    polygonBand(lmn$x,lmn$x,lmn$y,col=col,border=FALSE)


    head.pts <- matrix(c( 0.2,0,
                   0.1,-0.1,
                   0.075,-0.3,
                  -0.05,-0.6,
                  -0.03,-0.5,
                  -0.02,-0.45,
                   0.25,0.1),nrow=2)
    head.pts[1,] <- head.pts[1,]/x.scale
    head.pts[2,] <- head.pts[2,]/y.scale

    head.pts.rot <- rotMat %*% head.pts

    head.pts.rot <- matrix(as.vector(head.pts.rot),ncol=2,byrow=TRUE)

    #points(0.2+x.shift,0+y.shift,pch=16,col=col,cex=1.25)
    #points(0.1+x.shift,-0.1+y.shift,pch=16,col=col,cex=1.5)
    #points(0.075+x.shift,-0.3+y.shift,pch=16,col=col,cex=1.5)
    #points(-0.05+x.shift,-0.6+y.shift,pch=16,col=col,cex=1.4)
    #points(-0.03+x.shift,-0.5+y.shift,pch=16,col=col,cex=1.5)
    #points(-0.02+x.shift,-0.45+y.shift,pch=16,col=col,cex=1.5)
    #points(0.25+x.shift,0.1+y.shift,pch=16,col=col,cex=1.6)

    #for(i in 1:nrow(head.pts.rot))
        points(head.pts.rot[1,1]+x.shift,head.pts.rot[1,2]+y.shift,pch=16,col=col,cex=1.25)
        points(head.pts.rot[2,1]+x.shift,head.pts.rot[2,2]+y.shift,pch=16,col=col,cex=1.5)
        points(head.pts.rot[3,1]+x.shift,head.pts.rot[3,2]+y.shift,pch=16,col=col,cex=1.5)
        points(head.pts.rot[4,1]+x.shift,head.pts.rot[4,2]+y.shift,pch=16,col=col,cex=1.4)
        points(head.pts.rot[5,1]+x.shift,head.pts.rot[5,2]+y.shift,pch=16,col=col,cex=1.5)
        points(head.pts.rot[6,1]+x.shift,head.pts.rot[6,2]+y.shift,pch=16,col=col,cex=1.5)
        points(head.pts.rot[7,1]+x.shift,head.pts.rot[7,2]+y.shift,pch=16,col=col,cex=1.6)


    # legs
    # front
    seg.points <- matrix( c(
        0.15,-0.7,0.25,-2.5,
        0.25,-2.5,0.1,-4.5,

        0.15,-0.7,0.4,-3,
        0.4,-3,0.3,-4.5,

        0.45,-1.6,0.6,-3.25,
        0.6,-3.25,0.7,-3.4,
        0.7,-3.4,0.6,-5.5,

        0.55,-1.8,0.8,-3.25,
        0.8,-3.25,0.9,-3.4,
        0.9,-3.4,0.8,-5.5,

        0.95,-2.7,1.2,-4,
        1.2,-4,1.4,-4,
        1.4,-4,2.4,-5.5,

        0.95,-2.7,1.1,-4.3,
        1.1,-4.3,1.2,-4.4,
        1.2,-4.4,2,-5.8
    ),nrow=2)

    seg.points[1,] <- seg.points[1,]/x.scale
    seg.points[2,] <- seg.points[2,]/y.scale

    seg.points.rot <- rotMat %*% seg.points

    seg.points.rot <- matrix(as.vector(seg.points.rot),ncol=4,byrow=TRUE)

    for(i in 1:nrow(seg.points.rot))
        segments(x.scale*seg.points.rot[i,1]+x.shift, y.scale*seg.points.rot[i,2]+y.shift,
                 x.scale*seg.points.rot[i,3]+x.shift, y.scale*seg.points.rot[i,4]+y.shift, lwd=3,col=col)


    invisible(scale)
}




png("MosqTickFleaGlyphV2.png",res=600,height=11,width=8.5,units="in")

drawMosquito(scale=0.75,col=cdc.blue,bg="transparent")
addTickRot(a=2,b=5,x.shift=-3,y.shift=17,col=cdc.blue,theta=-pi/6)
addFleaRot(a=2,b=4,x.shift=3,y.shift=17,col=cdc.blue,bg="white",theta=pi/6)

dev.off()

library(magick)
img <- image_read("MosqTickFleaGlyphV2.png")
img <- image_background(img, color = "none")
img <- image_trim(img)
# img <- image_resize(img, geometry_size_percent(width=100))

library(hexSticker)


# CDC Colors
#s_width=0.9, s_height=1.25
plot(sticker(img,package="PooledInfRate",
             s_width=1.25,s_height=1.25,
             h_color=cdc.blue,h_fill="white",
             p_size=25,p_y=1.5, p_color=cdc.blue,
             s_x=1,s_y=0.75,
             h_size=1.2, # 1.2 is default
             filename="PIRV2.png",
             asp = TRUE,
             dpi=600,
             white_around_sticker = FALSE) )



