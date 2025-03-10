# adapted from 'GroupedCircularBoxplot' from Josh Berlinski, available at https://github.com/jdberlinski/GroupedCircularBoxplot
# modified by Fan Dai

# functions to obtain the Cartesian coordinates (x,y) from circular plot objects
# adapted from 'plotrix::draw.arc', 'draw.radial.line' and 'draw.circle' from 'plotrix' package by Jim Lemon et al.
# "par("usr")" set by plotting the baseline circle in 'GroupedCircularBoxplot.3D'
draw.arc.xy <- function(x=1, y=NULL, radius=1, angle1=deg1*pi/180,
                        angle2=deg2*pi/180, deg1=0, deg2=45, n=0.05,lwd=NA) {
  
  if (all(is.na(lwd)))
    lwd <- par("lwd")
  
  xylim<-par("usr") 
  ymult <- plotrix::getYmult()
  devunits <- dev.size("px")
  
  xy <- xy.coords(x, y)
  x <- xy$x
  y <- xy$y
  a1 <- pmin(angle1, angle2)
  a2 <- pmax(angle1, angle2)
  angle1 <- a1
  angle2 <- a2
  
  delta.angle <- (angle2 - angle1)
  if (n != as.integer(n)){
    n <- as.integer(1+delta.angle/n) 
  }
  delta.angle <- delta.angle/n
  angleS <- angle1 + seq(0, length=n) * delta.angle
  angleE <- c(angleS[-1], angle2)
  if (n > 1){
    half.lwd.user <- (lwd/2)*(xylim[2]-xylim[1])/devunits[1]
    adj.angle = delta.angle*half.lwd.user/(2*(radius+half.lwd.user))
    angleS[2:n] = angleS[2:n] - adj.angle
    angleE[1:(n-1)] = angleE[1:(n-1)] + adj.angle
  }
  p1x <- x + radius * cos(angleS)
  p1y <- y + radius * sin(angleS) * ymult
  p2x <- x + radius * cos(angleE)
  p2y <- y + radius * sin(angleE) * ymult
  return(list(n = n, x1 = p1x, y1 = p1y, x2 = p2x, y2 = p2y))
}

draw.radial.line.xy <- function(start, end, center=c(0, 0), angle=0, deg=NA, expand=FALSE) 
{
  if (is.na(deg))
    deg <- angle*180/pi
  angle <- deg*pi/180
  cosang <- cos(angle)
  sinang <- sin(angle)
  xylim <- par("usr")
  plotdim <- par("pin")
  ymult <- (xylim[4] - xylim[3])/(xylim[2] - xylim[1]) * plotdim[1]/plotdim[2]
  if (!expand || end <= 0)
  {
    px <- c(start, end)
    return(list(x = center[1]+px*cosang, y = center[2]+px*sinang*ymult))
  }
  else
  {
    devunits <- dev.size("px")
    lwdend <- (lwd/2)*(xylim[2]-xylim[1])/devunits[1]
    px <- c(start, start, end, end)
    lwdstart <- lwdend * start/end
    py <- c(lwdstart, -lwdstart, -lwdend, lwdend)
    pxt <- center[1] + px*cosang - py*sinang*ymult
    pyt <- center[2] + px*sinang + py*cosang*ymult
    return(list(x = pxt, y = pyt))
  }
}

draw.circle.xy <- function(x,y,radius,nv=100,border=NULL,density=NULL,angle=45) {
  
  xylim<-par("usr")
  plotdim<-par("pin")
  ymult<-plotrix::getYmult()
  angle.inc<-2*pi/nv
  angles<-seq(0,2*pi-angle.inc,by=angle.inc)
  
  for(circle in 1:length(radius)) {
    xv<-cos(angles)*radius[circle]+x
    yv<-sin(angles)*radius[circle]*ymult+y
  }
  return(list(x=xv,y=yv))
}

#' Create grouped circular boxplots around the 3D torus surface.
#'
#' @description Given a named list of circular objects, create a grouped circular boxplot.
#'
#' @param data_in List of circular objects to be plotted
#' @param template One of "degrees" (default), "radians" and "geographics" with eight directions
#' @param marg character specifying the baseline 2D circle region as either "large" (default) or "small"
#' @param shrink Numeric specifying the factor by which to scale the baseline circular plot. Numbers less than 1 will increase the size.
#' @param H Logical indicating if each data point should be drawn outside the hinges of the boxplot
#' @param constant Numeric specifying the multiplicative factor determining how far whiskers extend from box. A value of
#' "optimal" will choose values based on a von Mises distribution (see Buttarazzi et al. 2018)
#' @param lwd Scalar indicating the width of the quartile boxplot lines if \code{minimal = TRUE}
#' @param plot_cols Vector with the same length as `data_in`, specifying the color of the boxplot
#' @param line_cols Vector with the same length as `data_in`, specifying the color of the median lines
#' @param arrow_cols Vector with the same length as `data_in`, specifying the color of the arrows pointing to the
#' medians. Default value is the same as line_cols
#' @param draw_arrow Logical specifying if arrows pointing to each median should be drawn
#' @param minimal Logical. If true, a quartile boxplot is created with simple colored line and point will replace the box and median line. Radial lines
#' at fence points will also not be drawn
#' @param scale_widths Logical, should the width of each boxplot be scaled based on (the square root of) it's distance from the center?
#' @examples
#' library(circular)
#' library(GroupedCircularBoxplot)
#' set.seed(123)
#' data <- list(
#'     x = rvonmises(100, circular(pi), 5),
#'     y = rvonmises(100, circular(pi/2), 2.5),
#'     z = rvonmises(100, circular(7*pi/4), 8)
#' )
#' GroupedCircularBoxplot.3D(data)
#' @export
#' @author Josh Berlinski
#' @author Fan Dai
#' @importFrom graphics legend text par points
GroupedCircularBoxplot.3D <- function(
    data_in,
    template = "degrees",
    marg = "large",
    shrink = 1,
    H = FALSE,
    constant = "optimal",
    lwd = 2,
    plot_cols = RColorBrewer::brewer.pal(8, "Set2"),
    line_cols = RColorBrewer::brewer.pal(8, "Dark2"),
    arrow_cols = line_cols,
    draw_arrow = TRUE,
    minimal = FALSE,
    scale_widths = FALSE
) {
  
  # if only a circular vector is passed as data in, make it a list and don't plot a legend
  if (circular::is.circular(data_in)) {
    data_in <- list(a = data_in)
  }
  
  # error checking
  if (!is.list(data_in))
    stop("`data_in` must be a list of cicular vectors.")
  
  n_seq <- length(data_in)
  which_circ <- do.call(c, lapply(data_in, circular::is.circular))
  if (!all(which_circ))
    stop("`data_in` must be a list of circular vectors.")
  if (length(plot_cols) < n_seq)
    stop("Number of supplied colors is less than the number of sequences. Provide more colors.")
  
  # output list
  summary_output <- vector(mode = "list", length = n_seq)
  names(summary_output) <- names(data_in)
  
  # build torus mesh
  torus.R <- 0.8; torus.r <- 0.2
  n.ring <- n_seq + 1 # at least four vertices for three groups
  torus.M <- plot3D::mesh(seq(0, 2*pi,length.out=n.ring), seq(0, 2*pi,length.out=n.ring))
  alpha <- torus.M$x; beta <- torus.M$y
  torus.x <- (torus.R + torus.r*cos(alpha)) * cos(beta)
  torus.y <- (torus.R + torus.r*cos(alpha)) * sin(beta)
  torus.z <-  torus.r * sin(alpha) 
  tmp_seq = z_seq <- NULL
  for (i in 1:(n.ring-1)) {
    tmp_seq <- c(tmp_seq, sqrt(torus.x[i,1]^2+torus.y[i,1]^2))
    z_seq <- c(z_seq, torus.z[i,1])
  }
  
  size_val <- 1 / sqrt(tmp_seq)
  
  for (curr_seq in 1:n_seq) {
    A <- data_in[[curr_seq]]
    
    ##Check if Median is uniquely defined
    if(is.na(circular::median.circular(A)))
      stop("The median is not unique for this data set. \n The circular boxplot is not drawn.")
    
    if(constant=="optimal"){
      conc     <- circular::A1inv(circular::rho.circular(A))
      q1       <- circular::qvonmises(0.25, mu=circular::circular(0), kappa = conc)
      me       <- circular::qvonmises(0.5, mu=circular::circular(0), kappa = conc)
      q3       <- circular::qvonmises(0.75, mu=circular::circular(0), kappa = conc)
      box      <- range(c(q1,me,q3))
      q9965    <- circular::qvonmises(1-(0.007/2), mu=circular::circular(0),kappa = conc)
      q0035    <- circular::qvonmises((0.007/2), mu=circular::circular(0),kappa = conc)
      constant <- range(c(q9965,q3))/box
    }
    
    # TODO: find some way to reset the graphics parameters without nuking layouts
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    
    if (marg == "small"){
      par(oma=c(0,0,0,0))}else if (marg == "large"){
      par(mai=c(0.0,0.0,0,0))}
    
    # convert A to radians
    set1 <- circular::conversion.circular(A, units="radians", modulo="2pi", zero=0, rotation="counter")
    
    ## median and IQR
    x <- set1
    AM <- circular::median.circular(x) + pi
    x <- x[which(x != AM)]
    x2 <- as.matrix(sort(circular::circular(x - AM, modulo="2pi")))
    
    Combined <- cbind(seq_along(x2), rev(seq_along(x2)))
    Tukeyway <- apply(Combined, 1, min)
    
    data <- cbind(circular::circular(x2 + AM, modulo="2pi"), Tukeyway)
    colnames(data) <- c("observations",  "depth")
    
    CTM <- which(data[, 2] >= which.max(data[, 2]))
    CTM <- circular::mean.circular(circular::circular(data[CTM, 1], modulo = "2pi"))
    CTM <- circular::circular(CTM, modulo="2pi") # in case CTM > 2pi
    # CTM <- circular::circular(mean(circular::circular(data[CTM, 1], modulo = "2pi")), modulo="2pi")
    
    n <- length(x)
    depthofmedian <- floor((1+n) / 2)
    # depthofmedian <- round(((1+n)/2) - 0.1)
    depthofquartiles <- (1+depthofmedian)/2
    
    if (depthofquartiles %% 1 == 0) {
      quartiles <- which(data[,2] == depthofquartiles)
      # quartiles <- which(data[,2] == round(1+depthofmedian)/2)
      qA <-  circular::circular(data[quartiles[1], 1], modulo="2pi")
      qC <-  circular::circular(data[quartiles[2], 1], modulo="2pi")
    } else {
      depthq1 <- depthofquartiles + 0.5
      depthq2 <- depthofquartiles - 0.5
      q1 <- which(data[,2] == depthq1)
      q2 <- which(data[,2] == depthq2)
      qA  <- mean(circular::circular(data[c(q1[1], q2[1]), 1], modulo="2pi"))
      qC  <- mean(circular::circular(data[c(q1[2], q2[2]), 1], modulo="2pi"))
    }
    
    IQRdepth <- which(data[,2] >= depthofquartiles)
    IQR <- c(data[IQRdepth,1],qA,qC) # qA and qC should already be in here? unless they are an average
    IQRange <- range(c(qA,qC))
    
    set_1 <- set1
    fi <- CTM
    # fi <- circular::as.circular(CTM) # CTM Is already circular...
    
    # draw the template
    if (curr_seq == 1) { # plot the baseline circle with radius = torus.R + torus.r to set the coordinate parameters
        circular::plot.circular(
        circular::circular(NA, modulo = "2pi"),
        cex=0.5, axes=FALSE, shrink=shrink,
        template=NULL, control.circle=circular::circle.control(col="gray85", lty=2, lwd=0.5),
        xlim = c(-tmp_seq[curr_seq], tmp_seq[curr_seq])
      )
      
      pxy = draw.circle.xy(0,0,tmp_seq[curr_seq])
      rgl::plot3d(pxy$x, pxy$y, z_seq[curr_seq], 
                  type = "p", col = "gray85",pch = 20, size = 0.8,add = T)
      
      # label the directions
      lab_coord <- cbind(
        cos(circular::rad(circular::circular(c(0,45,90,135,180,225,270,315)))),
        sin(circular::rad(circular::circular(c(0,45,90,135,180,225,270,315))))
        )
        if(template=="degrees") {
          ax_labels <- c("0","45","90","135","180","225","270","315")
          tmult <- 1.2
        } else if(template=="radians") {
          ax_labels <- c(expression(0,frac(pi,4),frac(pi,2),frac(3*pi,4),
                                    pi,frac(5*pi,4),frac(3*pi,2),frac(7*pi,4)))
          tmult <- 1.25
        } else if(template=="geographics") {
          ax_labels <- c("E","NE","N","NW","W","SW","S","SE")
          tmult <- 1.2
        }
        rgl::text3d(
          tmult * circular::circular(lab_coord[,1], units = "radians"), 
          tmult * circular::circular(lab_coord[,2], units = "radians"),
          z_seq[curr_seq], texts = ax_labels,
        )
      }else {
        pxy = draw.circle.xy(0,0,tmp_seq[curr_seq])
        rgl::plot3d(pxy$x, pxy$y, z_seq[curr_seq], 
                    type = "p", col = "gray85",pch = 20, size = 0.8, add=TRUE)
      
    }
  
    
    # display data points
    if (H){
      set_1_circ = circular::circular(set_1)
      rgl::points3d(cos(set_1_circ)*tmp_seq[curr_seq], 
                    sin(set_1_circ)*tmp_seq[curr_seq], 
                    z_seq[curr_seq], 
                    cex=2, pch=8, add = T)
    }
    # is this bad?
    round_circ <- function(x) circular::rad(round(circular::deg(x))) # round to nearest degree
    
    # control wrap-around effect in case of median at pi (180?)
    if (round_circ(circular::circular(fi + pi, modulo="2pi")) == 0) {
      fi <- pi
      AM<- 2*pi
    } else {
      AM <- round_circ(circular::circular(fi + pi, modulo="2pi"))
    }
    
    # control wrap-around effect
    if (circular::range.circular(circular::circular(IQR, modulo = "2pi"))< ((2*pi)/(2*(constant + (1/2)))) ) {
      if (fi<pi) {
        setAnti <- subset(IQR, IQR>=fi & IQR<=AM)
        setClock<- subset(IQR, IQR<=fi | IQR>=AM)
        QAnti   <- round_circ(circular::circular(max(setAnti), modulo="2pi"))
        Qc      <- QAnti - round_circ(circular::range.circular(circular::circular(IQR, modulo = "2pi")))
        QClock  <- round_circ(circular::circular(Qc, modulo="2pi"))
        
        grid <- c(Qc, QAnti)
        
        astart <- QAnti
        aend <- Qc
        
        d <- round_circ(circular::range.circular(circular::circular(IQR, modulo = "2pi")))
        
        fA <- round_circ(QAnti + d*constant)
        fC <- round_circ(QClock - d*constant)
        
        semicircleClock <- subset(as.vector(set_1), as.vector(set_1) <= fi | as.vector(set_1) >= AM)
        semicircleAnti  <- subset(as.vector(set_1), as.vector(set_1) >= fi & as.vector(set_1) <= AM)
        semicircleClock <- c(semicircleClock, QClock)
        semicircleAnti <- c(semicircleAnti, QAnti)
        
        if (fC < 0) {
          swc <- subset(semicircleClock, semicircleClock >= round_circ(circular::circular(fC, modulo="2pi")) | semicircleClock <= QClock)
          swc <- c(swc, QClock)
          whiskerC <- circular::range.circular(circular::circular(swc))
          wC <- QClock - whiskerC
          faroutClock  <- subset(semicircleClock, semicircleClock >= AM & semicircleClock < round_circ(circular::circular(fC, modulo="2pi")))
        } else if (fC >= 0 & QClock >= pi){
          swc <- subset(semicircleClock, semicircleClock>=fC)
          swc <- c(swc, QClock)
          wC <- min(swc)
          faroutClock <- subset(semicircleClock, semicircleClock>=AM & semicircleClock<fC)
        } else if (fC >= 0 & QClock < pi){
          swc <- subset(semicircleClock, semicircleClock>=fC)
          swc <- c(swc, QClock)
          wC <- min(swc)
          faroutClock <- subset(semicircleClock, semicircleClock>=AM | semicircleClock<fC)
        }
        
        swa <- subset(semicircleAnti, semicircleAnti<=fA)
        swa <- c(swa, QAnti)
        wA <- max(swa)
        
        faroutAnti  <- subset(semicircleAnti, semicircleAnti>fA)
      } else if (fi == pi) {
        setAnti <- subset(IQR, IQR>=fi & IQR<=2*pi)
        setClock<- subset(IQR, IQR<=fi | IQR>=0)
        QAnti   <- round_circ(circular::circular(max(setAnti), modulo="2pi"))
        Qc      <- QAnti - round_circ(circular::range.circular(circular::circular(IQR, modulo = "2pi")))
        QClock  <- round_circ(circular::circular(Qc, modulo="2pi"))
        
        grid <- c(Qc, QAnti)
        
        astart <- QAnti
        aend <- Qc
        
        d <- round_circ(circular::range.circular(circular::circular(IQR, modulo = "2pi")))
        fA<- round_circ(QAnti + d*constant)
        fC<- round_circ(QClock - d*constant)
        
        semicircleClock <- subset(as.vector(set_1),as.vector(set_1)<=fi | as.vector(set_1)>= 0)
        semicircleAnti <- subset(as.vector(set_1),as.vector(set_1)>=fi & as.vector(set_1)<= 2*pi)
        
        semicircleClock <- c(semicircleClock, QClock)
        semicircleAnti <- c(semicircleAnti, QAnti)
        if (fC < 0) {
          swc <- subset(semicircleClock, semicircleClock >= round_circ(circular::circular(fC, modulo="2pi")) | semicircleClock <= QClock)
          swc <- c(swc, QClock)
          whiskerC <- circular::range.circular(circular::circular(swc))
          wC <- QClock - whiskerC
          faroutClock  <- subset(semicircleClock, semicircleClock>=0 & semicircleClock<round_circ(circular::circular(fC, modulo="2pi")))
        } else if (fC>=0){
          swc <- subset(semicircleClock, semicircleClock>=fC)
          swc <- c(swc, QClock)
          wC <- min(swc)
          faroutClock <- subset(semicircleClock, semicircleClock>=0 & semicircleClock<fC)
        }
        swa <- subset(semicircleAnti, semicircleAnti<=fA)
        swa <- c(swa, QAnti)
        wA <- max(swa)
        faroutAnti  <- subset(semicircleAnti, semicircleAnti>fA)
      } else if (fi>pi) {
        setAnti <- subset(IQR, IQR>=fi | IQR<=AM)
        setClock<- subset(IQR, IQR<=fi & IQR>=AM)
        QClock   <- min(setClock)
        Qa      <- QClock+circular::range.circular(circular::circular(IQR, modulo = "2pi"))
        QAnti  <- round_circ(circular::circular(Qa, modulo="2pi"))
        
        grid <- c(QClock, Qa)
        
        astart <- QClock
        aend <- Qa
        
        # define the whiskers
        d <- circular::range.circular(circular::circular(IQR, modulo = "2pi"))
        fC<- round_circ(QClock - d*constant)
        fA<- round_circ(QAnti + d*constant)
        semicircleClock <- subset(as.vector(set_1),as.vector(set_1)<=fi & as.vector(set_1)>=AM)
        semicircleAnti <- subset(as.vector(set_1),as.vector(set_1)>=fi | as.vector(set_1)<=AM)
        
        semicircleClock <- c(semicircleClock, QClock)
        semicircleAnti <- c(semicircleAnti, QAnti)
        swc <- subset(semicircleClock, semicircleClock>=fC)
        swc <- c(swc, QClock)
        wC <- min(swc)
        faroutClock <- subset(semicircleClock, semicircleClock<fC )
        
        if (fA>2*pi ) {
          swa <- subset(semicircleAnti, semicircleAnti<= round_circ(circular::circular(fA, modulo="2pi")) | semicircleAnti>= round_circ(circular::circular(QAnti, modulo="2pi")))
          swa <- c(swa, QAnti)
          whiskerA <- circular::range.circular(circular::circular(swa))
          wA <- QAnti+whiskerA
          faroutAnti <- subset(semicircleAnti, semicircleAnti> circular::circular(fA, modulo = "2pi") & semicircleAnti <= AM )
        } else if (fA<=2*pi & QAnti>=pi) {
          swa <- subset(semicircleAnti, semicircleAnti<=fA)
          swa <- c(swa, QAnti)
          wA <- max(swa)
          faroutAnti <- subset(semicircleAnti, semicircleAnti>fA | semicircleAnti<= AM  )
        } else if (fA<=2*pi & QAnti<pi) {
          swa <- subset(semicircleAnti, semicircleAnti<=fA)
          swa <- c(swa, QAnti)
          wA <- max(swa)
          faroutAnti <- subset(semicircleAnti, semicircleAnti>fA & semicircleAnti<= AM  )
        }
      }
      
      
      # plot and print far out values
      faroutvalues1 <- c(faroutClock, faroutAnti)
      compare <- set_1
      faroutvalues2 <- compare[compare %in% faroutvalues1]
      faroutvalues <- circular::as.circular(circular::circular(faroutvalues2), modulo="2pi")
      farout<- as.matrix(faroutvalues2)
      colnames(farout) <- c("Far out values")
      
      if (H) {
        rgl::points3d(cos(faroutvalues)*tmp_seq[curr_seq], 
                      sin(faroutvalues)*tmp_seq[curr_seq], 
                      z_seq[curr_seq], 
                      cex=2, col="white", add = T)
        
        rgl::points3d(cos(faroutvalues)*tmp_seq[curr_seq], 
                      sin(faroutvalues)*tmp_seq[curr_seq], 
                      z_seq[curr_seq], 
                      cex=2, pch=8, add = T)
        }else {
        if (minimal) {
          rgl::points3d(cos(faroutvalues)*tmp_seq[curr_seq], 
                        sin(faroutvalues)*tmp_seq[curr_seq], 
                        z_seq[curr_seq], 
                        cex=2, size=5.5, pch=8, col="gray30", add = T)
          }else {
            rgl::points3d(cos(faroutvalues)*tmp_seq[curr_seq], 
                   sin(faroutvalues)*tmp_seq[curr_seq], 
                   z_seq[curr_seq], 
                   cex=2, size=5.5, pch=8, col="gray30", add = T)
        }
      }
    } else {
      # case when range(box)>= (360/2(c+1/2))
      if (fi<=pi) {
        setAnti <- subset(IQR, IQR>=fi & IQR<=AM)
        setClock<- subset(IQR, IQR<=fi | IQR>=AM)
        QAnti   <- round_circ(circular::circular(max(setAnti), modulo="2pi"))
        Qc      <- QAnti-circular::range.circular(circular::circular(IQR, modulo = "2pi"))
        QClock  <- round_circ(circular::circular(Qc, modulo="2pi"))
        
        grid <- c(Qc, QAnti)
        
        astart <- QAnti
        aend <- Qc
        
        # define the whiskers
        semicircleClock <- subset(as.vector(set_1),as.vector(set_1)<=fi | as.vector(set_1)>=AM)
        semicircleAnti <- subset(as.vector(set_1),as.vector(set_1)>=fi & as.vector(set_1)<=AM)
        
        semicircleClock <- c(semicircleClock, QClock)
        semicircleAnti <- c(semicircleAnti, QAnti)
        if (QClock <= pi){
          swc <- subset(semicircleClock, semicircleClock >= AM | semicircleClock <= QClock)
        }
        else if (QClock > pi){
          swc <- subset(semicircleClock, semicircleClock >= AM & semicircleClock <= QClock)
        }
        
        whiskerC <- circular::range.circular(circular::circular(swc))
        wC <- QClock-whiskerC
        
        swa <- subset(semicircleAnti, semicircleAnti<=AM)
        wA <- max(swa)
      } else if (fi>=pi) {
        setAnti <- subset(IQR, IQR>=fi | IQR<=AM)
        setClock<- subset(IQR, IQR<=fi & IQR>=AM)
        QClock   <- min(setClock)
        Qa      <- QClock+circular::range.circular(circular::circular(IQR, modulo = "2pi"))
        QAnti  <- round_circ(circular::circular(Qa, modulo="2pi"))
        
        grid <- c(QClock, Qa)
        
        astart <- QClock
        aend <- Qa
        
        semicircleClock <- subset(as.vector(set_1),as.vector(set_1)<=fi & as.vector(set_1)>=AM)
        semicircleAnti <- subset(as.vector(set_1),as.vector(set_1)>=fi | as.vector(set_1)<=AM)
        semicircleClock <- c(semicircleClock, QClock)
        semicircleAnti <- c(semicircleAnti, QAnti)
        swc <- subset(semicircleClock, semicircleClock>=AM)
        wC <- min(swc)
        
        if(QAnti<pi){
          swa <- subset(semicircleAnti, semicircleAnti<=AM & semicircleAnti >= QAnti)
        }
        else if(QAnti>pi){
          swa <- subset(semicircleAnti, semicircleAnti<=AM | semicircleAnti >= QAnti)
        }
        whiskerA <- circular::range.circular(circular::circular(swa))
        wA <- QAnti+whiskerA
      }
    }
    
    # draw the plot
    if (!minimal) {
      size_ctrl <- ifelse(scale_widths, size_val[curr_seq], 1)
      tmp <- tmp_seq[curr_seq]
      
      # main box
      grid.seq <- seq(grid[1], grid[2], by=0.005)
      ngrid <- length(grid.seq)
      for(i in 1:ngrid){
        pxy = draw.radial.line.xy(tmp - 0.1*size_ctrl, tmp + 0.1*size_ctrl, 
                                  center=c(0,0), grid.seq[i])
        rgl::plot3d(pxy$x, pxy$y, z_seq[curr_seq], 
               type = "l", col = plot_cols[curr_seq],lwd = 2,
               cex = .3, pch = 20, add=T)
      }
      
      # median
      pxy = draw.radial.line.xy(tmp - 0.095*size_ctrl, tmp + 0.095*size_ctrl, 
                                center = c(0,0), CTM)
      rgl::plot3d(pxy$x, pxy$y, z_seq[curr_seq], 
             type = "l", col = line_cols[curr_seq],lwd = 4,
             cex = .3, pch = 20, add=T)
      
      # arc outlines of box
      pxy = draw.arc.xy(0, 0, tmp + 0.1*size_ctrl, grid[1], grid[2])
      start = c(pxy$x1[1], pxy$y1[1], z_seq[curr_seq])
      end = c(pxy$x2[pxy$n], pxy$y2[pxy$n], z_seq[curr_seq])
      rgl::arc3d(start, end, center = c(0, 0, z_seq[curr_seq]),
            radius = tmp + 0.1*size_ctrl, 
            lwd = 2, add = T)
      
      pxy = draw.arc.xy(0, 0, tmp - 0.1*size_ctrl, grid[1], grid[2])
      start = c(pxy$x1[1], pxy$y1[1], z_seq[curr_seq])
      end = c(pxy$x2[pxy$n], pxy$y2[pxy$n], z_seq[curr_seq])
      rgl::arc3d(start,end,center = c(0, 0, z_seq[curr_seq]),
            radius = tmp - 0.1*size_ctrl,
            lwd = 2, add = T)
      
      # endpoint outlines of box
      pxy = draw.radial.line.xy(tmp - 0.1*size_ctrl, tmp + 0.1*size_ctrl, 
                                center=c(0,0), grid[1])
      rgl::plot3d(pxy$x, pxy$y, z_seq[curr_seq], 
             type = "l", col = 1, lwd = 2, cex = .3, add=T)
      
      pxy = draw.radial.line.xy(tmp - 0.1*size_ctrl, tmp + 0.1*size_ctrl, 
                                center=c(0,0), grid[2])
      rgl::plot3d(pxy$x, pxy$y, z_seq[curr_seq], 
             type = "l", col = 1, lwd = 2, cex = .3, add=T)
      
      # whiskers
      pxy = draw.arc.xy(0,0, tmp_seq[curr_seq], wC,QClock)
      start = c(pxy$x1[1], pxy$y1[1], z_seq[curr_seq])
      end = c(pxy$x2[pxy$n], pxy$y2[pxy$n], z_seq[curr_seq])
      rgl::arc3d(start,end,center = c(0, 0, z_seq[curr_seq]),
            radius = tmp_seq[curr_seq],
            lwd = 2, add = T)
      
      # endpoints of whiskers
      pxy = draw.radial.line.xy(tmp_seq[curr_seq]-0.05, 0.05 + tmp_seq[curr_seq],
                                center=c(0,0),wC)
      rgl::plot3d(pxy$x, pxy$y, z_seq[curr_seq], 
             type = "l", col = 1, lwd = 2,
             cex = .3, add=T)
      
      pxy = draw.arc.xy(0, 0, tmp_seq[curr_seq], wA,QAnti)
      start = c(pxy$x1[1], pxy$y1[1], z_seq[curr_seq])
      end = c(pxy$x2[pxy$n], pxy$y2[pxy$n], z_seq[curr_seq])
      rgl::arc3d(start,end,center = c(0, 0, z_seq[curr_seq]),
            radius = tmp_seq[curr_seq],
            lwd = 2, add = T)
      
      pxy = draw.radial.line.xy(tmp_seq[curr_seq]-0.05, 0.05 + tmp_seq[curr_seq],
                                center=c(0,0),wA)
      rgl::plot3d(pxy$x, pxy$y, z_seq[curr_seq], 
             type = "l", col = 1, lwd = 2, cex = .3, add=T)
    }else {
      # draw the quartile boxplot
      
      # whiskers
      pxy = draw.arc.xy(0, 0, tmp_seq[curr_seq], wC,QClock)
      start = c(pxy$x1[1], pxy$y1[1], z_seq[curr_seq])
      end = c(pxy$x2[pxy$n], pxy$y2[pxy$n], z_seq[curr_seq])
      rgl::arc3d(start,end,center = c(0, 0, z_seq[curr_seq]),
            radius = tmp_seq[curr_seq],
            col="gray90", lty=1, lwd = 2, add = T)
      
      pxy = draw.arc.xy(0, 0, tmp_seq[curr_seq], wA,QAnti)
      start = c(pxy$x1[1], pxy$y1[1], z_seq[curr_seq])
      end = c(pxy$x2[pxy$n], pxy$y2[pxy$n], z_seq[curr_seq])
      rgl::arc3d(start,end,center = c(0, 0, z_seq[curr_seq]),
            radius = tmp_seq[curr_seq],
            col="gray90", lty=1, lwd = 2, add = T)
      
      # "box" as radial line
      pxy = draw.arc.xy(0, 0, tmp_seq[curr_seq], grid[1], grid[2])
      start = c(pxy$x1[1],pxy$y1[1],z_seq[curr_seq])
      end = c(pxy$x2[pxy$n],pxy$y2[pxy$n],z_seq[curr_seq])
      # check if the two points are opposite each other (arc is not well defined)
      if(sum(start[1:2] + end[1:2]) < 1e-10){
        end = c(pxy$x2[pxy$n-1],pxy$y2[pxy$n-1],z_seq[curr_seq])
      }
      rgl::arc3d(start,end,center = c(0,0,z_seq[curr_seq]),
                 radius = tmp_seq[curr_seq],
                 col=plot_cols[curr_seq], lwd = lwd, add = T)
      
      # median point
      rgl::points3d(cos(CTM)*tmp_seq[curr_seq], 
                    sin(CTM)*tmp_seq[curr_seq], 
                    z_seq[curr_seq], 
                    cex=2, col=line_cols[curr_seq], add = T)
    }
    
    gradi <- (as.matrix(circular::deg(data[,1])))
    output <- as.matrix(cbind(data,gradi))
    colnames(output) <- c("Obs.Radians", "Ranking", "Obs.Degrees")
    ##print(output)
    
    # draw an arrow indicating the median
    if (draw_arrow){
      arrow.end <- c(cos(as.numeric(CTM))*0.7*torus.R, 
                     sin(as.numeric(CTM))*0.7*torus.R, 
                     z_seq[curr_seq])
      rgl::arrow3d(c(0,0,0), arrow.end, 
              width = 0.15, s = 0.15, col = arrow_cols[curr_seq], add = T)
    }
    
    ## output object
    out = list()
    if(exists("faroutvalues")==TRUE){
      if(length(faroutvalues)!=0){out$farout = faroutvalues}
      else{out$farout = c("no far out values detected")}
    }else{out$farout = c("no far out values detected")}
    
    out$constant = constant
    summaryStatistics = data.frame(
      CircularMedian = as.numeric(circular::deg(circular::circular(fi, modulo="2pi"))),
      CounterClockwiseHinge = as.numeric(circular::deg(circular::circular(QAnti, modulo="2pi"))),
      ClockwiseHinge = as.numeric(circular::deg(circular::circular(QClock, modulo="2pi"))),
      CounterClockwiseWhisker = as.numeric(circular::deg(circular::circular(wA, modulo="2pi"))),
      ClockwiseWhisker = as.numeric(circular::deg(circular::circular(wC, modulo="2pi")))
    )
    out$statistics = summaryStatistics
    summary_output[[curr_seq]] <- out
    
  }
  
  return(invisible(summary_output))
  
}
