# adapted from CircularBoxplot from Buttarazzi D., Pandolfo G., Porzio G.C. (2018). A boxplot for circular data, Biometrics.
# the original function is also available in package 'bpDir' as 'CircularBoxplot'
# modified by josh berlinski

# TODO:
#   - return (invisible) list object similar to original function
#   - add example in documentation, potentially make a nicer description
#   - (low prio) make axis labels look nicer by default

#' Create grouped circular boxplot
#'
#' @description Given a named list of circular objects, create a grouped circular boxplot.
#'
#' @param data_in List of circular objects to be plotted
#' @param template One of "degrees", "radians", "geographics", or NULL
#' @param place If template is NULL, either "outside" or "inside" denoting where the axis should be drawn
#' @param units If template is NULL, units to use on the drawn axis
#' @param marg character specifying the plot region as either "large" (default) or "small".
#' @param shrink Numeric specifying the factor by which to scale the plot. Numbers less than 1 will increase the size.
#' @param H Logical indicating if each data point should be drawn outside the hinges of the boxplot
#' @param stack Logical indicating if drawn points should be stacked
#' @param constant Numeric specifying the multiplicative factor determining how far whiskers extend from box. A value of
#' "optimal" will choose values based on a von Mises distribution (see Buttarazzi et al. 2018)
#' @param rad_shift Scalar indicating the distance between each boxplot
#' @param lwd Scalar indicating the width of boxplot lines
#' @param plot_cols Vector with the same length as `data_in`, specifying the color of the boxplot
#' @param line_cols Vector with the same length as `data_in`, specifying the color of the median lines
#' @param arrow_cols Vector with the same length as `data_in`, specifying the color of the arrows pointing to the
#' medians. Default value is the same as line_cols
#' @param legend_pos String indicating where the legend should be drawn, or "none" for no legend
#' @param legend_title String indicating legend for title, if necessary
#' @param draw_arrow Logical specifying if arrows pointing to each median should be drawn
#' @param ordinal Logical, if `template` is NULL and `units` is "geographics", should the ordinal directions also be drawn?
#' @param template_options If template is "custom", a list containing named elements "ax_labels", "lab_coord", "shift",
#' and "gridlines." `ax_labels` is a character vector of axis labels to be plotted at the locations `lab_coord`, a
#' numeric vector in degrees counterclockwise starting from 0. `shift` is a numeric value that controls how far outside
#' the last boxplot the axis should be drawn, and `gridlines` is a logical indicating if dashed lines should be drawn at
#' the specified axis labels.
#' @param minimal Logical. If true, a simple colored line and point will replace the box and median line. Radial lines
#' at fence points will also not be drawn
#' @param scale_widths Logical, should the width of each boxplot be scaled based on it's distance from the center?
#' @export
#' @author Josh Berlinski
#' @author Davide Buttarazzi
#' @importFrom graphics legend text par points
GroupedCircularBoxplot <- function(
  data_in,
  template = "degrees",
  place = "none",
  units = "degrees",
  marg = "large",
  shrink = 1,
  H = FALSE,
  stack = FALSE,
  constant = "optimal",
  rad_shift = 0.5,
  lwd = 1,
  plot_cols = RColorBrewer::brewer.pal(8, "Set2"),
  line_cols = RColorBrewer::brewer.pal(8, "Dark2"),
  arrow_cols = line_cols,
  legend_pos = "topleft",
  legend_title = NULL,
  draw_arrow = TRUE,
  ordinal = FALSE,
  template_options = NULL,
  minimal = FALSE,
  scale_widths = FALSE
) {

  # if only a circular vector is passed as data in, make it a list and don't plot a legend
  if (circular::is.circular(data_in)) {
    data_in <- list(a = data_in)
    legend_pos <- "none"
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


  shift_val <- (0:(n_seq - 1)) * rad_shift
  size_val <- rev(1:n_seq * rad_shift)


  # let A be a list of vectors, each to be plotted

  for (curr_seq in 1:n_seq) {
    A <- data_in[[curr_seq]]
    delta <- shift_val[curr_seq]

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
    # oldpar <- par(no.readonly = TRUE)
    # on.exit(par(oldpar))

    if (marg == "small")
      par(oma=c(0,0,0,0))
    else if (marg == "large")
      par(mai=c(0.0,0.0,0,0))

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

    ##drawing the template
    if (curr_seq == 1) {
      circular::plot.circular(
        circular::circular(NA, modulo = "2pi"),
        cex=0.5, axes=FALSE, shrink=shrink,
        template=NULL, control.circle=circular::circle.control(col="gray60", lty=2, lwd=0.5),
        xlim = c(- 1 - max(shift_val), 1 + max(shift_val))
      )

      if (legend_pos != "none") {
        if (is.null(names(data_in)))
          names(data_in) <- 1:length(data_in)
        legend(
          legend_pos,
          legend = names(data_in),
          fill = plot_cols[1:length(data_in)],
          title = legend_title
        )
      }

      if (is.null(template)) {
        rad_grid <- seq(0, 315, by = 45)
        lab_coord <- cbind(
          cos(circular::rad(circular::circular(rad_grid))),
          sin(circular::rad(circular::circular(rad_grid)))
        )

        # this is awkward
        if (place == "outside" && units == "degrees") {
          plotrix::draw.arc(0, 0, 1.15, 0, 2 * pi, col = "burlywood4")
          bline <- c(1.09, 1.15)
          aline <- c(0, 1.08)
          ax_labels <- as.character(rad_grid)
          tmult <- 1.23
        } else if (place == "inside" && units == "degrees") {
          plotrix::draw.arc(0, 0, 0.7, 0, 2 * pi, col = "burlywood4")
          bline <- c(0.64, 0.7)
          aline <- c(0, 0.93)
          ax_labels <- as.character(rad_grid)
          tmult <- 0.5
        } else if (place == "outside" && units == "radians") {
          plotrix::draw.arc(0, 0, 1.15, 0, 2 * pi, col = "burlywood4")
          bline <- c(1.09, 1.15)
          aline <- c(0, 1.08)
          ax_labels <- c(expression(0, frac(pi, 4), frac(pi, 2), frac(3*pi, 4), pi, frac(5*pi, 4), frac(3*pi, 2), frac(7*pi, 4)))
          tmult <- 1.3
        } else if (place == "inside" && units == "radians") {
          plotrix::draw.arc(0, 0, 0.7, 0, 2 * pi, col = "burlywood4")
          bline <- c(0.64, 0.7)
          aline <- c(0, 0.93)
          ax_labels <- c(expression(0, frac(pi, 4), frac(pi, 2), frac(3*pi, 4), pi, frac(5*pi, 4), frac(3*pi, 2), frac(7*pi, 4)))
          tmult <- 0.45
        } else if (place == "outside" && units == "geographics") {
          if (ordinal) {
            ax_labels <- c("E", "NE", "N", "NW", "W", "SW", "S", "SE")
            rad_grid <- c(0, 45, 90, 135, 180, 225, 270, 315)
          } else {
            ax_labels <- c("E", "N", "W", "S")
            rad_grid <- c(0, 90, 180, 270)
          }
          tmult <- 1
          lab_coord <- cbind(
            cos(circular::rad(circular::circular(rad_grid))),
            sin(circular::rad(circular::circular(rad_grid)))
          )
          shift <- max(shift_val) + 1
          coord_offset <- cbind(
            shift * cos(circular::rad(circular::circular(rad_grid))),
            shift * sin(circular::rad(circular::circular(rad_grid)))
          )
          lab_coord <- lab_coord + coord_offset

          for (i in seq_along(rad_grid)) {
            if (i < length(rad_grid))
              plotrix::draw.arc(0, 0, shift + 1, deg1 = rad_grid[i] + 2, deg2 = rad_grid[i + 1] - 2, col = 1, lwd = 2, lty = 1)
            else
              plotrix::draw.arc(0, 0, shift + 1, deg1 = rad_grid[i] + 2, deg2 = 358, col = 1, lwd = 2, lty = 1)
          }
          text(
            tmult * circular::circular(lab_coord[, 1]), tmult * circular::circular(lab_coord[, 2]),
            labels = ax_labels, col = "black", cex = 1
          )
        }

        if (units != "geographics" && !is.null(units)) {
          for (i in rad_grid) {
            plotrix::draw.radial.line(aline[1], aline[2], center = c(0, 0), i, col = "azure2", lty = 2)
            plotrix::draw.radial.line(bline[1], bline[2], center = c(0, 0), i, col = "burlywood4", lty = 2)
          }
          text(
            tmult * circular::circular(lab_coord[, 1]), tmult * circular::circular(lab_coord[, 2]),
            labels = ax_labels, col = "burlywood4", cex = 0.6
          )
        }
      } else if (template == "custom") {
        ax_labels <- template_options$ax_labels
        rad_grid <- template_options$lab_coord
        tmult <- 1 #TODO: let there be an option
        shift_mod <- template_options$shift
        lab_coord <- cbind(
          cos(circular::rad(circular::circular(rad_grid))),
          sin(circular::rad(circular::circular(rad_grid)))
        )
        shift <- max(shift_val) + 1
        coord_offset <- cbind(
          (shift + shift_mod)* cos(circular::rad(circular::circular(rad_grid))),
          (shift + shift_mod) * sin(circular::rad(circular::circular(rad_grid)))
        )
        lab_coord <- lab_coord + coord_offset

        if (template_options$gridlines) {
          for (i in seq_along(rad_grid))
            plotrix::draw.radial.line(1, shift, center = c(0, 0), deg = rad_grid[i], col = "gray60", lty = 2, lwd = 0.5)
        }

        text(
          tmult * circular::circular(lab_coord[, 1]), tmult * circular::circular(lab_coord[, 2]),
          labels = ax_labels, col = "black", cex = 1
        )

      } else {
        lab_coord <- cbind(
          cos(circular::rad(circular::circular(c(0,90,180,270)))),
          sin(circular::rad(circular::circular(c(0,90,180,270))))
        )
        if(template=="degrees") {
          ax_labels <- c("0","90","180","270")
          tmult <- 0.82
        } else if(template=="radians") {
          ax_labels <- c(expression(0,frac(pi,2),pi,frac(3*pi,2)))
          tmult <- 0.65
        } else if(template=="geographics") {
          ax_labels <- c("E","N","W","S")
          tmult <- 0.82
        }

        text(
          tmult * circular::circular(lab_coord[,1], units = "radians"), tmult * circular::circular(lab_coord[,2], units = "radians"),
          labels = ax_labels, cex=0.6
        )
      }
    } else {
      plotrix::draw.circle(0,0,1 + delta,border="gray60", lty=2, lwd = 0.5)
    }

    ##drawing the plot
    if (H)
      circular::points.circular(circular::circular(set_1), cex=0.75, start.sep = delta)

    # TODO: what does this do, exactly?
    # it plots points over everything within the IQR, probably unnecessary,
    # gets drawn over by the box
    # points(circular::circular(IQR, modulo = "2pi"), cex=1.1, col="white", start.sep = delta)

    # is this bad?
    round_circ <- function(x) circular::rad(round(circular::deg(x))) # round to nearest degree

    ## controlling wrap-around effect in case of median at pi (180?)
    if (round_circ(circular::circular(fi + pi, modulo="2pi")) == 0) {
      fi <- pi
      AM<- 2*pi
    } else {
      AM <- round_circ(circular::circular(fi + pi, modulo="2pi"))
    }

    ## controlling wrap-around effect
    # im not touching this, its a (working) mess
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

        ## defining the whiskers
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


      ## plotting and printing far out values
      faroutvalues1 <- c(faroutClock, faroutAnti)
      compare <- set_1
      faroutvalues2 <- compare[compare %in% faroutvalues1]
      faroutvalues <- circular::as.circular(circular::circular(faroutvalues2), modulo="2pi")
      farout<- as.matrix(faroutvalues2)
      colnames(farout) <- c("Far out values")

      if (H) {
        circular::points.circular(faroutvalues, cex=0.8, col="white", start.sep = delta)
        circular::points.circular(faroutvalues, cex=0.8, pch=8, start.sep = delta)
      } else {
        if (stack)
          circular::points.circular(faroutvalues, cex=0.6, stack=stack, bins=500, sep=0.1, pch=8, start.sep = delta)
        else
          circular::points.circular(faroutvalues, cex=0.6, stack=stack, pch=8, start.sep = delta)
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

        ## defining the whiskers
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
    ### 1 draw box
    if (!minimal) {
      size_ctrl <- ifelse(scale_widths, sqrt(size_val[curr_seq]), 1)
      tmp <- 1 + delta

      plotrix::drawSectorAnnulus(
        angle1 = grid[1],
        angle2 = grid[2],
        radius1 = tmp - 0.1*size_ctrl,
        radius2 = tmp + 0.1*size_ctrl,
        col = plot_cols[curr_seq]
      )
      # draw median
      plotrix::draw.radial.line(tmp - 0.095*size_ctrl, tmp + 0.095*size_ctrl, center = c(0,0), CTM, col = line_cols[curr_seq], lwd = lwd)
      ### 2 draw arc outline of box
      plotrix::draw.arc(0, 0, tmp + 0.1*size_ctrl, grid[1], grid[2], col=1, lwd=lwd)
      plotrix::draw.arc(0, 0, tmp - 0.1*size_ctrl, grid[1], grid[2], col=1, lwd=lwd)
      ### 3
      plotrix::draw.radial.line(tmp - 0.1*size_ctrl, tmp + 0.1*size_ctrl, center=c(0,0), grid[1], col=1, lwd=lwd)
      plotrix::draw.radial.line(tmp - 0.1*size_ctrl, tmp + 0.1*size_ctrl, center=c(0,0), grid[2], col=1, lwd=lwd)
      ### 4
      plotrix::draw.arc(0,0,1 + delta,wC,QClock,col=1,lwd=lwd, lty=1)
      plotrix::draw.radial.line(0.95 + delta,1.05 + delta,center=c(0,0),wC,col=1,lwd=lwd)
      ### 5
      plotrix::draw.arc(0,0,1 + delta,wA,QAnti,col=1,lwd=lwd, lty=1)
      plotrix::draw.radial.line(0.95 + delta,1.05 + delta,center=c(0,0),wA,col=1,lwd=lwd)
    } else {
      plotrix::draw.arc(0,0,1 + delta,wC,QClock,col=1,lwd=0.25*lwd, lty=1)
      plotrix::draw.arc(0,0,1 + delta,wA,QAnti,col=1,lwd=0.25*lwd, lty=1)
      plotrix::draw.arc(0, 0, 1 + delta, grid[1], grid[2], col=plot_cols[curr_seq], lwd=lwd)
      circular::points.circular(CTM, cex=0.75, start.sep = delta, col=line_cols[curr_seq])
    }

    gradi <- (as.matrix(circular::deg(data[,1])))
    output <- as.matrix(cbind(data,gradi))
    colnames(output) <- c("Obs.Radians", "Ranking", "Obs.Degrees")
    ##print(output)

    ## drawing an arrow indicating the median
    # med_color <- ifelse(col_type == "fill", 1, line_cols[curr_seq])
    # arrows.circular(CTM, 0.78, col = line_cols[curr_seq], angle = 30, length = 0.1, )
    if (draw_arrow)
      circular::arrows.circular(CTM, 0.7, col = arrow_cols[curr_seq], angle = 30, length = 0.1, lwd = lwd)

    ## output object
    out = list()
    if(exists("faroutvalues")==TRUE){
      if(length(faroutvalues)!=0){out$farout = faroutvalues}
      else{out$farout = c("no far out values detected")}
    }
    else{out$farout = c("no far out values detected")}

    out$constant = constant

    summaryStatistics = data.frame(
      CircularMedian = as.numeric(circular::deg(circular::circular(fi, modulo="2pi"))),
      CounterClockwiseHinge = as.numeric(circular::deg(circular::circular(QAnti, modulo="2pi"))),
      ClockwiseHinge = as.numeric(circular::deg(circular::circular(QClock, modulo="2pi"))),
      CounterClockwiseWhisker = as.numeric(circular::deg(circular::circular(wA, modulo="2pi"))),
      ClockwiseWhisker = as.numeric(circular::deg(circular::circular(wC, modulo="2pi")))
    )

    out$statistics = summaryStatistics

  }
  points(0,0,pch=21, bg=1,  cex=1.1) # redraw center point, for fun
  return(invisible(out))

}
