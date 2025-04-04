# THIS CODE IS FULLY COMMENTED

#' @title plot_ellipse() function
#' @description Accepts parameters describing an ellipse and returns a simple plot
#' @param x The x coordinate of the ellipse's center point
#' @param y The y coordinate of the ellipse's center point
#' @param r The r value of the ellipse, which helps define its oblongness (x coordinate of the tilt vector)
#' @param s The s value of the ellipse, which helps define its oblongness (y coordinate of the tilt vector)
#' @param a The log area of the ellipse
#' @param z The height of the tilt vector, default value is 10
#' @export
plot_ellipse <- function(x,y,r,s,a,z=10) {

  # Type checking
  if (!("numeric" %in% class(x))) {stop("x must be numeric")}
  if (!("numeric" %in% class(y))) {stop("y must be numeric")}
  if (!("numeric" %in% class(r))) {stop("r must be numeric")}
  if (!("numeric" %in% class(s))) {stop("s must be numeric")}
  if (!("numeric" %in% class(a))) {stop("a must be numeric")}
  if (!("numeric" %in% class(z))) {stop("z must be numeric")}

  # Making plot coordinates
  ellipse <- make_ellipse_coords(x,y,r,s,a,z)

  # Plotting the ellipse
  ellipse_plot <- plot(ellipse$x,ellipse$y,xlab="",ylab="",asp=1,type="l")
  return(ellipse_plot)
}

#' @title plot_scenario() function
#' @description Accepts parameters describing an ellipse and cladogenetic scenario and returns a simple plot
#' @param x The x coordinate of the ancestral ellipse's center point
#' @param y The y coordinate of the ancestral ellipse's center point
#' @param r The r value of the ancestral ellipse, which helps define its oblongness (x coordinate of the tilt vector)
#' @param s The s value of the ancestral ellipse, which helps define its oblongness (y coordinate of the tilt vector)
#' @param a The log area of the ancestral ellipse
#' @param d The daughter configuration during cladogenesis (0,1), where 0 makes the left daughter D1 and 1 makes the right daughter D1
#' @param m The cladogenetic mode (0,1), where 0 is sympatric and 1 is allopatric
#' @param c The c value (not index) of the cladogenetic scenario's concentric circle --> ex. 0.5
#' @param h The h value (not index) of the cladogenetic scenario's direction line --> ex. pi/4, if required
#' @param z The height of the tilt vector, default value is 10
#' @param alpha The alpha parameter, default value is -3
#' @export
plot_scenario <- function(x,y,r,s,a,d,m,c,h,z=10,alpha=-3) {

  # Getting daughter info
  # We need x, y, and a for both the left and right daughters after cladogenesis
  # Based on the cladogenetic event input values
  x_left  = transform_node(char="x",daughter="left", d=d,m=m,c=c,h=h,x=x,y=y,r=r,s=s,a=a,alpha=alpha,z=z)
  y_left  = transform_node(char="y",daughter="left", d=d,m=m,c=c,h=h,x=x,y=y,r=r,s=s,a=a,alpha=alpha,z=z)
  a_left  = transform_node(char="a",daughter="left", d=d,m=m,c=c,h=h,x=x,y=y,r=r,s=s,a=a,alpha=alpha,z=z)
  x_right = transform_node(char="x",daughter="right",d=d,m=m,c=c,h=h,x=x,y=y,r=r,s=s,a=a,alpha=alpha,z=z)
  y_right = transform_node(char="y",daughter="right",d=d,m=m,c=c,h=h,x=x,y=y,r=r,s=s,a=a,alpha=alpha,z=z)
  a_right = transform_node(char="a",daughter="right",d=d,m=m,c=c,h=h,x=x,y=y,r=r,s=s,a=a,alpha=alpha,z=z)

  # Making plot
  # First, generating ellipse points for ancestor and both daughters
  ancestor <- make_ellipse_coords(x,y,r,s,a,z)
  left_daughter <- make_ellipse_coords(x_left,y_left,r,s,a_left,z)
  right_daughter <- make_ellipse_coords(x_right,y_right,r,s,a_right,z)

  # Determining graph boundaries based on the ellipses being plotted
  xlow <- min(min(ancestor$x),min(left_daughter$x),min(right_daughter$x))
  xhigh <- max(max(ancestor$x),max(left_daughter$x),max(right_daughter$x))
  ylow <- min(min(ancestor$y),min(left_daughter$y),min(right_daughter$y))
  yhigh <- max(max(ancestor$y),max(left_daughter$y),max(right_daughter$y))

  # Plotting all 3 ellipses, with the ancestor in black, the left daughter in blue, and the right daughter in red
  plot(ancestor$x,ancestor$y,xlab="",ylab="",asp=1,type="l",xlim=c(xlow,xhigh),ylim=c(ylow,yhigh))
  lines(left_daughter$x,left_daughter$y,col="blue")
  lines(right_daughter$x,right_daughter$y,col="red")
  ellipse_plot <- recordPlot()
  return(ellipse_plot)
}

#' @title plot_scenario_annotated() function
#' @description Accepts parameters describing an ellipse and cladogenetic scenario and returns an annotated plot using ggplot2
#' @param x The x coordinate of the ancestral ellipse's center point
#' @param y The y coordinate of the ancestral ellipse's center point
#' @param r The r value of the ancestral ellipse, which helps define its oblongness (x coordinate of the tilt vector)
#' @param s The s value of the ancestral ellipse, which helps define its oblongness (y coordinate of the tilt vector)
#' @param a The log area of the ancestral ellipse
#' @param d The daughter configuration during cladogenesis (0,1), where 0 makes the left daughter D1 and 1 makes the right daughter D1
#' @param m The cladogenetic mode (0,1), where 0 is sympatric and 1 is allopatric
#' @param c The c value (not index) of the cladogenetic scenario's concentric circle --> ex. 0.5
#' @param h The h value (not index) of the cladogenetic scenario's direction line --> ex. pi/4, if required
#' @param z The height of the tilt vector, default value is 10
#' @param alpha The alpha parameter, default value is -3
#' @export
plot_scenario_annotated <- function(x,y,r,s,a,d,m,c,h,z=10,alpha=-3) {

  # Getting daughter info -- see plot_scenario() for details
  x_left  = transform_node(char="x",daughter="left", d=d,m=m,c=c,h=h,x=x,y=y,r=r,s=s,a=a,alpha=alpha,z=z)
  y_left  = transform_node(char="y",daughter="left", d=d,m=m,c=c,h=h,x=x,y=y,r=r,s=s,a=a,alpha=alpha,z=z)
  a_left  = transform_node(char="a",daughter="left", d=d,m=m,c=c,h=h,x=x,y=y,r=r,s=s,a=a,alpha=alpha,z=z)
  x_right = transform_node(char="x",daughter="right",d=d,m=m,c=c,h=h,x=x,y=y,r=r,s=s,a=a,alpha=alpha,z=z)
  y_right = transform_node(char="y",daughter="right",d=d,m=m,c=c,h=h,x=x,y=y,r=r,s=s,a=a,alpha=alpha,z=z)
  a_right = transform_node(char="a",daughter="right",d=d,m=m,c=c,h=h,x=x,y=y,r=r,s=s,a=a,alpha=alpha,z=z)

  # Making plot
  ancestor <- make_ellipse_coords(x,y,r,s,a,z)
  left_daughter <- make_ellipse_coords(x_left,y_left,r,s,a_left,z)
  right_daughter <- make_ellipse_coords(x_right,y_right,r,s,a_right,z)
  # We also want an ellipse to represent the concentric circle used for the cladogenetic scenario
  c_ellipse <- make_ellipse_coords(x,y,r,s,a+log(c^2),z)

  # Determining graph boundaries based on the ellipses and the tilt point
  xlow <- min(min(ancestor$x),min(left_daughter$x),min(right_daughter$x),min(c_ellipse$x),r)
  xhigh <- max(max(ancestor$x),max(left_daughter$x),max(right_daughter$x),max(c_ellipse$x),r)
  ylow <- min(min(ancestor$y),min(left_daughter$y),min(right_daughter$y),min(c_ellipse$y),s)
  yhigh <- max(max(ancestor$y),max(left_daughter$y),max(right_daughter$y),max(c_ellipse$y),s)

  # Expanding the graph slightly, so annotations fit
  x_low_lim <- xlow - (xhigh - xlow) * 0.1
  x_high_lim <- xhigh + (xhigh - xlow) * 0.1
  y_low_lim <- ylow - (yhigh - ylow) * 0.1
  y_high_lim <- yhigh + (yhigh - ylow) * 0.1

  # d is the daughter configuration during cladogenesis (0,1), where 0 makes the left daughter D1 and 1 makes the right daughter D1
  # m is the cladogenetic mode (0,1), where 0 is sympatric and 1 is allopatric
  # c is the c value (not index) of the cladogenetic scenario's concentric circle --> ex. 0.5
  # h is the h value (not index) of the cladogenetic scenario's direction line --> ex. pi/4, if required

  # Gathering information for annotations
  if (d == 0) {left <- "d1"} else {left <- "d2"}
  if (d == 0) {right <- "d2"} else {right <- "d1"}
  if (m == 0) {mode <- paste0("sympatric, alpha=",alpha)} else {mode <- "allopatric"}

  # Annotation describing mode and daughter configuration
  top_text <- paste0(
    "Left Daughter: ",left," (blue)","\n",
    "Right Daughter: ",right," (red)","\n",
    "Mode: ",mode
  )

  # Annotation describing concentric circle, direction line, and log area of ancestor
  bottom_text <- paste0(
    "C: ",round(c,2),"\n",
    "H: ",round(h/pi,2)," x pi","\n",
    "Log Area (a): ",a
  )

  # Annotations for center point and tilt point, adjusted to fit in graph
  center_label <- paste0("Center: (",x,",",y,")")
  tilt_label <- paste0("Tilt Point: (",r,",",s,",",z,")")
  if (s > y) {
    vjust_center <- 1
    vjust_tilt <- 0}
  else {
    vjust_center <- 0
    vjust_tilt <- 1}
  # Direction line
  slope <- (y_right-y_left)/(x_right-x_left)
  intercept <- y_left - (x_left * slope)

  # Combining ellipse plots and annotations
  ellipse_plot <- ggplot() +
    lims(x=c(x_low_lim,x_high_lim),y=c(y_low_lim,y_high_lim)) +
    coord_fixed() +
    # Ancestor ellipse
    geom_polygon(data=ancestor,aes(x=x,y=y),color="black",fill=NA,linewidth=2) +
    # Direction line
    geom_abline(slope=slope,intercept=intercept,linetype="dashed") +
    # Left daughter ellipse
    geom_polygon(data=left_daughter,aes(x=x,y=y),color="blue",fill=NA) +
    # Right daughter ellipse
    geom_polygon(data=right_daughter,aes(x=x,y=y),color="red",fill=NA) +
    # Concentric circle
    geom_polygon(data=c_ellipse,aes(x=x,y=y),color="black",fill=NA,linetype="dashed") +
    # Tilt vector
    geom_segment(aes(x=x,xend=r,y=y,yend=s),linetype="dashed",color="black") +
    # Left daughter centroid
    geom_point(aes(x=x_left,y=y_left),color="blue") +
    # Right daughter centroid
    geom_point(aes(x=x_right,y=y_right),color="red") +
    # Labels
    annotate("label",x=x_low_lim,y=y_high_lim,label=top_text,hjust="inward",vjust="inward") +
    annotate("label",x=x_high_lim,y=y_low_lim,label=bottom_text,hjust="inward",vjust="inward") +
    annotate("label",x=x,y=y,label=center_label,vjust=vjust_center,hjust=0.5) +
    annotate("label",x=r,y=s,label=tilt_label,vjust=vjust_tilt,hjust=0.5) +
    # Ancestor centroid
    geom_point(aes(x=x,y=y)) +
    # Tilt point
    geom_point(aes(x=r,y=s)) +
    theme_bw()
  return(ellipse_plot)
}