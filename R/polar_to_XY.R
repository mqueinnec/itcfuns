#' Polar to XY coordianates
#' Convert polar coordinates to cartesian coordinates. Used for example if trees were mapped form plot center using distance and azimuth
#' @param azimuth - Numeric vector of azimuth from plot center
#' @param distance - Numeric vector of distance from plot center
#' @param xcenter - Cartesian X coordinate of plot center
#' @param ycenter - Cartesian Y coordinate of plot center
#' @param azimuth_offset - Optional. If the tree is leaning and measurement was from plot center to the stem, azimuth between stem and tree top.
#' @param distance_offset - Optional. If the tree is leaning and measurement was from plot center to the stem, distance between stem and tree top.
#'
#' @export



polar_to_XY <- function(azimuth,
                      distance,
                      xcenter,
                      ycenter,
                      azimuth_offset,
                      distance_offset) {
  #angle = azimuth - 90
  #angle[angle<0] <- 360 + angle[angle<0]
  angle = azimuth * pi/180
  #Convert to radians
  #angle = angle*pi/180
  angle = 2*pi - (angle - pi/2)
  x = xcenter + distance * cos(angle)
  y = ycenter + distance * sin(angle)

  # Apply offset
  if(!missing(azimuth_offset) & !missing(distance_offset)){
    angle_offset = azimuth_offset * pi/180
    angle_offset = 2*pi - (angle_offset - pi/2)

    x_offset = distance_offset * cos(angle_offset)
    y_offset = distance_offset * sin(angle_offset)
  }else{
    x_offset = 0
    y_offset = 0
  }

  data.frame(X = x + x_offset, Y = y + y_offset)
}
