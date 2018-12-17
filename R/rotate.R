# This is to remind you that ROMS/TOMS is formulated in general horizontal curvilinear coordinates. It is possible to formulate any application in Cartesian, spherical, or polar coordinates. You can activate a curvilnear application by turning on CPP option CURVILINEAR. Curvilinear applications include extra terms in the advection and any vector or tensor needs to be rotated to model generic coordinates (XI,ETA). 
# 
# In curvilinear applications, the coordinate rotation angle is found in the ROMS grid and history NetCDF files. See NetCDF variable angle. It is defined as the counterclockwise angle (radians) between the XI-axis and true EAST at RHO-points. 
# 
# Therefore, to transform between (XI,ETA) coordinates to (LON,LAT) coordinates, vectors and tensors need to be rotated according to 
# 
# Code:
#   u(LON,LAT)=u(XI,ETA)*cos(angle(i,j))-v(XI,ETA)*sin(angle(i,j))
# v(LON,LAT)=u(XI,ETA)*sin(angle(i,j))+v(XI,ETA)*cos(angle(i,j))
# 
# 
# conversely,
# 
# Code:
#   u(XI,ETA)=u(LON,LAT)*cos(angle(i,j))+v(LON,LAT)*sin(angle(i,j))
# v(XI,ETA)=u(LON,LAT)*sin(angle(i,j))-v(LON,LAT)*cos(angle(i,j))
#' @name romsrotate
#' @export
romsrotate <- function(uv, angle) {
  u <- uv[, 1] * cos(angle) - uv[,2] * sin(angle)
  v <- uv[, 1] * sin(angle) + uv[,2] * cos(angle)
  cbind(u, v)
}