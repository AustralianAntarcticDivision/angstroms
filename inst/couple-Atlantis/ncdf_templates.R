

## create NetCDF templates for Atlantis hydro
library(ncdf4)

### FOR TRANSPORT NC FILE
filename = sprintf("%s_transport.nc", tempfile())

model_title <- "Transport file [placeholder]"

transp_params <- ""




bgmfilepath <- bgmfiles::bgmfiles("antarctica_28")


nc_transp <- create_transport(filename, model_title = "Transport file Antarctica_28", bgmfilepath = bgmfilepath, 
                              bgmlevels = 1:9, time_steps = Sys.time() + (1:10) * 24 * 3600)



transport <- array(0, c(9, 90, 10))
## count:  
## If not specified and the variable does NOT have an unlimited dimension, the entire variable is written. 
## If the variable has an unlimited dimension, this argument must be specified.
ncvar_put(nc_transp, "transport", transport, count = dim(transport))
ncvar_get(nc_transp, "time")

nc_close(nc_transp)


### For T, S, Vertical Flux NC file
filename="RM_NEUS_variables_2010_20180325_fix.nc"

#define dimensions
timedim=ncdim_def("time", "", 1:length(t_tot), unlim=T, create_dimvar = F) #as.double(t_tot)
leveldim=ncdim_def("level", "", 1:nlevel, create_dimvar = F)
boxesdim=ncdim_def("boxes", "", 1:nboxes, create_dimvar = F)

#create variables
#NB!!!!!! Unlimited rec needs to be on the right - otherwise R complains!
#origMissVal_ex=0.0
var.time=ncvar_def("time","seconds since 1964-01-01 00:00:00 +10",timedim,prec="double")
var.box=ncvar_def("boxes", "", boxesdim, longname="Box IDs", prec='integer')
var.lev=ncvar_def("level","",leveldim,longname="layer index; 1=near surface; positice=down" ,prec="integer")
var.vertflux=ncvar_def("verticalflux","m3/s",list(leveldim, boxesdim, timedim),-999,longname="vertical flux averaged over floor of box",prec="float")
var.temp=ncvar_def("temperature","degree_C",list(leveldim, boxesdim, timedim),-999,longname="temperature volume averaged",prec="float")
var.salt=ncvar_def("salinity","psu",list(leveldim,boxesdim,timedim),-999,longname="salinity volume averaged",prec="float")

nc_varfile=nc_create(filename,list(var.time,var.box, var.lev, var.salt, var.temp, var.vertflux))

#assign global attributes to file
ncatt_put(nc_varfile,0,"title","Box averaged properties file, NEUS")
ncatt_put(nc_varfile,0,"geometry","neus_tmerc_RM.bgm")
ncatt_put(nc_varfile,0,"parameters","")

#assign attributes to variables
ncatt_put(nc_varfile,var.time,"dt",86400,prec="double")

#assign variables to file
ncvar_put(nc_varfile,var.vertflux,vertical_flux, count=c(nlevel,nboxes, ntimes))
ncvar_put(nc_varfile,var.time,t_tot)
ncvar_put(nc_varfile,var.lev,atl.level)
ncvar_put(nc_varfile,var.salt,salinity, count=c(nlevel,nboxes, ntimes))
ncvar_put(nc_varfile,var.temp,temperature, count=c(nlevel,nboxes, ntimes))
ncvar_put(nc_varfile,var.box,box.boxes)

nc_close(nc_varfile)

