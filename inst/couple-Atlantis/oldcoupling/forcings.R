temp <- rmatio::read.mat("inst/couple-Atlantis/Prydz_forcing/temp.mat")

str(temp)

## get a BGM and read it
bfile <- bgmfiles::bgmfiles("antarctica_28")
library(rbgm)
bgm <- bgmfile(bfile)
## we need the unsullied boxes to identify points inside them
box_bgm <- boxSpatial(bgm)
face_bgm <- faceSpatial(bgm)

temp <- rmatio::read.mat("inst/couple-Atlantis/Prydz_forcing/temp.mat")
u <- rmatio::read.mat("inst/couple-Atlantis/Prydz_forcing/u.mat")
t_len <- dim(temp$temp)[3]
#time=ncvar_get(nc_t,"ocean_time")
t_start = 1199145600  ## as.numeric(as.POSIXct("2008-01-01 00:00:00", tz = "GMT"))
dt <- 86400 

t_stop = t_start + t_len*86400
t_tot=seq(t_start,t_stop,dt)

label <- "antarctica_28"
#filename

f.temp = sprintf("%s_%s.nc", label,  "temp")
f.saln = sprintf("%s_%s.nc", label,  "salt")
nboxes <- nrow(box_bgm)
boxes <- seq_len(nboxes)
nlayers <- dim(temp$temp)[2]
library(ncdf4)
#define dimensions
dimb <- ncdim_def("b", "", boxes, create_dimvar = FALSE)
dimz <- ncdim_def("z", "" , seq_len(nlayers), create_dimvar=FALSE)
dimt <- ncdim_def("t1","", t_tot, unlim=TRUE)#,create_dimvar=FALSE)

#create variables
#NB!!!!!! Unlimited rec needs to be on the right - otherwise the program complains!
var.t <- ncvar_def("t", "seconds since 2008-01-01 00:00:00 +10",dimt,0,prec="double")
var.temp <-  ncvar_def("temperature", "DegC", list(dimz, dimb, dimt), 0, prec = "double")
var.saln <- ncvar_def("salinity","psu",list(dimz,dimb,dimt),0,prec="double")



#create file
nc_temp=nc_create(f.temp,list(var.t,var.temp))
nc_saln=nc_create(f.saln,list(var.t,var.saln))

#assign global attributes to temp file
ncatt_put(nc_temp,0,"title", sprintf("Temperature file %s", label))
ncatt_put(nc_temp,0,"geometry", basename(bfile))
ncatt_put(nc_temp,0,"parameters","")


#assign global attributes to saln file
ncatt_put(nc_saln,0,"title", sprintf("Salinity file %s", label))
ncatt_put(nc_saln,0,"geometry", basename(bfile))
ncatt_put(nc_saln,0,"parameters","")

#assign attributes to variables
ncatt_put(nc_temp,var.t,"dt",86400,prec="double")
ncatt_put(nc_saln,var.t,"dt",86400,prec="double")

#assign variables to file
ncvar_put(nc_temp, var.temp, temp$temp)
ncvar_put(nc_saln, var.saln, )

nc_close(nc_temp)
nc_close(nc_saln)
}