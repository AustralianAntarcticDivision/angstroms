

## create NetCDF templates for Atlantis hydro
library(ncdf4)

### FOR TRANSPORT NC FILE
transp_filename <- sprintf("%s_transport.nc", tempfile())
mass_filename  <- sprintf("%s_mass.nc", tempfile())
transp_params <- ""
mass_params <- ""

bgmfilepath <- bgmfiles::bgmfiles("antarctica_28")
library(angstroms)



create_transport(transp_filename, model_title = "Transport file Antarctica_28", bgmfilepath = bgmfilepath, 
                              bgmlevels = atlantis_depths, time_steps = time_steps)

create_mass(mass_filename, model_title = "Mass file Antarctica_28", bgmfilepath = bgmfilepath, 
                              bgmlevels = atlantis_depths, time_steps = time_steps)


nctransp <- ncdf4::nc_open(transp_filename, write = TRUE)
transport <- array(rnorm(9 * 90 * 10), c(9, 90, 10))
ncvar_put(nctransp, "transport", transport, count = dim(transport))
nc_close(nctransp)

ncmass <- ncdf4::nc_open(mass_filename, write = TRUE)

salinity <- temperature <- vertical_flux <- array(rnorm(9 * 28 * 10), c(9, 28, 10))
dm <- dim(vertical_flux)

#assign variables to file
ncvar_put(ncmass, "verticalflux", vertical_flux, count = dm)
ncvar_put(ncmass, "salinity", salinity, count = dm)
ncvar_put(ncmass, "temperature",temperature, count= dm)
nc_close(ncmass)

