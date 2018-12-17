exe <- "/home/shared/svn/HydroConstruct/trunk/HydroConstruct"
mass <- "/home/shared/data/Atlantis/HydroCoupling/MichaelSumner/file26522491e32_mass.nc"
uv <- "/home/shared/data/Atlantis/HydroCoupling/MichaelSumner/file265238d86d50_transport.nc"
prm <- here::here("inst/couple-Atlantis/HydroConstructRun/run-hydroconstruct.R")
cmd <- glue::glue("{exe} -f {uv} -t {mass} -s {mass} -r {prm}")

system(cmd)
