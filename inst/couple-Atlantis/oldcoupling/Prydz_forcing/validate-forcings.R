library(shinyrAtlantis)
fp <- "~/Git/angstroms/inst/couple-Atlantis/Prydz_forcing"
#aforce <- "/home/shared/data/Atlantis/EastAntarctica_Atlantis/run_models/EastAntarctica_2016-12-08/SOhydrodummy.nc"
salinity.file    <- file.path(fp, "SO28salt.nc"  )    
temperature.file <- file.path(fp, "SO28temp.nc"   )  
bgm.file         <- bgmfiles::bgmfiles("antarctica_28")
exchange.file <- file.path(fp, "SOhydrodummy.nc")
cum.depth <- c(0,5,10,20,50,100,200,3000)  # cumulative water layer depths

input.object <- make.sh.forcings.object(
  bgm.file         = bgm.file,
  exchange.file    = exchange.file,
  cum.depth        = cum.depth,
  temperature.file = temperature.file,
  salinity.file    = salinity.file
)

sh.forcings(input.object)
