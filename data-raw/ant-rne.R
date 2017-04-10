
antarctica <- subset(rnaturalearth::countries110, sovereignt == "Antarctica", select = "sov_a3")
devtools::use_data(antarctica)
