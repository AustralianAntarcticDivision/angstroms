library(shiny)
library(leaflet)

library(angstroms)

library(rworldmap)
library(tabularaster)
data(countriesLow)
longlat <- romscoords(roms_path)
bound <- boundary(longlat)
projection(bound) <- "+init=epsg:4326"
map <- raster::intersect(bound, countriesLow)
map <- romsmap(map, longlat)
roms_path <- file.path(getOption("default.datadir"), "data_local/acecrc.org.au/ROMS/mertz_sample/mer_his_1992_01.nc")

set_crs <- function(x, value) {
  projection(x) <- value
  x
}
epsg3857 <- "+init=epsg:3857 +proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null
+no_defs"
set_ext <- function(x) setExtent(x, extent(c(xmin(x), xmax(x), ymin(x), ymax(x)) * 10))
ll <- leaflet()
itime <- 0
vname <- "u"
ui <- fluidPage(
  leafletOutput("mymap"),
  p(),
  actionButton("increment_time", "Increment slice")
)

server <- function(input, output, session) {
  
  get_time <- eventReactive(input$increment_time, {
    itime <<- itime %% 31 + 1
    r0 <- romsdata(roms_path, varname = vname, slice = c(1, itime), transpose = TRUE)
    r1 <- romsdata(roms_path, varname = "v", slice = c(1, itime), transpose = TRUE)
    r <- sqrt(r0 * r0 + r1 * r1)
    set_ext(set_crs(r, epsg3857))
  }, ignoreNULL = FALSE)

  output$mymap <- renderLeaflet({
    leaflet() %>%
      #addMarkers(data = points())
    
      addRasterImage(get_time(), project = FALSE) ##%>% addPolygons(data = map, project = FALSE)
  })
}

shinyApp(ui, server)