#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(dplyr)
library(angstroms) ## devtools::install_github("mdsumner/angstroms")
library(ncdump) ## devtools::install_github("r-gris/ncdump")
f <- "/rdsi/PRIVATE/data_local/acecrc.org.au/ROMS/waom/sample/waom2_MinDepth20m_rx10.4_grd.nc"
nc <- NetCDF(f)
coord <- romscoords(f)
## put any raster into xy-index space (0, nc, 0, nr)
set_indextent <- function(x) {
  setExtent(x, extent(0, ncol(x), 0, nrow(x)))
}

vcolours <- viridis::viridis(100)
vartable <- nc$variable %>% filter(ndims > 1)
ui <- fluidPage(
   
   # Application title
   titlePanel("WAOM"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
         selectInput(inputId = "variable_name", label = "Variable", choices = vartable$name), 
         checkboxInput(inputId = "ice_contour", labe = "Ice contour", value = FALSE)
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("rasterPlot")
      )
   )
)

server <- function(input, output) {
   output$rasterPlot <- renderPlot({
      roms <- set_indextent(raster(f, varname = input$variable_name, ncdf = TRUE))
      afun <- NULL
      if (input$ice_contour) {
        ice <- raadtools::readice(latest = TRUE)
        dt <- getZ(ice)
        ice <- romsmap(rasterToContour(ice[[1]], levels = c(15, 50, 80)), coord)
        afun <- function() {
          plot(ice, add = TRUE)
          title(dt)
        }
      }
      plot(roms, col = vcolours, addfun = afun)
   })
}

# Run the application 
shinyApp(ui = ui, server = server, options = c(display.mode = "showcase"))

