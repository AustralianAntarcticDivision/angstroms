#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(angstroms)
library(raster)
library(rbgm)
library(bgmfiles)
romsfiles <- structure(list(fullname = c("/rdsi/PRIVATE/data_local/acecrc.org.au/ROMS/s_corney/cpolar/ocean_his_3101.nc", 
                                         "/rdsi/PRIVATE/data_local/acecrc.org.au/ROMS/s_corney/cpolar/ocean_his_3102.nc", 
                                         "/rdsi/PRIVATE/data_local/acecrc.org.au/ROMS/s_corney/cpolar/ocean_his_3103.nc", 
                                         "/rdsi/PRIVATE/data_local/acecrc.org.au/ROMS/s_corney/cpolar/ocean_his_3104.nc", 
                                         "/rdsi/PRIVATE/data_local/acecrc.org.au/ROMS/s_corney/cpolar/ocean_his_3105.nc", 
                                         "/rdsi/PRIVATE/data_local/acecrc.org.au/ROMS/s_corney/cpolar/ocean_his_3106.nc", 
                                         "/rdsi/PRIVATE/data_local/acecrc.org.au/ROMS/s_corney/cpolar/ocean_his_3107.nc", 
                                         "/rdsi/PRIVATE/data_local/acecrc.org.au/ROMS/s_corney/cpolar/ocean_his_3108.nc", 
                                         "/rdsi/PRIVATE/data_local/acecrc.org.au/ROMS/s_corney/cpolar/ocean_his_3109.nc", 
                                         "/rdsi/PRIVATE/data_local/acecrc.org.au/ROMS/s_corney/cpolar/ocean_his_3110.nc", 
                                         "/rdsi/PRIVATE/data_local/acecrc.org.au/ROMS/s_corney/cpolar/ocean_his_3111.nc", 
                                         "/rdsi/PRIVATE/data_local/acecrc.org.au/ROMS/s_corney/cpolar/ocean_his_3112.nc", 
                                         "/rdsi/PRIVATE/data_local/acecrc.org.au/ROMS/s_corney/cpolar/ocean_his_3201.nc", 
                                         "/rdsi/PRIVATE/data_local/acecrc.org.au/ROMS/s_corney/cpolar/ocean_his_3202.nc", 
                                         "/rdsi/PRIVATE/data_local/acecrc.org.au/ROMS/s_corney/cpolar/ocean_his_3203.nc", 
                                         "/rdsi/PRIVATE/data_local/acecrc.org.au/ROMS/s_corney/cpolar/ocean_his_3204.nc", 
                                         "/rdsi/PRIVATE/data_local/acecrc.org.au/ROMS/s_corney/cpolar/ocean_his_3205.nc", 
                                         "/rdsi/PRIVATE/data_local/acecrc.org.au/ROMS/s_corney/cpolar/ocean_his_3206.nc", 
                                         "/rdsi/PRIVATE/data_local/acecrc.org.au/ROMS/s_corney/cpolar/ocean_his_3207.nc", 
                                         "/rdsi/PRIVATE/data_local/acecrc.org.au/ROMS/s_corney/cpolar/ocean_his_3208.nc", 
                                         "/rdsi/PRIVATE/data_local/acecrc.org.au/ROMS/s_corney/cpolar/ocean_his_3209.nc", 
                                         "/rdsi/PRIVATE/data_local/acecrc.org.au/ROMS/s_corney/cpolar/ocean_his_3210.nc", 
                                         "/rdsi/PRIVATE/data_local/acecrc.org.au/ROMS/s_corney/cpolar/ocean_his_3211.nc", 
                                         "/rdsi/PRIVATE/data_local/acecrc.org.au/ROMS/s_corney/cpolar/ocean_his_3212.nc", 
                                         "/rdsi/PRIVATE/data_local/acecrc.org.au/ROMS/s_corney/cpolar/ocean_his_3301.nc", 
                                         "/rdsi/PRIVATE/data_local/acecrc.org.au/ROMS/s_corney/cpolar/ocean_his_3302.nc", 
                                         "/rdsi/PRIVATE/data_local/acecrc.org.au/ROMS/s_corney/cpolar/ocean_his_3303.nc", 
                                         "/rdsi/PRIVATE/data_local/acecrc.org.au/ROMS/s_corney/cpolar/ocean_his_3304.nc", 
                                         "/rdsi/PRIVATE/data_local/acecrc.org.au/ROMS/s_corney/cpolar/ocean_his_3305.nc", 
                                         "/rdsi/PRIVATE/data_local/acecrc.org.au/ROMS/s_corney/cpolar/ocean_his_3306.nc", 
                                         "/rdsi/PRIVATE/data_local/acecrc.org.au/ROMS/s_corney/cpolar/ocean_his_3307.nc", 
                                         "/rdsi/PRIVATE/data_local/acecrc.org.au/ROMS/s_corney/cpolar/ocean_his_3308.nc", 
                                         "/rdsi/PRIVATE/data_local/acecrc.org.au/ROMS/s_corney/cpolar/ocean_his_3309.nc", 
                                         "/rdsi/PRIVATE/data_local/acecrc.org.au/ROMS/s_corney/cpolar/ocean_his_3310.nc", 
                                         "/rdsi/PRIVATE/data_local/acecrc.org.au/ROMS/s_corney/cpolar/ocean_his_3311.nc", 
                                         "/rdsi/PRIVATE/data_local/acecrc.org.au/ROMS/s_corney/cpolar/ocean_his_3312.nc"
), date = structure(c(31496400, 34174800, 36594000, 39276000, 
                      41868000, 44546400, 47138400, 49816800, 52495200, 55087200, 57762000, 
                      60354000, 63032400, 65710800, 68220000, 70898400, 73490400, 76168800, 
                      78760800, 81439200, 84117600, 86709600, 89384400, 91976400, 94654800, 
                      97333200, 99752400, 102434400, 105026400, 107704800, 110296800, 
                      112975200, 115653600, 118245600, 120920400, 123512400), class = c("POSIXct", 
                                                                                        "POSIXt"), tzone = "")), class = "data.frame", .Names = c("fullname", 
                                                                                                                                                  "date"), 
row.names = c(NA, -36L))

library(dplyr)



file_db <- bind_rows(lapply(romsfiles$fullname, function(x) {
  nc <- ncdump::NetCDF(x)
  
  tlen <- filter(nc$dimension, name == "ocean_time")$len
  tibble(fullname = rep(x, tlen), band_level = seq_len(tlen))
}))
file_db$date <- seq(min(romsfiles$date), length = nrow(file_db), by = "1 day")
nc <- ncdump::NetCDF(romsfiles$fullname[1])
var4 <- nc$variable %>% filter(ndims == 4)
## get a BGM and read it
bfile <- bgmfiles::bgmfiles("antarctica_28")
bgm <- bgmfile(bfile)


priorityvar <- c("temp", "salt", "u", "v", "sqrt(u * u  + v * v)")
crazyvar <- c(priorityvar, setdiff(var4$name, priorityvar))

coords <- romscoords(romsfiles$fullname[1])
boxmap <- romsmap(boxSpatial(bgm), coords) 
facemap <- romsmap(faceSpatial(bgm), coords) 

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("ROMS and Atlantis"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      dateInput(inputId = "date", label = "Day of year", value = "1971-01-01", min = "1971-01-01", max = "1973-12-31"), 
      numericInput(inputId = "roms_level", label = "ROMS Level (31 is surface, 1 is sea floor)", value = 31, min = 1, max = 31, step = 1),
      selectInput(inputId = "box_id", label = "Atlantis box", selected = "<all>", choices = c("<all>", boxmap$label)), 
      selectInput(inputId = "variable", label = "Variable", choices = crazyvar ), 
      checkboxInput(inputId = "box_label", label = "Box Label", value = TRUE), 
      checkboxInput(inputId =  "box_zoom", label = "Zoom to Atlantis", value = TRUE), 
      checkboxInput(inputId = "extract", label = "Extract summary", value = FALSE), 
      checkboxInput(inputId = "mask", label = "Mask to Atlantis", value = FALSE)
      
      #selectInput(inputId = "atlantis", label = "Atlantis model", choices = c("antarctica_28", "antarctica_99", "VMPA_setas.bgm"))
      
      #,shiny::selectInput(inputId = "graticule", label = "Graticule", choices = c("graticule", "none"))
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("mapPlot")
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  output$mapPlot <- renderPlot({
    idx <- findInterval(as.Date(input$date), as.Date(file_db$date))
    titl <- sprintf("%s level: %i ", basename(file_db$fullname[idx]), input$roms_level)
 
    
    if (input$box_id == "<all>") {
      boxplot <- boxmap
      gcord <- coords# crop(coords, boxmap)
    } else {
      boxplot <- subset(boxmap, input$box_id == label)
      gcord <- coords # crop(coords, boxplot)
    }
    if (input$variable == "sqrt(u * u  + v * v)") {
      if (input$box_zoom) {
        roms <- crop(romsdata(file_db$fullname[idx], varname = "u", transpose = TRUE, slice = c(input$roms_level, file_db$band_level)), extent(boxplot) + 1)
        roms2 <- crop(romsdata(file_db$fullname[idx], varname = "v", transpose = TRUE, slice = c(input$roms_level, file_db$band_level)), extent(boxplot) + 1)
      } else {
        roms <- romsdata(file_db$fullname[idx], varname = "u", transpose = TRUE, slice = c(input$roms_level, file_db$band_level))
        roms2 <- romsdata(file_db$fullname[idx], varname = "v", transpose = TRUE, slice = c(input$roms_level, file_db$band_level))
        
      }
      roms <- sqrt(roms * roms + roms2 * roms2)
      
    } else {
      if (input$box_zoom) {
      # plot(1, main = input$date)
       roms <- crop(romsdata(file_db$fullname[idx], varname = input$variable, transpose = TRUE, slice = c(input$roms_level, file_db$band_level)), extent(boxplot) + 1)
      }  else {
        roms <- romsdata(file_db$fullname[idx], varname = input$variable, transpose = TRUE, slice = c(input$roms_level, file_db$band_level))
      }
    }
    if (input$mask) {
      roms <- mask(roms, boxplot)
    }
    if (input$extract) {
      par(mfrow = c(1, 2))
    } else {
      par(mfrow = c(1, 1))
    }
    plot(roms, col = viridis::viridis(100), asp = NA, addfun = function() {
      plot(boxplot, add = TRUE)
      if (input$box_label) text(boxplot, lab = boxplot$label, cex = 1.2, col = "hotpink")
      
      title(titl)
   
    })
    
    if (input$extract) {
      vals <- unlist(extract(roms, boxplot))
      plot(density(vals, na.rm = TRUE), main = sprintf("%s mean=%f", input$variable, mean(vals, na.rm = TRUE)))
      
    }
  })


}

# Run the application 
shinyApp(ui = ui, server = server)

