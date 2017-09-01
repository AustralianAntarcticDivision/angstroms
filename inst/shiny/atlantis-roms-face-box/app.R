
library(shiny)
library(angstroms)
library(raster)
library(rbgm)
library(bgmfiles)
library(tibble)
romsfiles <- tibble(fullname = "/mnt/temproms/ocean_his_31month_av.nc", 
                    date = as.POSIXct("2009-01-01", tz = "GMT"))
library(dplyr)
file_db <- bind_rows(lapply(romsfiles$fullname, function(x) {
  nc <- ncdump::NetCDF(x)
  tlen <- filter(nc$dimension, name == "ocean_time")$len
  tibble(fullname = rep(x, tlen), band_level = seq_len(tlen))
}))
file_db$date <- seq(min(romsfiles$date), length = nrow(file_db), by = "1 month")

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
      dateInput(inputId = "date", label = "Day of year", value = "2009-01-01", min = "2009-01-01", max = "2009-12-01"), 
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
      subface <- bind_rows(boxplot@data[, ".bx0"] %>% inner_join(facemap@data, c(".bx0" = "left")), 
                           boxplot@data[, ".bx0"] %>% inner_join(facemap@data, c(".bx0" = "right"))) %>% 
        select(.fx0) %>% distinct()
      faceplot <- facemap[match(subface$.fx0, facemap$.fx0), ]
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
      #par(mfrow = c(1, 2))
      op <- par(mfrow = c(1, 2))
    } else {
      op <- par(mfrow = c(1, 1))
    }
    plot(roms, col = viridis::viridis(100), asp = NA, addfun = function() {
      plot(boxplot, add = TRUE)
      if (input$box_label) text(boxplot, lab = boxplot$label, cex = 1.2, col = "hotpink")
      if (!input$box_id == "<all>") {
        text(coordinates(rgeos::gCentroid(faceplot, byid = TRUE)), label = gsub("ace", "", faceplot$label), col = "black", xpd = NA)
      }
      title(titl)
    
    })
    
    if (input$extract) {
      vals <- unlist(extract(roms, boxplot))
      plot(density(vals, na.rm = TRUE), main = sprintf("%s mean=%f", input$variable, mean(vals, na.rm = TRUE)))
      #facevals <- extract(roms, faceplot)
    } 
      par(op)
    
    
  })


}

# Run the application 
shinyApp(ui = ui, server = server)

