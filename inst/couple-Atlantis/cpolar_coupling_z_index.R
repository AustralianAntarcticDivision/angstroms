## ROMS Atlantis coupling MDSumner, Jessica Melbourne-Thomas April 2017

## NOTES

#https://confluence.csiro.au/display/Atlantis/Hydro+FAQ

#https://github.com/cecilieha/NoBA/blob/master/create_flux_ncdf_atlantis.R

##corners <- read.delim("https://raw.githubusercontent.com/cecilieha/NoBA/master/corners_neighbours_nordic.txt")

## dummy hydro
# sohydro <- "~/projects/AtlantisEastAntartica/SOhydrodummy.nc"
# library(ncdump)
# nc <- NetCDF(sohydro)

#' Convenience function to transform map projection . 
#'
#' Transform `x` to whatever the projection of `to` is. 
#' @param x  object to transform
#' @param to object with a map projection
#'
#' @return `x`, transformed
#' @export
#' @importFrom raster projection
#' @importFrom sp spTransform
project_to <- function(x, to) {
  spTransform(x, CRS(projection(to)))
}
## devtools::install_github(c("mdsumner/angstroms"))
library(angstroms)
library(rbgm) ## read BGM
library(bgmfiles) ## archive of BGM files
library(raadtools)
library(ncdump)
## get the Circumpolar files
cpolar <- raadtools:::cpolarfiles()
roms_file <- cpolar$fullname[2]

## get a BGM and read it
bfile <- bgmfiles::bgmfiles("antarctica_28")
bgm <- bgmfile(bfile)


## get the longitude/latitude arrays
roms_ll <- romscoords(roms_file, transpose = TRUE)

## we need the unsullied boxes to identify points inside them
boxes <- boxSpatial(bgm)
## which box does each point fall in
index_box <- function(box_sp, roms_ll) {
  ind <- sp::over(project_to(coords_points(roms_ll), box_sp) , as(box_sp, "SpatialPolygons"))
  tibble(box = ind, 
         cell = seq_len(ncell(roms_ll))) %>% 
    filter(!is.na(box))
}


## build the index for each box to the ROMS cells it contains
## and each face for the ROMS cells it traverses
library(dplyr)
box_roms_index <- index_box(boxes, roms_ll)
  
#' return the ramp of positive depths from surface down 
#' (so that the order is native to the NetCDF order)
roms_level <- function(Cs_r, h, cell) {
  extract(h, cell) *  -rev(Cs_r)
}
## put any raster into xy-index space (0, nc, 0, nr)
set_indextent <- function(x) {
  setExtent(x, extent(0, ncol(x), 0, nrow(x)))
}
h <- raster(roms_file, varname = "h")
## Cs_r is the S-coord stretching
Cs_r <- rawdata(roms_file, "Cs_r")


## build the level index between Atlantis and ROMS
box_roms_index$botz <- boxes$botz[box_roms_index$box]
list_nc_z_index <- vector('list', nrow(box_roms_index))
for (i in seq_len(nrow(box_roms_index))) {
  rl <- roms_level(Cs_r, h, box_roms_index$cell[i])
  ## implicit 0 at the surface, and implicit bottom based on ROMS
  list_nc_z_index[[i]] <- findInterval(rl, cumsum(c(0, rev(rbgm::build_dz(box_roms_index$botz[i])))))
if (i %% 1000 == 0) print(i)
  }

## join the box-xy-index to the level index
z_index <- bind_rows(lapply(list_nc_z_index, 
                 function(x) tibble(atlantis_level = x, roms_level = seq_along(x))), 
          .id = "cell_index") %>% 
  inner_join(mutate(box_roms_index, cell_index = as.character(row_number()))) %>% 
  select(-cell_index, -botz)

#saveRDS(z_index, file = "z_index.rds")
## driver function to loop over x raster by levels
## matching $cell for the right group
extract_at_level <- function(x, cell_level) {
  ulevel <- unique(cell_level$level)
  values <- numeric(nrow(cell_level))
  for (ul in seq_along(ulevel)) {
    asub <- cell_level$level == ulevel[ul]
    values[asub] <- extract(x[[ulevel[ul]]], 
            cell_level$cell[asub])
  }
  values
}

ru <- set_indextent(brick(roms_file, varname = "u", lvar = 4, level = 1))
z_index$u <- extract_at_level(ru, rename(z_index, level = roms_level, cell = cell))

rv <- set_indextent(brick(roms_file, varname = "v", lvar = 4, level = 1))
z_index$v <- extract_at_level(rv, rename(z_index, level = roms_level, cell = cell))


plot(romsmap(boxes[15, ], roms_ll), add = F, border = "red", lwd = 4)
plot(sqrt(ru[[10]]^2 + rv[[10]]^2), add = TRUE, 
     col = colorRampPalette(c("#132B43", high = "#56B1F7"))(10), zlim = c(0, 0.1))
plot(romsmap(boxes[15, ], roms_ll), add = TRUE, border = "red", lwd = 4)
ggplot(mutate(z_index %>% 
                filter(box == 15), 
              x = xFromCell(ru, cell), y = yFromCell(ru, cell)), 
       aes(x = x, y = y, fill = sqrt(u * u + v * v))) + geom_tile() + coord_equal() + 
  facet_wrap(~roms_level)


rm(ru, rv)


rdata <- set_indextent(brick(roms_file, varname = "temp", lvar = 4, level = 1))
z_index$temp <- extract_at_level(rdata, rename(z_index, level = roms_level, cell = cell))

rdata <- set_indextent(brick(roms_file, varname = "salt", lvar = 4, level = 1))
z_index$salt <- extract_at_level(rdata, rename(z_index, level = roms_level, cell = cell))
ggplot(z_index %>% mutate(x = xFromCell(rdata, cell), 
                          y = yFromCell(rdata, cell)), 
       aes(x = x, y = y, fill = temp)) + geom_tile() + 
  coord_equal()

ggplot(mutate(z_index %>% 
                filter(roms_level %% 3 == 0), 
              x = xFromCell(ru, cell), y = yFromCell(ru, cell)), 
       aes(x = x, y = y, fill = temp)) + geom_tile() + coord_equal() + 
  facet_wrap(~roms_level)
