## which box does each point fall in
index_box <- function(box_sp, roms_ll) {
  ind <- sp::over(project_to(coords_points(roms_ll), box_sp) , as(box_sp, "SpatialPolygons"))
  #  tibble(box = ind,
  tibble(box = box_sp$label[ind], 
         cell = seq_len(ncell(roms_ll))) %>% 
    filter(!is.na(box))
}

## 
index_face <- function(face_sp, roms_ll) {
  roms_face <- romsmap(project_to(face_sp, "+init=epsg:4326"), roms_ll)
  ind_face <- cellFromLine(roms_ll, roms_face)
  tibble(face = roms_face$label[rep(seq_len(nrow(roms_face)), lengths(ind_face))], 
         cell = unlist(ind_face))
}
