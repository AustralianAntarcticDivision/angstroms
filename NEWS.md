# angstroms dev

* Fix closing files, see #23. 

# angstroms 0.0.2.9005

* Now import reproj for coordinate reprojection. 

* Remove dependency on spbabel and rgdal. 

* Allow `rawdata()` to maintain original state of data values (`native` argument). 

* Add rectilinear support to `romscoords()`, by expanding the coordinates as for the curvilinear case.  This function now tries hard to find the coordinates, assuming they are 1D arrays if the 2D assumption fails. 

* Argument `transpose` now TRUE by default for `romsdata()` family and `romscoords()`.

* Argument `varname` is now the empty string by default to match raster behaviour, and benefit 
from messages abuot available variables. 

* Function `romscoords()` gains two new arguments `flip_y` to control latitude orientation, and `varname` in case the default raster doesn't match the coordinate requested in the rectlinear case. 

* Replace nabor import with FNN. 

* Remove rgeos dependency (for use on NCI). 

* Added 'create_transport' function to flesh out details of a hydrodynamic forcing file for Atlantis. 

* add better romshcoords logic, and romsdepth alias (simple arg provides original, 
probably wrong, implementation)

* add lvar argument to romsdata3d (needs more thought)

* update missing uses of `ncdf = TRUE` thanks to Ryan Morse in https://github.com/hypertidy/angstroms/issues/6

* added `databoundary` function, to trace around valid pixels

* more general approach to reading slices, with `romsdata3d` and `romsdata2d`

* S3 methods for rawdata

# angstroms 0.0.1

* added documentation for croproms and ROMS URL, thanks to CRAN feedback

* romsdata defaults to transpose

* updates and checks


# angstroms 0.0.0.9000

* use spbabel::sp not spFromTable

* Added ncraster function to read arbitrary slices

* Added romshcoords function to flesh out "h" by "Cs_r"

* Added a `NEWS.md` file to track changes to the package.



