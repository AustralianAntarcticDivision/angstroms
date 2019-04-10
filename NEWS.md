# angstroms dev

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



