
[![Travis-CI Build Status](https://travis-ci.org/hypertidy/angstroms.svg?branch=master)](https://travis-ci.org/hypertidy/angstroms) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/hypertidy/angstroms?branch=master&svg=true)](https://ci.appveyor.com/project/hypertidy/angstroms) [![Coverage Status](https://img.shields.io/codecov/c/github/hypertidy/angstroms/master.svg)](https://codecov.io/github/hypertidy/angstroms?branch=master) <!-- README.md is generated from README.Rmd. Please edit that file -->

angstroms
=========

The goal of angstroms is to provide easy access to Regional Ocean Modeling System (ROMS) output for R.

Installation
------------

You can install the development version of angstroms from github with:

``` r
# install.packages("devtools")
devtools::install_github("hypertidy/angstroms")
```

angstroms - R for ROMS
----------------------

Angstroms aims to make working with ROMS output as easy as possible in R. Rather than re-map explicitly the complex curvilinear grid in ROMS, the approach simplifies this by:

-   maintaining the internal index of ROMS as the default *georeferencing*
-   converting external data (maps, transects, points, etc.) into the native internal index space of ROMS
-   providing tools to read arbitrary slices from the grids (either 2D or 3D) as Raster objects
-   providing tools to recover the original full coordinates as needed

In combination these allow extraction and query from the ROMS output very easily.

The ability to deal with time series across multiple files is still in development, though can be used simply now with standard loops.

Some more examples:

<http://rpubs.com/cyclemumner/roms0>

<http://rpubs.com/cyclemumner/266770>

Contributing
------------

Here are some notes on lessons learned, all is work in progress.

-   calls to raster/brick/stack should be wrapped by an appropriate `romsdata()`, to ensure that `ncdf = TRUE` is consistent, and to reduce the amount of extent-setting and orienting code
-   romsdata() should have a convention to accept a raw matrix/array, see first point
-   use `set_indextent` on a brick/raster to ensure it is in index-extent, again to reduce code duplications
-   `ncdf=TRUE` is required, because on some systems we might fall back to rgdal which won't support NetCDF properly (can be done with SDS strings or VRT but that's a hassle)
-   use of raster (and so ncdf4) package could be better done more directly with RNetCDF / tidync
-   RNetCDF loses speed on ncdf4 for some metadata harvesting, but is otherwise a better choice right now (Nov 2017) especially with Thredds support for Win64 now on CRAN
-   orientation is an open question, but at the moment it looks like X-increasing, Y-decreasing, Z-decreasing, T-increasing - though it's not that simple - and for no good reason other than train-wreck-convention
-   romsdata() generally needs a bit more thought, should it take a raw array or only a source string, do 2d/3d variants make sense, how to orient specific slice conventions (roms\_xt etc.)

Conduct
-------

Please note that this project is released with a [Contributor Code of Conduct](CONDUCT.md). By participating in this project you agree to abide by its terms.
