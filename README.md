
<!-- badges: start -->

[![R build
status](https://github.com/AustralianAntarcticDivision/angstroms/workflows/R-CMD-check/badge.svg)](https://github.com/AustralianAntarcticDivision/angstroms/actions)
[![CRAN
status](http://www.r-pkg.org/badges/version/angstroms)](https://cran.r-project.org/package=angstroms)
![cranlogs](http://cranlogs.r-pkg.org./badges/angstroms)
![test-coverage](https://github.com/AustralianAntarcticDivision/angstroms/workflows/test-coverage/badge.svg)
![pkgdown](https://github.com/AustralianAntarcticDivision/angstroms/workflows/pkgdown/badge.svg)
[![R-CMD-check](https://github.com/AustralianAntarcticDivision/angstroms/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/AustralianAntarcticDivision/angstroms/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

<!-- README.md is generated from README.Rmd. Please edit that file -->

# angstroms

The goal of angstroms *was* ~~to provide easy access to Regional Ocean
Modeling System (ROMS)~~ but now it provides easy access and control
over *any* complex gridded output available to the
[raster](https://CRAN.R-project.org/package=raster) package in R.

The key features are:

- treat gridded data as existing index space
- handle coordinate arrays independently
- avoid remodelling data whenever possible
- re-map data into the grid space for extractions
- treat the grid as a mesh (which it always is).

This approach aligns with the plotting and interpretation tools in the
[quadmesh](https://CRAN.R-project.org/package=quadmesh) package,
particularly the
[mesh_plot()](https://hypertidy.github.io/quadmesh/reference/mesh_plot.html)
and
[quadmesh()](https://hypertidy.github.io/quadmesh/reference/quadmesh.html)
functions.

## Installation

Install the released version from [CRAN](https://CRAN.R-project.org/)
with:

``` r
install.packages("angstroms")
```

You can install the development version of angstroms from github with:

``` r
# install.packages("remotes")
remotes::install_github("AustralianAntarcticDivision/angstroms")
```

## angstroms - R for gridded model data (including ROMS)

Angstroms aims to make working with gridded output as easy as possible
in R. Rather than re-map explicitly the complex curvilinear grid in ROMS
or ACCESS (or many others), the approach simplifies this by:

- maintaining the internal index of a grid as the default
  *georeferencing*
- converting external data (maps, transects, points, etc.) into the
  native internal index space of ROMS
- providing tools to read arbitrary slices from the grids (either 2D or
  3D) as Raster objects
- providing tools to recover the original full coordinates as needed

In combination these allow extraction and query from the complex grids
output very easily.

The ability to deal with time series across multiple files is still in
development, though can be used simply now with standard loops.

## Contributing

Here are some notes on lessons learned, all is work in progress.

- calls to raster/brick/stack should be wrapped by an appropriate
  `romsdata()`, to ensure that `ncdf = TRUE` is consistent, and to
  reduce the amount of extent-setting and orienting code
- romsdata() should have a convention to accept a raw matrix/array, see
  first point
- use `set_indextent` on a brick/raster to ensure it is in index-extent,
  again to reduce code duplications
- `ncdf=TRUE` is required, because on some systems we might fall back to
  rgdal which won’t support NetCDF properly (can be done with SDS
  strings or VRT but that’s a hassle)
- use of raster (and so ncdf4) package could be better done more
  directly with RNetCDF / tidync
- RNetCDF loses speed on ncdf4 for some metadata harvesting, but is
  otherwise a better choice right now (Nov 2017) especially with Thredds
  support for Win64 now on CRAN
- orientation is an open question, but at the moment it looks like
  X-increasing, Y-decreasing, Z-decreasing, T-increasing - though it’s
  not that simple - and for no good reason other than
  train-wreck-convention
- romsdata() generally needs a bit more thought, should it take a raw
  array or only a source string, do 2d/3d variants make sense, how to
  orient specific slice conventions (roms_xt etc.)

## Conduct

Please note that this project is released with a [Contributor Code of
Conduct](CONDUCT.md). By participating in this project you agree to
abide by its terms.
