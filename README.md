
[![Travis-CI Build Status](https://travis-ci.org/mdsumner/angstroms.svg?branch=master)](https://travis-ci.org/mdsumner/angstroms) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/mdsumner/angstroms?branch=master&svg=true)](https://ci.appveyor.com/project/mdsumner/angstroms) [![Coverage Status](https://img.shields.io/codecov/c/github/mdsumner/angstroms/master.svg)](https://codecov.io/github/mdsumner/angstroms?branch=master) <!-- README.md is generated from README.Rmd. Please edit that file -->

angstroms
=========

The goal of angstroms is to provide easy access to Regional Ocean Modeling System (ROMS) output for R.

Installation
------------

You can install the development version of angstroms from github with:

``` r
# install.packages("devtools")
devtools::install_github("mdsumner/angstroms")
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

Please note that this project is released with a [Contributor Code of Conduct](CONDUCT.md). By participating in this project you agree to abide by its terms.
