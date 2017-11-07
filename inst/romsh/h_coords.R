# There are small differences (up to about 5m) which seem a little large to be
# just numerical/rounding errors. I haven't looked in detail yet as to what is
# going on, but the equation they use [along the lines of
# (hc*s_rho+Cs_r*h)/(hc+h)] does seem to be different to what you are using [
# just Cs_r*h  (?)]. I also suspect yours won't work if the Vtransform method is
# not 2, but that's a bit of a guess.


x <- path to roms file

h1 <- romsdepth(x,"rho") ## my function, just returns an array
h1 <- aperm(h1,c(2,1,3)) ## reorder dims and directions to match raster's
h1 <- h1[nrow(h1):1,,dim(h1)[3]:1]

h2 <- romshcoords(x) ## your version

hlev <- 10
h1r <- h2[[1]]
values(h1r) <- h1[,,hlev]
plot(h1r-h2[[hlev]])