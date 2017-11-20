#' Coordinates at depth
#'
#' Extract the multi-layer 'h'eight grid with S-coordinate stretching applied
#'
#' Compute ROMS grid depth from vertical stretched variables
#' Given a bathymetry (h), free-surface (zeta) and terrain-following parameters, this function computes the 3D depths for the requested C-grid location. If the free-surface is not provided, a zero value is assumed resulting in unperturb depths.  This function can be used when generating initial conditions or climatology data for an application. Check the following link for details: https://www.myroms.org/wiki/index.php/Vertical_S-coordinate
#' See https://github.com/dcherian/tools/blob/master/ROMS/arango/utility/set_depth.m
#' Original Matlab code by Hernan Arango.
#' @param x ROMS file name
#' @param grid_type string: "rho","psi","u","v","w"
#' @param slice integer: if non-missing, use this time slice to index into zeta (free-surface). Otherwise assume zeta is zero (and hence depth is time-independent)
#' @param ... dots
#' @param depth string: the name of the appropriate variable to use for depth
#' @param simple logical: if TRUE, use the old "simple" depth method, which may not be correct
#' @return RasterStack with a layer for every depth
#' @export
romshcoords <- function(x, grid_type = "rho", slice, ..., depth = "h", simple = FALSE){
  grid_type <- match.arg(tolower(grid_type),c("rho","psi","u","v","w"))
  S <- if (grid_type=="w") "Cs_w" else "Cs_r"
  h <- romsdata(x, varname = depth)
  Cs_r <- ncget(x, S)
  v <- values(h)
  if (simple) {
    ## simplistic, early version - probably should be defunct
    out <- set_indextent(brick(array(rep(rev(Cs_r), each = length(v)) * v,
                                     c(ncol(h), nrow(h), length(Cs_r))), transpose = TRUE))
  } else {

    Vtransform <- as.integer(ncget(x,"Vtransform"))
    if (!Vtransform %in% c(1,2)) stop("Vtransform must be 1 or 2")

    hc <- ncget(x,"hc")

    depth_grid <- if (grid_type=="w") "w" else "rho"

    zeta <- if (missing(slice)) 0 else stop("not coded yet")##angstroms::romsdata2d(x,"zeta",slice=slice,transpose=FALSE)
    N <- length(ncget(x,"Cs_r"))
    Np <- N+1

    h <- ncget(x,"h")
    hmin <- min(h)
    hmax <- max(h)

    Lp <- dim(h)[1]
    Mp <- dim(h)[2]
    L <- Lp-1
    M <- Mp-1

    z <- array(NA,dim=c(Lp,Mp,if (grid_type=="w") Np else N))

    ## Compute vertical stretching function, C(k):
    ##stretch <- stretching(x,depth_grid)
    if (depth_grid=="w") {
      stretch <- list(C=ncget(x,"Cs_w"),s=ncget(x,"s_w"))
    } else {
      stretch <- list(C=ncget(x,"Cs_r"),s=ncget(x,"s_rho"))
    }

    ## Average bathymetry and free-surface at requested C-grid type.
    if (grid_type=="rho") {
      hr <- h
      zetar <- zeta
    } else if (grid_type=="psi") {
      hp <- 0.25*(h[1:L,1:M]+h[2:Lp,1:M]+h[1:L,2:Mp]+h[2:Lp,2:Mp])
      zetap <- 0.25*(zeta[1:L,1:M]+zeta[2:Lp,1:M]+zeta[1:L,2:Mp]+zeta[2:Lp,2:Mp])
    } else if (grid_type=="u") {
      hu <- 0.5*(h[1:L,1:Mp]+h[2:Lp,1:Mp])
      zetau <- 0.5*(zeta[1:L,1:Mp]+zeta[2:Lp,1:Mp])
    } else if (grid_type=="v") {
      hv <- 0.5*(h[1:Lp,1:M]+h[1:Lp,2:Mp])
      zetav <- 0.5*(zeta[1:Lp,1:M]+zeta[1:Lp,2:Mp])
    } else if (grid_type=="w") {
      hr <- h
      zetar <- zeta
    } else {
      stop("unsupported grid_type: ",grid_type)
    }

    ## Compute depths (m) at requested C-grid location.

    if (Vtransform == 1) {
      if (grid_type=="rho") {
        for (k in seq_len(N)) {
          z0 <- (stretch$s[k]-stretch$C[k])*hc + stretch$C[k]*hr
          z[,,k] <- z0 + zetar*(1.0 + z0/hr)
        }
      } else if (grid_type=="psi") {
        for (k in seq_len(N)) {
          z0 <- (stretch$s[k]-stretch$C[k])*hc + stretch$C[k]*hp
          z[,,k] <- z0 + zetap*(1.0 + z0/hp)
        }
      } else if (grid_type=="u") {
        for (k in seq_len(N)) {
          z0 <- (stretch$s[k]-stretch$C[k])*hc + stretch$C[k]*hu
          z[,,k] <- z0 + zetau*(1.0 + z0/hu)
        }
      } else if (grid_type=="v") {
        for (k in seq_len(N)) {
          z0 <- (stretch$s[k]-stretch$C[k])*hc + stretch$C[k]*hv;
          z[,,k] <- z0 + zetav*(1.0 + z0/hv)
        }
      } else if (grid_type=="w") {
        z[,,1] <- -hr
        for (k in seq(from=2,to=Np,by=1)) {
          z0 <- (stretch$s[k]-stretch$C[k])*hc + stretch$C[k]*hr
          z[,,k] <- z0 + zetar*(1.0 + z0/hr)
        }
      } else {
        stop("unsupported grid_type: ",grid_type)
      }
    } else if (Vtransform == 2) {
      if (grid_type=="rho") {
        for (k in seq_len(N)) {
          z0 <- (hc*stretch$s[k]+stretch$C[k]*hr)/(hc+hr)
          z[,,k] <- zetar+(zeta+hr)*z0
        }
      } else if (grid_type=="psi") {
        for (k in seq_len(N)) {
          z0 <- (hc*stretch$s[k]+stretch$C[k]*hp)/(hc+hp)
          z[,,k] <- zetap+(zetap+hp)*z0
        }
      } else if (grid_type=="u") {
        for (k in seq_len(N)) {
          z0 <- (hc*stretch$s[k]+stretch$C[k]*hu)/(hc+hu)
          z[,,k] <- zetau+(zetau+hu)*z0
        }
      } else if (grid_type=="v") {
        for (k in seq_len(N)) {
          z0 <- (hc*stretch$s[k]+stretch$C[k]*hv)/(hc+hv)
          z[,,k] <- zetav+(zetav+hv)*z0
        }
      } else if (grid_type=="w") {
        for (k in seq_len(Np)) {
          z0 <- (hc*stretch$s[k]+stretch$C[k]*hr)/(hc+hr)
          z[,,k] <- zetar+(zetar+hr)*z0
        }
      } else {
        stop("unsupported grid_type: ",grid_type)
      }
    } else {
      stop("Vtransform must be 1 or 2")
    }
    ## FIXME all these flips and twirls can be applied more efficiently (or avoided)
    ## though should layers start at the surface and go down or ...
    out <- raster::flip(set_indextent(raster::brick(z, transpose = TRUE)), "y")
    out <- raster::subset(out, rev(seq_len(raster::nlayers(out))))

  }

 out
}

# FIXME: should be romsdepth the name anyway ...
#' @name romshcoords
#' @export
romsdepth <- function(x, ...) {...
 romshcoords(x, ...)
}


