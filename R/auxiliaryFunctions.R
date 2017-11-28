#' Calculate basis functions for new observations
#'
#' @param z elements where basis is going to be calculated. Assumes z
#'   has been scaled to [0, 1]
#' @param n_basis how many basis functions should be calculated
#' @param system Basis system to be used. Options are "cosine",
#'   "Fourier", "Haar", and "Daubechies"
#'
#' @return A matrix of dimension length(z) by n_basis with entries
#'   consisting of the first n_basis functions in the basis evaluated
#'   at the points z
calculateBasis <- function(z, n_basis, system = "Fourier") {
  system <- match.arg(system, c("cosine", "Fourier", "Haar", "Daubechies"))
  basis_function <- switch(system,
                           cosine = cosine_basis,
                           Fourier = fourier_basis,
                           Haar = haar_basis,
                           Daubechies = daubechies_basis)
  return(basis_function(z, n_basis))
}

#' Evaluates cosine basis for new observations
#' @inheritParams calculateBasis
#' @return A matrix of dimension length(z) by n_basis with entries
#'   consisting of the first n_basis cosine basis functions evaluated
#'   at the points z
cosine_basis <- function(z, n_basis) {
  basis <- matrix(NA, length(z), n_basis)
  basis[, 1] <- 1.0
  for (ii in 1:(n_basis - 1)) {
    basis[, ii + 1] <- sqrt(2) * cospi(ii * z)
  }
  attr(basis, "levels") <- seq_len(n_basis)
  return(basis)
}

#' Evaluates Fourier basis for new observations
#' @inheritParams calculateBasis
#' @return A matrix of dimension length(z) by n_basis with entries
#'   consisting of the first n_basis Fourier basis functions evaluated
#'   at the points z
fourier_basis <- function(z, n_basis) {
  basis <- matrix(NA, length(z), n_basis)
  basis[, 1] <- 1.0
  max_power <- (n_basis - 1) / 2
  for (ii in seq_len(max_power)) {
    basis[, 2 * ii] <- sqrt(2) * sinpi(2 * ii * z)
    basis[, 2 * ii + 1] <- sqrt(2) * cospi(2 * ii * z)
  }
  if (n_basis %% 2 == 0) {
    basis[, n_basis] <- sqrt(2) * sinpi(n_basis * z)
  }
  attr(basis, "levels") <- seq_len(n_basis)
  return(basis)
}


.pow2seq <- function(n) {
  seq <- rep(NA, n)
  pow <- 1
  for (ii in 1:n) {
    if (ii >= (2 ^ pow)) {
      pow <- pow + 1
    }
    seq[ii] <- pow
  }
  return(seq)
}


#' Evaluates Haar mother wavlet
#' @param x: float, value at which to evaluate the wavlet
#' @return float; the Haar mother wavlet evaluated at x
.haar_phi <- function(x) {
  if (0 <= x && x < 0.5) {
    return(1)
  } else if (0.5 <= x && x < 1) {
    return(-1)
  } else {
    return(0)
  }
}

#' Evaluates Haar basis for new observations
#' @inheritParams calculateBasis
#' @return A matrix of dimension length(z) by n_basis with entries
#'   consisting of the first n_basis Haar basis functions evaluated
#'   at the points z
haar_basis <- function(z, n_basis) {
  basis <- matrix(NA, length(z), n_basis)
  basis[, 1] <- 1.0

  kk <- 0
  jj <- 0

  for (ii in 2:n_basis) {
    if (jj == 2 ^ kk - 1) {
      kk <- kk + 1
      jj <- 0
    } else {
      jj <- jj + 1
    }
    basis[, ii] <- 2 ^ (kk / 2) * sapply(2 ^ kk * z - jj, .haar_phi)
  }

  attr(basis, "levels") <- .pow2seq(n_basis)
  return(basis)
}

#' Evaluates Daubechies basis for new observations
#' @inheritParams calculateBasis
#' @param filter_number: integer, controls smoothness
#' @param n_aux_basis: integer, the number of auxillary basis
#'   functions to compute
#' @param family: string, family of wavelet
#' @return A matrix of dimension length(z) by n_basis with entries
#'   consisting of the first n_basis Daubechies basis functions
#'   evaluated at the points z

daubechies_basis <- function(z, n_basis, n_aux_basis = max(n_basis, 2^12),
                             filter_number = 10, family = "DaubExPhase") {
  # Strategy is to project z onto a fine grid to use NN approximation
  aux <- wavethresh::GenW(n_aux_basis, filter.number = filter_number,
                          family = family)

  z_grid <- seq(0, 1, length.out = n_aux_basis)
  which_closest <- FNN::get.knnx(z_grid, as.matrix(z), k = 1)$nn.index

  aux[, 2:ncol(aux)] <- aux[, ncol(aux):2]
  aux <- aux[, 1:n_basis, drop = FALSE]
  basis <- aux[which_closest, ] / aux[1, 1]

  if (n_basis == 1) {
    basis <- as.matrix(basis)
  }

  attr(basis, "levels") <- .pow2seq(n_basis)
  return(basis)
}

#' Project onto density estimates
#' Uses binary search to determine the correct cutoff
normalize_density <- function(bin_size, estimates,
                              tolerance = 1e-3, max_iters = 500) {

  area <- bin_size * sum(pmax(estimates, 0.0))
  if (area < 1) {
    return(pmax(estimates / area, 0.0))
  }

  upper <- max(estimates)
  lower <- 0.0
  middle <- (upper + lower) / 2

  iter <- 1
  while (iter <= max_iters) {
    iter <- iter + 1

    density <- pmax(estimates - middle, 0.0)
    area <- bin_size * sum(density)

    if (abs(area - 1) < tolerance) {
      break
    }

    if (area > 1) {
      lower <- middle
    } else {
      upper <- middle
    }
    middle <- (upper + lower) / 2
  }

  return(pmax(estimates - middle, 0.0))
}

#' Removing bumps from density estimates
#'
#' @param binSize size of grid bins
#' @param estiamtes vector of density estimates on a grid
#' @param delta area threshold for removing a bump
#'
#' @return A vector of the density with bumps removed
remove_bumps <- function(binSize, estimates, delta) {
  runs <- rle(estimates > 0)
  n_runs <- length(runs$values)

  if (n_runs == 1) {
    return(estimates)
  }

  lower <- c(1, 1 + cumsum(runs$lengths))
  upper <- cumsum(runs$lengths)

  for (ii in 1:n_runs) {
    if (!runs$values[ii]) {
      next
    }
    area <- binSize * sum(estimates[lower[ii]:upper[ii]])
    if (area < delta) {
      estimates[lower[ii]:upper[ii]] <- 0.0
    }
  }

  return(estimates)
}

#' Sharpen estimated of a conditional density estimator
#'
#' @param binSize size of grid bins
#' @param estimates vector of density estimates on a grid
#' @param alpha sharpen coefficient
#'
#' @return A vector of the sharp density
sharpen <- function(binSize, estimates, alpha) {
  new_estimates=(estimates)^alpha
  # reescale so that it is a proper density
  new_estimates=(new_estimates/sum(new_estimates))/binSize

  return(new_estimates)
}


#' Post process density estimate to normalize and remove bumps
#'
#' @param binSize size of grid bins
#' @param estimates vector of density estimate on a grid
#' @param delta parameter for removing bumps; default is 0.0 implying
#'   no bump removal
#' @param threshold parameter for discarding small densities; default
#'   is 1e-6
#'
#' @return A vector of the density estimate
post_process <- function(binSize, estimates, delta = 0.0,
                         threshold = 1e-6,
                         alpha=1) {
  if (!is.vector(estimates)) {
    estimates <- as.vector(estimates)
  }

  estimates[estimates < threshold] <- 0.0

  if (all(estimates == 0)) {
    warning("Density has only negative values; replacing with uniform density")
    estimates[] <- 1.0
  }

  estimates <- normalize_density(binSize, estimates)
  if(delta!=0)
  {
    estimates <- remove_bumps(binSize, estimates, delta)
    estimates <- normalize_density(binSize, estimates)
  }
  if(alpha!=1)
  {
    estimates <- sharpen(binSize, estimates,alpha)
    estimates <- normalize_density(binSize, estimates)
  }

  return(estimates)
}


.findThresholdHPD=function(binSize,estimates,confidence)
{
  estimates=as.vector(estimates)
  maxDensity=max(estimates)
  minDensity=min(estimates)
  newCut=(maxDensity+minDensity)/2
  eps=1
  ii=1
  while(ii<=1000)
  {
    prob=sum(binSize*estimates*(estimates>newCut))
    eps=abs(confidence-prob)
    if(eps<0.0000001) break; # level found
    if(confidence>prob) maxDensity=newCut
    if(confidence<prob) minDensity=newCut
    newCut=(maxDensity+minDensity)/2
    ii=ii+1
  }
  return(newCut)
}

.findThresholdSymmetricMode=function(binSize,estimates,confidence)
{
  estimates=as.vector(estimates)
  mode=which.max(estimates)
  maxInverval=length(estimates)
  minInverval=1
  newCut=round(maxInverval/2)
  eps=1
  ii=1
  while(ii<=1000)
  {
    whichRegion=round(max(c(1,mode-newCut)):min(c(length(estimates),mode+newCut)))
    prob=sum(binSize*estimates[whichRegion])
    eps=abs(confidence-prob)
    if(eps<0.0000001) break; # level found
    if(confidence>prob) minInverval=newCut
    if(confidence<prob) maxInverval=newCut
    newCut=(maxInverval+minInverval)/2
    ii=ii+1
  }
  return(c(whichRegion[1],whichRegion[length(whichRegion)]))
}
