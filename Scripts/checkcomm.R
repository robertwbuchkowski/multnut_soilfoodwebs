#' A function check the community for errors before it is used in calculations.
#'
#' @param usin The community to check.
#' @param shuffleTL A Boolean stating whether the community should be sorted.
#' @param rmzeros A Boolean determining whether trophic species with zero biomass should be removed from the community.
#' @param verbose A Boolean. Do you want to see all the warnings?
#' @return The checked community.
#' @examples
#' checkcomm(intro_comm)
#' @export
checkcomm <- function(usin, shuffleTL = FALSE, rmzeros = TRUE, verbose = TRUE){


  # Remove zeros from the community
  if(rmzeros & any(usin$prop$B == 0)){

    usin$imat = usin$imat[usin$prop$B != 0,usin$prop$B != 0]

    usin$prop = subset(usin$prop, usin$prop$B !=0)

    if(verbose) warning("Removing zero biomass nodes.")
  }

  # Sort by trophic level
  if(shuffleTL) usin = TLsort(usin)

  # Calculate trophic level
  TL = TLcheddar(usin$imat)

  # Check for several errors in the community and return to user if any is true:

  if(!is.matrix(usin$imat)) stop("Feeding matrix must be of class matrix.")

  if(any(usin$prop$isDetritus >0)){
    # Make sure TL of detritus is 1
    if(!all(TL[usin$prop$isDetritus >0] ==1)) stop("Detirtus must have a trophic level position of 1.")


    # Detritus a and p are 1 and death is zero.
    if(!all(
      subset(usin$prop, usin$prop$isDetritus >0)$a == 1 & subset(usin$prop, usin$prop$isDetritus >0)$p == 1 & subset(usin$prop, usin$prop$isDetritus >0)$d == 0
    )){
      stop("Detritus must have a = 1, p = 1, and d = 0.")
    }

    # Make sure detritus recycling proportion sums to 1
    if(sum(usin$prop$DetritusRecycling) != 1 & verbose){
      warning("Rescaling Detritus Recycling to sum to 1.")

      usin$prop$DetritusRecycling = usin$prop$DetritusRecycling/sum(usin$prop$DetritusRecycling)
    }
  }

  # Properties and matrix have the same individuals
  if(!all(c(colnames(usin$imat) == rownames(usin$imat),colnames(usin$imat) == usin$prop$ID))) stop("Column names, row names, and property data frame IDs must all be the same and in the same order.")

  # Identify the cannibals and mutual predators after a reset in case it has changed since community was created.
  usin$prop$MutualPred = NULL
  if(!("MutualPred" %in% colnames(usin$prop))){
    usin = can_mutfeed(usin)
  }
  return(usin)
}
