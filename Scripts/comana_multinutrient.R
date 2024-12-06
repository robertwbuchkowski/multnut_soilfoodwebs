# Development of a new model with multiple elements:

# Creating a simple three node community:

# Create a data frame of feeding relationships:
feeding = data.frame(
Predator = c("Pred", "Pred"),
Prey = c("Prey1", "Prey2"),
Preference = c(1,1.2))

# Create a data frame of properties for each species:
if(F){
  properties = data.frame(
    ID = c(rep(rep(c("Pred", "Prey1", "Prey2"), each = 4), 4), rep(c("Pred", "Prey1", "Prey2"),each = 5)), # Name
    Element = c(rep(c("Carbon", "Nitrogen", "Phosphorus", "Calcium"), each = 3*4), rep("Carbon", 5*3)),
    Parameter = c(rep(c("a","E","Q", "canIMM"),3*4), rep(c("d", "B", "DetritusRecycling", "isDetritus", "isPlant"),3)))

  write.csv(properties, "Data/properties.csv")
}

# Read back in the example data:

properties = read.csv("Data/properties.csv")


# Build the foodweb:
source("Scripts/build_foodweb.R")
usin <- build_foodweb(feeding = feeding,
              properties = properties)

# Clean the environment:
rm(properties, feeding)

source("Scripts/comana.R")
comana(usin)


# Correct the diet to fix stoichiometry:

# Separate the imat and prop:
imat = usin$imat # row values of imat sets predator feeding preferences!
prop = usin$prop # properties of each trophic species
mineralization = comana(usin)$mineralization
consumption = comana(usin)$consumption
Nnodes = dim(imat)[1] # Number of nodes in the food web
AIJ = comana(usin)$AIJ

# just make one that is everything:
dietlimits = imat
dietlimits[1,3] = 1

#Identify the species that need correction by having negative mineralization and canIMM == 0 and more than 1 prey item
species = unname(which(apply(do.call("rbind", mineralization)* # This is the mineralizaiton
                               do.call("rbind",lapply(prop, function(x) (1-x$canIMM))), # This means that if canIMM == 1 the negative number is multiplied by zero and removed so that the test of needing correction fails. If canIMM ==0, then the numbers are left as is.
                             2, function(x) any(x < 0)) &
                         apply(imat > 0, 1, sum) > 1 # Species must have more than one food item
))

species


for(sp in species){

  while(TRUE){
    food = imat[sp,] > 0 # Get the list of food items
    numfood = sum(food)
    biomass = prop$Carbon$B[food] # Get the biomass of food items
    FT = unname(consumption[sp])
    limits = dietlimits[sp,food]

    # Try to see whether there is a solution where the diet can be optimized
    sol1 <- NULL

    # Build the constraints, because they change depending on how many elements are part of the model:

    AMAT = rbind(rep(1, numfood), # Constraint that consumption proportions sum to 1
                 -diag(nrow = numfood), # Diet limits
                 diag(nrow = numfood))

    for(i in 2:length(AIJ)){
      AMAT = rbind(AMAT, AIJ[[i]][sp, food]*FT)
    }

    # Transpose for the calculation below:
    AMAT = t(AMAT)

    # Calculate the bvec for the optimization:
    BVEC = c(1,
             -unname(limits),
             rep(0,numfood))

    for(i in 2:length(prop)){
      BVEC = c(BVEC, -prop[[i]]$E[sp]*prop$Carbon$B[sp])
    }

    try(
      sol1 <- quadprog::solve.QP(Dmat = diag(numfood), # We need the squared terms for each diet, use the identity matrix to get them f^T I f
                                 dvec = biomass/sum(biomass), # The diet is as close to the relative abundance as possible
                                 Amat = AMAT, # proportions sum to 1, Mineralization is zero, and none of the limits are exceeded
                                 bvec = BVEC, # proportions sum to 1, Mineralization is zero, and none of the limits are exceeded, and all values greater than zero
                                 meq = 1 # first position in Amat is an equality constraint)
      ),
      silent = T)

    #NEED TO CONFIRM THE ABOVE FUNCTION WITH AN EXAMPLE WHERE THERE IS A SOLUTION AND THEN IMPROVE THE LINEAR PROGRAM!


    # If there is no solution, then run a linear program to find the closest diet within the limits
    if(is.null(sol1)){
      sol1 <- lpSolve::lp(direction = "max",
                          objective.in = c(ai*FT),
                          const.mat = rbind(rep(1, length(ai)),-diag(nrow = length(ai)),diag(nrow = length(ai))),
                          const.dir = c("=", rep(">=", 2*length(ai))),
                          const.rhs = c(1, -limits, rep(0, length(ai)))
      )
    }

    # Confirm that solution is positive. Sometimes the solution produces a tiny negative number because of a rounding error (e.g. 1e-18). Set this to zero.
    solcheck = sol1$solution
    if(any(solcheck < -1e-10)) stop("Solution is returning a large negative value. There has been an error in the optimization.")
    if(any(solcheck < 0)){
      solcheck[solcheck < 0] = 0
    }

    if(any(solcheck < 1e-10)){
      warning(paste0("Diet correction removes an item from the diet of species ",sp,". May get strange model behavior. Check outputs for errors in diet proportions! This code just deletes the food item and distirbutes evenly across the other food items."))

      usin$imat[sp,usin$imat[sp,] > 0][which(solcheck < 1e-10)] = 0

    }else{
      break()
    }
  }

  FLIST = solcheck
  FLIST1 = FLIST/biomass
  FLIST1 = min(FLIST1[FLIST1 > 0])
  FLIST2 = (FLIST/biomass)/FLIST1

  usin$imat[sp,food] = FLIST2

  if(sum((comana(usin,shuffleTL = F)$fmat[sp,food]/FT - FLIST)^2) > 1e-10) stop("Check quadratic optimization. There is an issue.")

}
