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

feedingmatrix = matrix(0,
                       dimnames = list(unique(properties$ID), unique(properties$ID)),
                       ncol = length(unique(properties$ID)), nrow = length(unique(properties$ID)), byrow = TRUE)

#Run through the feeding list and add in the non-zero feeding links identified there.
for(i in 1:dim(feeding)[1]){
  feedingmatrix[feeding$Predator[i],feeding$Prey[i]] = feeding$Preference[i]
}

# Build the community
usin <- list(imat = feedingmatrix,
                  prop = properties)

# Clean the environment:
rm(feedingmatrix, properties, feeding)


# COPY OVER THE FUNCTION HERE TO CHECK THE COMMUNITY. DO THIS OUTSIDE OF COMANA.


# Community analysis function:
comana <- function(usin){
  # Separate the imat and prop:
  imat = usin$imat # row values of imat sets predator feeding preferences!
  prop = usin$prop # properties of each trophic species

  # Add in perfect production efficiency if p is not listed in the data frame:
  if(!any(prop$Parameter == "p")){
    tdf = unique(prop[,c("ID", "Element")])

    tdf[,c("Parameter")] = c("p")
    tdf[,c("Value")] = c(1)

    prop = rbind(prop, tdf)
    rm(tdf)
  }

  # Add in zero respiration if E is not listed in the data frame:
  if(!any(prop$Parameter == "E")){
    tdf = unique(prop[,c("ID", "Element")])

    tdf[,c("Parameter")] = c("E")
    tdf[,c("Value")] = c(0)

    prop = rbind(prop, tdf)
    rm(tdf)
  }

  # Test conflict between E and p:
  if(any(subset(prop, prop$Parameter == "E")$Value > 0) & any(subset(prop, prop$Parameter == "p")$Value < 1)) print("Error")

  Nnodes = dim(imat)[1] # Number of nodes in the food web

  # Break out the elements into a properties list for easier processing:
  element_list = unique(prop$Element)

  prop2 = vector(mode = "list", length = length(element_list))
  names(prop2) = element_list

  for(i in 1:length(prop2)){
    prop2[[i]] = subset(prop, prop$Element == element_list[i])

    prop2[[i]] = reshape(prop2[[i]][, c("ID", "Parameter", "Value")],
                         idvar = "ID",
                         timevar = "Parameter",
                         direction = "wide")

    names(prop2[[i]]) <- gsub("Value.", "", names(prop2[[i]]))
  }

  prop = prop2; rm(prop2)

  # Create a vector for the consumption rates
  temp_mat =
    -1*t(imat)*matrix(prop$Carbon$B, nrow = Nnodes, ncol = Nnodes)/matrix(rowSums(imat*matrix(prop$Carbon$B, nrow = Nnodes, ncol = Nnodes, byrow = T)), nrow = Nnodes, ncol = Nnodes, byrow = T)

  temp_mat[!is.finite(temp_mat)] = 0 # Replace non-finite values with 0 because total consumption was zero in this case

  diag(temp_mat) = prop$Carbon$a*prop$Carbon$p + diag(temp_mat) # Add in a*p term

  consumption = base::solve(temp_mat,(prop$Carbon$d*prop$Carbon$B + prop$Carbon$E*prop$Carbon$B))

  # Confirm that this solution is unique by showing Ax = 0 produces x = 0
  if(any(solve(temp_mat,rep(0, Nnodes)) != 0)){
    warning("Solution to the web is not unique!")
  }

  names(consumption) = colnames(imat) # Names match the trophic species names

  # Create an fmat vector
  fmat = vector(mode = "list", length = length(element_list))
  names(fmat) = element_list

  # Create a mineralization vector
  mineralization = fmat

  # Create a new matrix for feeding rates with the same dimensions as the food web matrix
  fmat[[1]] = (imat*matrix(prop$Carbon$B, nrow = Nnodes, ncol = Nnodes, byrow=TRUE)/rowSums(imat*matrix(prop$Carbon$B, nrow = Nnodes, ncol = Nnodes, byrow=TRUE)))*consumption

  fmat[[1]][!is.finite(fmat[[1]])] = 0 # Replace NaN with 0 for no feeding

  # Fix the detritus calculations: detritus receives dead material from all other trophic levels and so consumption is the losses minus the inputs it already gets from inefficient eating and dead biomass. This value can be negative indicating that inputs from outside the ecosystem are necessary to meet internal demand for C.
  if(any(prop$Carbon$isDetritus >0)){
    detritusPOS = which(prop$Carbon$isDetritus >0)

    for(i in detritusPOS){
      consumption[i] = sum(fmat[[1]][,i]) - prop$Carbon$DetritusRecycling[i]*(sum((1-prop$Carbon$a)*consumption) + sum(prop$Carbon$d*prop$B))
    }
  }

  # Calculate carbon mineralization using the production efficiency
  mineralization[[1]] = prop$Carbon$a*(1-prop$Carbon$p)*consumption + prop$Carbon$E*prop$Carbon$B

  # Calculate the mineralization rates of the various elements using the comparison to carbon:

  for(i in 2:length(element_list)) {
    current_element_properties = prop[[which(names(prop) == element_list[i])]]

    Qhat = prop$Carbon$Q/current_element_properties$Q

    mineralization[[i]] = (prop$Carbon$E*prop$Carbon$B - # Carbon mineralization rate based on a fixed proportion of biomass
                             rowSums((matrix(Qhat, nrow = Nnodes, ncol = Nnodes)* # Predator C:X ratio
                                        matrix(current_element_properties$a*current_element_properties$p, nrow = Nnodes, ncol = Nnodes)/ # Predator X assimilation and production rates
                                        matrix(Qhat, nrow = Nnodes, ncol = Nnodes, byrow = T) - # Prey C:X ratio
                                        matrix(prop$Carbon$a*prop$Carbon$p, nrow = Nnodes, ncol = Nnodes))* # Predator C assimilation and production rates
                                       fmat$Carbon))* # consumption rates
      Qhat* # multiply by C:X ratio to get back to units of X
      as.numeric(rowSums(imat)>0) # Make X mineralization zero for all nodes without prey items.
  }

 return(list(fmat = fmat, consumption = consumption, mineralization = mineralization))
}


# Correct the diet to fix stoichiometry:

# Separate the imat and prop:
imat = usin$imat # row values of imat sets predator feeding preferences!
prop = usin$prop # properties of each trophic species
mineralization = comana(usin)$mineralization
consumption = comana(usin)$consumption
Nnodes = dim(imat)[1] # Number of nodes in the food web

# just make one that is everything:
dietlimits = imat
dietlimits[1,3] = 1


# Break out the elements into a properties list for easier processing:
element_list = unique(prop$Element)

prop2 = vector(mode = "list", length = length(element_list))
names(prop2) = element_list

for(i in 1:length(prop2)){
  prop2[[i]] = subset(prop, prop$Element == element_list[i])

  prop2[[i]] = reshape(prop2[[i]][, c("ID", "Parameter", "Value")],
                       idvar = "ID",
                       timevar = "Parameter",
                       direction = "wide")

  names(prop2[[i]]) <- gsub("Value.", "", names(prop2[[i]]))
}

prop = prop2; rm(prop2)

# Build the matrices Aij for each element:

# NEED TO KEEP THESE ALL IN RELATIVE CARBON UNITS. DO NOT CONVERT TO OTHER ELEMENTS, FOR CLEAR COMPARISON.
AIJ = vector(mode = "list", length = length(element_list))

for(i in 2:length(element_list)) {
  current_element_properties = prop[[which(names(prop) == element_list[i])]]

  Qhat = prop$Carbon$Q/current_element_properties$Q

  AIJ[[i]] = (matrix(Qhat, nrow = Nnodes, ncol = Nnodes)* # Predator C:X ratio
                matrix(current_element_properties$a*current_element_properties$p, nrow = Nnodes, ncol = Nnodes)/ # Predator X assimilation and production rates
                matrix(Qhat, nrow = Nnodes, ncol = Nnodes, byrow = T) - # Prey C:X ratio
                matrix(prop$Carbon$a*prop$Carbon$p, nrow = Nnodes, ncol = Nnodes))

    # (prop$Carbon$E*prop$Carbon$B - # Carbon mineralization rate based on a fixed proportion of biomass
    #                        rowSums((matrix(Qhat, nrow = Nnodes, ncol = Nnodes)* # Predator C:X ratio
    #                                   matrix(current_element_properties$a*current_element_properties$p, nrow = Nnodes, ncol = Nnodes)/ # Predator X assimilation and production rates
    #                                   matrix(Qhat, nrow = Nnodes, ncol = Nnodes, byrow = T) - # Prey C:X ratio
    #                                   matrix(prop$Carbon$a*prop$Carbon$p, nrow = Nnodes, ncol = Nnodes))* # Predator C assimilation and production rates
    #                                  fmat$Carbon))* # consumption rates
    # Qhat* # multiply by C:X ratio to get back to units of X
    # as.numeric(rowSums(imat)>0) # Make X mineralization zero for all nodes without prey items.
}

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
    ai = AIJ[sp,food]
    biomass = prop$Carbon$B[food] # Get the biomass of food items
    FT = consumption[sp]
    limits = dietlimits[sp,food]

    # Try to see whether there is a solution where the diet can be optimized
    sol1 <- NULL
    try(

      # Build the constraints, because they change depending on how many elements are part of the model:


      # WORKING ABOVE TO GET THIS RIGHT BEFOREHAND: AIJ goes in the AMAT, while carbon respiration goes in BVEC. NEED TO KEEP ALL THESE IN RELATIVE CARBON UNITS!


      sol1 <- quadprog::solve.QP(Dmat = diag(length(ai)), # We need the squared terms for each diet, use the identity matrix to get them f^T I f
                                 dvec = biomass/sum(biomass), # The diet is as close to the relative abundance as possible
                                 Amat = t(
                                   rbind(rep(1, length(ai)), # Constraint that consumption proportions sum to 1
                                         ai*FT, # Nmin, but need one for each
                                         -diag(nrow = length(ai)), # Diet limits
                                         diag(nrow = length(ai)))), # proportions sum to 1, Nmin is zero, and none of the limits are exceeded
                                 bvec = c(1, 0, -limits, rep(0,length(ai))), # proportions sum to 1, Nmin is zero, and none of the limits are exceeded, and all values greater than zero
                                 meq = 1 # first position in Amat is an equality constraint)
      ),
      silent = T)

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
