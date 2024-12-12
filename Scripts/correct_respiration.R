#' A function to reduce growth by increasing respiration.
#'
#' @param usin The input community.
#' @param output_type Should the nutrient limitation be printed (TRUE) or included in the output as a second object (FALSE)?
#' @return The modified community with a higher respiration rate.
correct_respiration = function(usin, output_type = TRUE){

  # Produce a vector to print/output the nutrient limitation of each organism:
  nutlim <- rep(NA, dim(usin$imat)[1])

  #Identify the species that need correction by having negative mineralization and canIMM == 0 and more than 1 prey item
  species = unname(which(apply(do.call("rbind", comana(usin)$mineralization)* # This is the mineralizaiton
                                 do.call("rbind",lapply(usin$prop, function(x) (1-x$canIMM))), # This means that if canIMM == 1 the negative number is multiplied by zero and removed so that the test of needing correction fails. If canIMM ==0, then the numbers are left as is.
                               2, function(x) any(x < 0)) &
                           apply(usin$imat > 0, 1, sum) > 1 # Species must have more than one food item
  ))

  for(sp in species){
    # Separate the imat and prop:
    imat = usin$imat # row values of imat sets predator feeding preferences!
    prop = usin$prop # properties of each trophic species
    # mineralization = comana(usin)$mineralization
    Fij = comana(usin)$fmat$Carbon
    Nnodes = dim(imat)[1] # Number of nodes in the food web
    AIJ = comana(usin)$AIJ

    Ec_increased_to = lapply(AIJ, function(x) {-(1/prop$Carbon$B[sp])*sum(x[sp,]*Fij[sp,])})

    nutlim[sp] = names(AIJ)[which.max(Ec_increased_to)]

    usin$prop$Carbon$E[sp] = max(do.call("rbind",Ec_increased_to))

  }

  nutlim[is.na(nutlim)] = "Carbon"

  if(output_type){
    print(data.frame(ID = colnames(usin$imat),
                     `Limiting_nutrient` = nutlim))
    return(usin)
  }else{
    return(list(usin,output_type))
  }
}
