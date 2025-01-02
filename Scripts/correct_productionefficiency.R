#' A function to reduce growth by increasing respiration.
#'
#' @param usin The input community.
#' @param output_type Should the nutrient limitation be printed (TRUE) or included in the output as a second object (FALSE)?
#' @return The modified community with a higher respiration rate.
correct_productionefficiency = function(usin, output_type = TRUE){

  # Produce a vector to print/output the nutrient limitation of each organism:
  nutlim <- rep(NA, dim(usin$imat)[1])

  #Identify the species that need correction by having negative mineralization and canIMM == 0 and more than 1 prey item
  species = unname(which(apply(do.call("rbind", comana(usin)$mineralization)* # This is the mineralizaiton
                                 do.call("rbind",lapply(usin$prop, function(x) (1-x$canIMM))), # This means that if canIMM == 1 the negative number is multiplied by zero and removed so that the test of needing correction fails. If canIMM ==0, then the numbers are left as is.
                               2, function(x) any(x < 0)) &
                           apply(usin$imat > 0, 1, sum) > 1 # Species must have more than one food item
  ))


  # Produce a vector to print/output the nutrient limitation of each organism:
  nutlim <- rep("Carbon", dim(usin$imat)[1])

  element_list = names(usin$prop) # List of elements

  for(sp in species){

    id_limiting_element = rep(NA, length(element_list))


    # Separate the imat and prop:
    imat = usin$imat # row values of imat sets predator feeding preferences!
    prop = usin$prop # properties of each trophic species
    Nnodes = dim(imat)[1] # Number of nodes in the food web

    # Calculate the maximum production efficiency by element:
    for(i in 2:length(element_list)) {
      current_element_properties = prop[[which(names(prop) == element_list[i])]]

      Qhat = prop$Carbon$Q/current_element_properties$Q

      id_limiting_element[i] = rowSums( # Sum for each species
        (matrix(Qhat, nrow = Nnodes, ncol = Nnodes)* # Predator C:X ratio
        matrix(current_element_properties$a*current_element_properties$p, nrow = Nnodes, ncol = Nnodes)/ # Predator X assimilation and production rates
        matrix(Qhat, nrow = Nnodes, ncol = Nnodes, byrow = T))* # Prey C:X ratio
        comana(usin)$fmat$Carbon # Multiply by fmat carbon to get consumption rates
      )[sp] # Get the current species data
    }

    lim_element_id = which.min(id_limiting_element)

    nutlim[sp] = element_list[lim_element_id]

    p_temp = ((prop$Carbon$E[sp] + prop$Carbon$Ehat[sp])*prop$Carbon$B[sp] +
                min(id_limiting_element, na.rm = T))/
      (prop$Carbon$a[sp]* sum(comana(usin)$fmat$Carbon[sp,]))

    # Check if the p value is appropriate:
    if(p_temp < 0) stop("Error in p calculation, must be positive.")
    if(p_temp > 1) stop("Error in p calculation, must be smaller than 1.")
    if(p_temp > prop$Carbon$p[sp]) warning("Check in p calculation, should be smaller than existing p value.")

    # Replace p value in the new object:
    usin$prop$Carbon$p[sp] = p_temp
  }

  if(output_type){
    print(data.frame(ID = colnames(usin$imat),
                     `Limiting_nutrient` = nutlim))
    return(usin)
  }else{
    return(list(usin,output_type))
  }
}
