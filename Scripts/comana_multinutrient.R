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
    Parameter = c(rep(c("a","Rmin","Q", "canIMM"),3*4), rep(c("d", "B", "DetritusRecycling", "isDetritus", "isPlant"),3)))

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

# Separate the imat and prop:
imat = usin$imat # row values of imat sets predator feeding preferences!
prop = usin$prop # properties of each trophic species


Nnodes = dim(imat)[1] # Number of nodes in the food web



# Create a vector for the consumption rates
temp_mat =
  -1*t(imat)*matrix(prop$B, nrow = Nnodes, ncol = Nnodes)/matrix(rowSums(imat*matrix(prop$B, nrow = Nnodes, ncol = Nnodes, byrow = T)), nrow = Nnodes, ncol = Nnodes, byrow = T)

temp_mat[!is.finite(temp_mat)] = 0 # Replace non-finite values with 0 because total consumption was zero in this case

diag(temp_mat) = prop$a*prop$p + diag(temp_mat) # Add in a*p term

consumption = base::solve(temp_mat,(prop$d*prop$B))

# Confirm that this solution is unique by showing Ax = 0 produces x = 0
if(any(solve(temp_mat,rep(0, Nnodes)) != 0)){
  warning("Solution to the web is not unique!")
}

names(consumption) = colnames(imat) # Names match the trophic species names

# Create a new matrix for feeding rates with the same dimensions as the food web matrix
fmat = (imat*matrix(prop$B, nrow = Nnodes, ncol = Nnodes, byrow=TRUE)/rowSums(imat*matrix(prop$B, nrow = Nnodes, ncol = Nnodes, byrow=TRUE)))*consumption

fmat[!is.finite(fmat)] = 0 # Replace NaN with 0 for no feeding

# Fix the detritus calculations: detritus receives dead material from all other trophic levels and so consumption is the losses minus the inputs it already gets from inefficient eating and dead biomass. This value can be negative indicating that inputs from outside the ecosystem are necessary to meet internal demand for C.
if(any(prop$isDetritus >0)){
  detritusPOS = which(prop$isDetritus >0)

  for(i in detritusPOS){
    consumption[i] = sum(fmat[,i]) - prop$DetritusRecycling[i]*(sum((1-prop$a)*consumption) + sum(prop$d*prop$B))
  }
}

# Calculate carbon mineralization using the production efficiency
Cmin = prop$a*(1-prop$p)*consumption

# Calculate nitrogen mineralization as a matrix so that the nitrogen mineralized by each consumption interaction is recorded. Summing the rows of this matrix gives the nitrogen that is mineralized overall.
Nmin = (matrix(1/prop$CN, nrow = Nnodes, ncol = Nnodes, byrow = T) - matrix(prop$p/prop$CN, nrow=Nnodes, ncol = Nnodes))*matrix(prop$a, nrow=Nnodes, ncol = Nnodes)*fmat

# Calculate the nitrogen flux throughout the matrix as the amount leaving a node including the N that will be mineralized.
Nfmat = matrix(1/prop$CN, nrow = Nnodes, ncol = Nnodes, byrow = TRUE)*fmat
