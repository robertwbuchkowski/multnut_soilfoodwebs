# Development of a new model with multiple elements:

# Creating a simple three node community:

# Create a data frame of feeding relationships:
feeding = data.frame(
Predator = c("Pred", "Pred", "Prey1", "Prey2", "Prey2","Prey1", "Microbe1"),
Prey = c("Prey1", "Prey2", "Prey2", "Microbe1","Detritus","Detritus","Detritus"),
Preference = c(1,1.2,1,1,1,1,1))

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
source("Scripts/correct_diet.R")

usin2 = correct_diet(usin)

comana(usin2)

# Test the function with a simple two element community:

# Create a data frame of feeding relationships:
feeding = data.frame(
  Predator = c("Pred", "Pred", "Prey1", "Prey2", "Prey2","Prey1", "Microbe1"),
  Prey = c("Prey1", "Prey2", "Prey2", "Microbe1","Detritus","Detritus","Detritus"),
  Preference = c(1,1.2,1,1,1,1,1))

# Read back in the example data:
properties = read.csv("Data/properties2.csv")

# Build the foodweb:
source("Scripts/build_foodweb.R")
usin <- build_foodweb(feeding = feeding,
                      properties = properties)

# Clean the environment:
rm(properties, feeding)

# usin$prop$Nitrogen$a = 1

source("Scripts/comana.R")
comana(usin)

# Correct the diet to fix stoichiometry:
source("Scripts/correct_diet.R")

usin2 = correct_diet(usin)

comana(usin2)

# CORRECT DIET IS NOT WORKING, BECAUSE THERE SHOULD BE A ZERO MINERALIZATION POSSIBLE FOR THE LOWEST ELEMENT. ERROR IN THE PROGRAM!





# Test the function with a simple two element community:

# Create a data frame of feeding relationships:
feeding = data.frame(
  Predator = c("Pred", "Pred", "Prey1", "Prey2", "Prey2","Prey1", "Microbe1"),
  Prey = c("Prey1", "Prey2", "Prey2", "Microbe1","Detritus","Detritus","Detritus"),
  Preference = c(1,1.2,1,1,1,1,1))

# Read back in the example data:
properties = read.csv("Data/properties2.csv")

# Build the foodweb:
source("Scripts/build_foodweb.R")
usin <- build_foodweb(feeding = feeding,
                      properties = properties)

# Clean the environment:
rm(properties, feeding)

# usin$prop$Nitrogen$a = 1

source("Scripts/comana.R")
comana(usin)

# usin$imat[2,3] = 209.7189 # The correct rate
# comana(usin)$mineralization$Nitrogen

# Correct the diet to fix stoichiometry:
source("Scripts/correct_diet.R")

usin2 = correct_diet(usin)

comana(usin2)

# CORRECT DIET IS NOT WORKING, BECAUSE THERE SHOULD BE A ZERO MINERALIZATION POSSIBLE FOR THE LOWEST ELEMENT. ERROR IN THE PROGRAM!
