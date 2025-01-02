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

# Correct production efficiency:
source("Scripts/correct_productionefficiency.R")
debugonce(correct_productionefficiency)
usin3 = correct_productionefficiency(usin)

comana(usin3)$mineralization
comana(usin)$mineralization

comana(usin)$consumption
comana(usin3)$consumption


# Test p corrections for predator:

sptt = 3

usin_temp = usin
output = matrix(NA, nrow = 100, ncol = 5)
output[,1] = seq(0.01, 1, length = 100)


for(i in 1:100){

  usin_temp$prop$Carbon$p[sptt] = output[i,1]

  output[i,2:5] = do.call("c", lapply(comana(usin_temp)$mineralization, function(x) unname(x[sptt])))

}
output
usin3$prop$Carbon$p[sptt]

# Correct respiration:
source("Scripts/correct_respiration.R")
debugonce(correct_respiration)
usin3 = correct_respiration(usin)

comana(usin3)

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

# Correct respiration:
source("Scripts/correct_respiration.R")

debugonce(correct_respiration)
usin3 = correct_respiration(usin)

comana(usin3)

# The consumption rate changes, so the calculation in the correct_respiration function is not accurate. This is because consumption cannot be factored out when respiration rate maintains a second term in the equation...in other words, Fij is a function of E_[C,i] so both need to be changed together.

comana(usin)$consumption
comana(usin3)$consumption
