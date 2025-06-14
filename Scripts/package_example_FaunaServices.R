# Demonstrate MultnutFW code
# A: Robert Buchkowski

# Load the package:
library(multnutFW)
library(igraph)

## -------------------------------------------------------------
# Define edges: isopods eat detritus
edges <- c("detritus","isopod", "fungi", "isopod", "detritus", "fungi")
g <- make_graph(edges = edges, directed = TRUE)

# Define custom layout: isopods on top, detritus on bottom
layout_matrix <- matrix(c(
  -1, 0,   # detritus (bottom left)
   0, 1,   # isopod (top center)
   1, 0.25    # fungi (bottom right)
), ncol = 2, byrow = TRUE)

# Plot the graph
plot(g,
     layout = layout_matrix,
     vertex.color = "lightblue",
     vertex.size = 100,
     vertex.label.cex = 1.5,
     edge.arrow.size = 0.6,
     main = "Isopod-Detritus Network")


## -------------------------------------------------------------
# Create the tibble
interaction_table <- data.frame(
  Predator = c("isopod", "isopod", "fungi"),
  Prey = c("detritus", "fungi", "detritus"),
  Preference = c(1, 1, 1),
  aCarbon = c(0.3, 0.8, 0.8),
  aNitrogen = c(0.25, 0.75, 0.9),
  aPhosphorus = c(0.22, 0.73, 0.91)
)

# Display the tibble
structure(interaction_table)


## -------------------------------------------------------------
# Create the tibble
parameter_table <- data.frame(
  ID = c("isopod","fungi","detritus","isopod","fungi","detritus","isopod","fungi","detritus","isopod","fungi","detritus","isopod","fungi","detritus","isopod","fungi","detritus","isopod","fungi","detritus","isopod","fungi","detritus","isopod","fungi","detritus","isopod","fungi","detritus","isopod","fungi","detritus","isopod","fungi","detritus"),
  Element = c("Carbon", "Carbon", "Carbon", "Nitrogen","Nitrogen","Nitrogen", "Phosphorus","Phosphorus","Phosphorus","Carbon", "Carbon", "Carbon", "Nitrogen","Nitrogen","Nitrogen", "Phosphorus","Phosphorus","Phosphorus","Carbon", "Carbon", "Carbon","Carbon", "Carbon", "Carbon","Carbon", "Carbon", "Carbon","Carbon", "Carbon", "Carbon","Carbon", "Carbon", "Carbon","Carbon", "Carbon", "Carbon"),
  Parameter = c(rep("Q", 9), rep("canIMM", 9), rep("E",3), rep("d",3), rep("B",3),rep("DetritusRecycling",3), rep("isDetritus",3), rep("isPlant",3)),
  Value = c(0.5, 0.5, 0.5, 0.1, 0.1, 0.025, 0.015, 0.01, 0.008, 0,0,0,0,1,0,0,1,0, 0.2, 0.05, 0, 1, 0.2, 0, 5, 20, 100,0,0,1,0,0,1,0,0,0)
)

# Display the tibble
structure(parameter_table)


## ----echo = T-------------------------------------------------
# Assemble the food web:
foodweb = build_foodweb(feeding = as.data.frame(interaction_table), properties = as.data.frame(parameter_table))

structure(foodweb)


## -------------------------------------------------------------
trace(multnutFW::comana, quote(print(temp_mat)),at = 12)
output = comana(foodweb)
untrace(comana)


## -------------------------------------------------------------
trace(multnutFW::comana, quote(print((prop$Carbon$d * prop$Carbon$B +
        prop$Carbon$E * prop$Carbon$B + prop$Carbon$Ehat * prop$Carbon$B))),at = 12)
output = comana(foodweb)
untrace(comana)


## ----echo = T-------------------------------------------------
comana(foodweb)$consumption


## ----echo = T-------------------------------------------------
comana(foodweb)$fmat$Carbon


## ----echo = T-------------------------------------------------
comana(foodweb)$fmat$Nitrogen


## ----echo = T-------------------------------------------------
comana(foodweb)$AIJ$Nitrogen
comana(foodweb)$AIJ$Phosphorus


## ----echo = T-------------------------------------------------
comana(foodweb)$mineralization


## ----echo = T-------------------------------------------------
correct_diet(foodweb)$imat

## ----echo = T-------------------------------------------------
comana(correct_diet(foodweb))$mineralization


## ----echo = T-------------------------------------------------
foodweb2 = correct_respiration(foodweb)
foodweb2$prop$general$Carbon


## ----echo = T-------------------------------------------------
comana(foodweb2)$mineralization
