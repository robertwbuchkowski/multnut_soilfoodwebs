
# Load parameters:
ei = rbind(matrix(0, nrow = 4, ncol = 4), c(150, 150*0.05, 150*0.016, 150*0.018))
ii = c(0,10,0.5,0.01)
io = c(0,0.1, 0.1, 0.1)
intro_comm$prop$general$Nitrogen$canIMM[3] = 1
parameters = getPARAMS(intro_comm, externalinputs = ei, inorganicinputs = ii, inorganicloss = io, densitydependence = c(0,0,0,0,0))

(test1 = foodwebode(t = 1, y = parameters$yeqm, pars = parameters$parameters))

(test1 = foodwebode(t = 1, y = parameters$yeqm*c(rep(c(1,1,1,1,1.1), 4),1,1,1,1), pars = parameters$parameters))

test1[[1]][abs(test1[[1]]) < 1e-8] = 0

subset(test1[[1]], test1[[1]] != 0)


parameters2 = getPARAMS(intro_comm, externalinputs = ei, inorganicinputs = ii, inorganicloss = io, densitydependence = c(1,1,1,1,0))

(test2 = foodwebode(t = 1, y = parameters2$yeqm, pars = parameters2$parameters))

test2[[1]][abs(test2[[1]]) < 1e-8] = 0

subset(test2[[1]], test2[[1]] != 0)

parameters3 = getPARAMS(intro_comm, externalinputs = ei, inorganicinputs = ii, inorganicloss = io, functionalresponse = matrix(0.01, nrow = 5, ncol = 5))

(test3 = foodwebode(t = 1, y = parameters3$yeqm, pars = parameters3$parameters))

test3[[1]][abs(test3[[1]]) < 1e-8] = 0

subset(test3[[1]], test3[[1]] != 0)

# Test a few simulations:
library(deSolve)

m1 = ode(y = parameters$yeqm*c(rep(c(1,1,1,1,1.1), 4),1,1,1,1), times = 1:50, func =foodwebode, parms = parameters$parameters)

m2 = ode(y = parameters2$yeqm*c(rep(c(1,1,1,1,1.1), 4),1,1,1,1), times = 1:50, func =foodwebode, parms = parameters2$parameters)

m3 = ode(y = parameters3$yeqm*c(rep(c(1,1,1,1,1.1), 4),1,1,1,1), times = 1:50, func =foodwebode, parms = parameters3$parameters)

library(tidyverse)
tibble(data.frame(m1)) %>%
  mutate(Model = "A") %>%
  bind_rows(
    tibble(data.frame(m2)) %>%
      mutate(Model = "B")
  )%>%
  bind_rows(
    tibble(data.frame(m3)) %>%
      mutate(Model = "C")
  ) %>%
  pivot_longer(!time & !Model) %>%
  separate(name, into = c("Node","Element"), sep = "_") %>%
  filter(Model != "C") %>%
  ggplot(aes(x = time, y = value, color = Model, group = Model)) + geom_line() + facet_wrap(Node~Element, scales = "free")


m1a = ode(y = parameters$yeqm*c(rep(c(1,1,1,1,1.1), 4),1,1,1,0), times = 1:50, func =foodwebode, parms = parameters$parameters)



# List of matrices
matrices <- list(A, B, C)

expand_mat <- function(matrices){

  matrix_rows = nrow(matrices[[1]])
  matrix_cols = ncol(matrices[[1]])

  # Define the larger matrix
  larger_matrix <- matrix(0, nrow = length(matrices) * matrix_rows, ncol = length(matrices) * matrix_cols)

  # Place matrices on the diagonal of the larger matrix
  for (i in 1:length(matrices)) {
    start_row <- (i - 1) * matrix_rows + 1
    end_row <- i * matrix_rows
    start_col <- (i - 1) * matrix_cols + 1
    end_col <- i * matrix_cols
    larger_matrix[start_row:end_row, start_col:end_col] <- matrices[[i]]
  }

  return(larger_matrix)
}

expand_mat(list(matrix(data = 1, nrow = 2, ncol = 3),
                matrix(data = 2, nrow = 2, ncol = 3),
                matrix(data = 3, nrow = 2, ncol = 3)))




mineralization_list = vector('list', length = dim(Qmat)[1])

for(ii in 1:length(mineralization_list)){
  f.rhs.a = unname(c((netwithoutmineralization[,1]*Qmat - netwithoutmineralization)[ii,-1]))

  f.con.a = cbind(Qmat[ii,-1], diag(x=-1, nrow = dim(Qmat)[2]-1, ncol = dim(Qmat)[2]-1))

  f.dir.a = rep("=", dim(f.con.a)[1])

  f.rhs.b = rep(0, dim(Qmat)[2]) # where to add inorganic...

  f.con.b = diag(x=1, nrow = dim(Qmat)[2], ncol = dim(Qmat)[2])

  f.dir.b = rep(">=", dim(Qmat)[2])

  f.obj = rep(0, dim(Qmat)[2]); f.obj[1] = 1

  min_sol = lpSolve::lp(direction = "min", f.obj,
                        rbind(f.con.a, f.con.b),
                        c(f.dir.a, f.dir.b),
                        c(f.rhs.a, f.rhs.b))

  mineralization_list[[ii]] = min_sol$solution
}




f.rhs.a = c((netwithoutmineralization[,1]*Qmat - netwithoutmineralization)[,-1])

f.con.a = cbind(
  c(t(Qmat[,-1])),
  kronecker(matrix(1, nrow = dim(Qmat)[1], ncol = 1), diag(x=-1, nrow = dim(Qmat)[2]-1, ncol = dim(Qmat)[2]-1))
)

f.dir.a = rep("=", length(f.rhs.a))

f.rhs.b = rep(0, length(Qmat))

f.con.b = cbind(
  rep(0, length(Qmat)),
  kronecker(matrix(1, nrow = dim(Qmat)[1], ncol = 1), diag(x=1, nrow = dim(Qmat)[2], ncol = dim(Qmat)[2]))
)

f.rhs.b = rep(">=", length(Qmat))





# CODE FOR THE MINERALIZATION LINEAR PROGRAM:
# > library(multnutFW)
# > # Load parameters:
#   > ei = rbind(matrix(0, nrow = 4, ncol = 4), c(150, 150*0.05, 150*0.016, 150*0.018))
# > ii = c(0,10,0.5,0.01)
# > io = c(0,0.1, 0.1, 0.1)
# > parameters = getPARAMS(intro_comm, externalinputs = ei, inorganicinputs = ii, inorganicloss = io, densitydependence = c(0,0,0,0,0))
# > # Test a few simulations:
#   > library(deSolve)
# > m1 = ode(y = parameters$yeqm*c(rep(c(1,1,1,1,1.1), 4),1,1,1,1), times = 1:50, func =foodwebode, parms = parameters$parameters)
# Called from: func(time, state, parms, ...)
# Browse[1]> netwithoutmineralization
# Carbon     Nitrogen    Phosphorus      Calcium
# [1,]  2.289835e-15  0.006184426  5.724587e-17  0.001314754
# [2,]  2.170695e-01  0.011688356  1.259749e-02  0.026247480
# [3,]  1.823435e+00  0.141822731  5.267910e-01  0.592639846
# [4,]  4.857286e+00 -7.177002845  9.371657e-02  0.105431145
# [5,] -2.169079e+01 -1.854349809 -3.239330e-01 -0.363085669
# Browse[1]> Qmat
# [,1]  [,2]  [,3]  [,4]
# [1,]    1 0.222 0.030 0.020
# [2,]    1 0.208 0.020 0.022
# [3,]    1 0.200 0.016 0.018
# [4,]    1 0.200 0.016 0.018
# [5,]    1 0.050 0.016 0.018
# Browse[1]> netwithoutmineralization[,1]
# [1]  2.289835e-15  2.170695e-01  1.823435e+00  4.857286e+00
# [5] -2.169079e+01
# Browse[1]> netwithoutmineralization[,1]*Qmat
# [,1]          [,2]          [,3]          [,4]
# [1,]  2.289835e-15  5.083434e-16  6.869505e-17  4.579670e-17
# [2,]  2.170695e-01  4.515045e-02  4.341390e-03  4.775528e-03
# [3,]  1.823435e+00  3.646870e-01  2.917496e-02  3.282183e-02
# [4,]  4.857286e+00  9.714572e-01  7.771657e-02  8.743114e-02
# [5,] -2.169079e+01 -1.084540e+00 -3.470526e-01 -3.904342e-01
# Browse[1]> netwithoutmineralization[,1]*Qmat - netwithoutmineralization
# Carbon     Nitrogen    Phosphorus      Calcium
# [1,]      0 -0.006184426  1.144917e-17 -0.001314754
# [2,]      0  0.033462095 -8.256097e-03 -0.021471951
# [3,]      0  0.222864291 -4.976160e-01 -0.559818014
# [4,]      0  8.148460009 -1.600000e-02 -0.018000000
# [5,]      0  0.769810289 -2.311960e-02 -0.027348558
# Browse[1]> ?solve
# Browse[1]> ?lpSolve::lp
# Browse[1]> Q
# > f.obj = c(1,0,0)
# > f.con = c(0.222, -1, 0, 0,
#             +           0.030, 0, -1, 0,
#             +           0.020, 0, 0, -1)
# > f.obj = c(1,0,0,0)
# > f.dir = c("=", "=", "=")
# > f.con = matrix(c(0.222, -1, 0, 0,
#                    +           0.030, 0, -1, 0,
#                    +           0.020, 0, 0, -1), nrow = 3, byrow = T)
# > f.rhs = c(-0.006184426,1.144917e-17,-0.001314754)
# > f.con = matrix(c(0.222, -1, 0, 0,
#                    +           0.030, 0, -1, 0,
#                    +           0.020, 0, 0, -1,), nrow = 3, byrow = T)
# Error in c(0.222, -1, 0, 0, 0.03, 0, -1, 0, 0.02, 0, 0, -1, ) :
#   argument 13 is empty
# > f.con = matrix(c(0.222, -1, 0, 0,
#                    +           0.030, 0, -1, 0,
#                    +           0.020, 0, 0, -1,
#                    +           1,0,0,0,
#                    +           0,1,0,0,
#                    +           0,0,1,0,
#                    +           0,0,0,1), nrow = 7, byrow = T)
# > f.con
# [,1] [,2] [,3] [,4]
# [1,] 0.222   -1    0    0
# [2,] 0.030    0   -1    0
# [3,] 0.020    0    0   -1
# [4,] 1.000    0    0    0
# [5,] 0.000    1    0    0
# [6,] 0.000    0    1    0
# [7,] 0.000    0    0    1
# > f.dir = c("=", "=", "=", ">=", ">=", ">=", ">=")
# > f.rhs
# [1] -6.184426e-03  1.144917e-17 -1.314754e-03
# > f.rhs = c(-0.006184426,1.144917e-17,-0.001314754, 0,0,0,0)
# > lp ("max", f.obj, f.con, f.dir, f.rhs)
# Error in lp("max", f.obj, f.con, f.dir, f.rhs) :
#   could not find function "lp"
# > lpSolve::lp ("max", f.obj, f.con, f.dir, f.rhs)
# Error: status 3
# > f.obj
# [1] 1 0 0 0
# > f.con
# [,1] [,2] [,3] [,4]
# [1,] 0.222   -1    0    0
# [2,] 0.030    0   -1    0
# [3,] 0.020    0    0   -1
# [4,] 1.000    0    0    0
# [5,] 0.000    1    0    0
# [6,] 0.000    0    1    0
# [7,] 0.000    0    0    1
# > f.rhs
# [1] -6.184426e-03  1.144917e-17 -1.314754e-03  0.000000e+00
# [5]  0.000000e+00  0.000000e+00  0.000000e+00
# > f.dir
# [1] "="  "="  "="  ">=" ">=" ">=" ">="
# > lpSolve::lp ("max", f.obj, f.con, f.dir, f.rhs)
# Error: status 3
# > lpSolve::lp("max", f.obj, f.con, f.dir, f.rhs)
# Error: status 3
# > lpSolve::lp("min", f.obj, f.con, f.dir, f.rhs)
# Success: the objective function is 0
# > sol = lpSolve::lp("min", f.obj, f.con, f.dir, f.rhs)
# > sol
# Success: the objective function is 0
# > sol$objective
# [1] 1 0 0 0
# > sol$solution
# [1] 0.000000000 0.006184426 0.000000000 0.001314754
# > sol$solution*f.con
# [,1]         [,2]         [,3]         [,4]
# [1,] 0.0000000000 -0.001314754  0.000000000  0.000000000
# [2,] 0.0001855328  0.000000000 -0.001314754  0.000000000
# [3,] 0.0000000000  0.000000000  0.000000000 -0.001314754
# [4,] 0.0013147540  0.000000000  0.000000000  0.000000000
# [5,] 0.0000000000  0.001314754  0.000000000  0.000000000
# [6,] 0.0000000000  0.000000000  0.001314754  0.000000000
# [7,] 0.0000000000  0.000000000  0.000000000  0.001314754
# > rowSums(sol$solution*f.con)
# [1] -0.001314754 -0.001129221 -0.001314754  0.001314754  0.001314754
# [6]  0.001314754  0.001314754
# > f.rhs
# [1] -6.184426e-03  1.144917e-17 -1.314754e-03  0.000000e+00
# [5]  0.000000e+00  0.000000e+00  0.000000e+00
# > f.con
# [,1] [,2] [,3] [,4]
# [1,] 0.222   -1    0    0
# [2,] 0.030    0   -1    0
# [3,] 0.020    0    0   -1
# [4,] 1.000    0    0    0
# [5,] 0.000    1    0    0
# [6,] 0.000    0    1    0
# [7,] 0.000    0    0    1
# > 0.006184426*-1
# [1] -0.006184426
# > sol$solution%*%f.con
# Error in sol$solution %*% f.con : non-conformable arguments
# > f.con%*%sol$solution
# [,1]
# [1,] -0.006184426
# [2,]  0.000000000
# [3,] -0.001314754
# [4,]  0.000000000
# [5,]  0.006184426
# [6,]  0.000000000
# [7,]  0.001314754
# > sol = lpSolve::lp("min", f.obj, f.con, f.dir, f.rhs)
# > sol$sens.coef.to
# [1] 0
# > ?lp
# > sol = lpSolve::lp("min", f.obj, f.con, f.dir, f.rhs, compute.sens = T)
# > sol$compute.sens
# [1] 1
# > sol$sens.coef.from
# [1]   0.000000  -4.504505   0.000000 -50.000000
# > sol = lpSolve::lp("min", f.obj, f.con, f.dir, f.rhs)
