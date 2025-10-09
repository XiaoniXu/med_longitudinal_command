
library(tidyverse)
library(simstudy)
seed = 456

t_2 = 4
t_3 = 20

M_0 = 101.59790
M_1_E = 1.26326 
M_1_C = 3.97214


M_2_0 = 46.55165
M_2_E = 0.97056
M_2_C = -1.32312 
M_2_M_1 = 0.56713 

# hazard of D
D_1_0 = -4.867209182
D_1_E = .280396249
D_1_C =  .629431733
D_1_M_1 = -.001901315

D_2_0 = -7.411481
D_2_E =  .297685884
D_2_C =  .510771640 
D_2_M_2 = .007129399

# hazard of Y
Y_1_0 = -7.45413779
Y_1_E = 0.14343187
Y_1_C = 0.48480932
Y_1_M_1 = 0.02250126 


Y_2_0 = -7.798628
Y_2_E = .146895258
Y_2_C = .178517397
Y_2_M_1 = .009772939
Y_2_M_2 = .008042901

set.seed(seed)

## Baseline data definitions
# Why this variances?
def = defData(varname = "E", dist = "normal", formula=2, variance = 0.7) #formula is the mean
def = defData(def, varname = "C", dist = "binary", formula = 0.4)

# M1 definitions
def = defData(def, varname = "M1", dist = "normal", formula = "..M_0 + ..M_1_E * E + ..M_1_C * C", variance = 1)

# M2 definitions
def = defData(def, varname = "M2", dist = "normal", formula = "..M_2_0 + ..M_2_M_1 * M1 + ..M_2_E * E + ..M_2_C * C", variance = 1)

n = 100000 # individuals
dtSurv = genData(n, def)

sdef <- defSurv(varname = "Y", formula = "..Y_1_0 +..Y_1_E * E + ..Y_1_C * C + ..Y_1_M_1 * M1", shape = 1, scale = 1)
sdef <- defSurv(sdef, varname = "D", formula = "..D_1_0 + ..D_1_E * E + ..D_1_C * C + ..D_1_M_1 * M1", shape = 1, scale = 1)
dtSurv <- genSurv(dtSurv, sdef)

# Parallel foreach loop
for (i in 1:length(dtSurv$id)){
  curr_id <- dtSurv$id[i]
  row <- dtSurv[curr_id,]
  curr_Y <- row$Y
  curr_D <- row$D
  
  E <- row$E
  C <- row$C
  M1 <- row$M1
  M2 <- row$M2
  
  if (curr_D > t_2) {
    curr_D <- 0
    while (curr_D < t_2) {
      sdef <- defSurv(varname = "D", formula = "..D_2_0 + ..D_2_E * E + ..D_2_C * C + ..D_2_M_2 * M2", shape = 1, scale = 1)
      dtTem <- genData(1)
      dtTem <- genSurv(dtTem, sdef)
      row$D <- dtTem$D
      curr_D <- row$D
    }
  }
  
  if (curr_Y > t_2) {
    curr_Y <- 0
    while (curr_Y < t_2) {
      sdef <- defSurv(varname = "Y", formula = "..Y_2_0 + ..Y_2_E * E + ..Y_2_C * C + ..Y_2_M_1 * M1 + ..Y_2_M_2 * M2", shape = 1, scale = 1)
      dtTem <- genData(1)
      dtTem <- genSurv(dtTem, sdef)
      row$Y <- dtTem$Y
      curr_Y <- row$Y
    }
  } else {
    curr_Y <- curr_Y
  }
  dtSurv[curr_id,] = row
}

# define death indicator

dtSurv$D1 = ifelse((dtSurv$D <= t_2) & (dtSurv$D <= dtSurv$Y), 1, 0)
dtSurv$D2 = ifelse((dtSurv$D1 == 1) | ((dtSurv$D <= t_3) & (dtSurv$D <= dtSurv$Y)) , 1, 0)
dtSurv$nci = ifelse((dtSurv$D1 == 1) | (dtSurv$D2 == 1) , 1, 0)

# define missing

defM <- defMiss(varname = "M2", formula = "ifelse (Y < ..t_2, 1, 0)", logit.link = FALSE)
missMat <- genMiss(dtName = dtSurv, missDefs = defM, idvars = "id")
# Generate observed data from actual data and missing data matrix
dtObs <- genObs(dtSurv, missMat, idvars = "id")
defM <- defMiss(varname = "M2", formula = "ifelse(D1 == 1, 1, 0)", logit.link = FALSE)
missMat2 <- genMiss(dtName = dtObs, missDefs = defM, idvars = "id")
# Generate observed data from actual data and missing data matrix
dtObs <- genObs(dtObs, missMat2, idvars = "id")

dt_m = dtObs
dt_m$time_to_event = pmin(dt_m$Y, dt_m$D)
dt_m$eventHappened = ifelse(dt_m$Y < dt_m$D, 1, 0)

dt_m = dt_m[,-"D"]
dt_m = dt_m[,-"Y"]
simData = dt_m

# Restrict follow-up to 20 years
simData$eventHappened[simData$time_to_event > t_3] = 0
simData$time_to_event[simData$time_to_event > t_3] = t_3

simData$C = as.factor(simData$C)
simData$D1 = as.factor(simData$D1)
simData$D2 = as.factor(simData$D2)

filename <- paste0("simData.Y10D10.Population.RData")
save(simData, file = filename)


# Scenario with 10 % death and 35 % event
f <- function(x) {
  exp(-exp(-4.867209182 + D_1_E * E + D_1_C * C + D_1_M_1 * M1) * t_2) + 
    exp(-exp(x + D_2_E * E + D_2_C * C + D_2_M_2 * M2) * t_3) - 
    exp(-exp(x + D_2_E * E + D_2_C * C + D_2_M_2 * M2) * t_2) - 0.9
}
uniroot(f, c(-10, -4))$root

summary(simData$nci)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00000 0.00000 0.00000 0.09848 0.00000 1.00000 

f <- function(x) {
  exp(-exp(-7.45413779 + Y_1_E * E + Y_1_C * C + Y_1_M_1 * M1) * t_2) + 
    exp(-exp(x + D_2_E * E + D_2_C * C + D_2_M_2 * M2) * t_3) - 
    exp(-exp(x + D_2_E * E + D_2_C * C + D_2_M_2 * M2) * t_2) - 0.9
}
uniroot(f, c(-10, -4))$root