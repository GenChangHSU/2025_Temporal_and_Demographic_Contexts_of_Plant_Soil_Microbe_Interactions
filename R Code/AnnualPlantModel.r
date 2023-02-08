########################################################################################################################
#### Code to simulate Beverton-Holt annual plant model with seed bank
#### 2022.02.08 
########################################################################################################################



########################################################################################################################
#### Set working directory and load data files
library("tidyverse")



########################################################################################################################
#### Simulation function
BevertonHolt <- function(frame, parameters){
  
  Nsim = parameters$N
  Gsim = parameters$g
  Ssim = parameters$s
  Lsim = parameters$lambda
  Asim = parameters$A
  
  for(t in 2:dim(frame)[1]){
    x = frame[t-1, 2:(Nsim+1)]
    frame[t, 2:(Nsim+1)] = (1 - Gsim) * Ssim * x + Gsim * x * (Lsim / (1 + Asim %*% (Gsim * x)))
  }
  return(frame)
}



# ############################################################################################################################ Parameter setup
# #### Parameterize model (multiple species)
# # Number of species
# N = 5
# 
# # Seed germination rate
# g = rep(1, N)
# 
# # Survival rate of ungerminated seeds
# s = rep(0, N)
# 
# # Fecundity
# lambda = c(50, 50, 50, 50, 50)
# 
# # Interaction matrix
# A = matrix(c(1.0, 0.4, 0.5, 0.1, 0.5, 
#              0.3, 1.0, 0.2, 0.6, 0.5, 
#              0.5, 0.5, 1.0, 0.4, 0.5, 
#              0.5, 0.5, 0.3, 1.0, 0.5, 
#              0.1, 0.2, 0.3, 0.5, 1.0), N, N, byrow = TRUE)



############################################################################################################################ Parameter setup
#### Parameterize model (The coexist-pair FE-HO from Van Dyke et al., 2022, Nature)
#### N1 = FE (more abundant) & N2 = HO (less abundant)
# Number of species
N <- 2

# Seed germination rate
g <- c(0.752273, 0.666667)

# Survival rate of ungerminated seeds
s <- c(0.133750, 0.045000)

# Fecundity
lambda <- c(2129.949613, 736.667018)

# Interaction matrix
A <- matrix(c(0.588199, 1.410938,
              0.109412, 0.948399), N, N, byrow = TRUE)



########################################################################################################################
#### Simulation properties setup
# Aggregate parameters
parms <- list(N = N, 
              g = g, 
              s = s, 
              lambda = lambda, 
              A = A)

# Run time
Time <- 100

# Output matrix
out <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
colnames(out) <- c("Time", paste("N", c(1:N), sep = ""))

# Initial density
out[1, 2:(N+1)] <- rep(10, N)



########################################################################################################################
#### Run simulation
out <- BevertonHolt(out, parms)



# ########################################################################################################################
# #### Solve analytically (only works when coexistence is possible)
# b <- ((lambda * g)/(1 - (1 - g) * s)) - 1
# equilibrium <- data.frame(Time = rep(T, N), 
#                           Species = paste("N", c(1:N), sep=""), 
#                           Abundance = (solve(A, b) / g))



########################################################################################################################
#### Compile data and plot simulation time series
data <- as.data.frame(out) %>% 
  gather(Species, Abundance, -Time)

ggplot(data, aes(x = Time, y = Abundance, color = Species)) + 
  geom_line() 

# #### Relative abundance
# out %>% 
#   as.data.frame() %>%
#   mutate(RA.N2 = N2 / (N1 + N2)) %>%
#   ggplot(aes(x = Time, y = RA.N2)) + 
#   geom_line()
