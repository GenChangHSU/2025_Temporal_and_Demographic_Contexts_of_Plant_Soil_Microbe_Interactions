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


################################################################################
################################################################################
################################################################################
### Sensitivity analysis of parameters in Beverton-Holt model ##################
set.seed(123)

### A function for performing sensitivity analysis of parameters in the BH model
### Arguments:
### pchange_focal_par: percentage of increase in the focal parameter
### pchange_focal_par: percentage of increase/decrease in the non-focal parameter(s)
### n_generation: number of generation in each simulation
### n_sim: number of simulations for the analysis of each parameter
BH_sensitivity_one_par <- function(pchange_focal_par = 20, 
                                   pchange_nonfocal_par = 5,
                                   n_generation = 100,
                                   n_sim = 100){
  
  ### Beverton-Holt model 
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
  
  ### Default parameters (coexistence)
  # (1) Number of species
  N <- 2
  
  # (2) Seed germination rate
  g <- c(0.752273, 0.666667)
  
  # (3) Survival rate of ungerminated seeds
  s <- c(0.133750, 0.045000)
  
  # (4) Fecundity
  lambda <- c(2129.949613, 736.667018)
  
  # (5) Interaction matrix
  A <- matrix(c(0.588199, 1.410938,
                0.109412, 0.948399), N, N, byrow = TRUE)
  
  parms <- list(N = N, 
                g = g, 
                s = s, 
                lambda = lambda, 
                A = A)
  
  # Number of generations
  Time <- n_generation
  
  # A matrix to store the outputs
  out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
  colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
  
  # Initial abundance
  out_mat[1, 2:(N+1)] <- rep(10, N)
  
  # Model outputs
  out_mat <- BevertonHolt(out_mat, parms)
  
  # Extract the abundance in the last generation
  N2_equl_default <- out_mat[Time, N+1]
  
  
  ################################################################
  ### Sensitivity analysis of a single parameter #################
  ### Fix the parameters of N1 and change the parameters of N2 ###
  
  ### 1. Sensitivity of g
  g_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + pchange_focal_par/100))
    s_temp <- c(s[1], s[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100), 
                       A[2, 2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Extract the abundance in the last generation
    out_mat[Time, N+1]
    
  })

  ### 2. Sensitivity of s
  s_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    s_temp <- c(s[1], s[2]*(1 + pchange_focal_par/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100), 
                       A[2, 2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Extract the abundance in the last generation
    out_mat[Time, N+1]
    
  })
  
  ### 3. Sensitivity of lambda
  lambda_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    s_temp <- c(s[1], s[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + pchange_focal_par/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100), 
                       A[2, 2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Extract the abundance in the last generation
    out_mat[Time, N+1]
    
  })
  
  ### 4. Sensitivity of A_21
  A_21_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    s_temp <- c(s[1], s[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + pchange_focal_par/100), 
                       A[2, 2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Extract the abundance in the last generation
    out_mat[Time, N+1]
    
  })
  
  ### 5. Sensitivity of A_22
  A_22_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    s_temp <- c(s[1], s[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100), 
                       A[2, 2]*(1 + pchange_focal_par/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Extract the abundance in the last generation
    out_mat[Time, N+1]
    
  })
  
  ### Combine all the simulation results
  pars_sens_df <- data.frame(parameter = rep(c("g", "s", "lambda", "A_21", "A_22"), each = n_sim), 
                             N2_eql = c(g_sens, s_sens, lambda_sens, A_21_sens, A_22_sens))
  
  return(list(N2_equl_default = N2_equl_default, 
              pars_sens_df = pars_sens_df))
  
}


### A function for visualizing the simulation results
### Arguments:
### dat: outputs returned by BH_sensitivity_one_par()
### pchange_focal_par: percentage of increase in the focal parameter passed to BH_sensitivity_one_par()
### pchange_nonfocal_par: percentage of increase/decrease in the non-focal parameter(s) passed to BH_sensitivity_one_par()
plot_sim_one_par <- function(dat, pchange_focal_par, pchange_nonfocal_par){
  
  # Arrange the parameters by the average N2 abundance of the simulations
  pars_order <- dat$pars_sens_df %>% 
    group_by(parameter) %>% 
    summarise(mean = mean(N2_eql)) %>% 
    mutate(No = 1:5) %>% 
    arrange(desc(mean))
  
  # Reorder the parameter levels for plotting
  dat$pars_sens_df <- dat$pars_sens_df %>% 
    mutate(parameter = factor(parameter, level = pars_order$parameter, ordered = T))
  
  # A vector of x-axis labels
  x_labels <- c(expression(italic(alpha[21])),
                expression(italic(alpha[22])),
                expression(italic(g)),
                expression(italic(lambda)),
                expression(italic(s)))
  
  # Plot
  ggplot(data = dat$pars_sens_df) + 
    geom_hline(yintercept = dat$N2_equl_default, color = "grey60", linetype = "dashed", size = 1) + 
    geom_point(aes(x = parameter, y = N2_eql), position = position_jitter(width = 0.1), 
               color = "grey50", alpha = 0.75) + 
    stat_summary(aes(x = parameter, y = N2_eql), fun.y = mean,
                 fun.ymin = function(x) mean(x) - sd(x), 
                 fun.ymax = function(x) mean(x) + sd(x), 
                 geom = "pointrange", color = "red") + 
    labs(x = "Focal parameter", y = expression(N[2]), 
         title = glue::glue("Focal parameter: +{pchange_focal_par}% \n Non-focal parameters: ±{pchange_nonfocal_par}%")) +
    scale_x_discrete(labels = x_labels[pars_order$No]) + 
    scale_y_continuous(expand = c(0.25, 0)) + 
    theme_classic() + 
    theme(axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text = element_text(color = "black", size = 14),
          plot.title = element_text(hjust = 0.5, size = 16),
          axis.title.x = element_text(size = 15, margin = margin(t = 10)),
          axis.title.y = element_text(size = 15, margin = margin(r = 8)))
}

### Visualize the simulation results 
BH_sens_20_5 <- BH_sensitivity_one_par(pchange_focal_par = 20, 
                                       pchange_nonfocal_par = 5,
                                       n_generation = 100,
                                       n_sim = 100)
BH_sens_40_5 <- BH_sensitivity_one_par(pchange_focal_par = 40, 
                                       pchange_nonfocal_par = 5,
                                       n_generation = 100,
                                       n_sim = 100)
BH_sens_20_10 <- BH_sensitivity_one_par(pchange_focal_par = 20, 
                                       pchange_nonfocal_par = 10,
                                       n_generation = 100,
                                       n_sim = 100)
BH_sens_40_20 <- BH_sensitivity_one_par(pchange_focal_par = 40, 
                                       pchange_nonfocal_par = 20,
                                       n_generation = 100,
                                       n_sim = 100)

plot_sim_one_par(dat = BH_sens_20_5,
                 pchange_focal_par = 20, 
                 pchange_nonfocal_par = 5)
ggsave("./Outputs/one_par_sens_plot_20_5.tiff", width = 5, height = 4, dpi = 600, device = "tiff")
ggsave("./Outputs/one_par_sens_plot_20_5.pdf", width = 5, height = 4)

plot_sim_one_par(dat = BH_sens_40_5,
                 pchange_focal_par = 40, 
                 pchange_nonfocal_par = 5)
ggsave("./Outputs/one_par_sens_plot_40_5.tiff", width = 5, height = 4, dpi = 600, device = "tiff")

plot_sim_one_par(dat = BH_sens_20_10,
                 pchange_focal_par = 20, 
                 pchange_nonfocal_par = 10)
ggsave("./Outputs/one_par_sens_plot_20_10.tiff", width = 5, height = 4, dpi = 600, device = "tiff")

plot_sim_one_par(dat = BH_sens_40_20,
                 pchange_focal_par = 40, 
                 pchange_nonfocal_par = 20)
ggsave("./Outputs/one_par_sens_plot_40_20.tiff", width = 5, height = 4, dpi = 600, device = "tiff")



################################################################################
################################################################################
################################################################################
### Sensitivity analysis of two parameters

BH_sensitivity_two_pars <- function(pchange_focal_par = 20, 
                                    pchange_nonfocal_par = 5,
                                    n_generation = 100,
                                    n_sim = 100){
  
  ### Beverton-Holt model 
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
  
  ### Default parameters (coexistence)
  # (1) Number of species
  N <- 2
  
  # (2) Seed germination rate
  g <- c(0.752273, 0.666667)
  
  # (3) Survival rate of ungerminated seeds
  s <- c(0.133750, 0.045000)
  
  # (4) Fecundity
  lambda <- c(2129.949613, 736.667018)
  
  # (5) Interaction matrix
  A <- matrix(c(0.588199, 1.410938,
                0.109412, 0.948399), N, N, byrow = TRUE)
  
  parms <- list(N = N, 
                g = g, 
                s = s, 
                lambda = lambda, 
                A = A)
  
  # Number of generations
  Time <- n_generation
  
  # A matrix to store the outputs
  out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
  colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
  
  # Initial abundance
  out_mat[1, 2:(N+1)] <- rep(10, N)
  
  # Model outputs
  out_mat <- BevertonHolt(out_mat, parms)
  
  # Extract the abundance in the last generation
  N2_equl_default <- out_mat[Time, N+1]
  
  
  ################################################################
  ### Sensitivity analysis of two parameters #####################
  ### Fix the parameters of N1 and change the parameters of N2 ###
  
  ### 1. Sensitivity of g & s
  g_s_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + pchange_focal_par/100))
    s_temp <- c(s[1], s[2]*(1 + pchange_focal_par/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100), 
                       A[2, 2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Extract the abundance in the last generation
    out_mat[Time, N+1]
    
  })
  
  ### 2. Sensitivity of g & lambda
  g_lambda_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + pchange_focal_par/100))
    s_temp <- c(s[1], s[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + pchange_focal_par/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100), 
                       A[2, 2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Extract the abundance in the last generation
    out_mat[Time, N+1]
    
  })
  
  ### 3. Sensitivity of g & A_21
  g_A_21_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + pchange_focal_par/100))
    s_temp <- c(s[1], s[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + pchange_focal_par/100), 
                       A[2, 2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Extract the abundance in the last generation
    out_mat[Time, N+1]
    
  })
  
  ### 4. Sensitivity of g & A_22
  g_A_22_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + pchange_focal_par/100))
    s_temp <- c(s[1], s[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100), 
                       A[2, 2]*(1 + pchange_focal_par/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Extract the abundance in the last generation
    out_mat[Time, N+1]
    
  })
  
  ### 5. Sensitivity of s & lambda
  s_lambda_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    s_temp <- c(s[1], s[2]*(1 + pchange_focal_par/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + pchange_focal_par/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100), 
                       A[2, 2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Extract the abundance in the last generation
    out_mat[Time, N+1]
    
  })
  
  ### 6. Sensitivity of s & A_21
  s_A_21_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    s_temp <- c(s[1], s[2]*(1 + pchange_focal_par/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + pchange_focal_par/100), 
                       A[2, 2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Extract the abundance in the last generation
    out_mat[Time, N+1]
    
  })
  
  ### 7. Sensitivity of s & A_22
  s_A_22_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    s_temp <- c(s[1], s[2]*(1 + pchange_focal_par/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100), 
                       A[2, 2]*(1 + pchange_focal_par/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Extract the abundance in the last generation
    out_mat[Time, N+1]
    
  })
  
  ### 8. Sensitivity of lambda & A_21
  lambda_A_21_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    s_temp <- c(s[1], s[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + pchange_focal_par/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + pchange_focal_par/100), 
                       A[2, 2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Extract the abundance in the last generation
    out_mat[Time, N+1]
    
  })
  
  ### 9. Sensitivity of lambda & A_22
  lambda_A_22_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    s_temp <- c(s[1], s[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + pchange_focal_par/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100), 
                       A[2, 2]*(1 + pchange_focal_par/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Extract the abundance in the last generation
    out_mat[Time, N+1]
    
  })
  
  ### 10. Sensitivity of A_21 & A_22
  A_21_A_22_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    s_temp <- c(s[1], s[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + pchange_focal_par/100), 
                       A[2, 2]*(1 + pchange_focal_par/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Extract the abundance in the last generation
    out_mat[Time, N+1]
    
  })  
  
  ### Combine all the simulation results
  pars_sens_df <- data.frame(parameter = rep(c("g_s", "g_lambda", "g_A_21", "g_A_22", "s_lambda",
                                               "s_A_21", "s_A_22", "lambda_A_21", "lambda_A_22", "A_21_A_22"), 
                                             each = n_sim), 
                             N2_eql = c(g_s_sens, g_lambda_sens, g_A_21_sens, g_A_22_sens, s_lambda_sens,
                                        s_A_21_sens, s_A_22_sens, lambda_A_21_sens, lambda_A_22_sens, A_21_A_22_sens))
  
  return(list(N2_equl_default = N2_equl_default, 
              pars_sens_df = pars_sens_df))
  
}

plot_sim_two_pars <- function(dat, pchange_focal_par, pchange_nonfocal_par){
  
  # Arrange the parameters by the average N2 abundance of the simulations
  pars_order <- dat$pars_sens_df %>% 
    group_by(parameter) %>% 
    summarise(mean = mean(N2_eql)) %>% 
    mutate(No = 1:10) %>% 
    arrange(desc(mean))
  
  # Reorder the parameter levels for plotting
  dat$pars_sens_df <- dat$pars_sens_df %>% 
    mutate(parameter = factor(parameter, level = pars_order$parameter, ordered = T))
  
  # A vector of x-axis labels
  x_labels <- c(expression(italic(alpha[21]~+~alpha[22])),
                expression(italic(g~+~alpha[21])),
                expression(italic(g~+~alpha[22])),
                expression(italic(g~+~lambda)),
                expression(italic(g~+~s)),
                expression(italic(lambda~+~alpha[21])),
                expression(italic(lambda~+~alpha[22])),
                expression(italic(s~+~alpha[21])),
                expression(italic(s~+~alpha[22])),
                expression(italic(s~+~lambda)))
  
  # Plot
  ggplot(data = dat$pars_sens_df) + 
    geom_hline(yintercept = dat$N2_equl_default, color = "grey60", linetype = "dashed", size = 1) + 
    geom_point(aes(x = parameter, y = N2_eql), position = position_jitter(width = 0.1), 
               color = "grey50", alpha = 0.75) + 
    stat_summary(aes(x = parameter, y = N2_eql), fun.y = mean,
                 fun.ymin = function(x) mean(x) - sd(x), 
                 fun.ymax = function(x) mean(x) + sd(x), 
                 geom = "pointrange", color = "red") + 
    labs(x = "Focal parameters", y = expression(N[2]), 
         title = glue::glue("Focal parameters: +{pchange_focal_par}% \n Non-focal parameters: ±{pchange_nonfocal_par}%")) +
    scale_x_discrete(labels = x_labels[pars_order$No]) + 
    scale_y_continuous(expand = c(0.25, 0)) + 
    theme_classic() + 
    theme(axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text = element_text(color = "black", size = 14),
          plot.title = element_text(hjust = 0.5, size = 16),
          axis.title.x = element_text(size = 15, margin = margin(t = 10)),
          axis.title.y = element_text(size = 15, margin = margin(r = 8)))
}

### Visualize the simulation results 
BH_sens_20_5 <- BH_sensitivity_two_pars(pchange_focal_par = 20, 
                                        pchange_nonfocal_par = 5,
                                        n_generation = 100,
                                        n_sim = 100)
BH_sens_40_5 <- BH_sensitivity_two_pars(pchange_focal_par = 40, 
                                        pchange_nonfocal_par = 5,
                                        n_generation = 100,
                                        n_sim = 100)
BH_sens_20_10 <- BH_sensitivity_two_pars(pchange_focal_par = 20, 
                                         pchange_nonfocal_par = 10,
                                         n_generation = 100,
                                         n_sim = 100)
BH_sens_40_20 <- BH_sensitivity_two_pars(pchange_focal_par = 40, 
                                         pchange_nonfocal_par = 20,
                                         n_generation = 100,
                                         n_sim = 100)

plot_sim_two_pars(dat = BH_sens_20_5,
                  pchange_focal_par = 20, 
                  pchange_nonfocal_par = 5)
ggsave("./Outputs/two_pars_sens_plot_20_5.tiff", width = 10, height = 5, dpi = 600, device = "tiff")
ggsave("./Outputs/two_pars_sens_plot_20_5.pdf", width = 9, height = 5.5)

plot_sim_two_pars(dat = BH_sens_40_5,
                  pchange_focal_par = 40, 
                  pchange_nonfocal_par = 5)
ggsave("./Outputs/two_pars_sens_plot_40_5.tiff", width = 10, height = 5, dpi = 600, device = "tiff")

plot_sim_two_pars(dat = BH_sens_20_10,
                  pchange_focal_par = 20, 
                  pchange_nonfocal_par = 10)
ggsave("./Outputs/two_pars_sens_plot_20_10.tiff", width = 10, height = 5, dpi = 600, device = "tiff")

plot_sim_two_pars(dat = BH_sens_40_20,
                  pchange_focal_par = 40, 
                  pchange_nonfocal_par = 20)
ggsave("./Outputs/two_pars_sens_plot_40_20.tiff", width = 10, height = 5, dpi = 600, device = "tiff")



################################################################################
################################################################################
################################################################################
### Sensitivity analysis of three parameters

BH_sensitivity_three_pars <- function(pchange_focal_par = 20, 
                                      pchange_nonfocal_par = 5,
                                      n_generation = 100,
                                      n_sim = 100){
  
  ### Beverton-Holt model 
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
  
  ### Default parameters (coexistence)
  # (1) Number of species
  N <- 2
  
  # (2) Seed germination rate
  g <- c(0.752273, 0.666667)
  
  # (3) Survival rate of ungerminated seeds
  s <- c(0.133750, 0.045000)
  
  # (4) Fecundity
  lambda <- c(2129.949613, 736.667018)
  
  # (5) Interaction matrix
  A <- matrix(c(0.588199, 1.410938,
                0.109412, 0.948399), N, N, byrow = TRUE)
  
  parms <- list(N = N, 
                g = g, 
                s = s, 
                lambda = lambda, 
                A = A)
  
  # Number of generations
  Time <- n_generation
  
  # A matrix to store the outputs
  out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
  colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
  
  # Initial abundance
  out_mat[1, 2:(N+1)] <- rep(10, N)
  
  # Model outputs
  out_mat <- BevertonHolt(out_mat, parms)
  
  # Extract the abundance in the last generation
  N2_equl_default <- out_mat[Time, N+1]
  
  
  ################################################################
  ### Sensitivity analysis of three parameters ###################
  ### Fix the parameters of N1 and change the parameters of N2 ###
  
  ### 1. Sensitivity of g & s & lambda
  g_s_lambda_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + pchange_focal_par/100))
    s_temp <- c(s[1], s[2]*(1 + pchange_focal_par/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + pchange_focal_par/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100), 
                       A[2, 2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Extract the abundance in the last generation
    out_mat[Time, N+1]
    
  })
  
  ### 2. Sensitivity of g & s & A_21
  g_s_A_21_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + pchange_focal_par/100))
    s_temp <- c(s[1], s[2]*(1 + pchange_focal_par/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + pchange_focal_par/100), 
                       A[2, 2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Extract the abundance in the last generation
    out_mat[Time, N+1]
    
  })
  
  ### 3. Sensitivity of g & s & A_22
  g_s_A_22_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + pchange_focal_par/100))
    s_temp <- c(s[1], s[2]*(1 + pchange_focal_par/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100), 
                       A[2, 2]*(1 + pchange_focal_par/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Extract the abundance in the last generation
    out_mat[Time, N+1]
    
  })
  
  ### 4. Sensitivity of g & lambda & A_21
  g_lambda_A_21_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + pchange_focal_par/100))
    s_temp <- c(s[1], s[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + pchange_focal_par/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + pchange_focal_par/100), 
                       A[2, 2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Extract the abundance in the last generation
    out_mat[Time, N+1]
    
  })
  
  ### 5. Sensitivity of g & lambda & A_22
  g_lambda_A_22_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + pchange_focal_par/100))
    s_temp <- c(s[1], s[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + pchange_focal_par/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100), 
                       A[2, 2]*(1 + pchange_focal_par/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Extract the abundance in the last generation
    out_mat[Time, N+1]
    
  })
  
  ### 6. Sensitivity of g & A_21 & A_22
  g_A_21_A_22_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + pchange_focal_par/100))
    s_temp <- c(s[1], s[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + pchange_focal_par/100), 
                       A[2, 2]*(1 + pchange_focal_par/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Extract the abundance in the last generation
    out_mat[Time, N+1]
    
  })
  
  ### 7. Sensitivity of s & lambda & A_21
  s_lambda_A_21_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    s_temp <- c(s[1], s[2]*(1 + pchange_focal_par/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + pchange_focal_par/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + pchange_focal_par/100), 
                       A[2, 2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Extract the abundance in the last generation
    out_mat[Time, N+1]
    
  })
  
  ### 8. Sensitivity of s & lambda & A_22
  s_lambda_A_22_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    s_temp <- c(s[1], s[2]*(1 + pchange_focal_par/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + pchange_focal_par/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100), 
                       A[2, 2]*(1 + pchange_focal_par/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Extract the abundance in the last generation
    out_mat[Time, N+1]
    
  })
  
  ### 9. Sensitivity of s & A_21 & A_22
  s_A_21_A_22_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    s_temp <- c(s[1], s[2]*(1 + pchange_focal_par/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + pchange_focal_par/100), 
                       A[2, 2]*(1 + pchange_focal_par/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Extract the abundance in the last generation
    out_mat[Time, N+1]
    
  })
  
  ### 10. Sensitivity of lambda & A_21 & A_22
  lambda_A_21_A_22_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    s_temp <- c(s[1], s[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + pchange_focal_par/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + pchange_focal_par/100), 
                       A[2, 2]*(1 + pchange_focal_par/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Extract the abundance in the last generation
    out_mat[Time, N+1]
    
  })  
  
  ### Combine all the simulation results
  pars_sens_df <- data.frame(parameter = rep(c("g_s_lambda", "g_s_A_21", "g_s_A_22", "g_lambda_A_21", "g_lambda_A_22",
                                               "g_A_21_A_22", "s_lambda_A_21", "s_lambda_A_22", "s_A_21_A_22", "lambda_A_21_A_22"), 
                                             each = n_sim), 
                             N2_eql = c(g_s_lambda_sens, g_s_A_21_sens, g_s_A_22_sens, g_lambda_A_21_sens, g_lambda_A_22_sens,
                                        g_A_21_A_22_sens, s_lambda_A_21_sens, s_lambda_A_22_sens, s_A_21_A_22_sens, lambda_A_21_A_22_sens))
  
  return(list(N2_equl_default = N2_equl_default, 
              pars_sens_df = pars_sens_df))
  
}

plot_sim_three_pars <- function(dat, pchange_focal_par, pchange_nonfocal_par){
  
  # Arrange the parameters by the average N2 abundance of the simulations
  pars_order <- dat$pars_sens_df %>% 
    group_by(parameter) %>% 
    summarise(mean = mean(N2_eql)) %>% 
    mutate(No = 1:10) %>% 
    arrange(desc(mean))
  
  # Reorder the parameter levels for plotting
  dat$pars_sens_df <- dat$pars_sens_df %>% 
    mutate(parameter = factor(parameter, level = pars_order$parameter, ordered = T))
  
  # A vector of x-axis labels
  x_labels <- c(expression(italic(g~+~alpha[21]~+~alpha[22])),
                expression(italic(g~+~lambda~+~alpha[21])),
                expression(italic(g~+~lambda~+~alpha[22])),
                expression(italic(g~+~s~+~alpha[21])),
                expression(italic(g~+~s~+~alpha[22])),
                expression(italic(g~+~s~+~lambda)),
                expression(italic(lambda~+~alpha[21]~+~alpha[22])),
                expression(italic(s~+~alpha[21]~+~alpha[22])),
                expression(italic(s~+~lambda~+~alpha[21])),
                expression(italic(s~+~lambda~+~alpha[22])))
  
  # Plot
  ggplot(data = dat$pars_sens_df) + 
    geom_hline(yintercept = dat$N2_equl_default, color = "grey60", linetype = "dashed", size = 1) + 
    geom_point(aes(x = parameter, y = N2_eql), position = position_jitter(width = 0.1), 
               color = "grey50", alpha = 0.75) + 
    stat_summary(aes(x = parameter, y = N2_eql), fun.y = mean,
                 fun.ymin = function(x) mean(x) - sd(x), 
                 fun.ymax = function(x) mean(x) + sd(x), 
                 geom = "pointrange", color = "red") + 
    labs(x = "Focal parameters", y = expression(N[2]), 
         title = glue::glue("Focal parameters: +{pchange_focal_par}% \n Non-focal parameters: ±{pchange_nonfocal_par}%")) +
    scale_x_discrete(labels = x_labels[pars_order$No]) + 
    scale_y_continuous(expand = c(0.25, 0)) + 
    theme_classic() + 
    theme(axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(color = "black", size = 14, angle = 45),
          axis.text.y = element_text(color = "black", size = 14),
          plot.title = element_text(hjust = 0.5, size = 16),
          axis.title.x = element_text(size = 15, margin = margin(t = 0)),
          axis.title.y = element_text(size = 15, margin = margin(r = 8)))
}

### Visualize the simulation results 
BH_sens_20_5 <- BH_sensitivity_three_pars(pchange_focal_par = 20, 
                                          pchange_nonfocal_par = 5,
                                          n_generation = 100,
                                          n_sim = 100)
BH_sens_40_5 <- BH_sensitivity_three_pars(pchange_focal_par = 40, 
                                          pchange_nonfocal_par = 5,
                                          n_generation = 100,
                                          n_sim = 100)
BH_sens_20_10 <- BH_sensitivity_three_pars(pchange_focal_par = 20, 
                                           pchange_nonfocal_par = 10,
                                           n_generation = 100,
                                           n_sim = 100)
BH_sens_40_20 <- BH_sensitivity_three_pars(pchange_focal_par = 40, 
                                           pchange_nonfocal_par = 20,
                                           n_generation = 100,
                                           n_sim = 100)

plot_sim_three_pars(dat = BH_sens_20_5,
                    pchange_focal_par = 20, 
                    pchange_nonfocal_par = 5)
ggsave("./Outputs/three_pars_sens_plot_20_5.tiff", width = 10, height = 5, dpi = 600, device = "tiff")
ggsave("./Outputs/three_pars_sens_plot_20_5.pdf", width = 9, height = 5.5)

plot_sim_three_pars(dat = BH_sens_40_5,
                    pchange_focal_par = 40, 
                    pchange_nonfocal_par = 5)
ggsave("./Outputs/three_pars_sens_plot_40_5.tiff", width = 10, height = 5, dpi = 600, device = "tiff")

plot_sim_three_pars(dat = BH_sens_20_10,
                    pchange_focal_par = 20, 
                    pchange_nonfocal_par = 10)
ggsave("./Outputs/three_pars_sens_plot_20_10.tiff", width = 10, height = 5, dpi = 600, device = "tiff")

plot_sim_three_pars(dat = BH_sens_40_20,
                    pchange_focal_par = 40, 
                    pchange_nonfocal_par = 20)
ggsave("./Outputs/three_pars_sens_plot_40_20.tiff", width = 10, height = 5, dpi = 600, device = "tiff")



################################################################################
################################################################################
################################################################################
### Sensitivity analysis of four parameters

BH_sensitivity_four_pars <- function(pchange_focal_par = 20, 
                                     pchange_nonfocal_par = 5,
                                     n_generation = 100,
                                     n_sim = 100){
  
  ### Beverton-Holt model 
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
  
  ### Default parameters (coexistence)
  # (1) Number of species
  N <- 2
  
  # (2) Seed germination rate
  g <- c(0.752273, 0.666667)
  
  # (3) Survival rate of ungerminated seeds
  s <- c(0.133750, 0.045000)
  
  # (4) Fecundity
  lambda <- c(2129.949613, 736.667018)
  
  # (5) Interaction matrix
  A <- matrix(c(0.588199, 1.410938,
                0.109412, 0.948399), N, N, byrow = TRUE)
  
  parms <- list(N = N, 
                g = g, 
                s = s, 
                lambda = lambda, 
                A = A)
  
  # Number of generations
  Time <- n_generation
  
  # A matrix to store the outputs
  out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
  colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
  
  # Initial abundance
  out_mat[1, 2:(N+1)] <- rep(10, N)
  
  # Model outputs
  out_mat <- BevertonHolt(out_mat, parms)
  
  # Extract the abundance in the last generation
  N2_equl_default <- out_mat[Time, N+1]
  
  
  ################################################################
  ### Sensitivity analysis of four parameters ####################
  ### Fix the parameters of N1 and change the parameters of N2 ###
  
  ### 1. Sensitivity of g & s & lambda & A_21
  g_s_lambda_A_21_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + pchange_focal_par/100))
    s_temp <- c(s[1], s[2]*(1 + pchange_focal_par/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + pchange_focal_par/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + pchange_focal_par/100), 
                       A[2, 2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Extract the abundance in the last generation
    out_mat[Time, N+1]
    
  })
  
  ### 2. Sensitivity of g & s & A_21 & A_22
  g_s_A_21_A_22_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + pchange_focal_par/100))
    s_temp <- c(s[1], s[2]*(1 + pchange_focal_par/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + pchange_focal_par/100), 
                       A[2, 2]*(1 + pchange_focal_par/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Extract the abundance in the last generation
    out_mat[Time, N+1]
    
  })
  
  ### 3. Sensitivity of g & s & lambda & A_22
  g_s_lambda_A_22_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + pchange_focal_par/100))
    s_temp <- c(s[1], s[2]*(1 + pchange_focal_par/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + pchange_focal_par/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100), 
                       A[2, 2]*(1 + pchange_focal_par/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Extract the abundance in the last generation
    out_mat[Time, N+1]
    
  })
  
  ### 4. Sensitivity of g & lambda & A_21 & A_22
  g_lambda_A_21_A_22_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + pchange_focal_par/100))
    s_temp <- c(s[1], s[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + pchange_focal_par/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + pchange_focal_par/100), 
                       A[2, 2]*(1 + pchange_focal_par/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Extract the abundance in the last generation
    out_mat[Time, N+1]
    
  })
  
  ### 5. Sensitivity of s & lambda & A_21 & A_22
  s_lambda_A_21_A_22_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    s_temp <- c(s[1], s[2]*(1 + pchange_focal_par/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + pchange_focal_par/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + pchange_focal_par/100), 
                       A[2, 2]*(1 + pchange_focal_par/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Extract the abundance in the last generation
    out_mat[Time, N+1]
    
  })
  
  ### Combine all the simulation results
  pars_sens_df <- data.frame(parameter = rep(c("g_s_lambda_A_21", "g_s_A_21_A_22", "g_s_lambda_A_22", "g_lambda_A_21_A_22", "s_lambda_A_21_A_22"), each = n_sim), 
                             N2_eql = c(g_s_lambda_A_21_sens, g_s_A_21_A_22_sens, g_s_lambda_A_22_sens, g_lambda_A_21_A_22_sens, s_lambda_A_21_A_22_sens))
  
  return(list(N2_equl_default = N2_equl_default, 
              pars_sens_df = pars_sens_df))
  
}

plot_sim_four_pars <- function(dat, pchange_focal_par, pchange_nonfocal_par){
  
  # Arrange the parameters by the average N2 abundance of the simulations
  pars_order <- dat$pars_sens_df %>% 
    group_by(parameter) %>% 
    summarise(mean = mean(N2_eql)) %>% 
    mutate(No = 1:5) %>% 
    arrange(desc(mean))
  
  # Reorder the parameter levels for plotting
  dat$pars_sens_df <- dat$pars_sens_df %>% 
    mutate(parameter = factor(parameter, level = pars_order$parameter, ordered = T))
  
  # A vector of x-axis labels
  x_labels <- c(expression(italic(g~+~lambda~+~alpha[21]~+~alpha[22])),
                expression(italic(g~+~s~+~alpha[21]~+~alpha[22])),
                expression(italic(g~+~s~+~lambda~+~alpha[21])),
                expression(italic(g~+~s~+~lambda~+~alpha[22])),
                expression(italic(s~+~lambda~+~alpha[21]~+~alpha[22])))
  
  # Plot
  ggplot(data = dat$pars_sens_df) + 
    geom_hline(yintercept = dat$N2_equl_default, color = "grey60", linetype = "dashed", size = 1) + 
    geom_point(aes(x = parameter, y = N2_eql), position = position_jitter(width = 0.1), 
               color = "grey50", alpha = 0.75) + 
    stat_summary(aes(x = parameter, y = N2_eql), fun.y = mean,
                 fun.ymin = function(x) mean(x) - sd(x), 
                 fun.ymax = function(x) mean(x) + sd(x), 
                 geom = "pointrange", color = "red") + 
    labs(x = "Focal parameters", y = expression(N[2]), 
         title = glue::glue("Focal parameters: +{pchange_focal_par}% \n Non-focal parameter: ±{pchange_nonfocal_par}%")) +
    scale_x_discrete(labels = x_labels[pars_order$No]) + 
    scale_y_continuous(expand = c(0.25, 0)) + 
    theme_classic() + 
    theme(axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(color = "black", size = 12, angle = 45),
          axis.text.y = element_text(color = "black", size = 14),
          plot.title = element_text(hjust = 0.5, size = 16),
          axis.title.x = element_text(size = 15, margin = margin(t = -5)),
          axis.title.y = element_text(size = 15, margin = margin(r = 8)))
}

### Visualize the simulation results 
BH_sens_20_5 <- BH_sensitivity_four_pars(pchange_focal_par = 20, 
                                         pchange_nonfocal_par = 5,
                                         n_generation = 100,
                                         n_sim = 100)
BH_sens_40_5 <- BH_sensitivity_four_pars(pchange_focal_par = 40, 
                                         pchange_nonfocal_par = 5,
                                         n_generation = 100,
                                         n_sim = 100)
BH_sens_20_10 <- BH_sensitivity_four_pars(pchange_focal_par = 20, 
                                          pchange_nonfocal_par = 10,
                                          n_generation = 100,
                                          n_sim = 100)
BH_sens_40_20 <- BH_sensitivity_four_pars(pchange_focal_par = 40, 
                                          pchange_nonfocal_par = 20,
                                          n_generation = 100,
                                          n_sim = 100)

plot_sim_four_pars(dat = BH_sens_20_5,
                   pchange_focal_par = 20, 
                   pchange_nonfocal_par = 5)
ggsave("./Outputs/four_pars_sens_plot_20_5.tiff", width = 7, height = 5, dpi = 600, device = "tiff")
ggsave("./Outputs/four_pars_sens_plot_20_5.pdf", width = 7, height = 5)

plot_sim_four_pars(dat = BH_sens_40_5,
                   pchange_focal_par = 40, 
                   pchange_nonfocal_par = 5)
ggsave("./Outputs/four_pars_sens_plot_40_5.tiff", width = 7, height = 5, dpi = 600, device = "tiff")

plot_sim_four_pars(dat = BH_sens_20_10,
                   pchange_focal_par = 20, 
                   pchange_nonfocal_par = 10)
ggsave("./Outputs/four_pars_sens_plot_20_10.tiff", width = 7, height = 5, dpi = 600, device = "tiff")

plot_sim_four_pars(dat = BH_sens_40_20,
                   pchange_focal_par = 40, 
                   pchange_nonfocal_par = 20)
ggsave("./Outputs/four_pars_sens_plot_40_20.tiff", width = 7, height = 5, dpi = 600, device = "tiff")



################################################################################
################################################################################
################################################################################
### Sensitivity analysis using the relative abundance of N2 as the response ####
### One parameter
BH_sensitivity_one_par_relabd <- function(pchange_focal_par = 20, 
                                   pchange_nonfocal_par = 5,
                                   n_generation = 100,
                                   n_sim = 100){
  
  ### Beverton-Holt model 
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
  
  ### Default parameters (coexistence)
  # (1) Number of species
  N <- 2
  
  # (2) Seed germination rate
  g <- c(0.752273, 0.666667)
  
  # (3) Survival rate of ungerminated seeds
  s <- c(0.133750, 0.045000)
  
  # (4) Fecundity
  lambda <- c(2129.949613, 736.667018)
  
  # (5) Interaction matrix
  A <- matrix(c(0.588199, 1.410938,
                0.109412, 0.948399), N, N, byrow = TRUE)
  
  parms <- list(N = N, 
                g = g, 
                s = s, 
                lambda = lambda, 
                A = A)
  
  # Number of generations
  Time <- n_generation
  
  # A matrix to store the outputs
  out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
  colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
  
  # Initial abundance
  out_mat[1, 2:(N+1)] <- rep(10, N)
  
  # Model outputs
  out_mat <- BevertonHolt(out_mat, parms)
  
  # Relative abundance of N2 in the last generation
  N2_equl_default <- out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
  
  
  ################################################################
  ### Sensitivity analysis of a single parameter #################
  ### Fix the parameters of N1 and change the parameters of N2 ###
  
  ### 1. Sensitivity of g
  g_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + pchange_focal_par/100))
    s_temp <- c(s[1], s[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100), 
                       A[2, 2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Relative abundance of N2 in the last generation
    out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
    
  })
  
  ### 2. Sensitivity of s
  s_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    s_temp <- c(s[1], s[2]*(1 + pchange_focal_par/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100), 
                       A[2, 2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Relative abundance of N2 in the last generation
    out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
    
  })
  
  ### 3. Sensitivity of lambda
  lambda_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    s_temp <- c(s[1], s[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + pchange_focal_par/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100), 
                       A[2, 2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Relative abundance of N2 in the last generation
    out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
    
  })
  
  ### 4. Sensitivity of A_21
  A_21_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    s_temp <- c(s[1], s[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + pchange_focal_par/100), 
                       A[2, 2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Relative abundance of N2 in the last generation
    out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
    
  })
  
  ### 5. Sensitivity of A_22
  A_22_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    s_temp <- c(s[1], s[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100), 
                       A[2, 2]*(1 + pchange_focal_par/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Relative abundance of N2 in the last generation
    out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
    
  })
  
  ### Combine all the simulation results
  pars_sens_df <- data.frame(parameter = rep(c("g", "s", "lambda", "A_21", "A_22"), each = n_sim), 
                             N2_eql = c(g_sens, s_sens, lambda_sens, A_21_sens, A_22_sens))
  
  return(list(N2_equl_default = N2_equl_default, 
              pars_sens_df = pars_sens_df))
  
}

plot_sim_one_par_relabd <- function(dat, pchange_focal_par, pchange_nonfocal_par){
  
  # Arrange the parameters by the average N2 abundance of the simulations
  pars_order <- dat$pars_sens_df %>% 
    group_by(parameter) %>% 
    summarise(mean = mean(N2_eql)) %>% 
    mutate(No = 1:5) %>% 
    arrange(desc(mean))
  
  # Reorder the parameter levels for plotting
  dat$pars_sens_df <- dat$pars_sens_df %>% 
    mutate(parameter = factor(parameter, level = pars_order$parameter, ordered = T))
  
  # A vector of x-axis labels
  x_labels <- c(expression(italic(alpha[21])),
                expression(italic(alpha[22])),
                expression(italic(g)),
                expression(italic(lambda)),
                expression(italic(s)))
  
  y_line <- ifelse(max(dat$pars_sens_df$N2_eql > 0.99), 1, 0)
  y_lim <- ifelse(max(dat$pars_sens_df$N2_eql > 0.99), 1.1, NA)
  
  # Plot
  ggplot(data = dat$pars_sens_df) + 
    geom_hline(yintercept = dat$N2_equl_default, color = "grey60", linetype = "dashed", size = 1) + 
    geom_hline(yintercept = 0, color = "blue", linetype = "dashed", size = 1) + 
    geom_hline(yintercept = y_line, color = "blue", linetype = "dashed", size = 1) + 
    geom_point(aes(x = parameter, y = N2_eql), position = position_jitter(width = 0.1), 
               color = "grey50", alpha = 0.75) + 
    stat_summary(aes(x = parameter, y = N2_eql), fun.y = mean,
                 fun.ymin = function(x) mean(x) - sd(x), 
                 fun.ymax = function(x) mean(x) + sd(x), 
                 geom = "pointrange", color = "red") + 
    labs(x = "Focal parameter", y = expression(Relative~abundance~of~N[2]), 
         title = glue::glue("Focal parameter: +{pchange_focal_par}% \n Non-focal parameters: ±{pchange_nonfocal_par}%")) +
    scale_x_discrete(labels = x_labels[pars_order$No]) + 
    scale_y_continuous(expand = c(0.25, 0), limits = c(0, y_lim)) + 
    theme_classic() + 
    theme(axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(color = "black", size = 12),
          axis.text.y = element_text(color = "black", size = 14),
          plot.title = element_text(hjust = 0.5, size = 16),
          axis.title.x = element_text(size = 15, margin = margin(t = 10)),
          axis.title.y = element_text(size = 15, margin = margin(r = 8)))
}

### Visualize the simulation results 
BH_sens_20_5 <- BH_sensitivity_one_par_relabd(pchange_focal_par = 20, 
                                       pchange_nonfocal_par = 5,
                                       n_generation = 100,
                                       n_sim = 100)
BH_sens_40_5 <- BH_sensitivity_one_par_relabd(pchange_focal_par = 40, 
                                       pchange_nonfocal_par = 5,
                                       n_generation = 100,
                                       n_sim = 100)
BH_sens_20_10 <- BH_sensitivity_one_par_relabd(pchange_focal_par = 20, 
                                        pchange_nonfocal_par = 10,
                                        n_generation = 100,
                                        n_sim = 100)
BH_sens_40_20 <- BH_sensitivity_one_par_relabd(pchange_focal_par = 40, 
                                        pchange_nonfocal_par = 20,
                                        n_generation = 100,
                                        n_sim = 100)

plot_sim_one_par_relabd(dat = BH_sens_20_5,
                 pchange_focal_par = 20, 
                 pchange_nonfocal_par = 5)
ggsave("./Outputs/one_par_sens_plot_20_5_relabd.tiff", width = 5, height = 4, dpi = 600, device = "tiff")
ggsave("./Outputs/one_par_sens_plot_20_5_relabd.pdf", width = 5, height = 4)

plot_sim_one_par_relabd(dat = BH_sens_40_5,
                 pchange_focal_par = 40, 
                 pchange_nonfocal_par = 5)
ggsave("./Outputs/one_par_sens_plot_40_5_relabd.tiff", width = 5, height = 4, dpi = 600, device = "tiff")

plot_sim_one_par_relabd(dat = BH_sens_20_10,
                 pchange_focal_par = 20, 
                 pchange_nonfocal_par = 10)
ggsave("./Outputs/one_par_sens_plot_20_10_relabd.tiff", width = 5, height = 4, dpi = 600, device = "tiff")

plot_sim_one_par_relabd(dat = BH_sens_40_20,
                 pchange_focal_par = 40, 
                 pchange_nonfocal_par = 20)
ggsave("./Outputs/one_par_sens_plot_40_20_relabd.tiff", width = 5, height = 4, dpi = 600, device = "tiff")


### Two parameters
BH_sensitivity_two_pars_relabd <- function(pchange_focal_par = 20, 
                                    pchange_nonfocal_par = 5,
                                    n_generation = 100,
                                    n_sim = 100){
  
  ### Beverton-Holt model 
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
  
  ### Default parameters (coexistence)
  # (1) Number of species
  N <- 2
  
  # (2) Seed germination rate
  g <- c(0.752273, 0.666667)
  
  # (3) Survival rate of ungerminated seeds
  s <- c(0.133750, 0.045000)
  
  # (4) Fecundity
  lambda <- c(2129.949613, 736.667018)
  
  # (5) Interaction matrix
  A <- matrix(c(0.588199, 1.410938,
                0.109412, 0.948399), N, N, byrow = TRUE)
  
  parms <- list(N = N, 
                g = g, 
                s = s, 
                lambda = lambda, 
                A = A)
  
  # Number of generations
  Time <- n_generation
  
  # A matrix to store the outputs
  out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
  colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
  
  # Initial abundance
  out_mat[1, 2:(N+1)] <- rep(10, N)
  
  # Model outputs
  out_mat <- BevertonHolt(out_mat, parms)
  
  # Relative abundance of N2 in the last generation
  N2_equl_default <- out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
  
  
  ################################################################
  ### Sensitivity analysis of two parameters #####################
  ### Fix the parameters of N1 and change the parameters of N2 ###
  
  ### 1. Sensitivity of g & s
  g_s_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + pchange_focal_par/100))
    s_temp <- c(s[1], s[2]*(1 + pchange_focal_par/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100), 
                       A[2, 2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Relative abundance of N2 in the last generation
    out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
    
  })
  
  ### 2. Sensitivity of g & lambda
  g_lambda_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + pchange_focal_par/100))
    s_temp <- c(s[1], s[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + pchange_focal_par/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100), 
                       A[2, 2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Relative abundance of N2 in the last generation
    out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
    
  })
  
  ### 3. Sensitivity of g & A_21
  g_A_21_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + pchange_focal_par/100))
    s_temp <- c(s[1], s[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + pchange_focal_par/100), 
                       A[2, 2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Relative abundance of N2 in the last generation
    out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
    
  })
  
  ### 4. Sensitivity of g & A_22
  g_A_22_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + pchange_focal_par/100))
    s_temp <- c(s[1], s[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100), 
                       A[2, 2]*(1 + pchange_focal_par/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Relative abundance of N2 in the last generation
    out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
    
  })
  
  ### 5. Sensitivity of s & lambda
  s_lambda_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    s_temp <- c(s[1], s[2]*(1 + pchange_focal_par/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + pchange_focal_par/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100), 
                       A[2, 2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Relative abundance of N2 in the last generation
    out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
    
  })
  
  ### 6. Sensitivity of s & A_21
  s_A_21_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    s_temp <- c(s[1], s[2]*(1 + pchange_focal_par/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + pchange_focal_par/100), 
                       A[2, 2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Relative abundance of N2 in the last generation
    out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
    
  })
  
  ### 7. Sensitivity of s & A_22
  s_A_22_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    s_temp <- c(s[1], s[2]*(1 + pchange_focal_par/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100), 
                       A[2, 2]*(1 + pchange_focal_par/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Relative abundance of N2 in the last generation
    out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
    
  })
  
  ### 8. Sensitivity of lambda & A_21
  lambda_A_21_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    s_temp <- c(s[1], s[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + pchange_focal_par/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + pchange_focal_par/100), 
                       A[2, 2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Relative abundance of N2 in the last generation
    out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
    
  })
  
  ### 9. Sensitivity of lambda & A_22
  lambda_A_22_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    s_temp <- c(s[1], s[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + pchange_focal_par/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100), 
                       A[2, 2]*(1 + pchange_focal_par/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Relative abundance of N2 in the last generation
    out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
    
  })
  
  ### 10. Sensitivity of A_21 & A_22
  A_21_A_22_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    s_temp <- c(s[1], s[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + pchange_focal_par/100), 
                       A[2, 2]*(1 + pchange_focal_par/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Relative abundance of N2 in the last generation
    out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
    
  })  
  
  ### Combine all the simulation results
  pars_sens_df <- data.frame(parameter = rep(c("g_s", "g_lambda", "g_A_21", "g_A_22", "s_lambda",
                                               "s_A_21", "s_A_22", "lambda_A_21", "lambda_A_22", "A_21_A_22"), 
                                             each = n_sim), 
                             N2_eql = c(g_s_sens, g_lambda_sens, g_A_21_sens, g_A_22_sens, s_lambda_sens,
                                        s_A_21_sens, s_A_22_sens, lambda_A_21_sens, lambda_A_22_sens, A_21_A_22_sens))
  
  return(list(N2_equl_default = N2_equl_default, 
              pars_sens_df = pars_sens_df))
  
}

plot_sim_two_pars_relabd <- function(dat, pchange_focal_par, pchange_nonfocal_par){
  
  # Arrange the parameters by the average N2 abundance of the simulations
  pars_order <- dat$pars_sens_df %>% 
    group_by(parameter) %>% 
    summarise(mean = mean(N2_eql)) %>% 
    mutate(No = 1:10) %>% 
    arrange(desc(mean))
  
  # Reorder the parameter levels for plotting
  dat$pars_sens_df <- dat$pars_sens_df %>% 
    mutate(parameter = factor(parameter, level = pars_order$parameter, ordered = T))
  
  # A vector of x-axis labels
  x_labels <- c(expression(italic(alpha[21]~+~alpha[22])),
                expression(italic(g~+~alpha[21])),
                expression(italic(g~+~alpha[22])),
                expression(italic(g~+~lambda)),
                expression(italic(g~+~s)),
                expression(italic(lambda~+~alpha[21])),
                expression(italic(lambda~+~alpha[22])),
                expression(italic(s~+~alpha[21])),
                expression(italic(s~+~alpha[22])),
                expression(italic(s~+~lambda)))
  
  y_line <- ifelse(max(dat$pars_sens_df$N2_eql > 0.99), 1, 0)
  y_lim <- ifelse(max(dat$pars_sens_df$N2_eql > 0.99), 1.1, NA)
  
  # Plot
  ggplot(data = dat$pars_sens_df) + 
    geom_hline(yintercept = dat$N2_equl_default, color = "grey60", linetype = "dashed", size = 1) + 
    geom_hline(yintercept = 0, color = "blue", linetype = "dashed", size = 1) + 
    geom_hline(yintercept = y_line, color = "blue", linetype = "dashed", size = 1) + 
    geom_point(aes(x = parameter, y = N2_eql), position = position_jitter(width = 0.1), 
               color = "grey50", alpha = 0.75) + 
    stat_summary(aes(x = parameter, y = N2_eql), fun.y = mean,
                 fun.ymin = function(x) mean(x) - sd(x), 
                 fun.ymax = function(x) mean(x) + sd(x), 
                 geom = "pointrange", color = "red") + 
    labs(x = "Focal parameter", y = expression(Relative~abundance~of~N[2]), 
         title = glue::glue("Focal parameters: +{pchange_focal_par}% \n Non-focal parameters: ±{pchange_nonfocal_par}%")) +
    scale_x_discrete(labels = x_labels[pars_order$No]) + 
    scale_y_continuous(expand = c(0.25, 0), limits = c(0, y_lim)) + 
    theme_classic() + 
    theme(axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(color = "black", size = 12),
          axis.text.y = element_text(color = "black", size = 14),
          plot.title = element_text(hjust = 0.5, size = 16),
          axis.title.x = element_text(size = 15, margin = margin(t = 10)),
          axis.title.y = element_text(size = 15, margin = margin(r = 8)))
}

### Visualize the simulation results 
BH_sens_20_5 <- BH_sensitivity_two_pars_relabd(pchange_focal_par = 20, 
                                        pchange_nonfocal_par = 5,
                                        n_generation = 100,
                                        n_sim = 100)
BH_sens_40_5 <- BH_sensitivity_two_pars_relabd(pchange_focal_par = 40, 
                                        pchange_nonfocal_par = 5,
                                        n_generation = 100,
                                        n_sim = 100)
BH_sens_20_10 <- BH_sensitivity_two_pars_relabd(pchange_focal_par = 20, 
                                         pchange_nonfocal_par = 10,
                                         n_generation = 100,
                                         n_sim = 100)
BH_sens_40_20 <- BH_sensitivity_two_pars_relabd(pchange_focal_par = 40, 
                                         pchange_nonfocal_par = 20,
                                         n_generation = 100,
                                         n_sim = 100)

plot_sim_two_pars_relabd(dat = BH_sens_20_5,
                  pchange_focal_par = 20, 
                  pchange_nonfocal_par = 5)
ggsave("./Outputs/two_pars_sens_plot_20_5_relabd.tiff", width = 10, height = 5, dpi = 600, device = "tiff")
ggsave("./Outputs/two_pars_sens_plot_20_5_relabd.pdf", width = 9, height = 5.5)

plot_sim_two_pars_relabd(dat = BH_sens_40_5,
                  pchange_focal_par = 40, 
                  pchange_nonfocal_par = 5)
ggsave("./Outputs/two_pars_sens_plot_40_5_relabd.tiff", width = 10, height = 5, dpi = 600, device = "tiff")

plot_sim_two_pars_relabd(dat = BH_sens_20_10,
                  pchange_focal_par = 20, 
                  pchange_nonfocal_par = 10)
ggsave("./Outputs/two_pars_sens_plot_20_10_relabd.tiff", width = 10, height = 5, dpi = 600, device = "tiff")

plot_sim_two_pars_relabd(dat = BH_sens_40_20,
                  pchange_focal_par = 40, 
                  pchange_nonfocal_par = 20)
ggsave("./Outputs/two_pars_sens_plot_40_20_relabd.tiff", width = 10, height = 5, dpi = 600, device = "tiff")


### Three parameters
BH_sensitivity_three_pars_relabd <- function(pchange_focal_par = 20, 
                                      pchange_nonfocal_par = 5,
                                      n_generation = 100,
                                      n_sim = 100){
  
  ### Beverton-Holt model 
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
  
  ### Default parameters (coexistence)
  # (1) Number of species
  N <- 2
  
  # (2) Seed germination rate
  g <- c(0.752273, 0.666667)
  
  # (3) Survival rate of ungerminated seeds
  s <- c(0.133750, 0.045000)
  
  # (4) Fecundity
  lambda <- c(2129.949613, 736.667018)
  
  # (5) Interaction matrix
  A <- matrix(c(0.588199, 1.410938,
                0.109412, 0.948399), N, N, byrow = TRUE)
  
  parms <- list(N = N, 
                g = g, 
                s = s, 
                lambda = lambda, 
                A = A)
  
  # Number of generations
  Time <- n_generation
  
  # A matrix to store the outputs
  out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
  colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
  
  # Initial abundance
  out_mat[1, 2:(N+1)] <- rep(10, N)
  
  # Model outputs
  out_mat <- BevertonHolt(out_mat, parms)
  
  # Relative abundance of N2 in the last generation
  N2_equl_default <- out_mat[Time, N+1]/sum(out_mat[Time, N:(N+1)])
  
  
  ################################################################
  ### Sensitivity analysis of three parameters ###################
  ### Fix the parameters of N1 and change the parameters of N2 ###
  
  ### 1. Sensitivity of g & s & lambda
  g_s_lambda_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + pchange_focal_par/100))
    s_temp <- c(s[1], s[2]*(1 + pchange_focal_par/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + pchange_focal_par/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100), 
                       A[2, 2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Relative abundance of N2 in the last generation
    out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
    
  })
  
  ### 2. Sensitivity of g & s & A_21
  g_s_A_21_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + pchange_focal_par/100))
    s_temp <- c(s[1], s[2]*(1 + pchange_focal_par/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + pchange_focal_par/100), 
                       A[2, 2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Relative abundance of N2 in the last generation
    out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
    
  })
  
  ### 3. Sensitivity of g & s & A_22
  g_s_A_22_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + pchange_focal_par/100))
    s_temp <- c(s[1], s[2]*(1 + pchange_focal_par/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100), 
                       A[2, 2]*(1 + pchange_focal_par/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Relative abundance of N2 in the last generation
    out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
    
  })
  
  ### 4. Sensitivity of g & lambda & A_21
  g_lambda_A_21_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + pchange_focal_par/100))
    s_temp <- c(s[1], s[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + pchange_focal_par/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + pchange_focal_par/100), 
                       A[2, 2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Relative abundance of N2 in the last generation
    out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
    
  })
  
  ### 5. Sensitivity of g & lambda & A_22
  g_lambda_A_22_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + pchange_focal_par/100))
    s_temp <- c(s[1], s[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + pchange_focal_par/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100), 
                       A[2, 2]*(1 + pchange_focal_par/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Relative abundance of N2 in the last generation
    out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
    
  })
  
  ### 6. Sensitivity of g & A_21 & A_22
  g_A_21_A_22_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + pchange_focal_par/100))
    s_temp <- c(s[1], s[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + pchange_focal_par/100), 
                       A[2, 2]*(1 + pchange_focal_par/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Relative abundance of N2 in the last generation
    out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
    
  })
  
  ### 7. Sensitivity of s & lambda & A_21
  s_lambda_A_21_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    s_temp <- c(s[1], s[2]*(1 + pchange_focal_par/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + pchange_focal_par/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + pchange_focal_par/100), 
                       A[2, 2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Relative abundance of N2 in the last generation
    out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
    
  })
  
  ### 8. Sensitivity of s & lambda & A_22
  s_lambda_A_22_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    s_temp <- c(s[1], s[2]*(1 + pchange_focal_par/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + pchange_focal_par/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100), 
                       A[2, 2]*(1 + pchange_focal_par/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Relative abundance of N2 in the last generation
    out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
    
  })
  
  ### 9. Sensitivity of s & A_21 & A_22
  s_A_21_A_22_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    s_temp <- c(s[1], s[2]*(1 + pchange_focal_par/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + pchange_focal_par/100), 
                       A[2, 2]*(1 + pchange_focal_par/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Relative abundance of N2 in the last generation
    out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
    
  })
  
  ### 10. Sensitivity of lambda & A_21 & A_22
  lambda_A_21_A_22_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    s_temp <- c(s[1], s[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + pchange_focal_par/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + pchange_focal_par/100), 
                       A[2, 2]*(1 + pchange_focal_par/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Relative abundance of N2 in the last generation
    out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
    
  })  
  
  ### Combine all the simulation results
  pars_sens_df <- data.frame(parameter = rep(c("g_s_lambda", "g_s_A_21", "g_s_A_22", "g_lambda_A_21", "g_lambda_A_22",
                                               "g_A_21_A_22", "s_lambda_A_21", "s_lambda_A_22", "s_A_21_A_22", "lambda_A_21_A_22"), 
                                             each = n_sim), 
                             N2_eql = c(g_s_lambda_sens, g_s_A_21_sens, g_s_A_22_sens, g_lambda_A_21_sens, g_lambda_A_22_sens,
                                        g_A_21_A_22_sens, s_lambda_A_21_sens, s_lambda_A_22_sens, s_A_21_A_22_sens, lambda_A_21_A_22_sens))
  
  return(list(N2_equl_default = N2_equl_default, 
              pars_sens_df = pars_sens_df))
  
}

plot_sim_three_pars_relabd <- function(dat, pchange_focal_par, pchange_nonfocal_par){
  
  # Arrange the parameters by the average N2 abundance of the simulations
  pars_order <- dat$pars_sens_df %>% 
    group_by(parameter) %>% 
    summarise(mean = mean(N2_eql)) %>% 
    mutate(No = 1:10) %>% 
    arrange(desc(mean))
  
  # Reorder the parameter levels for plotting
  dat$pars_sens_df <- dat$pars_sens_df %>% 
    mutate(parameter = factor(parameter, level = pars_order$parameter, ordered = T))
  
  # A vector of x-axis labels
  x_labels <- c(expression(italic(g~+~alpha[21]~+~alpha[22])),
                expression(italic(g~+~lambda~+~alpha[21])),
                expression(italic(g~+~lambda~+~alpha[22])),
                expression(italic(g~+~s~+~alpha[21])),
                expression(italic(g~+~s~+~alpha[22])),
                expression(italic(g~+~s~+~lambda)),
                expression(italic(lambda~+~alpha[21]~+~alpha[22])),
                expression(italic(s~+~alpha[21]~+~alpha[22])),
                expression(italic(s~+~lambda~+~alpha[21])),
                expression(italic(s~+~lambda~+~alpha[22])))
  
  y_line <- ifelse(max(dat$pars_sens_df$N2_eql > 0.99), 1, 0)
  y_lim <- ifelse(max(dat$pars_sens_df$N2_eql > 0.99), 1.1, NA)
  
  # Plot
  ggplot(data = dat$pars_sens_df) + 
    geom_hline(yintercept = dat$N2_equl_default, color = "grey60", linetype = "dashed", size = 1) + 
    geom_hline(yintercept = 0, color = "blue", linetype = "dashed", size = 1) + 
    geom_hline(yintercept = y_line, color = "blue", linetype = "dashed", size = 1) + 
    geom_point(aes(x = parameter, y = N2_eql), position = position_jitter(width = 0.1), 
               color = "grey50", alpha = 0.75) + 
    stat_summary(aes(x = parameter, y = N2_eql), fun.y = mean,
                 fun.ymin = function(x) mean(x) - sd(x), 
                 fun.ymax = function(x) mean(x) + sd(x), 
                 geom = "pointrange", color = "red") + 
    labs(x = "Focal parameter", y = expression(Relative~abundance~of~N[2]), 
         title = glue::glue("Focal parameters: +{pchange_focal_par}% \n Non-focal parameters: ±{pchange_nonfocal_par}%")) +
    scale_x_discrete(labels = x_labels[pars_order$No]) + 
    scale_y_continuous(expand = c(0.25, 0), limits = c(0, y_lim)) + 
    theme_classic() + 
    theme(axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(color = "black", size = 12, angle = 45),
          axis.text.y = element_text(color = "black", size = 14),
          plot.title = element_text(hjust = 0.5, size = 16),
          axis.title.x = element_text(size = 15, margin = margin(t = -5)),
          axis.title.y = element_text(size = 15, margin = margin(r = 8)))
}

### Visualize the simulation results 
BH_sens_20_5 <- BH_sensitivity_three_pars_relabd(pchange_focal_par = 20, 
                                          pchange_nonfocal_par = 5,
                                          n_generation = 100,
                                          n_sim = 100)
BH_sens_40_5 <- BH_sensitivity_three_pars_relabd(pchange_focal_par = 40, 
                                          pchange_nonfocal_par = 5,
                                          n_generation = 100,
                                          n_sim = 100)
BH_sens_20_10 <- BH_sensitivity_three_pars_relabd(pchange_focal_par = 20, 
                                           pchange_nonfocal_par = 10,
                                           n_generation = 100,
                                           n_sim = 100)
BH_sens_40_20 <- BH_sensitivity_three_pars_relabd(pchange_focal_par = 40, 
                                           pchange_nonfocal_par = 20,
                                           n_generation = 100,
                                           n_sim = 100)

plot_sim_three_pars_relabd(dat = BH_sens_20_5,
                    pchange_focal_par = 20, 
                    pchange_nonfocal_par = 5)
ggsave("./Outputs/three_pars_sens_plot_20_5_relabd.tiff", width = 10, height = 5, dpi = 600, device = "tiff")
ggsave("./Outputs/three_pars_sens_plot_20_5_relabd.pdf", width = 9, height = 5.5)

plot_sim_three_pars_relabd(dat = BH_sens_40_5,
                    pchange_focal_par = 40, 
                    pchange_nonfocal_par = 5)
ggsave("./Outputs/three_pars_sens_plot_40_5_relabd.tiff", width = 10, height = 5, dpi = 600, device = "tiff")

plot_sim_three_pars_relabd(dat = BH_sens_20_10,
                    pchange_focal_par = 20, 
                    pchange_nonfocal_par = 10)
ggsave("./Outputs/three_pars_sens_plot_20_10_relabd.tiff", width = 10, height = 5, dpi = 600, device = "tiff")

plot_sim_three_pars_relabd(dat = BH_sens_40_20,
                    pchange_focal_par = 40, 
                    pchange_nonfocal_par = 20)
ggsave("./Outputs/three_pars_sens_plot_40_20_relabd.tiff", width = 10, height = 5, dpi = 600, device = "tiff")


### Four parameters
BH_sensitivity_four_pars_relabd <- function(pchange_focal_par = 20, 
                                     pchange_nonfocal_par = 5,
                                     n_generation = 100,
                                     n_sim = 100){
  
  ### Beverton-Holt model 
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
  
  ### Default parameters (coexistence)
  # (1) Number of species
  N <- 2
  
  # (2) Seed germination rate
  g <- c(0.752273, 0.666667)
  
  # (3) Survival rate of ungerminated seeds
  s <- c(0.133750, 0.045000)
  
  # (4) Fecundity
  lambda <- c(2129.949613, 736.667018)
  
  # (5) Interaction matrix
  A <- matrix(c(0.588199, 1.410938,
                0.109412, 0.948399), N, N, byrow = TRUE)
  
  parms <- list(N = N, 
                g = g, 
                s = s, 
                lambda = lambda, 
                A = A)
  
  # Number of generations
  Time <- n_generation
  
  # A matrix to store the outputs
  out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
  colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
  
  # Initial abundance
  out_mat[1, 2:(N+1)] <- rep(10, N)
  
  # Model outputs
  out_mat <- BevertonHolt(out_mat, parms)
  
  # Relative abundance of N2 in the last generation
  N2_equl_default <- out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
  
  
  ################################################################
  ### Sensitivity analysis of four parameters ####################
  ### Fix the parameters of N1 and change the parameters of N2 ###
  
  ### 1. Sensitivity of g & s & lambda & A_21
  g_s_lambda_A_21_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + pchange_focal_par/100))
    s_temp <- c(s[1], s[2]*(1 + pchange_focal_par/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + pchange_focal_par/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + pchange_focal_par/100), 
                       A[2, 2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Relative abundance of N2 in the last generation
    out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
    
  })
  
  ### 2. Sensitivity of g & s & A_21 & A_22
  g_s_A_21_A_22_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + pchange_focal_par/100))
    s_temp <- c(s[1], s[2]*(1 + pchange_focal_par/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + pchange_focal_par/100), 
                       A[2, 2]*(1 + pchange_focal_par/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Relative abundance of N2 in the last generation
    out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
    
  })
  
  ### 3. Sensitivity of g & s & lambda & A_22
  g_s_lambda_A_22_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + pchange_focal_par/100))
    s_temp <- c(s[1], s[2]*(1 + pchange_focal_par/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + pchange_focal_par/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100), 
                       A[2, 2]*(1 + pchange_focal_par/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Relative abundance of N2 in the last generation
    out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
    
  })
  
  ### 4. Sensitivity of g & lambda & A_21 & A_22
  g_lambda_A_21_A_22_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + pchange_focal_par/100))
    s_temp <- c(s[1], s[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + pchange_focal_par/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + pchange_focal_par/100), 
                       A[2, 2]*(1 + pchange_focal_par/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Relative abundance of N2 in the last generation
    out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
    
  })
  
  ### 5. Sensitivity of s & lambda & A_21 & A_22
  s_lambda_A_21_A_22_sens <- sapply(1:n_sim, function(x){
    
    # New parameters
    g_temp <- c(g[1], g[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    s_temp <- c(s[1], s[2]*(1 + pchange_focal_par/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + pchange_focal_par/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + pchange_focal_par/100), 
                       A[2, 2]*(1 + pchange_focal_par/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    # Number of generations
    Time <- n_generation
    
    # A matrix to store the outputs
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    # Initial abundance
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    # Model outputs
    out_mat <- BevertonHolt(out_mat, parms_temp)
    
    # Relative abundance of N2 in the last generation
    out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
    
  })
  
  ### Combine all the simulation results
  pars_sens_df <- data.frame(parameter = rep(c("g_s_lambda_A_21", "g_s_A_21_A_22", "g_s_lambda_A_22", "g_lambda_A_21_A_22", "s_lambda_A_21_A_22"), each = n_sim), 
                             N2_eql = c(g_s_lambda_A_21_sens, g_s_A_21_A_22_sens, g_s_lambda_A_22_sens, g_lambda_A_21_A_22_sens, s_lambda_A_21_A_22_sens))
  
  return(list(N2_equl_default = N2_equl_default, 
              pars_sens_df = pars_sens_df))
  
}

plot_sim_four_pars_relabd <- function(dat, pchange_focal_par, pchange_nonfocal_par){
  
  # Arrange the parameters by the average N2 abundance of the simulations
  pars_order <- dat$pars_sens_df %>% 
    group_by(parameter) %>% 
    summarise(mean = mean(N2_eql)) %>% 
    mutate(No = 1:5) %>% 
    arrange(desc(mean))
  
  # Reorder the parameter levels for plotting
  dat$pars_sens_df <- dat$pars_sens_df %>% 
    mutate(parameter = factor(parameter, level = pars_order$parameter, ordered = T))
  
  # A vector of x-axis labels
  x_labels <- c(expression(italic(g~+~lambda~+~alpha[21]~+~alpha[22])),
                expression(italic(g~+~s~+~alpha[21]~+~alpha[22])),
                expression(italic(g~+~s~+~lambda~+~alpha[21])),
                expression(italic(g~+~s~+~lambda~+~alpha[22])),
                expression(italic(s~+~lambda~+~alpha[21]~+~alpha[22])))
  
  y_line <- ifelse(max(dat$pars_sens_df$N2_eql > 0.99), 1, 0)
  y_lim <- ifelse(max(dat$pars_sens_df$N2_eql > 0.99), 1.1, NA)
  
  # Plot
  ggplot(data = dat$pars_sens_df) + 
    geom_hline(yintercept = dat$N2_equl_default, color = "grey60", linetype = "dashed", size = 1) + 
    geom_hline(yintercept = 0, color = "blue", linetype = "dashed", size = 1) + 
    geom_hline(yintercept = y_line, color = "blue", linetype = "dashed", size = 1) + 
    geom_point(aes(x = parameter, y = N2_eql), position = position_jitter(width = 0.1), 
               color = "grey50", alpha = 0.75) + 
    stat_summary(aes(x = parameter, y = N2_eql), fun.y = mean,
                 fun.ymin = function(x) mean(x) - sd(x), 
                 fun.ymax = function(x) mean(x) + sd(x), 
                 geom = "pointrange", color = "red") + 
    labs(x = "Focal parameter", y = expression(Relative~abundance~of~N[2]), 
         title = glue::glue("Focal parameters: +{pchange_focal_par}% \n Non-focal parameter: ±{pchange_nonfocal_par}%")) +
    scale_x_discrete(labels = x_labels[pars_order$No]) + 
    scale_y_continuous(expand = c(0.25, 0), limits = c(0, y_lim)) + 
    theme_classic() + 
    theme(axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(color = "black", size = 12, angle = 45),
          axis.text.y = element_text(color = "black", size = 14),
          plot.title = element_text(hjust = 0.5, size = 16),
          axis.title.x = element_text(size = 15, margin = margin(t = -5)),
          axis.title.y = element_text(size = 15, margin = margin(r = 8)))
}

### Visualize the simulation results 
BH_sens_20_5 <- BH_sensitivity_four_pars_relabd(pchange_focal_par = 20, 
                                         pchange_nonfocal_par = 5,
                                         n_generation = 100,
                                         n_sim = 100)
BH_sens_40_5 <- BH_sensitivity_four_pars_relabd(pchange_focal_par = 40, 
                                         pchange_nonfocal_par = 5,
                                         n_generation = 100,
                                         n_sim = 100)
BH_sens_20_10 <- BH_sensitivity_four_pars_relabd(pchange_focal_par = 20, 
                                          pchange_nonfocal_par = 10,
                                          n_generation = 100,
                                          n_sim = 100)
BH_sens_40_20 <- BH_sensitivity_four_pars_relabd(pchange_focal_par = 40, 
                                          pchange_nonfocal_par = 20,
                                          n_generation = 100,
                                          n_sim = 100)

plot_sim_four_pars_relabd(dat = BH_sens_20_5,
                   pchange_focal_par = 20, 
                   pchange_nonfocal_par = 5)
ggsave("./Outputs/four_pars_sens_plot_20_5_relabd.tiff", width = 7, height = 5, dpi = 600, device = "tiff")
ggsave("./Outputs/four_pars_sens_plot_20_5_relabd.pdf", width = 7, height = 5)

plot_sim_four_pars_relabd(dat = BH_sens_40_5,
                   pchange_focal_par = 40, 
                   pchange_nonfocal_par = 5)
ggsave("./Outputs/four_pars_sens_plot_40_5_relabd.tiff", width = 7, height = 5, dpi = 600, device = "tiff")

plot_sim_four_pars_relabd(dat = BH_sens_20_10,
                   pchange_focal_par = 20, 
                   pchange_nonfocal_par = 10)
ggsave("./Outputs/four_pars_sens_plot_20_10_relabd.tiff", width = 7, height = 5, dpi = 600, device = "tiff")

plot_sim_four_pars_relabd(dat = BH_sens_40_20,
                   pchange_focal_par = 40, 
                   pchange_nonfocal_par = 20)
ggsave("./Outputs/four_pars_sens_plot_40_20_relabd.tiff", width = 7, height = 5, dpi = 600, device = "tiff")



################################################################################
################################################################################
################################################################################
### One-parameter sensitivity analysis with different parameter settings #######

### Scenario 1. N1 = FE & N2 = HO
BH_sensitivity_FE_HO <- function(pchange_focal_par = 40, 
                                 pchange_nonfocal_par = 5,
                                 n_generation = 100,
                                 n_sim = 100){
  
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
  
  ### Default parameters
  N <- 2
  g <- c(0.752273, 0.666667)
  s <- c(0.133750, 0.045000)
  lambda <- c(2129.949613, 736.667018)
  A <- matrix(c(0.588199, 1.410938,
                0.109412, 0.948399), N, N, byrow = TRUE)
  
  parms <- list(N = N, 
                g = g, 
                s = s, 
                lambda = lambda, 
                A = A)
  
  Time <- n_generation
  
  out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
  colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
  
  out_mat[1, 2:(N+1)] <- rep(10, N)
  out_mat <- BevertonHolt(out_mat, parms)

  N2_equl_default <- out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
  
  
  ### All parameters increased
  g_temp <- c(g[1], g[2]*(1 + pchange_focal_par/100))
  s_temp <- c(s[1], s[2]*(1 + pchange_focal_par/100))
  lambda_temp <- c(lambda[1], lambda[2]*(1 + pchange_focal_par/100))
  A_temp <- matrix(c(A[1, 1], A[1, 2], 
                     A[2, 1]*(1 + pchange_focal_par/100), 
                     A[2, 2]*(1 + pchange_focal_par/100)), 
                   N, N, byrow = TRUE)
  
  parms_temp <- list(N = N, 
                     g = g_temp, 
                     s = s_temp, 
                     lambda = lambda_temp, 
                     A = A_temp)
  
  Time <- n_generation
  
  out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
  colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
  
  out_mat[1, 2:(N+1)] <- rep(10, N)
  out_mat <- BevertonHolt(out_mat, parms_temp)

  N2_equl_all_increase <- out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
  
  
  ### One parameter changed at once
  # (1) Sensitivity of g
  g_sens <- sapply(1:n_sim, function(x){
    
    g_temp <- c(g[1], g[2]*(1 + pchange_focal_par/100))
    s_temp <- c(s[1], s[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100), 
                       A[2, 2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    Time <- n_generation
    
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    out_mat <- BevertonHolt(out_mat, parms_temp)
    out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
    
  })
  
  # (2) Sensitivity of s
  s_sens <- sapply(1:n_sim, function(x){
    
    g_temp <- c(g[1], g[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    s_temp <- c(s[1], s[2]*(1 + pchange_focal_par/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100), 
                       A[2, 2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    Time <- n_generation
    
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    out_mat <- BevertonHolt(out_mat, parms_temp)
    out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
    
  })
  
  # (3) Sensitivity of lambda
  lambda_sens <- sapply(1:n_sim, function(x){
    
    g_temp <- c(g[1], g[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    s_temp <- c(s[1], s[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + pchange_focal_par/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100), 
                       A[2, 2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    Time <- n_generation

    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    out_mat <- BevertonHolt(out_mat, parms_temp)
    out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
    
  })
  
  # (4) Sensitivity of A_21
  A_21_sens <- sapply(1:n_sim, function(x){
    
    g_temp <- c(g[1], g[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    s_temp <- c(s[1], s[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + pchange_focal_par/100), 
                       A[2, 2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    Time <- n_generation
    
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    out_mat <- BevertonHolt(out_mat, parms_temp)
    out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
    
  })
  
  # (5) Sensitivity of A_22
  A_22_sens <- sapply(1:n_sim, function(x){
    
    g_temp <- c(g[1], g[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    s_temp <- c(s[1], s[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100), 
                       A[2, 2]*(1 + pchange_focal_par/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    Time <- n_generation
    
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    out_mat <- BevertonHolt(out_mat, parms_temp)
    out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
    
  })
  
  # Combine all the simulation results
  pars_sens_df <- data.frame(parameter = rep(c("g", "s", "lambda", "A_21", "A_22"), each = n_sim), 
                             N2_eql = c(g_sens, s_sens, lambda_sens, A_21_sens, A_22_sens))
  
  return(list(N2_equl_default = N2_equl_default,
              N2_equl_all_increase = N2_equl_all_increase,
              pars_sens_df = pars_sens_df))
}

### Plot function
plot_sim_one_par_relabd <- function(dat, pchange_focal_par, pchange_nonfocal_par){
  
  pars_order <- dat$pars_sens_df %>% 
    group_by(parameter) %>% 
    summarise(mean = mean(N2_eql)) %>% 
    mutate(No = 1:5) %>% 
    arrange(desc(mean))

  dat$pars_sens_df <- dat$pars_sens_df %>% 
    mutate(parameter = factor(parameter, level = pars_order$parameter, ordered = T))
  
  x_labels <- c(expression(italic(alpha[21])),
                expression(italic(alpha[22])),
                expression(italic(g)),
                expression(italic(lambda)),
                expression(italic(s)))
  
  ggplot(data = dat$pars_sens_df) + 
    geom_hline(yintercept = dat$N2_equl_default, color = "grey60", linetype = "dashed", size = 1) + 
    geom_hline(yintercept = dat$N2_equl_all_increase, color = "blue", linetype = "dashed", size = 1) + 
    geom_point(aes(x = parameter, y = N2_eql), position = position_jitter(width = 0.1), 
               color = "grey50", alpha = 0.75) + 
    stat_summary(aes(x = parameter, y = N2_eql), fun.y = mean,
                 fun.ymin = function(x) mean(x) - sd(x), 
                 fun.ymax = function(x) mean(x) + sd(x), 
                 geom = "pointrange", color = "red") + 
    labs(x = "Focal parameter", y = expression(Relative~abundance~of~N[2]), 
         title = glue::glue("Focal parameter: +{pchange_focal_par}% \n Non-focal parameters: ±{pchange_nonfocal_par}%")) +
    scale_x_discrete(labels = x_labels[pars_order$No]) + 
    scale_y_continuous(expand = c(0.25, 0)) + 
    theme_classic() + 
    theme(axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(color = "black", size = 12),
          axis.text.y = element_text(color = "black", size = 14),
          plot.title = element_text(hjust = 0.5, size = 16),
          axis.title.x = element_text(size = 15, margin = margin(t = 10)),
          axis.title.y = element_text(size = 15, margin = margin(r = 8)))
}

### Visualize the simulation results 
BH_sens_FE_HO <- BH_sensitivity_FE_HO()

plot_sim_one_par_relabd(dat = BH_sens_FE_HO,
                        pchange_focal_par = 40, 
                        pchange_nonfocal_par = 5)

ggsave("./Outputs/BH_sens_FE_HO.pdf", width = 5, height = 4)
ggsave("./Outputs/BH_sens_FE_HO.tiff", width = 5, height = 4)

### Scenario 2. N1 = FE & N2 = AC
BH_sensitivity_FE_AC <- function(pchange_focal_par = 40, 
                                 pchange_nonfocal_par = 5,
                                 n_generation = 100,
                                 n_sim = 100){
  
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
  
  ### Default parameters
  N <- 2
  g <- c(0.752273, 0.224031)
  s <- c(0.113750, 0.262546)
  lambda <- c(2129.949613, 2199.20)
  A <- matrix(c(0.588199, 0.252077,
                0.211761, 0.272905), N, N, byrow = TRUE)
  
  parms <- list(N = N, 
                g = g, 
                s = s, 
                lambda = lambda, 
                A = A)
  
  Time <- n_generation
  
  out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
  colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
  
  out_mat[1, 2:(N+1)] <- rep(10, N)
  out_mat <- BevertonHolt(out_mat, parms)
  
  N2_equl_default <- out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
  
  
  ### All parameters increased
  g_temp <- c(g[1], g[2]*(1 + pchange_focal_par/100))
  s_temp <- c(s[1], s[2]*(1 + pchange_focal_par/100))
  lambda_temp <- c(lambda[1], lambda[2]*(1 + pchange_focal_par/100))
  A_temp <- matrix(c(A[1, 1], A[1, 2], 
                     A[2, 1]*(1 + pchange_focal_par/100), 
                     A[2, 2]*(1 + pchange_focal_par/100)), 
                   N, N, byrow = TRUE)
  
  parms_temp <- list(N = N, 
                     g = g_temp, 
                     s = s_temp, 
                     lambda = lambda_temp, 
                     A = A_temp)
  
  Time <- n_generation
  
  out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
  colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
  
  out_mat[1, 2:(N+1)] <- rep(10, N)
  out_mat <- BevertonHolt(out_mat, parms_temp)
  
  N2_equl_all_increase <- out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
  
  
  ### One parameter changed at once
  # (1) Sensitivity of g
  g_sens <- sapply(1:n_sim, function(x){
    
    g_temp <- c(g[1], g[2]*(1 + pchange_focal_par/100))
    s_temp <- c(s[1], s[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100), 
                       A[2, 2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    Time <- n_generation
    
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    out_mat <- BevertonHolt(out_mat, parms_temp)
    out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
    
  })
  
  # (2) Sensitivity of s
  s_sens <- sapply(1:n_sim, function(x){
    
    g_temp <- c(g[1], g[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    s_temp <- c(s[1], s[2]*(1 + pchange_focal_par/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100), 
                       A[2, 2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    Time <- n_generation
    
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    out_mat <- BevertonHolt(out_mat, parms_temp)
    out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
    
  })
  
  # (3) Sensitivity of lambda
  lambda_sens <- sapply(1:n_sim, function(x){
    
    g_temp <- c(g[1], g[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    s_temp <- c(s[1], s[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + pchange_focal_par/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100), 
                       A[2, 2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    Time <- n_generation
    
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    out_mat <- BevertonHolt(out_mat, parms_temp)
    out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
    
  })
  
  # (4) Sensitivity of A_21
  A_21_sens <- sapply(1:n_sim, function(x){
    
    g_temp <- c(g[1], g[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    s_temp <- c(s[1], s[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + pchange_focal_par/100), 
                       A[2, 2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    Time <- n_generation
    
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    out_mat <- BevertonHolt(out_mat, parms_temp)
    out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
    
  })
  
  # (5) Sensitivity of A_22
  A_22_sens <- sapply(1:n_sim, function(x){
    
    g_temp <- c(g[1], g[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    s_temp <- c(s[1], s[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100), 
                       A[2, 2]*(1 + pchange_focal_par/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    Time <- n_generation
    
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    out_mat <- BevertonHolt(out_mat, parms_temp)
    out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
    
  })
  
  # Combine all the simulation results
  pars_sens_df <- data.frame(parameter = rep(c("g", "s", "lambda", "A_21", "A_22"), each = n_sim), 
                             N2_eql = c(g_sens, s_sens, lambda_sens, A_21_sens, A_22_sens))
  
  return(list(N2_equl_default = N2_equl_default,
              N2_equl_all_increase = N2_equl_all_increase,
              pars_sens_df = pars_sens_df))
}

### Plot function
plot_sim_one_par_relabd <- function(dat, pchange_focal_par, pchange_nonfocal_par){
  
  pars_order <- dat$pars_sens_df %>% 
    group_by(parameter) %>% 
    summarise(mean = mean(N2_eql)) %>% 
    mutate(No = 1:5) %>% 
    arrange(desc(mean))
  
  dat$pars_sens_df <- dat$pars_sens_df %>% 
    mutate(parameter = factor(parameter, level = pars_order$parameter, ordered = T))
  
  x_labels <- c(expression(italic(alpha[21])),
                expression(italic(alpha[22])),
                expression(italic(g)),
                expression(italic(lambda)),
                expression(italic(s)))
  
  ggplot(data = dat$pars_sens_df) + 
    geom_hline(yintercept = dat$N2_equl_default, color = "grey60", linetype = "dashed", size = 1) + 
    geom_hline(yintercept = dat$N2_equl_all_increase, color = "blue", linetype = "dashed", size = 1) + 
    geom_point(aes(x = parameter, y = N2_eql), position = position_jitter(width = 0.1), 
               color = "grey50", alpha = 0.75) + 
    stat_summary(aes(x = parameter, y = N2_eql), fun.y = mean,
                 fun.ymin = function(x) mean(x) - sd(x), 
                 fun.ymax = function(x) mean(x) + sd(x), 
                 geom = "pointrange", color = "red") + 
    labs(x = "Focal parameter", y = expression(Relative~abundance~of~N[2]), 
         title = glue::glue("Focal parameter: +{pchange_focal_par}% \n Non-focal parameters: ±{pchange_nonfocal_par}%")) +
    scale_x_discrete(labels = x_labels[pars_order$No]) + 
    scale_y_continuous(expand = c(0.25, 0)) + 
    theme_classic() + 
    theme(axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(color = "black", size = 12),
          axis.text.y = element_text(color = "black", size = 14),
          plot.title = element_text(hjust = 0.5, size = 16),
          axis.title.x = element_text(size = 15, margin = margin(t = 10)),
          axis.title.y = element_text(size = 15, margin = margin(r = 8)))
}

### Visualize the simulation results 
BH_sens_FE_AC <- BH_sensitivity_FE_AC()

plot_sim_one_par_relabd(dat = BH_sens_FE_AC,
                        pchange_focal_par = 40, 
                        pchange_nonfocal_par = 5)

ggsave("./Outputs/BH_sens_FE_AC.pdf", width = 5, height = 4)


### Scenario 3. N1 = AC & N2 = HO
BH_sensitivity_AC_HO <- function(pchange_focal_par = 40, 
                                 pchange_nonfocal_par = 5,
                                 n_generation = 100,
                                 n_sim = 100){
  
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
  
  ### Default parameters
  N <- 2
  g <- c(0.224031, 0.666667)
  s <- c(0.262546, 0.045000)
  lambda <- c(2199.20, 736.667018)
  A <- matrix(c(0.272905, 0.941539,
                0.137555, 0.948399), N, N, byrow = TRUE)
  
  parms <- list(N = N, 
                g = g, 
                s = s, 
                lambda = lambda, 
                A = A)
  
  Time <- n_generation
  
  out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
  colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
  
  out_mat[1, 2:(N+1)] <- rep(10, N)
  out_mat <- BevertonHolt(out_mat, parms)
  
  N2_equl_default <- out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
  
  
  ### All parameters increased
  g_temp <- c(g[1], g[2]*(1 + pchange_focal_par/100))
  s_temp <- c(s[1], s[2]*(1 + pchange_focal_par/100))
  lambda_temp <- c(lambda[1], lambda[2]*(1 + pchange_focal_par/100))
  A_temp <- matrix(c(A[1, 1], A[1, 2], 
                     A[2, 1]*(1 + pchange_focal_par/100), 
                     A[2, 2]*(1 + pchange_focal_par/100)), 
                   N, N, byrow = TRUE)
  
  parms_temp <- list(N = N, 
                     g = g_temp, 
                     s = s_temp, 
                     lambda = lambda_temp, 
                     A = A_temp)
  
  Time <- n_generation
  
  out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
  colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
  
  out_mat[1, 2:(N+1)] <- rep(10, N)
  out_mat <- BevertonHolt(out_mat, parms_temp)
  
  N2_equl_all_increase <- out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
  
  
  ### One parameter changed at once
  # (1) Sensitivity of g
  g_sens <- sapply(1:n_sim, function(x){
    
    g_temp <- c(g[1], g[2]*(1 + pchange_focal_par/100))
    s_temp <- c(s[1], s[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100), 
                       A[2, 2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    Time <- n_generation
    
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    out_mat <- BevertonHolt(out_mat, parms_temp)
    out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
    
  })
  
  # (2) Sensitivity of s
  s_sens <- sapply(1:n_sim, function(x){
    
    g_temp <- c(g[1], g[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    s_temp <- c(s[1], s[2]*(1 + pchange_focal_par/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100), 
                       A[2, 2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    Time <- n_generation
    
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    out_mat <- BevertonHolt(out_mat, parms_temp)
    out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
    
  })
  
  # (3) Sensitivity of lambda
  lambda_sens <- sapply(1:n_sim, function(x){
    
    g_temp <- c(g[1], g[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    s_temp <- c(s[1], s[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + pchange_focal_par/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100), 
                       A[2, 2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    Time <- n_generation
    
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    out_mat <- BevertonHolt(out_mat, parms_temp)
    out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
    
  })
  
  # (4) Sensitivity of A_21
  A_21_sens <- sapply(1:n_sim, function(x){
    
    g_temp <- c(g[1], g[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    s_temp <- c(s[1], s[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + pchange_focal_par/100), 
                       A[2, 2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    Time <- n_generation
    
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    out_mat <- BevertonHolt(out_mat, parms_temp)
    out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
    
  })
  
  # (5) Sensitivity of A_22
  A_22_sens <- sapply(1:n_sim, function(x){
    
    g_temp <- c(g[1], g[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    s_temp <- c(s[1], s[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100), 
                       A[2, 2]*(1 + pchange_focal_par/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    Time <- n_generation
    
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    out_mat <- BevertonHolt(out_mat, parms_temp)
    out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
    
  })
  
  # Combine all the simulation results
  pars_sens_df <- data.frame(parameter = rep(c("g", "s", "lambda", "A_21", "A_22"), each = n_sim), 
                             N2_eql = c(g_sens, s_sens, lambda_sens, A_21_sens, A_22_sens))
  
  return(list(N2_equl_default = N2_equl_default,
              N2_equl_all_increase = N2_equl_all_increase,
              pars_sens_df = pars_sens_df))
}

### Plot function
plot_sim_one_par_relabd <- function(dat, pchange_focal_par, pchange_nonfocal_par){
  
  pars_order <- dat$pars_sens_df %>% 
    group_by(parameter) %>% 
    summarise(mean = mean(N2_eql)) %>% 
    mutate(No = 1:5) %>% 
    arrange(desc(mean))
  
  dat$pars_sens_df <- dat$pars_sens_df %>% 
    mutate(parameter = factor(parameter, level = pars_order$parameter, ordered = T))
  
  x_labels <- c(expression(italic(alpha[21])),
                expression(italic(alpha[22])),
                expression(italic(g)),
                expression(italic(lambda)),
                expression(italic(s)))
  
  ggplot(data = dat$pars_sens_df) + 
    geom_hline(yintercept = dat$N2_equl_default, color = "grey60", linetype = "dashed", size = 1) + 
    geom_hline(yintercept = dat$N2_equl_all_increase, color = "blue", linetype = "dashed", size = 1) + 
    geom_point(aes(x = parameter, y = N2_eql), position = position_jitter(width = 0.1), 
               color = "grey50", alpha = 0.75) + 
    stat_summary(aes(x = parameter, y = N2_eql), fun.y = mean,
                 fun.ymin = function(x) mean(x) - sd(x), 
                 fun.ymax = function(x) mean(x) + sd(x), 
                 geom = "pointrange", color = "red") + 
    labs(x = "Focal parameter", y = expression(Relative~abundance~of~N[2]), 
         title = glue::glue("Focal parameter: +{pchange_focal_par}% \n Non-focal parameters: ±{pchange_nonfocal_par}%")) +
    scale_x_discrete(labels = x_labels[pars_order$No]) + 
    scale_y_continuous(expand = c(0.25, 0)) + 
    theme_classic() + 
    theme(axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(color = "black", size = 12),
          axis.text.y = element_text(color = "black", size = 14),
          plot.title = element_text(hjust = 0.5, size = 16),
          axis.title.x = element_text(size = 15, margin = margin(t = 10)),
          axis.title.y = element_text(size = 15, margin = margin(r = 8)))
}

### Visualize the simulation results 
BH_sens_AC_HO <- BH_sensitivity_AC_HO()

plot_sim_one_par_relabd(dat = BH_sens_AC_HO,
                        pchange_focal_par = 40, 
                        pchange_nonfocal_par = 5)

ggsave("./Outputs/BH_sens_AC_HO.pdf", width = 5, height = 4)


### Scenario 4. N1 = AC & N2 = SA
BH_sensitivity_AC_SA <- function(pchange_focal_par = 40, 
                                 pchange_nonfocal_par = 5,
                                 n_generation = 100,
                                 n_sim = 100){
  
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
  
  ### Default parameters
  N <- 2
  g <- c(0.224031, 0.349206)
  s <- c(0.262546, 0.338957)
  lambda <- c(2199.20, 1023.584)
  A <- matrix(c(0.272905, 0.017925,
                0.149194, 0.152118), N, N, byrow = TRUE)
  
  parms <- list(N = N, 
                g = g, 
                s = s, 
                lambda = lambda, 
                A = A)
  
  Time <- n_generation
  
  out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
  colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
  
  out_mat[1, 2:(N+1)] <- rep(10, N)
  out_mat <- BevertonHolt(out_mat, parms)
  
  N2_equl_default <- out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
  
  
  ### All parameters increased
  g_temp <- c(g[1], g[2]*(1 + pchange_focal_par/100))
  s_temp <- c(s[1], s[2]*(1 + pchange_focal_par/100))
  lambda_temp <- c(lambda[1], lambda[2]*(1 + pchange_focal_par/100))
  A_temp <- matrix(c(A[1, 1], A[1, 2], 
                     A[2, 1]*(1 + pchange_focal_par/100), 
                     A[2, 2]*(1 + pchange_focal_par/100)), 
                   N, N, byrow = TRUE)
  
  parms_temp <- list(N = N, 
                     g = g_temp, 
                     s = s_temp, 
                     lambda = lambda_temp, 
                     A = A_temp)
  
  Time <- n_generation
  
  out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
  colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
  
  out_mat[1, 2:(N+1)] <- rep(10, N)
  out_mat <- BevertonHolt(out_mat, parms_temp)
  
  N2_equl_all_increase <- out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
  
  
  ### One parameter changed at once
  # (1) Sensitivity of g
  g_sens <- sapply(1:n_sim, function(x){
    
    g_temp <- c(g[1], g[2]*(1 + pchange_focal_par/100))
    s_temp <- c(s[1], s[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100), 
                       A[2, 2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    Time <- n_generation
    
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    out_mat <- BevertonHolt(out_mat, parms_temp)
    out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
    
  })
  
  # (2) Sensitivity of s
  s_sens <- sapply(1:n_sim, function(x){
    
    g_temp <- c(g[1], g[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    s_temp <- c(s[1], s[2]*(1 + pchange_focal_par/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100), 
                       A[2, 2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    Time <- n_generation
    
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    out_mat <- BevertonHolt(out_mat, parms_temp)
    out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
    
  })
  
  # (3) Sensitivity of lambda
  lambda_sens <- sapply(1:n_sim, function(x){
    
    g_temp <- c(g[1], g[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    s_temp <- c(s[1], s[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + pchange_focal_par/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100), 
                       A[2, 2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    Time <- n_generation
    
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    out_mat <- BevertonHolt(out_mat, parms_temp)
    out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
    
  })
  
  # (4) Sensitivity of A_21
  A_21_sens <- sapply(1:n_sim, function(x){
    
    g_temp <- c(g[1], g[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    s_temp <- c(s[1], s[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + pchange_focal_par/100), 
                       A[2, 2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    Time <- n_generation
    
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    out_mat <- BevertonHolt(out_mat, parms_temp)
    out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
    
  })
  
  # (5) Sensitivity of A_22
  A_22_sens <- sapply(1:n_sim, function(x){
    
    g_temp <- c(g[1], g[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    s_temp <- c(s[1], s[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    lambda_temp <- c(lambda[1], lambda[2]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100))
    A_temp <- matrix(c(A[1, 1], A[1, 2], 
                       A[2, 1]*(1 + runif(n = 1, min = -pchange_nonfocal_par, max = pchange_nonfocal_par)/100), 
                       A[2, 2]*(1 + pchange_focal_par/100)), 
                     N, N, byrow = TRUE)
    
    parms_temp <- list(N = N, 
                       g = g_temp, 
                       s = s_temp, 
                       lambda = lambda_temp, 
                       A = A_temp)
    
    Time <- n_generation
    
    out_mat <- matrix(c(c(0:Time), rep(0, N*(Time+1))), (Time+1), N+1, byrow = F)
    colnames(out_mat) <- c("Time", paste("N", c(1:N), sep = ""))
    
    out_mat[1, 2:(N+1)] <- rep(10, N)
    
    out_mat <- BevertonHolt(out_mat, parms_temp)
    out_mat[Time, N+1]/sum(out_mat[Time, 2:(N+1)])
    
  })
  
  # Combine all the simulation results
  pars_sens_df <- data.frame(parameter = rep(c("g", "s", "lambda", "A_21", "A_22"), each = n_sim), 
                             N2_eql = c(g_sens, s_sens, lambda_sens, A_21_sens, A_22_sens))
  
  return(list(N2_equl_default = N2_equl_default,
              N2_equl_all_increase = N2_equl_all_increase,
              pars_sens_df = pars_sens_df))
}

### Plot function
plot_sim_one_par_relabd <- function(dat, pchange_focal_par, pchange_nonfocal_par){
  
  pars_order <- dat$pars_sens_df %>% 
    group_by(parameter) %>% 
    summarise(mean = mean(N2_eql)) %>% 
    mutate(No = 1:5) %>% 
    arrange(desc(mean))
  
  dat$pars_sens_df <- dat$pars_sens_df %>% 
    mutate(parameter = factor(parameter, level = pars_order$parameter, ordered = T))
  
  x_labels <- c(expression(italic(alpha[21])),
                expression(italic(alpha[22])),
                expression(italic(g)),
                expression(italic(lambda)),
                expression(italic(s)))
  
  ggplot(data = dat$pars_sens_df) + 
    geom_hline(yintercept = dat$N2_equl_default, color = "grey60", linetype = "dashed", size = 1) + 
    geom_hline(yintercept = dat$N2_equl_all_increase, color = "blue", linetype = "dashed", size = 1) + 
    geom_point(aes(x = parameter, y = N2_eql), position = position_jitter(width = 0.1), 
               color = "grey50", alpha = 0.75) + 
    stat_summary(aes(x = parameter, y = N2_eql), fun.y = mean,
                 fun.ymin = function(x) mean(x) - sd(x), 
                 fun.ymax = function(x) mean(x) + sd(x), 
                 geom = "pointrange", color = "red") + 
    labs(x = "Focal parameter", y = expression(Relative~abundance~of~N[2]), 
         title = glue::glue("Focal parameter: +{pchange_focal_par}% \n Non-focal parameters: ±{pchange_nonfocal_par}%")) +
    scale_x_discrete(labels = x_labels[pars_order$No]) + 
    scale_y_continuous(expand = c(0.25, 0)) + 
    theme_classic() + 
    theme(axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(color = "black", size = 12),
          axis.text.y = element_text(color = "black", size = 14),
          plot.title = element_text(hjust = 0.5, size = 16),
          axis.title.x = element_text(size = 15, margin = margin(t = 10)),
          axis.title.y = element_text(size = 15, margin = margin(r = 8)))
}

### Visualize the simulation results 
BH_sens_AC_SA <- BH_sensitivity_AC_SA()

plot_sim_one_par_relabd(dat = BH_sens_AC_SA,
                        pchange_focal_par = 40, 
                        pchange_nonfocal_par = 5)

ggsave("./Outputs/BH_sens_AC_SA.pdf", width = 5, height = 4)





