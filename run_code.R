## ---------------------------
##
## Script name: run_code.R
##
## Purpose of script:
##
## Author: Steve Soutar
##
## Date Created: 2024-02-08
##
## Email: Steven.Soutar@arcturis.data.com
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

library(rjags)
library(ggplot2)
library(ggpubr)
library(survival)

###########################
###########################
# Warm-up - no binary
# confounder.
###########################
###########################

N_ctrl <- 500 # Number in control arm.
N_trt <- 500 # Number in treatment arm.
N <- N_ctrl + N_trt
p_z <- 0.5 # Probability of being treated.

# Simulate treatment indicator:
z <- sample(c(0, 1), N, replace = TRUE, prob = c(1 - p_z, p_z))
lambda_base <- 0.1 # Baseline hazard.
beta_z <- log(0.5) # log(HR)

lambda <- lambda_base*exp(beta_z*z)
time <- rexp(N, rate = lambda)

# Apply administrative censoring:
study_length <- 30
delta <- ifelse(time < study_length, 1, 0)
time[delta == 0] <- study_length

# Time horizon for rMST:
tau <- 6

data_list_no_confounding <- list(N = N, z = z, time = time, 
                                 delta = delta, mu_log_lambda = mu_log_lambda, 
                                 tau_log_lambda = tau_log_lambda, 
                                 mu_beta_z = mu_beta_z, 
                                 tau_beta_z = tau_beta_z,
                                 tau = tau)

# Compile model:
model_no_confounding <- jags.model(file = textConnection(model_string_no_confounding), 
                                   data = data_list_no_confounding, n.chains = 1)

param_no_confounding <- c("lambda_base", "HR") # Parameters of interest.

N_iter <- 50000 # Total number of iterations.
burn_in <- 20000 # Burn-in period.

# Extract posterior sample:
posterior_sample_no_confounding <- coda.samples(model_no_confounding, 
                                                variable.names = param_no_confounding,
                                                n.iter = N_iter)[[1]]

posterior_sample_HR <- posterior_sample_no_confounding[-c(burn_in:N_iter), 1]
posterior_sample_lambda_base <- posterior_sample_no_confounding[-c(burn_in:N_iter), 2]

#posterior_sample_rMST_ctrl <- posterior_sample_no_confounding[-c(burn_in:N_iter), 3]
#posterior_sample_rMST_trt <- posterior_sample_no_confounding[-c(burn_in:N_iter), 4]

#####################
#####################
# Examine posterior
# sample:
#####################
#####################

hist(posterior_sample_HR)
hist(posterior_sample_lambda_base)

hist(posterior_sample_rMST_ctrl)
hist(posterior_sample_rMST_trt)

hist(posterior_sample_rMST_trt - posterior_sample_rMST_ctrl)

#####################################################################

####################
####################
# Binary confounder
####################
####################

N_ctrl <- 200
N_trt <- 200
N <- N_ctrl + N_trt
p_u <- 0.5 # = Pr(U = 1).

# Simulate confounder:
u <- sample(c(0, 1), N, prob = c(1 - p_u, p_u), replace = TRUE)

# Simulate treatment indicator:
p_z <- 1/(1 + exp(-(1 - 2*u))) # Pr(z|u).
z <- rbinom(N, 1, prob = p_z)
lambda_base <- 0.1 # Baseline hazard.
beta_z <- log(0.5) # log(HR)
beta_u <- log(2)
lambda <- lambda_base*exp(beta_z*z + beta_u*u)
time <- rexp(N, rate = lambda)

# Apply administrative censoring:
study_length <- 30
delta <- ifelse(time < study_length, 1, 0)
time[delta == 0] <- study_length

coxph(Surv(time, delta) ~ z) # Gives a biased estimate!!!

# Adjust for confounding:
coxph(Surv(time, delta) ~ z + u) # Recovers true parameter values.

# Hyperparameters:
mu_log_lambda <- 0 
tau_log_lambda <- 0.01

mu_beta_z <- 0
tau_beta_z <- 0.01

mu_beta_u <- 0
tau_beta_u <- 0.01

a_p <- 1
b_p <- 1

# Time horizon for rMST:
tau <- 6

data_list_confounding <- list(N = N, z = z, time = time, 
                              delta = delta, mu_log_lambda = mu_log_lambda, 
                              tau_log_lambda = tau_log_lambda, 
                              mu_beta_z = mu_beta_z, 
                              tau_beta_z = tau_beta_z, 
                              mu_beta_u = mu_beta_u, tau_beta_u = tau_beta_u, 
                              a_p = a_p, b_p = b_p, 
                              tau = tau)

model_confounding <- jags.model(file = textConnection(model_string_confounding), 
                                data = data_list_confounding, n.chains = 1)

param_confounding <- c("lambda_base", "HR", "beta_u", "p") 
posterior_sample_confounding <- coda.samples(model_confounding, variable.names = param_confounding,
                                             n.iter = N_iter)[[1]]

######################
######################
# Examine posterior
# sample:
######################
######################

# HR:
hist(posterior_sample_confounding[-c(burn_in:N_iter), 1])

# Baseline hazard:
hist(posterior_sample_confounding[-c(burn_in:N_iter), 3])

# Effect of u on survival:
hist(posterior_sample_confounding[-c(burn_in:N_iter), 2])

# Pr(u = 1):
hist(posterior_sample_confounding[-c(burn_in:N_iter), 4])


##################################################


model_confounding <- jags.model(file = textConnection(model), 
                                data = data_list_confounding, n.chains = 1)

param_confounding <- c("lambda_base", "HR", "beta_u", "p", "alpha") 
posterior_sample_confounding <- coda.samples(model_confounding, variable.names = param_confounding,
                                             n.iter = N_iter)[[1]]

# alpha:
posterior_sample_confounding[-c(burn_in:N_iter), 1]

posterior_sample_confounding[-c(burn_in:N_iter), 2]

# HR:
posterior_sample_confounding[-c(burn_in:N_iter), 3]

# beta_u:
posterior_sample_confounding[-c(burn_in:N_iter), 4]

# lambda_base:
posterior_sample_confounding[-c(burn_in:N_iter), 5]

# Pr(u = 1):
hist(posterior_sample_confounding[-c(burn_in:N_iter), 6])
