## ---------------------------
##
## Script name: simple_JAGS_ex.R
##
## Purpose of script: Simple example illustrating the use of JAGS for
##                    Bayesian inference.
##
## Author: Steve Soutar
##
## Date Created: 2024-02-13
##
## Email: Steven.Soutar@arcturis.data.com
##
## ---------------------------
##
## Notes: Inference for normal errors model.
##   
##
## ---------------------------  


model_string <- "

  ##################
  # Model block
  ##################
  
  model{
    
    #####################
    # Prior specification
    #####################

    mu ~ dnorm(0, 0.01)
    tau ~ dgamma(1, 1)
  
    # Derived variance parameter:
    sigma_sq <- 1/tau
    
    ###########################
    # Likelihood specification:
    ###########################
    
    for(i in 1:N){
  
      y[i] ~ dnorm(mu, tau) # dnorm is parametrised in terms of precision!
  
    }
  
  }

"
#########################
#########################
# Test on simulated data:
#########################
#########################

N <- 1000
mu <- 2
sigma_sq <- 0.5
y <- rnorm(N, mu, sqrt(sigma_sq))


# List of observed data:
data_list <- list(N = N, y = y)

# Compile model:
model <- jags.model(file = textConnection(model_string), 
           data = data_list, n.chains = 1)


param <- c("mu", "sigma_sq") # Parameters of interest.

N_iter <- 20000 # Total number of samples.
N_burn_in <- 10000 # Burn-in period.

# Extract posterior sample:
posterior_sample <- coda.samples(model, 
                      variable.names = param,
                      n.iter = N_iter)[[1]]

# mu:
hist(posterior_sample[-c(burn_in:N_iter), 1])
# sigma_sq:
hist(posterior_sample[-c(burn_in:N_iter), 2])

