## ---------------------------
##
## Script name: JAGS_code.R
##
## Purpose of script: JAGS code for QBA work.
##
## Author: Steve Soutar
##
## Date Created: 2024-02-07
##
## Email: Steven.Soutar@arcturis.data.com
##
## ---------------------------
##
## Notes: Two model specifications are presented here. Each model
##        specifies a proportional hazards model with exponential
##        baseline. However, they differ in the inclusion of a binary 
##        confounder. Interest lies in posterior estimation
##        of restricted mean survival for both arms of the trial.
##   
## ---------------------------  


##########################
##########################
# Model specification
# without confounder:
##########################
##########################

model_string_no_confounding <- "
  
  ###################
  # Data block:
  ###################
  
  data{
    
    C <- 1000000
    
    for(i in 1:N){
    
      zeros[i] <- 0
    
    }
  
  }
  
  ###############
  # Model block:
  ###############
  
  model{
    
    ######################
    # Prior specification:
    ######################
    
    # Baseline hazard (log scale):
    log_lambda_base ~ dnorm(mu_log_lambda, tau_log_lambda)
    lambda_base <- exp(log_lambda_base)
    
    # Treatment effect (log scale):
    beta_z ~ dnorm(mu_beta_z, tau_beta_z)
    
    # Hazard ratio:
    HR <- exp(beta_z)
    
    ###########################
    # Likelihood specification:
    ###########################
    for(i in 1:N){
      
      # Hazard for patient i (log scale):
      log_lambda[i] <- log_lambda_base + beta_z*z[i] 
      
      lambda[i] <- exp(log_lambda[i])
      
      # Cumulative hazard for patient i:
      H[i] <- lambda[i]*time[i]
        
      # log-survival for patient i:
      log_surv[i] <- -H[i]
      
      # Zeros trick:
      phi[i] <- C - delta[i]*log_lambda[i] - log_surv[i] # Phi must be kept positive!!
      zeros[i] ~ dpois(phi[i])
    
    }
  
    #########################
    # Compute restricted mean 
    # survival (rMST):
    # (Restricted to tau).
    #########################
    
    # ctrl arm:
    
    rMST_ctrl <- 1/lambda_base  -(1/lambda_base)*exp(-lambda_base*tau) 
    
    # trt arm:
    
    rMST_trt <- 1/(lambda_base*exp(beta_z)) - 1/(lambda_base*exp(beta_z))*exp(-lambda_base*exp(beta_z)*tau)
    
  }

"

##########################
##########################
# Model specification
# with binary confounder:
##########################
##########################

model_string_confounding <- "
  
  ###################
  # Data block:
  ###################
  
  data{
  
    C <- 1000000
    for(i in 1:N){
    
      zeros[i] <- 0
    
    }
  
  }
  
  ###############
  # Model block:
  ###############
  
  model{
    
    ######################
    # Prior specification:
    ######################
    
    # Baseline hazard (log scale):
    log_lambda_base ~ dnorm(mu_log_lambda, tau_log_lambda) # dnorm is parametrised in terms of precision!
    lambda_base <- exp(log_lambda_base)
    
    # Treatment effect (log scale):
    beta_z ~ dnorm(mu_beta_z, tau_beta_z)
    
    # Hazard ratio:
    HR <- exp(beta_z)
    
    # Binary confounder:
    for(i in 1:N){
    
      u[i] ~ dbern(p)
    
    }
    
    # Effect of u on survival (log scale):
    beta_u ~ dnorm(mu_beta_u, tau_beta_u)
    
    # Population probability of u = 1:
    p ~ dbeta(a_p, b_p) # p ~ U(0, 1) when a_p = b_p = 1.
    
    ###########################
    # Likelihood specification:
    ###########################
    
    for(i in 1:N){
      
      # Hazard for patient i (log scale):
      log_lambda[i] <- log_lambda_base + beta_z*z[i] + beta_u*u[i]
      
      lambda[i] <- exp(log_lambda[i])
      
      # Cumulative hazard for patient i:
      H[i] <- lambda[i]*time[i]
        
      # log-survival for patient i:
      log_surv[i] <- -H[i]
      
      # Zeros trick:
      phi[i] <- C - delta[i]*log_lambda[i] - log_surv[i] # Phi must be kept positive!!
      zeros[i] ~ dpois(phi[i])
    
    }
    
    #########################
    # Compute restricted mean 
    # survival (rMST):
    # (Restricted to tau).
    #########################
    
    # ctrl arm:
    
    #rMST_ctrl <- (-1/lambda_base)*exp(-lambda_base*tau) - (-1/lambda_base) 
    
    # trt arm:
    
    #rMST_trt <- (-1/lambda_base*HR)*exp(-lambda_base*HR*tau) - (-1/(lambda_base*exp(beta_z)))
    
  }

"

####################################################################
