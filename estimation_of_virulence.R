#### Estimation of virulence

# Libraries ------------------------------------------------------------------------------------------------------
library(foreign)

# Read data ------------------------------------------------------------------------------------------------------
dataset <- read.spss("/home/petr/Documents/Master_thesis/AgeEffectsVirulence.sav", to.data.frame = TRUE, use.value.labels = FALSE)
parasyte <- as.character(dataset$Pasteuria)
infecteds <- as.numeric(dataset$Infected)

# Subsetting w.r.t. age of infection ---------------------------------------------------------------------------
Infecteds_AaI_5 <- dataset[(dataset$InfectionAge == 5) & (dataset$Infected == 1) ,]
Exposed_AaI_5 <- dataset[(dataset$InfectionAge == 5) & (dataset$Infected == 0) ,]
Infecteds_AaI_15 <- dataset[(dataset$InfectionAge == 15) & (dataset$Infected == 1) ,]
Exposed_AaI_15 <- dataset[(dataset$InfectionAge == 15) & (dataset$Infected == 0) ,]
Infecteds_AaI_30 <- dataset[(dataset$InfectionAge == 30) & (dataset$Infected == 1) ,]
Exposed_AaI_30 <- dataset[(dataset$InfectionAge == 30) & (dataset$Infected == 0) ,]

Infecteds_AaI_longevities <- list(Infecteds_AaI_5$HostLongevity, Infecteds_AaI_15$HostLongevity, Infecteds_AaI_30$HostLongevity)
Exposed_AaI_longevities <- list(Exposed_AaI_5$HostLongevity, Exposed_AaI_15$HostLongevity, Exposed_AaI_30$HostLongevity)

# Choosing the lifetime threshold and initializing age axis ----------------------------------------------------
lifetime_threshold <- 150
age_axis <- 0:lifetime_threshold

# Computing the survival distributions for the subsets of interest ---------------------------------------------
survival_distribution <- function(longevity){
  survival_dist <- rep(0, lifetime_threshold +1)
  for (x in 0:lifetime_threshold) {
    survival_dist[x+1] <- sum(longevity >x)/ length(longevity)
  }
  return(survival_dist)
}


Survivals_Infecteds_AaI <- lapply(X=Infecteds_AaI_longevities, FUN=survival_distribution)
Survivals_Exposed_AaI <- lapply(X=Exposed_AaI_longevities, FUN=survival_distribution)

plot(age_axis, Survivals_Exposed_AaI[[2]])
# Delta Ansatzes to be optimized ------------------------------------------------------------------------------

delta <- function(parameters, age = age_axis){
  local_delta <- 0
  if(length(parameters) == 2){
    local_delta <- parameters[1]*age + parameters[2]
  }else if(length(parameters) == 3){
    local_delta <- parameters[1] * age^2 + parameters[2]*age + parameters[3]
  }else if(length(parameters) == 4){
    local_delta <- parameters[1] * age^3 + parameters[2]*age^2 + parameters[3]*age + parameters[4]
  }else if(length(parameters) == 5){
    local_delta <- parameters[1] * age^4 + parameters[2]*age^3 + parameters[3]*age^2 + parameters[4]*age + parameters[5]
  }else if(length(parameters) == 6){
    local_delta <- parameters[1] * age^5 + parameters[2]*age^4 + parameters[3]*age^3 + parameters[4]*age^2 + parameters[5]*age + parameters[6]
  }else (message("Length of argument is not between 1 and 6, it is:", length(parameters)))
  return (local_delta)
}

# Introducing -log likelihood function --------------------------------------------------------------------------
log_likelihood <- function(parameters, longevity, AaI){
  likelihood_local <- 0
  if(min(delta(parameters, age_axis))>=0){
     for (i in 1:length(longevity)) {
      likelihood_local <- likelihood_local - sum(delta(parameters, age_axis[AaI:longevity[i]])) + log(delta(parameters,longevity[i]))
    }
    
  }else{
    likelihood_local <- -Inf
  }
  
  return(-likelihood_local) #optim is minimizing, so we return sign flipped value to get maximization ## + log(factorial(length(longevity),
                                                                  ##consider adding the binomial terms to get the absolute likelihood,
                                                                  ##but it might differ from the factorial since we have different AaI now 
}

#### Optimization ------------------------------------------------------------------------------------------------

## Initializing variables ----------------------------------------------------------------------------------------
Optimal_likelihood_Infecteds_AaI_linear <- list(Inf,Inf,Inf)
Optimal_likelihood_Exposed_AaI_linear <- list(Inf,Inf,Inf)
Optimal_likelihood_Infecteds_AaI_quadratic <- list(Inf,Inf,Inf)
Optimal_likelihood_Exposed_AaI_quadratic <- list(Inf,Inf,Inf)
Optimal_likelihood_Infecteds_AaI_cubic <- list(Inf,Inf,Inf)
Optimal_likelihood_Exposed_AaI_cubic <- list(Inf,Inf,Inf)
Optimal_likelihood_Infecteds_AaI_quartic <- list(Inf,Inf,Inf)
Optimal_likelihood_Exposed_AaI_quartic <- list(Inf,Inf,Inf)

Optimal_likelihood_Infecteds_AaI_all_models <- list(Optimal_likelihood_Infecteds_AaI_linear, Optimal_likelihood_Infecteds_AaI_quadratic, Optimal_likelihood_Infecteds_AaI_cubic, Optimal_likelihood_Infecteds_AaI_quartic)
Optimal_likelihood_Exposed_AaI_all_models <- list(Optimal_likelihood_Exposed_AaI_linear, Optimal_likelihood_Exposed_AaI_quadratic, Optimal_likelihood_Exposed_AaI_cubic, Optimal_likelihood_Exposed_AaI_quartic)
Optimal_parameters_Infecteds_AaI_all_models <- list(list(c(0,0),c(0,0),c(0,0)),list(c(0,0,0),c(0,0,0),c(0,0,0)),list(c(0,0,0,0),c(0,0,0,0),c(0,0,0,0)), list(c(0,0,0,0,0),c(0,0,0,0,0),c(0,0,0,0,0)))
Optimal_parameters_Exposed_AaI_all_models <- list(list(c(0,0),c(0,0),c(0,0)),list(c(0,0,0),c(0,0,0),c(0,0,0)),list(c(0,0,0,0),c(0,0,0,0),c(0,0,0,0)), list(c(0,0,0,0,0),c(0,0,0,0,0),c(0,0,0,0,0)))

str(Optimal_parameters_Infecteds_AaI_all_models)
## Initial guesses for the parameters in the polynomial for delta -------------------------------------------------
##const_initial = c(0.00001)
##lin_initial = c(0.00001)
##quadrat_initial = c(0.000001)
##cubic_initial = c(0.00000001)
##quartic_initial = c(0.000000001)
const_initial <- c(0.1,0.001,0.00001)
lin_initial <- c(0.1,0.001,0.00001)
quadrat_initial <- c(0.01,0.0001,0.000001)
cubic_initial <- c(0.0001,0.000001,0.00000001)
quartic_initial <- c(0.00001,0.0000001,0.000000001) #maybe try even with the parameters from the original optimization (at least the orders of magnitude)

## Function to be utilized during the iteration over the different initial parameters -----------------------------
optim_subroutine <- function(parameters, Function=log_likelihood, infection, AaI_index){
  model_index <- (length(parameters)-1)
  Age_at_infection <- if (AaI_index == 1) 5 else if(AaI_index == 2) 15 else if(AaI_index == 3) 30 else 999
  if(infection == "I"){
    if(optim(par = parameters, fn=Function, longevity = Infecteds_AaI_longevities[[AaI_index]], AaI = Age_at_infection)$value < Optimal_likelihood_Infecteds_AaI_all_models[[model_index]][[AaI_index]]){
      return(optim(par = parameters, fn=Function, longevity = Infecteds_AaI_longevities[[AaI_index]], AaI = Age_at_infection))
      #message("Infecteds: Correction achieved for initial set of parameters: ", parameters)
    }
    
  }else if(infection == "E"){
    if(optim(par = parameters, fn=Function, AaI = Age_at_infection, longevity = Exposed_AaI_longevities[[AaI_index]])$value < Optimal_likelihood_Exposed_AaI_all_models[[model_index]][[AaI_index]]){
      return(optim(par = parameters, fn=Function, longevity = Exposed_AaI_longevities[[AaI_index]], AaI = Age_at_infection))
    }
      #message("Exposed: Correction achieved for initial set of parameters: ", parameters)
  }else{
    message("Invalid indication of the compartment Infected/Exposed through the parameter infection ")
  }
}

counter <-0
model_parameters <- 0
# Iteration over different initial parameters ----------------------------------------------------------------------------
for (i in const_initial) {
  
  message("i iteration began")
  for (j in lin_initial){
    message("j iteration began")
    counter <- counter + 1
    message("-------------------------------------------------------------",counter, "/9","-------------------------------------------------------------")
    model_parameters <- c(i,j)
    
    Optimal_parameters_Infecteds_AaI_all_models[[length(model_parameters)-1]][1] <- optim_subroutine(parameters = model_parameters, AaI_index = 1, infection = "I")$par
    Optimal_parameters_Infecteds_AaI_all_models[[length(model_parameters)-1]][2] <- optim_subroutine(parameters = model_parameters, AaI_index = 2, infection = "I")$par
    Optimal_parameters_Infecteds_AaI_all_models[[length(model_parameters)-1]][3] <- optim_subroutine(parameters = model_parameters, AaI_index = 3, infection = "I")$par
    
    Optimal_parameters_Exposed_AaI_all_models[[length(model_parameters)-1]][1] <- optim_subroutine(parameters = model_parameters, AaI_index = 1, infection = "E")$par
    Optimal_parameters_Exposed_AaI_all_models[[length(model_parameters)-1]][2] <- optim_subroutine(parameters = model_parameters, AaI_index = 2, infection = "E")$par
    Optimal_parameters_Exposed_AaI_all_models[[length(model_parameters)-1]][3] <- optim_subroutine(parameters = model_parameters, AaI_index = 3, infection = "E")$par
    
    Optimal_likelihood_Infecteds_AaI_all_models[[length(model_parameters)-1]][1] <- optim_subroutine(parameters = model_parameters, AaI_index = 1, infection = "I")$value
    Optimal_likelihood_Infecteds_AaI_all_models[[length(model_parameters)-1]][2] <- optim_subroutine(parameters = model_parameters, AaI_index = 2, infection = "I")$value
    Optimal_likelihood_Infecteds_AaI_all_models[[length(model_parameters)-1]][3] <- optim_subroutine(parameters = model_parameters, AaI_index = 3, infection = "I")$value
    
    Optimal_likelihood_Exposed_AaI_all_models[[length(model_parameters)-1]][1] <- optim_subroutine(parameters = model_parameters, AaI_index = 1, infection = "E")$value
    Optimal_likelihood_Exposed_AaI_all_models[[length(model_parameters)-1]][2] <- optim_subroutine(parameters = model_parameters, AaI_index = 2, infection = "E")$value
    Optimal_likelihood_Exposed_AaI_all_models[[length(model_parameters)-1]][3] <- optim_subroutine(parameters = model_parameters, AaI_index = 3, infection = "E")$value
    
    for (k in quadrat_initial) {
      
      message("k iteration began")
      model_parameters <- c(i,j,k)
      
      Optimal_parameters_Infecteds_AaI_all_models[[length(model_parameters)-1]][1] <- optim_subroutine(parameters = model_parameters, AaI_index = 1, infection = "I")$par
      Optimal_parameters_Infecteds_AaI_all_models[[length(model_parameters)-1]][2] <- optim_subroutine(parameters = model_parameters, AaI_index = 2, infection = "I")$par
      Optimal_parameters_Infecteds_AaI_all_models[[length(model_parameters)-1]][3] <- optim_subroutine(parameters = model_parameters, AaI_index = 3, infection = "I")$par
      
      Optimal_parameters_Exposed_AaI_all_models[[length(model_parameters)-1]][1] <- optim_subroutine(parameters = model_parameters, AaI_index = 1, infection = "E")$par
      Optimal_parameters_Exposed_AaI_all_models[[length(model_parameters)-1]][2] <- optim_subroutine(parameters = model_parameters, AaI_index = 2, infection = "E")$par
      Optimal_parameters_Exposed_AaI_all_models[[length(model_parameters)-1]][3] <- optim_subroutine(parameters = model_parameters, AaI_index = 3, infection = "E")$par
      
      Optimal_likelihood_Infecteds_AaI_all_models[[length(model_parameters)-1]][1] <- optim_subroutine(parameters = model_parameters, AaI_index = 1, infection = "I")$value
      Optimal_likelihood_Infecteds_AaI_all_models[[length(model_parameters)-1]][2] <- optim_subroutine(parameters = model_parameters, AaI_index = 2, infection = "I")$value
      Optimal_likelihood_Infecteds_AaI_all_models[[length(model_parameters)-1]][3] <- optim_subroutine(parameters = model_parameters, AaI_index = 3, infection = "I")$value
      
      Optimal_likelihood_Exposed_AaI_all_models[[length(model_parameters)-1]][1] <- optim_subroutine(parameters = model_parameters, AaI_index = 1, infection = "E")$value
      Optimal_likelihood_Exposed_AaI_all_models[[length(model_parameters)-1]][2] <- optim_subroutine(parameters = model_parameters, AaI_index = 2, infection = "E")$value
      Optimal_likelihood_Exposed_AaI_all_models[[length(model_parameters)-1]][3] <- optim_subroutine(parameters = model_parameters, AaI_index = 3, infection = "E")$value
      
      for (l in cubic_initial) {
        
          message("l iteration began")
          model_parameters <- c(i,j,k,l)
          
          Optimal_parameters_Infecteds_AaI_all_models[[length(model_parameters)-1]][1] <- optim_subroutine(parameters = model_parameters, AaI_index = 1, infection = "I")$par
          Optimal_parameters_Infecteds_AaI_all_models[[length(model_parameters)-1]][2] <- optim_subroutine(parameters = model_parameters, AaI_index = 2, infection = "I")$par
          Optimal_parameters_Infecteds_AaI_all_models[[length(model_parameters)-1]][3] <- optim_subroutine(parameters = model_parameters, AaI_index = 3, infection = "I")$par
          
          Optimal_parameters_Exposed_AaI_all_models[[length(model_parameters)-1]][1] <- optim_subroutine(parameters = model_parameters, AaI_index = 1, infection = "E")$par
          Optimal_parameters_Exposed_AaI_all_models[[length(model_parameters)-1]][2] <- optim_subroutine(parameters = model_parameters, AaI_index = 2, infection = "E")$par
          Optimal_parameters_Exposed_AaI_all_models[[length(model_parameters)-1]][3] <- optim_subroutine(parameters = model_parameters, AaI_index = 3, infection = "E")$par
          
          Optimal_likelihood_Infecteds_AaI_all_models[[length(model_parameters)-1]][1] <- optim_subroutine(parameters = model_parameters, AaI_index = 1, infection = "I")$value
          Optimal_likelihood_Infecteds_AaI_all_models[[length(model_parameters)-1]][2] <- optim_subroutine(parameters = model_parameters, AaI_index = 2, infection = "I")$value
          Optimal_likelihood_Infecteds_AaI_all_models[[length(model_parameters)-1]][3] <- optim_subroutine(parameters = model_parameters, AaI_index = 3, infection = "I")$value
          
          Optimal_likelihood_Exposed_AaI_all_models[[length(model_parameters)-1]][1] <- optim_subroutine(parameters = model_parameters, AaI_index = 1, infection = "E")$value
          Optimal_likelihood_Exposed_AaI_all_models[[length(model_parameters)-1]][2] <- optim_subroutine(parameters = model_parameters, AaI_index = 2, infection = "E")$value
          Optimal_likelihood_Exposed_AaI_all_models[[length(model_parameters)-1]][3] <- optim_subroutine(parameters = model_parameters, AaI_index = 3, infection = "E")$value
          
          for (m in quartic_initial) {
            
            message("m iteration began")
            # message("Structure of the variable at the beggining of the critical loop")
            # str(Optimal_likelihood_Infecteds_AaI_all_models)
            model_parameters <- c(i,j,k,l,m)
            Optimal_parameters_Infecteds_AaI_all_models[[length(model_parameters)-1]][1] <- optim_subroutine(parameters = model_parameters, AaI_index = 1, infection = "I")$par
            Optimal_parameters_Infecteds_AaI_all_models[[length(model_parameters)-1]][2] <- optim_subroutine(parameters = model_parameters, AaI_index = 2, infection = "I")$par
            Optimal_parameters_Infecteds_AaI_all_models[[length(model_parameters)-1]][3] <- optim_subroutine(parameters = model_parameters, AaI_index = 3, infection = "I")$par
            
            Optimal_parameters_Exposed_AaI_all_models[[length(model_parameters)-1]][1] <- optim_subroutine(parameters = model_parameters, AaI_index = 1, infection = "E")$par
            Optimal_parameters_Exposed_AaI_all_models[[length(model_parameters)-1]][2] <- optim_subroutine(parameters = model_parameters, AaI_index = 2, infection = "E")$par
            Optimal_parameters_Exposed_AaI_all_models[[length(model_parameters)-1]][3] <- optim_subroutine(parameters = model_parameters, AaI_index = 3, infection = "E")$par
            
            message("Structure of the variable before the critical paragraph, individual element")
            str(Optimal_likelihood_Infecteds_AaI_all_models[[4]][[3]])
            message("length(model_parameters)-1 = ",length(model_parameters)-1)
            message("optim(...)$Value = ",optim_subroutine(parameters = model_parameters, AaI_index = 3, infection = "I")$value)
            message("structure od optim()$value")
            str(optim_subroutine(parameters = model_parameters, AaI_index = 3, infection = "I")$value)
            
            
            Optimal_likelihood_Infecteds_AaI_all_models[[length(model_parameters)-1]][1] <- optim_subroutine(parameters = model_parameters, AaI_index = 1, infection = "I")$value
            Optimal_likelihood_Infecteds_AaI_all_models[[length(model_parameters)-1]][2] <- optim_subroutine(parameters = model_parameters, AaI_index = 2, infection = "I")$value
            #Optimal_likelihood_Infecteds_AaI_all_models[[length(model_parameters)-1]][3] <- optim_subroutine(parameters = model_parameters, AaI_index = 3, infection = "I")$value
            message("Structure of the variable after the critical paragraph")
            str(Optimal_likelihood_Infecteds_AaI_all_models)
            
            Optimal_likelihood_Exposed_AaI_all_models[[length(model_parameters)-1]][1] <- optim_subroutine(parameters = model_parameters, AaI_index = 1, infection = "E")$value
            Optimal_likelihood_Exposed_AaI_all_models[[length(model_parameters)-1]][2] <- optim_subroutine(parameters = model_parameters, AaI_index = 2, infection = "E")$value
            Optimal_likelihood_Exposed_AaI_all_models[[length(model_parameters)-1]][3] <- optim_subroutine(parameters = model_parameters, AaI_index = 3, infection = "E")$value
            message("4th paragraph ok")
        }
      }
    }
  }
}
