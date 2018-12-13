#Libraries
library(foreign)   #to read the spss file

#Daphnia cannot live longer than 200 days (in our data not even longer than 123)
lifetime_threshold = 150

#Initialize the variables ---------------------------------------------------------------
age_axis <- c(0:lifetime_threshold)
age_of_inf_5_axis <- c(0: (lifetime_threshold-5))
age_of_inf_15_axis <- c(0: (lifetime_threshold-15))
age_of_inf_30_axis <- c(0: (lifetime_threshold-30))

#Our polynomial Ansatz for delta---------------------------------------------------------
delta <- function(parameters, age){
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
  }
 return (local_delta)
}

#Reading the data -----------------------------------------------------------------------
dataset <- read.spss("/home/petr/Documents/Master_thesis/AgeEffectsVirulence.sav", to.data.frame = TRUE, use.value.labels = FALSE)
parasyte <- as.character(dataset$Pasteuria)
infecteds <- as.numeric(dataset$Infected)
aai <- as.numeric(dataset$InfectionAge)

#Subsetting --------------------------------

control_population <- dataset[parasyte == "Control ",] # Careful, the value of the variable contains the space after Control, possibly FIX the dataset later
host_longevity_control <- control_population$HostLongevity

infected_population <- dataset[infecteds == 1,]
longevity_aai_5_inf <- dataset[infecteds == 1 & aai == 5,]$HostLongevity
longevity_aai_15_inf <- dataset[infecteds == 1 & aai == 15,]$HostLongevity
longevity_aai_30_inf <- dataset[infecteds == 1 & aai == 30,]$HostLongevity

exposed_population <- dataset[(infecteds == 0) & (parasyte != "Control "),]
longevity_aai_5_exp <- dataset[infecteds == 0 & aai == 5,]$HostLongevity
longevity_aai_15_exp <- dataset[infecteds == 0 & aai == 15,]$HostLongevity
longevity_aai_30_exp <- dataset[infecteds == 0 & aai == 30,]$HostLongevity

longevity_uninfected <- c(host_longevity_control, longevity_aai_5_exp, longevity_aai_15_exp, longevity_aai_30_exp)
longevity_exposed <- c(longevity_aai_5_exp, longevity_aai_15_exp, longevity_aai_30_exp)

#Survivals of the control, infected and exposed populations ---------------------------------------
survivals_aai_5_inf <- rep(0,lifetime_threshold+1)
survivals_aai_15_inf <- rep(0,lifetime_threshold+1)
survivals_aai_30_inf <- rep(0,lifetime_threshold+1)
survivals_aai_5_exp <- rep(0,lifetime_threshold+1)
survivals_aai_15_exp <- rep(0,lifetime_threshold+1)
survivals_aai_30_exp <- rep(0,lifetime_threshold+1)
survivals_aai_5_control <- rep(0,lifetime_threshold+1)
survivals_aai_15_control <- rep(0,lifetime_threshold+1)
survivals_aai_30_control <- rep(0,lifetime_threshold+1)
survivals_uninfected <- rep(0,lifetime_threshold+1)
survivals_exposed <- rep(0,lifetime_threshold+1)

for (x in 0:lifetime_threshold) {
  survivals_aai_5_inf[x+1] <- sum(longevity_aai_5_inf>x)/ length(longevity_aai_5_inf)
  survivals_aai_15_inf[x+1] <- sum(longevity_aai_15_inf>x)/ length(longevity_aai_15_inf)
  survivals_aai_30_inf[x+1] <- sum(longevity_aai_30_inf>x)/ length(longevity_aai_30_inf)
  survivals_aai_5_exp[x+1] <- sum(longevity_aai_5_exp>x)/ length(longevity_aai_5_exp)
  survivals_aai_15_exp[x+1] <- sum(longevity_aai_15_exp>x)/ length(longevity_aai_15_exp)
  survivals_aai_30_exp[x+1] <- sum(longevity_aai_30_exp>x)/ length(longevity_aai_30_exp)
  survivals_uninfected[x+1] <- sum(longevity_uninfected>x)/ length(longevity_uninfected)
  survivals_exposed[x+1] <- sum(longevity_exposed>x)/length(longevity_exposed)
}

for (x in 4:lifetime_threshold){
  survivals_aai_5_control[x+1] <- sum(host_longevity_control>x)/ sum(host_longevity_control>5)
}

for (x in 14:lifetime_threshold){
  survivals_aai_15_control[x+1] <- sum(host_longevity_control>x)/ sum(host_longevity_control>15)
}

for (x in 29:lifetime_threshold) {
  survivals_aai_30_control[x+1] <- sum(host_longevity_control>x)/ sum(host_longevity_control>30)
}




# log likelihood------
#Likelihood function (log_likelihood) = log(Pruduct over individuals (Product over timesteps(probability of dying in that timestep))) -> log turns pruducts to sums
likelihood_general <- function(parameters, longevity, age_at_infection){
  likelihood_local <- 0
  if(min(delta(parameters, age_axis))>=0){
      #does the binomial factor play any role in the likelihood? It does not, right? Since under the logarithm it only plays role of additive constant (independent of delta)
      for (i in 1:length(longevity)) {
        likelihood_local <- likelihood_local -2*sum(delta(parameters, age_axis[age_at_infection:(longevity[i]-1)])) + log(1-exp(-delta(parameters,longevity[i]))) #instead of delta put in the right probability of dying at age longevity[i]
      }
    
    }else{
      likelihood_local <- -Inf
    }
  
  return(-likelihood_local) #optim is minimizing, so we return sign flipped value to get maximization #consider adding the binomial factor to have absolute likelihood,
                                                                                                      #keep in mind the different ages of infections which play role
}


######################## Optimization -------------------------------------------------------------------

const_initial = c(0.1,0.001,0.00001)
lin_initial = c(0.1,0.001,0.00001)
quadrat_initial = c(0.01,0.0001,0.000001)
cubic_initial = c(0.0001,0.000001,0.00000001)
quartic_initial = c(0.00001,0.0000001,0.000000001)

optimized_likelihood_aai_5_inf_linear <- Inf
optimized_likelihood_aai_15_inf_linear <- Inf
optimized_likelihood_aai_30_inf_linear <- Inf

optimized_likelihood_aai_5_exp_linear <- Inf
optimized_likelihood_aai_15_exp_linear <- Inf
optimized_likelihood_aai_30_exp_linear <- Inf

optimized_likelihood_aai_5_control_linear <- Inf
optimized_likelihood_aai_15_control_linear <- Inf
optimized_likelihood_aai_30_control_linear <- Inf

optimized_likelihood_uninfected_linear <- Inf
optimized_likelihood_exposed_linear <- Inf

optimized_likelihood_aai_5_inf_quadratic <- Inf
optimized_likelihood_aai_15_inf_quadratic<- Inf
optimized_likelihood_aai_30_inf_quadratic <- Inf

optimized_likelihood_aai_5_exp_quadratic <- Inf
optimized_likelihood_aai_15_exp_quadratic <- Inf
optimized_likelihood_aai_30_exp_quadratic <- Inf

optimized_likelihood_aai_5_control_quadratic <- Inf
optimized_likelihood_aai_15_control_quadratic <- Inf
optimized_likelihood_aai_30_control_quadratic <- Inf

optimized_likelihood_uninfected_quadratic <- Inf
optimized_likelihood_exposed_quadratic <- Inf

optimized_likelihood_aai_5_inf_cubic <- Inf
optimized_likelihood_aai_15_inf_cubic <- Inf
optimized_likelihood_aai_30_inf_cubic <- Inf

optimized_likelihood_aai_5_exp_cubic <- Inf
optimized_likelihood_aai_15_exp_cubic <- Inf
optimized_likelihood_aai_30_exp_cubic <- Inf

optimized_likelihood_aai_5_control_cubic <- Inf
optimized_likelihood_aai_15_control_cubic <- Inf
optimized_likelihood_aai_30_control_cubic <- Inf

optimized_likelihood_uninfected_cubic <- Inf
optimized_likelihood_exposed_cubic <- Inf

optimized_likelihood_aai_5_inf_quartic <- Inf
optimized_likelihood_aai_15_inf_quartic<- Inf
optimized_likelihood_aai_30_inf_quartic <- Inf

optimized_likelihood_aai_5_exp_quartic <- Inf
optimized_likelihood_aai_15_exp_quartic <- Inf
optimized_likelihood_aai_30_exp_quartic <- Inf

optimized_likelihood_aai_5_control_quartic <- Inf
optimized_likelihood_aai_15_control_quartic <- Inf
optimized_likelihood_aai_30_control_quartic <- Inf

optimized_likelihood_uninfected_quartic <- Inf
optimized_likelihood_exposed_quartic <- Inf

counter <- 0

for (i in const_initial) {
  
  for (j in lin_initial){
    counter <- counter + 1
    message("-------------------------------------------------------------",counter, "/9","-------------------------------------------------------------")
    
    if(optim(par = c(j, i), fn=likelihood_general, longevity = longevity_aai_5_inf, age_at_infection = 5)$value < optimized_likelihood_aai_5_inf_linear){
      optimized_likelihood_aai_5_inf_linear <- optim(par = c(j, i), fn=likelihood_general, longevity = longevity_aai_5_inf, age_at_infection = 5)$value
      optimized_parameters_aai_5_inf_linear <- optim(par = c(j, i), fn=likelihood_general, longevity = longevity_aai_5_inf, age_at_infection = 5)$par
      message("correction achieved for aai_5_inf population and i=",i," and j=", j)
    }

    if(optim(par = c(j, i), fn=likelihood_general, longevity = longevity_aai_15_inf, age_at_infection = 15)$value < optimized_likelihood_aai_15_inf_linear){
      optimized_likelihood_aai_15_inf_linear <- optim(par = c(j, i), fn=likelihood_general, longevity = longevity_aai_15_inf, age_at_infection = 15)$value
      optimized_parameters_aai_15_inf_linear <- optim(par = c(j, i), fn=likelihood_general, longevity = longevity_aai_15_inf, age_at_infection = 15)$par
      message("correction achieved for aai_15_inf population and i=",i," and j=", j)
    }

    if(optim(par = c(j, i), fn=likelihood_general, longevity = longevity_aai_30_inf, age_at_infection = 30)$value < optimized_likelihood_aai_30_inf_linear){
      optimized_likelihood_aai_30_inf_linear <- optim(par = c(j, i), fn=likelihood_general, longevity = longevity_aai_30_inf, age_at_infection = 30)$value
      optimized_parameters_aai_30_inf_linear <- optim(par = c(j, i), fn=likelihood_general, longevity = longevity_aai_30_inf, age_at_infection = 30)$par
      message("correction achieved for aai_30_inf population and i=",i," and j=", j)
    }

    if(optim(par = c(j, i), fn=likelihood_general, longevity = longevity_aai_5_exp, age_at_infection = 5)$value < optimized_likelihood_aai_5_exp_linear){
      optimized_likelihood_aai_5_exp_linear <- optim(par = c(j, i), fn=likelihood_general, longevity = longevity_aai_5_exp, age_at_infection = 5)$value
      optimized_parameters_aai_5_exp_linear <- optim(par = c(j, i), fn=likelihood_general, longevity = longevity_aai_5_exp, age_at_infection = 5)$par
      message("correction achieved for aai_5_exp population and i=",i," and j=", j)
    }

    if(optim(par = c(j, i), fn=likelihood_general, longevity = longevity_aai_15_exp, age_at_infection = 15)$value < optimized_likelihood_aai_15_exp_linear){
      optimized_likelihood_aai_15_exp_linear <- optim(par = c(j, i), fn=likelihood_general, longevity = longevity_aai_15_exp, age_at_infection = 15)$value
      optimized_parameters_aai_15_exp_linear <- optim(par = c(j, i), fn=likelihood_general, longevity = longevity_aai_15_exp, age_at_infection = 15)$par
      message("correction achieved for aai_15_exp population and i=",i," and j=", j)
    }

    if(optim(par = c(j, i), fn=likelihood_general, longevity = longevity_aai_30_exp, age_at_infection = 30)$value < optimized_likelihood_aai_30_exp_linear){
      optimized_likelihood_aai_30_exp_linear <- optim(par = c(j, i), fn=likelihood_general, longevity = longevity_aai_30_exp, age_at_infection = 30)$value
      optimized_parameters_aai_30_exp_linear <- optim(par = c(j, i), fn=likelihood_general, longevity = longevity_aai_30_exp, age_at_infection = 30)$par
      message("correction achieved for aai_30_exp population and i=",i," and j=", j)
    }
    ###control
    if(optim(par = c(j, i), fn=likelihood_general, longevity = host_longevity_control, age_at_infection = 5)$value < optimized_likelihood_aai_5_control_linear){
      optimized_likelihood_aai_5_control_linear <- optim(par = c(j, i), fn=likelihood_general, longevity = host_longevity_control, age_at_infection = 5)$value
      optimized_parameters_aai_5_control_linear <- optim(par = c(j, i), fn=likelihood_general, longevity = host_longevity_control, age_at_infection = 5)$par
      message("correction achieved for aai_5_control population and i=",i," and j=", j)
    }

    if(optim(par = c(j, i), fn=likelihood_general, longevity = host_longevity_control, age_at_infection = 15)$value < optimized_likelihood_aai_15_control_linear){
      optimized_likelihood_aai_15_control_linear <- optim(par = c(j, i), fn=likelihood_general, longevity = host_longevity_control, age_at_infection = 15)$value
      optimized_parameters_aai_15_control_linear <- optim(par = c(j, i), fn=likelihood_general, longevity = host_longevity_control, age_at_infection = 15)$par
      message("correction achieved for aai_15_control population and i=",i," and j=", j)
    }

    if(optim(par = c(j, i), fn=likelihood_general, longevity = host_longevity_control, age_at_infection = 30)$value < optimized_likelihood_aai_30_control_linear){
      optimized_likelihood_aai_30_control_linear <- optim(par = c(j, i), fn=likelihood_general, longevity = host_longevity_control, age_at_infection = 30)$value
      optimized_parameters_aai_30_control_linear <- optim(par = c(j, i), fn=likelihood_general, longevity = host_longevity_control, age_at_infection = 30)$par
      message("correction achieved for aai_30_control population and i=",i," and j=", j)
    }
    
    if(optim(par = c(j, i), fn=likelihood_general, longevity = longevity_uninfected, age_at_infection = 16)$value < optimized_likelihood_uninfected_linear){
      optimized_likelihood_uninfected_linear <- optim(par = c(j, i), fn=likelihood_general, longevity = longevity_uninfected, age_at_infection = 16)$value
      optimized_parameters_uninfected_linear <- optim(par = c(j, i), fn=likelihood_general, longevity = longevity_uninfected, age_at_infection = 16)$par
      message("correction achieved for uninfected population and i=",i," and j=", j)
    }
    
    if(optim(par = c(j, i), fn=likelihood_general, longevity = longevity_exposed, age_at_infection = 16)$value < optimized_likelihood_exposed_linear){
      optimized_likelihood_exposed_linear <- optim(par = c(j, i), fn=likelihood_general, longevity = longevity_exposed, age_at_infection = 16)$value
      optimized_parameters_exposed_linear <- optim(par = c(j, i), fn=likelihood_general, longevity = longevity_exposed, age_at_infection = 16)$par
      message("correction achieved for uninfected population and i=",i," and j=", j)
    }
    
    for (k in quadrat_initial) {
      
      if(optim(par = c(k, j, i), fn=likelihood_general, longevity = longevity_aai_5_inf, age_at_infection = 5)$value < optimized_likelihood_aai_5_inf_quadratic){
        optimized_likelihood_aai_5_inf_quadratic <- optim(par = c(k, j, i), fn=likelihood_general, longevity = longevity_aai_5_inf, age_at_infection = 5)$value
        optimized_parameters_aai_5_inf_quadratic <- optim(par = c(k, j, i), fn=likelihood_general, longevity = longevity_aai_5_inf, age_at_infection = 5)$par
        message("correction achieved for aai_5_inf population and i=",i," and j=", j, " and k=", k)
      }

      if(optim(par = c(k, j, i), fn=likelihood_general, longevity = longevity_aai_15_inf, age_at_infection = 15)$value < optimized_likelihood_aai_15_inf_quadratic){
        optimized_likelihood_aai_15_inf_quadratic <- optim(par = c(k, j, i), fn=likelihood_general, longevity = longevity_aai_15_inf, age_at_infection = 15)$value
        optimized_parameters_aai_15_inf_quadratic <- optim(par = c(k, j, i), fn=likelihood_general, longevity = longevity_aai_15_inf, age_at_infection = 15)$par
        message("correction achieved for aai_15_inf population and i=",i," and j=", j, " and k=", k)
      }

      if(optim(par = c(k, j, i), fn=likelihood_general, longevity = longevity_aai_30_inf, age_at_infection = 30)$value < optimized_likelihood_aai_30_inf_quadratic){
        optimized_likelihood_aai_30_inf_quadratic <- optim(par = c(k, j, i), fn=likelihood_general, longevity = longevity_aai_30_inf, age_at_infection = 30)$value
        optimized_parameters_aai_30_inf_quadratic <- optim(par = c(k, j, i), fn=likelihood_general, longevity = longevity_aai_30_inf, age_at_infection = 30)$par
        message("correction achieved for aai_30_inf population and i=",i," and j=", j, " and k=", k)
      }

      if(optim(par = c(k, j, i), fn=likelihood_general, longevity = longevity_aai_5_exp, age_at_infection = 5)$value < optimized_likelihood_aai_5_exp_quadratic){
        optimized_likelihood_aai_5_exp_quadratic <- optim(par = c(k, j, i), fn=likelihood_general, longevity = longevity_aai_5_exp, age_at_infection = 5)$value
        optimized_parameters_aai_5_exp_quadratic <- optim(par = c(k, j, i), fn=likelihood_general, longevity = longevity_aai_5_exp, age_at_infection = 5)$par
        message("correction achieved for aai_5_exp population and i=",i," and j=", j, " and k=", k)
      }

      if(optim(par = c(k, j, i), fn=likelihood_general, longevity = longevity_aai_15_exp, age_at_infection = 15)$value < optimized_likelihood_aai_15_exp_quadratic){
        optimized_likelihood_aai_15_exp_quadratic <- optim(par = c(k, j, i), fn=likelihood_general, longevity = longevity_aai_15_exp, age_at_infection = 15)$value
        optimized_parameters_aai_15_exp_quadratic <- optim(par = c(k, j, i), fn=likelihood_general, longevity = longevity_aai_15_exp, age_at_infection = 15)$par
        message("correction achieved for aai_15_exp population and i=",i," and j=", j, " and k=", k)
      }

      if(optim(par = c(k, j, i), fn=likelihood_general, longevity = longevity_aai_30_exp, age_at_infection = 30)$value < optimized_likelihood_aai_30_exp_quadratic){
        optimized_likelihood_aai_30_exp_quadratic <- optim(par = c(k, j, i), fn=likelihood_general, longevity = longevity_aai_30_exp, age_at_infection = 30)$value
        optimized_parameters_aai_30_exp_quadratic <- optim(par = c(k, j, i), fn=likelihood_general, longevity = longevity_aai_30_exp, age_at_infection = 30)$par
        message("correction achieved for aai_30_exp population and i=",i," and j=", j, " and k=", k)
      }
      ###control
      if(optim(par = c(k,j, i), fn=likelihood_general, longevity = host_longevity_control, age_at_infection = 5)$value < optimized_likelihood_aai_5_control_quadratic){
        optimized_likelihood_aai_5_control_quadratic <- optim(par = c(k,j, i), fn=likelihood_general, longevity = host_longevity_control, age_at_infection = 5)$value
        optimized_parameters_aai_5_control_quadratic <- optim(par = c(k,j, i), fn=likelihood_general, longevity = host_longevity_control, age_at_infection = 5)$par
        message("correction achieved for aai_5_control population and i=",i," and j=", j, " and k=",k)
      }

      if(optim(par = c(k,j, i), fn=likelihood_general, longevity = host_longevity_control, age_at_infection = 15)$value < optimized_likelihood_aai_15_control_quadratic){
        optimized_likelihood_aai_15_control_quadratic <- optim(par = c(k,j, i), fn=likelihood_general, longevity = host_longevity_control, age_at_infection = 15)$value
        optimized_parameters_aai_15_control_quadratic <- optim(par = c(k,j, i), fn=likelihood_general, longevity = host_longevity_control, age_at_infection = 15)$par
        message("correction achieved for aai_15_control population and i=",i," and j=", j, " and k=",k)
      }

      if(optim(par = c(k,j, i), fn=likelihood_general, longevity = host_longevity_control, age_at_infection = 30)$value < optimized_likelihood_aai_30_control_quadratic){
        optimized_likelihood_aai_30_control_quadratic <- optim(par = c(k,j, i), fn=likelihood_general, longevity = host_longevity_control, age_at_infection = 30)$value
        optimized_parameters_aai_30_control_quadratic <- optim(par = c(k,j, i), fn=likelihood_general, longevity = host_longevity_control, age_at_infection = 30)$par
        message("correction achieved for aai_30_control population and i=",i," and j=", j, " and k=",k)
      }
      
      if(optim(par = c(k,j, i), fn=likelihood_general, longevity = longevity_uninfected, age_at_infection = 16)$value < optimized_likelihood_uninfected_quadratic){
        optimized_likelihood_uninfected_quadratic <- optim(par = c(k,j, i), fn=likelihood_general, longevity = longevity_uninfected, age_at_infection = 16)$value
        optimized_parameters_uninfected_quadratic <- optim(par = c(k,j, i), fn=likelihood_general, longevity = longevity_uninfected, age_at_infection = 16)$par
        message("correction achieved for uninfected population and i=",i," and j=", j, "and k=",k)
      }
      if(optim(par = c(k,j, i), fn=likelihood_general, longevity = longevity_exposed, age_at_infection = 16)$value < optimized_likelihood_exposed_quadratic){
        optimized_likelihood_exposed_quadratic <- optim(par = c(k,j, i), fn=likelihood_general, longevity = longevity_exposed, age_at_infection = 16)$value
        optimized_parameters_exposed_quadratic <- optim(par = c(k,j, i), fn=likelihood_general, longevity = longevity_exposed, age_at_infection = 16)$par
        message("correction achieved for uninfected population and i=",i," and j=", j, "and k=",k)
      }
      for (l in cubic_initial) {
        if(optim(par = c(l, k, j, i), fn=likelihood_general, longevity = longevity_aai_5_inf, age_at_infection = 5)$value < optimized_likelihood_aai_5_inf_cubic){
          optimized_likelihood_aai_5_inf_cubic <- optim(par = c(l, k, j, i), fn=likelihood_general, longevity = longevity_aai_5_inf, age_at_infection = 5)$value
          optimized_parameters_aai_5_inf_cubic <- optim(par = c(l, k, j, i), fn=likelihood_general, longevity = longevity_aai_5_inf, age_at_infection = 5)$par
          message("correction achieved for aai_5_inf population and i=",i," and j=", j, " and k=", k, "and l=", l)
        }

        if(optim(par = c(l, k, j, i), fn=likelihood_general, longevity = longevity_aai_15_inf, age_at_infection = 15)$value < optimized_likelihood_aai_15_inf_cubic){
          optimized_likelihood_aai_15_inf_cubic <- optim(par = c(l, k, j, i), fn=likelihood_general, longevity = longevity_aai_15_inf, age_at_infection = 15)$value
          optimized_parameters_aai_15_inf_cubic <- optim(par = c(l, k, j, i), fn=likelihood_general, longevity = longevity_aai_15_inf, age_at_infection = 15)$par
          message("correction achieved for aai_15_inf population and i=",i," and j=", j, " and k=", k, "and l=", l)
        }

        if(optim(par = c(l, k, j, i), fn=likelihood_general, longevity = longevity_aai_30_inf, age_at_infection = 30)$value < optimized_likelihood_aai_30_inf_cubic){
          optimized_likelihood_aai_30_inf_cubic <- optim(par = c(l, k, j, i), fn=likelihood_general, longevity = longevity_aai_30_inf, age_at_infection = 30)$value
          optimized_parameters_aai_30_inf_cubic <- optim(par = c(l, k, j, i), fn=likelihood_general, longevity = longevity_aai_30_inf, age_at_infection = 30)$par
          message("correction achieved for aai_30_inf population and i=",i," and j=", j, " and k=", k, "and l=", l)
        }

        if(optim(par = c(l, k, j, i), fn=likelihood_general, longevity = longevity_aai_5_exp, age_at_infection = 5)$value < optimized_likelihood_aai_5_exp_cubic){
          optimized_likelihood_aai_5_exp_cubic <- optim(par = c(l, k, j, i), fn=likelihood_general, longevity = longevity_aai_5_exp, age_at_infection = 5)$value
          optimized_parameters_aai_5_exp_cubic <- optim(par = c(l, k, j, i), fn=likelihood_general, longevity = longevity_aai_5_exp, age_at_infection = 5)$par
          message("correction achieved for aai_5_exp population and i=",i," and j=", j, " and k=", k, "and l=", l)
        }

        if(optim(par = c(l, k, j, i), fn=likelihood_general, longevity = longevity_aai_15_exp, age_at_infection = 15)$value < optimized_likelihood_aai_15_exp_cubic){
          optimized_likelihood_aai_15_exp_cubic <- optim(par = c(l, k, j, i), fn=likelihood_general, longevity = longevity_aai_15_exp, age_at_infection = 15)$value
          optimized_parameters_aai_15_exp_cubic <- optim(par = c(l, k, j, i), fn=likelihood_general, longevity = longevity_aai_15_exp, age_at_infection = 15)$par
          message("correction achieved for aai_15_exp population and i=",i," and j=", j, " and k=", k, "and l=", l)
        }

        if(optim(par = c(l, k, j, i), fn=likelihood_general, longevity = longevity_aai_30_exp, age_at_infection = 30)$value < optimized_likelihood_aai_30_exp_cubic){
          optimized_likelihood_aai_30_exp_cubic <- optim(par = c(l, k, j, i), fn=likelihood_general, longevity = longevity_aai_30_exp, age_at_infection = 30)$value
          optimized_parameters_aai_30_exp_cubic <- optim(par = c(l, k, j, i), fn=likelihood_general, longevity = longevity_aai_30_exp, age_at_infection = 30)$par
          message("correction achieved for aai_30_exp population and i=",i," and j=", j, " and k=", k, "and l=", l)
        }

        ###control
        if(optim(par = c(l,k,j, i), fn=likelihood_general, longevity = host_longevity_control, age_at_infection = 5)$value < optimized_likelihood_aai_5_control_cubic){
          optimized_likelihood_aai_5_control_cubic <- optim(par = c(l,k,j, i), fn=likelihood_general, longevity = host_longevity_control, age_at_infection = 5)$value
          optimized_parameters_aai_5_control_cubic <- optim(par = c(l,k,j, i), fn=likelihood_general, longevity = host_longevity_control, age_at_infection = 5)$par
          message("correction achieved for aai_5_control population and i=",i," and j=", j, " and k=",k," and l=",l)
        }

        if(optim(par = c(l,k,j, i), fn=likelihood_general, longevity = host_longevity_control, age_at_infection = 15)$value < optimized_likelihood_aai_15_control_cubic){
          optimized_likelihood_aai_15_control_cubic <- optim(par = c(l,k,j, i), fn=likelihood_general, longevity = host_longevity_control, age_at_infection = 15)$value
          optimized_parameters_aai_15_control_cubic <- optim(par = c(l,k,j, i), fn=likelihood_general, longevity = host_longevity_control, age_at_infection = 15)$par
          message("correction achieved for aai_15_control population and i=",i," and j=", j, " and k=",k," and l=",l)
        }

        if(optim(par = c(l,k,j, i), fn=likelihood_general, longevity = host_longevity_control, age_at_infection = 30)$value < optimized_likelihood_aai_30_control_cubic){
          optimized_likelihood_aai_30_control_cubic <- optim(par = c(l,k,j, i), fn=likelihood_general, longevity = host_longevity_control, age_at_infection = 30)$value
          optimized_parameters_aai_30_control_cubic <- optim(par = c(l,k,j, i), fn=likelihood_general, longevity = host_longevity_control, age_at_infection = 30)$par
          message("correction achieved for aai_30_control population and i=",i," and j=", j, " and k=",k," and l=",l)
        }
        
        if(optim(par = c(l,k,j, i), fn=likelihood_general, longevity = longevity_uninfected, age_at_infection = 16)$value < optimized_likelihood_uninfected_cubic){
          optimized_likelihood_uninfected_cubic <- optim(par = c(l, k, j, i), fn=likelihood_general, longevity = longevity_uninfected, age_at_infection = 16)$value
          optimized_parameters_uninfected_cubic <- optim(par = c(l, k, j, i), fn=likelihood_general, longevity = longevity_uninfected, age_at_infection = 16)$par
          message("correction achieved for uninfected population and i=",i," and j=", j, "and k=",k," and l=",l)
        }
        if(optim(par = c(l,k,j, i), fn=likelihood_general, longevity = longevity_exposed, age_at_infection = 16)$value < optimized_likelihood_exposed_cubic){
          optimized_likelihood_exposed_cubic <- optim(par = c(l, k, j, i), fn=likelihood_general, longevity = longevity_exposed, age_at_infection = 16)$value
          optimized_parameters_exposed_cubic <- optim(par = c(l, k, j, i), fn=likelihood_general, longevity = longevity_exposed, age_at_infection = 16)$par
          message("correction achieved for uninfected population and i=",i," and j=", j, "and k=",k," and l=",l)
        }
        
        for (m in quartic_initial) {
          if(optim(par = c(m, l, k, j, i), fn=likelihood_general, longevity = longevity_aai_5_inf, age_at_infection = 5)$value < optimized_likelihood_aai_5_inf_quartic){
            optimized_likelihood_aai_5_inf_quartic <- optim(par = c(m, l, k, j, i), fn=likelihood_general, longevity = longevity_aai_5_inf, age_at_infection = 5)$value
            optimized_parameters_aai_5_inf_quartic <- optim(par = c(m, l, k, j, i), fn=likelihood_general, longevity = longevity_aai_5_inf, age_at_infection = 5)$par
            message("correction achieved for aai_5_inf population and i=",i," and j=", j, " and k=", k, "and l=", l, " and m=", m)
          }

          if(optim(par = c(m, l, k, j, i), fn=likelihood_general, longevity = longevity_aai_15_inf, age_at_infection = 15)$value < optimized_likelihood_aai_15_inf_quartic){
            optimized_likelihood_aai_15_inf_quartic <- optim(par = c(m, l, k, j, i), fn=likelihood_general, longevity = longevity_aai_15_inf, age_at_infection = 15)$value
            optimized_parameters_aai_15_inf_quartic <- optim(par = c(m, l, k, j, i), fn=likelihood_general, longevity = longevity_aai_15_inf, age_at_infection = 15)$par
            message("correction achieved for aai_15_inf population and i=",i," and j=", j, " and k=", k, "and l=", l, " and m=", m)
          }

          if(optim(par = c(m, l, k, j, i), fn=likelihood_general, longevity = longevity_aai_30_inf, age_at_infection = 30)$value < optimized_likelihood_aai_30_inf_quartic){
            optimized_likelihood_aai_30_inf_quartic <- optim(par = c(m, l, k, j, i), fn=likelihood_general, longevity = longevity_aai_30_inf, age_at_infection = 30)$value
            optimized_parameters_aai_30_inf_quartic <- optim(par = c(m, l, k, j, i), fn=likelihood_general, longevity = longevity_aai_30_inf, age_at_infection = 30)$par
            message("correction achieved for aai_30_inf population and i=",i," and j=", j, " and k=", k, "and l=", l, " and m=", m)
          }

          if(optim(par = c(m, l, k, j, i), fn=likelihood_general, longevity = longevity_aai_5_exp, age_at_infection = 5)$value < optimized_likelihood_aai_5_exp_quartic){
            optimized_likelihood_aai_5_exp_quartic <- optim(par = c(m, l, k, j, i), fn=likelihood_general, longevity = longevity_aai_5_exp, age_at_infection = 5)$value
            optimized_parameters_aai_5_exp_quartic <- optim(par = c(m, l, k, j, i), fn=likelihood_general, longevity = longevity_aai_5_exp, age_at_infection = 5)$par
            message("correction achieved for aai_5_exp population and i=",i," and j=", j, " and k=", k, "and l=", l, " and m=", m)
          }

          if(optim(par = c(m, l, k, j, i), fn=likelihood_general, longevity = longevity_aai_15_exp, age_at_infection = 15)$value < optimized_likelihood_aai_15_exp_quartic){
            optimized_likelihood_aai_15_exp_quartic <- optim(par = c(m, l, k, j, i), fn=likelihood_general, longevity = longevity_aai_15_exp, age_at_infection = 15)$value
            optimized_parameters_aai_15_exp_quartic <- optim(par = c(m, l, k, j, i), fn=likelihood_general, longevity = longevity_aai_15_exp, age_at_infection = 15)$par
            message("correction achieved for aai_15_exp population and i=",i," and j=", j, " and k=", k, "and l=", l, " and m=", m)
          }

          if(optim(par = c(m, l, k, j, i), fn=likelihood_general, longevity = longevity_aai_30_exp, age_at_infection = 30)$value < optimized_likelihood_aai_30_exp_quartic){
            optimized_likelihood_aai_30_exp_quartic <- optim(par = c(m, l, k, j, i), fn=likelihood_general, longevity = longevity_aai_30_exp, age_at_infection = 30)$value
            optimized_parameters_aai_30_exp_quartic<- optim(par = c(m, l, k, j, i), fn=likelihood_general, longevity = longevity_aai_30_exp, age_at_infection = 30)$par
            message("correction achieved for aai_30_exp population and i=",i," and j=", j, " and k=", k, "and l=", l, " and m=", m)
          }

          ###control
          if(optim(par = c(m,l,k,j, i), fn=likelihood_general, longevity = host_longevity_control, age_at_infection = 5)$value < optimized_likelihood_aai_5_control_quartic){
            optimized_likelihood_aai_5_control_quartic <- optim(par = c(m,l,k,j, i), fn=likelihood_general, longevity = host_longevity_control, age_at_infection = 5)$value
            optimized_parameters_aai_5_control_quartic <- optim(par = c(m,l,k,j, i), fn=likelihood_general, longevity = host_longevity_control, age_at_infection = 5)$par
            message("correction achieved for aai_5_control population and i=",i," and j=", j, " and k=",k," and l=",l," and m=",m)
          }

          if(optim(par = c(m,l,k,j, i), fn=likelihood_general, longevity = host_longevity_control, age_at_infection = 15)$value < optimized_likelihood_aai_15_control_quartic){
            optimized_likelihood_aai_15_control_quartic <- optim(par = c(m,l,k,j, i), fn=likelihood_general, longevity = host_longevity_control, age_at_infection = 15)$value
            optimized_parameters_aai_15_control_quartic <- optim(par = c(m,l,k,j, i), fn=likelihood_general, longevity = host_longevity_control, age_at_infection = 15)$par
            message("correction achieved for aai_15_control population and i=",i," and j=", j, " and k=",k," and l=",l," and m=",m)
          }

          if(optim(par = c(m,l,k,j, i), fn=likelihood_general, longevity = host_longevity_control, age_at_infection = 30)$value < optimized_likelihood_aai_30_control_quartic){
            optimized_likelihood_aai_30_control_quartic <- optim(par = c(m,l,k,j, i), fn=likelihood_general, longevity = host_longevity_control, age_at_infection = 30)$value
            optimized_parameters_aai_30_control_quartic <- optim(par = c(m,l,k,j, i), fn=likelihood_general, longevity = host_longevity_control, age_at_infection = 30)$par
            message("correction achieved for aai_30_control population and i=",i," and j=", j, " and k=",k," and l=",l," and m=",m)
          }
          
          if(optim(par = c(m,l,k,j, i), fn=likelihood_general, longevity = longevity_uninfected, age_at_infection = 16)$value < optimized_likelihood_uninfected_quartic){
            optimized_likelihood_uninfected_quartic <- optim(par = c(m,l, k, j, i), fn=likelihood_general, longevity = longevity_uninfected, age_at_infection = 16)$value
            optimized_parameters_uninfected_quartic <- optim(par = c(m,l, k, j, i), fn=likelihood_general, longevity = longevity_uninfected, age_at_infection = 16)$par
            message("correction achieved for uninfected population and i=",i," and j=", j, "and k=",k," and l=",l, " and m=",m)
          }
          if(optim(par = c(m,l,k,j, i), fn=likelihood_general, longevity = longevity_exposed, age_at_infection = 16)$value < optimized_likelihood_exposed_quartic){
            optimized_likelihood_exposed_quartic <- optim(par = c(m,l, k, j, i), fn=likelihood_general, longevity = longevity_exposed, age_at_infection = 16)$value
            optimized_parameters_exposed_quartic <- optim(par = c(m,l, k, j, i), fn=likelihood_general, longevity = longevity_exposed, age_at_infection = 16)$par
            message("correction achieved for uninfected population and i=",i," and j=", j, "and k=",k," and l=",l, " and m=",m)
          }
        }
      }
    }
  }
}
### deleted the rest below the optimization cycles, if needed check the original script for reference, only survival_model() kept
#### Survival model -------------------------------------------------------------
survival_model <- function(delta, x = age_axis){
  exponent <- c(rep(0,length(x)))
  for (a in x) {
    exponent[a+1] <- -sum(delta[1:(a+1)])
  }
  return(exp(exponent))
}
#messages ----------------------------------------------------------------------------------------------
message("parameters for aai_5_inf linear: ", optimized_parameters_aai_5_inf_linear)
message("parameters for aai_5_inf quadratic: ", optimized_parameters_aai_5_inf_quadratic)
message("parameters for aai_5_inf cubic: ", optimized_parameters_aai_5_inf_cubic)
message("parameters for aai_5_inf quartic: ", optimized_parameters_aai_5_inf_quartic)

message("parameters for aai_15_inf linear: ", optimized_parameters_aai_15_inf_linear)
message("parameters for aai_15_inf quadratic: ", optimized_parameters_aai_15_inf_quadratic)
message("parameters for aai_15_inf cubic: ", optimized_parameters_aai_15_inf_cubic)
message("parameters for aai_15_inf quartic: ", optimized_parameters_aai_15_inf_quartic)

message("parameters for aai_30_inf linear: ", optimized_parameters_aai_30_inf_linear)
message("parameters for aai_30_inf quadratic: ", optimized_parameters_aai_30_inf_quadratic)
message("parameters for aai_30_inf cubic: ", optimized_parameters_aai_30_inf_cubic)
message("parameters for aai_30_inf quartic: ", optimized_parameters_aai_30_inf_quartic)

message("parameters for aai_5_exp linear: ", optimized_parameters_aai_5_exp_linear)
message("parameters for aai_5_exp quadratic: ", optimized_parameters_aai_5_exp_quadratic)
message("parameters for aai_5_exp cubic: ", optimized_parameters_aai_5_exp_cubic)
message("parameters for aai_5_exp quartic: ", optimized_parameters_aai_5_exp_quartic)

message("parameters for aai_15_exp linear: ", optimized_parameters_aai_15_exp_linear)
message("parameters for aai_15_exp quadratic: ", optimized_parameters_aai_15_exp_quadratic)
message("parameters for aai_15_exp cubic: ", optimized_parameters_aai_15_exp_cubic)
message("parameters for aai_15_exp quartic: ", optimized_parameters_aai_15_exp_quartic)

message("parameters for aai_30_exp linear: ", optimized_parameters_aai_30_exp_linear)
message("parameters for aai_30_exp quadratic: ", optimized_parameters_aai_30_exp_quadratic)
message("parameters for aai_30_exp cubic: ", optimized_parameters_aai_30_exp_cubic)
message("parameters for aai_30_exp quartic: ", optimized_parameters_aai_30_exp_quartic)

message("parameters for aai_5_control linear: ", optimized_parameters_aai_5_control_linear)
message("parameters for aai_5_control quadratic: ", optimized_parameters_aai_5_control_quadratic)
message("parameters for aai_5_control cubic: ", optimized_parameters_aai_5_control_cubic)
message("parameters for aai_5_control quartic: ", optimized_parameters_aai_5_control_quartic)

message("parameters for aai_15_control linear: ", optimized_parameters_aai_15_control_linear)
message("parameters for aai_15_control quadratic: ", optimized_parameters_aai_15_control_quadratic)
message("parameters for aai_15_control cubic: ", optimized_parameters_aai_15_control_cubic)
message("parameters for aai_15_control quartic: ", optimized_parameters_aai_15_control_quartic)

message("parameters for aai_30_control linear: ", optimized_parameters_aai_30_control_linear)
message("parameters for aai_30_control quadratic: ", optimized_parameters_aai_30_control_quadratic)
message("parameters for aai_30_control cubic: ", optimized_parameters_aai_30_control_cubic)
message("parameters for aai_30_control quartic: ", optimized_parameters_aai_30_control_quartic)

#### LRTs --------------------------------------------------------------------------------------------------------------------------------
likelihoods_exp_linear <- c(optimized_likelihood_aai_5_exp_linear, optimized_likelihood_aai_15_exp_linear,optimized_likelihood_aai_30_exp_linear)
likelihoods_exp_quadratic <- c(optimized_likelihood_aai_5_exp_quadratic, optimized_likelihood_aai_15_exp_quadratic, optimized_likelihood_aai_30_exp_quadratic)
likelihoods_exp_cubic <- c(optimized_likelihood_aai_5_exp_cubic, optimized_likelihood_aai_15_exp_cubic, optimized_likelihood_aai_30_exp_cubic)
likelihoods_exp_quartic <- c(optimized_likelihood_aai_5_exp_quartic, optimized_likelihood_aai_15_exp_quartic, optimized_likelihood_aai_30_exp_quartic)

likelihoods_inf_linear <- c(optimized_likelihood_aai_5_inf_linear, optimized_likelihood_aai_15_inf_linear,optimized_likelihood_aai_30_inf_linear)
likelihoods_inf_quadratic <- c(optimized_likelihood_aai_5_inf_quadratic, optimized_likelihood_aai_15_inf_quadratic, optimized_likelihood_aai_30_inf_quadratic)
likelihoods_inf_cubic <- c(optimized_likelihood_aai_5_inf_cubic, optimized_likelihood_aai_15_inf_cubic, optimized_likelihood_aai_30_inf_cubic)
likelihoods_inf_quartic <- c(optimized_likelihood_aai_5_inf_quartic, optimized_likelihood_aai_15_inf_quartic, optimized_likelihood_aai_30_inf_quartic)

likelihoods_control_linear <- c(optimized_likelihood_aai_5_control_linear, optimized_likelihood_aai_15_control_linear,optimized_likelihood_aai_30_control_linear)
likelihoods_control_quadratic <- c(optimized_likelihood_aai_5_control_quadratic, optimized_likelihood_aai_15_control_quadratic, optimized_likelihood_aai_30_control_quadratic)
likelihoods_control_cubic <- c(optimized_likelihood_aai_5_control_cubic, optimized_likelihood_aai_15_control_cubic, optimized_likelihood_aai_30_control_cubic)
likelihoods_control_quartic <- c(optimized_likelihood_aai_5_control_quartic, optimized_likelihood_aai_15_control_quartic, optimized_likelihood_aai_30_control_quartic)

LRT_lin_quadratic_exp <- 2*(likelihoods_exp_linear - likelihoods_exp_quadratic)
LRT_quadratic_cub_exp <- 2*(likelihoods_exp_quadratic - likelihoods_exp_cubic)
LRT_cub_quartic_exp <- 2*(likelihoods_exp_cubic - likelihoods_exp_quartic)

LRT_lin_quadratic_inf <- 2*(likelihoods_inf_linear - likelihoods_inf_quadratic)
LRT_quadratic_cub_inf <- 2*(likelihoods_inf_quadratic - likelihoods_inf_cubic)
LRT_cub_quartic_inf <- 2*(likelihoods_inf_cubic - likelihoods_inf_quartic)

LRT_lin_quadratic_control<- 2*(likelihoods_control_linear - likelihoods_control_quadratic)
LRT_quadratic_cub_control <- 2*(likelihoods_control_quadratic - likelihoods_control_cubic)
LRT_cub_quartic_control <- 2*(likelihoods_control_cubic - likelihoods_control_quartic)

pchisq(LRT_lin_quadratic_exp, df = 1, lower.tail = FALSE)
pchisq(LRT_quadratic_cub_exp, df = 1, lower.tail = FALSE)
pchisq(LRT_cub_quartic_exp, df = 1, lower.tail = FALSE) # for age at infection 5, best fit is the cubic one, otherwise quartic

pchisq(LRT_lin_quadratic_inf, df = 1, lower.tail = FALSE)
pchisq(LRT_quadratic_cub_inf, df = 1, lower.tail = FALSE)
pchisq(LRT_cub_quartic_inf, df = 1, lower.tail = FALSE) # for infecteds quadratic is always the best fit

pchisq(LRT_lin_quadratic_control, df = 1, lower.tail = FALSE)
pchisq(LRT_quadratic_cub_control, df = 1, lower.tail = FALSE)
pchisq(LRT_cub_quartic_control, df = 1, lower.tail = FALSE) # for age at infection 30, best fit is the cubic one, otherwise quartic


#### Kolmogorov-Smirnov tests for survival distributions -----------------------------------
## Does Survival of control and exposed populations come from the same distribution?
ks.test(survivals_control, survivals_exposed)
ks.test(survivals_aai_5_exp[5:lifetime_threshold], survivals_aai_5_control[5:lifetime_threshold])
ks.test(survivals_aai_15_exp[15:lifetime_threshold], survivals_aai_15_control[15:lifetime_threshold]) ### Same distribution
ks.test(survivals_aai_30_exp[30:lifetime_threshold], survivals_aai_30_control[30:lifetime_threshold]) ### Same distribution

# The above Kolmogorov-Smirnov test is not completely right, since it is applied to the cummulative functions, which
# might not be completely independent

### plotting frequencies of lifetimes ---------------------------------------------
frequencies_deaths_exp_aai_5 <- rep(0, (lifetime_threshold-4))
frequencies_deaths_exp_aai_15 <- rep(0, (lifetime_threshold-14))
frequencies_deaths_exp_aai_30 <- rep(0, (lifetime_threshold-29))
frequencies_deaths_control_aai_5 <- rep(0, (lifetime_threshold-4))
frequencies_deaths_control_aai_15 <- rep(0, (lifetime_threshold-14))
frequencies_deaths_control_aai_30 <- rep(0, (lifetime_threshold-29))

for(i in 5:lifetime_threshold){
  frequencies_deaths_exp_aai_5[i-4] <- sum(longevity_aai_5_exp == i)
  frequencies_deaths_control_aai_5[i-4] <- sum(host_longevity_control == i)
}
for(i in 15:lifetime_threshold){
  frequencies_deaths_exp_aai_15[i-14] <- sum(longevity_aai_15_exp == i)
  frequencies_deaths_control_aai_15[i-14] <- sum(host_longevity_control == i)
}

for(i in 30:lifetime_threshold){
  frequencies_deaths_exp_aai_30[i-29] <- sum(longevity_aai_30_exp == i)
  frequencies_deaths_control_aai_30[i-29] <- sum(host_longevity_control == i)
}

frequencies_deaths_exp_aai_5 <- frequencies_deaths_exp_aai_5/sum(frequencies_deaths_exp_aai_5)
frequencies_deaths_exp_aai_15 <- frequencies_deaths_exp_aai_15/sum(frequencies_deaths_exp_aai_15)
frequencies_deaths_exp_aai_30 <- frequencies_deaths_exp_aai_30/sum(frequencies_deaths_exp_aai_30)
frequencies_deaths_control_aai_5 <- frequencies_deaths_control_aai_5/sum(frequencies_deaths_control_aai_5)
frequencies_deaths_control_aai_15 <- frequencies_deaths_control_aai_15/sum(frequencies_deaths_control_aai_15)
frequencies_deaths_control_aai_30 <- frequencies_deaths_control_aai_30/sum(frequencies_deaths_control_aai_30)

plot(age_of_inf_5_axis+5, frequencies_deaths_exp_aai_5, col = "blue")
lines(age_of_inf_5_axis+5, frequencies_deaths_control_aai_5, col = "red", type = "p")


plot(age_of_inf_15_axis+15, frequencies_deaths_exp_aai_15, col = "blue")
lines(age_of_inf_15_axis+15, frequencies_deaths_control_aai_15, col = "red", type = "p")

plot(age_of_inf_30_axis+30, frequencies_deaths_exp_aai_30, col = "blue")
lines(age_of_inf_30_axis+30, frequencies_deaths_control_aai_30, col = "red", type = "p")

### plotting histograms of lifetimes --------------------------------------------------------------

hist(longevity_aai_5_exp, freq = FALSE, col = rgb(1,0,0,0.7))
hist(host_longevity_control[host_longevity_control>5], freq = FALSE, add=T, col = rgb(0,1,0,0.6))

hist(longevity_aai_15_exp, freq = FALSE, col = rgb(1,0,0,0.7))
hist(host_longevity_control[host_longevity_control>15], freq = FALSE, add=T, col = rgb(0,1,0,0.6))

hist(longevity_aai_30_exp, freq = FALSE, col = rgb(1,0,0,0.7))
hist(host_longevity_control[host_longevity_control>30], freq = FALSE, add=T, col = rgb(0,1,0,0.6))


##### Plotting ------------------------------------------------------------------------------------
### Comparing exposed and control populations------------------------------------------------------
## Age at infection = 5 ---------------------------------------------------------------------------
plot(x = age_axis[5:lifetime_threshold], y = survivals_aai_5_exp[5:lifetime_threshold], xlab = "Age [Days]", ylab = "Survival_Exposed_At_Age=5", type = "s", lwd = 3)
lines(x = age_axis[5:lifetime_threshold], y = survivals_aai_5_control[5:lifetime_threshold], lwd = 3, type = "s", col = "red")
# Quartic model Control
lines(x = age_axis[5:lifetime_threshold], y = survival_model(delta = delta(age = age_axis[5:lifetime_threshold], parameters = optimized_parameters_aai_5_control_quartic), x = age_of_inf_5_axis), col = "red")
# Quartic model Exposed
lines(x = age_axis[5:lifetime_threshold], y = survival_model(delta = delta(age = age_axis[5:lifetime_threshold], parameters = optimized_parameters_aai_5_exp_quartic), x = age_of_inf_5_axis), col = "black")
legend("bottomleft", legend = c("Data exposed", "Data control", "Quartic model control", "Quartic model exposed"), col = c("black", "red", "black", "red"), lwd = c(3,3,1,1), cex = 0.8)

pdf("Survival_data_aai_5.pdf")
plot(x = age_axis[5:lifetime_threshold], y = survivals_aai_5_exp[5:lifetime_threshold], xlab = "Age [Days]", ylab = "Survival [fraction]", type = "s", lwd = 3)
lines(x = age_axis[5:lifetime_threshold], y = survivals_aai_5_control[5:lifetime_threshold], lwd = 3, col = "red", type = "s")
legend("bottomleft", legend = c("Data exposed (at age 5)", "Data control (since age 5)"), col = c("black", "red"), lwd = c(3,3), cex = 0.8)
dev.off()

plot(x = age_axis[5:lifetime_threshold], y = survival_model(delta = delta(age = age_axis[5:lifetime_threshold], parameters = optimized_parameters_aai_5_control_quartic), x = age_of_inf_5_axis), col = "red", xlab = "Age [Days]", ylab = "Survival_Exposed_At_Age=5", type = "l", lwd = 3)
# Quartic model Exposed
lines(x = age_axis[5:lifetime_threshold], y = survival_model(delta = delta(age = age_axis[5:lifetime_threshold], parameters = optimized_parameters_aai_5_exp_quartic), x = age_of_inf_5_axis), col = "black")
legend("bottomleft", legend = c("Quartic model exposed", "Quartic model control"), col = c("black", "red"), lwd = c(1,1), cex = 0.8)

## Age at infection = 15 -----------------------------------------------------------------------------
plot(x = age_axis[15:lifetime_threshold], y = survivals_aai_15_exp[15:lifetime_threshold], xlab = "Age [Days]", ylab = "Survival_Exposed_At_Age=15", type = "s", lwd = 3)
lines(x = age_axis[15:lifetime_threshold], y = survivals_aai_15_control[15:lifetime_threshold], lwd = 3, type = "s", col = "red")
# Quartic model Control
lines(x = age_axis[15:lifetime_threshold], y = survival_model(delta = delta(age = age_axis[15:lifetime_threshold], parameters = optimized_parameters_aai_15_control_quartic), x = age_of_inf_15_axis), col = "red")
# Quartic model Exposed
lines(x = age_axis[15:lifetime_threshold], y = survival_model(delta = delta(age = age_axis[15:lifetime_threshold], parameters = optimized_parameters_aai_15_exp_quartic), x = age_of_inf_15_axis), col = "black")
legend("bottomleft", legend = c("Data exposed", "Data control", "Quartic model control", "Quartic model exposed"), col = c("black", "red", "black", "red"), lwd = c(3,3,1,1), cex = 0.8)
str(survivals_aai_15_control[15:lifetime_threshold])

pdf("Survival_data_aai_15.pdf")
plot(x = age_axis[15:lifetime_threshold], y = survivals_aai_15_exp[15:lifetime_threshold], xlab = "Age [Days]", ylab = "Survival[fraction]", type = "s", lwd = 3)
lines(x = age_axis[15:lifetime_threshold], y = survivals_aai_15_control[15:lifetime_threshold], lwd = 3, type = "s", col = "red")
legend("bottomleft", legend = c("Data exposed (at age 15)", "Data control (since age 15)"), col = c("black", "red"), lwd = c(3,3), cex = 0.8)
dev.off()

plot(x = age_axis[15:lifetime_threshold], y = survival_model(delta = delta(age = age_axis[15:lifetime_threshold], parameters = optimized_parameters_aai_15_control_quartic), x = age_of_inf_15_axis), col = "red", xlab = "Age [Days]", ylab = "Survival_Exposed_At_Age=15", type = "l", lwd = 3)
# Quartic model Exposed
lines(x = age_axis[15:lifetime_threshold], y = survival_model(delta = delta(age = age_axis[15:lifetime_threshold], parameters = optimized_parameters_aai_15_exp_quartic), x = age_of_inf_15_axis), col = "black")
legend("bottomleft", legend = c("Quartic model exposed", "Quartic model control"), col = c("black", "red"), lwd = c(1,1), cex = 0.8)

## Age at infection = 30 ------------------------------------------------------------------------------
plot(x = age_axis[30:lifetime_threshold], y = survivals_aai_30_exp[30:lifetime_threshold], xlab = "Age [Days]", ylab = "Survival_Exposed_At_Age=30", type = "s", lwd = 3)
lines(x = age_axis[30:lifetime_threshold], y = survivals_aai_30_control[30:lifetime_threshold], lwd = 3, type = "s", col = "red")
# Quartic model Control
lines(x = age_axis[30:lifetime_threshold], y = survival_model(delta = delta(age = age_axis[30:lifetime_threshold], parameters = optimized_parameters_aai_30_control_quartic), x = age_of_inf_30_axis), col = "red")
# Quartic model Exposed
lines(x = age_axis[30:lifetime_threshold], y = survival_model(delta = delta(age = age_axis[30:lifetime_threshold], parameters = optimized_parameters_aai_30_exp_quartic), x = age_of_inf_30_axis), col = "black")
legend("bottomleft", legend = c("Data exposed", "Data control", "Quartic model control", "Quartic model exposed"), col = c("black", "red", "black", "red"), lwd = c(3,3,1,1), cex = 0.8)

pdf("Survival_data_aai_30.pdf")
plot(x = age_axis[30:lifetime_threshold], y = survivals_aai_30_exp[30:lifetime_threshold], xlab = "Age [Days]", ylab = "Survival [fraction]", type = "s", lwd = 3)
lines(x = age_axis[30:lifetime_threshold], y = survivals_aai_30_control[30:lifetime_threshold], lwd = 3, type = "s", col = "red")
legend("bottomleft", legend = c("Data exposed (at age 30)", "Data control (since age 30)"), col = c("black", "red"), lwd = c(3,3), cex = 0.8)
dev.off()

plot(x = age_axis[30:lifetime_threshold], y = survival_model(delta = delta(age = age_axis[30:lifetime_threshold], parameters = optimized_parameters_aai_30_control_quartic), x = age_of_inf_30_axis), col = "red", xlab = "Age [Days]", ylab = "Survival_Exposed_At_Age=30", type = "l", lwd = 3)
# Quartic model Exposed
lines(x = age_axis[30:lifetime_threshold], y = survival_model(delta = delta(age = age_axis[30:lifetime_threshold], parameters = optimized_parameters_aai_30_exp_quartic), x = age_of_inf_30_axis), col = "black")
legend("bottomleft", legend = c("Quartic model exposed", "Quartic model control"), col = c("black", "red"), lwd = c(1,1), cex = 0.8)


### Control plots ---------------------------------------------------------------------------------------------------------------------------
## AaI = 5, inf --------------------------------------------------------------------------------------------------
pdf("Fitting_delta_inf_at_age_5.pdf")
plot(x = age_axis[5:lifetime_threshold], y = survivals_aai_5_inf[5:lifetime_threshold], xlab = "Age [Days]", ylab = "Survival_Infected_At_Age=5")
# Linear model
lines(x = age_axis[5:lifetime_threshold], y = survival_model(delta = delta(age = age_axis[5:lifetime_threshold], parameters = optimized_parameters_aai_5_inf_linear), x = age_of_inf_5_axis), col = "yellow")
# Quadratic model
lines(x = age_axis[5:lifetime_threshold], y = survival_model(delta = delta(age = age_axis[5:lifetime_threshold], parameters = optimized_parameters_aai_5_inf_quadratic), x = age_of_inf_5_axis), col = "green")
# Cubic model
lines(x = age_axis[5:lifetime_threshold], y = survival_model(delta = delta(age = age_axis[5:lifetime_threshold], parameters = optimized_parameters_aai_5_inf_cubic), x = age_of_inf_5_axis), col = "blue")
# Quartic model
lines(x = age_axis[5:lifetime_threshold], y = survival_model(delta = delta(age = age_axis[5:lifetime_threshold], parameters = optimized_parameters_aai_5_inf_quartic), x = age_of_inf_5_axis), col = "red")
legend(90,1, legend = c("Data", 'Linear model', "Quadratic model", "Cubic model", "4th order polynomial model"), col = c("black", "yellow", "green", "blue", "red"), lty = 1, cex = 0.8)
dev.off()

## AaI = 15, inf ----------------------------------------------------------------------------------------------------
pdf("Fitting_delta_inf_at_age_15.pdf")
plot(x = age_axis[15:lifetime_threshold], y = survivals_aai_15_inf[15:lifetime_threshold], xlab = "Age [Days]", ylab = "Survival_Infected_At_Age=15")
# Linear model
lines(x = age_axis[15:lifetime_threshold], y = survival_model(delta = delta(age = age_axis[15:lifetime_threshold], parameters = optimized_parameters_aai_15_inf_linear), x = age_of_inf_15_axis), col = "yellow")
# Quadratic model
lines(x = age_axis[15:lifetime_threshold], y = survival_model(delta = delta(age = age_axis[15:lifetime_threshold], parameters = optimized_parameters_aai_15_inf_quadratic), x = age_of_inf_15_axis), col = "green")
# Cubic model
lines(x = age_axis[15:lifetime_threshold], y = survival_model(delta = delta(age = age_axis[15:lifetime_threshold], parameters = optimized_parameters_aai_15_inf_cubic), x = age_of_inf_15_axis), col = "blue")
# Quartic model
lines(x = age_axis[15:lifetime_threshold], y = survival_model(delta = delta(age = age_axis[15:lifetime_threshold], parameters = optimized_parameters_aai_15_inf_quartic), x = age_of_inf_15_axis), col = "red")
legend(95,1, legend = c("Data", 'Linear model', "Quadratic model", "Cubic model", "4th order polynomial model"), col = c("black", "yellow", "green", "blue", "red"), lty = 1, cex = 0.8)
dev.off()


## AaI = 30, inf ---------------------------------------------------------------------------------------------------
pdf("Fitting_delta_inf_at_age_30.pdf")
plot(x = age_axis[30:lifetime_threshold], y = survivals_aai_30_inf[30:lifetime_threshold], xlab = "Age [Days]", ylab = "Survival_Infected_At_Age=30")
# Linear model
lines(x = age_axis[30:lifetime_threshold], y = survival_model(delta = delta(age = age_axis[30:lifetime_threshold], parameters = optimized_parameters_aai_30_inf_linear), x = age_of_inf_30_axis), col = "yellow")
# Quadratic model
lines(x = age_axis[30:lifetime_threshold], y = survival_model(delta = delta(age = age_axis[30:lifetime_threshold], parameters = optimized_parameters_aai_30_inf_quadratic), x = age_of_inf_30_axis), col = "green")
# Cubic model
lines(x = age_axis[30:lifetime_threshold], y = survival_model(delta = delta(age = age_axis[30:lifetime_threshold], parameters = optimized_parameters_aai_30_inf_cubic), x = age_of_inf_30_axis), col = "blue")
# Quartic model
lines(x = age_axis[30:lifetime_threshold], y = survival_model(delta = delta(age = age_axis[30:lifetime_threshold], parameters = optimized_parameters_aai_30_inf_quartic), x = age_of_inf_30_axis), col = "red")
legend(100,1, legend = c("Data", 'Linear model', "Quadratic model", "Cubic model", "4th order polynomial model"), col = c("black", "yellow", "green", "blue", "red"), lty = 1, cex = 0.8)
dev.off()

## AaI = 5, exp -------------------------------------------------------------------------------------

pdf("Fitting_delta_exp_at_age_5.pdf")
plot(x = age_axis[5:lifetime_threshold], y = survivals_aai_5_exp[5:lifetime_threshold], xlab = "Age [Days]", ylab = "Survival_Exposed_At_Age=5", type = "s", lwd=2)
# Linear model
lines(x = age_axis[5:lifetime_threshold], y = survival_model(delta = delta(age = age_axis[5:lifetime_threshold], parameters = optimized_parameters_aai_5_exp_linear), x = age_of_inf_5_axis), col = "yellow")
# Quadratic model
lines(x = age_axis[5:lifetime_threshold], y = survival_model(delta = delta(age = age_axis[5:lifetime_threshold], parameters = optimized_parameters_aai_5_exp_quadratic), x = age_of_inf_5_axis), col = "green")
# Cubic model
lines(x = age_axis[5:lifetime_threshold], y = survival_model(delta = delta(age = age_axis[5:lifetime_threshold], parameters = optimized_parameters_aai_5_exp_cubic), x = age_of_inf_5_axis), col = "blue")
# Quartic model
lines(x = age_axis[5:lifetime_threshold], y = survival_model(delta = delta(age = age_axis[5:lifetime_threshold], parameters = optimized_parameters_aai_5_exp_quartic), x = age_of_inf_5_axis), col = "red")
legend(90,1, legend = c("Data", 'Linear model', "Quadratic model", "Cubic model", "4th order polynomial model"), col = c("black", "yellow", "green", "blue", "red"), lty = 1, cex = 0.8)
dev.off()

## AaI = 15, exp -------------------------------------------------------------------------------------

pdf("Fitting_delta_exp_at_age_15.pdf")
plot(x = age_axis[15:lifetime_threshold], y = survivals_aai_15_exp[15:lifetime_threshold], xlab = "Age [Days]", ylab = "Survival_Exposed_At_Age=15", type = "s", lwd = 2)
# Linear model
lines(x = age_axis[15:lifetime_threshold], y = survival_model(delta = delta(age = age_axis[15:lifetime_threshold], parameters = optimized_parameters_aai_15_exp_linear), x = age_of_inf_15_axis), col = "yellow")
# Quadratic model
lines(x = age_axis[15:lifetime_threshold], y = survival_model(delta = delta(age = age_axis[15:lifetime_threshold], parameters = optimized_parameters_aai_15_exp_quadratic), x = age_of_inf_15_axis), col = "green")
# Cubic model
lines(x = age_axis[15:lifetime_threshold], y = survival_model(delta = delta(age = age_axis[15:lifetime_threshold], parameters = optimized_parameters_aai_15_inf_cubic), x = age_of_inf_15_axis), col = "blue")
# Quartic model
lines(x = age_axis[15:lifetime_threshold], y = survival_model(delta = delta(age = age_axis[15:lifetime_threshold], parameters = optimized_parameters_aai_15_exp_quartic), x = age_of_inf_15_axis), col = "red")
legend(95,1, legend = c("Data", 'Linear model', "Quadratic model", "Cubic model", "4th order polynomial model"), col = c("black", "yellow", "green", "blue", "red"), lty = 1, cex = 0.8)
dev.off()

## AaI = 30, exp---------------------------------------------------------------------------------------
pdf("Fitting_delta_exp_at_age_30.pdf")
plot(x = age_axis[30:lifetime_threshold], y = survivals_aai_30_exp[30:lifetime_threshold], xlab = "Age [Days]", ylab = "Survival_Exposed_At_Age=30", type = "s", lwd=2)
# Linear model
lines(x = age_axis[30:lifetime_threshold], y = survival_model(delta = delta(age = age_axis[30:lifetime_threshold], parameters = optimized_parameters_aai_30_exp_linear), x = age_of_inf_30_axis), col = "yellow")
# Quadratic model
lines(x = age_axis[30:lifetime_threshold], y = survival_model(delta = delta(age = age_axis[30:lifetime_threshold], parameters = optimized_parameters_aai_30_exp_quadratic), x = age_of_inf_30_axis), col = "green")
# Cubic model
lines(x = age_axis[30:lifetime_threshold], y = survival_model(delta = delta(age = age_axis[30:lifetime_threshold], parameters = optimized_parameters_aai_30_exp_cubic), x = age_of_inf_30_axis), col = "blue")
# Quartic model
lines(x = age_axis[30:lifetime_threshold], y = survival_model(delta = delta(age = age_axis[30:lifetime_threshold], parameters = optimized_parameters_aai_30_exp_quartic), x = age_of_inf_30_axis), col = "red")
legend(100,1, legend = c("Data", 'Linear model', "Quadratic model", "Cubic model", "4th order polynomial model"), col = c("black", "yellow", "green", "blue", "red"), lty = 1, cex = 0.8)
dev.off()

## Uninfected population -------------------------------------------------

pdf("uninfected_population_survival_and_fit.pdf")
plot(x = age_axis[16:lifetime_threshold], y = survivals_uninfected[16:lifetime_threshold], xlab = "Age [Days]", ylab = "Survival_Uninfected", type = "s", lwd =2)
# Linear model
lines(x = age_axis[16:lifetime_threshold], y = survival_model(delta = delta(age = age_axis[16:lifetime_threshold], parameters = optimized_parameters_uninfected_linear), x = age_axis)[1:135], col = "yellow")
# Quadratic model
lines(x = age_axis[16:lifetime_threshold], y = survival_model(delta = delta(age = age_axis[16:lifetime_threshold], parameters = optimized_parameters_uninfected_quadratic), x = age_axis)[1:135], col = "green")
# Cubic model
lines(x = age_axis[16:lifetime_threshold], y = survival_model(delta = delta(age = age_axis[16:lifetime_threshold], parameters = optimized_parameters_uninfected_cubic), x = age_axis)[1:135], col = "blue")
# Quartic model
lines(x = age_axis[16:lifetime_threshold], y = survival_model(delta = delta(age = age_axis[16:lifetime_threshold], parameters = optimized_parameters_uninfected_quartic), x = age_axis)[1:135], col = "red")
legend(90,1, legend = c("Data", 'Linear model', "Quadratic model", "Cubic model", "4th order polynomial model"), col = c("black", "yellow", "green", "blue", "red"), lty = 1, cex = 0.8)
dev.off()

## Exposed population ------------------------------------------------------
pdf("exposed_population_survival_and_fit.pdf")
plot(x = age_axis[16:lifetime_threshold], y = survivals_exposed[16:lifetime_threshold], xlab = "Age [Days]", ylab = "Survival_Exposed", type = "s", lwd =2)
# Linear modelsurvival_model
lines(x = age_axis[16:lifetime_threshold], y = survival_model(delta = delta(age = age_axis[16:lifetime_threshold], parameters = optimized_parameters_exposed_linear), x = age_axis)[1:135], col = "yellow")
# Quadratic model
lines(x = age_axis[16:lifetime_threshold], y = survival_model(delta = delta(age = age_axis[16:lifetime_threshold], parameters = optimized_parameters_exposed_quadratic), x = age_axis)[1:135], col = "green")
# Cubic model
lines(x = age_axis[16:lifetime_threshold], y = survival_model(delta = delta(age = age_axis[16:lifetime_threshold], parameters = optimized_parameters_exposed_cubic), x = age_axis)[1:135], col = "blue")
# Quartic model
#par_trial <- c(9.2e-09, -1.577257e-06,  8.363635e-05, -1.362067e-03,  6.867880e-03)
lines(x = age_axis[16:lifetime_threshold], y = survival_model(delta = delta(age = age_axis[16:lifetime_threshold], parameters = optimized_parameters_exposed_quartic), x = age_axis)[1:135], col = "red")
legend(90,1, legend = c("Data", 'Linear model', "Quadratic model", "Cubic model", "4th order polynomial model"), col = c("black", "yellow", "green", "blue", "red"), lty = 1, cex = 0.8)
dev.off()

#optimized_likelihood_exposed_quartic
#likelihood_general(parameters = par_trial, longevity = longevity_exposed, age_at_infection = 16)
#optim_trial <- optim(par =optim_, fn = likelihood_general, longevity = longevity_exposed, age_at_infection = 16)
#optim_trial$value
# likelihood_general returns -log_likelihood, therefore we are searching for the minimum of the value, optimizer 
# basically succeeded, is our likelihood alright?



##TODO:check control plots and proceed with computing virulence
#### Virulence for different aai groups -----------------------------------------------------------------------------------------------------
virulence_aai_5 <- delta(age = age_axis[5:lifetime_threshold], parameters = optimized_parameters_aai_5_inf_quadratic) - delta(age = age_axis[5:lifetime_threshold], parameters = optimized_parameters_aai_5_exp_quartic)
virulence_aai_15 <- delta(age = age_axis[15:lifetime_threshold], parameters = optimized_parameters_aai_15_inf_quadratic) - delta(age = age_axis[15:lifetime_threshold], parameters = optimized_parameters_aai_15_exp_quartic)
virulence_aai_30 <- delta(age = age_axis[30:lifetime_threshold], parameters = optimized_parameters_aai_30_inf_quadratic) - delta(age = age_axis[30:lifetime_threshold], parameters = optimized_parameters_aai_30_exp_quartic)

pdf("Virulences_age_of_infection_structured_same_model_for_same_infection_status.pdf")
plot(x = age_axis[1:(lifetime_threshold-4)], y = virulence_aai_5, xlab = "Age of infection [Days]", ylab = "Virulence", col = "red", type = "l")
lines(x = age_axis[1:(lifetime_threshold-14)], y = virulence_aai_15, col = "green")
lines(x = age_axis[1:(lifetime_threshold-29)], y = virulence_aai_30, col = "blue")
points(x = c(89, 78, 68), y = c(0.224, 0.095, 0.06), col = "black", pch = 20, cex = 1.3)
legend(0,0.32, bty="n", legend = c("Age at infection = 5", "Age at infection = 15", "Age at infection = 30", "Points mark the last\navailable datapoint for fitting"), col = c("red", "green", "blue", "black"), lty =c(1,1,1,NA) , pch = c(NA,NA,NA,20), cex = 0.7, pt.cex = 1.3)
text(20.2,0.17, labels = "Sample size for exposed:", cex = 0.7)
text(45,0.17, labels = "16", cex = 0.7, col = "red")
text(51,0.17, labels = "56", cex = 0.7, col = "green")
text(57,0.17, labels = "64", cex = 0.7, col = "blue")
text(20.5,0.15, labels = "Sample size for infecteds:", cex = 0.7)
text(45,0.15, labels = "97", cex = 0.7, col = "red")
text(51,0.15, labels = "27", cex = 0.7, col = 'green')
text(57,0.15, labels = "14", cex = 0.7, col = 'blue')
dev.off()

### Plot susceptibility dependence on the age at infection --------------------------------------------------
y_sus <- c(length(longevity_aai_5_inf)/(length(longevity_aai_5_inf)+length(longevity_aai_5_exp)), length(longevity_aai_15_inf)/(length(longevity_aai_15_inf)+length(longevity_aai_15_exp)), length(longevity_aai_30_inf)/(length(longevity_aai_30_inf)+length(longevity_aai_30_exp)))
x_sus <- c(5,15,30)

plot(x_sus, y_sus, xlab = "Age at infection [Days]", ylab = "Ratio of infecteds in exposed population", type = "p", pch = 20, cex = 1.5, col = "blue" )



### Plotting histograms of lifespans for different ages of infection ------------------------------------------------------------
### seriously whats up with the high death rate around the age of 40 for aai_5_inf (does it relate to the step in the survival?)
#require("plotrix")
require("weights")
#longevities
inf_lifespan_aai_5 <- (longevity_aai_5_inf - 5)
inf_lifespan_aai_15 <- (longevity_aai_15_inf - 15)
inf_lifespan_aai_30 <- (longevity_aai_30_inf - 30)

#weights
virulence_proportion_aai_inf_5_aoi <-   virulence_aai_5[1:(max(longevity_aai_5_inf)-5)] / delta(age = age_axis[5:lifetime_threshold], parameters = optimized_parameters_aai_5_inf_quadratic)[1:(max(longevity_aai_5_inf)-5)]
virulence_proportion_aai_inf_15_aoi <-  virulence_aai_15[1:(max(longevity_aai_15_inf)-15)] / delta(age = age_axis[15:lifetime_threshold], parameters = optimized_parameters_aai_15_inf_quadratic)[1:(max(longevity_aai_15_inf)-15)]
virulence_proportion_aai_inf_30_aoi <-  virulence_aai_30[1:(max(longevity_aai_30_inf)-30)] / delta(age = age_axis[30:lifetime_threshold], parameters = optimized_parameters_aai_30_inf_quadratic)[1:(max(longevity_aai_30_inf)-30)]

prepare_weights <- function(proportion, longevity, aai){
  weights <- rep(0, length(longevity))
  for (i in 1:length(longevity)) {
    if(proportion[(longevity[i]-aai)]>0){
      weights[i] <- proportion[(longevity[i]-aai)]
    }else{
      weights[i] <- 0
    }
  }
  return(weights)
}

w_5 <- prepare_weights(proportion = virulence_proportion_aai_inf_5_aoi, longevity = inf_lifespan_aai_5, aai = 5)
w_15 <- prepare_weights(proportion = virulence_proportion_aai_inf_15_aoi, longevity = inf_lifespan_aai_15, aai = 15)
w_30 <- prepare_weights(proportion = virulence_proportion_aai_inf_30_aoi, longevity = inf_lifespan_aai_30, aai = 30)
length(virulence_proportion_aai_inf_30_aoi[inf_lifespan_aai_30[5]-30])

#histograms showing death distribution due too virulence only
#(significance of each death weighted by the ratio of virulence to total deathrate)

### histograms of death distributions for respective aai classes
pdf("Overall_distribution_of_deaths.pdf")
hist(inf_lifespan_aai_5, col = rgb(1,0,0,0.5), freq = FALSE, main = "Overall distribution of deaths", xlab = "Age of infection [Days]")
hist(inf_lifespan_aai_15, col = rgb(0,0,1,0.5), add = T, freq = FALSE)
hist(inf_lifespan_aai_30, col =  rgb(0,1,0,0.5), add = T, freq = FALSE)
legend("topright", bty="n", legend = c("Age at infection = 5", "Age at infection = 15", "Age at infection = 30"), col = c("red", "green", "blue"), lty =c(1,1,1) , cex = 0.9, lwd = 2)
dev.off()

pdf("Distribution_of_deaths_due_to_virulence.pdf")
wtd.hist(x = inf_lifespan_aai_5, weight = w_5, col = rgb(1,0,0,0.5), freq = FALSE, ylim = c(0,0.12), main = "Distributions of deaths due to infection", xlab = "Age of infection [days]", ylab = "Density")
wtd.hist(x = inf_lifespan_aai_15, breaks = seq(30,85, by = 5), weight = w_15, col = rgb(0,1,0,0.5), add = T, freq = FALSE)
wtd.hist(x = inf_lifespan_aai_30, breaks = seq(30,85, by = 5) , weight = w_30,col = rgb(0,0,1,0.5) , add = T, freq = FALSE)
legend(30,0.12, bty="n", legend = c("Age at infection = 5", "Age at infection = 15", "Age at infection = 30"), col = c("red", "green", "blue"), lty =c(1,1,1) , cex = 0.9, lwd = 2)
dev.off()


## Commparing the sample sizes with proper sizing when reflecting the deathrate
sample_sizes_for_aai_compartments <- c(length(longevity_aai_5_exp)+length(longevity_aai_5_inf),length(longevity_aai_15_exp)+length(longevity_aai_15_inf),length(longevity_aai_30_exp)+length(longevity_aai_30_inf))
expected_sample_sizes <- c(113, 113*sum(host_longevity_control>14)/sum(host_longevity_control>4), 113*sum(host_longevity_control>29)/sum(host_longevity_control>4))
expected_sample_sizes
#### TODO: ---------------------------------------------------------------------------------
### is the fit of aai_30_exp satisfactory? 
### Pick correct compartment for the computation of the natural deathrate (wait for answer by Roland)
### Including what we have in our model
### Relation between lifespan distribution and virulence?
### Look for better test then Kolmogorov-Smirnov (wait for Rolands reply)
### Trying to fit with Weibull distribution?
### As stated here: https://indico.cern.ch/event/217511/contributions/444928/attachments/349287/486926/compar.pdf
  # it should be fine using the KS test for cumulative distribution function (even supported in wikipedia, 
  # if 1-F(x) doesnt change the statistics(shouldnt!))
### Continue with improving the fit: Try moving the fit closer to the data manually,
  # if it looks better but likelihood decreases, we messed up in likelihood, if likelihood increases above maximum,
  # try better optimizer (bblme4,optimx)
### We could try fitting with Weibull, maybe each step in the survival individualy (so combination of two Weibul distributions),
  # or combination of two quadratic functions
### Likelihood hopefully corrected, try running the whole optimization and check the control plots
### (perhaps go back to optimizing even over the first 16 years with the uninfected and whole exposed population)
### Consider computing virulence since the time of castration and check that we dont include
  # in the exposed/infected population the individuals who died within 7 days after the exposure 
  # (infection status could not be verified so soon after exposure)