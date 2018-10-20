#Estimation of death rate of daphnias, depending on their age: delta(age)

#Libraries
library(foreign)   #to read the spss file

#Daphnia cannot live longer than 200 days (in our data not even longer than 123)
lifetime_threshold = 150

#Initialize the variables
age_axis <- c(0:lifetime_threshold)

#Assuming a linear relationship, we have the following Ansatz
delta <- function(parameters, age){
 return (parameters[1]*age + parameters[2])
}

#Let us optimize the regression parameters, to fit our data

#Reading the data
dataset <- read.spss("/home/petr/Documents/Master_thesis/AgeEffectsVirulence.sav", to.data.frame = TRUE, use.value.labels = FALSE)
parasyte <- as.character(dataset$Pasteuria)
infecteds <- as.numeric(dataset$Infected)

#Considering control population, infecteds and exposed (but uninfected) separately

control_population <- dataset[parasyte == "Control ",] # Careful, the value of the variable contains the space after Control, possibly FIX the dataset later
host_longevity_control <- control_population$HostLongevity

infected_population <- dataset[infecteds == 1,]
host_longevity_infecteds <-infected_population$HostLongevity

exposed_population <- dataset[(infecteds == 0) & (parasyte != "Control "),]
host_longevity_exposed <- exposed_population$HostLongevity


#Survivals of the control, infected and exposed populations
survivals_control <- c(rep(0,lifetime_threshold+1))
survivals_infecteds <- c(rep(0,lifetime_threshold+1))
survivals_exposed <- c(rep(0,lifetime_threshold+1))

for (x in 0:lifetime_threshold) {
  survivals_control[x] <- sum(host_longevity_control>x)/ length(host_longevity_control)
  
}

for (x in 0:lifetime_threshold) {
  survivals_infecteds[x] <- sum(host_longevity_infecteds>x)/ length(host_longevity_infecteds)
  
}

for (x in 0:lifetime_threshold) {
  survivals_exposed[x] <- sum(host_longevity_exposed>x)/ length(host_longevity_exposed)
  
}


#Likelihood function (log_likelihood) = log(Pruduct over individuals (Product over timesteps(probability of dying in that timestep))) -> log turns pruducts to sums
likelihood_control <- function(parameters, longevity = host_longevity_control){
  likelihood_local <- 0
  
  if((lifetime_threshold*parameters[1]+parameters[2])<1 & parameters[1]>0 & parameters[2]>0){
    
    #does the binomial factor play any role in the likelihood? It does not, right? Since under the logarithm it only plays role of additive constant (independent of delta)
    for (i in 1:length(longevity)) {
      for (a in 1:(longevity[i]-1)) {
        likelihood_local <- likelihood_local + log(1-delta(parameters, a))
      }
      
      likelihood_local <- likelihood_local + log(delta(parameters, longevity[i]))
    }

  }else{
    likelihood_local <- -99999999999
  } 
  return(-likelihood_local) #optim is minimizing, so we return sign flipped value to get maximization
}

likelihood_infecteds <- function(parameters, longevity = host_longevity_infecteds){
  likelihood_local <- 0
  
  if((lifetime_threshold*parameters[1]+parameters[2])<1 & parameters[1]>0 & parameters[2]>0){
    
    #does the binomial factor play any role in the likelihood? It does not, right? Since under the logarithm it only plays role of additive constant (independent of delta)
    for (i in 1:length(longevity)) {
      for (a in 1:(longevity[i]-1)) {
        likelihood_local <- likelihood_local + log(1-delta(parameters, a))
      }
      
      likelihood_local <- likelihood_local + log(delta(parameters, longevity[i]))
    }
    
  }else{
    likelihood_local <- -99999999
  } 
  return(-likelihood_local) #optim is minimizing, so we return sign flipped value to get maximization
}

likelihood_exposed <- function(parameters, longevity = host_longevity_exposed){
  likelihood_local <- 0
  
  if((lifetime_threshold*parameters[1]+parameters[2])<1 & parameters[1]>0 & parameters[2]>0){
    
    #does the binomial factor play any role in the likelihood? It does not, right? Since under the logarithm it only plays role of additive constant (independent of delta)
    for (i in 1:length(longevity)) {
      for (a in 1:(longevity[i]-1)) {
        likelihood_local <- likelihood_local + log(1-delta(parameters, a))
      }
      
      likelihood_local <- likelihood_local + log(delta(parameters, longevity[i]))
    }
    
  }else{
    likelihood_local <- -99999999
  } 
  return(-likelihood_local) #optim is minimizing, so we return sign flipped value to get maximization
}

#optimize likelihood w.r.t. alfa and beta
optimized_parameters_control <- optim(c(0.0001,0.001), fn=likelihood_control)$par
delta_values_control <- delta(optimized_parameters_control, age_axis)

optimized_parameters_infecteds <- optim(c(0.0001,0.001), fn=likelihood_infecteds)$par
delta_values_infecteds <- delta(optimized_parameters_infecteds, age_axis)

optimized_parameters_exposed <- optim(c(0.0001,0.001), fn=likelihood_exposed)$par
delta_values_exposed <- delta(optimized_parameters_exposed, age_axis)

#check overlap with survival data

survival_model_control <- function(delta = delta_values_control, x = age_axis){
  exponent <- c(rep(0,length(x)))
  for (a in x) {
    exponent[a+1] <- -sum(delta[1:(a+1)]*x[1:(a+1)])
  }
  return(exp(exponent))
}

survival_model_infecteds <- function(delta = delta_values_infecteds, x = age_axis){
  exponent <- c(rep(0,length(x)))
  for (a in x) {
    exponent[a+1] <- -sum(delta[1:(a+1)]*x[1:(a+1)])
  }
  return(exp(exponent))
}

survival_model_exposed <- function(delta = delta_values_exposed, x = age_axis){
  exponent <- c(rep(0,length(x)))
  for (a in x) {
    exponent[a+1] <- -sum(delta[1:(a+1)]*x[1:(a+1)])
  }
  return(exp(exponent))
}

#parameters for delta, after optimization

optimized_parameters_control
optimized_parameters_infecteds
optimized_parameters_exposed


#plotting the survival according to our liner model of delta

plot_control <- plot(age_axis, survivals_control, col="blue")
lines(age_axis, survival_model_control(), col="red")

plot_infecteds <- plot(age_axis, survivals_infecteds, col="blue")
lines(age_axis, survival_model_infecteds(), col="red")

plot_exposed <- plot(age_axis, survivals_exposed, col="blue")
lines(age_axis, survival_model_exposed(), col="red")

#plotting delta, for different populations
plot_deltas <- plot(age_axis, delta(optimized_parameters_control, age_axis), col = "blue", type = "l")
lines(age_axis, delta(optimized_parameters_infecteds, age_axis), col = "red")
lines(age_axis, delta(optimized_parameters_exposed, age_axis), col = "green")
