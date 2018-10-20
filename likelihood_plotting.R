#control script, to plot the likelihood and check the optimizer


#Libraries
library(foreign)   #to read the spss file
#library(ggplot2)

#read data
dataset <- read.spss("/home/petr/Documents/Master_thesis/AgeEffectsVirulence.sav", to.data.frame = TRUE, use.value.labels = FALSE)
parasyte <- as.character(dataset$Pasteuria)
infecteds <- as.numeric(dataset$Infected)

lifetime_threshold <- 150
age_axis <- c(0:lifetime_threshold)

#Assuming a linear relationship, we have the following Ansatz
delta <- function(parameters, age){
  local_delta <- 0
  local_delta <- parameters[1]*age + parameters[2]
  return (local_delta)
}


#try for control population
control_population <- dataset[parasyte == "Control ",] # Careful, the value of the variable contains the space after Control, possibly FIX the dataset later
host_longevity_control <- control_population$HostLongevity

likelihood_control <- function(parameters, longevity = host_longevity_control){
  likelihood_local <- 0
   if(max(delta(parameters, age_axis))<1 & min(delta(parameters, age_axis))>0){
     
     for (i in 1:length(longevity)) {
       for (a in 1:(longevity[i]-1)) {
         likelihood_local <- likelihood_local + log(1-delta(parameters, a))
       }
       
       likelihood_local <- likelihood_local + log(delta(parameters, longevity[i]))
     }
     
   }else{
     likelihood_local <- -Inf
   }
  
  return(-likelihood_local)
}

#linear boundaries
alfa_axis <- seq(from = 0.0001, to = 0.01, by = 0.0001)
beta_axis <- seq(from = 0.0001, to = 0.01, by = 0.0001)
likelihood_values <- matrix(nrow = length(alfa_axis), ncol = length(beta_axis))


for(i in 1:length(alfa_axis)){
  for (j in 1:length(beta_axis)) {
    likelihood_values[i,j] <- likelihood(c(alfa_axis[i], beta_axis[j]), longevity = host_longevity_control)
  }
}

max(likelihood_values)
