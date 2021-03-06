#Estimation of death rate of daphnias, depending on their age: delta(age)

#Libraries
library(foreign)   #to read the spss file

#Daphnia cannot live longer than 200 days (in our data not even longer than 123)
lifetime_threshold = 150

#Initialize the variables
age_axis <- c(0:lifetime_threshold)

#Our polynomial Ansatz for delta
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
  }
 return (local_delta)
}

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
  survivals_control[x+1] <- sum(host_longevity_control>x)/ length(host_longevity_control)
  
}

for (x in 0:lifetime_threshold) {
  survivals_infecteds[x+1] <- sum(host_longevity_infecteds>x)/ length(host_longevity_infecteds)
  
}

for (x in 0:lifetime_threshold) {
  survivals_exposed[x+1] <- sum(host_longevity_exposed>x)/ length(host_longevity_exposed)
  
}

#Likelihood function (log_likelihood) = log(Pruduct over individuals (Product over timesteps(probability of dying in that timestep))) -> log turns pruducts to sums
likelihood_control <- function(parameters, longevity = host_longevity_control){
  likelihood_local <- 0
  
      for (i in 1:length(longevity)) {
        for (a in 1:(longevity[i]-1)) {
          likelihood_local <- likelihood_local + log(1-delta(parameters, a))
        }
        
        likelihood_local <- likelihood_local + log(delta(parameters, longevity[i]))
      }

  return(-(likelihood_local + log(factorial(length(longevity))))) #optim is minimizing, so we return sign flipped value to get maximization
}

likelihood_infecteds <- function(parameters, longevity = host_longevity_infecteds){
  likelihood_local <- 0
  for (i in 1:length(longevity)) {
      for (a in 1:(longevity[i]-1)) {
        likelihood_local <- likelihood_local + log(1-delta(parameters, a))
      }
      likelihood_local <- likelihood_local + log(delta(parameters, longevity[i]))
    }
   
  return(-(likelihood_local + log(factorial(length(longevity))))) #optim is minimizing, so we return sign flipped value to get maximization
}

likelihood_exposed <- function(parameters, longevity = host_longevity_exposed){
  likelihood_local <- 0
  for (i in 1:length(longevity)) {
      for (a in 1:(longevity[i]-1)) {
        likelihood_local <- likelihood_local + log(1-delta(parameters, a))
      }
      
      likelihood_local <- likelihood_local + log(delta(parameters, longevity[i]))
    }
    
  return(-(likelihood_local + log(factorial(length(longevity)))))   #we minimize -likelihood => maximize likelihood; the term log(factorial()) is compensating for the binomial factors
}

#optimize likelihood w.r.t. coefficients of the polynomial, in delta
optimized_parameters_control_quartic <- optim(c(0.0000000001, 0.0000001, 0.000001, 0.001, 0.001), fn=likelihood_control)$par 
delta_values_control_quartic <- delta(optimized_parameters_control_quartic, age_axis)
optimal_likelihood_control_quartic <-optim(c(0.0000000001, 0.0000001, 0.000001, 0.001, 0.001), fn=likelihood_control)$value

optimized_parameters_infecteds_quartic <- optim(c(0.0000000001, 0.0000001, 0.000001, 0.001, 0.001), fn=likelihood_infecteds)$par
delta_values_infecteds_quartic <- delta(optimized_parameters_infecteds_quartic, age_axis)
optimal_likelihood_infecteds_quartic <-optim(c(0.0000000001, 0.0000001, 0.000001, 0.001, 0.001), fn=likelihood_infecteds)$value

optimized_parameters_exposed_quartic <- optim(c(0.0000000001, 0.0000001, 0.000001, 0.001,0.001), fn=likelihood_exposed)$par
delta_values_exposed_quartic <- delta(optimized_parameters_exposed_quartic, age_axis)
optimized_parameters_exposed_quartic <- optim(c(0.0000000001, 0.0000001, 0.000001, 0.001,0.001), fn=likelihood_exposed)$value

optimized_parameters_control_cubic <- optim(c(0.0000001, 0.000001, 0.001, 0.001), fn=likelihood_control)$par 
delta_values_control_cubic <- delta(optimized_parameters_control_cubic, age_axis)
optimized_parameters_control_cubic <- optim(c(0.0000001, 0.000001, 0.001, 0.001), fn=likelihood_control)$value

optimized_parameters_infecteds_cubic <- optim(c(0.0000001, 0.000001, 0.001, 0.001), fn=likelihood_infecteds)$par
delta_values_infecteds_cubic <- delta(optimized_parameters_infecteds_cubic, age_axis)
optimized_parameters_infecteds_cubic <- optim(c(0.0000001, 0.000001, 0.001, 0.001), fn=likelihood_infecteds)$value

optimized_parameters_exposed_cubic <- optim(c(0.0000001, 0.000001, 0.001,0.001), fn=likelihood_exposed)$par
delta_values_exposed_cubic <- delta(optimized_parameters_exposed_cubic, age_axis)
optimized_parameters_exposed_cubic <- optim(c(0.0000001, 0.000001, 0.001,0.001), fn=likelihood_exposed)$value

optimized_parameters_control_quadratic <- optim(c(0.000001, 0.001, 0.001), fn=likelihood_control)$par 
delta_values_control_quadratic <- delta(optimized_parameters_control_quadratic, age_axis)
optimized_parameters_control_quadratic <- optim(c(0.000001, 0.001, 0.001), fn=likelihood_control)$value 

optimized_parameters_infecteds_quadratic <- optim(c(0.000001, 0.001, 0.001), fn=likelihood_infecteds)$par
delta_values_infecteds_quadratic <- delta(optimized_parameters_infecteds_quadratic, age_axis)
optimized_parameters_infecteds_quadratic <- optim(c(0.000001, 0.001, 0.001), fn=likelihood_infecteds)$value

optimized_parameters_exposed_quadratic <- optim(c(0.000001, 0.001,0.001), fn=likelihood_exposed)$par
delta_values_exposed_quadratic <- delta(optimized_parameters_exposed_quadratic, age_axis)
optimized_parameters_exposed_quadratic <- optim(c(0.000001, 0.001,0.001), fn=likelihood_exposed)$value

optimized_parameters_control_linear <- optim(c(0.001, 0.001), fn=likelihood_control)$par
delta_values_control_linear <- delta(optimized_parameters_control_linear, age_axis)
optimized_parameters_control_linear <- optim(c(0.001, 0.001), fn=likelihood_control)$value

optimized_parameters_infecteds_linear <- optim(c(0.001, 0.001), fn=likelihood_infecteds)$par
delta_values_infecteds_linear <- delta(optimized_parameters_infecteds_linear, age_axis)
optimized_parameters_infecteds_linear <- optim(c(0.001, 0.001), fn=likelihood_infecteds)$value

optimized_parameters_exposed_linear <- optim(c(0.001,0.001), fn=likelihood_exposed)$par
delta_values_exposed_linear <- delta(optimized_parameters_exposed_linear, age_axis)
optimized_parameters_exposed_linear <- optim(c(0.001,0.001), fn=likelihood_exposed)$value

#check overlap with survival data

survival_model <- function(delta, x = age_axis){
  exponent <- c(rep(0,length(x)))
  for (a in x) {
    exponent[a+1] <- -sum(delta[1:(a+1)])
  }
  return(exp(exponent))
}


#parameters for delta, after optimization

optimized_parameters_control_quartic
optimized_parameters_infecteds_quartic
optimized_parameters_exposed_quartic

optimized_parameters_control_cubic
optimized_parameters_infecteds_cubic
optimized_parameters_exposed_cubic

optimized_parameters_control_quadratic
optimized_parameters_infecteds_quadratic
optimized_parameters_exposed_quadratic

optimized_parameters_control_linear
optimized_parameters_infecteds_linear
optimized_parameters_exposed_linear

#likelihood ratio test for control population
LR_lin_vs_quadratic_control <- likelihood_control(optimized_parameters_control_linear)/likelihood_control(optimized_parameters_control_quadratic)
LR_lin_vs_cubic_control <- likelihood_control(optimized_parameters_control_linear)/likelihood_control(optimized_parameters_control_cubic)
LR_lin_vs_quartic_control <- likelihood_control(optimized_parameters_control_linear)/likelihood_control(optimized_parameters_control_quartic)

LR_lin_vs_quadratic_control
LR_lin_vs_cubic_control
LR_lin_vs_quartic_control

#likelihood ratio test for exposed population
LR_lin_vs_quadratic_exposed <- likelihood_control(optimized_parameters_exposed_linear)/likelihood_control(optimized_parameters_exposed_quadratic)
LR_lin_vs_cubic_exposed <- likelihood_control(optimized_parameters_exposed_linear)/likelihood_control(optimized_parameters_exposed_cubic)
LR_lin_vs_quartic_exposed <- likelihood_control(optimized_parameters_exposed_linear)/likelihood_control(optimized_parameters_exposed_quartic)

LR_lin_vs_quadratic_exposed
LR_lin_vs_cubic_exposed
LR_lin_vs_quartic_exposed


#likelihood ratio test for infected population
LR_lin_vs_quadratic_infecteds <- likelihood_control(optimized_parameters_infecteds_linear)/likelihood_control(optimized_parameters_infecteds_quadratic)
LR_lin_vs_cubic_infecteds <- likelihood_control(optimized_parameters_infecteds_linear)/likelihood_control(optimized_parameters_infecteds_cubic)
LR_lin_vs_quartic_infecteds <- likelihood_control(optimized_parameters_infecteds_linear)/likelihood_control(optimized_parameters_infecteds_quartic)
#LR_lin_vs_quartic_infecteds <- likelihood_control(optimized_parameters_infecteds_linear)/likelihood_control(c(0,0,optimized_parameters_infecteds_quartic[3:5]))

LR_lin_vs_quadratic_infecteds
LR_lin_vs_cubic_infecteds
LR_lin_vs_quartic_infecteds

#from the LR test, it seems the linear model is the best for control and exposed populations and the cubic model is the best for the infected population


#plotting the survival compared to linear and quadratic model of delta

pdf("Control_population_survival.pdf")
plot_control <- plot(age_axis, survivals_control, col="blue", xlab = "Age [days]", ylab = "Survivals in control population [fraction]")
lines(age_axis, survival_model(delta = delta_values_control_linear), col="red")
lines(age_axis, survival_model(delta = delta_values_control_quadratic), col="green")
lines(age_axis, survival_model(delta = delta_values_control_cubic), col="yellow")
lines(age_axis, survival_model(delta = delta_values_control_quartic), col="orange")
legend(0,0.3, legend = c("Data", 'Linear model', "Quadratic model", "Cubic model", "4th order polynomial model"), col = c("blue", "red", "green", "yellow", "orange"), lty = 1)
dev.off()

pdf("Infected_population_survival.pdf")
plot_infecteds <- plot(age_axis, survivals_infecteds, col="blue", xlab = "Age [days]", ylab = "Survivals in infected population [fraction]")
lines(age_axis, survival_model(delta = delta_values_infecteds_linear), col="red")
lines(age_axis, survival_model(delta = delta_values_infecteds_quadratic), col="green")
lines(age_axis, survival_model(delta = delta_values_infecteds_cubic), col="yellow")
lines(age_axis, survival_model(delta = optimized_parameters_infecteds_quartic), col="orange")
legend(70,1, legend = c("Data", 'Linear model', "Quadratic model", "Cubic model", "4th order polynomial model"), col = c("blue", "red", "green", "yellow", "orange"), lty = 1)
dev.off()


pdf("Exposed_population_survival.pdf")
plot_exposed <- plot(age_axis, survivals_exposed, col="blue", xlab = "Age [days]", ylab = "Survivals in exposed population [fraction]")
lines(age_axis, survival_model(delta = delta_values_exposed_linear), col="red")
lines(age_axis, survival_model(delta = delta_values_exposed_quadratic), col="green")
lines(age_axis, survival_model(delta = delta_values_exposed_cubic), col="yellow")
lines(age_axis, survival_model(delta = delta_values_exposed_quartic), col="orange")
legend(0,0.3, legend = c("Data", 'Linear model', "Quadratic model", "Cubic model", "4th order polynomial model"), col = c("blue", "red", "green", "yellow", "orange"), lty = 1)
dev.off()

#plotting delta, for different populations

pdf("Quadratic_models.pdf")
plot_deltas <- plot(age_axis, delta(optimized_parameters_control_quadratic, age_axis), xlab = "Age [days]", ylab = "Delta under quadratic approx.", col = "blue", type = "l")
lines(age_axis, delta(optimized_parameters_infecteds_quadratic, age_axis), col = "red")
lines(age_axis, delta(optimized_parameters_exposed_quadratic, age_axis), col = "green")
legend(100,0.02, legend = c("Control", 'Infected', "Exposed"), col = c("blue", "red", "green"), lty = 1)
dev.off()

pdf("Linear_models.pdf")
plot_deltas <- plot(age_axis, delta(optimized_parameters_control_linear, age_axis), xlab = "Age [days]", ylab = "Delta under linear approx.", col = "blue", type = "l")
lines(age_axis, delta(optimized_parameters_infecteds_linear, age_axis), col = "red")
lines(age_axis, delta(optimized_parameters_exposed_linear, age_axis), col = "green")
legend(100,0.01, legend = c("Control", 'Infected', "Exposed"), col = c("blue", "red", "green"), lty = 1)
dev.off()

pdf("Cubic_models.pdf")
plot_deltas <- plot(age_axis, delta(optimized_parameters_control_cubic, age_axis), xlab = "Age [days]", ylab = "Delta under cubic approx.", col = "blue", type = "l")
lines(age_axis, delta(optimized_parameters_infecteds_cubic, age_axis), col = "red")
lines(age_axis, delta(optimized_parameters_exposed_cubic, age_axis), col = "green")
legend(0,0.30, legend = c("Control", 'Infected', "Exposed"), col = c("blue", "red", "green"), lty = 1)
dev.off()